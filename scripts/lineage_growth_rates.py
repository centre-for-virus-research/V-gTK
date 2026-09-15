#!/usr/bin/env python3
"""Estimate lineage growth rates for a protein's alleles in a V-gTK database.

A bolt-on: it reads a finished database - metadata, alignments, mutation calls
and a stored UShER tree - and writes a set of tables and a report. It changes
nothing unless `--write_db` is given.

What it estimates
-----------------

**Without the tree** (``--no-tree`` runs only this half). A tree is not needed
to measure a change in allele or lineage frequency, only dates:

* individual-level logistic regression of carrier frequency on time - the slope
  is the growth-rate advantage, a.k.a. the selection coefficient per year;
* the same fit on binned counts with a quasi-binomial dispersion, which is the
  honest version once you accept that genomic surveillance arrives in
  correlated batches;
* multinomial logistic regression across every residue competing at a site, so
  the estimates are mutually consistent rather than each against a different
  moving background;
* a Poisson log-linear fit to the counts, with and without a sampling-effort
  offset, separating "growing" from "growing faster than its competitors";
* the Frequency Increment Test, which asks whether the increments exceed drift
  without assuming the trajectory is logistic at all;
* Mann-Kendall with Sen's slope, and an early-versus-late Fisher test, as the
  assumption-free fallbacks.

**With the tree**, which adds the ancestry those cannot supply:

* parsimony reconstruction of how many times each allele arose independently,
  reported as the range over all most-parsimonious histories;
* per-origin growth rates, so a rise driven by one lucky introduction is
  distinguishable from the same rise repeated across many origins;
* for each origin: lineages-through-time slope, the Magallon-Sanderson
  method-of-moments rate over a range of assumed extinction fractions, the
  maximum-likelihood pure-birth rate, Pybus & Harvey's gamma, a coalescent
  skyline and the maximum-likelihood exponential-growth coalescent rate;
* the local branching index (Neher, Russell & Shraiman 2014), the standard
  tree-shape fitness proxy;
* a scan of every clade in the tree for logistic frequency growth, FDR
  corrected, which is where the fastest-growing lineages come from whether or
  not they carry a catalogued mutation.

Usage
-----
    # everything the database already knows about NS5A
    python scripts/lineage_growth_rates.py \\
        --db test_out/HCV_XML_full_plus_update/HCV_full_right_tree.db \\
        --protein NS5A --out dev/growth/NS5A

    # one genotype, genotypes taken from the curated reference list via the tree
    python scripts/lineage_growth_rates.py --db DB --protein NS5A \\
        --genotype 1 --genotype-reflist generic/hcv/ref_list_subtype_genotype.txt \\
        --out dev/growth/NS5A_g1

    # no tree, just the frequency time series
    python scripts/lineage_growth_rates.py --db DB --protein G --no-tree --out dev/growth/G
"""

import argparse
import json
import os
import sqlite3
import sys
import time

import numpy as np
import pandas as pd

SCRIPTS_DIR = os.path.dirname(os.path.abspath(__file__))
if SCRIPTS_DIR not in sys.path:
	sys.path.insert(0, SCRIPTS_DIR)

import growth_stats as gs  # noqa: E402
import protein_alleles as pal  # noqa: E402
import tree_growth_metrics as tgm  # noqa: E402

#: Extinction fractions the method-of-moments rate is reported at. mu/lambda
#: cannot be estimated from a clade's size and age, so the convention is to show
#: the rate under a low, a middling and a high assumption rather than pick one.
EPSILON_GRID = (0.0, 0.5, 0.9)

#: Default local-branching-index scale, as a fraction of the cohort's sampling
#: span. Neher et al. tune tau to the timescale over which fitness is being
#: predicted; a sixteenth of the span is the usual order of magnitude and is
#: exposed as --lbi-tau for anything else.
DEFAULT_LBI_TAU_FRACTION = 1.0 / 16.0

#: A calibrated protein start below this agreement with the catalogue calls is
#: refused. It is a deliberately high bar: genotyping a whole database off a
#: frame that is right three times in four would put a made-up residue on a
#: quarter of the sequences, and every growth rate downstream would be measuring
#: that instead of the virus.
MIN_CALIBRATION_AGREEMENT = 0.9

TREE_PREFERENCE = ('usher', 'iqtree')


def say(message):
	print('[%s] %s' % (time.strftime('%H:%M:%S'), message), flush=True)


# ---------------------------------------------------------------------------
# Inputs
# ---------------------------------------------------------------------------

def open_database(path):
	if not os.path.isfile(path):
		raise SystemExit('database not found: %s' % path)
	conn = sqlite3.connect('file:%s?mode=ro' % path, uri=True)
	conn.execute('PRAGMA query_only = ON')
	return conn


def choose_tree(conn, requested=None, segment=None):
	"""Pick the tree to work on. Returns ``(name, source, newick)``.

	UShER is preferred over IQ-TREE for the same reason clade assignment
	prefers it: the UShER tree holds every placed sample, while the IQ-TREE
	backbone only holds cluster representatives, and a growth rate measured on
	representatives is measuring the clustering.

	The segment is applied *before* either preference. A segmented build holds
	one tree per segment, and picking the longest newick instead hands a
	neuraminidase cohort the haemagglutinin tree. That happens to fail loudly
	when the two share no tips at all - but on a build where segments share
	accession names it would silently fit the wrong tree and report the result
	as if nothing were wrong.
	"""
	if not pal.table_exists(conn, 'trees'):
		raise SystemExit('this database has no trees table; re-run with --no-tree')
	rows = [row for row in conn.execute('SELECT name, source, newick, segment FROM trees')
			if row[2]]
	if not rows:
		raise SystemExit('the trees table is empty; re-run with --no-tree')
	if requested:
		for row in rows:
			if row[0] == requested:
				return row[0], row[1], row[2]
		raise SystemExit('no tree named %r (available: %s)'
						 % (requested, ', '.join(sorted(row[0] for row in rows))))

	if segment is not None:
		matching = [row for row in rows
					if str(row[3] or '').strip() == str(segment).strip()]
		if matching:
			rows = matching
		else:
			say('[warn] no tree in this database is labelled segment %s; falling back to '
				'the unlabelled trees, which may belong to another segment' % segment)

	for preferred in TREE_PREFERENCE:
		candidates = [row for row in rows if str(row[1]).strip().lower() == preferred]
		if candidates:
			# Within one segment and source, the longest newick is the most
			# complete tree.
			best = max(candidates, key=lambda row: len(row[2]))
			return best[0], best[1], best[2]
	best = max(rows, key=lambda row: len(row[2]))
	return best[0], best[1], best[2]


def apply_genotypes(conn, cohort, args, newick):
	"""Attach genotype labels, from the stored column or from the reference list."""
	source = 'meta_data.%s' % cohort.attrs.get('genotype_column') \
		if cohort.attrs.get('genotype_column') else None
	if args.genotype_reflist:
		if newick is None:
			raise SystemExit('--genotype-reflist needs a tree; drop --no-tree or '
							 'use a database with genotype columns')
		say('assigning genotypes from %s via the tree' % args.genotype_reflist)
		labels = pal.genotype_labels_from_reference_list(conn, args.genotype_reflist, newick)
		cohort['genotype'] = cohort['primary_accession'].map(
			lambda accession: labels.get(accession, {}).get('genotype', ''))
		cohort['subtype'] = cohort['primary_accession'].map(
			lambda accession: labels.get(accession, {}).get('subtype', ''))
		source = 'reference list %s via tree' % os.path.basename(args.genotype_reflist)
		resolved = int((cohort['genotype'].astype(str).str.strip() != '').sum())
		say('  %d of %d cohort sequences got a genotype' % (resolved, len(cohort)))
		if len(cohort) and resolved < 0.5 * len(cohort):
			say('  [warn] most sequences came back unlabelled. Tree-based clade assignment '
				'refuses to call a genotype when the enclosing clade spans more than half '
				'the tree, and on a heavily polytomous UShER tree that is the usual '
				'outcome - a query hangs off the root polytomy and has no local '
				'neighbourhood. Try --tree with the IQ-TREE backbone, which is properly '
				'resolved, or a database that already carries genotype columns.')
	cohort.attrs['genotype_source'] = source or 'none'

	if args.genotype:
		wanted = {value.strip() for value in args.genotype.split(',') if value.strip()}
		if source is None:
			raise SystemExit(
				'--genotype was given but this database carries no genotype column. '
				'Supply --genotype-reflist <curated reference list> to derive genotypes '
				'from the tree, or --genotype-column to name the column.')
		level = 'subtype' if args.genotype_level == 'subtype' else 'genotype'
		before = len(cohort)
		cohort = cohort[cohort[level].isin(wanted)].reset_index(drop=True)
		say('genotype filter %s in {%s}: %d of %d sequences'
			% (level, ', '.join(sorted(wanted)), len(cohort), before))
		if cohort.empty:
			raise SystemExit('no sequence matched the requested genotype')
	return cohort


# ---------------------------------------------------------------------------
# Alleles
# ---------------------------------------------------------------------------

class AlleleSet(object):
	"""Which sequences carry which residue at which site, and who was tested.

	`carriers` maps an allele id to a set of accessions. `evaluated` maps a site
	to the set of accessions that could have carried it - the denominator.
	Keeping the two apart is the difference between "3% of sequences carry this"
	and "3% of the sequences anyone looked at carry this".
	"""

	def __init__(self, protein, source, carriers, evaluated, residues, notes=None):
		self.protein = protein
		self.source = source
		self.carriers = carriers
		self.evaluated = evaluated
		self.residues = residues
		self.notes = notes or {}
		# Built once: scanning every allele for each site in turn is quadratic,
		# and a whole protein genotyped from the alignment has a few thousand
		# alleles over a few hundred sites.
		self._by_site = {}
		for allele, meta in residues.items():
			self._by_site.setdefault(meta['site'], []).append(allele)
		for alleles in self._by_site.values():
			alleles.sort()

	def sites(self):
		return sorted(self.evaluated)

	def alleles_at(self, site):
		return self._by_site.get(site, [])


def _allele_id(protein, site, residue):
	return '%s:%d%s' % (protein, int(site), residue)


def alleles_from_catalog(conn, protein, cohort_accessions):
	calls = pal.catalog_alleles(conn, protein)
	if calls.empty:
		raise SystemExit('no calls for protein %r in sequence_mutations' % protein)
	calls = calls[calls['primary_accession'].isin(cohort_accessions)]
	if calls.empty:
		raise SystemExit('no sequence in the cohort has a %r call; widen the date or '
						 'genotype filters' % protein)
	tested = set(calls['primary_accession'])
	carriers, residues, evaluated = {}, {}, {}
	for (site, residue), group in calls.groupby(['aa_position', 'alt_residue']):
		if not residue or residue == pal.UNKNOWN_RESIDUE:
			continue
		site = int(site)
		allele = _allele_id(protein, site, residue)
		carriers[allele] = set(group['primary_accession'])
		residues[allele] = {'site': site, 'residue': residue}
		evaluated.setdefault(site, set(tested))
	notes = {
		'denominator': 'sequences with at least one catalogued call for this protein',
		'caveat': ('sequence_mutations records the residues the annotator called, not '
				   'every sequence it looked at, so the denominator is the set of '
				   'sequences that produced any call for this protein (%d here). '
				   'Frequencies are relative to that set.' % len(tested)),
	}
	return AlleleSet(protein, 'catalog', carriers, evaluated, residues, notes)


def alleles_from_alignment(conn, protein, cohort_accessions, start, positions, master=None,
						   segment=None, coord_space='column'):
	frame, extraction = pal.alignment_residues(
		conn, positions, start, accessions=cohort_accessions, master=master,
		segment=segment, coord_space=coord_space)
	say('  read %d alignment rows for segment %s, used %d, skipped %d for a width other '
		'than the master\'s (%d)'
		% (extraction['rows_read'], segment, extraction['rows_used'],
		   extraction['rows_wrong_width'], extraction['master_width']))
	if frame.empty:
		raise SystemExit(
			'no alignment rows matched the cohort (master %s, segment %r, %d rows read, '
			'%d skipped for a different width). If the build stores one master per '
			'segment, pass --segment.'
			% (extraction['master'], segment, extraction['rows_read'],
			   extraction['rows_wrong_width']))
	carriers, residues, evaluated = {}, {}, {}
	for site in [column for column in frame.columns if column != 'primary_accession']:
		called = frame[['primary_accession', site]].copy()
		called.columns = ['primary_accession', 'residue']
		informative = called[~called['residue'].isin({pal.UNKNOWN_RESIDUE, ''})]
		if informative.empty:
			continue
		evaluated[int(site)] = set(informative['primary_accession'])
		for residue, group in informative.groupby('residue'):
			allele = _allele_id(protein, int(site), residue)
			carriers[allele] = set(group['primary_accession'])
			residues[allele] = {'site': int(site), 'residue': str(residue)}
	notes = {'denominator': 'every cohort sequence with a translatable codon at the site',
			 'protein_start': start, 'extraction': extraction}
	notes.update(frame_sanity(carriers, evaluated, residues,
							  skip_sites=extraction.get('master_gap_positions', ()),
							  terminal_site=extraction.get('last_position')))
	return AlleleSet(protein, 'alignment', carriers, evaluated, residues, notes)


#: Internal stop codons are what diagnose a wrong reading frame: a coding region
#: read correctly has essentially none, and one read in the wrong frame has them
#: every few dozen codons.
MAX_STOP_FRACTION = 0.01

#: Whole-codon deletions are *reported* above this fraction but are deliberately
#: NOT treated as a frame error. A protein aligned across subtypes has genuine
#: indel columns - influenza HA subtypes differ by insertions and deletions
#: against any single master, so an H3 sequence read against an H1 master really
#: is missing codons the master has. Treating that as a broken frame made the
#: check fire on every correctly-read influenza protein, which is worse than not
#: having the check at all.
MAX_DELETION_FRACTION = 0.25


def frame_sanity(carriers, evaluated, residues, skip_sites=(), terminal_site=None):
	"""Check a translated region looks like protein, and say so when it does not.

	This exists because `--protein-start` is an override the user is allowed to
	give and nothing else can verify. Getting it wrong does not fail - it
	silently produces hundreds of alleles with tiny p-values, every one of them
	describing a reading frame rather than a virus. Internal stop codons are the
	cheapest possible evidence that has happened.
	"""
	worst_stop = worst_deletion = 0.0
	stop_site = deletion_site = None
	skip_sites = set(skip_sites or ())
	if terminal_site is not None:
		# The last codon of a coding sequence is its own stop codon.
		skip_sites.add(int(terminal_site))
	by_site = {}
	for allele, meta in residues.items():
		by_site.setdefault(meta['site'], []).append(allele)
	for site, tested in evaluated.items():
		total = len(tested)
		if not total:
			continue
		# The last codon of a CDS is its stop codon, and a column where the
		# master is itself a gap is an insertion column. Counting either as
		# evidence of a broken frame makes the check fire on every correctly
		# read protein, which is worse than not having it.
		if site in skip_sites:
			continue
		for allele in by_site.get(site, ()):
			meta = residues[allele]
			fraction = len(carriers.get(allele, ())) / total
			if meta['residue'] == '*' and fraction > worst_stop:
				worst_stop, stop_site = fraction, site
			if meta['residue'] == pal.DELETION_RESIDUE and fraction > worst_deletion:
				worst_deletion, deletion_site = fraction, site
	notes = {'max_stop_fraction': worst_stop, 'max_deletion_fraction': worst_deletion,
			 'max_deletion_site': deletion_site, 'max_stop_site': stop_site}
	if worst_stop > MAX_STOP_FRACTION:
		notes['frame_warning'] = (
			'The translated region does not look like coding sequence: %.1f%% of '
			'sequences carry a stop codon at residue %s. The protein start is very '
			'probably wrong, and every allele frequency below describes a reading frame '
			'rather than the protein.' % (100 * worst_stop, stop_site))
		say('[warn] %s' % notes['frame_warning'])
	if worst_deletion > MAX_DELETION_FRACTION:
		notes['indel_note'] = (
			'%.1f%% of sequences carry a whole-codon deletion at residue %s. Across '
			'subtypes that is an indel column rather than a broken frame, so it is '
			'reported and not treated as an error - but alleles at such a site describe '
			'the alignment, not a substitution, and should be read as such.'
			% (100 * worst_deletion, deletion_site))
		say('[note] %s' % notes['indel_note'])
	return notes


def alleles_from_features(conn, protein, cohort_accessions, master=None, segment=None,
						  coord_space='auto', positions=None):
	"""Alignment route for a protein the database annotates as its own product.

	The master's own feature row gives the CDS bounds, so no calibration is
	needed and the coordinates are the ones the database already uses. What the
	feature row does *not* say is which of the two readings of those
	coordinates applies, so that is detected from the master's own translation
	and reported.
	"""
	master = master or pal.master_accession(conn, segment=segment)
	if master is None:
		raise SystemExit('no master reference for segment %r; this build has %d masters, '
						 'so --segment is required' % (segment, pal.count_masters(conn)))
	# A product carrying more than one feature row on the master is spliced: its
	# CDS is two or more exons, as influenza M2 and NEP both are. One
	# cds_start..cds_end span cannot represent that, and taking the first row
	# would silently translate whichever exon came back first - 26 nucleotides
	# of M2 exon 1, in the influenza build.
	spans = conn.execute('SELECT cds_start, cds_end FROM features WHERE accession = ? AND '
						 'LOWER(TRIM(product)) = LOWER(TRIM(?)) ORDER BY cds_start',
						 (master, protein)).fetchall()
	if len(spans) > 1:
		raise SystemExit(
			'%r is annotated on master %s as %d separate feature rows (%s): a spliced '
			'CDS. Reading it as one contiguous span would translate something that is '
			'not the protein, so it is refused here. Pick an unspliced product, or use '
			'--allele-source alignment with --protein-start and --positions for the '
			'exon you want.'
			% (protein, master, len(spans),
			   '; '.join('%s..%s' % (a, b) for a, b in spans)))

	start = pal.feature_cds_start(conn, protein, master)
	if start is None:
		raise SystemExit('the master (%s) has no %r feature to take coordinates from'
						 % (master, protein))
	row = conn.execute('SELECT cds_end FROM features WHERE accession = ? AND '
					   'LOWER(TRIM(product)) = LOWER(TRIM(?)) LIMIT 1',
					   (master, protein)).fetchone()
	try:
		end = int(float(row[0]))
	except (TypeError, ValueError, IndexError):
		raise SystemExit('the %r feature on %s has no usable end coordinate' % (protein, master))
	n_sites = max(0, (end - start + 1) // 3)
	if n_sites < 1:
		raise SystemExit('the %r feature on %s is shorter than one codon' % (protein, master))

	space, scores = pal.resolve_coord_space(conn, master, start, end, segment=segment,
											requested=coord_space)
	if scores:
		for name, entry in sorted(scores.items()):
			say('  coordinate reading %-9s stops=%d unresolved=%d  %s'
				% (name, entry['internal_stops'], entry['missing'] + entry['unknown'],
				   entry['protein_head']))
		say('  reading CDS %d..%d on master %s as %r' % (start, end, master, space))

	if positions is None:
		positions = range(1, n_sites + 1)
	allele_set = alleles_from_alignment(conn, protein, cohort_accessions, start, positions,
										master=master, segment=segment, coord_space=space)
	allele_set.notes['coord_space_scores'] = scores
	allele_set.notes['cds_bounds'] = (start, end)
	return allele_set


def build_alleles(conn, args, cohort):
	"""Resolve the requested protein into alleles by whichever route applies."""
	accessions = set(cohort['primary_accession'])
	catalogue_proteins = set(pal.list_proteins(conn).query("source == 'catalog'")['protein'])
	feature_products = {str(row[0]).strip().lower() for row
						in conn.execute('SELECT DISTINCT product FROM features')} \
		if pal.table_exists(conn, 'features') else set()

	source = args.allele_source
	if source == 'auto':
		if args.protein.strip().lower() in feature_products:
			source = 'features'
		elif args.protein in catalogue_proteins:
			source = 'catalog'
		else:
			raise SystemExit(
				'protein %r is neither an annotated product nor catalogued in this '
				'database. Known: %s' % (args.protein, ', '.join(sorted(
					catalogue_proteins | {p for p in feature_products}))))
		say('allele source resolved to %r' % source)

	masters = pal.count_masters(conn)
	if masters > 1 and not args.segment and source in ('features', 'alignment'):
		raise SystemExit(
			'this build has %d masters (one per segment) and --segment was not given. '
			'Coordinates resolved against the wrong segment\'s master produce a residue '
			'for every sequence and every one is wrong, so pick a segment explicitly.'
			% masters)

	requested_positions = None
	if args.positions:
		requested_positions = sorted({int(value) for value in args.positions.split(',')
									  if value.strip()})

	if source == 'catalog':
		return alleles_from_catalog(conn, args.protein, accessions)
	if source == 'features':
		return alleles_from_features(conn, args.protein, accessions, segment=args.segment,
									 coord_space=args.coord_space,
									 positions=requested_positions)

	start = args.protein_start
	calibration = None
	if start is None:
		if not args.calibrate_start:
			raise SystemExit('--allele-source alignment needs --protein-start, or '
							 '--calibrate-start to infer it from the catalogue calls')
		say('calibrating the start of %s against the catalogue calls' % args.protein)
		calibration = pal.calibrate_protein_start(conn, args.protein)
		say('  best start %s, agreement %.3f over %d scored calls (runner-up %.3f)'
			% (calibration['start'], calibration['agreement'], calibration['n_scored'],
			   calibration['runner_up_agreement']))
		if (calibration['start'] is None
				or calibration['agreement'] < MIN_CALIBRATION_AGREEMENT):
			raise SystemExit(
				'calibration agreed with only %.1f%% of the catalogue calls, below the '
				'%.0f%% required. The protein\'s numbering is probably not a fixed offset '
				'in this reference (genotype-relative numbering and indels both do this). '
				'Supply --protein-start explicitly, or use --allele-source catalog.'
				% (100 * calibration['agreement'], 100 * MIN_CALIBRATION_AGREEMENT))
		start = calibration['start']

	if args.positions:
		positions = sorted({int(value) for value in args.positions.split(',') if value.strip()})
	else:
		catalogued = pal.catalog_alleles(conn, args.protein)
		if catalogued.empty:
			raise SystemExit('no --positions given and no catalogued sites to fall back on')
		positions = sorted({int(value) for value in catalogued['aa_position']})
		say('genotyping the %d catalogued sites of %s' % (len(positions), args.protein))
	allele_set = alleles_from_alignment(conn, args.protein, accessions, start, positions,
										segment=args.segment,
										coord_space=(args.coord_space
													 if args.coord_space != 'auto' else 'column'))
	if calibration:
		allele_set.notes['calibration'] = calibration
	return allele_set


def site_labels(allele_set, site):
	"""``{accession: residue}`` at one site, dropping anything ambiguous.

	A sequence recorded as carrying two different residues at the same position
	is a mixed or miscalled base, not a member of two lineages, and it cannot be
	given a single label in a multinomial model.
	"""
	labels, ambiguous = {}, set()
	for allele in allele_set.alleles_at(site):
		residue = allele_set.residues[allele]['residue']
		for accession in allele_set.carriers[allele]:
			if labels.get(accession, residue) != residue:
				ambiguous.add(accession)
			labels[accession] = residue
	for accession in ambiguous:
		labels.pop(accession, None)
	return labels


# ---------------------------------------------------------------------------
# Tree-free estimators
# ---------------------------------------------------------------------------

def treefree_metrics(carrier_times, background_times, args, split_time=None):
	"""Every estimator that needs only dates and a carrier flag."""
	carrier_times = np.asarray(carrier_times, dtype=float)
	background_times = np.asarray(background_times, dtype=float)
	times = np.concatenate([carrier_times, background_times])
	flags = np.concatenate([np.ones_like(carrier_times), np.zeros_like(background_times)])

	individual = gs.logistic_growth(times, flags, firth='auto')
	row = {'n_evaluated': int(len(times)), 'n_carrier': int(len(carrier_times)),
		   'frequency': float(len(carrier_times)) / len(times) if len(times) else float('nan'),
		   'first_carrier_date': float(np.min(carrier_times)) if len(carrier_times) else float('nan'),
		   'last_carrier_date': float(np.max(carrier_times)) if len(carrier_times) else float('nan')}
	for key in ('growth_advantage_per_year', 'std_error', 'ci_low', 'ci_high',
				'wald_p', 'lrt_p', 'method', 'converged', 'doubling_time_years',
				'fitted_freq_start', 'fitted_freq_end'):
		row['logistic_' + key] = individual[key]
	row['selection_coefficient_per_generation'] = (
		individual['growth_advantage_per_year'] * args.generation_time)
	row['relative_reproduction_number'] = gs.reproduction_number(
		individual['growth_advantage_per_year'], args.generation_time, args.generation_sd)

	bins = gs.bin_trajectory(times, flags, args.bin_width, origin=args.bin_origin)
	usable = [entry for entry in bins if entry['n'] >= args.min_bin_size]
	row['n_time_bins'] = len(usable)
	if len(usable) >= 3:
		mid = np.array([entry['bin_mid'] for entry in usable])
		total = np.array([entry['n'] for entry in usable], dtype=float)
		positive = np.array([entry['n_carrier'] for entry in usable], dtype=float)
		binned = gs.logistic_growth(mid, positive / total, trials=total, quasi=True)
		for key in ('growth_advantage_per_year', 'std_error', 'ci_low', 'ci_high',
					'wald_p', 'dispersion', 'overdispersed'):
			row['binned_' + key] = binned[key]
		absolute = gs.poisson_growth(mid, positive)
		relative = gs.poisson_growth(mid, positive, exposure=total)
		row['count_growth_per_year'] = absolute['growth_rate_per_year']
		row['count_growth_p'] = absolute['p_value']
		row['count_doubling_time_years'] = absolute['doubling_time_years']
		row['count_growth_vs_sampling_per_year'] = relative['growth_rate_per_year']
		row['count_growth_vs_sampling_p'] = relative['p_value']
		fit = gs.frequency_increment_test(mid, positive / total, total)
		row['fit_statistic'] = fit['fit_statistic']
		row['fit_p'] = fit['p_value']
		trend = gs.mann_kendall(mid, positive / total)
		row['mk_tau'] = trend['tau']
		row['mk_p'] = trend['p_value']
		row['sen_slope_per_year'] = trend['sen_slope']
	else:
		for key in ('binned_growth_advantage_per_year', 'binned_std_error', 'binned_ci_low',
					'binned_ci_high', 'binned_wald_p', 'binned_dispersion',
					'count_growth_per_year', 'count_growth_p', 'count_doubling_time_years',
					'count_growth_vs_sampling_per_year', 'count_growth_vs_sampling_p',
					'fit_statistic', 'fit_p', 'mk_tau', 'mk_p', 'sen_slope_per_year'):
			row[key] = float('nan')
		row['binned_overdispersed'] = False

	if split_time is None:
		split_time = float(np.median(times)) if len(times) else float('nan')
	if np.isfinite(split_time):
		late_carrier = int(np.sum(carrier_times >= split_time))
		early_carrier = int(len(carrier_times) - late_carrier)
		late_total = int(np.sum(times >= split_time))
		early_total = int(len(times) - late_total)
		odds, p_value = gs.fisher_exact_2x2(late_carrier, late_total - late_carrier,
											early_carrier, early_total - early_carrier)
		row['late_vs_early_odds_ratio'] = odds
		row['late_vs_early_p'] = p_value
		row['split_date'] = split_time
	return row, bins


def analyse_alleles(allele_set, cohort, args):
	"""Tree-free growth estimates for every allele, plus the per-site multinomial."""
	dates = dict(zip(cohort['primary_accession'], cohort['decimal_date']))
	split = float(np.median(cohort['decimal_date'])) if len(cohort) else float('nan')
	rows, trajectories, multinomial = [], [], []

	for site in allele_set.sites():
		evaluated = [a for a in allele_set.evaluated[site] if a in dates]
		if len(evaluated) < args.min_allele_count:
			continue
		evaluated_times = np.array([dates[a] for a in evaluated])
		for allele in allele_set.alleles_at(site):
			carriers = [a for a in allele_set.carriers[allele] if a in dates]
			if len(carriers) < args.min_allele_count:
				continue
			carrier_set = set(carriers)
			carrier_times = np.array([dates[a] for a in carriers])
			background = np.array([dates[a] for a in evaluated if a not in carrier_set])
			if len(background) < args.min_allele_count:
				continue
			row, bins = treefree_metrics(carrier_times, background, args, split_time=split)
			row.update({'allele': allele, 'protein': allele_set.protein, 'site': site,
						'residue': allele_set.residues[allele]['residue'],
						'allele_source': allele_set.source})
			rows.append(row)
			for entry in bins:
				trajectories.append({'allele': allele, 'site': site, **entry})

		labels = site_labels(allele_set, site)
		labelled = [(dates[a], residue) for a, residue in labels.items() if a in dates]
		if len(labelled) >= 4 and len({residue for _, residue in labelled}) >= 2:
			records, info = gs.multinomial_growth([t for t, _ in labelled],
												  [r for _, r in labelled])
			for record in records:
				record.update({'protein': allele_set.protein, 'site': site,
							   'pivot': info['pivot'], 'converged': info['converged'],
							   'allele': _allele_id(allele_set.protein, site, record['category'])})
				multinomial.append(record)

	allele_frame = pd.DataFrame(rows)
	if not allele_frame.empty:
		for source, target in (('logistic_lrt_p', 'logistic_q'),
							   ('binned_wald_p', 'binned_q'),
							   ('fit_p', 'fit_q')):
			if source in allele_frame.columns:
				allele_frame[target] = gs.benjamini_hochberg(allele_frame[source])
		allele_frame = allele_frame.sort_values('logistic_growth_advantage_per_year',
												ascending=False).reset_index(drop=True)
	return allele_frame, pd.DataFrame(trajectories), pd.DataFrame(multinomial)


# ---------------------------------------------------------------------------
# Tree-based estimators
# ---------------------------------------------------------------------------

def clade_phylodynamics(tree, node, node_times, branch_times, label='', skyline_group=10):
	"""Every classical phylogenetic rate estimate for one clade.

	They are reported side by side deliberately. They make different
	assumptions - constant rates, complete sampling, no extinction, a
	Kingman genealogy - and where they disagree, the disagreement is the
	finding: a large positive coalescent rate next to a strongly negative
	gamma, for instance, says the clade grew and then stopped.
	"""
	nodes = tgm.subtree_nodes(tree, node)
	tips = [index for index in nodes if tree.is_tip[index] and np.isfinite(node_times[index])]
	row = {'node': int(node), 'node_label': label or tree.name[node],
		   'n_tips': len(tips), 'crown_age_years': float('nan'),
		   'stem_age_years': float('nan')}
	if len(tips) < 3:
		return row, []

	most_recent = float(np.max([node_times[index] for index in tips]))
	crown_time = node_times[node]
	row['first_date'] = float(np.min([node_times[index] for index in tips]))
	row['last_date'] = most_recent
	if np.isfinite(crown_time):
		row['crown_age_years'] = most_recent - float(crown_time)
	parent = tree.parent[node]
	if parent >= 0 and np.isfinite(node_times[parent]):
		row['stem_age_years'] = most_recent - float(node_times[parent])

	times, counts = tgm.lineages_through_time(tree, node_times, nodes)
	ltt = tgm.ltt_slope(times, counts)
	row['ltt_slope_per_year'] = ltt['slope_per_year']
	row['ltt_slope_se'] = ltt['std_error']
	row['ltt_r_squared'] = ltt['r_squared']
	row['ltt_p'] = ltt['p_value']

	for epsilon in EPSILON_GRID:
		tag = str(epsilon).replace('.', '')
		row['ms_crown_rate_eps%s' % tag] = tgm.magallon_sanderson(
			len(tips), row['crown_age_years'], epsilon, crown=True)
		row['ms_stem_rate_eps%s' % tag] = tgm.magallon_sanderson(
			len(tips), row['stem_age_years'], epsilon, crown=False)

	tree_length = float(np.nansum([branch_times[index] for index in nodes if index != node]))
	row['tree_length_years'] = tree_length
	row['yule_birth_rate_per_year'] = tgm.yule_birth_rate(len(tips), tree_length, crown=True)

	internal_ages = [most_recent - node_times[index] for index in nodes
					 if not tree.is_tip[index] and np.isfinite(node_times[index])
					 and len(tree.children[index]) > 1]
	row['gamma_statistic'] = tgm.gamma_statistic(internal_ages)

	intervals, coalescences, info = tgm.coalescent_intervals(tree, node_times, nodes)
	row['polytomy_fraction'] = info['polytomy_fraction']
	row['n_coalescences'] = info['n_coalescences']
	coalescent = tgm.coalescent_exponential_growth(intervals, coalescences)
	row['coalescent_growth_per_year'] = coalescent['growth_rate_per_year']
	row['coalescent_ci_low'] = coalescent['ci_low']
	row['coalescent_ci_high'] = coalescent['ci_high']
	row['coalescent_ne0'] = coalescent['n0']
	row['mean_terminal_branch_years'] = tgm.mean_terminal_branch_length(tree, branch_times, nodes)

	# The skyline is a trajectory, not a number, so it goes to its own table.
	# The single number worth carrying back is the slope of log(Ne) against
	# time, which is the growth rate the skyline implies - the same quantity the
	# exponential-growth coalescent estimates, but without assuming the rate was
	# constant, so the two disagreeing means the rate was not.
	points = tgm.skyline(intervals, coalescences, group=skyline_group)
	skyline_rows = [{'node': int(node), 'node_label': row['node_label'],
					 'tau_years': point['tau'], 'date': most_recent - point['tau'],
					 'ne': point['ne'], 'n_events': point.get('n_events', 1)}
					for point in points if point['ne'] > 0]
	if len(skyline_rows) >= 3:
		taus = np.array([entry['tau_years'] for entry in skyline_rows])
		sizes = np.log(np.array([entry['ne'] for entry in skyline_rows]))
		if np.ptp(taus) > 0:
			slope = float(np.polyfit(taus, sizes, 1)[0])
			row['skyline_growth_per_year'] = -slope
			row['skyline_points'] = len(skyline_rows)
	return row, skyline_rows


def analyse_tree(newick, cohort, allele_set, args):
	"""Build the working tree, then every tree-based estimate on it."""
	say('parsing the tree')
	full = tgm.parse_newick(newick)
	say('  %d nodes, %d tips' % (full.n_nodes, len(full.tips)))

	dates = dict(zip(cohort['primary_accession'], cohort['decimal_date']))
	known = set()
	for members in allele_set.evaluated.values():
		known |= members
	informative = set(dates) | known
	keep = [tip for tip in full.tips if full.name[tip] in informative]
	if len(keep) < 10:
		raise SystemExit('only %d tree tips are in the cohort; nothing to fit' % len(keep))
	say('  restricting to %d informative tips' % len(keep))
	tree, _ = tgm.induced_subtree(full, keep)
	del full

	tip_times = {tree.name[tip]: dates[tree.name[tip]] for tip in tree.tips
				 if tree.name[tip] in dates}
	node_times, dating = tgm.date_nodes(tree, tip_times, args.dating)
	if dating.get('fallback'):
		say('  dating fell back to mrca-bound: %s' % dating['fallback'])
	say('  dated with %s over %d tips' % (dating['method'], dating['n_dated_tips']))
	branch_times = tgm.time_branch_lengths(tree, node_times)

	tip_time = np.full(tree.n_nodes, np.nan)
	for tip in tree.tips:
		tip_time[tip] = dates.get(tree.name[tip], np.nan)

	span = float(np.nanmax(tip_time) - np.nanmin(tip_time)) if np.any(np.isfinite(tip_time)) else 1.0
	tau = args.lbi_tau if args.lbi_tau else max(span * DEFAULT_LBI_TAU_FRACTION, 1e-6)
	lbi_lengths = branch_times if args.lbi_units == 'time' else np.nan_to_num(tree.length, nan=0.0)
	_, lbi = tgm.local_branching_index(tree, lbi_lengths, tau)
	expansion = tgm.expansion_score(tree, tip_time, args.expansion_window)

	origin = args.bin_origin if args.bin_origin is not None else float(np.nanmin(tip_time))
	tip_bin = np.full(tree.n_nodes, -1, dtype=np.int64)
	finite = np.isfinite(tip_time)
	bin_index = np.floor((tip_time[finite] - origin) / args.bin_width).astype(int)
	tip_bin[np.flatnonzero(finite)] = bin_index
	n_bins = int(bin_index.max()) + 1 if bin_index.size else 0
	if n_bins < 3:
		raise SystemExit('the cohort spans fewer than three time bins; widen the date '
						 'range or shrink --bin-width')
	bin_mid = origin + (np.arange(n_bins) + 0.5) * args.bin_width

	say('  counting %d clades against %d time bins' % (tree.n_nodes, n_bins))
	bin_counts = tgm.clade_bin_counts(tree, tip_bin, n_bins)

	say('  scanning clades for frequency growth')
	scan = tgm.scan_clade_growth(tree, bin_counts, bin_mid,
								 min_clade_size=args.min_clade_size,
								 max_clade_fraction=args.max_clade_fraction,
								 max_nodes=args.max_nodes,
								 min_bin_total=args.min_bin_size)
	scan_frame = pd.DataFrame(scan)
	if not scan_frame.empty:
		scan_frame['q_value'] = gs.benjamini_hochberg(scan_frame['wald_p'])
		scan_frame['lbi'] = [lbi[int(node)] for node in scan_frame['node']]
		scan_frame['expansion_score'] = [expansion[int(node)] for node in scan_frame['node']]
		# One post-order pass over the whole tree, not one subtree walk per row:
		# with --max-nodes at its default this was the most expensive thing in
		# the scan, and quadratic on a ladder-shaped tree.
		earliest, latest = tgm.clade_time_bounds(tree, tip_time)
		scan_frame['first_date'] = [earliest[int(node)] for node in scan_frame['node']]
		scan_frame['last_date'] = [latest[int(node)] for node in scan_frame['node']]

	return {'tree': tree, 'node_times': node_times, 'branch_times': branch_times,
			'tip_time': tip_time, 'lbi': lbi, 'expansion': expansion,
			'bin_counts': bin_counts, 'bin_mid': bin_mid, 'dating': dating,
			'scan': scan_frame, 'lbi_tau': tau}


def analyse_origins(context, allele_set, cohort, allele_frame, args):
	"""Independent origins of each allele, and the growth of each one separately.

	This is the part a tree is genuinely needed for. A frequency rise measured
	across the whole database cannot tell a mutation that keeps arising and
	spreading from one that arose once inside a lineage which is spreading for
	its own reasons. Splitting the carriers by origin and fitting each
	separately is what separates the two.
	"""
	tree = context['tree']
	node_times = context['node_times']
	branch_times = context['branch_times']
	tip_time = context['tip_time']
	dates = dict(zip(cohort['primary_accession'], cohort['decimal_date']))

	ranked = allele_frame.copy()
	if not ranked.empty:
		ranked = ranked.reindex(ranked['logistic_growth_advantage_per_year'].abs()
								.sort_values(ascending=False).index)
	targets = list(ranked['allele'])[:args.max_tree_alleles] if not ranked.empty \
		else sorted(allele_set.carriers)[:args.max_tree_alleles]

	summaries, clusters, phylodynamics, skyline = [], [], [], []
	# Wild-type residues all resolve to the same ancestral lineage, so without a
	# cache the same 59,000-tip clade gets every phylodynamic rate fitted to it
	# once per allele.
	phylodynamic_cache = {}
	for allele in targets:
		site = allele_set.residues[allele]['site']
		carriers = allele_set.carriers[allele]
		evaluated = allele_set.evaluated[site]
		states = {}
		for tip in tree.tips:
			label = tree.name[tip]
			if label in carriers:
				states[label] = 1
			elif label in evaluated:
				states[label] = 0
		if sum(states.values()) < 2:
			continue
		result = tgm.reconstruct_origins(tree, states, root_state=args.ancestral_state,
										 resolution=args.origin_resolution)
		sizes = [len([t for t in cluster['tips'] if np.isfinite(tip_time[t])])
				 for cluster in result['clusters']]
		summary = {
			'allele': allele, 'protein': allele_set.protein, 'site': site,
			'residue': allele_set.residues[allele]['residue'],
			'n_carrier_tips': int(sum(states.values())),
			'ancestral_state': result['root_state'],
			'origins_min': result['origins_min'], 'origins_max': result['origins_max'],
			'origins_reconstructed': result['origins_realised'],
			'gains_min': result['gains_min'], 'gains_max': result['gains_max'],
			'losses_min': result['losses_min'], 'losses_max': result['losses_max'],
			'parsimony_cost': result['parsimony_cost'],
			'largest_origin_tips': max(sizes) if sizes else 0,
			'origins_with_min_tips': int(sum(1 for size in sizes if size >= args.min_cluster_size)),
			'carriers_in_largest_origin': (max(sizes) / sum(sizes)) if sum(sizes) else float('nan'),
		}
		# Near 1 means every carrier is its own origin: the residue keeps arising
		# and never goes anywhere, which is the signature of repeated selection
		# within hosts rather than a lineage transmitting. Near 0 means one
		# origin did all the spreading, and the frequency rise is that lineage's,
		# not necessarily the residue's.
		summary['origins_per_carrier'] = (float(summary['origins_min'])
										  / summary['n_carrier_tips']
										  if summary['n_carrier_tips'] else float('nan'))
		summaries.append(summary)

		ordered = sorted(result['clusters'],
						 key=lambda cluster: -len(cluster['tips']))[:args.max_clusters]
		for rank, cluster in enumerate(ordered, start=1):
			node = cluster['origin_node']
			dated = [tip for tip in cluster['tips'] if np.isfinite(tip_time[tip])]
			if len(dated) < args.min_cluster_size:
				continue
			member_times = np.array([tip_time[tip] for tip in dated])
			member_names = {tree.name[tip] for tip in dated}
			background = np.array([value for accession, value in dates.items()
								   if accession in evaluated and accession not in member_names])
			entry = {'allele': allele, 'origin_rank': rank, 'origin_node': int(node),
					 'node_label': tree.name[node], 'n_tips': len(dated),
					 'first_date': float(member_times.min()),
					 'last_date': float(member_times.max()),
					 'lbi': float(context['lbi'][node]),
					 'expansion_score': float(context['expansion'][node])}
			entry['date_span_years'] = float(member_times.max() - member_times.min())
			# A lineage sampled inside a single season has no trajectory. Fitting
			# one anyway returns a slope that only encodes *when* the lineage was
			# seen relative to the background, which is not a growth rate.
			if (len(background) >= args.min_allele_count
					and entry['date_span_years'] >= args.min_cluster_span):
				fitted, _ = treefree_metrics(member_times, background, args)
				for key in ('growth_advantage_per_year', 'ci_low', 'ci_high', 'lrt_p', 'method'):
					entry['logistic_' + key] = fitted['logistic_' + key]
			else:
				entry['logistic_skipped'] = ('span < %.2f years' % args.min_cluster_span
											 if entry['date_span_years'] < args.min_cluster_span
											 else 'too few background sequences')
			clusters.append(entry)
			if len(dated) >= args.min_phylodynamic_tips:
				if node not in phylodynamic_cache:
					fitted, points = clade_phylodynamics(
						tree, node, node_times, branch_times, label=tree.name[node],
						skyline_group=args.skyline_group)
					phylodynamic_cache[node] = fitted
					skyline.extend(points)
				record = dict(phylodynamic_cache[node])
				record.update({'allele': allele, 'origin_rank': rank,
							   'dating_method': context['dating']['method']})
				phylodynamics.append(record)

	return (pd.DataFrame(summaries), pd.DataFrame(clusters), pd.DataFrame(phylodynamics),
			pd.DataFrame(skyline))


# ---------------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------------

def _write(frame, directory, name):
	if frame is None or frame.empty:
		return None
	path = os.path.join(directory, name)
	frame.to_csv(path, sep='\t', index=False, float_format='%.6g')
	return path


def _format_table(frame, columns, limit=20):
	present = [column for column in columns if column in frame.columns]
	if frame.empty or not present:
		return '_nothing to report_\n'
	view = frame[present].head(limit).copy()
	for column in present:
		if pd.api.types.is_float_dtype(view[column]):
			view[column] = view[column].map(lambda value: '' if pd.isna(value) else '%.4g' % value)
	header = '| ' + ' | '.join(present) + ' |'
	rule = '|' + '|'.join(['---'] * len(present)) + '|'
	body = ['| ' + ' | '.join(str(value) for value in row) + ' |'
			for row in view.itertuples(index=False)]
	return '\n'.join([header, rule] + body) + '\n'


def write_report(path, args, cohort, allele_set, allele_frame, multinomial,
				 origins, clusters, phylodynamics, scan, context):
	span = (float(cohort['decimal_date'].min()), float(cohort['decimal_date'].max()))
	lines = ['# Lineage growth rates: %s' % allele_set.protein, '']
	lines.append('Database: `%s`  ' % os.path.abspath(args.db))
	lines.append('Generated: %s  ' % time.strftime('%Y-%m-%d %H:%M:%S'))
	lines.append('Genotype filter: %s (labels from %s)  '
				 % (args.genotype or 'none', cohort.attrs.get('genotype_source', 'none')))
	lines.append('Cohort: %d dated sequences, %.2f to %.2f  ' % (len(cohort), span[0], span[1]))
	lines.append('Allele source: `%s`  ' % allele_set.source)
	lines.append('Generation time assumed: %.3f years (sd %.3f)'
				 % (args.generation_time, args.generation_sd))
	lines.append('')
	if allele_set.notes.get('caveat'):
		lines += ['> **Denominator.** %s' % allele_set.notes['caveat'], '']
	extraction = allele_set.notes.get('extraction') or {}
	if extraction:
		bounds = allele_set.notes.get('cds_bounds')
		lines.append('Coordinates: master `%s`%s, CDS %s, read as **%s**; %d of %d requested '
					 'positions resolved.'
					 % (extraction.get('master'),
						(' (segment %s)' % extraction['segment']) if extraction.get('segment') else '',
						('%d..%d' % bounds) if bounds else 'from --protein-start',
						extraction.get('coord_space'),
						extraction.get('positions_resolved', 0),
						extraction.get('positions_requested', 0)))
		if extraction.get('rows_wrong_width'):
			lines.append('')
			lines.append('> **%d of %d alignment rows were skipped.** Their alignment is a '
						 'different width from the master\'s (%d), so they are in a '
						 'different coordinate space and reading fixed columns out of them '
						 'would invent residues. The frequencies below are over the %d rows '
						 'that share the master\'s frame.'
						 % (extraction['rows_wrong_width'], extraction['rows_read'],
							extraction['master_width'], extraction['rows_used']))
		lines.append('')
		lines.append('> **Site numbering** is the codon number counting from that CDS start '
					 'in the master\'s own frame. It is not a published numbering scheme. '
					 'For influenza HA in particular it is neither H1 nor H3 numbering: the '
					 'master frame is shared across subtypes that differ by indels, so the '
					 'numbers are comparable within this run and not to the literature.')
		lines.append('')
	if allele_set.notes.get('frame_warning'):
		lines += ['> **Reading frame check failed.** %s' % allele_set.notes['frame_warning'], '']
	if allele_set.notes.get('indel_note'):
		lines += ['> **Indel column.** %s' % allele_set.notes['indel_note'], '']
	if allele_set.notes.get('calibration'):
		calibration = allele_set.notes['calibration']
		lines += ['> **Protein start calibrated** to master coordinate %s, agreeing with '
				  '%.1f%% of %d catalogue calls.'
				  % (calibration['start'], 100 * calibration['agreement'],
					 calibration['n_scored']), '']
	if context:
		dating = context['dating']
		lines.append('Tree: %s, dated by `%s` over %d tips.'
					 % (args.tree or 'auto-selected', dating['method'], dating['n_dated_tips']))
		if dating['method'] == 'mrca-bound':
			lines.append('')
			lines.append('> **Node times are bounds, not estimates.** Each internal node is '
						 'dated to the earliest sample below it, so every clade looks younger '
						 'than it is and every age-based rate (Magallon-Sanderson, Yule, the '
						 'coalescent) is biased upwards. The frequency-based rates are '
						 'unaffected - they use sampling dates only.')
		if np.isfinite(dating.get('r_squared', float('nan'))):
			lines.append('Root-to-tip clock: rate %.3g, R2 %.3f.'
						 % (dating['clock_rate'], dating['r_squared']))
		lines.append('LBI scale tau = %.3f years.' % context['lbi_tau'])
		lines.append('')

	lines += ['## Alleles ranked by growth advantage', '',
			  'Slope of `logit(frequency) ~ time`: the growth-rate advantage of carriers '
			  'over everything else evaluated at that site, per year. `logistic_q` is the '
			  'likelihood-ratio p-value FDR corrected across alleles - the confidence '
			  'interval beside it is a Wald interval, so on a penalised fit the two can '
			  'disagree, and the likelihood ratio is the one to believe.', '']
	lines.append(_format_table(allele_frame, [
		'allele', 'n_carrier', 'n_evaluated', 'frequency',
		'logistic_growth_advantage_per_year', 'logistic_ci_low', 'logistic_ci_high',
		'logistic_q', 'logistic_doubling_time_years', 'binned_growth_advantage_per_year',
		'fit_p', 'mk_tau', 'late_vs_early_odds_ratio', 'logistic_method'], args.report_rows))

	if multinomial is not None and not multinomial.empty:
		lines += ['', '## Competing residues at each site (multinomial)', '',
				  'Jointly fitted, so the estimates at one site are mutually consistent. '
				  'The pivot residue is fixed at zero by construction.', '']
		lines.append(_format_table(multinomial.sort_values(
			'growth_advantage_per_year', ascending=False),
			['site', 'category', 'pivot', 'n', 'growth_advantage_per_year',
			 'ci_low', 'ci_high', 'p_value'], args.report_rows))

	if origins is not None and not origins.empty:
		lines += ['', '## Independent origins', '',
				  '`origins_min` and `origins_max` bracket the number of independent '
				  'acquisitions over *all* most-parsimonious histories. An allele whose '
				  'carriers sit in one origin is one observation however many sequences it '
				  'has; the same allele arising repeatedly is repeated evidence.', '']
		lines.append(_format_table(origins.sort_values('origins_min', ascending=False),
								   ['allele', 'n_carrier_tips', 'ancestral_state',
									'origins_min', 'origins_max', 'losses_min', 'losses_max',
									'origins_per_carrier', 'largest_origin_tips',
									'carriers_in_largest_origin', 'origins_with_min_tips'],
								   args.report_rows))
		lines.append('')
		lines.append('`ancestral_state = present` means the residue was reconstructed at the '
					 'root: it is the wild type, not a mutation that arose, and its '
					 'informative count is the losses. `origins_per_carrier` near 1 means '
					 'every carrier is its own origin - the residue keeps arising and does '
					 'not then spread - and near 0 means one origin did all the spreading, '
					 'so the frequency rise belongs to that lineage rather than to the '
					 'residue.')

	if clusters is not None and not clusters.empty:
		lines += ['', '## Growth of each independent origin', '']
		lines.append(_format_table(clusters.sort_values('n_tips', ascending=False),
								   ['allele', 'origin_rank', 'n_tips', 'first_date',
									'last_date', 'logistic_growth_advantage_per_year',
									'logistic_ci_low', 'logistic_ci_high', 'logistic_lrt_p',
									'lbi', 'expansion_score'], args.report_rows))

	if phylodynamics is not None and not phylodynamics.empty:
		lines += ['', '## Phylodynamic rates per origin lineage', '',
				  'Independent of the frequency fits above: these read a rate off the shape '
				  'of the genealogy rather than off how common the lineage is.', '']
		if context and context['dating']['method'] == 'mrca-bound':
			lines += ['> **These four columns are not trustworthy under `mrca-bound` '
					  'dating.** Dating an internal node to the earliest sample below it '
					  'puts every deep node at the start of the window, so the '
					  'lineages-through-time curve starts at its maximum and only falls, '
					  'and the LTT slope, the Magallon-Sanderson and Yule rates and gamma '
					  'all describe that instead of the virus - a gamma of -50 is the '
					  'dating, not a real slowdown. Use a tree with informative branch '
					  'lengths and `--dating root-to-tip`, or read only the '
					  'frequency-based rates.', '']
		deduped = phylodynamics.drop_duplicates(subset=['node'])
		lines.append(_format_table(deduped.sort_values('n_tips', ascending=False),
								   ['allele', 'node_label', 'n_tips', 'crown_age_years',
									'ltt_slope_per_year', 'ms_crown_rate_eps00',
									'ms_crown_rate_eps09', 'yule_birth_rate_per_year',
									'gamma_statistic', 'coalescent_growth_per_year',
									'coalescent_ci_low', 'coalescent_ci_high',
									'skyline_growth_per_year', 'polytomy_fraction'],
								   args.report_rows))

	if scan is not None and not scan.empty:
		lines += ['', '## Fastest-growing clades in the tree', '',
				  'Every clade of at least %d cohort sequences, fitted against the whole '
				  'cohort and FDR corrected. Independent of the protein - this is where a '
				  'growing lineage shows up whether or not it carries a catalogued '
				  'mutation.' % args.min_clade_size, '']
		lines.append(_format_table(scan[scan['q_value'] < args.report_q]
								   if 'q_value' in scan.columns else scan,
								   ['node_label', 'clade_size', 'clade_fraction', 'first_date',
									'n_time_bins_occupied', 'growth_advantage_per_year',
									'ci_low', 'ci_high', 'q_value', 'lbi', 'expansion_score',
									'dispersion'],
								   args.report_rows))

	lines += ['', '## Reading these numbers', '',
			  '* A growth advantage is *relative*. It says the lineage is outcompeting the '
			  'sequences it is being compared against, not that it is growing in absolute '
			  'terms - compare `count_growth_per_year` for that.',
			  '* Sampling is not random. Frequencies are frequencies **among sequenced '
			  'isolates**, and anything that changes who gets sequenced (a treatment '
			  'cohort, an outbreak investigation, a new surveillance programme) moves them '
			  'without any change in the virus.',
			  '* `logistic_method = firth` marks a fit that was penalised because the '
			  'allele is nearly perfectly separated in time. The estimate is finite and '
			  'usable; the plain maximum-likelihood one would not have been.',
			  '* Overdispersion (`dispersion` far above 1) means the sequences are '
			  'correlated within time bins. The binned intervals already account for it; '
			  'the individual-level ones do not.',
			  '* A clade first sampled at the very end of the window will always come out '
			  'near the top of the clade scan with a tiny q-value. The model is being '
			  'asked about the whole window, and over that window the clade did go from '
			  'absent to present. Check `first_date` and `n_time_bins_occupied`, and use '
			  '`--fit-window` to ask the question over a period the cohort actually '
			  'covers.',
			  '']
	with open(path, 'w', encoding='utf-8') as handle:
		handle.write('\n'.join(lines))
	return path


def write_plots(directory, allele_frame, trajectories, context, args):
	"""Frequency trajectories, and an LTT plot for the biggest clade."""
	import matplotlib
	matplotlib.use('Agg')
	import matplotlib.pyplot as plt

	written = []
	if trajectories is not None and not trajectories.empty and not allele_frame.empty:
		top = list(allele_frame['allele'].head(6))
		figure, axis = plt.subplots(figsize=(8, 5))
		for allele in top:
			series = trajectories[trajectories['allele'] == allele].sort_values('bin_mid')
			series = series[series['n'] >= args.min_bin_size]
			if series.empty:
				continue
			axis.plot(series['bin_mid'], series['frequency'], marker='o', label=allele)
			axis.fill_between(series['bin_mid'], series['ci_low'], series['ci_high'], alpha=0.15)
		axis.set_xlabel('year')
		axis.set_ylabel('frequency among evaluated sequences')
		axis.set_title('Allele frequency trajectories')
		axis.legend(fontsize='small')
		figure.tight_layout()
		path = os.path.join(directory, 'frequency_trajectories.png')
		figure.savefig(path, dpi=150)
		plt.close(figure)
		written.append(path)

	if context and context['scan'] is not None and not context['scan'].empty:
		tree = context['tree']
		node = int(context['scan'].iloc[0]['node'])
		times, counts = tgm.lineages_through_time(tree, context['node_times'],
												  tgm.subtree_nodes(tree, node))
		if len(times) > 3:
			figure, axis = plt.subplots(figsize=(7, 4.5))
			axis.step(times, counts, where='post')
			axis.set_yscale('log')
			axis.set_xlabel('year')
			axis.set_ylabel('lineages')
			axis.set_title('Lineages through time: %s' % (tree.name[node] or node))
			figure.tight_layout()
			path = os.path.join(directory, 'ltt_fastest_clade.png')
			figure.savefig(path, dpi=150)
			plt.close(figure)
			written.append(path)
	return written


def write_back(db_path, run_id, args, frames):
	"""Append the results into the database under a run id."""
	conn = sqlite3.connect(db_path)
	try:
		summary = pd.DataFrame([{
			'run_id': run_id, 'created_at': time.strftime('%Y-%m-%d %H:%M:%S'),
			'protein': args.protein, 'genotype': args.genotype or '',
			'tree': args.tree or 'auto', 'dating': args.dating,
			'parameters': json.dumps({key: value for key, value in vars(args).items()
									  if key not in {'db'}}, default=str),
		}])
		summary.to_sql('lineage_growth_runs', conn, if_exists='append', index=False)
		for table, frame in frames.items():
			if frame is None or frame.empty:
				continue
			stamped = frame.copy()
			stamped.insert(0, 'run_id', run_id)
			stamped.to_sql(table, conn, if_exists='append', index=False)
		conn.commit()
	finally:
		conn.close()


# ---------------------------------------------------------------------------
# Command line
# ---------------------------------------------------------------------------

def parse_args(argv=None):
	parser = argparse.ArgumentParser(
		description='Estimate lineage growth rates for a protein in a V-gTK database.',
		formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument('--db', required=True, help='finished V-gTK sqlite database')
	parser.add_argument('--protein', help='protein to analyse (see --list-proteins)')
	parser.add_argument('--out', help='output directory')
	parser.add_argument('--list-proteins', action='store_true',
						help='print the proteins this database can be asked about and exit')

	cohort = parser.add_argument_group('cohort')
	cohort.add_argument('--genotype', help='comma-separated genotypes to keep')
	cohort.add_argument('--genotype-level', choices=('genotype', 'subtype'), default='genotype')
	cohort.add_argument('--genotype-column', help='meta_data column holding the genotype')
	cohort.add_argument('--genotype-reflist',
						help='curated reference list; genotypes are inherited from tree '
							 'neighbours, for a database with no genotype columns')
	cohort.add_argument('--min-year', type=int, help='earliest collection year to keep')
	cohort.add_argument('--max-year', type=int, help='latest collection year to keep')
	cohort.add_argument('--min-precision', choices=('year', 'month', 'day'), default='year',
						help='discard dates less precise than this')
	cohort.add_argument('--segment', help='restrict to one segment')
	cohort.add_argument('--include-excluded', action='store_true',
						help='keep sequences the database flags with exclusion_status. '
							 'They are dropped by default, which is right when the flag '
							 'means failed QC - but a build can carry exclusions that are '
							 'not about sequence quality at all, and then the default '
							 'silently removes most of the data. Check what the flag means '
							 'in your build before using this')
	cohort.add_argument('--fit-window', type=float,
						help='keep only the last N years of the cohort. Worth setting: a '
							 'database spanning a century puts most of its time axis in '
							 'decades holding a handful of sequences, and a clade first '
							 'seen at the very end of such a window always looks like the '
							 'fastest-growing thing in the tree')

	alleles = parser.add_argument_group('alleles')
	alleles.add_argument('--allele-source', choices=('auto', 'catalog', 'features', 'alignment'),
						 default='auto')
	alleles.add_argument('--protein-start', type=int,
						 help='first nucleotide of the protein in master coordinates '
							  '(--allele-source alignment)')
	alleles.add_argument('--calibrate-start', action='store_true',
						 help='infer --protein-start from the catalogue calls, and refuse '
							  'if the agreement is poor')
	alleles.add_argument('--positions', help='comma-separated residue positions to genotype')
	alleles.add_argument('--coord-space', choices=pal.COORD_SPACES, default='auto',
						 help="how to read features.cds_start against the stored alignment. "
							  "'column' treats it as an alignment column, 'ungapped' as a "
							  "position in the master's own sequence. They differ only when "
							  "the master carries gaps - the influenza HA master has 132, "
							  "and reading it the wrong way puts the frame 18 codons out. "
							  "'auto' translates the master both ways and takes the reading "
							  "with no internal stop codons")
	alleles.add_argument('--min-allele-count', type=int, default=10,
						 help='skip alleles with fewer carriers, or fewer background '
							  'sequences, than this')

	timing = parser.add_argument_group('time')
	timing.add_argument('--bin-width', type=float, default=1.0, help='time bin width in years')
	timing.add_argument('--bin-origin', type=float, help='left edge of the first bin')
	timing.add_argument('--min-bin-size', type=int, default=5,
						help='ignore time bins holding fewer sequences than this')
	timing.add_argument('--generation-time', type=float, default=1.0,
						help='mean generation interval in years, for the per-generation '
							 'selection coefficient and the reproduction number')
	timing.add_argument('--generation-sd', type=float, default=0.0,
						help='standard deviation of the generation interval; 0 assumes a '
							 'fixed interval, which overstates R')

	tree = parser.add_argument_group('tree')
	tree.add_argument('--no-tree', action='store_true', help='frequency estimators only')
	tree.add_argument('--tree', help='name of the row in the trees table to use')
	tree.add_argument('--tree-file',
						help='a newick file to use instead of a tree stored in the '
							 'database. This is what to reach for when a tree has been '
							 'rebuilt but not yet written back: the growth rates then '
							 'describe the new topology rather than the stale one the '
							 'database still holds')
	tree.add_argument('--dating', choices=tgm.DATING_METHODS, default='mrca-bound')
	tree.add_argument('--lbi-tau', type=float, help='local branching index scale in years')
	tree.add_argument('--lbi-units', choices=('time', 'divergence'), default='time')
	tree.add_argument('--expansion-window', type=float, default=5.0,
						help='window in years for the recent-versus-previous expansion score')
	tree.add_argument('--min-clade-size', type=int, default=20)
	tree.add_argument('--max-clade-fraction', type=float, default=0.9)
	tree.add_argument('--max-nodes', type=int, default=20000,
						help='largest number of clades to fit in the scan')
	tree.add_argument('--max-tree-alleles', type=int, default=25,
						help='alleles to reconstruct origins for, strongest signal first')
	tree.add_argument('--max-clusters', type=int, default=20,
						help='independent origins per allele to report in detail')
	tree.add_argument('--min-cluster-size', type=int, default=3)
	tree.add_argument('--min-cluster-span', type=float, default=1.0,
						help='shortest sampling span, in years, an origin lineage must '
							 'cover before a growth rate is fitted to it')
	tree.add_argument('--ancestral-state', choices=('auto', 'absent', 'present'),
						default='auto',
						help="state at the root. 'auto' lets parsimony choose, which is "
							 "what keeps a wild-type residue from being counted as an "
							 "independent origin in every sequence that carries it")
	tree.add_argument('--min-phylodynamic-tips', type=int, default=10,
						help='smallest origin lineage to fit phylodynamic rates to')
	tree.add_argument('--skyline-group', type=int, default=10,
						help='coalescent events averaged per skyline point. The classic '
							 'ungrouped estimator has enormous variance per point and is '
							 'mostly noise')
	tree.add_argument('--origin-resolution', choices=('deltran', 'acctran'), default='deltran',
						help='which most-parsimonious history to cut into clusters; the '
							 'reported origin range covers both')

	output = parser.add_argument_group('output')
	output.add_argument('--report-rows', type=int, default=20)
	output.add_argument('--report-q', type=float, default=0.05)
	output.add_argument('--plots', action='store_true')
	output.add_argument('--write_db', action='store_true',
						help='append the results into the database (it is opened read-only '
							 'for everything else)')
	return parser.parse_args(argv)


def main(argv=None):
	args = parse_args(argv)
	conn = open_database(args.db)

	if args.list_proteins:
		frame = pal.list_proteins(conn)
		print(frame.to_string(index=False) if not frame.empty
			  else 'this database records no proteins')
		return 0
	if not args.protein or not args.out:
		raise SystemExit('--protein and --out are required unless --list-proteins is given')
	os.makedirs(args.out, exist_ok=True)

	newick = None
	if not args.no_tree:
		if args.tree_file:
			if not os.path.isfile(args.tree_file):
				raise SystemExit('tree file not found: %s' % args.tree_file)
			with open(args.tree_file, 'r', encoding='utf-8') as handle:
				newick = handle.read().strip()
			if not newick:
				raise SystemExit('tree file is empty: %s' % args.tree_file)
			args.tree = os.path.basename(args.tree_file)
			say('using tree from file %s (%.1f MB)'
				% (args.tree_file, os.path.getsize(args.tree_file) / 1e6))
		else:
			name, source, newick = choose_tree(conn, args.tree, segment=args.segment)
			args.tree = name
			say('using tree %r (source %s)' % (name, source))

	say('loading the cohort')
	cohort = pal.load_cohort(conn, genotype_column=args.genotype_column,
							 min_precision=args.min_precision, min_year=args.min_year,
							 max_year=args.max_year, segment=args.segment,
							 drop_excluded=not args.include_excluded)
	if args.include_excluded:
		say('keeping sequences flagged exclusion_status - check what that flag means here')
	cohort = apply_genotypes(conn, cohort, args, newick)
	if args.fit_window and not cohort.empty:
		cutoff = float(cohort['decimal_date'].max()) - args.fit_window
		before = len(cohort)
		cohort = cohort[cohort['decimal_date'] >= cutoff].reset_index(drop=True)
		say('fit window: kept %d of %d sequences from %.2f onwards'
			% (len(cohort), before, cutoff))
	if cohort.empty:
		raise SystemExit('the cohort is empty after filtering')
	say('cohort: %d sequences, %.2f to %.2f'
		% (len(cohort), cohort['decimal_date'].min(), cohort['decimal_date'].max()))
	if args.bin_origin is None:
		args.bin_origin = float(cohort['decimal_date'].min())

	say('resolving alleles for %s' % args.protein)
	allele_set = build_alleles(conn, args, cohort)
	say('  %d alleles over %d sites' % (len(allele_set.carriers), len(allele_set.evaluated)))

	say('fitting frequency growth models')
	allele_frame, trajectories, multinomial = analyse_alleles(allele_set, cohort, args)
	say('  %d alleles passed the count thresholds' % len(allele_frame))

	context = None
	origins = clusters = phylodynamics = skyline = scan = pd.DataFrame()
	if not args.no_tree:
		context = analyse_tree(newick, cohort, allele_set, args)
		scan = context['scan']
		say('reconstructing independent origins')
		origins, clusters, phylodynamics, skyline = analyse_origins(
			context, allele_set, cohort, allele_frame, args)
		say('  %d alleles reconstructed, %d origin lineages detailed'
			% (len(origins), len(clusters)))

	frames = {'lineage_growth_alleles': allele_frame,
			  'lineage_growth_trajectories': trajectories,
			  'lineage_growth_multinomial': multinomial,
			  'lineage_growth_origins': origins,
			  'lineage_growth_origin_lineages': clusters,
			  'lineage_growth_phylodynamics': phylodynamics,
			  'lineage_growth_skyline': skyline,
			  'lineage_growth_clades': scan}
	for table, frame in frames.items():
		_write(frame, args.out, table.replace('lineage_growth_', '') + '.tsv')

	report = write_report(os.path.join(args.out, 'report.md'), args, cohort, allele_set,
						  allele_frame, multinomial, origins, clusters, phylodynamics,
						  scan, context)
	say('report written to %s' % report)

	if args.plots:
		for path in write_plots(args.out, allele_frame, trajectories, context, args):
			say('plot written to %s' % path)

	if args.write_db:
		run_id = '%s_%s_%s' % (args.protein, (args.genotype or 'all').replace(',', '-'),
							   time.strftime('%Y%m%dT%H%M%S'))
		conn.close()
		write_back(args.db, run_id, args, frames)
		say('results appended to %s under run_id %s' % (args.db, run_id))
	else:
		conn.close()
	return 0


if __name__ == '__main__':
	sys.exit(main())
