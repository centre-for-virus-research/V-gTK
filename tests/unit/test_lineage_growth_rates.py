import os
import sqlite3

import pandas as pd
import pytest

import lineage_growth_rates as lgr
import protein_alleles as pal


# A four-codon "genome": ATG <site 2> CCC GGG. Site 2 is AAA (K) in carriers and
# GAA (E) in everything else, so the frame contains no stop codon and the
# reading-frame check should pass.
CARRIER_CODON = 'AAA'
BACKGROUND_CODON = 'GAA'


def _sequences():
	# Ten sequences a year from 2000 to 2019; carriers start in 2010 and rise
	# from one a year to ten, which is a clean logistic sweep from 0 to 1.
	rows = []
	for index, year in enumerate(range(2000, 2020)):
		carriers = max(0, year - 2009)
		for slot in range(10):
			rows.append({
				'accession': 'Q%04d' % len(rows),
				'year': year,
				'month': 1 + (slot % 12),
				'carrier': slot < carriers,
			})
	return rows


def _alignment(carrier):
	return 'ATG' + (CARRIER_CODON if carrier else BACKGROUND_CODON) + 'CCC' + 'GGG'


@pytest.fixture
def database(tmp_path):
	path = tmp_path / 'synthetic.db'
	conn = sqlite3.connect(str(path))
	rows = _sequences()

	meta = [{'primary_accession': 'MASTER', 'collection_year': '2000',
			 'collection_mon': '1', 'collection_day': '1', 'accession_type': 'master',
			 'segment': '', 'exclusion_status': '', 'host': '', 'country': ''}]
	alignments = [{'primary_accession': 'MASTER', 'alignment': _alignment(False),
				   'alignment_name': 'MASTER', 'segment': '', 'sequence_id': 'MASTER'}]
	features = [{'accession': 'MASTER', 'master_ref_accession': 'MASTER',
				 'reference_accession': 'MASTER', 'aln_start': '1', 'aln_end': '12',
				 'cds_start': '1', 'cds_end': '12', 'product': 'G'}]
	mutations = []
	for row in rows:
		meta.append({'primary_accession': row['accession'],
					 'collection_year': str(row['year']),
					 'collection_mon': str(row['month']), 'collection_day': '15',
					 'accession_type': 'query', 'segment': '', 'exclusion_status': '',
					 'host': '', 'country': ''})
		alignments.append({'primary_accession': row['accession'],
						   'alignment': _alignment(row['carrier']),
						   'alignment_name': 'MASTER', 'segment': '',
						   'sequence_id': row['accession']})
		features.append({'accession': row['accession'], 'master_ref_accession': 'MASTER',
						 'reference_accession': 'MASTER', 'aln_start': '1', 'aln_end': '12',
						 'cds_start': '1', 'cds_end': '12', 'product': 'G'})
		residue = 'K' if row['carrier'] else 'E'
		mutations.append({'primary_accession': row['accession'], 'mutation_id': 'G:2%s' % residue,
						  'protein_name': 'G', 'segment': '', 'aa_position': 2,
						  'alt_residue': residue, 'combination_id': None})
		# A duplicate row under another drug combination, which must not double
		# the sequence's weight in any frequency.
		mutations.append({'primary_accession': row['accession'], 'mutation_id': 'G:2%s' % residue,
						  'protein_name': 'G', 'segment': '', 'aa_position': 2,
						  'alt_residue': residue, 'combination_id': 'G:2%s+9Z' % residue})

	carriers = [row['accession'] for row in rows if row['carrier']]
	others = [row['accession'] for row in rows if not row['carrier']]
	newick = '((%s)fast,(%s)slow,MASTER)root;' % (','.join(carriers), ','.join(others))

	pd.DataFrame(meta).to_sql('meta_data', conn, index=False)
	pd.DataFrame(alignments).to_sql('sequence_alignment', conn, index=False)
	pd.DataFrame(features).to_sql('features', conn, index=False)
	pd.DataFrame(mutations).to_sql('sequence_mutations', conn, index=False)
	pd.DataFrame([{'name': 'usher', 'source': 'usher', 'segment_key': None,
				   'segment': None, 'newick': newick, 'created_at': '2020-01-01'}]
				 ).to_sql('trees', conn, index=False)
	pd.DataFrame([{'primary_accession': 'Q0000', 'reason': 'test exclusion'}]
				 ).to_sql('excluded_accessions', conn, index=False)
	conn.commit()
	conn.close()
	return str(path)


def _read(directory, name):
	return pd.read_csv(os.path.join(directory, name), sep='\t')


# --- dates ---------------------------------------------------------------

def test_partial_dates_land_mid_period_not_at_the_start():
	# Defaulting to 1 January would push every imprecise date earlier than the
	# precise ones, and imprecision is correlated with time.
	assert pal.decimal_date(2020, None, None) == (2020.5, 'year')
	month_value, precision = pal.decimal_date(2020, 7, None)
	assert precision == 'month' and 2020.4 < month_value < 2020.6
	day_value, precision = pal.decimal_date(2020, 3, 15)
	assert precision == 'day' and abs(day_value - (2020 + 74.5 / 366)) < 1e-9
	assert pal.decimal_date('', None, None)[1] == 'none'


def test_cohort_drops_excluded_accessions(database):
	conn = sqlite3.connect(database)
	cohort = pal.load_cohort(conn)
	assert 'Q0000' not in set(cohort['primary_accession'])
	assert len(cohort) == 200          # 200 queries + master - 1 excluded
	conn.close()


# --- end to end ----------------------------------------------------------

def test_tree_free_run_recovers_the_planted_growth(database, tmp_path):
	out = tmp_path / 'notree'
	lgr.main(['--db', database, '--protein', 'G', '--out', str(out),
			  '--allele-source', 'features', '--no-tree'])
	alleles = _read(str(out), 'alleles.tsv').set_index('allele')
	assert alleles.loc['G:2K', 'logistic_growth_advantage_per_year'] > 0.5
	assert alleles.loc['G:2K', 'logistic_q'] < 1e-6
	# The competing residue must fall at the same rate it rises.
	assert alleles.loc['G:2E', 'logistic_growth_advantage_per_year'] == pytest.approx(
		-alleles.loc['G:2K', 'logistic_growth_advantage_per_year'], rel=1e-6)
	assert alleles.loc['G:2K', 'fit_p'] < 0.05
	# Ten leading bins at frequency zero are ties, so tau cannot reach 1 even
	# though the trajectory never falls.
	assert alleles.loc['G:2K', 'mk_tau'] > 0.7
	assert alleles.loc['G:2K', 'late_vs_early_odds_ratio'] > 1
	assert os.path.exists(os.path.join(str(out), 'report.md'))


def test_catalog_source_deduplicates_combination_rows(database, tmp_path):
	# Every sequence appears twice in sequence_mutations; the carrier count must
	# still be the number of sequences.
	out = tmp_path / 'catalog'
	lgr.main(['--db', database, '--protein', 'G', '--out', str(out),
			  '--allele-source', 'catalog', '--no-tree'])
	alleles = _read(str(out), 'alleles.tsv').set_index('allele')
	# The master carries no catalogue call, so the denominator is the 199
	# queries that were evaluated - not the whole cohort.
	assert alleles.loc['G:2K', 'n_evaluated'] == 199
	assert alleles.loc['G:2K', 'n_carrier'] == 55


def test_multinomial_agrees_with_the_pairwise_fit_for_two_residues(database, tmp_path):
	out = tmp_path / 'mlr'
	lgr.main(['--db', database, '--protein', 'G', '--out', str(out),
			  '--allele-source', 'features', '--no-tree'])
	alleles = _read(str(out), 'alleles.tsv').set_index('allele')
	multinomial = _read(str(out), 'multinomial.tsv').set_index('allele')
	assert multinomial.loc['G:2E', 'growth_advantage_per_year'] == 0.0   # the pivot
	assert multinomial.loc['G:2K', 'growth_advantage_per_year'] == pytest.approx(
		alleles.loc['G:2K', 'logistic_growth_advantage_per_year'], rel=0.05)


def test_tree_run_finds_one_origin_and_the_growing_clade(database, tmp_path):
	out = tmp_path / 'tree'
	lgr.main(['--db', database, '--protein', 'G', '--out', str(out),
			  '--allele-source', 'features'])
	origins = _read(str(out), 'origins.tsv').set_index('allele')
	assert origins.loc['G:2K', 'ancestral_state'] == 'absent'
	assert origins.loc['G:2K', 'origins_min'] == 1
	assert origins.loc['G:2K', 'origins_max'] == 1
	# The background residue is the ancestral one, so its count is its losses.
	assert origins.loc['G:2E', 'ancestral_state'] == 'present'
	assert origins.loc['G:2E', 'losses_min'] >= 1

	clades = _read(str(out), 'clades.tsv')
	fastest = clades.sort_values('growth_advantage_per_year', ascending=False).iloc[0]
	assert fastest['node_label'] == 'fast'
	assert fastest['growth_advantage_per_year'] > 0.5
	assert fastest['q_value'] < 0.01


def test_fit_window_restricts_the_cohort(database, tmp_path):
	out = tmp_path / 'window'
	lgr.main(['--db', database, '--protein', 'G', '--out', str(out),
			  '--allele-source', 'features', '--no-tree', '--fit-window', '5'])
	alleles = _read(str(out), 'alleles.tsv').set_index('allele')
	assert alleles.loc['G:2K', 'n_evaluated'] < 100


def test_results_can_be_written_back_under_a_run_id(database, tmp_path):
	out = tmp_path / 'writeback'
	lgr.main(['--db', database, '--protein', 'G', '--out', str(out),
			  '--allele-source', 'features', '--no-tree', '--write_db'])
	conn = sqlite3.connect(database)
	runs = pd.read_sql_query('SELECT * FROM lineage_growth_runs', conn)
	alleles = pd.read_sql_query('SELECT * FROM lineage_growth_alleles', conn)
	conn.close()
	assert len(runs) == 1 and runs.iloc[0]['protein'] == 'G'
	assert set(alleles['run_id']) == {runs.iloc[0]['run_id']}
	assert 'G:2K' in set(alleles['allele'])


def test_a_wrong_reading_frame_is_reported_rather_than_analysed_silently(database, tmp_path):
	out = tmp_path / 'badframe'
	# Starting one base late puts the whole region out of frame.
	lgr.main(['--db', database, '--protein', 'G', '--out', str(out),
			  '--allele-source', 'alignment', '--protein-start', '2',
			  '--positions', '1,2,3', '--no-tree', '--min-allele-count', '5'])
	report = open(os.path.join(str(out), 'report.md'), encoding='utf-8').read()
	assert 'Reading frame check failed' in report


def test_a_missing_protein_names_what_is_available(database, tmp_path):
	with pytest.raises(SystemExit) as failure:
		lgr.main(['--db', database, '--protein', 'NOPE', '--out', str(tmp_path / 'x'),
				  '--no-tree'])
	assert 'NOPE' in str(failure.value) and 'G' in str(failure.value)


def test_requesting_a_genotype_without_one_in_the_database_explains_the_options(database, tmp_path):
	with pytest.raises(SystemExit) as failure:
		lgr.main(['--db', database, '--protein', 'G', '--out', str(tmp_path / 'y'),
				  '--no-tree', '--genotype', '1'])
	assert '--genotype-reflist' in str(failure.value)


# --- coordinate conventions and segmented builds -------------------------

def test_master_is_selected_per_segment(database):
	# The synthetic build has one master; the influenza builds have eight, one
	# per segment, and picking the wrong one is silent.
	conn = sqlite3.connect(database)
	try:
		assert pal.count_masters(conn) == 1
		assert pal.master_accession(conn) == 'MASTER'
		assert pal.master_accession(conn, segment='') == 'MASTER'
		assert pal.master_accession(conn, segment='4') is None
	finally:
		conn.close()


def test_the_two_coordinate_readings_agree_only_on_a_gapless_master():
	gapless = 'ATG' * 10
	assert (pal.codon_columns(gapless, 1, [1, 2, 3], 'column')
			== pal.codon_columns(gapless, 1, [1, 2, 3], 'ungapped'))
	# Two leading gaps shift the ungapped reading off the column reading.
	gapped = '--' + 'ATG' * 10
	assert (pal.codon_columns(gapped, 3, [1], 'column')
			!= pal.codon_columns(gapped, 3, [1], 'ungapped'))


def test_coordinate_reading_is_detected_from_the_master_translation():
	# A master with leading gaps: read by column the CDS is a clean protein,
	# read as ungapped it runs off the end of the sequence.
	coding = 'ATG' + 'AAA' * 20 + 'TAA'
	master = '-' * 12 + coding
	start = 13
	end = start + len(coding) - 1
	space, scores = pal.detect_coord_space(master, start, end)
	assert space == 'column'
	assert scores['column']['internal_stops'] == 0
	assert scores['column']['missing'] == 0
	assert scores['ungapped']['missing'] > 0


def test_multiple_masters_require_an_explicit_segment(database, tmp_path):
	import shutil
	copy = str(tmp_path / 'two_masters.db')
	shutil.copy(database, copy)
	conn = sqlite3.connect(copy)
	conn.execute("INSERT INTO meta_data (primary_accession, collection_year, collection_mon, "
				 "collection_day, accession_type, segment) "
				 "VALUES ('MASTER2','2001','1','1','master','7')")
	conn.commit()
	conn.close()
	with pytest.raises(SystemExit) as failure:
		lgr.main(['--db', copy, '--protein', 'G', '--out', str(tmp_path / 'z'),
				  '--allele-source', 'features', '--no-tree'])
	assert '--segment' in str(failure.value)


def test_tree_is_chosen_for_the_requested_segment(database, tmp_path):
	# A segmented build holds one tree per segment. Choosing by "longest newick"
	# hands a segment-6 cohort the segment-4 tree, whose tips it shares none of.
	import shutil
	copy = str(tmp_path / 'two_trees.db')
	shutil.copy(database, copy)
	conn = sqlite3.connect(copy)
	try:
		conn.execute("UPDATE trees SET segment = '4'")
		conn.execute("INSERT INTO trees (name, source, segment_key, segment, newick, "
					 "created_at) VALUES ('usher_seg6','usher','s6','6','(X,Y);',"
					 "'2020-01-01')")
		conn.commit()
		assert lgr.choose_tree(conn, None, segment='6')[0] == 'usher_seg6'
		# The bigger segment-4 tree still wins when segment 4 is what was asked for.
		assert lgr.choose_tree(conn, None, segment='4')[0] == 'usher'
		# An explicit --tree name still overrides the segment.
		assert lgr.choose_tree(conn, 'usher_seg6', segment='4')[0] == 'usher_seg6'
	finally:
		conn.close()


def test_a_spliced_product_is_refused_rather_than_read_as_one_span(database, tmp_path):
	# Influenza M2 and NEP are annotated as two feature rows, one per exon.
	# Taking the first row silently translates 26 nucleotides of exon 1.
	import shutil
	copy = str(tmp_path / 'spliced.db')
	shutil.copy(database, copy)
	conn = sqlite3.connect(copy)
	try:
		conn.execute("INSERT INTO features (accession, master_ref_accession, "
					 "reference_accession, aln_start, aln_end, cds_start, cds_end, product) "
					 "VALUES ('MASTER','MASTER','MASTER','1','12','7','12','G')")
		conn.commit()
	finally:
		conn.close()
	with pytest.raises(SystemExit) as failure:
		lgr.main(['--db', copy, '--protein', 'G', '--out', str(tmp_path / 'spl'),
				  '--allele-source', 'features', '--no-tree'])
	message = str(failure.value)
	assert 'spliced' in message and 'two' not in message.lower().split()[0:1]
	assert '--protein-start' in message


def test_a_tree_file_can_stand_in_for_the_database_tree(database, tmp_path):
	# A rebuilt tree is usually not in the database yet. --tree-file lets the
	# growth rates describe the new topology instead of the stale stored one.
	conn = sqlite3.connect(database)
	try:
		newick = conn.execute("SELECT newick FROM trees WHERE name='usher'").fetchone()[0]
	finally:
		conn.close()
	path = tmp_path / 'external.nwk'
	path.write_text(newick, encoding='utf-8')

	out = tmp_path / 'fromfile'
	lgr.main(['--db', database, '--protein', 'G', '--out', str(out),
			  '--allele-source', 'features', '--tree-file', str(path)])
	assert os.path.exists(os.path.join(str(out), 'origins.tsv'))
	report = open(os.path.join(str(out), 'report.md'), encoding='utf-8').read()
	assert 'external.nwk' in report


def test_a_missing_tree_file_is_refused_rather_than_ignored(database, tmp_path):
	# Silently falling back to the database tree would answer a different
	# question from the one that was asked.
	with pytest.raises(SystemExit) as failure:
		lgr.main(['--db', database, '--protein', 'G', '--out', str(tmp_path / 'x'),
				  '--allele-source', 'features',
				  '--tree-file', str(tmp_path / 'nope.nwk')])
	assert 'tree file not found' in str(failure.value)


def test_excluded_sequences_are_dropped_by_default_and_kept_on_request(database, tmp_path):
	# A build can carry exclusion flags that are not about sequence quality. The
	# default is still to honour them, but silently losing most of the data is
	# worse than being able to ask for it back.
	import shutil
	copy = str(tmp_path / 'flagged.db')
	shutil.copy(database, copy)
	conn = sqlite3.connect(copy)
	try:
		conn.execute("UPDATE meta_data SET exclusion_status='1' "
					 "WHERE primary_accession LIKE 'Q00%'")
		flagged = conn.execute("SELECT COUNT(*) FROM meta_data "
							   "WHERE exclusion_status='1'").fetchone()[0]
		conn.commit()
	finally:
		conn.close()
	assert flagged > 0

	default = pal.load_cohort(sqlite3.connect(copy))
	kept = pal.load_cohort(sqlite3.connect(copy), drop_excluded=False)
	assert len(kept) - len(default) == flagged

	# And the flag is reachable from the command line.
	out = tmp_path / 'incl'
	lgr.main(['--db', copy, '--protein', 'G', '--out', str(out),
			  '--allele-source', 'features', '--no-tree', '--include-excluded'])
	report = open(os.path.join(str(out), 'report.md'), encoding='utf-8').read()
	assert 'Cohort:' in report


def test_alignment_route_defaults_to_all_sites_when_no_catalog(database, tmp_path):
	out = tmp_path / 'aln_auto_sites'
	lgr.main(['--db', database, '--protein', 'G', '--out', str(out),
			  '--allele-source', 'alignment', '--protein-start', '1', '--no-tree'])
	report = open(os.path.join(str(out), 'report.md'), encoding='utf-8').read()
	assert 'Cohort:' in report
