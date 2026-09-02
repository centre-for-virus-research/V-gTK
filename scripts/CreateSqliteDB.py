#!/usr/bin/env python3
import os
import re
import csv
import shutil
import sqlite3
import sys
import pandas as pd
from Bio import SeqIO
from os.path import join, normpath
from datetime import datetime
from argparse import ArgumentParser
from ExportRefListFromUpdateDb import load_reference_file_table
from clade_from_tree import assign_labels_from_tree
from accession_utils import normalise_accession
import segment_utils


class CreateSqliteDB:
	def __init__(
		self,
		meta_data,
		features,
		pad_aln,
		gene_info,
		m49_countries,
		m49_interm_region,
		m49_regions,
		m49_sub_regions,
		proj_settings,
		fasta_sequence_file,
		insertions,
		host_taxa_file,
		base_dir,
		output_dir,
		db_name,
		db_status,
		tree_file=None,
		iqtree_file=None,
		usher_tree=None,
		cluster_tsv=None,
		cluster_min_seq_id=None,
		filtered_ids_file=None,
		filtered_details_file=None,
		tree_manifest=None,
		reference_tsv=None,
		clade_assignments=None,
		update=False,
		update_db=None,
		batch_id=None,
		is_segmented=None,
	):
		# None = infer from the data (legacy behaviour). True/False = authoritative.
		self.is_segmented = is_segmented
		self.meta_data = meta_data
		self.features = features
		self.pad_aln = pad_aln
		self.gene_info = gene_info if gene_info not in {None, "", "None", "null"} else None
		self.m49_countries = m49_countries
		self.m49_interm_region = m49_interm_region
		self.m49_regions = m49_regions
		self.m49_sub_regions = m49_sub_regions
		self.proj_settings = proj_settings
		self.fasta_sequence_file = fasta_sequence_file
		self.insertions = insertions
		self.host_taxa_file = host_taxa_file
		self.base_dir = base_dir
		self.output_dir = output_dir
		self.db_name = db_name
		self.db_status = db_status
		self.tree_file = tree_file
		self.iqtree_file = iqtree_file
		self.usher_tree = usher_tree
		self.cluster_tsv = cluster_tsv
		self.cluster_min_seq_id = cluster_min_seq_id
		self.filtered_ids_file = filtered_ids_file
		self.filtered_details_file = filtered_details_file
		self.tree_manifest = tree_manifest
		self.reference_tsv = reference_tsv
		self.clade_assignments = clade_assignments
		self.update = bool(update)
		self.update_db = update_db
		self.batch_id = batch_id or datetime.now().strftime("batch_%Y%m%d_%H%M%S")
		# Columns fabricated for rows that do not exist yet (today: the cluster
		# placeholder). They are written on INSERT and deliberately left out of the
		# ON CONFLICT ... DO UPDATE SET list so a re-supplied accession keeps the
		# value an earlier run computed for it. {table: [column, ...]}
		self._insert_only_columns = {}

	@staticmethod
	def _normalize_segment_value(value):
		"""Delegates to :mod:`segment_utils` - the single normalisation authority.

		This used to scrape every digit out of the string: ``4.0`` became segment 40,
		and the polymerase names inverted - ``PB2`` (segment 1) became ``2`` and ``PB1``
		(segment 2) became ``1``. Missing still maps to ``""`` for this call site.
		"""
		normalised = segment_utils.normalise_segment(value)
		if normalised is None or normalised.casefold() in segment_utils.PANDAS_NULL_TOKENS:
			return ""
		return normalised

	@staticmethod
	def _read_tree_file(tree_path):
		if not tree_path:
			return None
		try:
			with open(tree_path, "r", encoding="utf-8") as handle:
				return handle.read().strip()
		except FileNotFoundError:
			return None

	@staticmethod
	def _load_tree_manifest(manifest_path):
		if not manifest_path or not os.path.isfile(manifest_path):
			return []
		rows = []
		with open(manifest_path, "r", encoding="utf-8") as f:
			reader = csv.DictReader(f, delimiter='\t')
			for row in reader:
				path_val = (row.get("path") or "").strip()
				if not path_val:
					continue
				rows.append({
					"source": (row.get("source") or "").strip() or "unknown",
					"name": (row.get("name") or "").strip() or None,
					"segment_key": (row.get("segment_key") or "").strip() or None,
					"path": path_val,
				})
		return rows

	@staticmethod
	def _segment_from_key(segment_key):
		if not segment_key:
			return None
		key = str(segment_key).strip()
		if not key:
			return None
		if key.isdigit():
			return key

		patterns = [
			r"(?:^|[_-])segment[_-]?(\d+)(?:$|[_-])",
			r"(?:^|[_-])seg[_-]?(\d+)(?:$|[_-])",
			r"(?:^|[_-])refset[_-]?(\d+)(?:$|[_-])",
			r"(?:^|[_-])(\d+)(?:$|[_-])",
		]
		for pat in patterns:
			m = re.search(pat, key, flags=re.IGNORECASE)
			if m:
				return m.group(1)
		return None

	def _load_filtered_ids(self):
		if not self.filtered_ids_file:
			return set()
		try:
			with open(self.filtered_ids_file, "r", encoding="utf-8") as f:
				raw = {line.strip() for line in f if line.strip()}
		except FileNotFoundError:
			return set()
		# These are matched against meta_data.primary_accession, which
		# _normalise_identity_columns has already made bare. Normalising only one
		# side of that comparison would silently stop the QC exclusion applying -
		# strictly worse than the un-normalised state it replaced. Carrying BOTH
		# spellings keeps the set a superset of what it used to be, so this can
		# only ever exclude more, never less.
		return raw | self._bare_variants(raw)

	@staticmethod
	def _bare_variants(ids):
		"""The canonical spelling of every id that has one.

		Composite ``<accession>_<segment>`` ids are left alone: an accession may
		itself contain an underscore (``NC_001542``), so splitting one apart is
		guesswork. They stay comparable because the composite is rebuilt from the
		already-normalised accession on the other side.
		"""
		bare = set()
		for value in ids:
			normalised = normalise_accession(value)
			if normalised and normalised != value:
				bare.add(normalised)
		return bare

	@staticmethod
	def _require_file(path, label):
		if not path or not os.path.isfile(path):
			raise FileNotFoundError(f"{label} file not found: {path}")

	@staticmethod
	def _read_tsv_required(path, required_columns, label, dtype=None, keep_na_strings=False):
		"""Read a TSV, optionally keeping NA-like strings as literal text.

		`keep_na_strings` is opt-in per call site rather than global on purpose.
		Some columns use 'NA' as a deliberate "not applicable" encoding - a real
		HCV matrix carries 274k host_validated='NA' cells that are meant to reach
		SQL as NULL - while others use it as a real name. Influenza's neuraminidase
		gene and segment are both literally called "NA", so for those inputs the
		default sentinel handling erases the value.
		"""
		df = pd.read_csv(path, sep="\t", dtype=dtype, keep_default_na=not keep_na_strings)
		missing = [c for c in required_columns if c not in df.columns]
		if missing:
			raise ValueError(f"{label} is missing required columns: {', '.join(missing)}")
		return df

	@staticmethod
	def _read_csv_required(path, required_columns, label, dtype=None):
		df = pd.read_csv(path, sep=",", dtype=dtype)
		missing = [c for c in required_columns if c not in df.columns]
		if missing:
			raise ValueError(f"{label} is missing required columns: {', '.join(missing)}")
		return df

	@staticmethod
	def _ensure_primary_accession(df, label, aliases=None):
		aliases = aliases or []
		if "primary_accession" in df.columns:
			return df
		for alias in aliases:
			if alias in df.columns:
				df["primary_accession"] = df[alias]
				return df
		raise ValueError(f"{label} is missing required columns: primary_accession")

	@staticmethod
	def _normalize_alignment_columns(df, label):
		if "sequence_id" in df.columns and "primary_accession" in df.columns:
			# Enforce sequence_id as the primary query key. 
			# Preserve the reference target name in alignment_name if not already present.
			if "alignment_name" not in df.columns:
				df["alignment_name"] = df["primary_accession"]
			df["primary_accession"] = df["sequence_id"]
		elif "primary_accession" not in df.columns and "sequence_id" in df.columns:
			df["primary_accession"] = df["sequence_id"]
		elif "primary_accession" in df.columns and "sequence_id" not in df.columns:
			df["sequence_id"] = df["primary_accession"]
		else:
			missing = [c for c in ["primary_accession", "sequence_id"] if c not in df.columns]
			if missing:
				raise ValueError(f"{label} is missing required columns: {', '.join(missing)}")
		return df

	@staticmethod
	def _normalize_filtered_reason(reason):
		text = "" if reason is None else str(reason).strip()
		if not text:
			return ""
		prefix = "alignment_filtering:"
		if text.startswith(prefix):
			detail = text[len(prefix):].strip()
		else:
			detail = text
		detail = re.sub(r"^In sequence #\d+\s+['\"][^'\"]+['\"]:\s*", "", detail, flags=re.IGNORECASE)
		if not detail:
			return "alignment_filtering"
		return f"alignment_filtering: {detail}"

	@classmethod
	def _normalize_summary_reason(cls, reason):
		text = "" if reason is None else str(reason).strip()
		if not text:
			return ""
		lower_text = text.lower()
		if (
			text.startswith("alignment_filtering:")
			or lower_text.startswith("in sequence #")
			or "unable to align" in lower_text
			or "number of seed matches" in lower_text
			or lower_text == "not significant blast hit"
			or "reference not present in master-projected reference_aln" in lower_text
			or lower_text.startswith("filtered against reference")
		):
			return cls._normalize_filtered_reason(text)
		return text

	def _load_filtered_details(self):
		reasons = {}
		if not self.filtered_details_file or not os.path.isfile(self.filtered_details_file):
			return reasons
		try:
			df = pd.read_csv(self.filtered_details_file, sep="\t", dtype=str).fillna("")
		except Exception:
			return reasons
		if "seq_name" not in df.columns:
			return reasons
		for _, row in df.iterrows():
			seq = str(row.get("seq_name", "")).strip()
			if not seq:
				continue
			err = str(row.get("error", "")).strip()
			ref = str(row.get("reference", "")).strip()
			if err:
				reason = self._normalize_filtered_reason(err)
			elif ref:
				reason = f"alignment_filtering: filtered against reference {ref}"
			else:
				reason = "alignment_filtering"
			reasons[seq] = reason
		# Same asymmetry as _load_filtered_ids: these keys are looked up with an
		# already-normalised accession, so the bare spelling has to be present too.
		for key, reason in list(reasons.items()):
			bare = normalise_accession(key)
			if bare and bare not in reasons:
				reasons[bare] = reason
		return reasons

	def _add_cluster_column(self, df_meta_data):
		if not self.cluster_tsv:
			return df_meta_data
		try:
			cluster_df = pd.read_csv(self.cluster_tsv, sep="\t", header=None, dtype=str)
		except FileNotFoundError:
			return df_meta_data
		if cluster_df.shape[1] < 2:
			return df_meta_data
		cluster_df = cluster_df.iloc[:, :2]
		cluster_df.columns = ["cluster_rep", "member"]
		cluster_map = dict(zip(cluster_df["member"], cluster_df["cluster_rep"]))
		# Mapped onto the already-normalised primary_accession below, so a cluster
		# TSV written with versioned members would otherwise map to nothing and
		# quietly drop every sequence's cluster assignment.
		for member, rep in list(cluster_map.items()):
			bare = normalise_accession(member)
			if bare and bare not in cluster_map:
				cluster_map[bare] = rep
		try:
			min_id = float(self.cluster_min_seq_id) if self.cluster_min_seq_id is not None else None
		except (TypeError, ValueError):
			min_id = None
		col_name = f"cluster_{int(round(min_id * 100))}pct" if min_id is not None else "cluster"
		if "primary_accession" in df_meta_data.columns:
			df_meta_data[col_name] = df_meta_data["primary_accession"].map(cluster_map)
		return df_meta_data

	@staticmethod
	def _normalize_reference_field(value):
		text = "" if value is None else str(value).strip()
		if text.lower() in {"", "na", "n/a", "nan", "none", "null"}:
			return ""
		return text

	def _load_reference_lookup(self):
		if not self.reference_tsv:
			return {}
		try:
			reference_df = load_reference_file_table(self.reference_tsv)
		except FileNotFoundError:
			return {}
		if "genotype" not in reference_df.columns:
			return {}
		reference_df = reference_df[["primary_accession", "genotype", "subtype"]].copy()
		if reference_df.empty:
			return {}
		for col in ["primary_accession", "genotype", "subtype"]:
			reference_df[col] = reference_df[col].map(self._normalize_reference_field)
		reference_df = reference_df[reference_df["primary_accession"] != ""]
		reference_df = reference_df.drop_duplicates(subset=["primary_accession"], keep="first")
		return {
			row["primary_accession"]: {
				"nearest_reference_genotype": row["genotype"],
				"nearest_reference_subtype": row["subtype"],
			}
			for _, row in reference_df.iterrows()
		}

	def _candidate_assignment_trees(self):
		"""(origin, newick) pairs that can carry queries + labelled references, most
		trustworthy first. The UShER tree holds every placed sample, so it is
		preferred; the IQ-TREE backbone only holds cluster representatives.

		`origin` is the provenance token stored in meta_data.genotype_origin. The
		agreed vocabulary has exactly two tree tokens, so UShER trees are
		'tree_usher' and every other tree - IQ-TREE, VeryFastTree, an unrecognised
		manifest source - is recorded under the generic 'tree_iqtree'.
		"""
		candidates = []
		manifest = self._load_tree_manifest(self.tree_manifest)

		def _add(path, origin):
			newick = self._read_tree_file(path)
			if newick:
				candidates.append((origin, newick))

		_add(self.usher_tree, "tree_usher")
		for entry in manifest:
			if entry.get("source") == "usher":
				_add(entry.get("path"), "tree_usher")
		_add(self.iqtree_file, "tree_iqtree")
		for entry in manifest:
			if entry.get("source") == "iqtree":
				_add(entry.get("path"), "tree_iqtree")
		_add(self.tree_file, "tree_iqtree")
		for entry in manifest:
			if entry.get("source") not in {"usher", "iqtree"}:
				_add(entry.get("path"), "tree_iqtree")
		return candidates

	def _tree_based_reference_labels(self, lookup):
		"""Assign genotype/subtype to query tips from their nearest labelled
		reference in the phylogenetic tree.

		Several trees routinely contain the same accession (one UShER tree per
		segment plus an IQ-TREE backbone). The first tree that says ANYTHING about
		an accession supplies BOTH of its fields, and later trees are not consulted
		for it. Merging the two fields independently - which is what this did -
		produced chimeric labels: a genotype from the UShER neighbourhood and a
		subtype from the IQ-TREE neighbourhood, a pair no reference in either tree
		carries. 69 of the HCV references have a genotype and no subtype, so "this
		tree gives a genotype but no subtype" is the normal case, not a corner.
		"""
		candidates = self._candidate_assignment_trees()
		if not candidates:
			return {}
		reference_labels = {
			acc: {
				"genotype": labels.get("nearest_reference_genotype", ""),
				"subtype": labels.get("nearest_reference_subtype", ""),
			}
			for acc, labels in lookup.items()
		}
		merged = {}
		for origin, newick in candidates:
			try:
				assignments = assign_labels_from_tree(newick, reference_labels)
			except Exception as exc:  # never let a malformed tree break the DB build
				print(f"[warn] Tree-based clade assignment skipped for one tree: {exc}")
				continue
			for acc, labels in assignments.items():
				if acc in merged:
					continue
				genotype = str(labels.get("genotype") or "").strip()
				subtype = str(labels.get("subtype") or "").strip()
				if not genotype and not subtype:
					continue
				merged[acc] = {"genotype": genotype, "subtype": subtype, "origin": origin}
		return merged

	def _load_clade_assignments(self):
		"""Precomputed assignments from CladeAssignment.py (the EPA-ng fallback,
		used when the run produced no usable tree). Columns: primary_accession,
		genotype, subtype."""
		path = self.clade_assignments
		if not path or not os.path.isfile(path):
			return {}
		assignments = {}
		with open(path, "r", newline="", encoding="utf-8") as handle:
			for row in csv.DictReader(handle, delimiter="\t"):
				accession = (row.get("primary_accession") or "").strip()
				if not accession:
					continue
				genotype = (row.get("genotype") or "").strip()
				subtype = (row.get("subtype") or "").strip()
				if genotype or subtype:
					assignments[accession] = {"genotype": genotype, "subtype": subtype}
		return assignments

	#: meta_data.genotype_origin / meta_data.subtype_origin vocabulary, in
	#: precedence order: the first source that supplies a value wins, and the value
	#: it supplied is the one stored. A curated reference-list entry is never
	#: overwritten by inference.
	GENOTYPE_ORIGINS = (
		"curated_reflist",
		"tree_usher",
		"tree_iqtree",
		"epa_placement",
		"blast_tophit",
		"gisaid_declared",
		"ncbi_declared",
	)
	GENOTYPE_ORIGIN_UNRESOLVED = "unresolved"

	@staticmethod
	def _vendor_declared_label_columns(df_meta_data):
		"""meta_data columns holding a genotype/subtype the *vendor* declared.

		Only columns that literally name a genotype or a subtype are read. Notably
		NOT `serotype`: influenza's H1N1 is a serotype, not a subtype in the sense
		nearest_reference_subtype carries, and quietly copying it in would relabel
		half a database.
		"""
		gisaid, ncbi = {}, {}
		for col in df_meta_data.columns:
			name = str(col).strip().lower()
			if name.startswith("nearest_reference"):
				continue
			if name.endswith("genotype"):
				field = "genotype"
			elif name.endswith("subtype"):
				field = "subtype"
			else:
				continue
			if name.startswith("gisaid"):
				gisaid.setdefault(field, col)
			elif name in {"genotype", "subtype"} or name.startswith(("ncbi", "genbank")):
				ncbi.setdefault(field, col)
		return {"gisaid_declared": gisaid, "ncbi_declared": ncbi}

	def _add_reference_columns(self, df_meta_data, df_aln):
		"""Resolve nearest_reference_genotype/_subtype and record WHERE each came from.

		Four columns are written: the two labels and, beside each, the source that
		produced it (see GENOTYPE_ORIGINS). Without the origin columns a stored
		'6a' on a reference curated as "genotype 6, subtype not assigned" is
		indistinguishable from a curated one - the subtype letter is inferred from
		whatever the sequence aligned against, and nothing in the database said so.

		The two fields are taken from the SAME source wherever that source has
		both, because they used to be resolved by two independent climbs and could
		be assembled out of two different neighbourhoods. A source that carries a
		genotype but no subtype still falls through for the subtype alone - that is
		deliberate (a curated subtype of "NA" means "not assigned", and inference
		may fill it) and the origin columns are what make the fall-through visible.
		"""
		lookup = self._load_reference_lookup()
		if not lookup or "primary_accession" not in df_meta_data.columns:
			return df_meta_data

		meta_accessions = df_meta_data["primary_accession"].fillna("").astype(str).str.strip()

		# Preferred source for queries: the phylogenetic tree neighbourhood.
		tree_labels = self._tree_based_reference_labels(lookup)

		# EPA-ng placement results, only produced when the run had no usable tree.
		epa_labels = self._load_clade_assignments()

		# Fallback for queries missing from the tree: the best BLAST hit, i.e. the
		# reference each query was aligned against (alignment_name).
		nearest_ref_map = {}
		if df_aln is not None and not df_aln.empty and {"primary_accession", "alignment_name"}.issubset(df_aln.columns):
			aln_map_df = df_aln[["primary_accession", "alignment_name"]].copy()
			aln_map_df["primary_accession"] = aln_map_df["primary_accession"].fillna("").astype(str).str.strip()
			aln_map_df["alignment_name"] = aln_map_df["alignment_name"].fillna("").astype(str).str.strip()
			aln_map_df = aln_map_df[(aln_map_df["primary_accession"] != "") & (aln_map_df["alignment_name"] != "")]
			aln_map_df = aln_map_df.drop_duplicates(subset=["primary_accession"], keep="first")
			nearest_ref_map = dict(zip(aln_map_df["primary_accession"], aln_map_df["alignment_name"]))

		# Vendor-declared labels, read positionally alongside the accessions.
		vendor_columns = self._vendor_declared_label_columns(df_meta_data)
		vendor_values = {}
		for origin, fields in vendor_columns.items():
			if not fields:
				continue
			vendor_values[origin] = {
				field: df_meta_data[col].fillna("").astype(str).str.strip().tolist()
				for field, col in fields.items()
			}

		def _clean(value):
			return self._normalize_reference_field(value)

		genotypes, subtypes, genotype_origins, subtype_origins = [], [], [], []
		for position, accession in enumerate(meta_accessions):
			candidates = []

			direct = lookup.get(accession)
			if direct:
				candidates.append((
					"curated_reflist",
					_clean(direct.get("nearest_reference_genotype", "")),
					_clean(direct.get("nearest_reference_subtype", "")),
				))

			tree = tree_labels.get(accession)
			if tree:
				candidates.append((tree.get("origin", "tree_iqtree"), _clean(tree.get("genotype")), _clean(tree.get("subtype"))))

			epa = epa_labels.get(accession)
			if epa:
				candidates.append(("epa_placement", _clean(epa.get("genotype")), _clean(epa.get("subtype"))))

			blast_ref = nearest_ref_map.get(accession, "")
			blast = lookup.get(blast_ref) if blast_ref else None
			if blast:
				candidates.append((
					"blast_tophit",
					_clean(blast.get("nearest_reference_genotype", "")),
					_clean(blast.get("nearest_reference_subtype", "")),
				))

			for origin in ("gisaid_declared", "ncbi_declared"):
				fields = vendor_values.get(origin)
				if not fields:
					continue
				candidates.append((
					origin,
					_clean(fields.get("genotype", [""] * len(meta_accessions))[position] if "genotype" in fields else ""),
					_clean(fields.get("subtype", [""] * len(meta_accessions))[position] if "subtype" in fields else ""),
				))

			genotype = subtype = ""
			genotype_origin = subtype_origin = self.GENOTYPE_ORIGIN_UNRESOLVED
			for origin, candidate_genotype, candidate_subtype in candidates:
				if not genotype and candidate_genotype:
					genotype, genotype_origin = candidate_genotype, origin
				if not subtype and candidate_subtype:
					subtype, subtype_origin = candidate_subtype, origin
				if genotype and subtype:
					break

			genotypes.append(genotype)
			subtypes.append(subtype)
			genotype_origins.append(genotype_origin)
			subtype_origins.append(subtype_origin)

		index = df_meta_data.index
		df_meta_data["nearest_reference_genotype"] = pd.Series(genotypes, index=index, dtype=object)
		df_meta_data["nearest_reference_subtype"] = pd.Series(subtypes, index=index, dtype=object)
		df_meta_data["genotype_origin"] = pd.Series(genotype_origins, index=index, dtype=object)
		df_meta_data["subtype_origin"] = pd.Series(subtype_origins, index=index, dtype=object)
		return df_meta_data

	def load_fasta(self):
		fasta_data = []
		for record in SeqIO.parse(self.fasta_sequence_file, "fasta"):
			fasta_data.append({"header": record.id, "sequence": str(record.seq)})
		return pd.DataFrame(fasta_data)

	@staticmethod
	def _normalize_db_status(db_status):
		s = (db_status or "").strip().lower()
		if s in {"new", "new db", "create", "created", "fresh"}:
			return "new db"
		if s in {"modified", "update", "updated", "changed", "last modified", "last updated"}:
			return "last updated"
		return db_status

	def _resolve_creation_type(self):
		"""What info.creation_type should say about THIS run.

		The pipeline never passed -ds, and the CLI default is the literal string
		"new db", so `if db_status` was always true and every --update run stamped
		itself as a fresh build: the shipped update-mode database carries two rows
		that both say 'new db' and nothing in it distinguishes "built once" from
		"built then updated". The run mode is the authority when the declared
		status contradicts it, and the contradiction is reported rather than
		absorbed.
		"""
		accurate = "last updated" if self.update else "new db"
		declared = self._normalize_db_status(self.db_status) if self.db_status else ""
		if not declared:
			return accurate
		if self.update and declared == "new db":
			print(
				"[CreateSqliteDB][warn] --update was requested but --db_status says 'new db'. "
				"Recording info.creation_type='last updated' instead: this run modified an "
				"existing database. Pass -ds 'last updated' from the caller to silence this."
			)
			return accurate
		return declared

	def _check_update_history(self, conn):
		"""Refuse or warn when the target's recorded history contradicts --update.

		info.creation_type is the only human-facing "how was this database built?"
		record, and update mode is the one path that runs against years of curated
		data. Pointing it at something that is not a database this pipeline built is
		not recoverable afterwards, so it is refused; a database with no recorded
		history (built before info was populated, or by hand) is allowed through
		with a loud warning because that is a legitimate legacy state.
		"""
		if not self.update:
			return
		if not self._table_exists(conn, "meta_data"):
			raise ValueError(
				f"--update was pointed at '{self.update_db}', which has no meta_data table. "
				"That is not a database this pipeline built, so there is nothing to update "
				"incrementally - run without --update to build it, or pass the right --update_db."
			)
		history = []
		if self._table_exists(conn, "info"):
			cols = self._table_columns(conn, "info")
			if "creation_type" in cols:
				history = [
					str(row[0] or "").strip()
					for row in conn.execute("SELECT creation_type FROM info").fetchall()
				]
		if not history:
			print(
				"[CreateSqliteDB][warn] --update target has no recorded build history "
				"(info.creation_type is empty or absent). Proceeding, but this database "
				"cannot say how it was built; every run from now on records itself."
			)
			return
		print(
			f"[CreateSqliteDB] --update target records {len(history)} previous build(s): "
			f"{history[:5]}{'...' if len(history) > 5 else ''}"
		)

	def _db_path(self):
		if self.update and not self.update_db:
			raise ValueError("--update requires --update_db path")
		return join(self.base_dir, self.output_dir, self.db_name + ".db")

	@staticmethod
	def _paths_equivalent(path_a, path_b):
		if not path_a or not path_b:
			return False
		return os.path.realpath(os.path.abspath(path_a)) == os.path.realpath(os.path.abspath(path_b))

	def _prepare_update_target_db(self, db_path):
		if not self.update:
			return
		if not self.update_db:
			raise ValueError("--update requires --update_db path")
		if not os.path.isfile(self.update_db):
			raise FileNotFoundError(f"update_db file not found: {self.update_db}")
		if self._paths_equivalent(db_path, self.update_db):
			return
		os.makedirs(os.path.dirname(db_path), exist_ok=True)
		shutil.copyfile(self.update_db, db_path)

	@staticmethod
	def _table_exists(conn, table):
		row = conn.execute("SELECT 1 FROM sqlite_master WHERE type='table' AND name=?", (table,)).fetchone()
		return row is not None

	@staticmethod
	def _table_columns(conn, table):
		if not CreateSqliteDB._table_exists(conn, table):
			return []
		rows = conn.execute(f"PRAGMA table_info({table})").fetchall()
		return [r[1] for r in rows]

	@staticmethod
	def _normalize_key_series(s):
		return s.fillna("").astype(str).str.strip()

	@staticmethod
	def _cluster_placeholder_columns(columns):
		return [c for c in columns if re.fullmatch(r"cluster_\d+pct", str(c or "").strip())]

	@staticmethod
	def _collect_nonblank_segments_from_df(df):
		if df is None or "segment" not in df.columns:
			return set()
		values = (
			df["segment"]
			.fillna("")
			.astype(str)
			.str.strip()
		)
		return {value for value in values.tolist() if value != ""}

	def _collect_nonblank_segments_from_table(self, conn, table):
		if not self._table_exists(conn, table):
			return set()
		if "segment" not in self._table_columns(conn, table):
			return set()
		rows = conn.execute(
			f"SELECT DISTINCT TRIM(COALESCE(CAST(segment AS TEXT), '')) FROM {table} "
			"WHERE TRIM(COALESCE(CAST(segment AS TEXT), '')) != ''"
		).fetchall()
		return {str(row[0]).strip() for row in rows if row and str(row[0]).strip()}

	def _should_force_unsegmented_segment_one(self, conn, dfs):
		"""Should blank segments be filled in with '1'?

		For a non-segmented virus every row is segment 1 by definition, so filling
		blanks is right. For a segmented one it is a fabrication: a blank means the
		segment could not be determined, and stamping '1' on it silently asserts
		PB2 for influenza.

		The two cases are indistinguishable from the data alone - a segmented build
		whose rows happen to all be segment 1 (a small test subset, or one where
		every other segment was excluded) looks exactly like a non-segmented one.
		So when the caller tells us, we believe them; inference is only the
		fallback for callers that do not pass the flag.
		"""
		# getattr: callers that build the object with __new__ (several tests do)
		# never run __init__, so the attribute may not exist.
		declared = getattr(self, 'is_segmented', None)
		if declared is not None:
			return not declared

		observed = set()
		for df in dfs:
			observed.update(self._collect_nonblank_segments_from_df(df))
		if conn is not None:
			for table in ["meta_data", "features", "sequence_alignment", "insertions"]:
				observed.update(self._collect_nonblank_segments_from_table(conn, table))
		return observed.issubset({"1"})

	@staticmethod
	def _force_segment_one_df(df, create_if_missing=False):
		if df is None:
			return df
		if "segment" not in df.columns:
			if not create_if_missing:
				return df
			df["segment"] = "1"
			return df
		segment = df["segment"].fillna("").astype(str).str.strip()
		df["segment"] = segment.mask(segment == "", "1")
		return df

	def _backfill_segment_one_in_db(self, conn, tables):
		for table in tables:
			if not self._table_exists(conn, table):
				continue
			columns = self._table_columns(conn, table)
			if "segment" not in columns:
				continue
			# For a non-segmented virus a blank segment and segment '1' are the same
			# row, so stamping '1' on the blank one can collide with a row that is
			# already there - and now that every upsert table carries a UNIQUE index
			# over its key, that collision is an IntegrityError that aborts the build
			# rather than a silent duplicate. The blank row is the stale one (it is
			# the state a database has before it learns about segments), so drop it
			# in favour of the row this run just wrote.
			key_cols = self._infer_key_cols(table, pd.DataFrame(columns=columns))
			other_keys = [c for c in key_cols if c != "segment" and c in columns]
			if other_keys:
				match = " AND ".join(
					f"TRIM(COALESCE(CAST(other.{self._quote_identifier(c)} AS TEXT), '')) = "
					f"TRIM(COALESCE(CAST({table}.{self._quote_identifier(c)} AS TEXT), ''))"
					for c in other_keys
				)
				cursor = conn.execute(
					f"DELETE FROM {table} WHERE TRIM(COALESCE(CAST(segment AS TEXT), '')) = '' "
					f"AND EXISTS (SELECT 1 FROM {table} AS other WHERE other.rowid <> {table}.rowid "
					f"AND TRIM(COALESCE(CAST(other.segment AS TEXT), '')) = '1' AND {match})"
				)
				if cursor.rowcount:
					print(
						f"[CreateSqliteDB] Dropped {cursor.rowcount} blank-segment row(s) from '{table}' "
						f"superseded by the segment '1' row with the same key {other_keys}"
					)
			conn.execute(
				f"UPDATE {table} SET segment='1' WHERE TRIM(COALESCE(CAST(segment AS TEXT), '')) = ''"
			)

	def _fetch_existing_keys(self, conn, table, key_cols):
		if not self._table_exists(conn, table):
			return set()
		df = pd.read_sql_query(f"SELECT {', '.join(key_cols)} FROM {table}", conn)
		for c in key_cols:
			if c not in df.columns:
				raise ValueError(f"DB table '{table}' is missing expected key column '{c}'")
			df[c] = self._normalize_key_series(df[c])
		if len(key_cols) == 1:
			return set(df[key_cols[0]].tolist())
		return set(map(tuple, df[key_cols].itertuples(index=False, name=None)))

	@staticmethod
	def _dedupe_incoming_df(df, key_cols):
		before = len(df)
		df2 = df.drop_duplicates(subset=key_cols, keep="first")
		return df2, (before - len(df2))

	def _align_df_to_existing_schema(self, conn, table, df):
		existing_cols = self._table_columns(conn, table)
		if not existing_cols:
			return df

		incoming_cols = list(df.columns)
		extra_cols = [c for c in incoming_cols if c not in existing_cols]
		if extra_cols:
			print(
				f"[CreateSqliteDB][warn] Incoming '{table}' has columns not present in existing DB schema; "
				f"dropping: {extra_cols}"
			)
			df = df.drop(columns=extra_cols)

		# A dropped column that differs from a stored one only by case is not an
		# extra column at all - SQLite identifiers are case-insensitive, so the
		# incoming values were headed for that stored column and have just been
		# thrown away. This has to raise BEFORE the tolerant path below, which
		# would otherwise read 'Segment' as "absent from this batch, keep what is
		# stored" and write the row with its segment silently unchanged.
		existing_by_case = {str(c).casefold(): c for c in existing_cols}
		case_collisions = [
			(c, existing_by_case[str(c).casefold()])
			for c in extra_cols
			if str(c).casefold() in existing_by_case
		]
		if case_collisions:
			pairs = "; ".join(f"incoming {inc!r} vs stored {db!r}" for inc, db in case_collisions)
			raise ValueError(
				f"Incoming '{table}' dataframe is missing columns required by existing DB schema: "
				f"{[db for _, db in case_collisions]}. These differ from the incoming columns only by "
				f"case ({pairs}), and SQLite identifiers are case-insensitive, so the incoming values "
				"would be dropped rather than stored. Rename one column of each pair - in the incoming "
				"table or in the database - so the two spellings match exactly."
			)

		missing_cols = [c for c in existing_cols if c not in df.columns]
		if self.update and table in self.UPSERT_TABLES:
			# A column the incoming batch does not carry must LEAVE THE STORED VALUE
			# ALONE. Fabricating it here and handing it to INSERT OR REPLACE is how a
			# run without --cluster_tsv wrote the literal string 'NA- see tree' over
			# the cluster representative an earlier run had computed (100 such rows in
			# test_out/update_test/rabv-jul0425-update-test.db). The upsert now omits
			# absent columns from the UPDATE half instead, so they survive untouched.
			insert_only = []
			for cluster_col in self._cluster_placeholder_columns(missing_cols):
				# Still fabricated - but only for rows being INSERTed, i.e. accessions
				# the database has never seen. An existing row keeps its own value.
				df[cluster_col] = "NA- see tree"
				insert_only.append(cluster_col)
			self._record_insert_only_columns(table, insert_only)
			preserved = [c for c in existing_cols if c not in df.columns]
			if preserved:
				shown = preserved[:10]
				suffix = "" if len(preserved) == len(shown) else f" (+{len(preserved) - len(shown)} more)"
				print(
					f"[CreateSqliteDB] Incoming '{table}' does not carry {len(preserved)} column(s) "
					f"present in the DB; their stored values are preserved: {shown}{suffix}"
				)
			return df[[c for c in existing_cols if c in df.columns]]
		if missing_cols:
			# External columns are now namespaced at the merge step, so a DB built
			# before that change holds their old raw spellings. Without this hint the
			# operator reads the bare list as "the export lost columns" rather than
			# "the merge renamed them", and retries against a DB that _ensure_update_columns
			# has already half-altered.
			from merge_into_gB_matrix import NormalizeAndMerge

			incoming_slugs = {NormalizeAndMerge.slugify_column(c): c for c in df.columns}
			# The namespaced-column branch used to index incoming_slugs with a slug it
			# had only matched against the SUFFIX of some other slug, so the lookup
			# raised KeyError and the ValueError this hint exists to explain never
			# reached the operator.
			renamed = {}
			for old_name in missing_cols:
				slug = NormalizeAndMerge.slugify_column(old_name)
				match = incoming_slugs.get(slug)
				if match is None:
					match = next(
						(value for key, value in incoming_slugs.items() if slug == key.split('_', 1)[-1]),
						None,
					)
				if match is not None:
					renamed[old_name] = match
			hint = ""
			if renamed or any('_' in c for c in df.columns):
				hint = (
					" These look like external columns that are now namespaced at the merge "
					"step (e.g. 'Lineage' is produced as '<source>_lineage'). A database built "
					"before that change cannot be updated in place - rebuild it, or rename the "
					"columns in the existing DB to match."
				)
			raise ValueError(
				f"Incoming '{table}' dataframe is missing columns required by existing DB schema: {missing_cols}.{hint}"
			)

		return df[existing_cols]

	@staticmethod
	def _quote_identifier(name):
		"""Quote a SQL identifier. Never interpolate a bare column name.

		`ALTER TABLE t ADD COLUMN HA INSDC_Upload TEXT` does not fail - SQLite
		parses `HA` as the name and `INSDC_Upload TEXT` as the type, silently
		creating the wrong column. GISAID ships ten such headers.
		"""
		return '"' + str(name).replace('"', '""') + '"'

	@staticmethod
	def assert_no_case_insensitive_duplicates(columns, label):
		"""Fail loudly, and by name, before SQLite fails cryptically.

		SQLite column identifiers are case-insensitive while pandas labels are
		case-sensitive, so a frame carrying both `segment` and `Segment` is legal
		in pandas and illegal in SQLite. Left to sqlite the build dies with
		`duplicate column name: Segment` and no indication of which two columns
		or which input produced them.
		"""
		seen, collisions = {}, []
		for position, col in enumerate(columns, start=1):
			key = str(col).casefold()
			if key in seen:
				collisions.append(f"{seen[key][0]!r} (col {seen[key][1]}) vs {col!r} (col {position})")
			else:
				seen[key] = (col, position)
		if collisions:
			raise ValueError(
				f"'{label}' has columns that differ only by case, which SQLite cannot "
				f"represent: {'; '.join(collisions)}. Map or rename one of each pair - "
				f"external columns should be namespaced at the merge step."
			)

	#: meta_data columns whose value may legitimately BE the string "NA".
	#: Influenza segment 6 is neuraminidase, so `segment_name` is literally 'NA'
	#: and `segment`/`segment_declared` can be too. Everywhere else in meta_data
	#: 'NA' is a deliberate "not applicable" encoding that must reach SQL as NULL
	#: - a real HCV matrix carries 274k such host_validated cells - which is why
	#: this is a per-column exemption rather than a file-wide keep_default_na.
	NA_IS_A_VALUE_COLUMNS = frozenset({
		"segment", "segment_name", "segment_declared", "segment_validated",
	})

	def _read_meta_data_tsv(self):
		"""Read meta_data losslessly, then restore null semantics per column.

		Reading with pandas' defaults would erase every neuraminidase row's
		`segment_name` - the exact bug this change set exists to fix, reintroduced
		by the columns it adds.
		"""
		frame = self._read_tsv_required(
			self.meta_data, ["primary_accession"], "meta_data",
			dtype=str, keep_na_strings=True,
		)
		sentinels = {
			"", "#N/A", "#N/A N/A", "#NA", "-1.#IND", "-1.#QNAN", "-NaN", "-nan",
			"1.#IND", "1.#QNAN", "<NA>", "N/A", "NA", "NULL", "NaN", "None",
			"n/a", "nan", "null",
		}
		# One vectorised mask over the whole frame rather than a per-column loop.
		# .replace() would route through pandas' deprecated silent-downcasting
		# path, and assigning column by column fragments a 101-column x 3.9M-row
		# frame badly enough that pandas warns about it.
		#
		# For the segment columns only a blank is missing, so that influenza's
		# neuraminidase - genuinely named "NA" - survives. Everywhere else the
		# full sentinel set applies, because there 'NA' means "not applicable"
		# and is meant to reach SQL as NULL.
		protected = frame.columns.isin(list(self.NA_IS_A_VALUE_COLUMNS))
		condition = frame.isin(sentinels)
		if protected.any():
			condition.loc[:, protected] = frame.loc[:, protected].eq("")
		return frame.mask(condition)

	def _ensure_update_columns(self, conn, table, df):
		if not self.update or not self._table_exists(conn, table):
			return
		self.assert_no_case_insensitive_duplicates(df.columns, f"incoming {table}")
		# Track additions as we go: reading the schema once and never updating it
		# meant a repeated label issued the same ALTER twice, and the second raised.
		existing = {str(c).casefold() for c in self._table_columns(conn, table)}
		for col in df.columns:
			if str(col).casefold() in existing:
				continue
			conn.execute(f"ALTER TABLE {table} ADD COLUMN {self._quote_identifier(col)} TEXT")
			existing.add(str(col).casefold())

	def _resolve_key_cols_for_existing_schema(self, conn, table, key_cols, incoming_cols):
		existing_cols = set(self._table_columns(conn, table))
		if not existing_cols:
			return key_cols

		resolved = [c for c in key_cols if c in existing_cols and c in incoming_cols]
		if resolved != key_cols:
			print(
				f"[CreateSqliteDB][warn] Adjusted key columns for '{table}' to match existing DB schema: "
				f"requested={key_cols}, resolved={resolved}"
			)

		if resolved:
			return resolved

		fallback_candidates = [
			"primary_accession", "accession", "header", "name", "m49_code", "id"
		]
		for c in fallback_candidates:
			if c in existing_cols and c in incoming_cols:
				print(
					f"[CreateSqliteDB][warn] Falling back to key column '{c}' for table '{table}' "
					"due to schema mismatch."
				)
				return [c]

		raise ValueError(
			f"Cannot resolve compatible key columns for table '{table}'. "
			f"Requested keys={key_cols}; existing columns={sorted(existing_cols)}; "
			f"incoming columns={sorted(list(incoming_cols))}"
		)

	def _infer_key_cols(self, table, df):
		cols = list(df.columns)
		if table == "meta_data":
			return ["primary_accession", "segment"] if "segment" in cols else ["primary_accession"]
		if table == "sequences":
			return ["header", "segment"] if "segment" in cols else ["header"]
		if table == "sequence_alignment":
			key = ["primary_accession"]
			if "alignment_name" in cols:
				key.append("alignment_name")
			if "segment" in cols:
				key.append("segment")
			return key
		if table == "features":
			key = ["accession"]
			if "cds_start_OG_seq" in cols and "cds_end_OG_seq" in cols:
				key.extend(["cds_start_OG_seq", "cds_end_OG_seq"])
			else:
				key.extend(["cds_start", "cds_end"])
			key.append("product")
			if "segment" in cols:
				key.append("segment")
			missing = [c for c in key if c not in cols]
			if not missing:
				return key
			if "primary_accession" in cols and "segment" in cols:
				return ["primary_accession", "segment"]
			if "primary_accession" in cols:
				return ["primary_accession"]
			if "accession" in cols and "segment" in cols:
				return ["accession", "segment"]
			if "accession" in cols:
				return ["accession"]
			return cols[:1]
		if table == "insertions":
			key = [c for c in ["primary_accession", "reference", "insertion"] if c in cols]
			if "segment" in cols:
				key.append("segment")
			return key if key else cols[:1]
		if table == "host_taxa":
			for c in ["taxonomy_id", "host_taxa_id", "tax_id", "id"]:
				if c in cols:
					return [c]
			return [cols[0]] if cols else []
		if table == "genes":
			for c in ["name", "gene_name", "description"]:
				if c in cols:
					return [c]
			return cols[:1]
		if table in ["m49_country", "m49_intermediate", "m49_regions", "m49_sub_regions", "project_settings"]:
			return cols[:1] if cols else []
		return cols[:1] if cols else []

	#: The columns whose value IS the pipeline's identity key, per table. Every one
	#: of them is joined - directly or through valid_pairs - against
	#: meta_data.primary_accession, so all of them have to carry the same (bare)
	#: spelling of a GenBank name.
	#:
	#: meta_data.accession_version is deliberately absent, and must stay absent. It
	#: is the ONE column whose job is to carry the version, so that GenBankFetcher
	#: can notice on an update run that a record has been revised (PV547761.1 ->
	#: PV547761.2) and re-fetch it. Normalising it would make every record look
	#: unrevised for ever.
	#:
	#: sequence_alignment.sequence_id is listed because _normalize_alignment_columns
	#: makes primary_accession a verbatim copy of it; normalising one and not the
	#: other would leave a single row disagreeing with itself about its own name.
	#: 'primary_accession' appears under features because the features TSV is
	#: allowed to name its key column either way (see feat_acc_col in create_db).
	IDENTITY_COLUMNS = {
		"meta_data": ("primary_accession", "locus", "gi_number"),
		"sequence_alignment": ("primary_accession", "sequence_id", "alignment_name"),
		"features": ("accession", "primary_accession", "master_ref_accession", "reference_accession"),
		"sequences": ("header",),
		"insertions": ("primary_accession",),
	}

	#: How many offending values to name in the warning before trailing off. Enough
	#: to identify the producer, few enough not to bury the message under a batch
	#: of 3.9M rows.
	IDENTITY_WARNING_SAMPLE = 5

	def _normalise_identity_columns(self, df, table, columns=None):
		"""Force this table's identity columns to the bare accession spelling.

		This class used to do NO accession normalisation at all: whatever spelling
		arrived was stored verbatim after .str.strip(). That was correct only by
		luck - every upstream producer happens to emit bare accessions - and it
		failed silently in three different ways the moment one of them emitted
		'PV547761.1' instead of 'PV547761':

		* a versioned accession in features.tsv or the padded alignment, against a
		  bare meta_data, is not in valid_pairs, so the row is filtered out with no
		  message. For features that can empty the table;
		* sequences.header comes straight from Bio.SeqIO's record.id and keeps
		  whatever the FASTA carried, so the sequence row vanished while its
		  meta_data row stayed - a record with metadata and no sequence;
		* on an --update run the UNIQUE indexes are over the raw column, so bare and
		  versioned are two distinct keys: INSERT ... ON CONFLICT found no conflict
		  and simply inserted, giving a DUPLICATE row per upsert table instead of an
		  overwrite. Nothing logged it and the run exited 0.

		So the value is repaired here, at ingest, where the table is read - and the
		repair is REPORTED. A version leak upstream is a bug in the producer; being
		quietly absorbed here is how it would survive to the next release. Only a
		real version strip warns: leading/trailing whitespace is tidied silently,
		as the surrounding .str.strip() calls already do.

		Normalisation goes through accession_utils.normalise_accession, so it is
		conservative: only a string that is entirely an accession followed by
		'.<digits>' loses its suffix. A strain name, a GISAID id, a cluster label or
		a segment group name containing a dot passes through untouched, which
		split('.')[0] would not have done.
		"""
		if df is None or df.empty:
			return df

		columns = self.IDENTITY_COLUMNS.get(table, ()) if columns is None else columns
		for column in columns:
			if column not in df.columns:
				continue

			values = df[column].tolist()
			new_values = []
			version_strips = []
			touched = False
			for value in values:
				# A missing value stays missing. Coercing it to '' here would turn a
				# NULL locus into an empty string in every existing database.
				if value is None or (not isinstance(value, str) and pd.isna(value)):
					new_values.append(value)
					continue
				text = str(value)
				stripped = text.strip()
				bare = normalise_accession(text)
				replacement = stripped if bare is None else bare
				if replacement != text:
					touched = True
				if replacement != stripped:
					version_strips.append((stripped, replacement))
				new_values.append(replacement)

			if touched:
				df[column] = pd.Series(new_values, index=df.index, dtype=object)

			if version_strips:
				sample = ", ".join(
					f"{raw} -> {new}" for raw, new in version_strips[:self.IDENTITY_WARNING_SAMPLE]
				)
				if len(version_strips) > self.IDENTITY_WARNING_SAMPLE:
					sample += ", ..."
				print(
					f"[CreateSqliteDB][warn] {len(version_strips)} value(s) in {table}.{column} "
					f"arrived carrying a GenBank version suffix and were normalised to the bare "
					f"accession this pipeline joins on: {sample}. Fix the producer of that column - "
					"a versioned identity key joins against nothing, so the row would have been "
					"dropped on a fresh build and duplicated on an --update run. "
					"meta_data.accession_version is the only column allowed to carry a version."
				)
		return df

	#: Tables written with an upsert rather than an append. Every one of them needs
	#: a UNIQUE index over its key or the "replace" half never fires.
	UPSERT_TABLES = frozenset({"meta_data", "sequence_alignment", "features", "insertions", "sequences"})

	#: Index names are pinned per table so a database that already carries the
	#: original three is not given a second, identical index under a new name.
	#: sequences/insertions had none at all: they are in UPSERT_TABLES and were
	#: written with INSERT OR REPLACE, which with nothing to conflict on degrades
	#: to a plain INSERT - so replaying an identical update grew both tables by one
	#: row per record.
	UPSERT_INDEX_NAMES = {
		"features": "idx_features_upsert",
		"sequence_alignment": "idx_seq_alignment_upsert",
		"meta_data": "idx_metadata_upsert",
		"sequences": "idx_sequences_upsert",
		"insertions": "idx_insertions_upsert",
	}

	def _record_insert_only_columns(self, table, columns):
		store = getattr(self, "_insert_only_columns", None)
		if store is None:
			store = {}
			self._insert_only_columns = store
		store[table] = list(columns)

	def _insert_only_columns_for(self, table):
		return list((getattr(self, "_insert_only_columns", None) or {}).get(table, []))

	def _blank_out_null_key_columns(self, conn, table, key_cols):
		"""Make stored NULL key values reachable by the upsert.

		SQLite treats NULLs as DISTINCT in a UNIQUE index, so a pre-existing row
		with a NULL key column can never conflict with anything: INSERT OR REPLACE
		appends beside it and the stale row survives forever. Incoming keys are
		normalised to '' by _normalize_key_series, so normalising stored NULLs the
		same way is what lets the two meet - and it matches how the rest of this
		file compares keys (TRIM(COALESCE(col, ''))).

		A row whose stored key is NULL where the batch carries a real value (a
		database that predates the `segment` column, say) is still a different key
		afterwards and is NOT silently merged; it is reported so an operator can
		see that those rows will not be superseded.
		"""
		for col in key_cols:
			quoted = self._quote_identifier(col)
			row = conn.execute(f"SELECT COUNT(*) FROM {table} WHERE {quoted} IS NULL").fetchone()
			null_count = int(row[0]) if row else 0
			if not null_count:
				continue
			conn.execute(f"UPDATE {table} SET {quoted} = '' WHERE {quoted} IS NULL")
			print(
				f"[CreateSqliteDB][warn] '{table}'.{col} held {null_count} NULL value(s) in the "
				f"upsert key {key_cols}; normalised to '' so the UNIQUE index can see them. "
				"Rows whose key is blank are only superseded by an incoming row with a blank key."
			)

	def _describe_duplicate_keys(self, conn, table, key_cols, limit=5):
		cols_sql = ", ".join(self._quote_identifier(c) for c in key_cols)
		try:
			rows = conn.execute(
				f"SELECT {cols_sql}, COUNT(*) FROM {table} GROUP BY {cols_sql} "
				f"HAVING COUNT(*) > 1 ORDER BY COUNT(*) DESC LIMIT {int(limit)}"
			).fetchall()
		except sqlite3.Error:
			return []
		return rows

	def _create_table_unique_indexes(self, conn, table, key_cols):
		"""Create the UNIQUE index that makes the upsert an upsert.

		`key_cols` must be the key the merge actually resolved (see _infer_key_cols
		and _resolve_key_cols_for_existing_schema) - an index over a different set
		would not be a usable ON CONFLICT target, and the write would silently fall
		back to appending.

		Returns the index name, or None when no index could be built.
		"""
		index_name = self.UPSERT_INDEX_NAMES.get(table)
		if not index_name:
			return None
		existing_cols = self._table_columns(conn, table)
		if not existing_cols:
			return None
		existing_lower = {str(c).casefold(): c for c in existing_cols}
		resolved = [existing_lower[str(c).casefold()] for c in (key_cols or []) if str(c).casefold() in existing_lower]
		if not resolved:
			return None

		self._blank_out_null_key_columns(conn, table, resolved)
		cols_sql = ", ".join(self._quote_identifier(c) for c in resolved)
		try:
			conn.execute(
				f"CREATE UNIQUE INDEX IF NOT EXISTS {self._quote_identifier(index_name)} "
				f"ON {table} ({cols_sql})"
			)
		except sqlite3.IntegrityError as exc:
			# These indexes were only ever created in update mode, so a database
			# built before this change has none and the first update is the first
			# time the constraint exists. Say exactly which rows block it instead of
			# letting "UNIQUE constraint failed" escape with no context.
			examples = self._describe_duplicate_keys(conn, table, resolved)
			rendered = "; ".join(
				f"{tuple(row[:-1])} x{row[-1]}" for row in examples
			) or "(could not be listed)"
			raise ValueError(
				f"Cannot create UNIQUE index '{index_name}' on {table}({', '.join(resolved)}): the "
				f"existing database already holds rows that share a key, e.g. {rendered}. Update mode "
				"needs this index to replace rows instead of appending them, so the duplicates have to "
				"go first - de-duplicate the table (keep one row per key) and re-run, or rebuild the "
				f"database from scratch. Underlying error: {exc}"
			) from exc
		return index_name

	@staticmethod
	def _unique_index_for_keys(conn, table, key_cols):
		"""Find a non-partial UNIQUE index whose columns are exactly `key_cols`.

		Returned in the index's own column order, which is what ON CONFLICT needs.
		"""
		wanted = {str(c).casefold() for c in key_cols}
		if not wanted:
			return None, []
		try:
			index_rows = conn.execute(f'PRAGMA index_list("{table}")').fetchall()
		except sqlite3.Error:
			return None, []
		for row in index_rows:
			name, is_unique = row[1], row[2]
			partial = row[4] if len(row) > 4 else 0
			if not is_unique or partial:
				continue
			info = conn.execute(f'PRAGMA index_info("{name}")').fetchall()
			cols = [r[2] for r in info]
			if any(c is None for c in cols):
				continue  # index over an expression
			if {str(c).casefold() for c in cols} == wanted:
				return name, cols
		return None, []

	#: Columns whose blank value is a statement, not a gap: an update must be able
	#: to CLEAR them. Everywhere else a blank incoming cell is treated as "this
	#: batch says nothing" and the stored value is kept - a thin batch must not
	#: erase curated metadata (a real HCV matrix carries 274k host_validated='NA'
	#: cells, all of which reach this point as NULL).
	CLEARABLE_ON_UPDATE_COLUMNS = frozenset({"exclusion_status", "exclusion_criteria"})

	def _upsert_dataframe(self, conn, table, df, key_cols):
		"""Column-wise upsert: replace what this batch supplies, keep the rest.

		INSERT OR REPLACE deletes the conflicting row and inserts the incoming one
		wholesale, so any column that is absent or blank in this batch silently
		erases whatever an earlier run had stored. ON CONFLICT ... DO UPDATE SET
		touches only the columns present in the batch, and only when they carry a
		value.
		"""
		columns = list(df.columns)
		placeholders = ", ".join(["?"] * len(columns))
		sql_cols = ", ".join(self._quote_identifier(c) for c in columns)
		payload = [tuple(x) for x in df.values]

		index_name, index_cols = self._unique_index_for_keys(conn, table, key_cols)
		if not index_name:
			# No usable conflict target: keep the historical behaviour rather than
			# failing the run, but say so - this is the state in which an update
			# duplicates rows instead of replacing them.
			print(
				f"[CreateSqliteDB][warn] No UNIQUE index over {key_cols} on '{table}'; falling back to "
				"INSERT OR REPLACE, which cannot preserve columns this batch does not carry."
			)
			conn.executemany(
				f"INSERT OR REPLACE INTO {table} ({sql_cols}) VALUES ({placeholders})", payload
			)
			return len(df)

		key_lower = {str(c).casefold() for c in index_cols}
		insert_only = {str(c).casefold() for c in self._insert_only_columns_for(table)}
		clearable = {c.casefold() for c in self.CLEARABLE_ON_UPDATE_COLUMNS}

		assignments = []
		for col in columns:
			lower = str(col).casefold()
			if lower in key_lower or lower in insert_only:
				continue
			quoted = self._quote_identifier(col)
			if lower in clearable:
				assignments.append(f"{quoted} = excluded.{quoted}")
			else:
				assignments.append(
					f"{quoted} = CASE WHEN excluded.{quoted} IS NULL "
					f"OR TRIM(CAST(excluded.{quoted} AS TEXT)) = '' "
					f"THEN {self._quote_identifier(table)}.{quoted} ELSE excluded.{quoted} END"
				)
		conflict_cols = ", ".join(self._quote_identifier(c) for c in index_cols)
		if assignments:
			sql = (
				f"INSERT INTO {table} ({sql_cols}) VALUES ({placeholders}) "
				f"ON CONFLICT ({conflict_cols}) DO UPDATE SET {', '.join(assignments)}"
			)
		else:
			sql = (
				f"INSERT INTO {table} ({sql_cols}) VALUES ({placeholders}) "
				f"ON CONFLICT ({conflict_cols}) DO NOTHING"
			)
		conn.executemany(sql, payload)
		return len(df)

	def merge_table_append_nonredundant(self, conn, df, table, key_cols=None, update_exclusions=None):
		if df is None:
			return 0
		df = df.copy()
		key_cols = key_cols if key_cols is not None else self._infer_key_cols(table, df)
		
		if self.update and self._table_exists(conn, table):
			self._ensure_update_columns(conn, table, df)
			df = self._align_df_to_existing_schema(conn, table, df)
			key_cols = self._resolve_key_cols_for_existing_schema(conn, table, key_cols, set(df.columns))
			if table in self.UPSERT_TABLES:
				self._create_table_unique_indexes(conn, table, key_cols)

		for c in key_cols:
			if c not in df.columns:
				raise ValueError(f"Incoming '{table}' dataframe missing key column '{c}'")
			df[c] = self._normalize_key_series(df[c])
			
		df, dropped_internal = self._dedupe_incoming_df(df, key_cols)
		if dropped_internal:
			print(f"[CreateSqliteDB] Dropped {dropped_internal} duplicate incoming rows in '{table}' by key {key_cols}")

		# Both write paths below go through SQLite's case-insensitive identifiers,
		# so check once here rather than letting either fail opaquely.
		self.assert_no_case_insensitive_duplicates(df.columns, f"incoming {table}")

		if not self.update:
			df.to_sql(table, conn, if_exists="replace", index=False)
			# Build the index on a fresh build too. to_sql(if_exists='replace') drops
			# the table and its indexes, so this has to run after the write - and
			# without it the very first --update against a new database is the first
			# time the constraint exists, which is precisely when a pre-existing
			# violation would surface as an unexplained CREATE failure.
			if table in self.UPSERT_TABLES:
				self._create_table_unique_indexes(conn, table, key_cols)
			return len(df)

		if not self._table_exists(conn, table):
			df.to_sql(table, conn, if_exists="replace", index=False)
			if table in self.UPSERT_TABLES:
				self._create_table_unique_indexes(conn, table, key_cols)
			print(f"[CreateSqliteDB] Created table '{table}' with {len(df)} rows (table did not exist)")
			return len(df)

		if table in self.UPSERT_TABLES:
			self._upsert_dataframe(conn, table, df, key_cols)
			print(f"[CreateSqliteDB] Atomically upserted {len(df)} rows into '{table}' by key {key_cols}")
			return len(df)

		existing_keys = self._fetch_existing_keys(conn, table, key_cols)
		if len(key_cols) == 1:
			key = key_cols[0]
			new_mask = ~df[key].isin(existing_keys)
			dup_keys = df.loc[~new_mask, key].tolist()
		else:
			incoming_keys = list(map(tuple, df[key_cols].itertuples(index=False, name=None)))
			keep_flags = [k not in existing_keys for k in incoming_keys]
			new_mask = pd.Series(keep_flags, index=df.index)
			dup_keys = [k for k, keep in zip(incoming_keys, keep_flags) if not keep]
			
		df_new = df.loc[new_mask].copy()
		if not df_new.empty:
			df_new.to_sql(table, conn, if_exists="append", index=False)
			print(f"[CreateSqliteDB] Appended {len(df_new)} new rows into '{table}' (non-redundant)")
		else:
			print(f"[CreateSqliteDB] No new rows to append into '{table}' (all duplicates)")
			
		if update_exclusions is not None and dup_keys:
			now_str = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
			for k in dup_keys:
				update_exclusions.append({
					"batch_id": self.batch_id, 
					"table_name": table, 
					"key": str(k), 
					"reason": "duplicate_key_in_db", 
					"date": now_str
				})
		return len(df_new)

	@staticmethod
	def _table_row_count(conn, table):
		if not CreateSqliteDB._table_exists(conn, table):
			return 0
		row = conn.execute(f"SELECT COUNT(*) FROM {table}").fetchone()
		return int(row[0]) if row else 0

	@staticmethod
	def _count_distinct_nonexcluded_query_accessions(conn):
		if not CreateSqliteDB._table_exists(conn, "meta_data"):
			return 0
		cols = [row[1] for row in conn.execute("PRAGMA table_info(meta_data)").fetchall()]
		if "primary_accession" not in cols:
			return 0

		where = ["primary_accession IS NOT NULL", "TRIM(CAST(primary_accession AS TEXT)) != ''"]
		if "accession_type" in cols:
			where.append("LOWER(TRIM(COALESCE(CAST(accession_type AS TEXT), ''))) NOT IN ('reference', 'master')")
		if "exclusion_status" in cols:
			where.append(
				"LOWER(TRIM(COALESCE(CAST(exclusion_status AS TEXT), ''))) IN ('', '0', 'false', 'no', 'na', 'none', 'nan')"
			)

		sql = f"SELECT COUNT(DISTINCT primary_accession) FROM meta_data WHERE {' AND '.join(where)}"
		row = conn.execute(sql).fetchone()
		return int(row[0]) if row else 0

	def create_db(self):
		self._require_file(self.meta_data, "meta_data")
		self._require_file(self.features, "features")
		self._require_file(self.pad_aln, "pad_aln")
		if self.gene_info is not None:
			self._require_file(self.gene_info, "gene_info")
		self._require_file(self.m49_countries, "m49_countries")
		self._require_file(self.m49_interm_region, "m49_interm_region")
		self._require_file(self.m49_regions, "m49_regions")
		self._require_file(self.m49_sub_regions, "m49_sub_regions")
		self._require_file(self.proj_settings, "proj_settings")
		self._require_file(self.fasta_sequence_file, "fasta_sequences")
		self._require_file(self.insertions, "insertions")
		self._require_file(self.host_taxa_file, "host_taxa_file")

		excluded_records = []
		update_exclusions = []
		filtered_ids = self._load_filtered_ids()
		filtered_details = self._load_filtered_details()

		df_meta_data = self._read_meta_data_tsv()
		if "segment" not in df_meta_data.columns:
			df_meta_data["segment"] = ""

		# Strict baseline normalisation for metadata fields
		df_meta_data["primary_accession"] = df_meta_data["primary_accession"].fillna("").astype(str).str.strip()
		# ...and the identity spelling, which .str.strip() never touched. accession_version
		# is not in IDENTITY_COLUMNS and must not be: it is the revision marker.
		df_meta_data = self._normalise_identity_columns(df_meta_data, "meta_data")
		df_meta_data["segment"] = df_meta_data["segment"].fillna("").astype(str).str.strip().map(self._normalize_segment_value)

		# Build reference segment mapping from reference_tsv to backfill missing metadata segments
		ref_to_seg = {}
		if self.reference_tsv:
			try:
				ref_df = load_reference_file_table(self.reference_tsv)
				if "primary_accession" in ref_df.columns and "segment" in ref_df.columns:
					for _, row in ref_df.iterrows():
						acc = str(row["primary_accession"]).strip()
						# The same normaliser meta_data.segment went through above.
						# This was an inline digit-scrape, which disagreed with it in
						# three ways: '01' stayed '01' and so matched nothing, '1.0'
						# became segment 10, and a gene name like 'PB1' became segment
						# 1 (PB1 is segment 2). A mismatch here is silent - the value
						# is only used to backfill, and a wrong backfill then drops the
						# row from valid_pairs.
						normalized_seg = self._normalize_segment_value(row["segment"])
						if acc and normalized_seg:
							ref_to_seg[acc] = normalized_seg
			except Exception as e:
				print(f"[warn] Failed to load segments from reference_tsv: {e}")

		if ref_to_seg:
			missing_seg_mask = df_meta_data["segment"] == ""
			if missing_seg_mask.any():
				mapped_segs = df_meta_data.loc[missing_seg_mask, "primary_accession"].map(ref_to_seg).fillna("")
				df_meta_data.loc[missing_seg_mask, "segment"] = mapped_segs

		# Build absolute reference mapping dictionary based on metadata
		acc_to_seg_map = dict(zip(df_meta_data["primary_accession"], df_meta_data["segment"]))

		if "exclusion_status" not in df_meta_data.columns:
			df_meta_data["exclusion_status"] = ""
		if "exclusion_criteria" not in df_meta_data.columns:
			df_meta_data["exclusion_criteria"] = ""
		if "exclusion" in df_meta_data.columns:
			raw_exclusion = df_meta_data["exclusion"].fillna("").astype(str).str.strip()
			existing_criteria = df_meta_data["exclusion_criteria"].fillna("").astype(str).str.strip()
			df_meta_data["exclusion_criteria"] = existing_criteria.mask(existing_criteria == "", raw_exclusion)
		
		acc_type_col = "accession_type" if "accession_type" in df_meta_data.columns else None
		if acc_type_col:
			acc_type_norm = df_meta_data[acc_type_col].fillna("").str.strip().str.lower()
			is_ref_or_master = acc_type_norm.isin(["reference", "master"])
		else:
			is_ref_or_master = pd.Series(False, index=df_meta_data.index)

		if filtered_ids:
			print(f"[CreateSqliteDB] Processing {len(filtered_ids)} alignment filtering rules...")
			acc_series = df_meta_data["primary_accession"]
			seg_series = df_meta_data["segment"]
			composite_ids = acc_series + "_" + seg_series
			
			filtered_mask = (acc_series.isin(filtered_ids) | composite_ids.isin(filtered_ids)) & (~is_ref_or_master)
			
			if filtered_mask.any():
				filtered_acc = df_meta_data.loc[filtered_mask, "primary_accession"]
				filtered_seg = df_meta_data.loc[filtered_mask, "segment"]
				
				def get_reason(acc, seg):
					return str(filtered_details.get(f"{acc}_{seg}", filtered_details.get(acc, "alignment_filtering"))).strip()
				
				reason_map = pd.Series([get_reason(a, s) for a, s in zip(filtered_acc, filtered_seg)], index=filtered_acc.index)
				reason_map = reason_map.mask(reason_map == "", "alignment_filtering")
				
				existing_criteria = df_meta_data.loc[filtered_mask, "exclusion_criteria"].fillna("").astype(str).str.strip()
				df_meta_data.loc[filtered_mask, "exclusion_criteria"] = existing_criteria.where(
					existing_criteria != "",
					reason_map,
				)
				df_meta_data.loc[filtered_mask, "exclusion_status"] = "1"
				print(f"[CreateSqliteDB] Marked {int(filtered_mask.sum())} filtered segments as excluded in meta_data")

		is_ref_or_master = is_ref_or_master.reindex(df_meta_data.index, fill_value=False)
		exclusion_reason = pd.Series("", index=df_meta_data.index, dtype=str)
		for col in ["exclusion", "exclusion_criteria"]:
			if col in df_meta_data.columns:
				col_values = df_meta_data[col].fillna("").astype(str).str.strip()
				exclusion_reason = exclusion_reason.where(exclusion_reason != "", col_values)

		status_mask = pd.Series(False, index=df_meta_data.index)
		if "exclusion_status" in df_meta_data.columns:
			status_values = df_meta_data["exclusion_status"].fillna("").astype(str).str.strip()
			status_mask = ~status_values.str.lower().isin(["", "0", "false", "no", "na", "none", "nan"])
			fallback_reason = "metadata_exclusion"
			exclusion_reason = exclusion_reason.mask((exclusion_reason == "") & status_mask, fallback_reason)
		df_meta_data.loc[status_mask, "exclusion_status"] = "1"
		df_meta_data.loc[(exclusion_reason != "") & (~status_mask), "exclusion_status"] = "1"
		df_meta_data.loc[df_meta_data["exclusion_criteria"].fillna("").astype(str).str.strip() == "", "exclusion_criteria"] = exclusion_reason

		exclusion_mask = exclusion_reason.astype(str).str.strip() != ""
		incoming_query_accessions = set()
		passed_qc_query_accessions = set()
		if "primary_accession" in df_meta_data.columns:
			incoming_query_accessions = set(df_meta_data.loc[~is_ref_or_master, "primary_accession"])
			incoming_query_accessions.discard("")
			passed_qc_query_accessions = set(df_meta_data.loc[(~is_ref_or_master) & (~exclusion_mask), "primary_accession"])
			passed_qc_query_accessions.discard("")
		excluded_rows = df_meta_data[exclusion_mask]
		if not excluded_rows.empty:
			print(f"[CreateSqliteDB] Found {len(excluded_rows)} segment rows with exclusions in meta_data")
			for idx, row in excluded_rows.iterrows():
				acc = str(row.get("primary_accession", "")).strip()
				reason = str(exclusion_reason.loc[idx]).strip()
				if acc:
					excluded_records.append({"primary_accession": acc, "reason": reason})

		# Load alignment data with strict string type casting
		df_aln = self._read_tsv_required(self.pad_aln, [], "pad_aln", dtype=str)
		df_aln = self._normalize_alignment_columns(df_aln, "pad_aln")
		df_aln["primary_accession"] = df_aln["primary_accession"].fillna("").astype(str).str.strip()
		df_aln["sequence_id"] = df_aln["sequence_id"].fillna("").astype(str).str.strip()
		# Before valid_pairs is built from meta_data: a versioned accession here is
		# simply absent from it, and the alignment row is filtered away in silence.
		df_aln = self._normalise_identity_columns(df_aln, "sequence_alignment")

		# Ingest and normalise/backfill segment values safely within df_aln
		if "segment" not in df_aln.columns:
			df_aln["segment"] = df_aln["primary_accession"].map(acc_to_seg_map).fillna("")
		else:
			df_aln["segment"] = df_aln["segment"].fillna("").astype(str).str.strip()
			df_aln["segment"] = df_aln["segment"].replace("", float("NaN")).fillna(df_aln["primary_accession"].map(acc_to_seg_map)).fillna("")
		df_aln["segment"] = df_aln["segment"].map(self._normalize_segment_value)

		df_meta_data = self._add_cluster_column(df_meta_data)
		df_meta_data = self._add_reference_columns(df_meta_data, df_aln)
		
		# Establish reference sets using pure python string primitives
		valid_df = df_meta_data.loc[(~exclusion_mask) | is_ref_or_master, ["primary_accession", "segment"]].dropna()
		valid_pairs = set(zip(
			valid_df["primary_accession"].astype(str).str.strip(),
			valid_df["segment"].astype(str).str.strip()
		))

		acc_to_seg_map = dict(zip(
			df_meta_data["primary_accession"].astype(str).str.strip(),
			df_meta_data["segment"].astype(str).str.strip()
		))

		df_features = self._read_tsv_required(self.features, [], "features", dtype=str)
		feat_acc_col = "accession" if "accession" in df_features.columns else "primary_accession"
		if feat_acc_col in df_features.columns:
			df_features[feat_acc_col] = df_features[feat_acc_col].fillna("").astype(str).str.strip()
		# master_ref_accession and reference_accession are normalised here too: they
		# are joined against meta_data by VerifyMutations, not by this file, so a
		# version leak in either survives the valid_pairs filter and only shows up
		# much later as a feature map with no reference.
		df_features = self._normalise_identity_columns(df_features, "features")

		# Ingest and normalise/backfill segment values safely within df_features
		if "segment" not in df_features.columns:
			df_features["segment"] = df_features[feat_acc_col].map(acc_to_seg_map).fillna("") if feat_acc_col in df_features.columns else ""
		else:
			df_features["segment"] = df_features["segment"].fillna("").astype(str).str.strip()
			if feat_acc_col in df_features.columns:
				df_features["segment"] = df_features["segment"].replace("", float("NaN")).fillna(df_features[feat_acc_col].map(acc_to_seg_map)).fillna("")
		df_features["segment"] = df_features["segment"].map(self._normalize_segment_value)

		# Apply high-performance structural list comprehension filters against valid_pairs
	
		if valid_pairs:
			if feat_acc_col in df_features.columns:
				df_features = df_features[[t in valid_pairs for t in zip(df_features[feat_acc_col], df_features["segment"])]]
			
			# FIX: Filter using the clean primary_accession column to match valid_pairs structural layout
			df_aln = df_aln[[t in valid_pairs for t in zip(df_aln["primary_accession"], df_aln["segment"])]]

		if self.gene_info is not None:
			# keep_na_strings: influenza's neuraminidase gene is named "NA". Without
			# this the shipped IAV database stores it as ('Neuraminidase', None, '',
			# 'whole_genome') - a gene with no name.
			df_gene = self._read_tsv_required(self.gene_info, [], "gene_info", keep_na_strings=True)
		else:
			unique_products = []
			if "product" in df_features.columns:
				unique_products = df_features["product"].dropna().unique()
			elif "feature" in df_features.columns:
				unique_products = df_features["feature"].dropna().unique()
			unique_products = [str(p).strip() for p in unique_products if str(p).strip()]
			gene_records = []
			for p in unique_products:
				gene_records.append({
					"description": p,
					"display_name": p,
					"name": p,
					"parent_name": "whole_genome"
				})
			gene_records.append({
				"description": "Whole genome",
				"display_name": "Whole genome",
				"name": "whole_genome",
				"parent_name": "NULL"
			})
			df_gene = pd.DataFrame(gene_records)
			
		df_m49_country = self._read_csv_required(self.m49_countries, ["m49_code"], "m49_countries", dtype={"m49_code": str})
		df_m49_interm = self._read_csv_required(self.m49_interm_region, [], "m49_interm_region")
		df_m49_region = self._read_csv_required(self.m49_regions, [], "m49_regions")
		df_m49_sub_region = self._read_csv_required(self.m49_sub_regions, [], "m49_sub_regions")
		df_proj_setting = self._read_tsv_required(self.proj_settings, [], "proj_settings")
		
		df_insertions = self._read_tsv_required(self.insertions, [], "insertions", dtype=str)
		df_insertions = self._ensure_primary_accession(df_insertions, "insertions", aliases=["accession", "sequence_id"])
		# Filtered against valid_pairs below with the same membership idiom as
		# features, so it needs the same canonical spelling.
		df_insertions = self._normalise_identity_columns(df_insertions, "insertions")
		df_insertions["primary_accession"] = df_insertions["primary_accession"].fillna("").astype(str).str.strip()
		if "segment" not in df_insertions.columns:
			df_insertions["segment"] = df_insertions["primary_accession"].map(acc_to_seg_map).fillna("")
		else:
			df_insertions["segment"] = df_insertions["segment"].fillna("").astype(str).str.strip()
			df_insertions["segment"] = df_insertions["segment"].replace("", float("NaN")).fillna(df_insertions["primary_accession"].map(acc_to_seg_map)).fillna("")
		df_insertions["segment"] = df_insertions["segment"].map(self._normalize_segment_value)
			
		if valid_pairs and "primary_accession" in df_insertions.columns:
			df_insertions = df_insertions[[t in valid_pairs for t in zip(df_insertions["primary_accession"], df_insertions["segment"])]]
		
		df_host_taxa = self._read_tsv_required(self.host_taxa_file, [], "host_taxa_file", dtype=str)
		host_meta_col = next((c for c in ["host_taxa_id", "taxonomy_id", "host_tax_id"] if c in df_meta_data.columns), None)
		if host_meta_col and not df_host_taxa.empty:
			valid_taxa = set(df_meta_data[host_meta_col].dropna().astype(str).str.strip())
			taxa_col = next((c for c in ["taxa_id", "taxonomy_id", "host_taxa_id", "tax_id", "id"] if c in df_host_taxa.columns), None)
			if taxa_col:
				df_host_taxa = df_host_taxa[df_host_taxa[taxa_col].astype(str).str.strip().isin(valid_taxa)]
		
		df_fasta_sequences = self.load_fasta()
		df_fasta_sequences["header"] = df_fasta_sequences["header"].fillna("").astype(str).str.strip()
		# record.id keeps whatever the FASTA carried, so this is the one identity
		# column the pipeline does not control the spelling of at all.
		df_fasta_sequences = self._normalise_identity_columns(df_fasta_sequences, "sequences")
		if "segment" not in df_fasta_sequences.columns:
			df_fasta_sequences["segment"] = df_fasta_sequences["header"].map(acc_to_seg_map).fillna("")
		else:
			df_fasta_sequences["segment"] = df_fasta_sequences["segment"].fillna("").astype(str).str.strip()
			df_fasta_sequences["segment"] = df_fasta_sequences["segment"].replace("", float("NaN")).fillna(df_fasta_sequences["header"].map(acc_to_seg_map)).fillna("")
		df_fasta_sequences["segment"] = df_fasta_sequences["segment"].map(self._normalize_segment_value)

		if valid_pairs and "header" in df_fasta_sequences.columns:
			df_fasta_sequences = df_fasta_sequences[[t in valid_pairs for t in zip(df_fasta_sequences["header"], df_fasta_sequences["segment"])]]

		db_path = self._db_path()
		os.makedirs(os.path.dirname(db_path), exist_ok=True)
		self._prepare_update_target_db(db_path)
		conn = sqlite3.connect(db_path)
		cursor = conn.cursor()
		cursor.execute("PRAGMA foreign_keys = ON;")
		self._check_update_history(conn)
		cursor.execute("DROP TABLE IF EXISTS excluded_accessions;")
		cursor.execute("CREATE TABLE IF NOT EXISTS trees (name TEXT, source TEXT, segment_key TEXT, segment TEXT, newick TEXT, created_at TEXT);")
		cursor.execute("CREATE TABLE IF NOT EXISTS info (creation_type TEXT, date TEXT);")
		cursor.execute("CREATE TABLE IF NOT EXISTS update_exclusions (batch_id TEXT, table_name TEXT, key TEXT, reason TEXT, date TEXT);")
		cursor.execute("CREATE TABLE IF NOT EXISTS update_batches (batch_id TEXT PRIMARY KEY, started_at TEXT, finished_at TEXT, update_db TEXT, mode TEXT);")
		cursor.execute("CREATE TABLE IF NOT EXISTS update_table_deltas (batch_id TEXT, table_name TEXT, before_count INTEGER, after_count INTEGER, delta INTEGER);")

		now_str = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

		cursor.execute("INSERT OR REPLACE INTO update_batches (batch_id, started_at, finished_at, update_db, mode) VALUES (?, ?, ?, ?, ?)", (self.batch_id, now_str, None, self.update_db if self.update else None, "update" if self.update else "create"))

		if self._should_force_unsegmented_segment_one(conn, [df_meta_data, df_features, df_aln, df_insertions]):
			df_meta_data = self._force_segment_one_df(df_meta_data, create_if_missing=True)
			df_features = self._force_segment_one_df(df_features, create_if_missing=True)
			df_aln = self._force_segment_one_df(df_aln, create_if_missing=True)
			df_insertions = self._force_segment_one_df(df_insertions, create_if_missing=True)
			df_fasta_sequences = self._force_segment_one_df(df_fasta_sequences, create_if_missing=True)

		tables_for_delta = ["meta_data", "features", "sequence_alignment", "genes", "sequences", "insertions", "host_taxa"]
		before_counts = {t: self._table_row_count(conn, t) for t in tables_for_delta}
		before_query_count = self._count_distinct_nonexcluded_query_accessions(conn)

		self.merge_table_append_nonredundant(conn, df_meta_data, "meta_data", None, update_exclusions)
		self.merge_table_append_nonredundant(conn, df_features, "features", None, update_exclusions)
		self.merge_table_append_nonredundant(conn, df_aln, "sequence_alignment", None, update_exclusions)
		self.merge_table_append_nonredundant(conn, df_gene, "genes", None, update_exclusions)
		self.merge_table_append_nonredundant(conn, df_m49_country, "m49_country", ["m49_code"], update_exclusions)
		self.merge_table_append_nonredundant(conn, df_m49_interm, "m49_intermediate", None, update_exclusions)
		self.merge_table_append_nonredundant(conn, df_m49_region, "m49_regions", None, update_exclusions)
		self.merge_table_append_nonredundant(conn, df_m49_sub_region, "m49_sub_regions", None, update_exclusions)
		self.merge_table_append_nonredundant(conn, df_proj_setting, "project_settings", None, update_exclusions)
		self.merge_table_append_nonredundant(conn, df_fasta_sequences, "sequences", None, update_exclusions)
		self.merge_table_append_nonredundant(conn, df_insertions, "insertions", None, update_exclusions)
		self.merge_table_append_nonredundant(conn, df_host_taxa, "host_taxa", None, update_exclusions)
		
		if self._should_force_unsegmented_segment_one(conn, [df_meta_data, df_features, df_aln, df_insertions]):
			# 'sequences' belongs here: it carries a segment column and is upserted on
			# (header, segment), so a row left blank is a row no later batch can reach.
			self._backfill_segment_one_in_db(conn, ["meta_data", "features", "sequence_alignment", "insertions", "sequences"])

		df_excluded = pd.DataFrame(excluded_records, columns=["primary_accession", "reason"])
		if not df_excluded.empty:
			df_excluded = df_excluded.drop_duplicates(subset=["primary_accession"])
			print(f"[CreateSqliteDB] Marked {len(df_excluded)} excluded accessions directly on meta_data")

		if update_exclusions:
			pd.DataFrame(update_exclusions).to_sql("update_exclusions", conn, if_exists="append", index=False)
			print(f"[CreateSqliteDB] Logged {len(update_exclusions)} update-mode duplicate keys into update_exclusions")

		accession_to_segment = {}
		if "primary_accession" in df_meta_data.columns and "segment" in df_meta_data.columns:
			for _, row in df_meta_data[["primary_accession", "segment"]].dropna().iterrows():
				acc = str(row["primary_accession"]).strip()
				seg = str(row["segment"]).strip()
				if acc and seg and acc not in accession_to_segment:
					accession_to_segment[acc] = seg

		is_single_segment = self._should_force_unsegmented_segment_one(conn, [df_meta_data, df_features, df_aln, df_insertions])

		tree_records = []
		for name, source, tree_path in [("veryfasttree", "veryfasttree", self.tree_file), ("iqtree", "iqtree", self.iqtree_file), ("usher", "usher", self.usher_tree)]:
			newick = self._read_tree_file(tree_path)
			if newick:
				if is_single_segment and source == "usher":
					continue
				tree_records.append({"name": name, "source": source, "segment_key": None, "segment": None, "newick": newick})

		for entry in self._load_tree_manifest(self.tree_manifest):
			newick = self._read_tree_file(entry["path"])
			if not newick:
				continue
			seg_key = entry.get("segment_key")
			seg_num = accession_to_segment.get(seg_key) if seg_key else None
			if seg_num is None:
				seg_num = self._segment_from_key(seg_key)
			name = entry.get("name") or f"{entry['source']}_{seg_key or 'tree'}"
			if is_single_segment and entry.get("source") == "usher":
				continue
			tree_records.append({"name": name, "source": entry["source"], "segment_key": seg_key, "segment": seg_num, "newick": newick})

		if is_single_segment:
			usher_newick = None
			if self.usher_tree:
				usher_newick = self._read_tree_file(self.usher_tree)
			if not usher_newick and self.tree_manifest:
				for entry in self._load_tree_manifest(self.tree_manifest):
					if entry.get("source") == "usher":
						usher_newick = self._read_tree_file(entry["path"])
						if usher_newick:
							break
			if usher_newick:
				tree_records.append({
					"name": "Usher_tree_full_segment_1",
					"source": "usher",
					"segment_key": None,
					"segment": "1",
					"newick": usher_newick
				})

		if tree_records:
			df_tree = pd.DataFrame(tree_records)
			df_tree["created_at"] = now_str
			if not self.update:
				df_tree.to_sql("trees", conn, if_exists="replace", index=False)
			else:
				if self._table_exists(conn, "trees"):
					if is_single_segment:
						conn.execute("DELETE FROM trees WHERE source = 'usher'")
					for _, tr in df_tree.iterrows():
						if is_single_segment and tr.get("source") == "usher":
							continue
						conn.execute(
							"DELETE FROM trees WHERE COALESCE(source, '')=? AND COALESCE(name, '')=? AND COALESCE(segment_key, '')=? AND COALESCE(segment, '')=?",
							(
								str(tr.get("source") or "").strip(),
								str(tr.get("name") or "").strip(),
								str(tr.get("segment_key") or "").strip(),
								str(tr.get("segment") or "").strip(),
							),
						)
				df_tree.to_sql("trees", conn, if_exists="append", index=False)

		creation_type = self._resolve_creation_type()
		# Two columns, appended - the info table's shape is unchanged, so anything
		# reading (creation_type, date) keeps working; only the value is now honest.
		pd.DataFrame([{"creation_type": creation_type, "date": now_str}]).to_sql("info", conn, if_exists="append", index=False)

		after_counts = {t: self._table_row_count(conn, t) for t in tables_for_delta}
		after_query_count = self._count_distinct_nonexcluded_query_accessions(conn)
		deltas = []
		for table_name in tables_for_delta:
			before = before_counts.get(table_name, 0)
			after = after_counts.get(table_name, 0)
			deltas.append({"batch_id": self.batch_id, "table_name": table_name, "before_count": int(before), "after_count": int(after), "delta": int(after - before)})
		pd.DataFrame(deltas).to_sql("update_table_deltas", conn, if_exists="append", index=False)
		cursor.execute("UPDATE update_batches SET finished_at=? WHERE batch_id=?", (datetime.now().strftime("%Y-%m-%d %H:%M:%S"), self.batch_id))

		reason_counts = pd.Series(dtype=int)
		detailed_counts = pd.Series(dtype=int)
		if not df_excluded.empty:
			normalized_reason = df_excluded['reason'].map(self._normalize_summary_reason)
			normalized_reason = normalized_reason.mask(normalized_reason == "", "alignment_filtering")
			reason_counts = normalized_reason.value_counts()
			detailed_reason = normalized_reason
			detailed_reason = detailed_reason.mask(detailed_reason == "", "alignment_filtering")
			detailed_counts = detailed_reason.value_counts().head(10)

		delta_by_table = {row["table_name"]: int(row["delta"]) for row in deltas}
		summary_lines = []
		summary_lines.append("\n" + "="*80)
		summary_lines.append("[CreateSqliteDB] DB Construction Summary:")
		summary_lines.append("="*80)
		summary_lines.append(f"Total sequences loaded      : {len(df_meta_data)}")
		summary_lines.append(f"  (includes references/masters and query sequences)")
		summary_lines.append(f"Query sequences passing QC  : {len(passed_qc_query_accessions)}")
		summary_lines.append(f"Sequences marked as excluded: {len(df_excluded)}")
		summary_lines.append(f"  (references/masters + failed QC)")

		segment_summary = None
		if "segment" in df_meta_data.columns:
			query_rows = df_meta_data.loc[~is_ref_or_master].copy()
			query_exclusion_mask = exclusion_mask.loc[query_rows.index]
			segment_series = query_rows["segment"].fillna("").astype(str).str.strip()
			nonblank_segment_values = {seg for seg in segment_series.unique() if seg}
			if len(nonblank_segment_values) > 1:
				normalized_segment = segment_series.mask(segment_series == "", "NA")
				segment_frame = pd.DataFrame({
					"segment": normalized_segment,
					"included": (~query_exclusion_mask).astype(int),
					"excluded": query_exclusion_mask.astype(int),
				})
				segment_summary = (
					segment_frame.groupby("segment", as_index=False)[["included", "excluded"]]
					.sum()
				)
				segment_summary["_segment_num"] = pd.to_numeric(segment_summary["segment"], errors="coerce")
				segment_summary = (
					segment_summary.sort_values(["_segment_num", "segment"], na_position="last")
					.drop(columns=["_segment_num"])
				)

		if segment_summary is not None and not segment_summary.empty:
			summary_lines.append("-" * 80)
			summary_lines.append("[CreateSqliteDB] Query Sequences by Segment (passed/failed QC):")
			summary_lines.append(f"{'Segment':<45} | {'Passed QC':<10} | {'Failed QC':<10}")
			summary_lines.append("-" * 80)
			for _, row in segment_summary.iterrows():
				summary_lines.append(
					f"{str(row['segment']):<45} | {int(row['included']):<10} | {int(row['excluded']):<10}"
				)
		if not reason_counts.empty:
			summary_lines.append("-" * 80)
			summary_lines.append(f"{'Category':<45} | {'Count':<10}")
			summary_lines.append("-" * 80)
			for reason, count in reason_counts.items():
				summary_lines.append(f"{reason:<45} | {count:<10}")
			summary_lines.append("\n[CreateSqliteDB] Detailed Exclusion Reasons (Top 10):")
			for reason, count in detailed_counts.items():
				summary_lines.append(f"  - {count:<4} : {reason}")

		if self.update:
			summary_lines.append("\n[CreateSqliteDB] Update Summary:")
			summary_lines.append(f"Sequences added to DB       : {delta_by_table.get('meta_data', 0)}")
			summary_lines.append(f"Input query sequences       : {len(incoming_query_accessions)}")
			summary_lines.append(f"Input query sequences passed QC : {len(passed_qc_query_accessions)}")
			summary_lines.append(f"Input query sequences failed QC : {len(incoming_query_accessions - passed_qc_query_accessions)}")
			summary_lines.append(f"Total passed-QC queries in DB   : {after_query_count} (delta: {after_query_count - before_query_count})")

		summary_lines.append("="*80 + "\n")
		summary_str = "\n".join(summary_lines)
		print(summary_str)

		summary_path = os.path.join(self.base_dir, self.output_dir, "db_summary.txt")
		with open(summary_path, "w", encoding="utf-8") as sf:
			sf.write(summary_str)

		conn.commit()
		conn.close()


def process(args):
	db_creator = CreateSqliteDB(
		args.meta_data,
		args.features,
		args.pad_aln,
		args.gene_info,
		args.m49_countries,
		args.m49_interm_region,
		args.m49_regions,
		args.m49_sub_regions,
		args.proj_settings,
		args.fasta_sequences,
		args.insertion_file,
		args.host_taxa_file,
		args.base_dir,
		args.output_dir,
		args.db_name,
		args.db_status,
		args.tree_file,
		args.iqtree_file,
		args.usher_tree,
		args.cluster_tsv,
		args.cluster_min_seq_id,
		args.filtered_ids,
		args.filtered_details,
		args.tree_manifest,
		args.reference_tsv,
		clade_assignments=args.clade_assignments,
		update=args.update,
		update_db=args.update_db,
		batch_id=args.batch_id,
		is_segmented=(None if args.is_segmented is None else args.is_segmented == 'Y'),
	)
	db_creator.create_db()


if __name__ == "__main__":
	parser = ArgumentParser(description='Creating sqlite DB')
	parser.add_argument('-m', '--meta_data', help='Meta data table', default="tmp/GenBank-matrix/gB_matrix_raw.tsv")
	parser.add_argument('-b', '--base_dir', help='Base directory', default="tmp")
	parser.add_argument('-o', '--output_dir', help='tmp directory where the database is stored', default="SqliteDB")
	parser.add_argument('-rf', '--features', help='Features table', default="tmp/Tables/features.tsv")
	parser.add_argument('-p', '--pad_aln', help='Padded alignment file', default="tmp/Tables/sequence_alignment.tsv")
	parser.add_argument('-g', '--gene_info', help='Gene table', default=None)
	parser.add_argument('--is_segmented', choices=['Y', 'N'], default=None,
		help="Whether the virus is segmented. When given it is authoritative for the 'fill blank segments with 1' decision, which otherwise has to be inferred from the data and cannot tell a non-segmented build from a segmented one that happens to hold only segment 1.")
	parser.add_argument('-mc', '--m49_countries', help='M49 countries', default="assets/m49_country.csv")
	parser.add_argument('-mir', '--m49_interm_region', help='M49 intermediate regions', default="assets/m49_intermediate_region.csv")
	parser.add_argument('-mr', '--m49_regions', help='M49 regions', default="assets/m49_region.csv")
	parser.add_argument('-msr', '--m49_sub_regions', help='M49 sub-regions', default="assets/m49_sub_region.csv")
	parser.add_argument('-s', '--proj_settings', help='Project settings', default="tmp/Software_info/software_info.tsv")
	parser.add_argument('-fa', '--fasta_sequences', help='Fasta sequences', default="tmp/GenBank-matrix/sequences.fa")
	parser.add_argument('-i', '--insertion_file', help='Nextalign insertion file', default="tmp/Tables/insertions.tsv")
	parser.add_argument('-ht', '--host_taxa_file', help='Host Taxanomy file', default="tmp/HostTaxa/Host_taxa.tsv")
	parser.add_argument('-d', '--db_name', help='Name of the Sqlite database', default="gdb")
	parser.add_argument('-ds', '--db_status', default=None,
		help='Database status recorded in info.creation_type: "new db" or "last updated". '
			 'Defaults to the run mode ("last updated" with --update, "new db" without), and '
			 'the run mode wins if the two contradict each other.')
	parser.add_argument('-t', '--tree_file', help='VeryFastTree Newick file', default=None)
	parser.add_argument('-it', '--iqtree_file', help='IQ-TREE Newick file', default=None)
	parser.add_argument('-ut', '--usher_tree', help='UShER output Newick file', default=None)
	parser.add_argument('--tree_manifest', help='TSV manifest with columns: source, name, segment_key, path', default=None)
	parser.add_argument('-ct', '--cluster_tsv', help='MMseqs clustering TSV (rep\tmember)', default=None)
	parser.add_argument('-ci', '--cluster_min_seq_id', help='MMseqs min sequence identity used for clustering', default=None)
	parser.add_argument('-fi', '--filtered_ids', help='File with filtered sequence IDs (one per line) to exclude from DB', default=None)
	parser.add_argument('-fd', '--filtered_details', help='TSV with filtered sequence details (seq_name, reference, error, warnings)', default=None)
	parser.add_argument('--update', action='store_true', help='Append/update into an existing DB instead of replacing tables')
	parser.add_argument('--update_db', default=None, help='Path to existing DB file for update mode')
	parser.add_argument('--batch_id', default=None, help='Optional update batch identifier for audit logging')
	parser.add_argument('--reference_tsv', help='Optional reference TSV with columns: primary_accession, status,segment,genotype,subtype, to help resolve segment info for tree records', default=None)
	parser.add_argument('--clade_assignments', help='Optional TSV (primary_accession, genotype, subtype) from CladeAssignment.py, used when no tree is available', default=None)
	args = parser.parse_args()
	if args.update_db:
		args.update_db = normpath(args.update_db)
	try:
		process(args)
	except Exception as exc:
		print(f"ERROR: {exc}", file=sys.stderr)
		sys.exit(2)