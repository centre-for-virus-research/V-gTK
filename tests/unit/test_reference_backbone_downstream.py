"""Downstream consequences of a backbone that keeps reference insertion columns.

Once the alignment is wider than the master and the master row carries gaps:

* MMseqs completeness must be measured against the master's length, not width;
* UShER's faToVcf reference is pinned to the master when it can be;
* reference insertions already held as columns must not reach the insertions
  table a second time;
* both pipeline entry points must hand the same backbone to every step.
"""

from pathlib import Path

import MMseqsClustering
from GenerateTables import GenerateTables
from UsherPlacement import UsherPlacement

ROOT = Path(__file__).resolve().parents[2]


def _write(path: Path, text: str) -> str:
	path.parent.mkdir(parents=True, exist_ok=True)
	path.write_text(text, encoding="utf-8")
	return str(path)


# ---------------------------------------------------------------------------
# MMseqs completeness
# ---------------------------------------------------------------------------

def test_reference_length_is_the_masters_ungapped_length(tmp_path):
	# Versioned on one side, bare on the other: both normalise to NC_000001.
	aln = _write(tmp_path / "aln.fasta", ">Q1\nACGTACGT--\n>NC_000001.1\nACGT--ACGT\n")
	assert MMseqsClustering.reference_length_from_alignment(aln, ["NC_000001"]) == 8
	assert MMseqsClustering.reference_length_from_alignment(aln, ["ABSENT"]) is None
	assert MMseqsClustering.reference_length_from_alignment(aln, []) is None


def test_completeness_against_master_length_keeps_complete_genomes_in_step_one(tmp_path):
	# Width 20, master 10 nt: a complete 10 nt genome is 50% of the width but 100%
	# of the genome.
	aln = _write(tmp_path / "aln.fasta",
				 ">M\nACGTACGTAC----------\n>Q1\n----------ACGTACGTAC\n>FRAG\nACG-----------------\n")
	complete, remainder = tmp_path / "c.seq", tmp_path / "r.seq"

	by_width = MMseqsClustering.split_by_completeness(aln, 0.9, str(complete), str(remainder))
	assert by_width == (0, 3)

	by_master = MMseqsClustering.split_by_completeness(aln, 0.9, str(complete), str(remainder),
													   reference_length=10)
	assert by_master == (2, 1)
	assert ">FRAG" in remainder.read_text()


# ---------------------------------------------------------------------------
# UShER reference
# ---------------------------------------------------------------------------

def _usher(tmp_path, reference_ids):
	msa = _write(tmp_path / "msa.fasta", ">R1\nACGT\n>NC_000001\nACGA\n>Q\nACGG\n")
	return UsherPlacement(padded_aln=msa, output_dir=str(tmp_path / "out"), reference_ids=reference_ids), msa


def test_usher_reference_is_the_master_when_it_is_a_representative(tmp_path):
	processor, msa = _usher(tmp_path, ["NC_000001"])
	reps = _write(tmp_path / "reps.fasta", ">R1\nACGT\n>NC_000001\nACGA\n")
	assert processor.resolve_reference_id(cluster_rep=reps, alignment_fasta=msa) == "NC_000001"


def test_usher_reference_falls_back_when_the_master_is_not_a_representative(tmp_path):
	processor, msa = _usher(tmp_path, ["NC_000001"])
	reps = _write(tmp_path / "reps.fasta", ">R1\nACGT\n")
	assert processor.resolve_reference_id(cluster_rep=reps, alignment_fasta=msa) == "R1"


def test_usher_reference_without_masters_keeps_the_first_representative(tmp_path):
	processor, msa = _usher(tmp_path, [])
	reps = _write(tmp_path / "reps.fasta", ">R1\nACGT\n>NC_000001\nACGA\n")
	assert processor.resolve_reference_id(cluster_rep=reps, alignment_fasta=msa) == "R1"


def test_usher_reference_in_update_mode_prefers_the_master_in_the_alignment(tmp_path):
	# The reference list may carry the versioned spelling; the alignment is bare.
	processor, msa = _usher(tmp_path, ["NC_000001.1"])
	assert processor.resolve_reference_id(alignment_fasta=msa) == "NC_000001"


# ---------------------------------------------------------------------------
# Insertions table
# ---------------------------------------------------------------------------

def _nextalign_tree(root: Path):
	_write(root / "query_aln" / "R1" / "R1.insertions.csv",
		   "seqName,insertions,aaInsertions\nQ1,100:ACG,\nQ2,,\n")
	_write(root / "reference_aln" / "M" / "M.insertions.csv",
		   "seqName,insertions,aaInsertions\nR1,50:TTTTTT,\n")


def _insertions(tmp_path, reference_insertions):
	nextalign = tmp_path / "Nextalign"
	_nextalign_tree(nextalign)
	tables = GenerateTables("gb.tsv", str(tmp_path), "Tables", "hits.tsv", [], "host.tsv",
							str(nextalign), "x@example.com", None, reference_insertions)
	tables.create_insertion_table({})
	lines = (tmp_path / "Tables" / "insertions.tsv").read_text().splitlines()
	return lines[0], sorted(lines[1:])


def test_insertions_without_a_built_backbone_are_unchanged(tmp_path):
	header, rows = _insertions(tmp_path, None)
	assert header == "primary_accession\treference\tinsertion\tsegment"
	assert rows == ["Q1\tR1\t100:ACG\t", "R1\tM\t50:TTTTTT\t"]


def test_insertions_with_a_built_backbone_use_its_dropped_bases_for_references(tmp_path):
	dropped = _write(tmp_path / "ref_set_aligned" / "dropped_insertions.tsv",
					 "primary_accession\treference\tinsertion\tsegment\nR1\tM\t812:AAACCC\t\n")
	_, rows = _insertions(tmp_path, dropped)
	assert rows == ["Q1\tR1\t100:ACG\t", "R1\tM\t812:AAACCC\t"]


# ---------------------------------------------------------------------------
# Workflow wiring (both entry points)
# ---------------------------------------------------------------------------

def test_nextflow_builds_the_backbone_only_for_fresh_builds_without_a_curated_one():
	text = (ROOT / "vgtk-init.nf").read_text(encoding="utf-8")
	assert "process BUILD_REFERENCE_ALIGNMENT" in text
	assert ("def BUILD_REF_ALIGNMENT = !UPDATE_MODE && ref_backbone_dir == 'UNSET' "
			"&& params.build_ref_alignment.toString().toBoolean()") in text
	assert "BUILD_REFERENCE_ALIGNMENT(BLAST_ALIGNMENT.out.ref_seqs_fasta," in text
	assert '"build_ref_alignment", "ref_aln_min_insertion_support"' in text
	assert "check_cmd mafft --version" in text


def test_nextflow_hands_one_backbone_to_filtering_padding_and_tables():
	text = (ROOT / "vgtk-init.nf").read_text(encoding="utf-8")
	# COLLECT_FILTERED_SEQUENCES and PAD_ALIGNMENT disagreeing about the backbone is
	# what deleted 633,988 influenza sequences (MISSING_H3_report.md).
	assert "effective_ref_list,\n                               ref_backbone_ch)" in text
	assert "ref_backbone_ch,\n                  COLLECT_FILTERED_SEQUENCES.out.filtered_ids)" in text
	assert "params.ref_list,\n                    ref_backbone_ch)" in text
	assert "--reference_insertions !{ref_set_aligned_dir}/dropped_insertions.tsv" in text
	assert text.count('--ref_list "!{params.ref_list}"') == 2  # MMseqs and UShER


def test_config_declares_the_backbone_params():
	text = (ROOT / "nextflow.config").read_text(encoding="utf-8")
	assert "build_ref_alignment = true" in text
	assert "ref_aln_min_insertion_support = 2" in text


def test_bash_pipeline_mirrors_the_backbone_wiring():
	text = (ROOT / "vgtk-rabv.sh").read_text(encoding="utf-8")
	assert 'python "${SCRIPTS}/BuildReferenceAlignment.py"' in text
	assert '--precomputed_ref_dir "${REF_BACKBONE_DIR}" --ref_list "${REF_LIST}"' in text
	assert '--precomputed_ref_dir "${REF_BACKBONE_DIR:-UNSET}"' in text
	assert '--reference_insertions "${REF_BACKBONE_DIR}/dropped_insertions.tsv"' in text
	assert text.count('--ref_list "${REF_LIST}" \\') >= 3  # MMseqs and both UShER calls
	assert "_check_dep mafft" in text
