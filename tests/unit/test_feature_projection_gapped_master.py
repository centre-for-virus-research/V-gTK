"""Feature projection on a master row that carries gaps.

With a guide or built backbone, a row can start or end inside an insertion
column (the master is gapped there). CalcAlignmentCord clamps a feature to the
master bases the row covers, so its ``cds_*`` columns end on a master base while
``aln_*`` ends on the row's last base. The validator must compare like with like.

Found on real runs: 13 HCV reference rows (HCV_test, built backbone) and 3
influenza NA rows (segmented_xml_test, curated backbone) were flagged although
their projection was correct.
"""

import sqlite3

import BuildReferenceAlignment as bra
from ValidateDbTree import validate_feature_projection_integrity


def _db(tmp_path, master_row, rows):
	conn = sqlite3.connect(str(tmp_path / "t.db"))
	conn.execute("CREATE TABLE sequence_alignment (primary_accession TEXT, alignment_name TEXT, alignment TEXT, segment TEXT)")
	conn.execute("CREATE TABLE features (accession TEXT, master_ref_accession TEXT, product TEXT, "
				 "aln_start TEXT, aln_end TEXT, cds_start TEXT, cds_end TEXT)")
	conn.execute("INSERT INTO sequence_alignment VALUES ('M', 'M', ?, '1')", (master_row,))
	for row in rows:
		conn.execute("INSERT INTO features VALUES (?, ?, ?, ?, ?, ?, ?)", row)
	conn.commit()
	return conn


# Master: 12 bases with a 2-column gap at columns 11-12, gene over columns 1-10.
MASTER = "ACGTACGTAC--GT"


def test_row_ending_in_an_insertion_column_is_not_flagged(tmp_path):
	conn = _db(tmp_path, MASTER, [
		("M", "M", "gene", "1", "14", "1", "10"),
		# Q's last base is column 12 (master gap); its last covered master base is column 10.
		("Q", "M", "gene", "1", "12", "1", "10"),
	])
	result = validate_feature_projection_integrity(conn)
	assert result["ok"], result["examples"]


def test_a_genuinely_wrong_cds_end_is_still_flagged(tmp_path):
	conn = _db(tmp_path, MASTER, [
		("M", "M", "gene", "1", "14", "1", "10"),
		("Q", "M", "gene", "1", "12", "1", "8"),
	])
	result = validate_feature_projection_integrity(conn)
	assert not result["ok"]
	assert result["examples"][0]["accession"] == "Q"


def test_without_a_stored_master_row_the_old_rule_applies(tmp_path):
	conn = sqlite3.connect(str(tmp_path / "t.db"))
	conn.execute("CREATE TABLE features (accession TEXT, master_ref_accession TEXT, product TEXT, "
				 "aln_start TEXT, aln_end TEXT, cds_start TEXT, cds_end TEXT)")
	conn.execute("INSERT INTO features VALUES ('M', 'M', 'gene', '1', '14', '1', '10')")
	conn.execute("INSERT INTO features VALUES ('Q', 'M', 'gene', '1', '12', '1', '10')")
	conn.commit()
	assert validate_feature_projection_integrity(conn)["ok"]
	conn.execute("UPDATE features SET cds_end='9' WHERE accession='Q'")
	conn.commit()
	assert not validate_feature_projection_integrity(conn)["ok"]


# ---------------------------------------------------------------------------
# The builder side: a partial genome's last bases stay beside its codons
# ---------------------------------------------------------------------------

def test_trailing_partial_codon_goes_into_the_next_empty_codon_slot():
	block = "ATGAAA---------"
	left, new_block, right = bra.tuck_partial_codons("", block, "GC")
	assert (left, right) == ("", "")
	assert new_block == "ATGAAAGC-------"


def test_leading_partial_codon_goes_into_the_previous_empty_codon_slot():
	block = "------ATGAAA"
	left, new_block, right = bra.tuck_partial_codons("TG", block, "")
	assert (left, right) == ("", "")
	assert new_block == "----TGATGAAA"


def test_partial_codons_stay_in_edge_columns_when_the_row_reaches_the_block_edge():
	block = "ATGAAACCC"
	assert bra.tuck_partial_codons("A", block, "GC") == ("A", block, "GC")


def test_tucking_never_changes_the_row_bases():
	block = "---ATG---AAA------"
	left, new_block, right = bra.tuck_partial_codons("C", block, "TT")
	assert (left + new_block + right).replace("-", "") == "CATGAAATT"
	assert len(new_block) == len(block)
