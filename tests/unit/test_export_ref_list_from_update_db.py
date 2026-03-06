from pathlib import Path

from ExportRefListFromUpdateDb import export_ref_list


def test_export_ref_list_from_update_db_writes_reference_rows(tmp_path: Path, basic_update_db: Path):
	out = tmp_path / "ref_list.tsv"
	export_ref_list(str(basic_update_db), str(out))

	text = out.read_text(encoding="utf-8")
	assert "REF1\tmaster\t1" in text
	assert "REF2\treference\t2" in text
	assert "Q_OLD" not in text