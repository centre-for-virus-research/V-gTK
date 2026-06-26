import csv
from pathlib import Path

from NextalignAlignment import NextalignAlignment


def _write_matrix(path: Path):
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(["gi_number", "exclusion_status", "exclusion_criteria"])
        writer.writerow(["ACC_OK", "0", ""])
        writer.writerow(["ACC_FAIL", "0", ""])


def _read_matrix(path: Path):
    with path.open("r", newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def test_update_gb_matrix_marks_failed_accessions(tmp_path: Path):
    matrix = tmp_path / "gB_matrix_raw.tsv"
    with matrix.open("w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(["gi_number", "accession_type", "exclusion_status", "exclusion_criteria"])
        writer.writerow(["ACC_OK", "query", "0", ""])
        writer.writerow(["ACC_FAIL", "query", "0", ""])

    aln_dir = tmp_path / "Nextalign" / "query_aln" / "REF1"
    aln_dir.mkdir(parents=True, exist_ok=True)
    (aln_dir / "REF1.errors.csv").write_text(
        'seqName,errors,warnings\nACC_FAIL,"In sequence #2 ''ACC_FAIL'': Unable to align, not enough matches.",\n',
        encoding="utf-8",
    )

    processor = NextalignAlignment(
        gb_matrix=str(matrix),
        query_dir=str(tmp_path / "query"),
        ref_dir=str(tmp_path / "ref"),
        ref_fa_file=str(tmp_path / "ref.fa"),
        master_seq_dir=str(tmp_path / "master"),
        tmp_dir=str(tmp_path),
        master_ref="MASTER1",
        nextalign_dir="Nextalign",
        reference_alignment=None,
    )

    processor.update_gb_matrix([str(tmp_path / "Nextalign" / "query_aln")], str(matrix))

    rows = {r["gi_number"]: r for r in _read_matrix(matrix)}
    assert rows["ACC_FAIL"]["exclusion_status"] == "1"
    assert rows["ACC_FAIL"]["exclusion_criteria"] == "In sequence #2 ACC_FAIL: Unable to align, not enough matches."
    assert rows["ACC_OK"]["exclusion_status"] == "0"


def test_update_gb_matrix_ignores_reference_alignment_errors(tmp_path: Path):
    matrix = tmp_path / "gB_matrix_raw.tsv"
    with matrix.open("w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(["gi_number", "accession_type", "exclusion_status", "exclusion_criteria"])
        writer.writerow(["REF_FAIL", "reference", "0", ""])
        writer.writerow(["Q_FAIL", "query", "0", ""])

    query_dir = tmp_path / "Nextalign" / "query_aln" / "REF1"
    query_dir.mkdir(parents=True, exist_ok=True)
    (query_dir / "REF1.errors.csv").write_text(
        "seqName,errors,warnings\nQ_FAIL,query failure,\n",
        encoding="utf-8",
    )

    ref_dir = tmp_path / "Nextalign" / "reference_aln" / "MASTER1"
    ref_dir.mkdir(parents=True, exist_ok=True)
    (ref_dir / "MASTER1.errors.csv").write_text(
        "seqName,errors,warnings\nREF_FAIL,reference failure,\n",
        encoding="utf-8",
    )

    processor = NextalignAlignment(
        gb_matrix=str(matrix),
        query_dir=str(tmp_path / "query"),
        ref_dir=str(tmp_path / "ref"),
        ref_fa_file=str(tmp_path / "ref.fa"),
        master_seq_dir=str(tmp_path / "master"),
        tmp_dir=str(tmp_path),
        master_ref="MASTER1",
        nextalign_dir="Nextalign",
        reference_alignment=None,
    )

    processor.update_gb_matrix(
        [str(tmp_path / "Nextalign" / "query_aln"), str(tmp_path / "Nextalign" / "reference_aln")],
        str(matrix),
    )

    rows = {r["gi_number"]: r for r in _read_matrix(matrix)}
    assert rows["Q_FAIL"]["exclusion_status"] == "1"
    assert rows["Q_FAIL"]["exclusion_criteria"] == "query failure"
    assert rows["REF_FAIL"]["exclusion_status"] == "0"
    assert rows["REF_FAIL"]["exclusion_criteria"] == ""


def test_process_reference_alignment_mode_runs_query_only(tmp_path: Path, monkeypatch):
    matrix = tmp_path / "gB_matrix_raw.tsv"
    _write_matrix(matrix)

    query_dir = tmp_path / "query"
    ref_dir = tmp_path / "ref"
    query_dir.mkdir(parents=True, exist_ok=True)
    ref_dir.mkdir(parents=True, exist_ok=True)
    (query_dir / "REF1.fa").write_text(">Q\nATGC\n", encoding="utf-8")
    (ref_dir / "REF1.fa").write_text(">R\nATGC\n", encoding="utf-8")

    processor = NextalignAlignment(
        gb_matrix=str(matrix),
        query_dir=str(query_dir),
        ref_dir=str(ref_dir),
        ref_fa_file=str(tmp_path / "ref.fa"),
        master_seq_dir=str(tmp_path / "master"),
        tmp_dir=str(tmp_path),
        master_ref="MASTER1",
        nextalign_dir="Nextalign",
        reference_alignment="provided_alignment.fa",
    )

    calls = {"query": 0, "master": 0, "update": 0}

    def fake_query(*_args, **_kwargs):
        calls["query"] += 1

    def fake_master(*_args, **_kwargs):
        calls["master"] += 1

    def fake_update(*_args, **_kwargs):
        calls["update"] += 1

    monkeypatch.setattr(processor, "nextalign_query", fake_query)
    monkeypatch.setattr(processor, "nextalign_master", fake_master)
    monkeypatch.setattr(processor, "update_gb_matrix", fake_update)

    processor.process()

    assert calls == {"query": 1, "master": 0, "update": 1}


def test_process_non_reference_alignment_runs_master_path(tmp_path: Path, monkeypatch):
    matrix = tmp_path / "gB_matrix_raw.tsv"
    _write_matrix(matrix)

    query_dir = tmp_path / "query"
    ref_dir = tmp_path / "ref"
    master_dir = tmp_path / "master"
    query_dir.mkdir(parents=True, exist_ok=True)
    ref_dir.mkdir(parents=True, exist_ok=True)
    master_dir.mkdir(parents=True, exist_ok=True)

    (query_dir / "REF1.fa").write_text(">Q\nATGC\n", encoding="utf-8")
    (ref_dir / "REF1.fa").write_text(">R\nATGC\n", encoding="utf-8")
    (master_dir / "MASTER1.fasta").write_text(">MASTER1\nATGC\n", encoding="utf-8")

    processor = NextalignAlignment(
        gb_matrix=str(matrix),
        query_dir=str(query_dir),
        ref_dir=str(ref_dir),
        ref_fa_file=str(tmp_path / "ref.fa"),
        master_seq_dir=str(master_dir),
        tmp_dir=str(tmp_path),
        master_ref="MASTER1",
        nextalign_dir="Nextalign",
        reference_alignment=None,
    )

    calls = {"query": 0, "master": 0, "update": 0, "dedup": 0}

    def fake_query(*_args, **_kwargs):
        calls["query"] += 1

    def fake_master(query_acc_path, ref_acc_path, query_aln_op):
        calls["master"] += 1
        master = Path(ref_acc_path).stem
        out_dir = Path(query_aln_op) / master
        out_dir.mkdir(parents=True, exist_ok=True)
        (out_dir / f"{master}.aligned.fasta").write_text(">MASTER1\nATGC\n", encoding="utf-8")

    class _FakeRemoveRedundantSequence:
        def __init__(self, _input_seq, _output_seq):
            pass

        def remove_redundant_fasta(self):
            calls["dedup"] += 1

    def fake_update(*_args, **_kwargs):
        calls["update"] += 1

    monkeypatch.setattr(processor, "nextalign_query", fake_query)
    monkeypatch.setattr(processor, "nextalign_master", fake_master)
    monkeypatch.setattr(processor, "update_gb_matrix", fake_update)
    monkeypatch.setattr("NextalignAlignment.RemoveRedundantSequence", _FakeRemoveRedundantSequence)
    processor.process()

    assert calls["query"] == 1
    assert calls["master"] == 1
    assert calls["dedup"] == 1
    assert calls["update"] == 1


def test_validate_alignment_integrity_checks_all_sequences(tmp_path: Path):
    processor = NextalignAlignment(
        gb_matrix=str(tmp_path / "gB_matrix.tsv"),
        query_dir=str(tmp_path / "query"),
        ref_dir=str(tmp_path / "ref"),
        ref_fa_file=str(tmp_path / "ref.fa"),
        master_seq_dir=str(tmp_path / "master"),
        tmp_dir=str(tmp_path),
        master_ref="MASTER1",
        nextalign_dir="Nextalign",
        reference_alignment=None,
    )

    # 10 columns total. Conserved region is final 30% (columns 7 to 9, 0-indexed).
    # Master sequence:
    # "A T G C G T A A C G" (indices 0 to 9)
    # Downstream conserved (indices 7, 8, 9): A, C, G
    
    # CASE 1: Both queries align perfectly in conserved region
    clean_fasta = tmp_path / "clean.fasta"
    clean_fasta.write_text(
        ">MASTER\nATGCGTAACG\n"
        ">QUERY1\nATGCGTAACG\n"  # matches ACG -> 100% identity
        ">QUERY2\nATGCGTAACG\n", # matches ACG -> 100% identity
        encoding="utf-8",
    )
    passed, min_id = processor._validate_alignment_integrity(str(clean_fasta))
    assert passed is True
    assert min_id == 1.0

    # CASE 2: One query aligns perfectly, second query is misaligned/mismatched in conserved region
    dirty_fasta = tmp_path / "dirty.fasta"
    dirty_fasta.write_text(
        ">MASTER\nATGCGTAACG\n"
        ">QUERY1\nATGCGTAACG\n"  # matches ACG -> 100% identity
        ">QUERY2\nATGCGTGGGG\n", # matches GGG instead of AACG -> 33% identity in conserved region
        encoding="utf-8",
    )
    passed, min_id = processor._validate_alignment_integrity(str(dirty_fasta))
    assert passed is False
    assert abs(min_id - 0.333333) < 1e-5

