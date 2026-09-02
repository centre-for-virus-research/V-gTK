"""Unit and adversarial coverage for scripts/NextalignAlignment.py.

The happy paths live here - a well-formed errors.csv, a well-formed gB matrix
and two ideal alignments - next to the boundaries this stage is actually handed
by NCBI, by the filesystem and by nextalign itself: paths carrying shell
metacharacters or spaces, errors.csv files whose header is not the one the
parser expects, ragged TSV rows, alignments whose records are not all the same
length, and 3' ends full of ``N``.

Most of the boundary tests were written against a defect that was first
reproduced on the unfixed module. They now pass, and each keeps its docstring
saying what it caught, so a regression re-breaks a test with an explanation
attached rather than an anonymous assertion.

The module no longer goes through a shell, so the harness intercepts
``subprocess.run`` and inspects the **argv list**. Two tests deliberately go
further and execute a stub ``nextalign`` on PATH, because the difference between
"a list" and "a string a shell re-parses" cannot be proven by inspecting the
list alone.

Everything here is synthetic and writes only into tmp_path.
"""

import csv
import os
import subprocess
import sys
from pathlib import Path

import pytest

import NextalignAlignment as na_module
from NextalignAlignment import NextalignAlignment, NextalignNotAvailable


MATRIX_HEADER = ["gi_number", "accession_type", "exclusion_status", "exclusion_criteria"]


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------


def make(tmp_path: Path, **kwargs) -> NextalignAlignment:
    defaults = dict(
        gb_matrix=str(tmp_path / "gB_matrix_raw.tsv"),
        query_dir=str(tmp_path / "query"),
        ref_dir=str(tmp_path / "ref"),
        ref_fa_file=str(tmp_path / "ref.fa"),
        master_seq_dir=str(tmp_path / "master"),
        tmp_dir=str(tmp_path),
        master_ref="MASTER1",
        nextalign_dir="Nextalign",
        reference_alignment=None,
    )
    defaults.update(kwargs)
    return NextalignAlignment(**defaults)


def write_matrix(path: Path, rows, header=None):
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(header or MATRIX_HEADER)
        writer.writerows(rows)


def read_matrix(path: Path):
    with path.open("r", newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_errors_csv(aln_root: Path, group: str, text: str) -> str:
    """Lay out a ``<aln_root>/query_aln/<group>/<group>.errors.csv`` tree."""
    group_dir = aln_root / "query_aln" / group
    group_dir.mkdir(parents=True, exist_ok=True)
    (group_dir / f"{group}.errors.csv").write_text(text, encoding="utf-8")
    return str(aln_root / "query_aln")


def write_fasta(path: Path, records):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        for seq_id, sequence in records:
            handle.write(f">{seq_id}\n{sequence}\n")


class Recorder(list):
    """Captured subprocess.run calls, newest-friendly accessors."""

    def argv(self, index=0):
        return self[index]["argv"]

    def kwargs(self, index=0):
        return self[index]["kwargs"]

    def arg_after(self, flag, index=0):
        argv = self.argv(index)
        return argv[argv.index(flag) + 1]


def capture_commands(monkeypatch, return_code=0, side_effect=None):
    """Intercept the argv list handed to subprocess.run.

    Deliberately records the LIST, not a joined string: the whole point of the
    fix is that the arguments never become a string a shell can re-parse.
    """
    seen = Recorder()

    class _Completed:
        def __init__(self, code):
            self.returncode = code

    def fake_run(command, **kwargs):
        seen.append({"argv": list(command), "kwargs": kwargs})
        if side_effect is not None:
            side_effect(list(command))
        return _Completed(return_code)

    monkeypatch.setattr(na_module.subprocess, "run", fake_run)
    return seen


def stub_nextalign(tmp_path: Path, monkeypatch, body="exit 0"):
    """Put a controllable real `nextalign` on PATH.

    Used only where the test must prove something about actual process
    execution rather than about the argv list.
    """
    bin_dir = tmp_path / "stub_bin"
    bin_dir.mkdir(parents=True, exist_ok=True)
    stub = bin_dir / "nextalign"
    stub.write_text(f"#!/bin/sh\n{body}\n", encoding="utf-8")
    stub.chmod(0o755)
    monkeypatch.setenv("PATH", str(bin_dir) + os.pathsep + os.environ["PATH"])
    return stub


# ===========================================================================
# 1. Command execution: argv is never flattened into a shell string
# ===========================================================================
#
# Both entry points used to build a proper argv list and then throw the safety
# away with `os.system(" ".join(command))`. Everything in this section was a
# consequence of that one line.


def test_command_is_never_handed_to_a_shell(tmp_path: Path, monkeypatch):
    """subprocess.run receives a list, and shell is never enabled."""
    seen = capture_commands(monkeypatch)

    make(tmp_path).nextalign_query(
        str(tmp_path / "ACC1.fa"), str(tmp_path / "ref.fa"), str(tmp_path / "out")
    )

    assert isinstance(seen.argv(), list)
    assert seen.kwargs().get("shell", False) is False


def test_paths_containing_spaces_survive_as_single_arguments(tmp_path: Path, monkeypatch):
    """A data directory called 'RABV sequences' used to be word-split by
    /bin/sh, so nextalign ran against a path that did not exist."""
    seen = capture_commands(monkeypatch)
    query = tmp_path / "RABV sequences" / "GROUP 1.fa"
    write_fasta(query, [("Q", "ATGC")])
    reference = tmp_path / "RABV sequences" / "REF 1.fa"
    write_fasta(reference, [("R", "ATGC")])

    make(tmp_path).nextalign_query(str(query), str(reference), str(tmp_path / "out"))

    argv = seen.argv()
    assert argv[-1] == str(query)
    assert seen.arg_after("--input-ref") == str(reference)


def test_shell_metacharacters_in_a_path_are_not_executed(tmp_path: Path, monkeypatch):
    """Proven end to end against a real process, not just the argv list: a
    directory named 'grp; touch FILE #' used to create FILE."""
    stub_nextalign(tmp_path, monkeypatch)
    sentinel = tmp_path / "INJECTED"

    hostile_dir = tmp_path / f"grp; touch {sentinel} #"
    hostile_dir.mkdir(parents=True)
    query = hostile_dir / "ACC1.fa"
    write_fasta(query, [("Q", "ATGC")])

    make(tmp_path).nextalign_query(str(query), str(tmp_path / "ref.fa"), str(tmp_path / "out"))

    assert not sentinel.exists(), "a path fragment was executed as a shell command"


def test_shell_chaining_cannot_mask_a_failed_alignment(tmp_path: Path, monkeypatch):
    """os.system returned the status of the LAST command in the string, so a
    path containing '; true' made a failed alignment report success."""
    stub_nextalign(tmp_path, monkeypatch, body="exit 1")

    hostile_dir = tmp_path / "grp; true #"
    hostile_dir.mkdir(parents=True)
    query = hostile_dir / "ACC1.fa"
    write_fasta(query, [("Q", "ATGC")])

    aligned = make(tmp_path).nextalign_query(
        str(query), str(tmp_path / "ref.fa"), str(tmp_path / "out")
    )

    assert aligned is False


def test_query_failure_reports_the_real_exit_code(tmp_path: Path, monkeypatch, capsys):
    """os.system returned a wait status, so exit 3 was printed as 'code 768'."""
    stub_nextalign(tmp_path, monkeypatch, body="exit 3")

    aligned = make(tmp_path).nextalign_query(
        str(tmp_path / "ACC1.fa"), str(tmp_path / "ref.fa"), str(tmp_path / "out")
    )

    assert aligned is False
    assert "code 3" in capsys.readouterr().out


def test_empty_accession_stays_an_argument(tmp_path: Path, monkeypatch):
    """An empty --output-basename used to vanish from the joined string, so
    nextalign silently took the next flag as the basename value."""
    seen = capture_commands(monkeypatch)
    monkeypatch.setattr(NextalignAlignment, "path_to_basename", staticmethod(lambda _p: ""))

    make(tmp_path).nextalign_query(
        str(tmp_path / "weird"), str(tmp_path / "ref.fa"), str(tmp_path / "out")
    )

    assert seen.arg_after("--output-basename") == ""


def test_query_command_carries_the_baseline_profile(tmp_path: Path, monkeypatch):
    """Queries use profile 0 and never relax."""
    seen = capture_commands(monkeypatch)
    processor = make(tmp_path, max_threads=14)

    processor.nextalign_query(
        str(tmp_path / "ACC1.fa"), str(tmp_path / "ref.fa"), str(tmp_path / "out")
    )

    assert len(seen) == 1
    baseline = processor.relaxation_profiles[0]
    assert seen.arg_after("--min-seeds") == str(baseline["min_seeds"])
    assert seen.arg_after("--seed-spacing") == str(baseline["seed_spacing"])
    assert seen.arg_after("--min-match-rate") == str(baseline["min_match_rate"])
    assert seen.arg_after("--jobs") == "14"


def test_missing_binary_raises_immediately_for_a_query(tmp_path: Path, monkeypatch):
    """A missing nextalign is an environment fault, not sequence data."""
    def boom(_command, **_kwargs):
        raise FileNotFoundError(2, "No such file or directory", "nextalign")

    monkeypatch.setattr(na_module.subprocess, "run", boom)

    with pytest.raises(NextalignNotAvailable):
        make(tmp_path).nextalign_query(
            str(tmp_path / "ACC1.fa"), str(tmp_path / "ref.fa"), str(tmp_path / "out")
        )


def test_missing_binary_is_not_retried_or_blamed_on_divergence(tmp_path: Path, monkeypatch, capsys):
    """The shell turned a missing binary into status 127, which the relaxation
    loop retried five times and then reported as 'extreme sequence
    divergence' - pointing the operator at the data, not the conda env."""
    calls = []

    def boom(command, **_kwargs):
        calls.append(list(command))
        raise FileNotFoundError(2, "No such file or directory", "nextalign")

    monkeypatch.setattr(na_module.subprocess, "run", boom)

    with pytest.raises(NextalignNotAvailable):
        make(tmp_path).nextalign_master(
            str(tmp_path / "ref.fa"), str(tmp_path / "MASTER1.fa"), str(tmp_path / "out")
        )

    assert len(calls) == 1, "a missing binary must not be retried through every profile"
    assert "divergence" not in capsys.readouterr().out


# ===========================================================================
# 2. _read_nextalign_errors: the positional fallback used to be unreachable
# ===========================================================================
#
# `handle.seek(0)` sat OUTSIDE the `with open(...)` block, so every errors.csv
# without exactly the seqName+errors header raised ValueError on a closed file.


def test_headerless_errors_csv_falls_back_to_positional_columns(tmp_path: Path):
    errors = tmp_path / "REF1.errors.csv"
    errors.write_text("ACC_FAIL,not enough matches\n", encoding="utf-8")

    failed = NextalignAlignment._read_nextalign_errors(str(errors), "REF1")

    assert failed == {"ACC_FAIL": ["not enough matches"]}


def test_unexpected_errors_csv_header_does_not_crash(tmp_path: Path):
    """A nextalign build naming the column 'error' rather than 'errors' used to
    take the whole pipeline down instead of degrading to the positional reader."""
    errors = tmp_path / "REF1.errors.csv"
    errors.write_text("seqName,error\nACC_FAIL,not enough matches\n", encoding="utf-8")

    failed = NextalignAlignment._read_nextalign_errors(str(errors), "REF1")

    assert failed == {"ACC_FAIL": ["not enough matches"]}


def test_unrecognised_header_row_is_not_read_as_a_sequence(tmp_path: Path):
    """The positional fallback must not turn the header itself into a failure."""
    errors = tmp_path / "REF1.errors.csv"
    errors.write_text("seqName,error\nACC_FAIL,boom\n", encoding="utf-8")

    assert "seqName" not in NextalignAlignment._read_nextalign_errors(str(errors), "REF1")


def test_empty_errors_csv_reports_no_failures(tmp_path: Path):
    """A zero-byte errors.csv - a killed or out-of-disk nextalign."""
    errors = tmp_path / "REF1.errors.csv"
    errors.write_text("", encoding="utf-8")

    assert NextalignAlignment._read_nextalign_errors(str(errors), "REF1") == {}


def test_header_only_errors_csv_reports_no_failures(tmp_path: Path):
    """The ordinary 'nothing failed' output - nextalign still writes the header."""
    errors = tmp_path / "REF1.errors.csv"
    errors.write_text("seqName,errors,warnings\n", encoding="utf-8")

    assert NextalignAlignment._read_nextalign_errors(str(errors), "REF1") == {}


def test_repeated_errors_for_one_accession_accumulate(tmp_path: Path):
    errors = tmp_path / "REF1.errors.csv"
    errors.write_text(
        "seqName,errors,warnings\nACC_FAIL,first,\nACC_FAIL,second,\n", encoding="utf-8"
    )

    failed = NextalignAlignment._read_nextalign_errors(str(errors), "REF1")

    assert failed == {"ACC_FAIL": ["first", "second"]}


def test_reference_own_errors_are_filtered_out(tmp_path: Path):
    """--include-reference puts the reference in the alignment; its own error
    rows must not be attributed to a query."""
    errors = tmp_path / "REF1.errors.csv"
    errors.write_text(
        "seqName,errors,warnings\nREF1,reference noise,\nACC_FAIL,real failure,\n",
        encoding="utf-8",
    )

    failed = NextalignAlignment._read_nextalign_errors(str(errors), "REF1")

    assert failed == {"ACC_FAIL": ["real failure"]}


def test_reference_filtering_also_applies_to_the_positional_fallback(tmp_path: Path):
    errors = tmp_path / "REF1.errors.csv"
    errors.write_text("REF1,reference noise\nACC_FAIL,real failure\n", encoding="utf-8")

    failed = NextalignAlignment._read_nextalign_errors(str(errors), "REF1")

    assert failed == {"ACC_FAIL": ["real failure"]}


def test_blank_error_text_is_not_a_failure(tmp_path: Path):
    errors = tmp_path / "REF1.errors.csv"
    errors.write_text("seqName,errors,warnings\nACC_OK,,some warning\n", encoding="utf-8")

    assert NextalignAlignment._read_nextalign_errors(str(errors), "REF1") == {}


def test_missing_error_column_values_do_not_crash(tmp_path: Path):
    """A short row gives DictReader a None value, which str(None) used to turn
    into the literal string 'None' and record as a failure reason."""
    errors = tmp_path / "REF1.errors.csv"
    errors.write_text("seqName,errors,warnings\nACC_OK\n", encoding="utf-8")

    assert NextalignAlignment._read_nextalign_errors(str(errors), "REF1") == {}


# ===========================================================================
# 3. update_gb_matrix: the pipeline's central table is rewritten atomically
# ===========================================================================


def test_ragged_row_does_not_destroy_the_matrix(tmp_path: Path):
    """DATA LOSS regression. The matrix used to be truncated by open(..., 'w')
    before the rows were serialised: a row with more fields than the header
    raised part-way through writerows, leaving only the rows written before the
    raise. Reproduced on the unfixed code as 3 rows in, 1 row on disk."""
    matrix = tmp_path / "gB_matrix_raw.tsv"
    write_matrix(
        matrix,
        [
            ["ACC_OK", "query", "0", ""],
            ["ACC_RAGGED", "query", "0", "", "leftover"],
            ["ACC_LAST", "query", "0", ""],
        ],
    )
    before = matrix.read_bytes()

    with pytest.raises(ValueError, match="more fields than"):
        make(tmp_path).update_gb_matrix([], str(matrix))

    assert matrix.read_bytes() == before, "the matrix was modified despite the failure"


def test_write_failure_leaves_the_previous_matrix_intact(tmp_path: Path, monkeypatch):
    """The general case of the above: any exception during serialisation must
    not cost the caller their matrix."""
    matrix = tmp_path / "gB_matrix_raw.tsv"
    write_matrix(matrix, [["ACC_OK", "query", "0", ""]])
    before = matrix.read_bytes()

    class ExplodingWriter:
        def __init__(self, *_a, **_k):
            pass

        def writeheader(self):
            pass

        def writerows(self, _rows):
            raise OSError("No space left on device")

    monkeypatch.setattr(na_module.csv, "DictWriter", ExplodingWriter)

    with pytest.raises(OSError):
        make(tmp_path).update_gb_matrix([], str(matrix))

    assert matrix.read_bytes() == before


def test_failed_write_leaves_no_temp_file_behind(tmp_path: Path, monkeypatch):
    matrix = tmp_path / "gB_matrix_raw.tsv"
    write_matrix(matrix, [["ACC_OK", "query", "0", ""]])

    class ExplodingWriter:
        def __init__(self, *_a, **_k):
            pass

        def writeheader(self):
            raise OSError("boom")

    monkeypatch.setattr(na_module.csv, "DictWriter", ExplodingWriter)

    with pytest.raises(OSError):
        make(tmp_path).update_gb_matrix([], str(matrix))

    assert not (tmp_path / "gB_matrix_raw.tsv.tmp").exists()
    assert sorted(p.name for p in tmp_path.iterdir()) == ["gB_matrix_raw.tsv"]


def test_missing_exclusion_columns_default_to_zero(tmp_path: Path):
    """reader.fieldnames used to be appended to IN PLACE while the DictReader
    was still parsing, so the new keys were materialised on every row with
    restval None; row.get('exclusion_status', '0') then returned None - the key
    exists - and every row was written blank instead of '0'."""
    matrix = tmp_path / "gB_matrix_raw.tsv"
    write_matrix(matrix, [["ACC_OK", "query"]], header=["gi_number", "accession_type"])

    make(tmp_path).update_gb_matrix([], str(matrix))

    row = read_matrix(matrix)[0]
    assert row["exclusion_status"] == "0"
    assert row["exclusion_criteria"] == ""


def test_blank_exclusion_status_is_normalised_to_zero(tmp_path: Path):
    matrix = tmp_path / "gB_matrix_raw.tsv"
    write_matrix(matrix, [["ACC_OK", "query", "", ""]])

    make(tmp_path).update_gb_matrix([], str(matrix))

    assert read_matrix(matrix)[0]["exclusion_status"] == "0"


def test_empty_matrix_file_reports_a_usable_error(tmp_path: Path):
    """Used to raise TypeError: argument of type 'NoneType' is not iterable,
    which says nothing about which file is malformed."""
    matrix = tmp_path / "gB_matrix_raw.tsv"
    matrix.write_text("", encoding="utf-8")

    with pytest.raises(ValueError, match="no header row"):
        make(tmp_path).update_gb_matrix([], str(matrix))


def test_update_is_idempotent(tmp_path: Path):
    """Criteria used to be appended unconditionally, so a resumed nextflow run
    produced 'not enough matches; not enough matches' and grew the field on
    every pass."""
    matrix = tmp_path / "gB_matrix_raw.tsv"
    write_matrix(matrix, [["ACC_FAIL", "query", "0", ""]])
    aln = write_errors_csv(
        tmp_path / "Nextalign", "REF1",
        "seqName,errors,warnings\nACC_FAIL,not enough matches,\n",
    )
    processor = make(tmp_path)

    processor.update_gb_matrix([aln], str(matrix))
    first = read_matrix(matrix)[0]["exclusion_criteria"]
    processor.update_gb_matrix([aln], str(matrix))
    second = read_matrix(matrix)[0]["exclusion_criteria"]

    assert second == first == "not enough matches"


def test_distinct_criteria_still_accumulate(tmp_path: Path):
    """Idempotency must not collapse genuinely different reasons."""
    matrix = tmp_path / "gB_matrix_raw.tsv"
    write_matrix(matrix, [["ACC_FAIL", "query", "0", ""]])
    aln = write_errors_csv(
        tmp_path / "Nextalign", "REF1",
        "seqName,errors,warnings\nACC_FAIL,not enough matches,\nACC_FAIL,too divergent,\n",
    )

    make(tmp_path).update_gb_matrix([aln], str(matrix))

    assert read_matrix(matrix)[0]["exclusion_criteria"] == "not enough matches; too divergent"


@pytest.mark.parametrize(
    "raw",
    ['"left\tright"', '"line one\nline two"', '"carriage\rreturn"'],
    ids=["tab", "newline", "carriage-return"],
)
def test_whitespace_in_error_text_keeps_the_tsv_rectangular(tmp_path: Path, raw):
    """The matrix is read downstream both with csv.DictReader and with a bare
    line.strip().split('\\t') (ValidateMatrix.py:107), so quoting an embedded
    tab satisfies the first reader while silently breaking the second."""
    matrix = tmp_path / "gB_matrix_raw.tsv"
    write_matrix(matrix, [["ACC_FAIL", "query", "0", ""]])
    aln = write_errors_csv(
        tmp_path / "Nextalign", "REF1", f"seqName,errors,warnings\nACC_FAIL,{raw},\n",
    )

    make(tmp_path).update_gb_matrix([aln], str(matrix))

    lines = matrix.read_text(encoding="utf-8").splitlines()
    assert len(lines) == 2, "the criteria text introduced a line break"
    assert len(lines[1].split("\t")) == len(lines[0].split("\t"))


def test_utf8_matrix_survives_a_c_locale(tmp_path: Path):
    """Both open() calls used to omit encoding=, so under LC_ALL=C - the default
    in many schedulers and containers - a strain name with an accent raised
    UnicodeDecodeError."""
    matrix = tmp_path / "gB_matrix_raw.tsv"
    matrix.write_text(
        "gi_number\taccession_type\texclusion_status\texclusion_criteria\tstrain\n"
        "ACC1\tquery\t0\t\tRABV/Côte-d'Ivoire/2021\n",
        encoding="utf-8",
    )
    scripts_dir = str(Path(__file__).resolve().parents[2] / "scripts")
    program = (
        "import sys\n"
        f"sys.path.insert(0, {scripts_dir!r})\n"
        "import locale\n"
        "if locale.getpreferredencoding(False).lower().replace('-', '') == 'utf8':\n"
        "    raise SystemExit(88)\n"
        "from NextalignAlignment import NextalignAlignment as C\n"
        f"C({str(matrix)!r}, 'q', 'r', 'f', 'ms', {str(tmp_path)!r}, 'M', 'N', None)"
        f".update_gb_matrix([], {str(matrix)!r})\n"
    )
    result = subprocess.run(
        [sys.executable, "-X", "utf8=0", "-c", program],
        env={**os.environ, "LC_ALL": "C", "LANG": "C", "PYTHONUTF8": "0",
             "PYTHONCOERCECLOCALE": "0"},
        capture_output=True, text=True,
    )
    if result.returncode == 88:
        pytest.skip("this interpreter cannot be forced out of UTF-8 mode")

    assert result.returncode == 0, result.stderr
    assert "Côte" in matrix.read_text(encoding="utf-8")


def test_master_and_exclusion_list_rows_are_never_excluded(tmp_path: Path):
    """Extends the existing 'reference' test to the other two protected types,
    and to the case/whitespace variants the .strip().lower() guard implies."""
    matrix = tmp_path / "gB_matrix_raw.tsv"
    write_matrix(
        matrix,
        [
            ["ACC_MASTER", " Master ", "0", ""],
            ["ACC_EXCL", "EXCLUSION_LIST", "0", ""],
            ["ACC_REF", "Reference", "0", ""],
            ["ACC_QUERY", "query", "0", ""],
        ],
    )
    aln = write_errors_csv(
        tmp_path / "Nextalign", "REF1",
        "seqName,errors,warnings\n"
        "ACC_MASTER,boom,\nACC_EXCL,boom,\nACC_REF,boom,\nACC_QUERY,boom,\n",
    )

    make(tmp_path).update_gb_matrix([aln], str(matrix))

    rows = {r["gi_number"]: r for r in read_matrix(matrix)}
    assert rows["ACC_MASTER"]["exclusion_status"] == "0"
    assert rows["ACC_EXCL"]["exclusion_status"] == "0"
    assert rows["ACC_REF"]["exclusion_status"] == "0"
    assert rows["ACC_QUERY"]["exclusion_status"] == "1"


def test_existing_criteria_are_preserved_and_appended(tmp_path: Path):
    matrix = tmp_path / "gB_matrix_raw.tsv"
    write_matrix(matrix, [["ACC_FAIL", "query", "1", "too short"]])
    aln = write_errors_csv(
        tmp_path / "Nextalign", "REF1",
        "seqName,errors,warnings\nACC_FAIL,not enough matches,\n",
    )

    make(tmp_path).update_gb_matrix([aln], str(matrix))

    assert read_matrix(matrix)[0]["exclusion_criteria"] == "too short; not enough matches"


def test_upstream_exclusions_are_not_reset(tmp_path: Path):
    """A sequence excluded by an earlier stage must stay excluded even though
    nextalign has nothing to say about it."""
    matrix = tmp_path / "gB_matrix_raw.tsv"
    write_matrix(matrix, [["ACC_OLD", "query", "1", "excluded upstream"]])

    make(tmp_path).update_gb_matrix([], str(matrix))

    row = read_matrix(matrix)[0]
    assert row["exclusion_status"] == "1"
    assert row["exclusion_criteria"] == "excluded upstream"


def test_unrelated_columns_are_preserved_verbatim(tmp_path: Path):
    matrix = tmp_path / "gB_matrix_raw.tsv"
    write_matrix(
        matrix,
        [["ACC1", "query", "0", "", "Rabies lyssavirus", "2021-06-01"]],
        header=MATRIX_HEADER + ["organism", "collection_date"],
    )

    make(tmp_path).update_gb_matrix([], str(matrix))

    row = read_matrix(matrix)[0]
    assert row["organism"] == "Rabies lyssavirus"
    assert row["collection_date"] == "2021-06-01"


def test_only_query_aln_directories_are_scanned(tmp_path: Path):
    """process() passes reference_aln in the list but the loop skips anything
    not named query_aln. Pinned because the signature invites the opposite
    reading, and because the guard is what keeps master failures out."""
    matrix = tmp_path / "gB_matrix_raw.tsv"
    write_matrix(matrix, [["ACC_A", "query", "0", ""], ["ACC_B", "query", "0", ""]])

    write_errors_csv(tmp_path / "Nextalign", "REF1", "seqName,errors,warnings\nACC_A,boom,\n")
    ref_group = tmp_path / "Nextalign" / "reference_aln" / "MASTER1"
    ref_group.mkdir(parents=True)
    (ref_group / "MASTER1.errors.csv").write_text(
        "seqName,errors,warnings\nACC_B,boom,\n", encoding="utf-8"
    )

    make(tmp_path).update_gb_matrix(
        [str(tmp_path / "Nextalign" / "query_aln"), str(tmp_path / "Nextalign" / "reference_aln")],
        str(matrix),
    )

    rows = {r["gi_number"]: r for r in read_matrix(matrix)}
    assert rows["ACC_A"]["exclusion_status"] == "1"
    assert rows["ACC_B"]["exclusion_status"] == "0"


def test_a_missing_alignment_directory_is_not_fatal(tmp_path: Path):
    """process() creates these, but a resumed or hand-driven run may not have."""
    matrix = tmp_path / "gB_matrix_raw.tsv"
    write_matrix(matrix, [["ACC1", "query", "0", ""]])

    make(tmp_path).update_gb_matrix([str(tmp_path / "nope" / "query_aln")], str(matrix))

    assert read_matrix(matrix)[0]["exclusion_status"] == "0"


def test_extra_failures_are_applied_alongside_errors_csv(tmp_path: Path):
    """Groups that never produced an errors.csv reach the matrix through the
    extra_failures channel."""
    matrix = tmp_path / "gB_matrix_raw.tsv"
    write_matrix(matrix, [["ACC_CSV", "query", "0", ""], ["ACC_NORUN", "query", "0", ""]])
    aln = write_errors_csv(
        tmp_path / "Nextalign", "REF1", "seqName,errors,warnings\nACC_CSV,boom,\n"
    )

    make(tmp_path).update_gb_matrix(
        [aln], str(matrix), {"ACC_NORUN": ["nextalign failed to run"]}
    )

    rows = {r["gi_number"]: r for r in read_matrix(matrix)}
    assert rows["ACC_CSV"]["exclusion_status"] == "1"
    assert rows["ACC_NORUN"]["exclusion_status"] == "1"
    assert rows["ACC_NORUN"]["exclusion_criteria"] == "nextalign failed to run"


def test_accession_matching_is_exact(tmp_path: Path):
    """Documents a coupling worth knowing about: seqName must equal gi_number
    character for character. If nextalign ever emits the full FASTA header
    (accession + description) the failure is dropped in silence rather than
    flagged - there is no diagnostic for 'error reported for an unknown
    accession'."""
    matrix = tmp_path / "gB_matrix_raw.tsv"
    write_matrix(matrix, [["ACC1", "query", "0", ""]])
    aln = write_errors_csv(
        tmp_path / "Nextalign", "REF1",
        "seqName,errors,warnings\nACC1 Rabies lyssavirus strain X,boom,\n",
    )

    make(tmp_path).update_gb_matrix([aln], str(matrix))

    assert read_matrix(matrix)[0]["exclusion_status"] == "0"

def test_update_gb_matrix_marks_failed_accessions(tmp_path: Path):
    """The end-to-end happy path, with a verbatim nextalign error line: the
    named accession is excluded and carries the reason, its neighbour is left
    alone."""
    matrix = tmp_path / "gB_matrix_raw.tsv"
    write_matrix(matrix, [["ACC_OK", "query", "0", ""], ["ACC_FAIL", "query", "0", ""]])

    aln = write_errors_csv(
        tmp_path / "Nextalign", "REF1",
        'seqName,errors,warnings\n'
        'ACC_FAIL,"In sequence #2 ACC_FAIL: Unable to align, not enough matches.",\n',
    )

    make(tmp_path).update_gb_matrix([aln], str(matrix))

    rows = {r["gi_number"]: r for r in read_matrix(matrix)}
    assert rows["ACC_FAIL"]["exclusion_status"] == "1"
    assert rows["ACC_FAIL"]["exclusion_criteria"] == (
        "In sequence #2 ACC_FAIL: Unable to align, not enough matches."
    )
    assert rows["ACC_OK"]["exclusion_status"] == "0"


def test_update_gb_matrix_ignores_reference_alignment_errors(tmp_path: Path):
    """Overlaps test_only_query_aln_directories_are_scanned above, and is kept
    because it additionally pins the reference row's exclusion_criteria as
    still empty rather than merely unset."""
    matrix = tmp_path / "gB_matrix_raw.tsv"
    write_matrix(matrix, [["REF_FAIL", "reference", "0", ""], ["Q_FAIL", "query", "0", ""]])

    query_aln = write_errors_csv(
        tmp_path / "Nextalign", "REF1", "seqName,errors,warnings\nQ_FAIL,query failure,\n"
    )
    ref_group = tmp_path / "Nextalign" / "reference_aln" / "MASTER1"
    ref_group.mkdir(parents=True, exist_ok=True)
    (ref_group / "MASTER1.errors.csv").write_text(
        "seqName,errors,warnings\nREF_FAIL,reference failure,\n", encoding="utf-8"
    )

    make(tmp_path).update_gb_matrix(
        [query_aln, str(tmp_path / "Nextalign" / "reference_aln")], str(matrix)
    )

    rows = {r["gi_number"]: r for r in read_matrix(matrix)}
    assert rows["Q_FAIL"]["exclusion_status"] == "1"
    assert rows["Q_FAIL"]["exclusion_criteria"] == "query failure"
    assert rows["REF_FAIL"]["exclusion_status"] == "0"
    assert rows["REF_FAIL"]["exclusion_criteria"] == ""


# ===========================================================================
# 4. _validate_alignment_integrity: the heuristic that can abort the pipeline
# ===========================================================================
#
# A False here makes nextalign_master relax and, after five profiles, raise
# RuntimeError and kill the run. Every false negative is a pipeline outage, not
# a dropped sequence.


def test_ragged_alignment_fails_closed_instead_of_crashing(tmp_path: Path):
    """The column loop was bounded by len(reference) but indexed the query, so
    any record shorter than the reference raised IndexError. An alignment is
    rectangular by definition; a record of a different length means nextalign
    produced something malformed, so this fails closed."""
    aln = tmp_path / "ragged.fasta"
    write_fasta(aln, [("REF", "ATGCGTAACG"), ("Q1", "ATGCG")])

    passed, identity = make(tmp_path)._validate_alignment_integrity(str(aln))

    assert passed is False
    assert identity == 0.0


def test_over_long_record_also_fails_closed(tmp_path: Path):
    aln = tmp_path / "long.fasta"
    write_fasta(aln, [("REF", "ATGCGTAACG"), ("Q1", "ATGCGTAACGTTTT")])

    assert make(tmp_path)._validate_alignment_integrity(str(aln))[0] is False


def test_ambiguity_codes_are_not_counted_as_mismatches(tmp_path: Path):
    """N was compared as a literal character, so a sequence whose called bases
    all match but whose 3' end is a run of N - routine in GenBank - scored 0%
    and, via nextalign_master, aborted the whole run."""
    aln = tmp_path / "ntail.fasta"
    write_fasta(
        aln,
        [("REF", "ATGCGTAACGTTGCACGATC"), ("Q1", "ATGCGTAACGTTGCNNNNNN")],
    )

    passed, identity = make(tmp_path)._validate_alignment_integrity(str(aln))

    assert passed is True, f"a perfectly-matching sequence with an N tail scored {identity:.2%}"


@pytest.mark.parametrize("code", list("RYSWKMBDHVN"))
def test_every_iupac_ambiguity_code_is_skipped(tmp_path: Path, code):
    """The full IUPAC set, not just N - GenBank records carry all of them."""
    aln = tmp_path / f"iupac_{code}.fasta"
    write_fasta(
        aln,
        [("REF", "ATGCGTAACGTTGCACGATC"), ("Q1", "ATGCGTAACGTTGC" + code * 6)],
    )

    assert make(tmp_path)._validate_alignment_integrity(str(aln))[0] is True


def test_ambiguity_does_not_hide_a_real_mismatch(tmp_path: Path):
    """Skipping N must not become a way for a frame-shifted sequence to pass:
    the called bases still have to agree."""
    aln = tmp_path / "mixed_n.fasta"
    write_fasta(
        aln,
        [("REF", "ATGCGTAACGTTGCACGATC"), ("Q1", "ATGCGTAACGTTGCNTTTTT")],
    )

    passed, identity = make(tmp_path)._validate_alignment_integrity(str(aln))

    # One N is skipped; the remaining five called bases are scored, and only
    # one of them agrees with the reference.
    assert passed is False
    assert identity == pytest.approx(1 / 5)


def test_no_comparable_columns_is_not_reported_as_total_mismatch(tmp_path: Path):
    """When every reference column in the window is a gap - a query extending
    past the reference's 3' end - valid_positions is 0. 'Nothing comparable'
    used to be returned as identity 0.0, which is the difference between
    continuing and aborting the run."""
    aln = tmp_path / "overhang.fasta"
    write_fasta(
        aln,
        [("REF", "ATGCGTAACGTTGCAC" + "-" * 8), ("Q1", "ATGCGTAACGTTGCAC" + "AAAAAAAA")],
    )

    passed, _identity = make(tmp_path)._validate_alignment_integrity(str(aln))

    assert passed is True


def test_unvalidatable_sequence_is_announced(tmp_path: Path, capsys):
    """Passing an unvalidatable sequence is only acceptable if it is visible."""
    aln = tmp_path / "overhang.fasta"
    write_fasta(
        aln,
        [("REF", "ATGCGTAACGTTGCAC" + "-" * 8), ("Q1", "ATGCGTAACGTTGCAC" + "AAAAAAAA")],
    )

    make(tmp_path)._validate_alignment_integrity(str(aln))

    assert "not validated" in capsys.readouterr().out


def test_master_does_not_abort_on_a_valid_alignment_with_an_n_tail(tmp_path: Path, monkeypatch):
    """The consequence at the call site: nextalign_master used to burn all five
    relaxation profiles on a valid alignment and then take the pipeline down
    because one master reference had an N tail."""
    out_root = tmp_path / "reference_aln"

    def write_alignment(_argv):
        write_fasta(
            out_root / "MASTER1" / "MASTER1.aligned.fasta",
            [("MASTER1", "ATGCGTAACGTTGCACGATC"), ("Q1", "ATGCGTAACGTTGCNNNNNN")],
        )

    seen = capture_commands(monkeypatch, side_effect=write_alignment)

    make(tmp_path).nextalign_master(
        str(tmp_path / "ref.fa"), str(tmp_path / "MASTER1.fasta"), str(out_root)
    )

    assert len(seen) == 1, "a valid alignment should not need a relaxation profile"


def test_gap_only_query_scores_zero(tmp_path: Path):
    aln = tmp_path / "allgap.fasta"
    write_fasta(aln, [("REF", "ATGCGTAACG"), ("Q1", "-" * 10)])

    passed, identity = make(tmp_path)._validate_alignment_integrity(str(aln))

    assert passed is False
    assert identity == 0.0


def test_single_record_alignment_fails_closed(tmp_path: Path):
    """Every query failed to align, so only the reference is left."""
    aln = tmp_path / "lonely.fasta"
    write_fasta(aln, [("REF", "ATGCGTAACG")])

    assert make(tmp_path)._validate_alignment_integrity(str(aln)) == (False, 0.0)


def test_missing_alignment_file_fails_closed(tmp_path: Path):
    assert make(tmp_path)._validate_alignment_integrity(str(tmp_path / "nope.fasta")) == (False, 0.0)


def test_comparison_is_case_insensitive(tmp_path: Path):
    """Soft-masked (lowercase) output must not read as 100% mismatch."""
    aln = tmp_path / "softmask.fasta"
    write_fasta(aln, [("REF", "ATGCGTAACG"), ("Q1", "atgcgtaacg")])

    passed, identity = make(tmp_path)._validate_alignment_integrity(str(aln))

    assert passed is True
    assert identity == 1.0


def test_min_identity_reports_the_worst_query(tmp_path: Path):
    aln = tmp_path / "mixed.fasta"
    write_fasta(
        aln,
        [
            ("REF", "ATGCGTAACG"),
            ("GOOD", "ATGCGTAACG"),
            ("BAD", "ATGCGTGGGG"),
            ("ALSO_GOOD", "ATGCGTAACG"),
        ],
    )

    passed, identity = make(tmp_path)._validate_alignment_integrity(str(aln))

    assert passed is False
    assert identity == pytest.approx(1 / 3)


def test_one_bad_query_does_not_hide_behind_good_ones(tmp_path: Path, capsys):
    aln = tmp_path / "mixed.fasta"
    write_fasta(
        aln, [("REF", "ATGCGTAACG"), ("GOOD", "ATGCGTAACG"), ("BAD", "ATGCGTGGGG")]
    )

    make(tmp_path)._validate_alignment_integrity(str(aln))

    assert "BAD" in capsys.readouterr().out


def test_first_record_is_taken_as_the_reference(tmp_path: Path):
    """Undocumented coupling to nextalign's --include-reference output order:
    record 0 is trusted to be the reference. If that order ever changes, every
    query is scored against another query and this stage silently stops
    validating anything."""
    aln = tmp_path / "reordered.fasta"
    write_fasta(aln, [("Q1", "ATGCGTGGGG"), ("REF", "ATGCGTAACG")])

    passed, _identity = make(tmp_path)._validate_alignment_integrity(str(aln))

    assert passed is False  # scored REF against Q1, not Q1 against REF


def test_master_stops_at_the_first_profile_that_validates(tmp_path: Path, monkeypatch):
    out_root = tmp_path / "reference_aln"

    def write_alignment(_argv):
        write_fasta(
            out_root / "MASTER1" / "MASTER1.aligned.fasta",
            [("MASTER1", "ATGCGTAACG"), ("Q1", "ATGCGTAACG")],
        )

    seen = capture_commands(monkeypatch, side_effect=write_alignment)
    processor = make(tmp_path)

    processor.nextalign_master(str(tmp_path / "ref.fa"), str(tmp_path / "MASTER1.fa"), str(out_root))

    assert len(seen) == 1
    assert seen.arg_after("--min-seeds") == str(processor.relaxation_profiles[0]["min_seeds"])


def test_master_walks_every_profile_before_aborting(tmp_path: Path, monkeypatch):
    """Succeeds but writes nothing, so validation fails on every profile."""
    seen = capture_commands(monkeypatch)
    processor = make(tmp_path)

    with pytest.raises(RuntimeError, match="Alignment Frame Failure"):
        processor.nextalign_master(
            str(tmp_path / "ref.fa"), str(tmp_path / "MASTER1.fa"), str(tmp_path / "out")
        )

    assert [seen.arg_after("--min-seeds", i) for i in range(len(seen))] == [
        str(p["min_seeds"]) for p in processor.relaxation_profiles
    ]

def test_validate_alignment_integrity_checks_all_sequences(tmp_path: Path):
    """Both halves of the contract on one ten-column alignment: a clean file
    scores 1.0, and a single query that diverges in the conserved window drags
    the reported minimum down and fails the check."""
    processor = make(tmp_path)

    # 10 columns total. Conserved region is the final 30% (columns 7 to 9,
    # 0-indexed). Master "A T G C G T A A C G" -> downstream bases A, C, G.
    clean_fasta = tmp_path / "clean.fasta"
    clean_fasta.write_text(
        ">MASTER\nATGCGTAACG\n"
        ">QUERY1\nATGCGTAACG\n"   # matches ACG -> 100% identity
        ">QUERY2\nATGCGTAACG\n",  # matches ACG -> 100% identity
        encoding="utf-8",
    )
    passed, min_id = processor._validate_alignment_integrity(str(clean_fasta))
    assert passed is True
    assert min_id == 1.0

    dirty_fasta = tmp_path / "dirty.fasta"
    dirty_fasta.write_text(
        ">MASTER\nATGCGTAACG\n"
        ">QUERY1\nATGCGTAACG\n"   # matches ACG -> 100% identity
        ">QUERY2\nATGCGTGGGG\n",  # matches GGG -> 33% identity in the window
        encoding="utf-8",
    )
    passed, min_id = processor._validate_alignment_integrity(str(dirty_fasta))
    assert passed is False
    assert abs(min_id - 0.333333) < 1e-5


# ===========================================================================
# 5. Accession identity: path_to_basename and the master lookup
# ===========================================================================
#
# The pipeline keys on the BARE accession (meta_data.primary_accession,
# locus, FASTA headers, newick tips, filenames) and keeps the versioned form
# only in meta_data.accession_version, where GenBankFetcher uses it to notice a
# revised record. path_to_basename must produce the bare form for every spelling
# a filename might carry.


def test_basename_canonicalises_to_the_bare_accession():
    """`split('.')[0]` happened to give the right answer for 'NC_001542.fasta'
    and the wrong one for everything else."""
    assert NextalignAlignment.path_to_basename("/tmp/x/NC_001542.fasta") == "NC_001542"
    assert NextalignAlignment.path_to_basename("NC_001542.1.fasta") == "NC_001542"
    assert NextalignAlignment.path_to_basename("NC_001542.1.aligned.fasta") == "NC_001542"
    assert NextalignAlignment.path_to_basename("KX148218") == "KX148218"


def test_dot_leading_filename_does_not_produce_an_empty_accession():
    """An empty accession became a directory name and a command-line value."""
    assert NextalignAlignment.path_to_basename(".snakemake_timestamp") != ""


def test_basename_leaves_a_non_accession_group_name_alone():
    """Segmented runs name groups '<ref>_<segment>_<strand>', and influenza
    strain names contain dots. Neither is an accession and neither may be
    truncated."""
    assert NextalignAlignment.path_to_basename("NC_002016_4_plus.fa") == "NC_002016_4_plus"
    assert NextalignAlignment.path_to_basename("EPI_ISL_402124.fasta") == "EPI_ISL_402124"


@pytest.mark.parametrize("accession", ["NC_001542", "NC_001542.1", "KX148218", "PV547761.2"])
def test_basename_round_trips_the_master_lookup(tmp_path: Path, accession):
    """THE invariant that makes redundancy removal work. process() looks for
    <reference_aln>/<master>/<master>.aligned.fasta, where <master> comes from
    get_master_list() and the directory was named by path_to_basename(). Both
    sides must agree for every spelling either side might use."""
    processor = make(tmp_path, master_ref=accession)

    assert NextalignAlignment.path_to_basename(f"{accession}.fasta") == processor.get_master_list()[0]


def test_dedup_runs_for_a_versioned_master_accession(tmp_path: Path, monkeypatch):
    """The end-to-end consequence. A master list carrying 'NC_001542.1' - which
    is what an update DB or a hand-written reference list can hold - used to be
    compared against a directory named 'NC_001542', so redundancy removal was
    skipped for every master, silently and with no log line."""
    for name in ("query", "ref", "master"):
        (tmp_path / name).mkdir()
    write_fasta(tmp_path / "query" / "G1.fa", [("Q", "ATGC")])
    write_fasta(tmp_path / "ref" / "G1.fa", [("R", "ATGC")])
    write_fasta(tmp_path / "master" / "NC_001542.1.fasta", [("NC_001542.1", "ATGC")])

    processor = make(tmp_path, master_ref="NC_001542.1")
    deduped = []

    def fake_master(_query_path, ref_acc_path, out_root):
        accession = NextalignAlignment.path_to_basename(ref_acc_path)
        write_fasta(
            Path(out_root) / accession / f"{accession}.aligned.fasta",
            [("A", "ATGC"), ("B", "ATGC")],
        )

    class FakeDedup:
        def __init__(self, input_seq, _output_seq):
            self.input_seq = input_seq

        def remove_redundant_fasta(self):
            deduped.append(self.input_seq)

    monkeypatch.setattr(processor, "nextalign_query", lambda *_a: True)
    monkeypatch.setattr(processor, "nextalign_master", fake_master)
    monkeypatch.setattr(processor, "update_gb_matrix", lambda *_a: None)
    monkeypatch.setattr(na_module, "RemoveRedundantSequence", FakeDedup)

    processor.process()

    assert len(deduped) == 1, "redundancy removal never ran for the master reference"


def test_a_master_with_no_alignment_is_announced(tmp_path: Path, monkeypatch, capsys):
    """The silent skip is what hid the version mismatch for so long."""
    for name in ("query", "ref", "master"):
        (tmp_path / name).mkdir()
    processor = make(tmp_path, master_ref="NC_001542")
    monkeypatch.setattr(processor, "update_gb_matrix", lambda *_a: None)

    processor.process()

    assert "redundancy removal skipped for NC_001542" in capsys.readouterr().out


# ===========================================================================
# 6. get_master_list
# ===========================================================================


def test_keyboard_interrupt_is_not_swallowed(tmp_path: Path, monkeypatch):
    """The bare `except:` caught BaseException, so Ctrl-C during the master-list
    load left the run continuing with an empty master list."""
    master_file = tmp_path / "ref_list.tsv"
    master_file.write_text("primary_accession\taccession_type\n", encoding="utf-8")

    def boom(_path):
        raise KeyboardInterrupt()

    monkeypatch.setattr(na_module, "load_master_accessions_from_file", boom)

    with pytest.raises(KeyboardInterrupt):
        make(tmp_path, master_ref=str(master_file)).get_master_list()


def test_a_parse_failure_is_reported_not_silent(tmp_path: Path, monkeypatch, capsys):
    master_file = tmp_path / "ref_list.tsv"
    master_file.write_text("primary_accession\n", encoding="utf-8")

    def boom(_path):
        raise ValueError("malformed reference table")

    monkeypatch.setattr(na_module, "load_master_accessions_from_file", boom)

    assert make(tmp_path, master_ref=str(master_file)).get_master_list() == []
    assert "malformed reference table" in capsys.readouterr().out


def test_comma_separated_master_ref_is_split_and_stripped(tmp_path: Path):
    assert make(tmp_path, master_ref="A_1234, B_5678 ,,C_9012").get_master_list() == [
        "A_1234", "B_5678", "C_9012"
    ]


def test_update_db_takes_precedence_over_master_ref(tmp_path: Path, monkeypatch):
    monkeypatch.setattr(na_module, "load_master_accessions", lambda _db: ["NC_001542.1"])

    assert make(tmp_path, master_ref="KX148218", update_db="x.db").get_master_list() == ["NC_001542"]


def test_empty_update_db_result_falls_back_to_master_ref(tmp_path: Path, monkeypatch):
    monkeypatch.setattr(na_module, "load_master_accessions", lambda _db: [])

    assert make(tmp_path, master_ref="KX148218", update_db="x.db").get_master_list() == ["KX148218"]


def test_blank_master_entries_are_dropped(tmp_path: Path, monkeypatch):
    monkeypatch.setattr(na_module, "load_master_accessions", lambda _db: ["NC_001542", "", None, "  "])

    assert make(tmp_path, update_db="x.db").get_master_list() == ["NC_001542"]


# ===========================================================================
# 7. process(): pairing queries with references
# ===========================================================================


def test_query_group_without_a_reference_is_not_aligned_silently(tmp_path: Path, monkeypatch, capsys):
    """process() assumed ref_file == each_query_file and never checked that the
    reference existed, so a group with no reference was handed to nextalign
    with a path that was not there and vanished downstream without a log line."""
    query_dir = tmp_path / "query"
    ref_dir = tmp_path / "ref"
    query_dir.mkdir()
    ref_dir.mkdir()
    write_fasta(query_dir / "G1.fa", [("ACC1", "ATGC")])
    write_fasta(query_dir / "G2.fa", [("ACC2", "ATGC")])
    write_fasta(ref_dir / "G1.fa", [("R", "ATGC")])

    processor = make(tmp_path, query_dir=str(query_dir), ref_dir=str(ref_dir),
                     reference_alignment="provided.fa")
    pairs = []
    monkeypatch.setattr(processor, "nextalign_query",
                        lambda q, r, _o: pairs.append((Path(q).name, os.path.exists(r))) or True)
    monkeypatch.setattr(processor, "update_gb_matrix", lambda *_a: None)

    processor.process()

    assert [name for name, exists in pairs if not exists] == []
    assert "No reference sequence file for query group G2.fa" in capsys.readouterr().out


def test_an_unreferenced_group_is_excluded_in_the_matrix(tmp_path: Path, monkeypatch):
    """Not aligning it is not enough - those accessions must not pass into the
    tree build as though they had aligned."""
    query_dir = tmp_path / "query"
    ref_dir = tmp_path / "ref"
    query_dir.mkdir()
    ref_dir.mkdir()
    write_fasta(query_dir / "G2.fa", [("ACC2", "ATGC"), ("ACC3", "ATGC")])

    matrix = tmp_path / "gB_matrix_raw.tsv"
    write_matrix(matrix, [["ACC2", "query", "0", ""], ["ACC3", "query", "0", ""]])

    processor = make(tmp_path, gb_matrix=str(matrix), query_dir=str(query_dir),
                     ref_dir=str(ref_dir), reference_alignment="provided.fa")
    capture_commands(monkeypatch)

    processor.process()

    rows = {r["gi_number"]: r for r in read_matrix(matrix)}
    assert rows["ACC2"]["exclusion_status"] == "1"
    assert rows["ACC3"]["exclusion_status"] == "1"
    assert "No reference sequence file" in rows["ACC2"]["exclusion_criteria"]


def test_hidden_files_in_the_query_dir_are_skipped(tmp_path: Path, monkeypatch):
    """os.listdir returns editor swap files, .DS_Store and nextflow/snakemake
    bookkeeping entries, each of which was launched as its own alignment job."""
    query_dir = tmp_path / "query"
    ref_dir = tmp_path / "ref"
    query_dir.mkdir()
    ref_dir.mkdir()
    write_fasta(query_dir / "G1.fa", [("Q", "ATGC")])
    (query_dir / ".DS_Store").write_text("junk", encoding="utf-8")
    (query_dir / ".snakemake_timestamp").write_text("", encoding="utf-8")
    write_fasta(ref_dir / "G1.fa", [("R", "ATGC")])

    processor = make(tmp_path, query_dir=str(query_dir), ref_dir=str(ref_dir),
                     reference_alignment="provided.fa")
    aligned = []
    monkeypatch.setattr(processor, "nextalign_query", lambda q, _r, _o: aligned.append(Path(q).name) or True)
    monkeypatch.setattr(processor, "update_gb_matrix", lambda *_a: None)

    processor.process()

    assert aligned == ["G1.fa"]


def test_subdirectories_in_the_query_dir_are_skipped(tmp_path: Path, monkeypatch):
    query_dir = tmp_path / "query"
    ref_dir = tmp_path / "ref"
    (query_dir / "leftovers").mkdir(parents=True)
    ref_dir.mkdir()
    write_fasta(query_dir / "G1.fa", [("Q", "ATGC")])
    write_fasta(ref_dir / "G1.fa", [("R", "ATGC")])

    processor = make(tmp_path, query_dir=str(query_dir), ref_dir=str(ref_dir),
                     reference_alignment="provided.fa")
    aligned = []
    monkeypatch.setattr(processor, "nextalign_query", lambda q, _r, _o: aligned.append(Path(q).name) or True)
    monkeypatch.setattr(processor, "update_gb_matrix", lambda *_a: None)

    processor.process()

    assert aligned == ["G1.fa"]


def test_a_failed_query_alignment_does_not_pass_silently(tmp_path: Path, monkeypatch):
    """nextalign_query only printed on failure and returned None either way, so
    process() could not tell an aligned group from a failed one. When nextalign
    dies before writing errors.csv the sequences were never marked excluded and
    passed into the tree build as if they had aligned."""
    query_dir = tmp_path / "query"
    ref_dir = tmp_path / "ref"
    query_dir.mkdir()
    ref_dir.mkdir()
    write_fasta(query_dir / "G1.fa", [("ACC1", "ATGC")])
    write_fasta(ref_dir / "G1.fa", [("R", "ATGC")])

    matrix = tmp_path / "gB_matrix_raw.tsv"
    write_matrix(matrix, [["ACC1", "query", "0", ""]])

    capture_commands(monkeypatch, return_code=1)

    processor = make(tmp_path, gb_matrix=str(matrix), query_dir=str(query_dir),
                     ref_dir=str(ref_dir), reference_alignment="provided.fa")
    processor.process()

    assert read_matrix(matrix)[0]["exclusion_status"] == "1"


def test_a_successful_run_excludes_nothing(tmp_path: Path, monkeypatch):
    """The other half of the contract: a clean run must not mark anything."""
    query_dir = tmp_path / "query"
    ref_dir = tmp_path / "ref"
    query_dir.mkdir()
    ref_dir.mkdir()
    write_fasta(query_dir / "G1.fa", [("ACC1", "ATGC")])
    write_fasta(ref_dir / "G1.fa", [("R", "ATGC")])

    matrix = tmp_path / "gB_matrix_raw.tsv"
    write_matrix(matrix, [["ACC1", "query", "0", ""]])

    capture_commands(monkeypatch, return_code=0)

    processor = make(tmp_path, gb_matrix=str(matrix), query_dir=str(query_dir),
                     ref_dir=str(ref_dir), reference_alignment="provided.fa")
    processor.process()

    assert read_matrix(matrix)[0]["exclusion_status"] == "0"


def test_process_creates_both_alignment_output_directories(tmp_path: Path, monkeypatch):
    query_dir = tmp_path / "query"
    query_dir.mkdir()
    processor = make(tmp_path, query_dir=str(query_dir), reference_alignment="provided.fa")
    monkeypatch.setattr(processor, "update_gb_matrix", lambda *_a: None)

    processor.process()

    assert (tmp_path / "Nextalign" / "query_aln").is_dir()
    assert (tmp_path / "Nextalign" / "reference_aln").is_dir()


def test_a_missing_query_directory_is_reported_not_crashed(tmp_path: Path, monkeypatch):
    """os.listdir on an absent directory used to raise FileNotFoundError deep
    in process(); an empty upstream stage should be diagnosable."""
    processor = make(tmp_path, query_dir=str(tmp_path / "absent"),
                     reference_alignment="provided.fa")
    monkeypatch.setattr(processor, "update_gb_matrix", lambda *_a: None)

    processor.process()  # must not raise


def test_query_groups_are_aligned_in_a_stable_order(tmp_path: Path, monkeypatch):
    """os.listdir order is arbitrary; a reproducible run needs a fixed order."""
    query_dir = tmp_path / "query"
    ref_dir = tmp_path / "ref"
    query_dir.mkdir()
    ref_dir.mkdir()
    for name in ("G3.fa", "G1.fa", "G2.fa"):
        write_fasta(query_dir / name, [("Q", "ATGC")])
        write_fasta(ref_dir / name, [("R", "ATGC")])

    processor = make(tmp_path, query_dir=str(query_dir), ref_dir=str(ref_dir),
                     reference_alignment="provided.fa")
    aligned = []
    monkeypatch.setattr(processor, "nextalign_query", lambda q, _r, _o: aligned.append(Path(q).name) or True)
    monkeypatch.setattr(processor, "update_gb_matrix", lambda *_a: None)

    processor.process()

    assert aligned == ["G1.fa", "G2.fa", "G3.fa"]


def test_process_reference_alignment_mode_runs_query_only(tmp_path: Path, monkeypatch):
    """With a supplied reference alignment the master leg is skipped entirely,
    and the matrix is still updated."""
    matrix = tmp_path / "gB_matrix_raw.tsv"
    write_matrix(matrix, [["ACC_OK", "query", "0", ""], ["ACC_FAIL", "query", "0", ""]])

    query_dir = tmp_path / "query"
    ref_dir = tmp_path / "ref"
    write_fasta(query_dir / "REF1.fa", [("Q", "ATGC")])
    write_fasta(ref_dir / "REF1.fa", [("R", "ATGC")])

    processor = make(tmp_path, query_dir=str(query_dir), ref_dir=str(ref_dir),
                     reference_alignment="provided_alignment.fa")

    calls = {"query": 0, "master": 0, "update": 0}
    monkeypatch.setattr(processor, "nextalign_query",
                        lambda *_a, **_k: calls.__setitem__("query", calls["query"] + 1))
    monkeypatch.setattr(processor, "nextalign_master",
                        lambda *_a, **_k: calls.__setitem__("master", calls["master"] + 1))
    monkeypatch.setattr(processor, "update_gb_matrix",
                        lambda *_a, **_k: calls.__setitem__("update", calls["update"] + 1))

    processor.process()

    assert calls == {"query": 1, "master": 0, "update": 1}


def test_process_non_reference_alignment_runs_master_path(tmp_path: Path, monkeypatch):
    """Without a supplied alignment the master is aligned and deduplicated as
    well as the queries."""
    matrix = tmp_path / "gB_matrix_raw.tsv"
    write_matrix(matrix, [["ACC_OK", "query", "0", ""], ["ACC_FAIL", "query", "0", ""]])

    query_dir = tmp_path / "query"
    ref_dir = tmp_path / "ref"
    master_dir = tmp_path / "master"
    write_fasta(query_dir / "REF1.fa", [("Q", "ATGC")])
    write_fasta(ref_dir / "REF1.fa", [("R", "ATGC")])
    write_fasta(master_dir / "MASTER1.fasta", [("MASTER1", "ATGC")])

    processor = make(tmp_path, query_dir=str(query_dir), ref_dir=str(ref_dir),
                     master_seq_dir=str(master_dir))

    calls = {"query": 0, "master": 0, "update": 0, "dedup": 0}

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

    monkeypatch.setattr(processor, "nextalign_query",
                        lambda *_a, **_k: calls.__setitem__("query", calls["query"] + 1))
    monkeypatch.setattr(processor, "nextalign_master", fake_master)
    monkeypatch.setattr(processor, "update_gb_matrix",
                        lambda *_a, **_k: calls.__setitem__("update", calls["update"] + 1))
    monkeypatch.setattr(na_module, "RemoveRedundantSequence", _FakeRemoveRedundantSequence)

    processor.process()

    assert calls == {"query": 1, "master": 1, "update": 1, "dedup": 1}
