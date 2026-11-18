import csv
from pathlib import Path

from SoftwareVersion import SoftwareVersionChecker


class FakeResult:
    def __init__(self, stdout="", stderr=""):
        self.stdout = stdout
        self.stderr = stderr


def test_get_software_version_and_tsv_write(tmp_path: Path, monkeypatch):
    def fake_run(command, stdout, stderr, text):
        cmd = " ".join(command)
        if "nextalign" in cmd:
            return FakeResult(stdout="nextalign version 2.0.0\n")
        if "blastn" in cmd:
            return FakeResult(stdout="blastn: 2.15.0+\n")
        if "mafft" in cmd:
            return FakeResult(stderr="mafft v7.520\n")
        if "iqtree2" in cmd:
            return FakeResult(stdout="IQ-TREE multicore version 2.3.5\n")
        if "FastTree" in cmd:
            return FakeResult(stdout="FastTree Version 2.1\n")
        if "python" in cmd:
            return FakeResult(stdout="Python 3.11.9\n")
        return FakeResult(stdout="unknown\n")

    monkeypatch.setattr("SoftwareVersion.subprocess.run", fake_run)

    checker = SoftwareVersionChecker(str(tmp_path), "Software_info", "software_info.tsv")
    out_file = tmp_path / "Software_info" / "software_info.tsv"

    assert out_file.exists()

    with out_file.open("r", encoding="utf-8", newline="") as f:
        rows = list(csv.reader(f, delimiter="\t"))

    assert rows[0] == ["Software", "Version"]
    names = [r[0] for r in rows[1:]]
    assert "Nextalign" in names
    assert "BLAST" in names
    assert "Python" in names
    assert "Time of creation" in names
    assert "vgtk version" in names
    assert dict(rows[1:])["vgtk version"] == "v.1.1.0"

    assert checker.software_versions["MAFFT"] == "mafft v7.520"


def test_get_software_version_not_installed(monkeypatch):
    def fake_run(*args, **kwargs):
        raise FileNotFoundError()

    checker = SoftwareVersionChecker.__new__(SoftwareVersionChecker)
    monkeypatch.setattr("SoftwareVersion.subprocess.run", fake_run)

    assert checker.get_software_version("Whatever", ["whatever", "--version"]) == "Not Installed"
