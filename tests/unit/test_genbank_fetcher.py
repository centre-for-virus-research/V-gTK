import csv
import sqlite3
from pathlib import Path

import pytest

import GenBankFetcher as gb_module
from GenBankFetcher import GenBankFetcher


class _FakeResponse:
    def __init__(self, json_data=None, text_data="", status_code=200):
        self._json_data = json_data
        self.text = text_data
        self.status_code = status_code

    def json(self):
        return self._json_data

    def raise_for_status(self):
        if self.status_code >= 400:
            raise Exception(f"HTTP {self.status_code}")


def _write_db(path: Path, meta_col: str, meta_values, excluded_values=None):
    excluded_values = excluded_values or []
    conn = sqlite3.connect(str(path))
    try:
        cur = conn.cursor()
        cur.execute(f"CREATE TABLE meta_data ({meta_col} TEXT, exclusion_status TEXT, exclusion_criteria TEXT)")
        cur.executemany(
            f"INSERT INTO meta_data({meta_col}, exclusion_status, exclusion_criteria) VALUES (?, '', '')",
            [(v,) for v in meta_values],
        )
        cur.executemany(
            f"INSERT INTO meta_data({meta_col}, exclusion_status, exclusion_criteria) VALUES (?, '1', 'excluded')",
            [(v,) for v in excluded_values],
        )
        conn.commit()
    finally:
        conn.close()


def test_fetch_ids_paginates_and_strips_versions(tmp_path: Path, monkeypatch):
    fetcher = GenBankFetcher(
        taxid="11292",
        base_url="https://example/",
        email="x@y.com",
        output_dir=str(tmp_path),
        batch_size=2,
        sleep_time=0,
        base_dir="GenBank-XML",
        update_file=None,
    )

    def fake_get(url):
        if "retmax=0" in url:
            return _FakeResponse({"esearchresult": {"count": "3", "webenv": "W", "querykey": "1"}})
        if "retstart=0" in url:
            return _FakeResponse({"esearchresult": {"idlist": ["A.1", "B.2"]}})
        if "retstart=2" in url:
            return _FakeResponse({"esearchresult": {"idlist": ["C.9"]}})
        raise AssertionError(f"Unexpected URL: {url}")

    monkeypatch.setattr("GenBankFetcher.requests.get", fake_get)
    ids = fetcher.fetch_ids()
    assert ids == ["A", "B", "C"]


def test_fetch_genbank_data_adds_ref_list_and_saves_batches(tmp_path: Path, monkeypatch):
    ref_file = tmp_path / "refs.tsv"
    ref_file.write_text("REF1\tmaster\n", encoding="utf-8")

    fetcher = GenBankFetcher(
        taxid="11292",
        base_url="https://example/",
        email="x@y.com",
        output_dir=str(tmp_path),
        batch_size=100,
        sleep_time=0,
        base_dir="GenBank-XML",
        update_file=None,
        ref_list=str(ref_file),
    )
    fetcher.efetch_batch_size = 2

    calls = []

    def fake_get(url):
        calls.append(url)
        return _FakeResponse(json_data={}, text_data="<GBSet></GBSet>")

    saved = []

    def fake_save(data, marker):
        saved.append((data, marker))

    monkeypatch.setattr("GenBankFetcher.requests.get", fake_get)
    monkeypatch.setattr(fetcher, "save_data", fake_save)

    fetcher.fetch_genbank_data(["A", "A", "B"])

    assert len(saved) >= 1
    assert all(data == "<GBSet></GBSet>" for data, _ in saved)
    assert any("efetch.fcgi" in c for c in calls)


def test_fetch_genbank_data_adds_headered_ref_list(tmp_path: Path, monkeypatch):
    ref_file = tmp_path / "refs.tsv"
    ref_file.write_text(
        "primary_accession\tstatus\tsegment\nREF1\tmaster\t1\nREF2\treference\t1\n",
        encoding="utf-8",
    )

    fetcher = GenBankFetcher(
        taxid="11292",
        base_url="https://example/",
        email="x@y.com",
        output_dir=str(tmp_path),
        batch_size=100,
        sleep_time=0,
        base_dir="GenBank-XML",
        update_file=None,
        ref_list=str(ref_file),
    )
    fetcher.efetch_batch_size = 10

    called_urls = []

    def fake_get(url):
        called_urls.append(url)
        return _FakeResponse(json_data={}, text_data="<GBSet></GBSet>")

    monkeypatch.setattr("GenBankFetcher.requests.get", fake_get)
    monkeypatch.setattr(fetcher, "save_data", lambda *_: None)

    fetcher.fetch_genbank_data(["A"])

    assert any("REF1" in url and "REF2" in url for url in called_urls)


def test_update_requires_meta_data_table(tmp_path: Path):
    bad_db = tmp_path / "bad.db"
    conn = sqlite3.connect(str(bad_db))
    try:
        conn.execute("CREATE TABLE wrong_table (x TEXT)")
        conn.commit()
    finally:
        conn.close()

    fetcher = GenBankFetcher(
        taxid="11292",
        base_url="https://example/",
        email="x@y.com",
        output_dir=str(tmp_path),
        batch_size=100,
        sleep_time=0,
        base_dir="GenBank-XML",
        update_file=None,
    )

    with pytest.raises(ValueError, match="meta_data"):
        fetcher.update(str(bad_db))


def test_update_from_db_fetches_only_updated_and_new_accessions(tmp_path: Path, monkeypatch):
    db_path = tmp_path / "prev.db"
    _write_db(
        db_path,
        meta_col="accession_version",
        meta_values=["A.1", "B.2", "REF1.1"],
        excluded_values=["X1"],
    )

    fetcher = GenBankFetcher(
        taxid="11292",
        base_url="https://example/",
        email="x@y.com",
        output_dir=str(tmp_path),
        batch_size=100,
        sleep_time=0,
        base_dir="GenBank-XML",
        update_file=None,
    )

    monkeypatch.setattr(fetcher, "fetch_accs", lambda: ["A.1", "A.2", "B.2", "C.1", "X1.2"])

    captured = {}

    def fake_fetch_genbank_data(ids):
        captured["ids"] = ids

    monkeypatch.setattr(fetcher, "fetch_genbank_data", fake_fetch_genbank_data)

    fetcher.update(str(db_path))

    assert captured["ids"] == ["A.2", "C.1"]

    log_path = tmp_path / "update_accessions.tsv"
    assert log_path.exists()
    rows = list(csv.DictReader(log_path.open("r", encoding="utf-8"), delimiter="\t"))
    assert {tuple((r["old_accession_version"], r["new_accession_version"])) for r in rows} == {
        ("A.1", "A.2"),
        ("NA", "C.1"),
    }


def test_update_from_db_fallbacks_to_primary_accession_column(tmp_path: Path, monkeypatch):
    db_path = tmp_path / "prev_primary.db"
    _write_db(
        db_path,
        meta_col="primary_accession",
        meta_values=["A", "B"],
        excluded_values=[],
    )

    fetcher = GenBankFetcher(
        taxid="11292",
        base_url="https://example/",
        email="x@y.com",
        output_dir=str(tmp_path),
        batch_size=100,
        sleep_time=0,
        base_dir="GenBank-XML",
        update_file=None,
    )

    monkeypatch.setattr(fetcher, "fetch_accs", lambda: ["A.1", "B.1", "C.1"])
    captured = {}
    monkeypatch.setattr(fetcher, "fetch_genbank_data", lambda ids: captured.setdefault("ids", ids))

    fetcher.update(str(db_path))

    assert captured["ids"] == ["C.1"]


def test_split_accession_version_handles_edge_cases(tmp_path: Path):
    fetcher = GenBankFetcher(
        taxid="11292",
        base_url="https://example/",
        email="x@y.com",
        output_dir=str(tmp_path),
        batch_size=100,
        sleep_time=0,
        base_dir="GenBank-XML",
        update_file=None,
    )

    assert fetcher._split_accession_version(None) == (None, None)
    assert fetcher._split_accession_version("") == (None, None)
    assert fetcher._split_accession_version("ABC123") == ("ABC123", None)
    assert fetcher._split_accession_version("ABC123.7") == ("ABC123", 7)
    assert fetcher._split_accession_version("ABC123.X") == ("ABC123", None)


def test_detect_meta_data_acc_col_raises_when_no_supported_columns(tmp_path: Path):
    db_path = tmp_path / "no_supported_cols.db"
    conn = sqlite3.connect(str(db_path))
    try:
        conn.execute("CREATE TABLE meta_data (foo TEXT, bar TEXT)")
        conn.commit()

        fetcher = GenBankFetcher(
            taxid="11292",
            base_url="https://example/",
            email="x@y.com",
            output_dir=str(tmp_path),
            batch_size=100,
            sleep_time=0,
            base_dir="GenBank-XML",
            update_file=None,
        )

        with pytest.raises(ValueError, match="Could not find an accession column"):
            fetcher._detect_meta_data_acc_col(conn, preferred="accession_version")
    finally:
        conn.close()


def test_compute_missing_ids_skips_excluded_and_unversioned_existing(tmp_path: Path):
    fetcher = GenBankFetcher(
        taxid="11292",
        base_url="https://example/",
        email="x@y.com",
        output_dir=str(tmp_path),
        batch_size=100,
        sleep_time=0,
        base_dir="GenBank-XML",
        update_file=None,
    )

    ncbi_ids = ["A.1", "A.2", "B.1", "EXC.2", "NEW.1"]
    meta_accessions = ["A", "B.1"]
    excluded_primary = {"EXC"}

    missing_ids, updated_versions, new_accessions = fetcher._compute_missing_ids(
        ncbi_ids,
        meta_accessions,
        excluded_primary,
    )

    assert missing_ids == ["NEW.1"]
    assert updated_versions == []
    assert new_accessions == ["NEW.1"]


def test_save_data_uses_incrementing_suffix_when_file_exists(tmp_path: Path):
    fetcher = GenBankFetcher(
        taxid="11292",
        base_url="https://example/",
        email="x@y.com",
        output_dir=str(tmp_path),
        batch_size=100,
        sleep_time=0,
        base_dir="GenBank-XML",
        update_file=None,
    )

    xml_dir = tmp_path / "GenBank-XML"
    xml_dir.mkdir(parents=True, exist_ok=True)
    first_path = xml_dir / "batch-100.xml"
    first_path.write_text("old", encoding="utf-8")

    fetcher.save_data("<GBSet></GBSet>", 100)

    second_path = xml_dir / "batch-100_1.xml"
    assert second_path.exists()
    assert second_path.read_text(encoding="utf-8") == "<GBSet></GBSet>"


def test_update_with_no_missing_ids_does_not_download(tmp_path: Path, monkeypatch):
    db_path = tmp_path / "prev_all_up_to_date.db"
    _write_db(
        db_path,
        meta_col="accession_version",
        meta_values=["A.1", "B.2"],
        excluded_values=[],
    )

    fetcher = GenBankFetcher(
        taxid="11292",
        base_url="https://example/",
        email="x@y.com",
        output_dir=str(tmp_path),
        batch_size=100,
        sleep_time=0,
        base_dir="GenBank-XML",
        update_file=None,
    )

    monkeypatch.setattr(fetcher, "fetch_accs", lambda: ["A.1", "B.2"])

    called = {"download": False}

    def fake_fetch(ids):
        called["download"] = True

    monkeypatch.setattr(fetcher, "fetch_genbank_data", fake_fetch)

    fetcher.update(str(db_path))

    assert called["download"] is False


def test_fetch_ids_retries_on_connection_error_then_succeeds(tmp_path: Path, monkeypatch):
    fetcher = GenBankFetcher(
        taxid="11292",
        base_url="https://example/",
        email="x@y.com",
        output_dir=str(tmp_path),
        batch_size=2,
        sleep_time=1,
        base_dir="GenBank-XML",
        update_file=None,
    )

    calls = {"retstart0": 0}
    waits = []

    def fake_sleep(seconds):
        waits.append(seconds)

    def fake_get(url):
        if "retmax=0" in url:
            return _FakeResponse({"esearchresult": {"count": "2", "webenv": "W", "querykey": "1"}})
        if "retstart=0" in url:
            calls["retstart0"] += 1
            if calls["retstart0"] == 1:
                raise gb_module.requests.exceptions.ConnectionError("transient")
            return _FakeResponse({"esearchresult": {"idlist": ["A.1", "B.2"]}})
        raise AssertionError(f"Unexpected URL: {url}")

    monkeypatch.setattr("GenBankFetcher.requests.get", fake_get)
    monkeypatch.setattr("GenBankFetcher.sleep", fake_sleep)

    ids = fetcher.fetch_ids()
    assert ids == ["A", "B"]
    assert waits == [1]


def test_fetch_ids_raises_after_max_retries(tmp_path: Path, monkeypatch):
    fetcher = GenBankFetcher(
        taxid="11292",
        base_url="https://example/",
        email="x@y.com",
        output_dir=str(tmp_path),
        batch_size=2,
        sleep_time=0,
        base_dir="GenBank-XML",
        update_file=None,
    )

    def fake_get(url):
        if "retmax=0" in url:
            return _FakeResponse({"esearchresult": {"count": "2", "webenv": "W", "querykey": "1"}})
        if "retstart=0" in url:
            raise gb_module.requests.exceptions.ConnectionError("always fails")
        raise AssertionError(f"Unexpected URL: {url}")

    monkeypatch.setattr("GenBankFetcher.requests.get", fake_get)
    monkeypatch.setattr("GenBankFetcher.sleep", lambda *_: None)

    with pytest.raises(gb_module.requests.exceptions.ConnectionError):
        fetcher.fetch_ids()


def test_fetch_accs_retries_on_incomplete_json_then_succeeds(tmp_path: Path, monkeypatch):
    fetcher = GenBankFetcher(
        taxid="11292",
        base_url="https://example/",
        email="x@y.com",
        output_dir=str(tmp_path),
        batch_size=2,
        sleep_time=0,
        base_dir="GenBank-XML",
        update_file=None,
    )

    calls = {"retstart0": 0}

    def fake_get(url):
        if "retmax=0" in url:
            return _FakeResponse({"esearchresult": {"count": "2", "webenv": "W", "querykey": "1"}})
        if "retstart=0" in url:
            calls["retstart0"] += 1
            if calls["retstart0"] == 1:
                return _FakeResponse({"wrong": {}})
            return _FakeResponse({"esearchresult": {"idlist": ["A.1", "B.2"]}})
        raise AssertionError(f"Unexpected URL: {url}")

    monkeypatch.setattr("GenBankFetcher.requests.get", fake_get)
    monkeypatch.setattr("GenBankFetcher.sleep", lambda *_: None)

    ids = fetcher.fetch_accs()
    assert ids == ["A.1", "B.2"]


def test_fetch_genbank_data_429_retry_then_success(tmp_path: Path, monkeypatch):
    fetcher = GenBankFetcher(
        taxid="11292",
        base_url="https://example/",
        email="x@y.com",
        output_dir=str(tmp_path),
        batch_size=100,
        sleep_time=1,
        base_dir="GenBank-XML",
        update_file=None,
    )
    fetcher.efetch_batch_size = 2

    class _Resp429:
        status_code = 429

    class _Transient429Response:
        text = ""

        def raise_for_status(self):
            raise gb_module.requests.exceptions.HTTPError("429", response=_Resp429())

    attempts = {"count": 0}
    waits = []
    saved = []

    def fake_sleep(seconds):
        waits.append(seconds)

    def fake_get(_url):
        attempts["count"] += 1
        if attempts["count"] == 1:
            return _Transient429Response()
        return _FakeResponse(json_data={}, text_data="<GBSet></GBSet>")

    monkeypatch.setattr("GenBankFetcher.requests.get", fake_get)
    monkeypatch.setattr("GenBankFetcher.sleep", fake_sleep)
    monkeypatch.setattr(fetcher, "save_data", lambda data, marker: saved.append((data, marker)))

    fetcher.fetch_genbank_data(["A", "B"])

    assert len(saved) == 1
    assert waits[0] == 11


def test_download_test_run_fallbacks_to_empty_ids_on_error(tmp_path: Path, monkeypatch):
    fetcher = GenBankFetcher(
        taxid="11292",
        base_url="https://example/",
        email="x@y.com",
        output_dir=str(tmp_path),
        batch_size=100,
        sleep_time=0,
        base_dir="GenBank-XML",
        update_file=None,
        test_run=True,
    )

    monkeypatch.setattr("GenBankFetcher.requests.get", lambda _url: (_ for _ in ()).throw(Exception("boom")))
    captured = {}
    monkeypatch.setattr(fetcher, "fetch_genbank_data", lambda ids: captured.setdefault("ids", ids))

    fetcher.download()

    assert captured["ids"] == []


def test_download_test_run_uses_configured_test_fetch_count(tmp_path: Path, monkeypatch):
    fetcher = GenBankFetcher(
        taxid="11292",
        base_url="https://example/",
        email="x@y.com",
        output_dir=str(tmp_path),
        batch_size=100,
        sleep_time=0,
        base_dir="GenBank-XML",
        update_file=None,
        test_run=True,
        test_fetch_count=100,
    )

    requested_urls = []
    sampled = {}
    captured = {}

    def fake_get(url):
        requested_urls.append(url)
        return _FakeResponse({"esearchresult": {"idlist": [f"ACC{i}" for i in range(150)]}})

    def fake_sample(population, k):
        sampled["population_size"] = len(population)
        sampled["k"] = k
        return list(population)[:k]

    monkeypatch.setattr("GenBankFetcher.requests.get", fake_get)
    monkeypatch.setattr(gb_module.random, "sample", fake_sample)
    monkeypatch.setattr(fetcher, "save_data", lambda data, marker: captured.setdefault("markers", []).append(marker))

    fetcher.download()

    assert any("retmax=100" in url for url in requested_urls)
    assert sampled == {"population_size": 150, "k": 100}
    assert captured["markers"] == [10, 20, 30, 40, 50, 60, 70, 80, 90, 100]


def test_update_does_not_fetch_refs_from_ref_list_in_db_backed_mode(tmp_path: Path, monkeypatch):
    db_path = tmp_path / "prev_refs.db"
    _write_db(
        db_path,
        meta_col="accession_version",
        meta_values=["REF1.1", "Q1.1"],
        excluded_values=[],
    )

    ref_file = tmp_path / "refs.tsv"
    ref_file.write_text("REF1\tmaster\nREF2\treference\n", encoding="utf-8")

    fetcher = GenBankFetcher(
        taxid="11292",
        base_url="https://example/",
        email="x@y.com",
        output_dir=str(tmp_path),
        batch_size=100,
        sleep_time=0,
        base_dir="GenBank-XML",
        update_file=str(db_path),
        ref_list=str(ref_file),
    )

    monkeypatch.setattr(fetcher, "fetch_accs", lambda: ["Q1.1"])
    captured = {}
    monkeypatch.setattr(fetcher, "fetch_genbank_data", lambda ids: captured.setdefault("ids", ids))

    fetcher.update(str(db_path))

    assert "ids" not in captured


def test_update_ref_list_parser_ignores_malformed_and_excluded_refs(tmp_path: Path, monkeypatch, basic_update_db: Path):
    ref_file = tmp_path / "refs.tsv"
    ref_file.write_text(
        "REF1\tmaster\n"
        "REF_NEW\treference\n"
        "BROKEN_ONLY_ONE_COL\n"
        "Q_EXCL\treference\n",
        encoding="utf-8",
    )

    fetcher = GenBankFetcher(
        taxid="11292",
        base_url="https://example/",
        email="x@y.com",
        output_dir=str(tmp_path),
        batch_size=100,
        sleep_time=0,
        base_dir="GenBank-XML",
        update_file=str(basic_update_db),
        ref_list=str(ref_file),
    )

    monkeypatch.setattr(fetcher, "fetch_accs", lambda: ["Q_OLD.1"])  # no NCBI deltas
    captured = {}
    monkeypatch.setattr(fetcher, "fetch_genbank_data", lambda ids: captured.setdefault("ids", ids))

    fetcher.update(str(basic_update_db))

    # Update mode is DB-backed for references: ref_list should not force fetching REF_NEW.
    assert "ids" not in captured


def test_update_test_run_samples_brand_new_accessions_first(tmp_path: Path, monkeypatch):
    db_path = tmp_path / "prev_test_update.db"
    _write_db(
        db_path,
        meta_col="accession_version",
        meta_values=["A.1", "B.1"],
        excluded_values=[],
    )

    fetcher = GenBankFetcher(
        taxid="11292",
        base_url="https://example/",
        email="x@y.com",
        output_dir=str(tmp_path),
        batch_size=100,
        sleep_time=0,
        base_dir="GenBank-XML",
        update_file=str(db_path),
        test_run=True,
    )

    monkeypatch.setattr(fetcher, "iter_accs", lambda: iter([["A.2", "NEW1.1", "NEW2.1", "NEW3.1"]]))

    captured = {}
    monkeypatch.setattr(fetcher, "fetch_genbank_data", lambda ids: captured.setdefault("ids", ids))

    fetcher.update(str(db_path))

    # Should select only brand-new primary accessions, excluding A.2 update candidate
    assert set(captured["ids"]) == {"NEW1.1", "NEW2.1", "NEW3.1"}


def test_update_test_run_stops_after_first_100_missing_accessions(tmp_path: Path, monkeypatch):
    db_path = tmp_path / "prev_test_update.db"
    _write_db(
        db_path,
        meta_col="accession_version",
        meta_values=["OLD.1"],
        excluded_values=[],
    )

    fetcher = GenBankFetcher(
        taxid="11292",
        base_url="https://example/",
        email="x@y.com",
        output_dir=str(tmp_path),
        batch_size=50,
        sleep_time=0,
        base_dir="GenBank-XML",
        update_file=str(db_path),
        test_run=True,
    )

    pages_seen = []

    def fake_iter_accs():
        for page_idx in range(3):
            page = [f"NEW{page_idx}_{i}.1" for i in range(50)]
            pages_seen.append(page_idx)
            yield page

    monkeypatch.setattr(fetcher, "iter_accs", fake_iter_accs)

    captured = {}
    monkeypatch.setattr(fetcher, "fetch_genbank_data", lambda ids: captured.setdefault("ids", ids))

    fetcher.update(str(db_path))

    assert pages_seen == [0, 1]
    assert len(captured["ids"]) == 100
    assert captured["ids"][0] == "NEW0_0.1"
    assert captured["ids"][-1] == "NEW1_49.1"


def test_update_test_run_respects_configured_test_fetch_count(tmp_path: Path, monkeypatch):
    db_path = tmp_path / "prev_test_update.db"
    _write_db(
        db_path,
        meta_col="accession_version",
        meta_values=["OLD.1"],
        excluded_values=[],
    )

    fetcher = GenBankFetcher(
        taxid="11292",
        base_url="https://example/",
        email="x@y.com",
        output_dir=str(tmp_path),
        batch_size=50,
        sleep_time=0,
        base_dir="GenBank-XML",
        update_file=str(db_path),
        test_run=True,
        test_fetch_count=25,
    )

    pages_seen = []

    def fake_iter_accs():
        for page_idx in range(2):
            page = [f"NEW{page_idx}_{i}.1" for i in range(50)]
            pages_seen.append(page_idx)
            yield page

    monkeypatch.setattr(fetcher, "iter_accs", fake_iter_accs)

    captured = {}
    monkeypatch.setattr(fetcher, "fetch_genbank_data", lambda ids: captured.setdefault("ids", ids))

    fetcher.update(str(db_path))

    assert pages_seen == [0]
    assert len(captured["ids"]) == 25
    assert captured["ids"][0] == "NEW0_0.1"
    assert captured["ids"][-1] == "NEW0_24.1"


def test_fetch_genbank_data_does_not_test_sample_when_update_mode(tmp_path: Path, monkeypatch):
    fetcher = GenBankFetcher(
        taxid="11292",
        base_url="https://example/",
        email="x@y.com",
        output_dir=str(tmp_path),
        batch_size=100,
        sleep_time=0,
        base_dir="GenBank-XML",
        update_file=str(tmp_path / "dummy.db"),
        test_run=True,
    )
    fetcher.efetch_batch_size = 100

    sample_called = {"called": False}

    def fake_sample(_population, _k):
        sample_called["called"] = True
        return []

    monkeypatch.setattr(gb_module.random, "sample", fake_sample)
    monkeypatch.setattr("GenBankFetcher.requests.get", lambda _url: _FakeResponse(json_data={}, text_data="<GBSet></GBSet>"))

    seen = []
    monkeypatch.setattr(fetcher, "save_data", lambda data, marker: seen.append((data, marker)))

    fetcher.fetch_genbank_data(["X1.1", "X2.1"])

    assert sample_called["called"] is False
    assert len(seen) == 1
