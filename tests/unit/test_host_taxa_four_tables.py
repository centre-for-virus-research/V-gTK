"""Extensive TDD unit and integration tests for Workstream 3:
All 4 host taxa SQLite tables across every virus (fresh & --update_db),
composite update keys (preserving synonyms), post-GISAID ValidateMatrix
host resolution, semicolon qualifier stripping, and generic 'Duck' -> Anatidae (8830).
"""

from __future__ import annotations

import csv
import sqlite3
import subprocess
import sys
from pathlib import Path

import pytest

from CreateSqliteDB import CreateSqliteDB
from HostTaxaTable import HostTaxaTable
from ValidateMatrix import ValidateMatrix

REPO_ROOT = Path(__file__).resolve().parents[2]
VALIDATE_DB_TREE = REPO_ROOT / "scripts" / "ValidateDbTree.py"
INFLUENZA_HOST_ALIASES = REPO_ROOT / "generic" / "influenza" / "host-aliases.tsv"


def _write_tsv(path: Path, header, rows) -> Path:
    with path.open("w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(header)
        w.writerows(rows)
    return path


def _write_mini_taxdump(taxa_dir: Path) -> tuple[Path, Path]:
    taxa_dir.mkdir(parents=True, exist_ok=True)
    names = taxa_dir / "names.dmp"
    nodes = taxa_dir / "nodes.dmp"
    names.write_text(
        "\n".join(
            [
                "1\t|\troot\t|\t\t|\tscientific name\t|",
                "2759\t|\tEukaryota\t|\t\t|\tscientific name\t|",
                "7711\t|\tChordata\t|\t\t|\tscientific name\t|",
                "40674\t|\tMammalia\t|\t\t|\tscientific name\t|",
                "9443\t|\tPrimates\t|\t\t|\tscientific name\t|",
                "9604\t|\tHominidae\t|\t\t|\tscientific name\t|",
                "9605\t|\tHomo\t|\t\t|\tscientific name\t|",
                "9606\t|\tHomo sapiens\t|\t\t|\tscientific name\t|",
                "9606\t|\thuman\t|\t\t|\tgenbank common name\t|",
                "9606\t|\tman\t|\t\t|\tcommon name\t|",
                "9913\t|\tBos taurus\t|\t\t|\tscientific name\t|",
                "9913\t|\tdairy cow\t|\t\t|\tcommon name\t|",
                "9823\t|\tSus scrofa\t|\t\t|\tscientific name\t|",
                "9823\t|\tswine\t|\t\t|\tcommon name\t|",
                "9681\t|\tFelidae\t|\t\t|\tscientific name\t|",
                "8782\t|\tAves\t|\t\t|\tscientific name\t|",
                "8782\t|\tavian\t|\t\t|\tcommon name\t|",
                "8825\t|\tAnseriformes\t|\t\t|\tscientific name\t|",
                "8830\t|\tAnatidae\t|\t\t|\tscientific name\t|",
                "8835\t|\tAnas\t|\t\t|\tscientific name\t|",
                "8839\t|\tAnas platyrhynchos\t|\t\t|\tscientific name\t|",
                "8839\t|\tmallard\t|\t\t|\tgenbank common name\t|",
                "8839\t|\tduck\t|\t\t|\tcommon name\t|",
                "9031\t|\tGallus gallus\t|\t\t|\tscientific name\t|",
                "9031\t|\tchicken\t|\t\t|\tgenbank common name\t|",
                "54971\t|\tArenaria interpres\t|\t\t|\tscientific name\t|",
                "54971\t|\truddy turnstone\t|\t\t|\tgenbank common name\t|",
                "33554\t|\tCarnivora\t|\t\t|\tscientific name\t|",
                "30419\t|\tAccipitridae\t|\t\t|\tscientific name\t|",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    nodes.write_text(
        "\n".join(
            [
                "1\t|\t1\t|\tno rank\t|",
                "2759\t|\t1\t|\tsuperkingdom\t|",
                "7711\t|\t2759\t|\tphylum\t|",
                "40674\t|\t7711\t|\tclass\t|",
                "9443\t|\t40674\t|\torder\t|",
                "9604\t|\t9443\t|\tfamily\t|",
                "9605\t|\t9604\t|\tgenus\t|",
                "9606\t|\t9605\t|\tspecies\t|",
                "9913\t|\t40674\t|\tspecies\t|",
                "9823\t|\t40674\t|\tspecies\t|",
                "33554\t|\t40674\t|\torder\t|",
                "9681\t|\t33554\t|\tfamily\t|",
                "8782\t|\t7711\t|\tclass\t|",
                "8825\t|\t8782\t|\torder\t|",
                "8830\t|\t8825\t|\tfamily\t|",
                "8835\t|\t8830\t|\tgenus\t|",
                "8839\t|\t8835\t|\tspecies\t|",
                "9031\t|\t8782\t|\tspecies\t|",
                "54971\t|\t8782\t|\tspecies\t|",
                "30419\t|\t8782\t|\tfamily\t|",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    return names, nodes


class TestValidateMatrixGisaidHostsAndDuckAnatidae:
    def test_extra_host_aliases_semicolon_stripping_and_duck_to_anatidae(self, tmp_path: Path):
        taxa_dir = tmp_path / "Taxa"
        _write_mini_taxdump(taxa_dir)
        base_host_map = _write_tsv(tmp_path / "host_mapping.tsv", ["raw", "mapped"], [["Homo sapiens", "Homo sapiens"]])
        country_map = _write_tsv(tmp_path / "country_mapping.tsv", ["raw", "mapped"], [])
        country_file = _write_tsv(tmp_path / "m49.tsv", ["display_name", "m49_code"], [["United Kingdom", "826"]])

        gb_matrix = _write_tsv(
            tmp_path / "gB_matrix_merged.tsv",
            ["primary_accession", "host", "country", "collection_date", "host_validated", "host_taxa_id", "host_scientific_name"],
            [
                # Already validated GenBank row — must stay untouched in fill_missing_only mode
                ["GB1", "Mallard", "United Kingdom", "2024-01-01", "Anas platyrhynchos", "8839", "Anas platyrhynchos"],
                # Unvalidated GISAID rows
                ["GI_HUMAN", "Human", "United Kingdom", "2024-01-01", "", "", ""],
                ["GI_DUCK", "Duck", "United Kingdom", "2024-01-01", "", "", ""],
                ["GI_MALLARD_TAIL", "mallard; gender M; age hatch year", "United Kingdom", "2024-01-01", "", "", ""],
                ["GI_TURNSTONE", "ruddy turnstone; After hatch year", "United Kingdom", "2024-01-01", "", "", ""],
                ["GI_SWINE_TAIL", "swine; Juvenile", "United Kingdom", "2024-01-01", "", "", ""],
                ["GI_OTHER_AVIAN", "Other avian", "United Kingdom", "2024-01-01", "", "", ""],
                ["GI_OTHER_MAMMAL", "Other mammals", "United Kingdom", "2024-01-01", "", "", ""],
                ["GI_FELINE", "Feline", "United Kingdom", "2024-01-01", "", "", ""],
                ["GI_EAGLE", "Eagle", "United Kingdom", "2024-01-01", "", "", ""],
                ["GI_COW", "Dairy cow", "United Kingdom", "2024-01-01", "", "", ""],
                ["GI_ENV", "Environment", "United Kingdom", "2024-01-01", "", "", ""],
                ["GI_PETFOOD", "Pet Food", "United Kingdom", "2024-01-01", "", "", ""],
            ],
        )

        vm = ValidateMatrix(
            url="https://example.invalid/taxdump.tar.gz",
            taxa_path="Taxa",
            base_dir=str(tmp_path),
            output_dir="out",
            gb_matrix=str(gb_matrix),
            country_file=str(country_file),
            assets="",
            host_map=str(base_host_map),
            country_map=str(country_map),
            extra_host_map=str(INFLUENZA_HOST_ALIASES),
            fill_missing_only=True,
        )
        vm.validate_matrix()

        out_path = tmp_path / "out" / "gB_matrix_validated.tsv"
        with out_path.open("r", encoding="utf-8") as fh:
            rows = {r["primary_accession"]: r for r in csv.DictReader(fh, delimiter="\t")}

        assert rows["GB1"]["host_taxa_id"] == "8839"
        assert rows["GI_HUMAN"]["host_taxa_id"] == "9606"
        # Generic 'Duck' MUST map to Anatidae (8830), NOT Anas platyrhynchos (8839)!
        assert rows["GI_DUCK"]["host_taxa_id"] == "8830"
        assert rows["GI_DUCK"]["host_scientific_name"] == "Anatidae"
        # Semicolon tails stripped cleanly
        assert rows["GI_MALLARD_TAIL"]["host_taxa_id"] == "8839"
        assert rows["GI_TURNSTONE"]["host_taxa_id"] == "54971"
        assert rows["GI_SWINE_TAIL"]["host_taxa_id"] == "9823"
        # High-count influenza GISAID aliases resolved
        assert rows["GI_OTHER_AVIAN"]["host_taxa_id"] == "8782"
        assert rows["GI_OTHER_MAMMAL"]["host_taxa_id"] == "40674"
        assert rows["GI_FELINE"]["host_taxa_id"] == "9681"
        assert rows["GI_EAGLE"]["host_taxa_id"] == "30419"
        assert rows["GI_COW"]["host_taxa_id"] == "9913"
        # Non-biological strings stay unmapped
        assert rows["GI_ENV"]["host_taxa_id"] == ""
        assert rows["GI_PETFOOD"]["host_taxa_id"] == ""


class TestFourHostTaxaTablesInSqlite:
    def _build_minimal_inputs(self, work: Path, taxa_ids: list[str], db_name: str, names: Path, nodes: Path):
        work.mkdir(parents=True, exist_ok=True)
        gb_tsv = _write_tsv(
            work / "meta.tsv",
            ["primary_accession", "accession_type", "dataset_source", "host_taxa_id", "segment"],
            [[f"ACC_{i}_{tid}", "query", "genbank", tid, "1"] for i, tid in enumerate(taxa_ids)],
        )
        ht_dir = work / "HostTaxa"
        ht = HostTaxaTable(
            gb_matrix=str(gb_tsv),
            output_dir="HostTaxa",
            names_file=str(names),
            nodes_file=str(nodes),
            base_dir=str(work),
            host_output_file="Host_taxa.tsv",
            child_output_file="Host_taxa_children.tsv",
            lineage_output_file="Host_taxa_lineage.tsv",
            lineage_lookup_output_file="Host_taxa_lineage_lookup.tsv",
        )
        ht.write_tables()
        return gb_tsv, ht_dir

    def test_create_and_update_load_all_four_host_tables_preserving_synonyms(self, tmp_path: Path):
        names, nodes = _write_mini_taxdump(tmp_path / "Taxa")
        gb1, ht1 = self._build_minimal_inputs(tmp_path / "run1", ["9606", "8830"], "rabv_like", names, nodes)

        # Prepare minimal companion tables for CreateSqliteDB
        def _write_csv(path: Path, header, rows) -> Path:
            with path.open("w", newline="", encoding="utf-8") as fh:
                w = csv.writer(fh, delimiter=",")
                w.writerow(header)
                w.writerows(rows)
            return path

        def _prep_support(d: Path, accs: list[str]):
            m49_c = _write_csv(d / "m49_c.csv", ["m49_code", "country"], [["826", "UK"]])
            m49_i = _write_csv(d / "m49_i.csv", ["code", "name"], [])
            m49_r = _write_csv(d / "m49_r.csv", ["code", "name"], [])
            m49_s = _write_csv(d / "m49_s.csv", ["code", "name"], [])
            proj = _write_tsv(d / "proj.tsv", ["setting", "value"], [["version", "test"]])
            genes = _write_tsv(d / "genes.tsv", ["name", "display_name"], [["N", "N"]])
            feats = _write_tsv(d / "feats.tsv", ["accession", "product", "cds_start", "cds_end"], [[a, "N", "1", "12"] for a in accs])
            ins = _write_tsv(d / "ins.tsv", ["primary_accession", "pos"], [])
            aln = _write_tsv(d / "aln.tsv", ["primary_accession", "alignment_name", "segment", "alignment"], [[a, "aln1", "1", "ATGCGTACGTAA"] for a in accs])
            fa = d / "seqs.fa"
            fa.write_text("".join(f">{a}\nATGCGTACGTAA\n" for a in accs), encoding="utf-8")
            ref_list = d / "refs.txt"
            ref_list.write_text(f"{accs[0]}\tmaster\t1\n", encoding="utf-8")
            return m49_c, m49_i, m49_r, m49_s, proj, genes, feats, ins, aln, fa, ref_list

        d1 = tmp_path / "run1"
        m49_c, m49_i, m49_r, m49_s, proj, genes, feats, ins, aln, fa, ref_list = _prep_support(d1, ["ACC_0_9606", "ACC_1_8830"])
        c1 = CreateSqliteDB(
            meta_data=str(gb1),
            features=str(feats),
            pad_aln=str(aln),
            gene_info=str(genes),
            m49_countries=str(m49_c),
            m49_interm_region=str(m49_i),
            m49_regions=str(m49_r),
            m49_sub_regions=str(m49_s),
            proj_settings=str(proj),
            fasta_sequence_file=str(fa),
            insertions=str(ins),
            host_taxa_file=str(ht1 / "Host_taxa.tsv"),
            host_taxa_children_file=str(ht1 / "Host_taxa_children.tsv"),
            host_taxa_lineage_file=str(ht1 / "Host_taxa_lineage.tsv"),
            host_taxa_lineage_lookup_file=str(ht1 / "Host_taxa_lineage_lookup.tsv"),
            base_dir=str(d1),
            output_dir="SqliteDB",
            db_name="rabv_like",
            db_status="new db",
            reference_tsv=str(ref_list),
        )
        c1.create_db()

        db1 = d1 / "SqliteDB" / "rabv_like.db"
        conn = sqlite3.connect(db1)
        tables = {r[0] for r in conn.execute("SELECT name FROM sqlite_master WHERE type='table'")}
        for expected in ("host_taxa", "host_taxa_children", "host_taxa_lineage", "host_taxa_lineage_lookup"):
            assert expected in tables, f"Missing table {expected}"
        # 9606 has 3 rows in host_taxa (scientific, common 'human', common 'man')
        ht_rows_1 = conn.execute("SELECT taxa_id, name, name_type, common_name FROM host_taxa WHERE CAST(taxa_id AS TEXT)='9606'").fetchall()
        assert len(ht_rows_1) >= 2
        assert all(r[3] == "human" for r in ht_rows_1)
        assert conn.execute("SELECT COUNT(*) FROM host_taxa_children").fetchone()[0] > 0
        assert conn.execute("SELECT COUNT(*) FROM host_taxa_lineage WHERE CAST(taxa_id AS TEXT)='9606'").fetchone()[0] == 1
        assert conn.execute("SELECT COUNT(*) FROM host_taxa_lineage_lookup WHERE CAST(desc_taxa_id AS TEXT)='9606'").fetchone()[0] > 0
        conn.close()

        # Now run --update_db adding taxid 8839 (mallard, multiple synonym rows) + existing 9606!
        d2 = tmp_path / "run2"
        gb2, ht2 = self._build_minimal_inputs(d2, ["9606", "8839"], "rabv_like", names, nodes)
        m49_c2, m49_i2, m49_r2, m49_s2, proj2, genes2, feats2, ins2, aln2, fa2, ref_list2 = _prep_support(d2, ["ACC_0_9606", "ACC_1_8839"])
        c2 = CreateSqliteDB(
            meta_data=str(gb2),
            features=str(feats2),
            pad_aln=str(aln2),
            gene_info=str(genes2),
            m49_countries=str(m49_c2),
            m49_interm_region=str(m49_i2),
            m49_regions=str(m49_r2),
            m49_sub_regions=str(m49_s2),
            proj_settings=str(proj2),
            fasta_sequence_file=str(fa2),
            insertions=str(ins2),
            host_taxa_file=str(ht2 / "Host_taxa.tsv"),
            host_taxa_children_file=str(ht2 / "Host_taxa_children.tsv"),
            host_taxa_lineage_file=str(ht2 / "Host_taxa_lineage.tsv"),
            host_taxa_lineage_lookup_file=str(ht2 / "Host_taxa_lineage_lookup.tsv"),
            base_dir=str(d2),
            output_dir="SqliteDB",
            db_name="rabv_like",
            db_status="last updated",
            reference_tsv=str(ref_list2),
            update=True,
            update_db=str(db1),
        )
        c2.create_db()

        db2 = d2 / "SqliteDB" / "rabv_like.db"
        conn = sqlite3.connect(db2)
        # Synonyms for BOTH 9606 and 8839 MUST survive --update without duplication or dropping!
        ht_8839 = conn.execute("SELECT taxa_id, name, name_type, common_name FROM host_taxa WHERE CAST(taxa_id AS TEXT)='8839'").fetchall()
        assert len(ht_8839) >= 2, f"Expected multiple synonym rows for 8839, got {ht_8839}"
        ht_9606 = conn.execute("SELECT taxa_id, name, name_type FROM host_taxa WHERE CAST(taxa_id AS TEXT)='9606'").fetchall()
        assert len(ht_9606) == len(set(ht_9606)) and len(ht_9606) == len(ht_rows_1)
        assert conn.execute("SELECT COUNT(*) FROM host_taxa_lineage WHERE CAST(taxa_id AS TEXT) IN ('9606','8839')").fetchone()[0] == 2
        assert conn.execute("SELECT COUNT(*) FROM host_taxa_lineage_lookup WHERE CAST(desc_taxa_id AS TEXT)='8839'").fetchone()[0] > 0
        # No duplicate rows in children or lineage_lookup
        c_all = conn.execute("SELECT child_taxa_id, parent_taxa_id FROM host_taxa_children").fetchall()
        assert len(c_all) == len(set(c_all))
        ll_all = conn.execute("SELECT lineage_taxa_id, desc_taxa_id FROM host_taxa_lineage_lookup").fetchall()
        assert len(ll_all) == len(set(ll_all))
        conn.close()

    def test_workflow_glue_wires_all_four_host_tables_in_nf_and_sh(self):
        nf = (REPO_ROOT / "vgtk-init.nf").read_text(encoding="utf-8")
        sh = (REPO_ROOT / "vgtk-rabv.sh").read_text(encoding="utf-8")
        for fname in ("Host_taxa_children.tsv", "Host_taxa_lineage.tsv", "Host_taxa_lineage_lookup.tsv"):
            assert fname in nf, f"vgtk-init.nf missing {fname}"
            assert fname in sh, f"vgtk-rabv.sh missing {fname}"
        for flag in ("-htc", "-htl", "-htll"):
            assert flag in nf, f"vgtk-init.nf missing {flag}"
            assert flag in sh, f"vgtk-rabv.sh missing {flag}"
