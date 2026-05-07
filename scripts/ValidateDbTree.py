import argparse
import os
import sqlite3
from io import StringIO

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from Bio import Phylo

ACCESSION_COLUMNS = {
    "primary_accession",
    "accession",
    "accession_id",
    "accession_version",
    "sequence_id",
    "sequence_name",
    "seq_name",
    "sample",
    "query",
    "reference",
    "header",
}

CLUSTER_PLACEHOLDER_VALUES = {
    "na-see-tree",
    "na- see tree",
    "na - see tree",
}

NA_SET = {"", "NA", "Na", "N/A", "na", "n/a", "-", None}
HOST_TAXA_META_COLUMNS = ["host_taxa_id", "taxonomy_id", "host_tax_id"]
HOST_TAXA_TABLE_COLUMNS = ["taxa_id", "taxonomy_id", "host_taxa_id", "tax_id", "id"]
MUTATION_REQUIRED_COLUMNS = ["primary_accession", "mutation_id", "protein_name", "aa_position", "alt_residue"]
DRUG_RESISTANCE_REQUIRED_COLUMNS = ["primary_accession", "combination_id", "combination_status", "mutations_detected", "mutations_required"]
COMPACT_MUTATION_SUMMARY_COLUMNS = ["primary_accession", "relevant_mutations_present", "total_relevant_mutation_count"]
COMPLETED_SIGNATURE_COLUMNS = ["primary_accession", "signature_id", "signature_kind"]


def get_table_columns(conn, table_name):
    cursor = conn.cursor()
    cursor.execute(f"PRAGMA table_info({table_name})")
    return [row[1] for row in cursor.fetchall()]


def collect_table_values(conn, table_name, columns):
    values = set()
    cursor = conn.cursor()
    for col in columns:
        cursor.execute(f"SELECT DISTINCT {col} FROM {table_name}")
        values.update([r[0] for r in cursor.fetchall() if r[0]])
    return values


def find_cluster_column(columns):
    for col in columns:
        if col == "cluster" or col.startswith("cluster_"):
            return col
    return None


def is_cluster_placeholder(value):
    if value is None:
        return True
    norm = str(value).strip().lower()
    if not norm:
        return True
    return norm in CLUSTER_PLACEHOLDER_VALUES


def fetch_tree_newick(conn):
    rows = fetch_trees(conn)
    if not rows:
        return None, None, None

    ranked_rows = []
    for row in rows:
        tree = Phylo.read(StringIO(row["newick"]), "newick")
        terminal_count = len([tip for tip in tree.get_terminals() if tip.name])
        source_rank = 2
        if row["source"] == "usher":
            source_rank = 0
        elif row["source"] == "iqtree":
            source_rank = 1
        ranked_rows.append((source_rank, -terminal_count, row["name"] or "", row))

    best_row = sorted(ranked_rows, key=lambda item: (item[0], item[1], item[2]))[0][3]
    return best_row["name"], best_row["source"], best_row["newick"]


def fetch_trees(conn, source=None):
    cols = get_table_columns(conn, "trees")
    select_cols = ["name", "source", "newick"]
    if "segment_key" in cols:
        select_cols.append("segment_key")
    if "segment" in cols:
        select_cols.append("segment")

    cursor = conn.cursor()
    query = f"SELECT {', '.join(select_cols)} FROM trees WHERE newick IS NOT NULL AND length(newick) > 0"
    params = []
    if source:
        query += " AND source = ?"
        params.append(source)
    cursor.execute(query, params)

    rows = []
    for record in cursor.fetchall():
        row = {k: v for k, v in zip(select_cols, record)}
        row.setdefault("segment_key", None)
        row.setdefault("segment", None)
        rows.append(row)
    return rows


def fetch_table_column(conn, table_name, column_name):
    cursor = conn.cursor()
    cursor.execute(f"SELECT {column_name} FROM {table_name}")
    return [r[0] for r in cursor.fetchall()]


def table_exists(conn, table_name):
    cursor = conn.cursor()
    cursor.execute(
        "SELECT 1 FROM sqlite_master WHERE type='table' AND name=? LIMIT 1",
        (table_name,),
    )
    return cursor.fetchone() is not None


def fetch_distinct_values(conn, table, column, where_sql=None, params=()):
    cursor = conn.cursor()
    sql = f"SELECT DISTINCT {column} AS v FROM {table}"
    if where_sql:
        sql += f" WHERE {where_sql}"
    cursor.execute(sql, params)
    out = set()
    for row in cursor.fetchall():
        value = row[0]
        if value is None:
            continue
        text = str(value).strip()
        if text:
            out.add(text)
    return out


def fetch_count(conn, table, where_sql=None, params=()):
    cursor = conn.cursor()
    sql = f"SELECT COUNT(*) FROM {table}"
    if where_sql:
        sql += f" WHERE {where_sql}"
    cursor.execute(sql, params)
    row = cursor.fetchone()
    return int(row[0]) if row else 0


def compare_sets(expected, observed):
    missing = sorted(expected - observed)
    extra = sorted(observed - expected)
    return {
        "ok": len(missing) == 0 and len(extra) == 0,
        "expected_count": len(expected),
        "observed_count": len(observed),
        "missing": missing,
        "extra": extra,
    }


def find_first_present(columns, candidates):
    return next((col for col in candidates if col in columns), None)


def get_meta_nonexcluded_filter(meta_cols, exclusion_column, exclude_value):
    if exclusion_column not in meta_cols:
        return None, ()
    where_sql = f"({exclusion_column} IS NULL OR CAST({exclusion_column} AS TEXT) != ?)"
    return where_sql, (str(exclude_value),)


def get_meta_excluded_filter(meta_cols, exclusion_column, exclude_value):
    if exclusion_column not in meta_cols:
        return None, ()
    where_sql = f"NOT ({exclusion_column} IS NULL OR CAST({exclusion_column} AS TEXT) != ?)"
    return where_sql, (str(exclude_value),)


def get_expected_accessions(conn, accession_column="primary_accession", exclusion_column="exclusion_status", exclude_value="1"):
    if not table_exists(conn, "meta_data"):
        raise RuntimeError("Table 'meta_data' not found in DB.")

    meta_cols = get_table_columns(conn, "meta_data")
    if accession_column not in meta_cols:
        raise RuntimeError(
            f"meta_data does not contain accession column '{accession_column}'. Available columns: {meta_cols}"
        )

    where_sql, params = get_meta_nonexcluded_filter(meta_cols, exclusion_column, exclude_value)
    excluded_where_sql, excluded_params = get_meta_excluded_filter(meta_cols, exclusion_column, exclude_value)
    raw_expected = fetch_distinct_values(conn, "meta_data", accession_column, where_sql=where_sql, params=params)
    meta_excluded_accessions = set()
    if excluded_where_sql:
        meta_excluded_accessions = fetch_distinct_values(
            conn,
            "meta_data",
            accession_column,
            where_sql=excluded_where_sql,
            params=excluded_params,
        )
    legacy_excluded_accessions = set()
    if table_exists(conn, "excluded_accessions"):
        excluded_cols = get_table_columns(conn, "excluded_accessions")
        excluded_col = find_first_present(excluded_cols, [accession_column, "primary_accession", "accession"])
        if excluded_col:
            legacy_excluded_accessions = fetch_distinct_values(conn, "excluded_accessions", excluded_col)

    expected = raw_expected - legacy_excluded_accessions
    total_rows = fetch_count(conn, "meta_data")
    nonexcluded_rows = fetch_count(conn, "meta_data", where_sql=where_sql, params=params) if where_sql else total_rows
    accepted_accessions = set(expected)
    filtered_accessions = set(meta_excluded_accessions) | set(legacy_excluded_accessions)

    return {
        "expected_accessions": expected,
        "accepted_accessions": accepted_accessions,
        "filtered_accessions": filtered_accessions,
        "accepted_filtered_union": accepted_accessions | filtered_accessions,
        "meta_total_rows": total_rows,
        "meta_nonexcluded_rows": nonexcluded_rows,
        "meta_has_exclusion": exclusion_column in meta_cols,
        "excluded_accessions": filtered_accessions,
        "meta_columns": meta_cols,
    }


def validate_accession_table(conn, expected_accessions, table_name, candidate_columns, label):
    if not table_exists(conn, table_name):
        return {"ok": False, "error": f"Table '{table_name}' not found.", "table": table_name, "label": label}

    cols = get_table_columns(conn, table_name)
    col = find_first_present(cols, candidate_columns)
    if col is None:
        return {
            "ok": False,
            "error": f"Table '{table_name}' missing any of columns {candidate_columns}. Columns present: {cols}",
            "table": table_name,
            "label": label,
        }

    observed = fetch_distinct_values(conn, table_name, col)
    result = compare_sets(expected_accessions, observed)
    result.update({"table": table_name, "column": col, "label": label})
    return result


def validate_host_taxa(conn, meta_cols, where_sql=None, params=()):
    if not table_exists(conn, "host_taxa"):
        return {"ok": False, "error": "Table 'host_taxa' not found.", "table": "host_taxa"}

    meta_host_col = find_first_present(meta_cols, HOST_TAXA_META_COLUMNS)
    host_cols = get_table_columns(conn, "host_taxa")
    host_col = find_first_present(host_cols, HOST_TAXA_TABLE_COLUMNS)

    if meta_host_col is None and host_col is None:
        return {
            "ok": True,
            "skipped": True,
            "table": "host_taxa",
            "column": None,
            "reason": "No host taxonomy ID columns were present in either meta_data or host_taxa.",
        }

    if meta_host_col is None:
        return {"ok": False, "error": f"meta_data missing any host taxonomy column {HOST_TAXA_META_COLUMNS}.", "table": "host_taxa"}
    if host_col is None:
        return {
            "ok": False,
            "error": f"host_taxa missing any taxonomy identifier column {HOST_TAXA_TABLE_COLUMNS}. Columns present: {host_cols}",
            "table": "host_taxa",
        }

    expected = fetch_distinct_values(conn, "meta_data", meta_host_col, where_sql=where_sql, params=params)
    expected = {value for value in expected if value not in NA_SET}
    observed = fetch_distinct_values(conn, "host_taxa", host_col)

    result = compare_sets(expected, observed)
    result.update(
        {
            "table": "host_taxa",
            "column": host_col,
            "label": f"host_taxa vs meta_data.{meta_host_col}",
            "meta_column": meta_host_col,
        }
    )
    return result


def format_consistency_result(title, result, show_n=25):
    lines = [f"[{title}]"]
    if result.get("skipped"):
        lines.append("  OK: True (skipped)")
        lines.append(f"  Reason: {result.get('reason', 'not provided')}")
        return "\n".join(lines)
    if "error" in result:
        lines.append("  OK: False")
        lines.append(f"  ERROR: {result['error']}")
        return "\n".join(lines)

    lines.append(f"  OK: {result.get('ok')}")
    lines.append(f"  Table: {result.get('table')}  Column: {result.get('column')}")
    lines.append(f"  Expected: {result.get('expected_count')}  Observed: {result.get('observed_count')}")
    missing = result.get("missing", [])
    extra = result.get("extra", [])
    lines.append(f"  Missing: {len(missing)}  Extra: {len(extra)}")

    if missing:
        lines.append(f"  Missing examples (up to {show_n}):")
        for value in missing[:show_n]:
            lines.append(f"    - {value}")
    if extra:
        lines.append(f"  Extra examples (up to {show_n}):")
        for value in extra[:show_n]:
            lines.append(f"    - {value}")
    return "\n".join(lines)


def result_failed(result):
    return not result.get("ok", False) and not result.get("skipped", False)


def validate_mutation_tables(conn, expected_accessions):
    results = []

    has_legacy_mutation_tables = table_exists(conn, "sequence_mutations") or table_exists(conn, "sequence_drug_resistance")
    has_compact_mutation_tables = table_exists(conn, "sequence_relevant_mutation_summary") or table_exists(conn, "completed_signatures_only")

    if not has_legacy_mutation_tables and not has_compact_mutation_tables:
        return [
            {
                "title": "mutation tables",
                "ok": True,
                "skipped": True,
                "reason": "No mutation tables present in DB.",
            }
        ]

    if table_exists(conn, "sequence_relevant_mutation_summary"):
        summary_cols = get_table_columns(conn, "sequence_relevant_mutation_summary")
        missing_cols = [col for col in COMPACT_MUTATION_SUMMARY_COLUMNS if col not in summary_cols]
        if missing_cols:
            results.append(
                {
                    "title": "sequence_relevant_mutation_summary schema",
                    "ok": False,
                    "error": f"Missing required columns: {missing_cols}. Present columns: {summary_cols}",
                }
            )
        else:
            cursor = conn.cursor()
            cursor.execute(
                """
                SELECT COUNT(*)
                FROM sequence_relevant_mutation_summary
                WHERE primary_accession IS NULL OR TRIM(CAST(primary_accession AS TEXT)) = ''
                   OR relevant_mutations_present IS NULL OR TRIM(CAST(relevant_mutations_present AS TEXT)) = ''
                   OR total_relevant_mutation_count IS NULL
                """
            )
            blank_required = int(cursor.fetchone()[0])

            cursor.execute(
                """
                SELECT COUNT(*)
                FROM (
                    SELECT primary_accession, COUNT(*) AS c
                    FROM sequence_relevant_mutation_summary
                    GROUP BY 1
                    HAVING c > 1
                )
                """
            )
            duplicate_rows = int(cursor.fetchone()[0])

            observed_accessions = fetch_distinct_values(conn, "sequence_relevant_mutation_summary", "primary_accession")
            orphan_accessions = sorted(observed_accessions - expected_accessions)

            cursor.execute(
                """
                SELECT primary_accession, relevant_mutations_present, total_relevant_mutation_count
                FROM sequence_relevant_mutation_summary
                """
            )
            count_mismatches = 0
            for primary_accession, mutation_list, total_count in cursor.fetchall():
                mutation_ids = [value for value in str(mutation_list or '').split(';') if value]
                try:
                    parsed_total = int(total_count)
                except (TypeError, ValueError):
                    count_mismatches += 1
                    continue
                if parsed_total != len(mutation_ids) or parsed_total < 1:
                    count_mismatches += 1

            results.append(
                {
                    "title": "sequence_relevant_mutation_summary integrity",
                    "ok": blank_required == 0 and duplicate_rows == 0 and len(orphan_accessions) == 0 and count_mismatches == 0,
                    "table": "sequence_relevant_mutation_summary",
                    "blank_required_rows": blank_required,
                    "duplicate_keys": duplicate_rows,
                    "orphan_accessions": orphan_accessions,
                    "count_mismatches": count_mismatches,
                    "observed_accessions": len(observed_accessions),
                }
            )

    if table_exists(conn, "completed_signatures_only"):
        completed_cols = get_table_columns(conn, "completed_signatures_only")
        missing_cols = [col for col in COMPLETED_SIGNATURE_COLUMNS if col not in completed_cols]
        if missing_cols:
            results.append(
                {
                    "title": "completed_signatures_only schema",
                    "ok": False,
                    "error": f"Missing required columns: {missing_cols}. Present columns: {completed_cols}",
                }
            )
        else:
            cursor = conn.cursor()
            cursor.execute(
                """
                SELECT COUNT(*)
                FROM completed_signatures_only
                WHERE primary_accession IS NULL OR TRIM(CAST(primary_accession AS TEXT)) = ''
                   OR signature_id IS NULL OR TRIM(CAST(signature_id AS TEXT)) = ''
                   OR signature_kind IS NULL OR TRIM(CAST(signature_kind AS TEXT)) = ''
                """
            )
            blank_required = int(cursor.fetchone()[0])

            cursor.execute(
                """
                SELECT COUNT(*)
                FROM (
                    SELECT primary_accession, signature_id, COUNT(*) AS c
                    FROM completed_signatures_only
                    GROUP BY 1,2
                    HAVING c > 1
                )
                """
            )
            duplicate_rows = int(cursor.fetchone()[0])

            observed_accessions = fetch_distinct_values(conn, "completed_signatures_only", "primary_accession")
            orphan_accessions = sorted(observed_accessions - expected_accessions)

            cursor.execute(
                """
                SELECT COUNT(*)
                FROM completed_signatures_only
                WHERE signature_kind NOT IN ('single', 'combination')
                """
            )
            invalid_signature_kinds = int(cursor.fetchone()[0])

            results.append(
                {
                    "title": "completed_signatures_only integrity",
                    "ok": blank_required == 0 and duplicate_rows == 0 and len(orphan_accessions) == 0 and invalid_signature_kinds == 0,
                    "table": "completed_signatures_only",
                    "blank_required_rows": blank_required,
                    "duplicate_keys": duplicate_rows,
                    "orphan_accessions": orphan_accessions,
                    "invalid_signature_kinds": invalid_signature_kinds,
                    "observed_accessions": len(observed_accessions),
                }
            )

    if table_exists(conn, "sequence_mutations"):
        mut_cols = get_table_columns(conn, "sequence_mutations")
        missing_cols = [col for col in MUTATION_REQUIRED_COLUMNS if col not in mut_cols]
        if missing_cols:
            results.append(
                {
                    "title": "sequence_mutations schema",
                    "ok": False,
                    "error": f"Missing required columns: {missing_cols}. Present columns: {mut_cols}",
                }
            )
        else:
            cursor = conn.cursor()
            cursor.execute(
                """
                SELECT COUNT(*)
                FROM sequence_mutations
                WHERE primary_accession IS NULL OR TRIM(CAST(primary_accession AS TEXT)) = ''
                   OR mutation_id IS NULL OR TRIM(CAST(mutation_id AS TEXT)) = ''
                   OR protein_name IS NULL OR TRIM(CAST(protein_name AS TEXT)) = ''
                   OR aa_position IS NULL OR TRIM(CAST(aa_position AS TEXT)) = ''
                   OR alt_residue IS NULL OR TRIM(CAST(alt_residue AS TEXT)) = ''
                """
            )
            blank_required = int(cursor.fetchone()[0])

            observed_accessions = fetch_distinct_values(conn, "sequence_mutations", "primary_accession")
            orphan_accessions = sorted(observed_accessions - expected_accessions)

            combination_col = "combination_id" if "combination_id" in mut_cols else None
            if combination_col:
                duplicate_sql = """
                    SELECT COUNT(*)
                    FROM (
                        SELECT primary_accession, mutation_id, protein_name,
                               COALESCE(TRIM(CAST(segment AS TEXT)), ''),
                               COALESCE(TRIM(CAST(aa_position AS TEXT)), ''),
                               alt_residue,
                               COALESCE(TRIM(CAST(combination_id AS TEXT)), ''),
                               COUNT(*) AS c
                        FROM sequence_mutations
                        GROUP BY 1,2,3,4,5,6,7
                        HAVING c > 1
                    )
                """
            else:
                duplicate_sql = """
                    SELECT COUNT(*)
                    FROM (
                        SELECT primary_accession, mutation_id, protein_name,
                               COALESCE(TRIM(CAST(segment AS TEXT)), ''),
                               COALESCE(TRIM(CAST(aa_position AS TEXT)), ''),
                               alt_residue,
                               COUNT(*) AS c
                        FROM sequence_mutations
                        GROUP BY 1,2,3,4,5,6
                        HAVING c > 1
                    )
                """
            cursor.execute(duplicate_sql)
            duplicate_rows = int(cursor.fetchone()[0])

            results.append(
                {
                    "title": "sequence_mutations integrity",
                    "ok": blank_required == 0 and duplicate_rows == 0 and len(orphan_accessions) == 0,
                    "table": "sequence_mutations",
                    "blank_required_rows": blank_required,
                    "duplicate_keys": duplicate_rows,
                    "orphan_accessions": orphan_accessions,
                    "observed_accessions": len(observed_accessions),
                }
            )
    elif table_exists(conn, "sequence_drug_resistance"):
        results.append(
            {
                "title": "sequence_mutations table",
                "ok": False,
                "error": "sequence_drug_resistance exists but sequence_mutations table is missing.",
            }
        )

    if table_exists(conn, "sequence_drug_resistance"):
        dr_cols = get_table_columns(conn, "sequence_drug_resistance")
        missing_cols = [col for col in DRUG_RESISTANCE_REQUIRED_COLUMNS if col not in dr_cols]
        if missing_cols:
            results.append(
                {
                    "title": "sequence_drug_resistance schema",
                    "ok": False,
                    "error": f"Missing required columns: {missing_cols}. Present columns: {dr_cols}",
                }
            )
        else:
            cursor = conn.cursor()
            cursor.execute(
                """
                SELECT COUNT(*)
                FROM sequence_drug_resistance
                WHERE primary_accession IS NULL OR TRIM(CAST(primary_accession AS TEXT)) = ''
                   OR combination_id IS NULL OR TRIM(CAST(combination_id AS TEXT)) = ''
                   OR combination_status IS NULL OR TRIM(CAST(combination_status AS TEXT)) = ''
                   OR mutations_detected IS NULL
                   OR mutations_required IS NULL
                """
            )
            blank_required = int(cursor.fetchone()[0])

            cursor.execute(
                """
                SELECT COUNT(*)
                FROM (
                    SELECT primary_accession, combination_id, COUNT(*) AS c
                    FROM sequence_drug_resistance
                    GROUP BY 1,2
                    HAVING c > 1
                )
                """
            )
            duplicate_rows = int(cursor.fetchone()[0])

            observed_accessions = fetch_distinct_values(conn, "sequence_drug_resistance", "primary_accession")
            orphan_accessions = sorted(observed_accessions - expected_accessions)

            cursor.execute(
                """
                SELECT COUNT(*)
                FROM sequence_drug_resistance
                WHERE CAST(mutations_detected AS INTEGER) < 1
                   OR CAST(mutations_required AS INTEGER) < 1
                """
            )
            non_positive_counts = int(cursor.fetchone()[0])

            cursor.execute(
                """
                SELECT COUNT(*)
                FROM sequence_drug_resistance
                WHERE combination_status NOT IN ('complete', 'partial')
                """
            )
            invalid_status_values = int(cursor.fetchone()[0])

            combo_key_mismatches = 0
            detected_count_mismatches = 0
            status_mismatches = 0
            if table_exists(conn, "sequence_mutations"):
                cursor.execute(
                    """
                    WITH mutation_combo_counts AS (
                        SELECT primary_accession, combination_id, COUNT(*) AS mutation_count
                        FROM sequence_mutations
                        WHERE combination_id IS NOT NULL AND TRIM(CAST(combination_id AS TEXT)) != ''
                        GROUP BY 1,2
                    ),
                    resistance_keys AS (
                        SELECT primary_accession, combination_id, mutations_detected, mutations_required, combination_status
                        FROM sequence_drug_resistance
                    )
                    SELECT COUNT(*)
                    FROM (
                        SELECT primary_accession, combination_id FROM mutation_combo_counts
                        EXCEPT
                        SELECT primary_accession, combination_id FROM resistance_keys
                    )
                    UNION ALL
                    SELECT COUNT(*)
                    FROM (
                        SELECT primary_accession, combination_id FROM resistance_keys
                        EXCEPT
                        SELECT primary_accession, combination_id FROM mutation_combo_counts
                    )
                    """
                )
                combo_key_mismatch_rows = cursor.fetchall()
                combo_key_mismatches = sum(int(row[0]) for row in combo_key_mismatch_rows)

                cursor.execute(
                    """
                    WITH mutation_combo_counts AS (
                        SELECT primary_accession, combination_id, COUNT(*) AS mutation_count
                        FROM sequence_mutations
                        WHERE combination_id IS NOT NULL AND TRIM(CAST(combination_id AS TEXT)) != ''
                        GROUP BY 1,2
                    )
                    SELECT COUNT(*)
                    FROM sequence_drug_resistance sdr
                    JOIN mutation_combo_counts mcc
                      ON sdr.primary_accession = mcc.primary_accession
                     AND sdr.combination_id = mcc.combination_id
                    WHERE CAST(sdr.mutations_detected AS INTEGER) != mcc.mutation_count
                    """
                )
                detected_count_mismatches = int(cursor.fetchone()[0])

                cursor.execute(
                    """
                    SELECT COUNT(*)
                    FROM sequence_drug_resistance
                    WHERE combination_status != CASE
                        WHEN CAST(mutations_detected AS INTEGER) >= CAST(mutations_required AS INTEGER) THEN 'complete'
                        ELSE 'partial'
                    END
                    """
                )
                status_mismatches = int(cursor.fetchone()[0])

            results.append(
                {
                    "title": "sequence_drug_resistance integrity",
                    "ok": (
                        blank_required == 0
                        and duplicate_rows == 0
                        and len(orphan_accessions) == 0
                        and non_positive_counts == 0
                        and invalid_status_values == 0
                        and combo_key_mismatches == 0
                        and detected_count_mismatches == 0
                        and status_mismatches == 0
                    ),
                    "table": "sequence_drug_resistance",
                    "blank_required_rows": blank_required,
                    "duplicate_keys": duplicate_rows,
                    "orphan_accessions": orphan_accessions,
                    "non_positive_counts": non_positive_counts,
                    "invalid_status_values": invalid_status_values,
                    "combo_key_mismatches": combo_key_mismatches,
                    "detected_count_mismatches": detected_count_mismatches,
                    "status_mismatches": status_mismatches,
                    "observed_accessions": len(observed_accessions),
                }
            )

    return results


def _try_parse_int(value):
    text = str(value).strip() if value is not None else ""
    if text in NA_SET:
        return None
    try:
        return int(text)
    except ValueError:
        return None


def validate_feature_projection_integrity(conn):
    title = "feature projection integrity"
    if not table_exists(conn, "features"):
        return {
            "title": title,
            "ok": True,
            "skipped": True,
            "reason": "Table 'features' not present in DB.",
        }

    feature_cols = get_table_columns(conn, "features")
    required_cols = ["accession", "master_ref_accession", "product", "aln_start", "aln_end", "cds_start", "cds_end"]
    missing_cols = [col for col in required_cols if col not in feature_cols]
    if missing_cols:
        return {
            "title": title,
            "ok": True,
            "skipped": True,
            "reason": (
                "features table does not contain projection-detail columns needed for this check: "
                f"missing {missing_cols}"
            ),
        }

    cursor = conn.cursor()
    cursor.execute(
        "SELECT accession, master_ref_accession, product, aln_start, aln_end, cds_start, cds_end FROM features"
    )
    feature_rows = cursor.fetchall()

    master_spans = {}
    for accession, _, product, _, _, cds_start, cds_end in feature_rows:
        accession_text = str(accession).strip() if accession is not None else ""
        product_text = str(product).strip() if product is not None else ""
        start = _try_parse_int(cds_start)
        end = _try_parse_int(cds_end)
        if not accession_text or not product_text or start is None or end is None:
            continue
        master_spans.setdefault((accession_text, product_text), []).append((start, end))

    offending_rows = []
    unresolved_master_rows = 0
    invalid_aln_rows = 0
    for accession, master_ref_accession, product, aln_start, aln_end, _, _ in feature_rows:
        accession_text = str(accession).strip() if accession is not None else ""
        master_text = str(master_ref_accession).strip() if master_ref_accession is not None else ""
        product_text = str(product).strip() if product is not None else ""
        if not accession_text or not master_text or not product_text or accession_text == master_text:
            continue

        aln_start_int = _try_parse_int(aln_start)
        aln_end_int = _try_parse_int(aln_end)
        if aln_start_int is None or aln_end_int is None:
            invalid_aln_rows += 1
            continue

        master_product_spans = master_spans.get((master_text, product_text))
        if not master_product_spans:
            unresolved_master_rows += 1
            continue

        overlaps_master_feature = any(
            aln_end_int >= master_start and aln_start_int <= master_end
            for master_start, master_end in master_product_spans
        )
        if not overlaps_master_feature:
            offending_rows.append(
                {
                    "accession": accession_text,
                    "master_ref_accession": master_text,
                    "product": product_text,
                    "aln_start": str(aln_start),
                    "aln_end": str(aln_end),
                }
            )

    return {
        "title": title,
        "ok": len(offending_rows) == 0 and invalid_aln_rows == 0,
        "table": "features",
        "offending_rows": len(offending_rows),
        "invalid_aln_rows": invalid_aln_rows,
        "unresolved_master_rows": unresolved_master_rows,
        "examples": offending_rows[:25],
    }


def format_mutation_integrity_result(result, show_n=25):
    lines = [f"[{result['title']}]"]
    if result.get("skipped"):
        lines.append("  OK: True (skipped)")
        lines.append(f"  Reason: {result.get('reason', 'not provided')}")
        return "\n".join(lines)
    if "error" in result:
        lines.append("  OK: False")
        lines.append(f"  ERROR: {result['error']}")
        return "\n".join(lines)

    lines.append(f"  OK: {result.get('ok')}")
    if result.get("table"):
        lines.append(f"  Table: {result['table']}")
    for key in [
        "blank_required_rows",
        "duplicate_keys",
        "non_positive_counts",
        "invalid_status_values",
        "combo_key_mismatches",
        "detected_count_mismatches",
        "status_mismatches",
        "observed_accessions",
    ]:
        if key in result:
            lines.append(f"  {key}: {result[key]}")
    orphan_accessions = result.get("orphan_accessions", [])
    if orphan_accessions:
        lines.append(f"  Orphan accession examples (up to {show_n}):")
        for value in orphan_accessions[:show_n]:
            lines.append(f"    - {value}")
    return "\n".join(lines)


def format_feature_integrity_result(result, show_n=25):
    lines = [f"[{result['title']}]"]
    if result.get("skipped"):
        lines.append("  OK: True (skipped)")
        lines.append(f"  Reason: {result.get('reason', 'not provided')}")
        return "\n".join(lines)
    if "error" in result:
        lines.append("  OK: False")
        lines.append(f"  ERROR: {result['error']}")
        return "\n".join(lines)

    lines.append(f"  OK: {result.get('ok')}")
    if result.get("table"):
        lines.append(f"  Table: {result['table']}")
    lines.append(f"  offending_rows: {result.get('offending_rows', 0)}")
    lines.append(f"  invalid_aln_rows: {result.get('invalid_aln_rows', 0)}")
    lines.append(f"  unresolved_master_rows: {result.get('unresolved_master_rows', 0)}")
    examples = result.get("examples", [])
    if examples:
        lines.append(f"  Offending examples (up to {show_n}):")
        for example in examples[:show_n]:
            lines.append(
                "    - accession={accession} product={product} master_ref_accession={master_ref_accession} aln_start={aln_start} aln_end={aln_end}".format(
                    **example
                )
            )
    return "\n".join(lines)


def main(argv=None):
    parser = argparse.ArgumentParser(description="Validate SQLite DB contents against tree coverage and table consistency")
    parser.add_argument("--db", required=True, help="Path to SQLite DB")
    parser.add_argument("--outdir", required=True, help="Output directory for report and plot")
    parser.add_argument(
        "--test-mode",
        action="store_true",
        help="Relax strict failures for known test-only edge cases (e.g., no placeable queries)",
    )
    parser.add_argument(
        "--expect-segment-trees",
        action="store_true",
        help="In segmented test mode, require at least one USHER tree per segment in meta_data",
    )
    parser.add_argument(
        "--check-update-integrity",
        action="store_true",
        help="Validate update audit tables and update-mode integrity constraints",
    )
    parser.add_argument(
        "--allow-no-trees",
        action="store_true",
        help="Allow validation to pass when the DB intentionally contains no trees",
    )
    parser.add_argument(
        "--accession-column",
        default="primary_accession",
        help="Accession column in meta_data used as the expected accession set",
    )
    parser.add_argument(
        "--exclusion-column",
        default="exclusion_status",
        help="Exclusion column in meta_data used to filter the expected accession set",
    )
    parser.add_argument(
        "--exclude-value",
        default="1",
        help="Value in the exclusion column that marks a row as excluded",
    )
    args = parser.parse_args(argv)

    os.makedirs(args.outdir, exist_ok=True)

    conn = sqlite3.connect(args.db)
    try:
        meta_info = get_expected_accessions(
            conn,
            accession_column=args.accession_column,
            exclusion_column=args.exclusion_column,
            exclude_value=args.exclude_value,
        )
        expected_accessions = meta_info["expected_accessions"]
        meta_columns = meta_info["meta_columns"]
        where_sql, where_params = get_meta_nonexcluded_filter(meta_columns, args.exclusion_column, args.exclude_value)

        tree_name, tree_source, newick = fetch_tree_newick(conn)
        if not newick and not args.allow_no_trees:
            raise SystemExit("No tree found in DB (trees table is empty)")

        tree_available = bool(newick)
        tree = None
        tree_terminals = set()
        if tree_available:
            tree = Phylo.read(StringIO(newick), "newick")
            tree_terminals = {t.name for t in tree.get_terminals() if t.name}

        meta_accessions = [x for x in fetch_table_column(conn, "meta_data", args.accession_column) if x]
        seq_headers = [x for x in fetch_table_column(conn, "sequences", "header") if x]

        meta_set = {str(x).strip() for x in meta_accessions if str(x).strip()}
        seq_set = {str(x).strip() for x in seq_headers if str(x).strip()}
        excluded_accessions = set(meta_info["excluded_accessions"])
        accepted_accessions = set(meta_info["accepted_accessions"])
        filtered_accessions = set(meta_info["filtered_accessions"])
        accepted_filtered_union = set(meta_info["accepted_filtered_union"])
        expected_meta_set = set(expected_accessions)
        missing_in_sequences = sorted(expected_meta_set - seq_set)
        missing_in_meta = sorted(tree_terminals - meta_set) if tree_available else []

        cluster_col = find_cluster_column(meta_columns)
        centroid_set = set()
        missing_centroids_in_tree = []
        extra_in_tree = []
        if cluster_col:
            cursor = conn.cursor()
            cursor.execute(f"SELECT {args.accession_column}, {cluster_col} FROM meta_data")
            cluster_rows = cursor.fetchall()
            centroid_set = {
                str(cluster_val).strip()
                for _, cluster_val in cluster_rows
                if not is_cluster_placeholder(cluster_val)
            }
            if tree_available:
                missing_centroids_in_tree = sorted(centroid_set - tree_terminals)
                extra_in_tree = sorted(tree_terminals - centroid_set)

        missing_in_tree = sorted(expected_meta_set - tree_terminals) if tree_available else []

        consistency_results = {
            "sequence_alignment vs meta_data": validate_accession_table(
                conn,
                expected_meta_set,
                "sequence_alignment",
                ["primary_accession", "accession", "sequence_id"],
                "sequence_alignment vs meta_data",
            ),
            "features vs meta_data": validate_accession_table(
                conn,
                expected_meta_set,
                "features",
                ["accession", "primary_accession"],
                "features vs meta_data",
            ),
            "host_taxa vs meta_data": validate_host_taxa(
                conn,
                meta_columns,
                where_sql=where_sql,
                params=where_params,
            ),
        }
        mutation_integrity_results = validate_mutation_tables(conn, expected_meta_set)
        feature_integrity_result = validate_feature_projection_integrity(conn) if args.check_update_integrity else {
            "title": "feature projection integrity",
            "ok": True,
            "skipped": True,
            "reason": "Update integrity checks not requested.",
        }

        cursor = conn.cursor()
        cursor.execute("SELECT COUNT(*) FROM meta_data")
        meta_count = cursor.fetchone()[0]
        cursor.execute("SELECT COUNT(*) FROM sequences")
        seq_count = cursor.fetchone()[0]
        cursor.execute("SELECT COUNT(*) FROM trees WHERE source='usher'")
        usher_tree_count = cursor.fetchone()[0]

        segment_count = None
        if "segment" in meta_columns:
            segment_values = [
                x for x in fetch_table_column(conn, "meta_data", "segment")
                if x is not None and str(x).strip() != ""
            ]
            segment_count = len(set(segment_values))

        segmented_validation_ok = False
        if tree_available and args.expect_segment_trees and segment_count is not None and segment_count > 1:
            usher_trees = fetch_trees(conn, source="usher")

            cursor.execute(f"SELECT {args.accession_column}, segment FROM meta_data")
            meta_rows = cursor.fetchall()

            accession_to_segment = {}
            accessions_by_segment = {}
            segments_in_meta = set()
            for acc, seg in meta_rows:
                if acc:
                    accession_to_segment[str(acc)] = str(seg).strip() if seg is not None else None
                if seg is None:
                    continue
                seg_str = str(seg).strip()
                if not seg_str:
                    continue
                segments_in_meta.add(seg_str)
                if acc and str(acc) not in excluded_accessions:
                    accessions_by_segment.setdefault(seg_str, set()).add(str(acc))

            trees_by_segment = {}
            for tr in usher_trees:
                seg = tr.get("segment")
                seg_key = tr.get("segment_key")
                if (seg is None or str(seg).strip() == "") and seg_key:
                    seg = accession_to_segment.get(str(seg_key))
                if seg is None or str(seg).strip() == "":
                    continue
                seg_str = str(seg).strip()
                t = Phylo.read(StringIO(tr["newick"]), "newick")
                terms = {x.name for x in t.get_terminals() if x.name}
                trees_by_segment.setdefault(seg_str, set()).update(terms)

            missing_segment_trees = sorted(s for s in segments_in_meta if s not in trees_by_segment)
            per_segment_missing = {}
            per_segment_extra = {}
            for seg in sorted(segments_in_meta):
                expected_segment_accessions = accessions_by_segment.get(seg, set())
                tree_terms = trees_by_segment.get(seg, set())
                missing_segment_accessions = sorted(expected_segment_accessions - tree_terms)
                if missing_segment_accessions:
                    per_segment_missing[seg] = missing_segment_accessions
                extra_segment_accessions = sorted(tree_terms - expected_segment_accessions)
                if extra_segment_accessions:
                    per_segment_extra[seg] = extra_segment_accessions

            print(f"[info] Segmented validation: meta segments={len(segments_in_meta)} tree segments={len(trees_by_segment)}")
            if missing_segment_trees:
                print(f"[info] Missing tree segments: {', '.join(missing_segment_trees[:10])}")
            if per_segment_missing:
                first_seg = sorted(per_segment_missing.keys())[0]
                print(f"[info] Segment {first_seg} missing accessions in tree (first 10): {', '.join(per_segment_missing[first_seg][:10])}")
            if per_segment_extra:
                first_seg = sorted(per_segment_extra.keys())[0]
                print(f"[info] Segment {first_seg} extra accessions in tree (first 10): {', '.join(per_segment_extra[first_seg][:10])}")

            if missing_segment_trees:
                raise SystemExit(
                    "Validation failed: segmented run expects at least one UShER tree per segment "
                    f"(missing segments={', '.join(missing_segment_trees[:10])})"
                )
            if per_segment_missing and not args.test_mode:
                raise SystemExit("Validation failed: per-segment UShER tree is missing unfiltered accessions")
            if per_segment_missing and args.test_mode:
                print("[info] Test mode: allowing per-segment UShER trees to be subsampled relative to meta_data accessions")

            segmented_validation_ok = True

        report_path = os.path.join(args.outdir, "db_tree_validation.txt")
        missing_tree_path = os.path.join(args.outdir, "missing_in_tree.txt")
        missing_seq_path = os.path.join(args.outdir, "missing_in_sequences.txt")
        missing_meta_path = os.path.join(args.outdir, "missing_in_meta.txt")
        missing_centroids_path = os.path.join(args.outdir, "missing_centroids_in_tree.txt")
        extra_in_tree_path = os.path.join(args.outdir, "extra_in_tree.txt")

        with open(report_path, "w", encoding="utf-8") as report:
            report.write("=== ValidateDB Summary ===\n")
            report.write(f"DB: {args.db}\n")
            report.write(f"meta_data total rows: {meta_info['meta_total_rows']}\n")
            if meta_info["meta_has_exclusion"]:
                report.write(
                    f"meta_data non-excluded rows: {meta_info['meta_nonexcluded_rows']} "
                    f"(excluding {args.exclusion_column} == {args.exclude_value})\n"
                )
            else:
                report.write("meta_data exclusion column not found; using all rows.\n")
            report.write(f"Accepted distinct accessions: {len(accepted_accessions)}\n")
            report.write(f"Filtered distinct accessions: {len(filtered_accessions)}\n")
            report.write(f"Accepted + filtered distinct accessions: {len(accepted_filtered_union)}\n")
            report.write(f"Excluded accessions considered: {len(excluded_accessions)}\n")
            report.write(f"Expected distinct accessions (from meta_data.{args.accession_column}): {len(expected_meta_set)}\n")
            report.write("\n")

            for title, result in consistency_results.items():
                report.write(format_consistency_result(title, result))
                report.write("\n\n")

            report.write("=== Tree Validation Summary ===\n")
            if tree_available:
                report.write(f"Tree status: present\n")
                report.write(f"Tree source: {tree_source or tree_name}\n")
                report.write(f"Tree name: {tree_name or ''}\n")
            else:
                report.write("Tree status: not present (allowed by --allow-no-trees)\n")
            report.write(f"Meta rows: {meta_count}\n")
            report.write(f"Sequence rows: {seq_count}\n")
            report.write(f"Tree terminals: {len(tree_terminals)}\n")
            report.write(f"UShER tree rows: {usher_tree_count}\n")
            if segment_count is not None:
                report.write(f"Segment count in meta_data: {segment_count}\n")
            if cluster_col:
                report.write(f"Cluster column: {cluster_col}\n")
                report.write(f"Centroid count: {len(centroid_set)}\n")
                report.write(f"Missing centroids in tree: {len(missing_centroids_in_tree)}\n")
                report.write(f"Extra nodes in tree: {len(extra_in_tree)}\n")
            report.write(f"Missing in tree (meta_data -> tree): {len(missing_in_tree)}\n")
            report.write(f"Missing in sequences (meta_data -> sequences): {len(missing_in_sequences)}\n")
            report.write(f"Missing in meta_data (tree -> meta_data): {len(missing_in_meta)}\n")

            if missing_in_tree:
                report.write("\nMissing in tree (first 50):\n")
                report.write("\n".join(missing_in_tree[:50]) + "\n")
            if missing_in_sequences:
                report.write("\nMissing in sequences (first 50):\n")
                report.write("\n".join(missing_in_sequences[:50]) + "\n")
            if missing_in_meta:
                report.write("\nMissing in meta_data (first 50):\n")
                report.write("\n".join(missing_in_meta[:50]) + "\n")

            if cluster_col and missing_centroids_in_tree:
                report.write("\nMissing centroids in tree (first 50):\n")
                report.write("\n".join(missing_centroids_in_tree[:50]) + "\n")
            if cluster_col and extra_in_tree:
                report.write("\nExtra nodes in tree (first 50):\n")
                report.write("\n".join(extra_in_tree[:50]) + "\n")

            report.write("\nTable-level accession coverage:\n")
            table_list = ["meta_data", "sequences", "sequence_alignment", "features", "insertions", "host_taxa"]
            if table_exists(conn, "sequence_mutations"):
                table_list.append("sequence_mutations")
            if table_exists(conn, "sequence_drug_resistance"):
                table_list.append("sequence_drug_resistance")
            for table_name in table_list:
                if not table_exists(conn, table_name):
                    report.write(f"- {table_name}: table missing\n")
                    continue
                cols = get_table_columns(conn, table_name)
                acc_cols = [c for c in cols if c in ACCESSION_COLUMNS]
                if not acc_cols:
                    report.write(f"- {table_name}: no accession-like columns found\n")
                    continue
                values = collect_table_values(conn, table_name, acc_cols)
                if tree_available:
                    missing_from_table = sorted(tree_terminals - values)
                    report.write(
                        f"- {table_name}: columns={','.join(acc_cols)} values={len(values)} missing_from_table={len(missing_from_table)}\n"
                    )
                else:
                    report.write(
                        f"- {table_name}: columns={','.join(acc_cols)} values={len(values)} tree_coverage_skipped=1\n"
                    )

            if args.check_update_integrity:
                report.write("\nUpdate integrity checks:\n")
                if table_exists(conn, "update_batches") and table_exists(conn, "update_table_deltas"):
                    cursor.execute("SELECT batch_id, mode, started_at, finished_at FROM update_batches ORDER BY rowid DESC LIMIT 1")
                    last_batch = cursor.fetchone()
                    if last_batch:
                        batch_id, mode, started_at, finished_at = last_batch
                        report.write(f"- last_batch_id: {batch_id}\n")
                        report.write(f"- last_batch_mode: {mode}\n")
                        report.write(f"- last_batch_started: {started_at}\n")
                        report.write(f"- last_batch_finished: {finished_at}\n")
                        cursor.execute(
                            "SELECT table_name, before_count, after_count, delta FROM update_table_deltas WHERE batch_id=? ORDER BY table_name",
                            (batch_id,),
                        )
                        rows = cursor.fetchall()
                        for table_name, before_count, after_count, delta in rows:
                            report.write(f"  - {table_name}: before={before_count} after={after_count} delta={delta}\n")
                else:
                    report.write("- update audit tables missing (update_batches/update_table_deltas)\n")

                if table_exists(conn, "features") and table_exists(conn, "meta_data"):
                    feature_cols = get_table_columns(conn, "features")
                    if "segment" in feature_cols and "segment" in meta_columns and "accession" in feature_cols:
                        cursor.execute(
                            """
                            SELECT COUNT(*)
                            FROM features f
                            JOIN meta_data m ON m.primary_accession = f.accession
                            WHERE COALESCE(TRIM(f.segment), '') != COALESCE(TRIM(m.segment), '')
                            """
                        )
                        seg_mismatch_count = cursor.fetchone()[0]
                        report.write(f"- feature_segment_mismatch: {seg_mismatch_count}\n")
                    report.write(format_feature_integrity_result(feature_integrity_result))
                    report.write("\n")

            if mutation_integrity_results:
                report.write("\nMutation integrity checks:\n")
                for result in mutation_integrity_results:
                    report.write(format_mutation_integrity_result(result))
                    report.write("\n\n")

        if missing_in_tree:
            with open(missing_tree_path, "w", encoding="utf-8") as handle:
                handle.write("\n".join(missing_in_tree) + "\n")
        if missing_in_sequences:
            with open(missing_seq_path, "w", encoding="utf-8") as handle:
                handle.write("\n".join(missing_in_sequences) + "\n")
        if missing_in_meta:
            with open(missing_meta_path, "w", encoding="utf-8") as handle:
                handle.write("\n".join(missing_in_meta) + "\n")

        if cluster_col and missing_centroids_in_tree:
            with open(missing_centroids_path, "w", encoding="utf-8") as handle:
                handle.write("\n".join(missing_centroids_in_tree) + "\n")
        if cluster_col and extra_in_tree:
            with open(extra_in_tree_path, "w", encoding="utf-8") as handle:
                handle.write("\n".join(extra_in_tree) + "\n")

        if tree_available:
            print(f"[info] Tree source: {tree_source or tree_name}")
            print(f"[info] Tree name: {tree_name or ''}")
        else:
            print("[info] Tree status: not present (allowed by --allow-no-trees)")
        print(f"[info] Accepted distinct accessions: {len(accepted_accessions)}")
        print(f"[info] Filtered distinct accessions: {len(filtered_accessions)}")
        print(f"[info] Accepted + filtered distinct accessions: {len(accepted_filtered_union)}")
        print(f"[info] Meta rows: {meta_count}, Sequence rows: {seq_count}, Tree terminals: {len(tree_terminals)}")
        print(f"[info] UShER tree rows: {usher_tree_count}")
        if segment_count is not None:
            print(f"[info] Segment count in meta_data: {segment_count}")
        for title, result in consistency_results.items():
            if result.get("skipped"):
                print(f"[info] {title}: skipped ({result.get('reason')})")
            elif "error" in result:
                print(f"[info] {title}: ERROR: {result['error']}")
            else:
                print(
                    f"[info] {title}: ok={result['ok']} missing={len(result.get('missing', []))} extra={len(result.get('extra', []))}"
                )
        for result in mutation_integrity_results:
            if result.get("skipped"):
                print(f"[info] {result['title']}: skipped ({result.get('reason')})")
            elif "error" in result:
                print(f"[info] {result['title']}: ERROR: {result['error']}")
            else:
                print(
                    f"[info] {result['title']}: ok={result['ok']} "
                    f"blank_required_rows={result.get('blank_required_rows', 0)} "
                    f"duplicate_keys={result.get('duplicate_keys', 0)} "
                    f"orphan_accessions={len(result.get('orphan_accessions', []))}"
                )
        if args.check_update_integrity:
            if feature_integrity_result.get("skipped"):
                print(f"[info] {feature_integrity_result['title']}: skipped ({feature_integrity_result.get('reason')})")
            elif "error" in feature_integrity_result:
                print(f"[info] {feature_integrity_result['title']}: ERROR: {feature_integrity_result['error']}")
            else:
                print(
                    f"[info] {feature_integrity_result['title']}: ok={feature_integrity_result.get('ok')} "
                    f"offending_rows={feature_integrity_result.get('offending_rows', 0)} "
                    f"invalid_aln_rows={feature_integrity_result.get('invalid_aln_rows', 0)} "
                    f"unresolved_master_rows={feature_integrity_result.get('unresolved_master_rows', 0)}"
                )
        if cluster_col:
            print(f"[info] Cluster column: {cluster_col}, Centroid count: {len(centroid_set)}")
            print(f"[info] Missing centroids in tree: {len(missing_centroids_in_tree)}")
            print(f"[info] Extra nodes in tree: {len(extra_in_tree)}")
        print(f"[info] Missing in tree: {len(missing_in_tree)}")
        print(f"[info] Missing in sequences: {len(missing_in_sequences)}")
        print(f"[info] Missing in meta_data: {len(missing_in_meta)}")
        if missing_in_tree:
            print(f"[info] Missing in tree (first 10): {', '.join(missing_in_tree[:10])}")
        if missing_in_sequences:
            print(f"[info] Missing in sequences (first 10): {', '.join(missing_in_sequences[:10])}")
        if missing_in_meta:
            print(f"[info] Missing in meta_data (first 10): {', '.join(missing_in_meta[:10])}")
        if cluster_col and missing_centroids_in_tree:
            print(f"[info] Missing centroids in tree (first 10): {', '.join(missing_centroids_in_tree[:10])}")
        if cluster_col and extra_in_tree:
            print(f"[info] Extra nodes in tree (first 10): {', '.join(extra_in_tree[:10])}")

        fig = plt.figure(figsize=(12, 18))
        ax = fig.add_subplot(1, 1, 1)
        if tree_available:
            Phylo.draw(tree, do_show=False, axes=ax)
        else:
            ax.text(0.5, 0.5, "No tree present in DB", ha="center", va="center", fontsize=16)
            ax.set_axis_off()
        plot_path = os.path.join(args.outdir, "db_tree.png")
        plt.tight_layout()
        fig.savefig(plot_path, dpi=150)

        overlap_count = len(meta_set & tree_terminals)
        disjoint_sets = overlap_count == 0 and len(meta_set) > 0 and len(tree_terminals) > 0

        if args.test_mode and tree_source == "usher" and disjoint_sets:
            print(
                "[warn] Test mode: tree terminals and meta_data are disjoint; "
                "treating as no-placeable-queries scenario and not failing validation."
            )
            return

        if tree_available and args.expect_segment_trees and segment_count is not None and segment_count > 1 and not segmented_validation_ok:
            if usher_tree_count < segment_count:
                raise SystemExit(
                    "Validation failed: segmented run expects at least one UShER tree per segment "
                    f"(segments={segment_count}, usher_trees={usher_tree_count})"
                )

        if args.check_update_integrity:
            if not table_exists(conn, "update_batches") or not table_exists(conn, "update_table_deltas"):
                raise SystemExit("Validation failed: update audit tables are missing")
            cursor.execute("SELECT batch_id FROM update_batches ORDER BY rowid DESC LIMIT 1")
            row = cursor.fetchone()
            if row:
                batch_id = row[0]
                cursor.execute("SELECT COUNT(*) FROM update_table_deltas WHERE batch_id=?", (batch_id,))
                if cursor.fetchone()[0] == 0:
                    raise SystemExit("Validation failed: no update_table_deltas rows found for latest batch")
            if table_exists(conn, "features") and table_exists(conn, "meta_data"):
                feature_cols = get_table_columns(conn, "features")
                if "segment" in feature_cols and "segment" in meta_columns and "accession" in feature_cols:
                    cursor.execute(
                        """
                        SELECT COUNT(*)
                        FROM features f
                        JOIN meta_data m ON m.primary_accession = f.accession
                        WHERE COALESCE(TRIM(f.segment), '') != COALESCE(TRIM(m.segment), '')
                        """
                    )
                    if cursor.fetchone()[0] > 0:
                        raise SystemExit("Validation failed: segment contamination detected in features table")
            if result_failed(feature_integrity_result):
                raise SystemExit("Validation failed: feature projection integrity check failed")

        failed_consistency_checks = [title for title, result in consistency_results.items() if result_failed(result)]
        failed_mutation_checks = [result["title"] for result in mutation_integrity_results if result_failed(result)]
        if failed_consistency_checks or failed_mutation_checks:
            if args.test_mode and args.expect_segment_trees and segment_count is not None and segment_count > 1 and not args.check_update_integrity:
                print("[info] Test mode: skipping strict DB consistency failures for segmented validation run")
            else:
                raise SystemExit(
                    "Validation failed: DB consistency checks failed for " + ", ".join(failed_consistency_checks + failed_mutation_checks)
                )

        if not tree_available:
            return

        if tree_source == "usher":
            if segmented_validation_ok:
                return
            if missing_in_tree or missing_in_sequences or missing_in_meta:
                raise SystemExit("Validation failed: UShER tree does not match accessions")
            if cluster_col and extra_in_tree:
                print(
                    "[info] UShER tree contains non-centroid terminals "
                    f"({len(extra_in_tree)} extra vs centroid context); full accession coverage passed."
                )
        elif cluster_col:
            if missing_centroids_in_tree or extra_in_tree:
                raise SystemExit("Validation failed: IQ-TREE does not match centroid set")
        else:
            if missing_in_tree or missing_in_sequences or missing_in_meta:
                raise SystemExit("Validation failed: missing accessions detected")

    finally:
        conn.close()


if __name__ == "__main__":
    main()
