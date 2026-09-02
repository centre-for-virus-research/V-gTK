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



def describe_missing_accessions(conn, missing, meta_cols, limit=2000):
    """Print the accession_type / exclusion_status breakdown of missing rows.

    A consistency failure lists accessions; what matters is what *kind* of row
    they are. An `exclusion_list` reference with exclusion_status != 1 is a
    different bug from a query that genuinely failed to align, and the accession
    alone does not distinguish them.
    """
    if not missing or "accession_type" not in meta_cols:
        return

    has_status = "exclusion_status" in meta_cols
    sample = [str(a) for a in missing[:limit]]
    marks = ",".join("?" * len(sample))
    columns = "accession_type" + (", exclusion_status" if has_status else "")

    try:
        rows = conn.execute(
            f"SELECT {columns} FROM meta_data WHERE primary_accession IN ({marks})",
            sample,
        ).fetchall()
    except sqlite3.Error:
        return

    tally = {}
    for row in rows:
        kind = str(row[0] or "").strip().lower() or "(blank)"
        status = str(row[1] or "").strip() if has_status else "?"
        tally[(kind, status)] = tally.get((kind, status), 0) + 1

    for (kind, status), count in sorted(tally.items(), key=lambda kv: -kv[1]):
        print(
            f"[info]     -> {count} of them: accession_type={kind!r}"
            + (f" exclusion_status={status!r}" if has_status else "")
        )
        if kind == "exclusion_list" and status != "1":
            print(
                "[warn]        an 'exclusion_list' row with exclusion_status != 1 "
                "is the bug: it is a reference the list says to exclude, so it "
                "will never have an alignment, but it is still being required to "
                "have one. Check FilterAndExtractSequences classified it."
            )


def validate_sequence_alignment_vs_meta(conn, expected_accessions, accession_column="primary_accession", exclusion_column="exclusion_status", exclude_value="1"):
    table_name = "sequence_alignment"
    label = "sequence_alignment vs meta_data"
    if not table_exists(conn, table_name):
        return {"ok": False, "error": f"Table '{table_name}' not found.", "table": table_name, "label": label}

    cols = get_table_columns(conn, table_name)
    col = find_first_present(cols, ["primary_accession", "accession", "sequence_id"])
    if col is None:
        return {
            "ok": False,
            "error": f"Table '{table_name}' missing any of columns ['primary_accession', 'accession', 'sequence_id']. Columns present: {cols}",
            "table": table_name,
            "label": label,
        }

    meta_cols = get_table_columns(conn, "meta_data") if table_exists(conn, "meta_data") else []
    if accession_column not in meta_cols:
        observed = fetch_distinct_values(conn, table_name, col)
    elif exclusion_column in meta_cols and "accession_type" in meta_cols:
        cursor = conn.cursor()
        cursor.execute(
            f"""
            SELECT DISTINCT s.{col} AS v
            FROM {table_name} s
            LEFT JOIN meta_data m
              ON m.{accession_column} = s.{col}
            WHERE s.{col} IS NOT NULL AND s.{col} != ''
              AND (
                    m.{accession_column} IS NULL
                    OR m.{exclusion_column} != ?
                    OR LOWER(COALESCE(m.accession_type, '')) NOT IN ('reference', 'master')
                  )
            """,
            (str(exclude_value),),
        )
        observed = {str(row[0]).strip() for row in cursor.fetchall() if row and row[0]}
    else:
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
    mismatched_cds_rows = 0
    for accession, master_ref_accession, product, aln_start, aln_end, cds_start, cds_end in feature_rows:
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
        cds_start_int = _try_parse_int(cds_start)
        cds_end_int = _try_parse_int(cds_end)
        # Clamp the projected CDS coordinates to the sequence's actual alignment window
        if cds_start_int is not None and cds_end_int is not None:
            cds_start_int = max(cds_start_int, aln_start_int)
            cds_end_int = min(cds_end_int, aln_end_int)
        
        expected_matches = False
        for master_start, master_end in master_product_spans:
            expected_start = max(master_start, aln_start_int)
            expected_end = min(master_end, aln_end_int)
            if expected_start <= expected_end:
                if cds_start_int == expected_start and cds_end_int == expected_end:
                    expected_matches = True
                    break

        if not expected_matches:
            mismatched_cds_rows += 1
            offending_rows.append(
                {
                    "accession": accession_text,
                    "master_ref_accession": master_text,
                    "product": product_text,
                    "aln_start": str(aln_start),
                    "aln_end": str(aln_end),
                    "cds_start": str(cds_start),
                    "cds_end": str(cds_end),
                    "reason": "cds_span_mismatch",
                }
            )
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
                    "cds_start": str(cds_start),
                    "cds_end": str(cds_end),
                    "reason": "aln_span_outside_master_feature",
                }
            )

    return {
        "title": title,
        "ok": len(offending_rows) == 0 and invalid_aln_rows == 0,
        "table": "features",
        "offending_rows": len(offending_rows),
        "invalid_aln_rows": invalid_aln_rows,
        "unresolved_master_rows": unresolved_master_rows,
        "mismatched_cds_rows": mismatched_cds_rows,
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
    lines.append(f"  mismatched_cds_rows: {result.get('mismatched_cds_rows', 0)}")
    examples = result.get("examples", [])
    if examples:
        lines.append(f"  Offending examples (up to {show_n}):")
        for example in examples[:show_n]:
            lines.append(
                "    - accession={accession} product={product} master_ref_accession={master_ref_accession} aln_start={aln_start} aln_end={aln_end} cds_start={cds_start} cds_end={cds_end} reason={reason}".format(
                    **example
                )
            )
    return "\n".join(lines)
# ---------------------------------------------------------------------------
# Extended DB invariants
#
# The checks below answer "is this finished database internally coherent?" as
# opposed to "did the pipeline crash?".  Every failure mode they cover is one
# that produces a database that opens fine, loads fine, and is quietly wrong:
# a column that stopped being written, a gene whose name was eaten by a bad
# delimiter, an alignment that grew residues the submitted sequence never had,
# a tree leaf that resolves to nothing.  They are reported as warnings by
# default (see --strict-invariants) because a legitimate build can trip some of
# them (e.g. a reference list with no genotype labels).
# ---------------------------------------------------------------------------

NUCLEOTIDE_ALPHABET = "ACGTURYSWKMBDHVN"
ALIGNMENT_EXTRA_CHARS = "-."
GENE_PARENT_SENTINELS = {"", "na", "n/a", "null", "none", "nan", "-"}

# meta_data columns that a finished DB should essentially always have populated
# for at least one row.  A column here that is 100% blank means the writer for
# it silently stopped producing values (the influenza neuraminidase "NA" case
# erased a whole column exactly this way).
CORE_META_COLUMNS = [
    "primary_accession",
    "accession_version",
    "organism",
    "taxonomy",
    "accession_type",
    "length",
    "real_length",
    "segment",
    "collection_date",
    "collection_year",
    "country",
    "host",
    "host_taxa_id",
    "host_scientific_name",
    "country_validated",
]

# Columns whose emptiness means "the labelling mechanism never ran".  Kept
# separate so the report can say precisely what is missing.
REFERENCE_LABEL_META_COLUMNS = ["nearest_reference_genotype", "nearest_reference_subtype"]

MONTH_ABBREVIATIONS = {
    "jan": 1, "feb": 2, "mar": 3, "apr": 4, "may": 5, "jun": 6,
    "jul": 7, "aug": 8, "sep": 9, "oct": 10, "nov": 11, "dec": 12,
}


def quote_ident(name):
    """Quote a SQLite identifier. Column names come from PRAGMA table_info, but
    quoting keeps a column literally called "segment"/"Segment"/"index" safe."""
    return '"' + str(name).replace('"', '""') + '"'


def invariant_result(title, findings=None, skipped=False, reason=None, error=None, severity="warning"):
    result = {"title": title, "severity": severity, "findings": findings or []}
    if error is not None:
        result["ok"] = False
        result["error"] = error
        return result
    if skipped:
        result["ok"] = True
        result["skipped"] = True
        result["reason"] = reason or "not applicable to this DB"
        return result
    result["ok"] = not result["findings"]
    return result


def add_finding(findings, code, count, detail, examples=None):
    if count:
        findings.append(
            {
                "code": code,
                "count": int(count),
                "detail": detail,
                "examples": [str(x) for x in (examples or [])][:25],
            }
        )


def _blank_sql(column):
    col = quote_ident(column)
    return f"({col} IS NULL OR TRIM(CAST({col} AS TEXT)) = '')"


def validate_meta_column_population(conn, meta_cols, where_sql=None, params=()):
    """Flag meta_data columns that are 100% NULL/blank.

    Real trigger: pandas read the influenza gene name "NA" as a null and the
    shipped DB ended up with a nameless gene; the same class of accident makes a
    whole metadata column silently empty, and nothing else in the pipeline
    notices because empty is a legal value everywhere."""
    title = "meta_data column population"
    if not table_exists(conn, "meta_data"):
        return invariant_result(title, skipped=True, reason="Table 'meta_data' not present.")

    cursor = conn.cursor()
    where_clause = f" WHERE {where_sql}" if where_sql else ""
    empty_core = []
    empty_other = []
    empty_labels = []

    # One pass over the table, not one per column. A query per column is a full
    # table scan per column, and meta_data is wide: on a real influenza build
    # that is 106 scans of a 28 GB table - about 3 TB of reads, and it ran for
    # half an hour without finishing. Rolled into a single SELECT of N
    # conditional sums it is one scan. Chunked because SQLite caps an expression
    # tree at SQLITE_MAX_COLUMN (2000 by default), and to keep the SQL readable.
    populated = {}
    chunk_size = 64
    for start in range(0, len(meta_cols), chunk_size):
        chunk = meta_cols[start:start + chunk_size]
        sums = ", ".join(
            f"SUM(CASE WHEN NOT {_blank_sql(col)} THEN 1 ELSE 0 END)" for col in chunk
        )
        cursor.execute(f"SELECT {sums} FROM meta_data{where_clause}", params)
        row = cursor.fetchone() or ()
        for col, value in zip(chunk, row):
            populated[col] = int(value or 0)

    for col in meta_cols:
        if populated.get(col):
            continue
        if col in REFERENCE_LABEL_META_COLUMNS:
            empty_labels.append(col)
        elif col in CORE_META_COLUMNS:
            empty_core.append(col)
        else:
            empty_other.append(col)

    findings = []
    add_finding(
        findings,
        "core_meta_column_entirely_empty",
        len(empty_core),
        "meta_data columns that a finished DB should populate are 100% NULL/blank",
        empty_core,
    )
    add_finding(
        findings,
        "reference_label_columns_unpopulated",
        len(empty_labels),
        (
            "genotype/subtype label columns exist but not one row carries a value; "
            "if the reference list ships genotype labels the labelling step did not run"
        ),
        empty_labels,
    )
    result = invariant_result(title, findings)
    result["empty_other_columns"] = empty_other
    result["checked_columns"] = len(meta_cols)
    return result


def validate_segment_labels(conn, meta_cols, exclusion_column="exclusion_status", exclude_value="1"):
    """Segment must be a segment label, present, and agree across every table.

    Real trigger: a segment/Segment case collision plus ten divergent segment
    normalisers left most rows with an empty segment and some carrying a gene
    name instead of a segment number.  An empty or gene-named segment does not
    raise anything - it just silently partitions the DB wrongly."""
    title = "segment label integrity"
    if not table_exists(conn, "meta_data") or "segment" not in meta_cols:
        return invariant_result(title, skipped=True, reason="meta_data has no 'segment' column.")

    cursor = conn.cursor()
    findings = []

    where_sql, params = get_meta_nonexcluded_filter(meta_cols, exclusion_column, exclude_value)
    suffix = f" AND ({where_sql})" if where_sql else ""
    cursor.execute(f"SELECT COUNT(*) FROM meta_data WHERE {_blank_sql('segment')}{suffix}", params)
    blank_segments = int(cursor.fetchone()[0])
    add_finding(
        findings,
        "blank_segment_in_meta_data",
        blank_segments,
        "non-excluded meta_data rows carry no segment label",
    )

    meta_segments = fetch_distinct_values(conn, "meta_data", "segment")

    if len(meta_segments) == 1:
        only = next(iter(meta_segments))
        if only != "1":
            add_finding(
                findings,
                "single_segment_not_labelled_1",
                1,
                "a non-segmented DB should label every row segment='1'",
                [only],
            )

    for table_name, _acc_col in (
        ("sequences", "header"),
        ("sequence_alignment", "primary_accession"),
        ("features", "accession"),
        ("insertions", "primary_accession"),
        ("trees", "name"),
    ):
        if not table_exists(conn, table_name):
            continue
        if "segment" not in get_table_columns(conn, table_name):
            continue
        observed = fetch_distinct_values(conn, table_name, "segment")
        unknown = sorted(observed - meta_segments)
        add_finding(
            findings,
            f"segment_absent_from_meta_data:{table_name}",
            len(unknown),
            f"{table_name}.segment values that no meta_data row uses",
            unknown,
        )

    gene_like = set()
    if table_exists(conn, "genes"):
        gene_cols = get_table_columns(conn, "genes")
        for col in ("name", "display_name"):
            if col in gene_cols:
                gene_like |= {v.lower() for v in fetch_distinct_values(conn, "genes", col)}
    if table_exists(conn, "features") and "product" in get_table_columns(conn, "features"):
        gene_like |= {v.lower() for v in fetch_distinct_values(conn, "features", "product")}
    gene_like.discard("whole_genome")
    segment_gene_names = sorted(v for v in meta_segments if v.lower() in gene_like)
    add_finding(
        findings,
        "segment_holds_a_gene_name",
        len(segment_gene_names),
        "meta_data.segment values that are gene/product names rather than segment labels",
        segment_gene_names,
    )

    result = invariant_result(title, findings)
    result["meta_segments"] = sorted(meta_segments)
    return result


def validate_referential_integrity(conn, expected_accessions, accession_column="primary_accession"):
    """Every child-table row must point at a meta_data accession, and every
    non-excluded accession must own a raw sequence.

    Real trigger: an update run that re-wrote meta_data but not its children (or
    vice versa) leaves rows addressed to accessions the DB no longer describes.
    Those rows are invisible to every join the toolkit does, so nothing errors -
    the data is simply unreachable."""
    title = "cross-table referential integrity"
    if not table_exists(conn, "meta_data"):
        return invariant_result(title, skipped=True, reason="Table 'meta_data' not present.")

    cursor = conn.cursor()
    meta_accession_set = fetch_distinct_values(conn, "meta_data", accession_column)
    findings = []
    children = [
        ("sequences", "header", True),
        ("sequence_alignment", "primary_accession", False),
        ("insertions", "primary_accession", False),
        ("features", "accession", False),
    ]
    for table_name, child_col, check_reverse in children:
        if not table_exists(conn, table_name):
            continue
        cols = get_table_columns(conn, table_name)
        col = find_first_present(cols, [child_col, "primary_accession", "accession", "header"])
        if col is None:
            continue
        # Set difference in Python rather than a LEFT JOIN on TRIM(CAST(...)):
        # the expression wrapper stops SQLite building an automatic index, and
        # the join degrades to O(rows x rows) - minutes on a production DB for
        # an answer two hash sets give in a second.
        observed_ids = fetch_distinct_values(conn, table_name, col)
        orphans = sorted(observed_ids - meta_accession_set)
        add_finding(
            findings,
            f"orphan_rows:{table_name}",
            len(orphans),
            f"{table_name}.{col} values with no matching meta_data.{accession_column}",
            orphans,
        )
        if check_reverse:
            observed = fetch_distinct_values(conn, table_name, col)
            missing = sorted(set(expected_accessions) - observed)
            add_finding(
                findings,
                f"missing_from_child:{table_name}",
                len(missing),
                f"non-excluded meta_data accessions with no row in {table_name}",
                missing,
            )

    # Pointer columns: each names the reference an accession was aligned or
    # projected against.  sequence_alignment.alignment_name in particular is how
    # CreateSqliteDB._add_reference_columns resolves a query's genotype from its
    # best BLAST hit - an unresolvable value there does not raise, it just makes
    # nearest_reference_genotype come back empty for that query.
    pointer_columns = [
        ("sequence_alignment", "alignment_name"),
        ("insertions", "reference"),
        ("features", "master_ref_accession"),
        ("features", "reference_accession"),
    ]
    for table_name, pointer_col in pointer_columns:
        if not table_exists(conn, table_name):
            continue
        if pointer_col not in get_table_columns(conn, table_name):
            continue
        observed = fetch_distinct_values(conn, table_name, pointer_col)
        unresolved = sorted(observed - meta_accession_set)
        add_finding(
            findings,
            f"reference_pointer_unresolved:{table_name}.{pointer_col}",
            len(unresolved),
            f"{table_name}.{pointer_col} values that name no meta_data.{accession_column}",
            unresolved,
        )
    return invariant_result(title, findings)


def _alphabet_offenders(text, allowed):
    """Return True when `text` uses a character outside `allowed`.

    Done in Python on a streamed row rather than in SQL: a GLOB negated class
    (`*[^ACGT...]*`) backtracks from every offset and costs ~150s over a
    production alignment table, and a nested-REPLACE chain copies the column
    once per allowed letter.  One set() per row is an order of magnitude
    cheaper than either."""
    return bool(set(text) - allowed)


def validate_sequence_content(conn, accession_column="primary_accession"):
    """Stored sequences must be non-empty, nucleotide-alphabet, and must agree
    with the lengths meta_data advertises for them.

    Real trigger: meta_data.length / real_length / a,t,g,c,n are computed by a
    different step than the one that stores the sequence.  If either side is
    written for the wrong accession, or a sequence is truncated on the way in,
    both tables still look plausible on their own; only the cross-check shows
    it.  meta_data.length must equal the stored sequence length and real_length
    must equal a+t+g+c - both hold exactly on every shipped DB."""
    title = "sequence content and length coherence"
    if not table_exists(conn, "sequences"):
        return invariant_result(title, skipped=True, reason="Table 'sequences' not present.")

    seq_cols = get_table_columns(conn, "sequences")
    header_col = find_first_present(seq_cols, ["header", "primary_accession", "accession"])
    if header_col is None or "sequence" not in seq_cols:
        return invariant_result(title, skipped=True, reason=f"sequences table lacks header/sequence columns: {seq_cols}")

    findings = []
    qh = quote_ident(header_col)
    allowed = set(NUCLEOTIDE_ALPHABET) | set(NUCLEOTIDE_ALPHABET.lower())

    # One streaming pass over the sequence text: emptiness, alphabet and stored
    # length all come out of it, so the column is read once rather than once per
    # check (it is over a gigabyte on a production HCV database).
    empty_rows = 0
    bad_alphabet = []
    bad_alphabet_count = 0
    stored_lengths = {}
    for header, sequence in conn.execute(f"SELECT {qh}, sequence FROM sequences"):
        if sequence is None or not str(sequence).strip():
            empty_rows += 1
            continue
        text = str(sequence)
        key = str(header).strip() if header is not None else ""
        if key:
            stored_lengths[key] = len(text)
        if _alphabet_offenders(text, allowed):
            bad_alphabet_count += 1
            if len(bad_alphabet) < 25:
                bad_alphabet.append(key)

    add_finding(findings, "empty_sequence_rows", empty_rows, "sequences rows with a NULL/zero-length sequence")
    add_finding(
        findings,
        "sequence_outside_nucleotide_alphabet",
        bad_alphabet_count,
        f"sequences containing characters outside {NUCLEOTIDE_ALPHABET}",
        bad_alphabet,
    )

    meta_cols = get_table_columns(conn, "meta_data") if table_exists(conn, "meta_data") else []
    qacc = quote_ident(accession_column)
    if accession_column in meta_cols and "length" in meta_cols:
        rows = []
        for accession, declared in conn.execute(f'SELECT {qacc}, "length" FROM meta_data'):
            if accession is None:
                continue
            key = str(accession).strip()
            declared_text = "" if declared is None else str(declared).strip()
            if not declared_text or key not in stored_lengths:
                continue
            try:
                declared_int = int(float(declared_text))
            except ValueError:
                continue
            if declared_int != stored_lengths[key]:
                rows.append((key, declared_text, stored_lengths[key]))
        add_finding(
            findings,
            "meta_length_disagrees_with_stored_sequence",
            len(rows),
            "meta_data.length differs from LENGTH(sequences.sequence)",
            [f"{r[0]} meta={r[1]} stored={r[2]}" for r in rows],
        )

    if accession_column in meta_cols and {"real_length", "a", "t", "g", "c"}.issubset(set(meta_cols)):
        cursor = conn.cursor()
        cursor.execute(
            f"""
            SELECT m.{qacc}, m."real_length",
                   (CAST(m."a" AS INTEGER)+CAST(m."t" AS INTEGER)+CAST(m."g" AS INTEGER)+CAST(m."c" AS INTEGER))
            FROM meta_data m
            WHERE NOT {_blank_sql('real_length')}
              AND CAST(m."real_length" AS INTEGER) !=
                  (CAST(m."a" AS INTEGER)+CAST(m."t" AS INTEGER)+CAST(m."g" AS INTEGER)+CAST(m."c" AS INTEGER))
            """
        )
        rows = cursor.fetchall()
        add_finding(
            findings,
            "real_length_disagrees_with_base_counts",
            len(rows),
            "meta_data.real_length differs from a+t+g+c",
            [f"{r[0]} real_length={r[1]} atgc={r[2]}" for r in rows],
        )
    return invariant_result(title, findings)


def validate_alignment_geometry(conn, accession_column="primary_accession"):
    """An alignment row is the stored sequence with gaps inserted, so stripping
    the gaps can never make it longer than the raw sequence, and every row in
    one alignment must have the same width.

    Real trigger: the padding step splices reference flanks onto short queries.
    When it over-reaches it appends residues the submitted sequence never had -
    the DB then serves an aligned sequence that disagrees with its own raw
    sequence, with nothing logged.  A ragged alignment width means one row was
    padded against a different reference than the rest of its segment."""
    title = "alignment geometry"
    if not table_exists(conn, "sequence_alignment"):
        return invariant_result(title, skipped=True, reason="Table 'sequence_alignment' not present.")

    aln_cols = get_table_columns(conn, "sequence_alignment")
    acc_col = find_first_present(aln_cols, ["primary_accession", "accession", "sequence_id"])
    if acc_col is None or "alignment" not in aln_cols:
        return invariant_result(title, skipped=True, reason=f"sequence_alignment lacks accession/alignment columns: {aln_cols}")

    findings = []
    qa = quote_ident(acc_col)
    allowed = set(NUCLEOTIDE_ALPHABET + ALIGNMENT_EXTRA_CHARS)
    allowed |= {c.lower() for c in allowed}

    raw_lengths = {}
    if table_exists(conn, "sequences"):
        seq_cols = get_table_columns(conn, "sequences")
        header_col = find_first_present(seq_cols, ["header", "primary_accession", "accession"])
        if header_col and "sequence" in seq_cols:
            qh = quote_ident(header_col)
            # LENGTH() is evaluated inside SQLite so no sequence text is
            # transferred; the accession match then happens in a dict.  A SQL
            # join on TRIM(CAST(...)) cannot use an index and runs quadratically.
            for header, seq_len in conn.execute(f"SELECT {qh}, LENGTH(sequence) FROM sequences WHERE sequence IS NOT NULL"):
                if header is not None:
                    raw_lengths[str(header).strip()] = seq_len

    group_col = "segment" if "segment" in aln_cols else ("alignment_name" if "alignment_name" in aln_cols else None)
    select_cols = [qa, "alignment"] + ([quote_ident(group_col)] if group_col else [])

    empty_rows = 0
    all_gap = []
    all_gap_count = 0
    bad_alphabet = []
    bad_alphabet_count = 0
    over_long = []
    widths = {}
    # Single streaming pass: emptiness, all-gap rows, alphabet, ungapped length
    # and per-group width all come from the same read of the alignment column.
    for row in conn.execute(f"SELECT {', '.join(select_cols)} FROM sequence_alignment"):
        accession = str(row[0]).strip() if row[0] is not None else ""
        alignment = row[1]
        group_value = str(row[2]).strip() if group_col and row[2] is not None else ""
        if alignment is None or not str(alignment).strip():
            empty_rows += 1
            continue
        text = str(alignment)
        if group_col:
            widths.setdefault(group_value, set()).add(len(text))
        ungapped = len(text) - text.count("-") - text.count(".")
        if ungapped == 0:
            all_gap_count += 1
            if len(all_gap) < 25:
                all_gap.append(accession)
        if _alphabet_offenders(text, allowed):
            bad_alphabet_count += 1
            if len(bad_alphabet) < 25:
                bad_alphabet.append(accession)
        raw_len = raw_lengths.get(accession)
        if raw_len is not None and ungapped > raw_len:
            over_long.append((accession, ungapped, raw_len))

    add_finding(findings, "empty_alignment_rows", empty_rows, "sequence_alignment rows with no alignment string")
    add_finding(findings, "all_gap_alignment_rows", all_gap_count, "alignment rows that are nothing but gaps", all_gap)
    add_finding(
        findings,
        "alignment_outside_nucleotide_alphabet",
        bad_alphabet_count,
        f"alignments containing characters outside {NUCLEOTIDE_ALPHABET}{ALIGNMENT_EXTRA_CHARS}",
        bad_alphabet,
    )
    add_finding(
        findings,
        "alignment_longer_than_raw_sequence",
        len(over_long),
        "ungapped alignment has more residues than the stored raw sequence (padding invented bases)",
        [f"{r[0]} ungapped={r[1]} raw={r[2]}" for r in over_long],
    )
    ragged = sorted((g, lens) for g, lens in widths.items() if len(lens) > 1)
    add_finding(
        findings,
        "ragged_alignment_width",
        len(ragged),
        f"alignment rows sharing a {group_col} disagree on alignment length" if group_col else "alignment rows disagree on length",
        [f"{group_col}={g} distinct_lengths={len(lens)} min={min(lens)} max={max(lens)}" for g, lens in ragged],
    )
    return invariant_result(title, findings)


def validate_duplicate_records(conn, accession_column="primary_accession"):
    """One accession may appear once per segment, never twice.

    Real trigger: an incremental --update that re-inserts instead of upserting
    doubles rows.  Counts still look sane, joins still work, but every
    downstream aggregate silently double-counts the duplicated accessions."""
    title = "duplicate record keys"
    findings = []
    cursor = conn.cursor()
    targets = [
        ("meta_data", accession_column),
        ("sequences", "header"),
        ("sequence_alignment", "primary_accession"),
    ]
    for table_name, key_col in targets:
        if not table_exists(conn, table_name):
            continue
        cols = get_table_columns(conn, table_name)
        col = find_first_present(cols, [key_col, "primary_accession", "accession", "header"])
        if col is None:
            continue
        group_cols = [quote_ident(col)]
        if "segment" in cols:
            group_cols.append("COALESCE(TRIM(CAST(\"segment\" AS TEXT)), '')")
        group_sql = ", ".join(group_cols)
        # NB: HAVING must reference COUNT(*) directly. meta_data has a column
        # literally named "n" (the ambiguous-base count), so `HAVING n > 1`
        # silently resolves to that column instead of the alias and reports
        # every row with more than one N as a duplicate.
        cursor.execute(
            f"""
            SELECT {quote_ident(col)}, COUNT(*)
            FROM {quote_ident(table_name)}
            WHERE {quote_ident(col)} IS NOT NULL AND TRIM(CAST({quote_ident(col)} AS TEXT)) != ''
            GROUP BY {group_sql}
            HAVING COUNT(*) > 1
            """
        )
        rows = cursor.fetchall()
        add_finding(
            findings,
            f"duplicate_key:{table_name}",
            len(rows),
            f"{table_name} has repeated ({col}, segment) keys",
            [f"{r[0]} x{r[1]}" for r in rows],
        )
    return invariant_result(title, findings)


def validate_gene_table(conn):
    """Every gene needs a name, and its parent must exist.

    Real trigger: generic/rabv/Tables/gene_info.csv and
    generic/other/Tables/gene_info.csv end their last row with spaces / a comma
    instead of a tab, so the whole_genome row is loaded as a single field named
    'whole_genome    NULL'.  Every other gene points at parent 'whole_genome',
    which then does not exist - the gene hierarchy is silently broken.  The
    influenza gene 'NA' being read as a null produced a nameless gene the same
    way."""
    title = "genes table integrity"
    if not table_exists(conn, "genes"):
        return invariant_result(title, skipped=True, reason="Table 'genes' not present.")

    gene_cols = get_table_columns(conn, "genes")
    if "name" not in gene_cols:
        return invariant_result(title, skipped=True, reason=f"genes table has no 'name' column: {gene_cols}")

    cursor = conn.cursor()
    findings = []

    cursor.execute(f"SELECT COUNT(*) FROM genes WHERE {_blank_sql('name')}")
    blank_names = cursor.fetchone()[0]
    cursor.execute(f"SELECT description FROM genes WHERE {_blank_sql('name')} LIMIT 25")
    add_finding(
        findings,
        "gene_with_no_name",
        blank_names,
        "genes rows with a blank name (a gene name that looks like a null, e.g. influenza 'NA', was erased)",
        [str(r[0]) for r in cursor.fetchall()],
    )

    if "display_name" in gene_cols:
        cursor.execute(f"SELECT COUNT(*) FROM genes WHERE {_blank_sql('display_name')}")
        add_finding(findings, "gene_with_no_display_name", cursor.fetchone()[0], "genes rows with a blank display_name")

    cursor.execute("SELECT name FROM genes WHERE name IS NOT NULL")
    names = [str(r[0]) for r in cursor.fetchall()]
    merged = [n for n in names if "\t" in n or "  " in n.strip() or "," in n]
    add_finding(
        findings,
        "gene_name_looks_like_a_merged_field",
        len(merged),
        "gene names containing a tab, a comma or a run of spaces - the source row was split on the wrong delimiter",
        merged,
    )

    name_set = {n.strip() for n in names if n.strip()}
    if "parent_name" in gene_cols:
        cursor.execute("SELECT name, parent_name FROM genes")
        dangling = []
        for name, parent in cursor.fetchall():
            parent_text = str(parent).strip() if parent is not None else ""
            if parent_text.lower() in GENE_PARENT_SENTINELS:
                continue
            if parent_text not in name_set:
                dangling.append(f"{name} -> {parent_text}")
        add_finding(
            findings,
            "gene_parent_does_not_resolve",
            len(dangling),
            "genes.parent_name values that match no genes.name",
            dangling,
        )

    cursor.execute(
        f"SELECT name, COUNT(*) FROM genes WHERE NOT {_blank_sql('name')} GROUP BY name HAVING COUNT(*) > 1"
    )
    rows = cursor.fetchall()
    add_finding(findings, "duplicate_gene_name", len(rows), "genes table has repeated names", [f"{r[0]} x{r[1]}" for r in rows])
    return invariant_result(title, findings)


def validate_tree_tip_resolution(conn, accession_column="primary_accession"):
    """Every leaf label in every stored tree must name a meta_data accession.

    Real trigger: the reference alignment carries backbone genomes that are not
    part of the curated reference list, so they end up as leaves in the segment
    trees but have no metadata, no sequence and no alignment anywhere in the DB.
    Clicking one in the toolkit resolves to nothing.  The existing validator
    only cross-checks the single 'best' tree, and in segmented mode it returns
    before that check runs at all."""
    title = "tree tip resolution"
    if not table_exists(conn, "trees"):
        return invariant_result(title, skipped=True, reason="Table 'trees' not present.")

    cursor = conn.cursor()
    cursor.execute("SELECT COUNT(*) FROM trees WHERE newick IS NOT NULL AND LENGTH(TRIM(newick)) > 0")
    if cursor.fetchone()[0] == 0:
        return invariant_result(title, skipped=True, reason="No non-empty trees stored (tree-free build).")

    meta_accessions = set()
    if table_exists(conn, "meta_data"):
        meta_accessions = fetch_distinct_values(conn, "meta_data", accession_column)

    findings = []
    cursor.execute("SELECT COUNT(*) FROM trees WHERE newick IS NULL OR LENGTH(TRIM(newick)) = 0")
    add_finding(findings, "empty_tree_rows", cursor.fetchone()[0], "trees rows with a NULL/empty newick")

    unresolved = set()
    unparsable = []
    per_tree = {}
    trees_seen = 0
    for row in fetch_trees(conn):
        trees_seen += 1
        try:
            tree = Phylo.read(StringIO(row["newick"]), "newick")
        except Exception as exc:  # a malformed newick is itself a finding
            unparsable.append(f"{row['name']}: {exc}")
            continue
        tips = {tip.name.strip() for tip in tree.get_terminals() if tip.name}
        if not tips:
            unparsable.append(f"{row['name']}: tree has no labelled tips")
            continue
        missing = tips - meta_accessions
        if missing:
            per_tree[row["name"]] = len(missing)
            unresolved |= missing

    add_finding(findings, "unparsable_tree", len(unparsable), "stored newick strings that could not be parsed", unparsable)
    add_finding(
        findings,
        "tree_tip_absent_from_meta_data",
        len(unresolved),
        "distinct tree leaf labels that match no meta_data accession",
        sorted(unresolved),
    )

    result = invariant_result(title, findings)
    result["trees_checked"] = trees_seen
    result["trees_with_unresolved_tips"] = per_tree
    return result


def validate_segment_tree_coverage(conn, meta_cols):
    """Each segment carried in meta_data should have at least one tree.

    Real trigger: a segmented build where one segment's tree job failed still
    writes a complete-looking DB; the segment simply has no phylogeny and the
    toolkit shows an empty tree panel for it."""
    title = "segment tree coverage"
    if not table_exists(conn, "trees") or "segment" not in meta_cols:
        return invariant_result(title, skipped=True, reason="No trees table or no meta_data.segment column.")
    tree_cols = get_table_columns(conn, "trees")
    if "segment" not in tree_cols:
        return invariant_result(title, skipped=True, reason="trees table has no 'segment' column.")

    tree_segments = fetch_distinct_values(conn, "trees", "segment", where_sql="newick IS NOT NULL AND LENGTH(TRIM(newick)) > 0")
    if not tree_segments:
        return invariant_result(title, skipped=True, reason="No segment-labelled trees stored.")
    meta_segments = fetch_distinct_values(conn, "meta_data", "segment")
    missing = sorted(meta_segments - tree_segments)
    findings = []
    add_finding(findings, "segment_without_tree", len(missing), "meta_data segments with no tree row", missing)
    return invariant_result(title, findings)


def _parse_collection_date(text):
    """Return (year, month, day) with None for parts the string does not carry.
    Handles the GenBank 'DD-Mon-YYYY'/'Mon-YYYY'/'YYYY' forms and the ISO
    'YYYY-MM-DD'/'YYYY-MM' forms that GISAID and newer records use."""
    raw = "" if text is None else str(text).strip()
    if not raw or raw.lower() in {"na", "n/a", "-", "unknown", "none", "null"}:
        return None, None, None
    parts = raw.split("-")
    try:
        if len(parts) == 3 and len(parts[0]) == 4 and parts[0].isdigit():
            return int(parts[0]), int(parts[1]), int(parts[2])
        if len(parts) == 3 and parts[0].isdigit():
            return int(parts[2]), MONTH_ABBREVIATIONS.get(parts[1][:3].lower()), int(parts[0])
        if len(parts) == 2 and len(parts[0]) == 4 and parts[0].isdigit():
            return int(parts[0]), int(parts[1]), None
        if len(parts) == 2:
            return int(parts[1]), MONTH_ABBREVIATIONS.get(parts[0][:3].lower()), None
        if len(parts) == 1 and parts[0].isdigit() and len(parts[0]) == 4:
            return int(parts[0]), None, None
    except (ValueError, IndexError):
        return None, None, None
    return None, None, None


def validate_collection_dates(conn, meta_cols, today=None):
    """Collection dates must be in the past, plausible, and must not lose their
    day/month on the way into the split columns.

    Real trigger: the day/month splitter only understands GenBank's
    'DD-Mon-YYYY'.  Every record whose collection_date arrives as ISO
    'YYYY-MM-DD' - which is what GISAID and recent GenBank submissions use -
    keeps its year and silently loses day and month, so date-resolution filters
    quietly drop those samples.  A future collection_date means a two-digit year
    was expanded into the wrong century."""
    title = "collection date sanity"
    if not table_exists(conn, "meta_data") or "collection_date" not in meta_cols:
        return invariant_result(title, skipped=True, reason="meta_data has no 'collection_date' column.")

    import datetime as _datetime

    today = today or _datetime.date.today()
    cursor = conn.cursor()
    select_cols = ["primary_accession" if "primary_accession" in meta_cols else "rowid", "collection_date"]
    for optional in ("collection_year", "collection_mon", "collection_day"):
        if optional in meta_cols:
            select_cols.append(optional)
    cursor.execute(f"SELECT {', '.join(quote_ident(c) for c in select_cols)} FROM meta_data")

    future_dates = []
    implausible = []
    year_disagreements = []
    lost_precision = []
    # Stream the cursor rather than fetchall(): a production DB carries millions
    # of meta_data rows and this check does not need them all in memory at once.
    for row in cursor:
        record = dict(zip(select_cols, row))
        accession = str(record[select_cols[0]])
        year, month, day = _parse_collection_date(record.get("collection_date"))
        if year is None:
            continue
        if year > today.year or (
            month and day and _datetime.date(year, month, day) > today
        ):
            future_dates.append(f"{accession} collection_date={record.get('collection_date')}")
        if year < 1900:
            implausible.append(f"{accession} collection_date={record.get('collection_date')}")
        stored_year = str(record.get("collection_year") or "").strip()
        if stored_year and stored_year.isdigit() and int(stored_year) != year:
            year_disagreements.append(
                f"{accession} collection_date={record.get('collection_date')} collection_year={stored_year}"
            )
        blank_parts = []
        if month is not None and "collection_mon" in select_cols and not str(record.get("collection_mon") or "").strip():
            blank_parts.append("collection_mon")
        if day is not None and "collection_day" in select_cols and not str(record.get("collection_day") or "").strip():
            blank_parts.append("collection_day")
        if blank_parts:
            lost_precision.append(
                f"{accession} collection_date={record.get('collection_date')} blank={'+'.join(blank_parts)}"
            )

    findings = []
    add_finding(findings, "collection_date_in_the_future", len(future_dates), "collection dates later than today", future_dates)
    add_finding(findings, "collection_year_before_1900", len(implausible), "collection years that predate viral sequencing", implausible)
    add_finding(
        findings,
        "collection_year_disagrees_with_collection_date",
        len(year_disagreements),
        "collection_year does not match the year in collection_date",
        year_disagreements,
    )
    add_finding(
        findings,
        "date_precision_lost_in_split_columns",
        len(lost_precision),
        "collection_date carries a month/day that collection_mon/collection_day did not receive",
        lost_precision,
    )
    return invariant_result(title, findings)


def run_extended_invariant_checks(
    conn,
    meta_cols,
    expected_accessions,
    accession_column="primary_accession",
    exclusion_column="exclusion_status",
    exclude_value="1",
    where_sql=None,
    where_params=(),
):
    """Run every extended invariant check and return the ordered result list."""
    return [
        validate_meta_column_population(conn, meta_cols, where_sql=where_sql, params=where_params),
        validate_segment_labels(conn, meta_cols, exclusion_column=exclusion_column, exclude_value=exclude_value),
        validate_referential_integrity(conn, expected_accessions, accession_column=accession_column),
        validate_duplicate_records(conn, accession_column=accession_column),
        validate_sequence_content(conn, accession_column=accession_column),
        validate_alignment_geometry(conn, accession_column=accession_column),
        validate_gene_table(conn),
        validate_tree_tip_resolution(conn, accession_column=accession_column),
        validate_segment_tree_coverage(conn, meta_cols),
        validate_collection_dates(conn, meta_cols),
    ]


def format_invariant_result(result, show_n=10):
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
    other_empty = result.get("empty_other_columns")
    if other_empty:
        # Informational only: optional columns (GISAID merge artefacts, run-mode
        # specific fields) that are legitimately empty in many builds.  Listed so
        # a reviewer can spot one that should NOT have been empty.
        lines.append(f"  (also entirely empty, not treated as a finding: {', '.join(other_empty)})")
    if not result.get("findings"):
        return "\n".join(lines)
    for finding in result["findings"]:
        lines.append(f"  - {finding['code']}: {finding['count']} ({finding['detail']})")
        for example in finding.get("examples", [])[:show_n]:
            lines.append(f"      * {example}")
    return "\n".join(lines)


def get_extra_tree_label(tree_source):
    if tree_source == "usher":
        return "Non-centroid UShER leaves"
    return "Extra nodes in tree"


def validation_fail(message):
    print("[info] Validation status: FAIL")
    raise SystemExit(message)


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
        "--segment-tree-source",
        choices=["usher", "iqtree"],
        default="usher",
        help="Tree source to use for segmented validation",
    )
    parser.add_argument(
        "--check-update-integrity",
        action="store_true",
        help="Validate update audit tables and update-mode integrity constraints",
    )
    parser.add_argument(
        "--allow-empty-update",
        action="store_true",
        help="Allow validation to pass even if the latest update batch made no database changes",
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
    parser.add_argument(
        "--skip-invariant-checks",
        action="store_true",
        help="Skip the extended DB invariant checks (column population, segment labels, referential integrity, sequence/alignment geometry, genes, tree tips, dates)",
    )
    parser.add_argument(
        "--strict-invariants",
        action="store_true",
        help="Fail validation when an extended DB invariant check reports a finding (default: report as warnings)",
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
            validation_fail("No tree found in DB (trees table is empty)")

        tree_available = bool(newick)
        tree = None
        tree_terminals = set()
        if tree_available:
            tree = Phylo.read(StringIO(newick), "newick")
            all_trees = fetch_trees(conn, source=tree_source)
            for tr in all_trees:
                t = Phylo.read(StringIO(tr["newick"]), "newick")
                tree_terminals.update(tip.name for tip in t.get_terminals() if tip.name)

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

        # Collected before anything else that needs it: the reference/master set
        # is required whether or not the database has a cluster column.
        reference_master_set = set()
        reference_master_set_by_segment = {}
        if "accession_type" in meta_columns:
            # `segment` is optional: non-segmented databases have no such column,
            # and selecting it unconditionally raises OperationalError. This ran
            # inside the cluster-column branch before, which happened to hide the
            # problem because those databases were never reached.
            has_segment = "segment" in meta_columns
            cursor = conn.cursor()
            select_cols = f"{args.accession_column}, accession_type"
            if has_segment:
                select_cols += ", segment"
            cursor.execute(f"SELECT {select_cols} FROM meta_data")
            for row in cursor.fetchall():
                acc, accession_type = row[0], row[1]
                seg = row[2] if has_segment else None
                if not acc:
                    continue
                accession_text = str(acc).strip()
                if not accession_text:
                    continue
                accession_type_text = str(accession_type).strip().lower() if accession_type is not None else ""
                if accession_type_text not in {"reference", "master"}:
                    continue
                reference_master_set.add(accession_text)
                seg_str = str(seg).strip() if seg is not None and str(seg).strip() else None
                if seg_str:
                    reference_master_set_by_segment.setdefault(seg_str, set()).add(accession_text)

        # Every reference and the master must appear in a tree. They are the
        # backbone - cluster seeds, tree topology, and the anchors genotype calls
        # are made against - so one going missing is a defect, never a sampling
        # artefact. Checked against the union of ALL trees regardless of source,
        # since a reference may sit in one segment's tree and not another's.
        #
        # This is the check that catches a test-mode subsample dropping
        # references. Measured across the reference databases: HCV_full 238/238,
        # rabv_update 28/28, IAV 396/396 all satisfy it.
        # Only the UShER tree is expected to hold every reference. UShER places
        # every eligible sequence, so a reference absent from it was dropped
        # upstream. IQ-TREE is built from MMseqs cluster representatives, so a
        # non-centroid reference is legitimately absent from it - measured,
        # rabv_update's IQ-TREE holds 13 of 28 references and HCV_full's holds
        # 219 of 238, while their UShER trees hold 28/28 and 238/238. Applying
        # this check to IQ-TREE would fail correct base_tree_only runs.
        references_missing_from_tree = []
        usher_tree_present = False
        if tree_available and reference_master_set:
            usher_terminals = set()
            for tr in fetch_trees(conn, source="usher"):
                try:
                    parsed = Phylo.read(StringIO(tr["newick"]), "newick")
                except Exception:
                    continue
                usher_tree_present = True
                usher_terminals.update(
                    tip.name for tip in parsed.get_terminals() if tip.name
                )
            if usher_tree_present:
                # Excluded rows are context only - an update run carries
                # reference/master rows marked exclusion_status=1 that are
                # deliberately not placed. Requiring them in the tree would fail
                # a correct update.
                references_missing_from_tree = sorted(
                    reference_master_set - usher_terminals - excluded_accessions
                )

        cluster_col = find_cluster_column(meta_columns)
        cluster_rows = []
        unclustered_accessions = set()
        centroid_set = set()
        centroid_set_by_segment = {}
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
            if args.segment_tree_source == "iqtree" and "segment" in meta_columns:
                cursor.execute(f"SELECT {args.accession_column}, segment, {cluster_col} FROM meta_data")
                for acc, seg, cluster_val in cursor.fetchall():
                    if is_cluster_placeholder(cluster_val):
                        continue
                    seg_str = str(seg).strip() if seg is not None and str(seg).strip() else None
                    if not seg_str:
                        continue
                    centroid_set_by_segment.setdefault(seg_str, set()).add(str(cluster_val).strip())
            if tree_available:
                missing_centroids_in_tree = sorted(centroid_set - tree_terminals)
                extra_in_tree = sorted(t for t in (tree_terminals - centroid_set) if t not in reference_master_set)

        if cluster_col:
            # An accession with no cluster assignment was never offered to the
            # clustering step, so it cannot appear in a tree built from cluster
            # output. Requiring it there is not a weaker check - it is the only
            # correct one.
            #
            # This is not hypothetical. In test mode TEST_SUBSAMPLE_CLUSTER_INPUT
            # caps the clustering input at params.test_max_cluster_seqs (120 for
            # the HCV_OM_test profile); everything past the cap is never
            # clustered, never placed, and can never reach the tree. Measured on
            # both a test DB and a 274,606-sequence production DB, every single
            # accession absent from the UShER tree was unclustered, and no
            # clustered accession was ever missing from it.
            unclustered_accessions = {
                str(acc).strip()
                for acc, cluster_val in cluster_rows
                if acc and str(acc).strip() and is_cluster_placeholder(cluster_val)
            }

        missing_in_tree = (
            sorted(expected_meta_set - tree_terminals - unclustered_accessions)
            if tree_available
            else []
        )

        consistency_results = {
            "sequence_alignment vs meta_data": validate_sequence_alignment_vs_meta(
                conn,
                expected_meta_set,
                accession_column=args.accession_column,
                exclusion_column=args.exclusion_column,
                exclude_value=args.exclude_value,
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
            ),
        }
        mutation_integrity_results = validate_mutation_tables(conn, expected_meta_set)
        if args.skip_invariant_checks:
            invariant_results = []
        else:
            invariant_results = run_extended_invariant_checks(
                conn,
                meta_columns,
                expected_meta_set,
                accession_column=args.accession_column,
                exclusion_column=args.exclusion_column,
                exclude_value=args.exclude_value,
                where_sql=where_sql,
                where_params=where_params,
            )
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
        cursor.execute("SELECT COUNT(*) FROM trees WHERE source=?", (args.segment_tree_source,))
        segment_tree_count = cursor.fetchone()[0]

        segment_count = None
        if "segment" in meta_columns:
            segment_values = [
                x for x in fetch_table_column(conn, "meta_data", "segment")
                if x is not None and str(x).strip() != ""
            ]
            segment_count = len(set(segment_values))

        segmented_validation_ok = False
        iqtree_segment_report_lines = []
        if tree_available and args.expect_segment_trees and segment_count is not None and segment_count > 1:
            segment_trees = fetch_trees(conn, source=args.segment_tree_source)
            segment_tree_label = "UShER" if args.segment_tree_source == "usher" else "IQ-TREE"

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
            for tr in segment_trees:
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
                if args.segment_tree_source == "iqtree" and centroid_set_by_segment:
                    expected_segment_accessions = centroid_set_by_segment.get(seg, set())
                else:
                    expected_segment_accessions = accessions_by_segment.get(seg, set())
                tree_terms = trees_by_segment.get(seg, set())
                missing_segment_accessions = sorted(expected_segment_accessions - tree_terms)
                if missing_segment_accessions:
                    per_segment_missing[seg] = missing_segment_accessions
                segment_reference_master_set = reference_master_set_by_segment.get(seg, set())
                extra_segment_accessions = sorted(t for t in (tree_terms - expected_segment_accessions) if t not in segment_reference_master_set)
                if extra_segment_accessions:
                    per_segment_extra[seg] = extra_segment_accessions
                if args.segment_tree_source == "iqtree":
                    iqtree_segment_report_lines.append(
                        f"[info] IQ-TREE segment {seg}: centroid_nodes_expected={len(expected_segment_accessions)} "
                        f"centroid_nodes_observed={len(tree_terms)} missing_centroids={len(missing_segment_accessions)} "
                        f"extra_tree_nodes={len(extra_segment_accessions)}"
                    )

            print(f"[info] Segmented validation: meta segments={len(segments_in_meta)} tree segments={len(trees_by_segment)}")
            if args.segment_tree_source == "iqtree" and centroid_set_by_segment:
                print("[info] IQ-TREE segmented validation: comparing per-segment centroid nodes only")
                for line in iqtree_segment_report_lines:
                    print(line)
            if missing_segment_trees:
                print(f"[info] Missing tree segments: {', '.join(missing_segment_trees[:10])}")
            if per_segment_missing:
                first_seg = sorted(per_segment_missing.keys())[0]
                print(f"[info] Segment {first_seg} missing accessions in tree (first 10): {', '.join(per_segment_missing[first_seg][:10])}")
            if per_segment_extra:
                first_seg = sorted(per_segment_extra.keys())[0]
                print(f"[info] Segment {first_seg} extra accessions in tree (first 10): {', '.join(per_segment_extra[first_seg][:10])}")

            if missing_segment_trees:
                validation_fail(
                    f"Validation failed: segmented run expects at least one {segment_tree_label} tree per segment "
                    f"(missing segments={', '.join(missing_segment_trees[:10])})"
                )
            if per_segment_missing and not args.test_mode:
                validation_fail(f"Validation failed: per-segment {segment_tree_label} tree is missing unfiltered accessions")
            if per_segment_missing and args.test_mode:
                print(f"[info] Test mode: allowing per-segment {segment_tree_label} trees to be subsampled relative to meta_data accessions")

            segmented_validation_ok = True

        report_path = os.path.join(args.outdir, "db_tree_validation.txt")
        missing_tree_path = os.path.join(args.outdir, "missing_in_tree.txt")
        missing_seq_path = os.path.join(args.outdir, "missing_in_sequences.txt")
        missing_meta_path = os.path.join(args.outdir, "missing_in_meta.txt")
        missing_centroids_path = os.path.join(args.outdir, "missing_centroids_in_tree.txt")
        extra_in_tree_path = os.path.join(args.outdir, "extra_in_tree.txt")
        extra_tree_label = get_extra_tree_label(tree_source)

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
            report.write(f"{('UShER' if args.segment_tree_source == 'usher' else 'IQ-TREE')} tree rows: {segment_tree_count}\n")
            if segment_count is not None:
                report.write(f"Segment count in meta_data: {segment_count}\n")
            if cluster_col:
                report.write(f"Cluster column: {cluster_col}\n")
                report.write(f"Centroid count: {len(centroid_set)}\n")
                report.write(f"Missing centroids in tree: {len(missing_centroids_in_tree)}\n")
                report.write(f"{extra_tree_label}: {len(extra_in_tree)}\n")
            if args.segment_tree_source == "iqtree" and cluster_col and segment_count is not None and segment_count > 1:
                report.write("\n=== IQ-TREE Segment Validation ===\n")
                report.write("Per-segment centroid-node checks only\n")
                for line in iqtree_segment_report_lines:
                    report.write(line + "\n")
            else:
                report.write(
                    f"References/master missing from the UShER tree: "
                    f"{len(references_missing_from_tree) if usher_tree_present else 'not checked (no UShER tree)'}\n"
                )
                if references_missing_from_tree:
                    report.write("\nReferences missing from the UShER tree (first 50):\n")
                    report.write("\n".join(references_missing_from_tree[:50]) + "\n")
                report.write(
                    f"Unclustered accessions (never offered to clustering, so not "
                    f"eligible for the tree): {len(unclustered_accessions)}\n"
                )
                report.write(
                    f"Missing in tree (clustered meta_data -> tree): {len(missing_in_tree)}\n"
                )
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
                    report.write(f"\n{extra_tree_label} (first 50):\n")
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
                        changed_rows = [row for row in rows if int(row[3]) != 0]
                        meta_delta = next((int(delta) for table_name, _before, _after, delta in rows if table_name == "meta_data"), 0)
                        report.write(f"- changed_tables: {len(changed_rows)}\n")
                        report.write(f"- sequences_added_to_db: {meta_delta}\n")
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

            if invariant_results:
                report.write("\nExtended DB invariant checks")
                report.write(" (fatal)\n" if args.strict_invariants else " (warnings only)\n")
                for result in invariant_results:
                    report.write(format_invariant_result(result))
                    report.write("\n\n")

            report.write("\nValidation status: PASS\n")

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
        if args.segment_tree_source == "iqtree" and cluster_col and segment_count is not None and segment_count > 1:
            print(
                f"[info] IQ-TREE segmented validation: metadata rows={meta_count}, sequences={seq_count}, tree terminals={len(tree_terminals)}, iqtree rows={segment_tree_count}"
            )
        else:
            print(
                f"[info] Sequence counts by DB storage/origin location:\n"
                f"[info]   - metadata entries (meta_data table): {meta_count} (includes query sequences, references, and excluded rows)\n"
                f"[info]   - raw sequence records (sequences table): {seq_count} (raw sequences for non-excluded query/reference entries)\n"
                f"[info]   - tree leaves (trees table terminals): {len(tree_terminals)}"
            )
            print(f"[info] {('UShER' if args.segment_tree_source == 'usher' else 'IQ-TREE')} tree rows: {segment_tree_count}")
        if segment_count is not None:
            print(f"[info] Segment count in meta_data: {segment_count}")
        print("[info] Running DB consistency checks across tables:")
        print("[info]   - sequence_alignment vs meta_data: verifying that every non-excluded accession in meta_data has an aligned sequence")
        print("[info]   - features vs meta_data: verifying that every non-excluded accession has projected coordinate features")
        print("[info]   - host_taxa vs meta_data: verifying that host taxonomy lookup records exist for each host_taxa_id in meta_data")
        for title, result in consistency_results.items():
            if result.get("skipped"):
                print(f"[info] {title}: skipped ({result.get('reason')})")
            elif "error" in result:
                print(f"[info] {title}: ERROR: {result['error']}")
            else:
                missing = result.get('missing', [])
                extra = result.get('extra', [])
                print(
                    f"[info] {title}: ok={result['ok']} missing={len(missing)} extra={len(extra)}"
                )
                if missing:
                    print(f"[info]   {title} missing (present in meta_data but missing in {result.get('table', 'table')}, first 10): {', '.join(missing[:10])}")
                    # Say WHAT the missing rows are, not just that they are
                    # missing. A CI failure on this check reported 22 bare
                    # accessions and nothing else; they turned out to be
                    # reference-list `exclusion_list` entries (influenza B
                    # references, which cannot align to influenza A) that had not
                    # been marked exclusion_status=1, so they were wrongly
                    # expected to have alignments. That breakdown would have
                    # named the cause immediately.
                    describe_missing_accessions(conn, missing, meta_columns)
                if extra:
                    print(f"[info]   {title} extra (present in {result.get('table', 'table')} but missing/excluded in meta_data, first 10): {', '.join(extra[:10])}")
        for result in invariant_results:
            if result.get("skipped"):
                print(f"[info] {result['title']}: skipped ({result.get('reason')})")
            elif "error" in result:
                print(f"[warn] {result['title']}: ERROR: {result['error']}")
            elif result.get("ok"):
                print(f"[info] {result['title']}: ok=True")
            else:
                print(f"[warn] {result['title']}: ok=False")
                for finding in result.get("findings", []):
                    print(f"[warn]   {finding['code']}: {finding['count']} - {finding['detail']}")
                    for example in finding.get("examples", [])[:10]:
                        print(f"[warn]     * {example}")
        for result in mutation_integrity_results:
            if result.get("skipped"):
                print(f"[info] {result['title']}: skipped ({result.get('reason')})")
            elif "error" in result:
                print(f"[info] {result['title']}: ERROR: {result['error']}")
            else:
                orphan_accessions = result.get('orphan_accessions', [])
                print(
                    f"[info] {result['title']}: ok={result['ok']} "
                    f"blank_required_rows={result.get('blank_required_rows', 0)} "
                    f"duplicate_keys={result.get('duplicate_keys', 0)} "
                    f"orphan_accessions={len(orphan_accessions)}"
                )
                if orphan_accessions:
                    print(f"[info]   {result['title']} orphan accessions (first 10): {', '.join(orphan_accessions[:10])}")
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
                    f"unresolved_master_rows={feature_integrity_result.get('unresolved_master_rows', 0)} "
                    f"mismatched_cds_rows={feature_integrity_result.get('mismatched_cds_rows', 0)}"
                )
                if not feature_integrity_result.get('ok'):
                    examples = feature_integrity_result.get('examples', [])
                    if examples:
                        print(f"[info] feature projection integrity: offending row details (up to 25):")
                        for ex in examples:
                            reason = ex.get('reason', 'unknown')
                            print(
                                f"[info]   accession={ex['accession']} master={ex['master_ref_accession']} "
                                f"product={ex['product']} aln=({ex['aln_start']},{ex['aln_end']}) "
                                f"cds=({ex['cds_start']},{ex['cds_end']}) reason={reason}"
                            )
        if cluster_col:
            if args.segment_tree_source == "iqtree" and segment_count is not None and segment_count > 1:
                print("[info] IQ-TREE segmented validation: per-segment centroid-node checks only")
                for line in iqtree_segment_report_lines:
                    print(line)
            else:
                print(f"[info] Cluster column: {cluster_col}, Centroid count: {len(centroid_set)}")
                print(f"[info] Missing centroids in tree: {len(missing_centroids_in_tree)}")
                print(f"[info] {extra_tree_label}: {len(extra_in_tree)}")
                if usher_tree_present:
                    print(f"[info] References/master missing from the UShER tree: "
                          f"{len(references_missing_from_tree)}")
                else:
                    print("[info] References/master vs tree: not checked "
                          "(no UShER tree; IQ-TREE holds centroids, not every reference)")
                if references_missing_from_tree:
                    print(f"[info] References missing (first 10): "
                          f"{', '.join(references_missing_from_tree[:10])}")
                print(f"[info] Unclustered accessions (not eligible for the tree): "
                      f"{len(unclustered_accessions)}")
                print(f"[info] Missing in tree: {len(missing_in_tree)} "
                      f"(clustered accessions only; unclustered ones are excluded by construction)")
                print(f"[info] Missing in sequences: {len(missing_in_sequences)}")
                print(f"[info] Missing in meta_data: {len(missing_in_meta)}")
                if missing_in_tree:
                    print(f"[info] Missing in tree (first 10): {', '.join(missing_in_tree[:10])}")
                if missing_in_sequences:
                    print(f"[info] Missing in sequences (first 10): {', '.join(missing_in_sequences[:10])}")
                if missing_in_meta:
                    print(f"[info] Missing in meta_data (first 10): {', '.join(missing_in_meta[:10])}")
                if missing_centroids_in_tree:
                    print(f"[info] Missing centroids in tree (first 10): {', '.join(missing_centroids_in_tree[:10])}")
                if extra_in_tree:
                    print(f"[info] {extra_tree_label} (first 10): {', '.join(extra_in_tree[:10])}")

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

        if args.test_mode and tree_source == args.segment_tree_source and disjoint_sets:
            print(
                f"[warn] Test mode: tree terminals and meta_data are disjoint for {args.segment_tree_source}; "
                "treating as no-placeable-queries scenario and not failing validation."
            )
            return

        if tree_available and args.expect_segment_trees and segment_count is not None and segment_count > 1 and not segmented_validation_ok:
            if segment_tree_count < segment_count:
                validation_fail(
                    f"Validation failed: segmented run expects at least one {('UShER' if args.segment_tree_source == 'usher' else 'IQ-TREE')} tree per segment "
                    f"(segments={segment_count}, trees={segment_tree_count})"
                )

        if args.check_update_integrity:
            if not table_exists(conn, "update_batches") or not table_exists(conn, "update_table_deltas"):
                validation_fail("Validation failed: update audit tables are missing")
            update_batches_cols = get_table_columns(conn, "update_batches")
            if "mode" in update_batches_cols:
                cursor.execute("SELECT batch_id, mode FROM update_batches ORDER BY rowid DESC LIMIT 1")
                row = cursor.fetchone()
                batch_id = row[0] if row else None
                mode = row[1] if row else None
            else:
                cursor.execute("SELECT batch_id FROM update_batches ORDER BY rowid DESC LIMIT 1")
                row = cursor.fetchone()
                batch_id = row[0] if row else None
                mode = "update"
            if batch_id:
                cursor.execute("SELECT table_name, delta FROM update_table_deltas WHERE batch_id=?", (batch_id,))
                delta_rows = cursor.fetchall()
                if len(delta_rows) == 0:
                    validation_fail("Validation failed: no update_table_deltas rows found for latest batch")
                if str(mode).strip().lower() == "update" and not args.allow_empty_update:
                    if all(int(delta) == 0 for _table_name, delta in delta_rows):
                        validation_fail("Validation failed: latest update batch made no DB changes")
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
                        validation_fail("Validation failed: segment contamination detected in features table")
            if result_failed(feature_integrity_result):
                validation_fail("Validation failed: feature projection integrity check failed")

        failed_invariant_checks = [result["title"] for result in invariant_results if result_failed(result)]
        if failed_invariant_checks and args.strict_invariants:
            validation_fail(
                "Validation failed: DB invariant checks failed for " + ", ".join(failed_invariant_checks)
            )

        failed_consistency_checks = [title for title, result in consistency_results.items() if result_failed(result)]
        failed_mutation_checks = [result["title"] for result in mutation_integrity_results if result_failed(result)]
        if failed_consistency_checks or failed_mutation_checks:
            if args.test_mode and args.expect_segment_trees and segment_count is not None and segment_count > 1 and not args.check_update_integrity:
                print("[info] Test mode: skipping strict DB consistency failures for segmented validation run")
            else:
                validation_fail(
                    "Validation failed: DB consistency checks failed for " + ", ".join(failed_consistency_checks + failed_mutation_checks)
                )

        if not tree_available:
            return

        if references_missing_from_tree:
            # References are the backbone - cluster seeds, tree topology, and the
            # anchors genotype calls are made against. A missing one is never
            # acceptable, in test mode or otherwise. Checked on every path, not
            # just the segment-tree branch.
            validation_fail(
                f"Validation failed: {len(references_missing_from_tree)} reference/master "
                f"accessions are absent from the UShER tree "
                f"(first 10: {', '.join(references_missing_from_tree[:10])})"
            )

        if tree_source == args.segment_tree_source:
            if segmented_validation_ok:
                return
            if missing_in_tree or missing_in_sequences or missing_in_meta:
                validation_fail(f"Validation failed: {('UShER' if tree_source == 'usher' else 'IQ-TREE')} tree does not match accessions")
            if cluster_col and extra_in_tree:
                print(
                    f"[info] {('UShER' if tree_source == 'usher' else 'IQ-TREE')} tree contains non-centroid leaves beyond the MMseqs centroid context "
                    f"({len(extra_in_tree)} leaves); full accession coverage passed."
                )
        elif cluster_col:
            if missing_centroids_in_tree or extra_in_tree:
                validation_fail("Validation failed: IQ-TREE does not match centroid set")
        else:
            if missing_in_tree or missing_in_sequences or missing_in_meta:
                validation_fail("Validation failed: missing accessions detected")

        print("[info] Validation status: PASS")

    finally:
        conn.close()


if __name__ == "__main__":
    main()
