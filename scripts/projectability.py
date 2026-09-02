#!/usr/bin/env python3
"""Single source of truth for "can PadAlignment project this reference?".

Two scripts need to answer that question and they used to answer it differently:

* :mod:`PadAlignment` answers it by *doing* it - for each master it resolves a
  per-segment backbone MSA (``refset_<segment>_aln.fasta`` from
  ``--precomputed_ref_dir``, or a backbone exported from ``--update_db``) and
  projects every ``query_aln/<ref>/`` sub-alignment whose ``<ref>`` is a row of
  that backbone.
* :mod:`CollectFilteredSequences` answered it by reading
  ``Nextalign/reference_aln/*/*.aligned.fasta`` - the output of aligning all 418
  references against the 8 *master* sequences - and deleting the queries of every
  reference missing from it.

Those are not the same set, and the second one is catastrophically smaller for a
segment whose references span deep antigenic diversity. Influenza segment 4's
master is ``AB573800`` (H1N1) and segment 6's is ``AB472016`` (H5N2); nextalign
cannot seed-match an H3/H7/H9 HA against an H1 master, so only 30 of 117 HA and
29 of 83 NA references survived ``reference_aln``. Their queries - 633,988
sequences, including 242,686 of the 242,732 H3 - were written to
``filtered_sequences_ids.txt`` and skipped by ``PadAlignment --skip_ids``, even
though ``refset_4_aln.fasta`` holds all 117 HA references and projects every one
of them correctly. The database shipped with 0.02% of its H3.

This module holds the predicate once so the two callers cannot drift again. It
mirrors :meth:`PadAlignment.PadAlignment.process_all_masters`'s resolution order
exactly; if that order ever changes, change it here.

Deliberate design notes
-----------------------
* The reference set is read from the **backbone FASTA headers**, never from
  ``ref_list_refmast.txt``. 22 accessions present in the influenza refsets
  (``NC_0073xx``, ``EU6366xx``, ``CY015122``, ``GU052512``) are absent from that
  list, and its 22 ``exclusion_list`` rows carry ``B``/``C``/``D`` in the segment
  column (influenza B/C/D), which is not a segment number. The ref list is used
  only to learn *which masters exist and which segment each one covers*.
* ``generic/influenza/ref_set_aligned/refset_*_aln.fasta`` are CRLF-terminated.
  Header parsing strips ``\\r`` explicitly - a bare ``line[1:].rstrip("\\n")``
  yields ids with a trailing carriage return that silently match nothing.
* Absence of a precomputed/DB backbone is not an error here. The non-segmented
  (RABV) and HCV paths pass no ``--precomputed_ref_dir``, and for them the
  ``reference_aln`` fallback *is* the correct answer, because they are
  single-master and their references are near-identical to that master.
"""

import os
import re

import accession_utils
import segment_utils


#: Spellings of "no value" that ``nextflow.config`` and ``vgtk-init.nf`` use for
#: an unset path parameter before it reaches argparse.
_UNSET_TOKENS = frozenset({"", "null", "none", "unset"})


def normalise_optional_path(path_value):
    """Return ``path_value`` as a real path, or ``None`` for an unset sentinel.

    Mirrors :meth:`PadAlignment.PadAlignment._normalize_optional_path`. Nextflow
    interpolates a missing param as the literal string ``null`` and
    ``vgtk-init.nf`` normalises some to ``UNSET``; both must read as absent.
    """
    if path_value is None:
        return None
    text = str(path_value).strip()
    if text.casefold() in _UNSET_TOKENS:
        return None
    return text


def read_fasta_ids(fasta_path):
    """Return the canonical accession of every record in ``fasta_path``.

    Headers are normalised the way the rest of the pipeline normalises them:
    first whitespace-delimited token, then :mod:`accession_utils`. CRLF-safe.
    Returns an empty set rather than raising if the file cannot be read - the
    caller decides whether an empty backbone is fatal.

    The version used to come off with ``.split(".")[0]``, which cuts at the
    first dot and so truncates any header id that is not a plain accession -
    exactly the ids this set is then membership-tested against by
    ``CollectFilteredSequences`` (a ``query_aln/<ref>`` directory name) and by
    ``PadAlignment.find_orphan_references``. A backbone row this set fails to
    spell the way its query directory is spelled marks every query in that
    group unprojectable, which is the shape of the incident this module was
    written to prevent.
    """
    ids = set()
    if not fasta_path or not os.path.isfile(fasta_path):
        return ids
    try:
        with open(fasta_path, "r", encoding="utf-8", errors="replace") as handle:
            for line in handle:
                if not line.startswith(">"):
                    continue
                token = line[1:].strip().split()
                if not token:
                    continue
                ids.add(accession_utils.normalise_accession(token[0]))
    except OSError as exc:
        print(f"[warn] Could not read backbone alignment {fasta_path}: {exc}")
    return ids


def find_precomputed_reference_alignment(precomputed_ref_dir, segment_value):
    """Resolve the backbone MSA for ``segment_value`` inside ``precomputed_ref_dir``.

    Canonical name is ``refset_<segment>_aln.fasta``. Falls back to a
    single-digit-run filename match for custom naming, and - only when no segment
    is known at all - to the sole FASTA in the directory. Returns ``None`` when
    nothing resolves, which is the caller's signal to fall back to
    ``reference_aln``.

    This is the function :meth:`PadAlignment.PadAlignment.find_precomputed_reference_alignment`
    delegates to, so the two can never disagree about which file gets opened.
    """
    if not precomputed_ref_dir or not os.path.isdir(precomputed_ref_dir):
        return None

    segment = segment_utils.normalise_segment(segment_value)
    if segment is not None and segment.casefold() in segment_utils.PANDAS_NULL_TOKENS:
        segment = None

    if not segment:
        fasta_files = sorted(
            os.path.join(precomputed_ref_dir, fname)
            for fname in os.listdir(precomputed_ref_dir)
            if fname.lower().endswith((".fasta", ".fa"))
        )
        if len(fasta_files) == 1:
            print(
                f"[warn] No segment value supplied; using sole precomputed alignment "
                f"{os.path.basename(fasta_files[0])}."
            )
            return fasta_files[0]
        return None

    preferred = os.path.join(precomputed_ref_dir, f"refset_{segment}_aln.fasta")
    if os.path.exists(preferred):
        return preferred

    fallback_matches = []
    for fname in os.listdir(precomputed_ref_dir):
        if not fname.lower().endswith((".fasta", ".fa")):
            continue
        seg_match = re.search(r"(\d+)", fname)
        if seg_match and seg_match.group(1) == segment:
            fallback_matches.append(os.path.join(precomputed_ref_dir, fname))

    if len(fallback_matches) == 1:
        return fallback_matches[0]
    if len(fallback_matches) > 1:
        chosen = sorted(fallback_matches)[0]
        print(
            f"[warn] Multiple precomputed alignment files matched segment {segment}: "
            f"{', '.join(os.path.basename(x) for x in sorted(fallback_matches))}. "
            f"Using {os.path.basename(chosen)}."
        )
        return chosen

    return None


def backbone_files(nextalign_dir=None, precomputed_ref_dir=None, master_segment_map=None):
    """Return the backbone MSA paths PadAlignment will actually project through.

    Resolution order per master, matching
    :meth:`PadAlignment.PadAlignment.process_all_masters`:

    1. ``precomputed_ref_dir/refset_<segment>_aln.fasta`` for that master's segment
    2. ``nextalign_dir/reference_aln/<master>/<master>.aligned.fasta``

    When ``master_segment_map`` is empty (non-segmented viruses, or a ref list
    with no ``master`` rows) every ``reference_aln`` output is taken, which is
    what the single-master RABV/HCV builds want.

    Returns ``(paths, source)`` where ``source`` is ``"precomputed"``,
    ``"reference_aln"`` or ``"mixed"``, for logging.
    """
    precomputed_ref_dir = normalise_optional_path(precomputed_ref_dir)
    nextalign_dir = normalise_optional_path(nextalign_dir)

    paths = []
    used_precomputed = False
    used_reference_aln = False

    if master_segment_map:
        for master, segment in sorted(master_segment_map.items()):
            resolved = find_precomputed_reference_alignment(precomputed_ref_dir, segment)
            if resolved:
                paths.append(resolved)
                used_precomputed = True
                continue
            if nextalign_dir:
                fallback = os.path.join(
                    nextalign_dir, "reference_aln", master, f"{master}.aligned.fasta"
                )
                if os.path.isfile(fallback):
                    paths.append(fallback)
                    used_reference_aln = True
    elif precomputed_ref_dir and os.path.isdir(precomputed_ref_dir):
        # No master->segment map, but a backbone directory exists: every FASTA in
        # it is a backbone PadAlignment could open.
        for fname in sorted(os.listdir(precomputed_ref_dir)):
            if fname.lower().endswith((".fasta", ".fa")):
                paths.append(os.path.join(precomputed_ref_dir, fname))
                used_precomputed = True

    if not paths and nextalign_dir:
        ref_aln_root = os.path.join(nextalign_dir, "reference_aln")
        if os.path.isdir(ref_aln_root):
            for master in sorted(os.listdir(ref_aln_root)):
                candidate = os.path.join(ref_aln_root, master, f"{master}.aligned.fasta")
                if os.path.isfile(candidate):
                    paths.append(candidate)
                    used_reference_aln = True

    if used_precomputed and used_reference_aln:
        source = "mixed"
    elif used_precomputed:
        source = "precomputed"
    else:
        source = "reference_aln"

    # dict.fromkeys keeps insertion order while de-duplicating; two masters on the
    # same segment legitimately resolve to the same backbone file.
    return list(dict.fromkeys(paths)), source


def projectable_reference_ids(nextalign_dir=None, precomputed_ref_dir=None,
                              master_segment_map=None):
    """Return the set of reference accessions PadAlignment can project.

    A ``query_aln/<ref>/`` sub-alignment reaches the merged MSA if and only if
    ``<ref>`` is a row of one of the backbones returned by :func:`backbone_files`.
    Anything else is genuinely unprojectable and its queries are legitimately
    dropped.

    Returns ``(ids, source)``. An empty ``ids`` with a non-empty backbone list is
    a real failure and callers must treat it as one - it is exactly the state
    that silently deleted 16% of the influenza database.
    """
    paths, source = backbone_files(
        nextalign_dir=nextalign_dir,
        precomputed_ref_dir=precomputed_ref_dir,
        master_segment_map=master_segment_map,
    )
    ids = set()
    for path in paths:
        ids |= read_fasta_ids(path)
    return ids, source


def load_master_segment_map(ref_list_path):
    """Return ``{master_accession: segment}`` from a ref-list TSV.

    Used only to decide which backbone file to open per master - never as the
    reference set itself. Returns ``{}`` for a missing/2-column/unsegmented ref
    list, which callers must read as "not segmented", not as "no masters".
    """
    ref_list_path = normalise_optional_path(ref_list_path)
    if not ref_list_path or not os.path.isfile(ref_list_path):
        return {}
    try:
        from ExportRefListFromUpdateDb import load_reference_file_table
        df = load_reference_file_table(ref_list_path)
    except Exception as exc:  # pragma: no cover - defensive, matches PadAlignment
        print(f"[warn] Could not parse ref list {ref_list_path}: {exc}")
        return {}

    if "segment" not in df.columns or "accession_type" not in df.columns:
        return {}
    if not df["segment"].astype(str).str.strip().ne("").any():
        return {}

    masters = df[df["accession_type"].astype(str).str.strip().str.lower() == "master"]
    mapping = {}
    for _, row in masters.iterrows():
        acc = str(row["primary_accession"]).strip()
        segment = segment_utils.normalise_segment(row.get("segment"))
        if acc:
            mapping[acc] = segment
    return mapping
