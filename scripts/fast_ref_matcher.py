#!/usr/bin/env python3
"""High-throughput reference top-hit matcher for massive intraspecies viral builds (e.g. 14M+ SARS-CoV-2).

Why this exists:
  Running monolithic `blastn -num_threads 48` on 14.5M full-length 30 kb SARS-CoV-2
  queries against 5,836 99.8%-identical Pango references exhausts kernel VMA maps
  (`/proc/sys/vm/max_map_count = 65530`, crashing with `CThread::Run() -- error creating thread`)
  and requires >35 days of gapped seed extensions across all 5,836 near-identical references.

How this works (when all references align cleanly to the single master via `nextalign`):
  1. Aligns all references once to the master sequence with `nextalign` (~2.5s for 5,836 refs),
     building an exact sparse allele-difference inverted index + prefix-sum array + insertion
     signature table relative to the master.
  2. Streams `query_seq.fa` in 100,000-sequence batches:
       - Detects minus-strand queries via conserved master 25-mers and reverse-complements
         them prior to alignment (recording `strand = 'minus'` for downstream handling).
       - Aligns the batch to the master with `nextalign --jobs <threads>` (~20s per 100k seqs).
       - Computes exact pairwise mismatches + shared-insertion bonuses across ALL references
         simultaneously via forked multiprocessing workers (~2s per 100k seqs).
       - Routes any sequence that `nextalign` failed to align into `unaligned_fasta_out` so
         the caller runs standard chunked `blastn` fallback with zero loss of sensitivity.
"""

import csv
from concurrent.futures import ProcessPoolExecutor
import os
import shutil
import subprocess
import tempfile
import time
import numpy as np

_RC_TRANS = str.maketrans(
    "ACGTRYMKWSBDHVNacgtrymkwsbdhvn",
    "TGCAYRKMWSVHDBNtgcayrkmwsvhdbn",
)

# Shared read-only globals inherited via fork COW by worker processes
_G_REF_IDS = None
_G_MASTER = None
_G_REF_STARTS = None
_G_REF_ENDS = None
_G_PREF_DIFF = None
_G_DIFF_REFS = None
_G_PAD_OR_N_REFS = None
_G_SAME_ALLELE = None
_G_INS_REFS = None
_G_REF_INS_PENALTY = None


def reverse_complement(seq):
    return seq.translate(_RC_TRANS)[::-1]


def iter_fasta_blocks(path, chunk_bytes=64 * 1024 * 1024):
    """Yield (id_token, seq_upper) from a FASTA file using 64 MiB block reads."""
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        carry = ""
        first_chunk = True
        while True:
            chunk = fh.read(chunk_bytes)
            if not chunk:
                break
            data = carry + chunk
            if first_chunk:
                gt = data.find(">")
                if gt == -1:
                    carry = ""
                    continue
                data = data[gt + 1:]
                first_chunk = False
            parts = data.split("\n>")
            carry = parts[-1]
            for blk in parts[:-1]:
                hdr, sep_nl, body = blk.partition("\n")
                if sep_nl or hdr:
                    h = hdr.strip()
                    qid = h.split(None, 1)[0] if h else ""
                    if qid:
                        yield qid, body.replace("\n", "").replace("\r", "").replace(" ", "").upper()
        if carry and not first_chunk:
            hdr, sep_nl, body = carry.partition("\n")
            if sep_nl or hdr:
                h = hdr.strip()
                qid = h.split(None, 1)[0] if h else ""
                if qid:
                    yield qid, body.replace("\n", "").replace("\r", "").replace(" ", "").upper()


def _load_aligned_matrix(fasta_path):
    ids = []
    byte_rows = []
    expected_len = None
    for qid, seq in iter_fasta_blocks(fasta_path):
        b = seq.encode("ascii", errors="replace")
        if expected_len is None:
            expected_len = len(b)
        elif len(b) != expected_len:
            return None, None
        ids.append(qid)
        byte_rows.append(b)
    if not ids:
        return [], np.empty((0, 0), dtype=np.uint8)
    arr = np.frombuffer(b"".join(byte_rows), dtype=np.uint8).reshape(len(ids), expected_len)
    return ids, arr


def _load_insertions_csv(csv_path):
    """Return {seqName: tuple_of_(insertion_token, bonus_weight)} from nextalign .insertions.csv."""
    out = {}
    if not csv_path or not os.path.isfile(csv_path):
        return out
    with open(csv_path, "r", encoding="utf-8", errors="replace", newline="") as fh:
        reader = csv.reader(fh)
        header = next(reader, None)
        if not header:
            return out
        for row in reader:
            if len(row) < 2 or not row[1]:
                continue
            acc = row[0].strip().split()[0]
            tokens = []
            for tok in row[1].split(";"):
                tok = tok.strip()
                if not tok:
                    continue
                ins_seq = tok.split(":", 1)[-1] if ":" in tok else tok
                tokens.append((tok, max(4, len(ins_seq) * 2)))
            if tokens:
                out[acc] = tuple(tokens)
    return out


def _build_strand_probes(master_seq, k=25, n_probes=16):
    acgt = set("ACGT")
    probes_plus = []
    step = max(1, (len(master_seq) - k) // (n_probes + 2))
    for idx in range(1, n_probes + 1):
        start = idx * step
        kmer = master_seq[start:start + k]
        if len(kmer) == k and set(kmer) <= acgt:
            probes_plus.append(kmer)
    probes_minus = [reverse_complement(km) for km in probes_plus]
    return tuple(probes_plus), tuple(probes_minus)


def _score_worker_slice(args):
    q_ids, q_bytes_blob, width, q_strands, q_ins_list = args
    n_q = len(q_ids)
    if n_q == 0:
        return []
    q_arr = np.frombuffer(q_bytes_blob, dtype=np.uint8).reshape(n_q, width)
    M = _G_MASTER
    ref_ids = _G_REF_IDS
    ref_starts = _G_REF_STARTS
    ref_ends = _G_REF_ENDS
    pref_diff = _G_PREF_DIFF
    diff_refs_at_pos = _G_DIFF_REFS
    pad_or_n_at_pos = _G_PAD_OR_N_REFS
    same_allele_refs = _G_SAME_ALLELE
    ins_refs_map = _G_INS_REFS
    ref_ins_penalty = _G_REF_INS_PENALTY

    lines = []
    for qi in range(n_q):
        q = q_arr[qi]
        ng = np.where(q != 45)[0]
        if len(ng) == 0:
            continue
        qs = int(ng[0])
        qe = int(ng[-1])
        dist = (pref_diff[:, qe + 1] - pref_diff[:, qs]).astype(np.int32)
        offset = 0
        q_slice = q[qs:qe + 1]
        q_diffs = np.where(q_slice != M[qs:qe + 1])[0] + qs
        for c in q_diffs:
            b = int(q[c])
            if b == 78:  # 'N' in query: neutral across all references
                dr = diff_refs_at_pos[c]
                if len(dr):
                    dist[dr] -= 1
                continue
            offset += 1
            pn = pad_or_n_at_pos[c]
            if len(pn):
                dist[pn] -= 1
            dr = diff_refs_at_pos[c]
            if len(dr):
                dist[dr] -= 1
            sr = same_allele_refs.get((int(c), b))
            if sr is not None:
                dist[sr] -= 1

        tot_mm = dist + offset
        ov_start = np.maximum(qs, ref_starts)
        ov_end = np.minimum(qe, ref_ends)
        ov_len = np.maximum(1, ov_end - ov_start + 1)
        score = ov_len - (4 * tot_mm) - ref_ins_penalty

        q_ins = q_ins_list[qi]
        if q_ins:
            for tok, weight in q_ins:
                r_idxs = ins_refs_map.get(tok)
                if r_idxs is not None:
                    score[r_idxs] += weight + 2

        bi = int(np.argmax(score))
        mm_best = max(0, int(tot_mm[bi]))
        span_best = int(ov_len[bi])
        pident = max(0.0, min(100.0, 100.0 * (span_best - mm_best) / max(1, span_best)))
        lines.append(f"{q_ids[qi]}\t{ref_ids[bi]}\t{pident:.3f}\t{q_strands[qi]}\n")
    return lines


def try_fast_nextalign_tophits(
    query_fasta,
    db_fasta,
    master_acc,
    tophits_out,
    unaligned_fasta_out,
    threads=8,
    batch_size=250000,
):
    """Attempt fast exact master-anchored top-hit matching.

    Returns True when reference indexing succeeded (writing matched hits to
    `tophits_out` and any Nextalign-unaligned fallback queries to
    `unaligned_fasta_out`), or False when preconditions are not met so the
    caller should run full chunked `blastn` instead.
    """
    global _G_REF_IDS, _G_MASTER, _G_REF_STARTS, _G_REF_ENDS
    global _G_PREF_DIFF, _G_DIFF_REFS, _G_PAD_OR_N_REFS, _G_SAME_ALLELE, _G_INS_REFS, _G_REF_INS_PENALTY

    if shutil.which("nextalign") is None:
        return False

    ref_seqs = dict(iter_fasta_blocks(db_fasta))
    if not ref_seqs or master_acc not in ref_seqs:
        return False

    master_seq = ref_seqs[master_acc]
    probes_plus, probes_minus = _build_strand_probes(master_seq)

    with tempfile.TemporaryDirectory(prefix="fast_ref_match_") as tmpdir:
        master_fa = os.path.join(tmpdir, "master.fa")
        with open(master_fa, "w", encoding="utf-8") as mfh:
            mfh.write(f">{master_acc}\n{master_seq}\n")

        ref_out_dir = os.path.join(tmpdir, "ref_aln")
        os.makedirs(ref_out_dir, exist_ok=True)
        cmd_ref = [
            "nextalign", "run",
            "--jobs", str(max(1, int(threads))),
            "--input-ref", master_fa,
            "--output-all", ref_out_dir,
            "--output-basename", "refs",
            db_fasta,
        ]
        res = subprocess.run(cmd_ref, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        ref_aln_fa = os.path.join(ref_out_dir, "refs.aligned.fasta")
        if res.returncode != 0 or not os.path.isfile(ref_aln_fa):
            return False

        ref_ids, ref_arr = _load_aligned_matrix(ref_aln_fa)
        if ref_ids is None or len(ref_ids) != len(ref_seqs) or master_acc not in ref_ids:
            # Some references could not be aligned cleanly to the single master (divergent virus)
            return False

        t0 = time.time()
        master_idx = ref_ids.index(master_acc)
        M = ref_arr[master_idx].copy()
        width = len(M)

        ref_clean = ref_arr.copy()
        ref_starts = np.zeros(len(ref_ids), dtype=np.int32)
        ref_ends = np.zeros(len(ref_ids), dtype=np.int32)
        for ri in range(len(ref_ids)):
            ng = np.where(ref_clean[ri] != 45)[0]
            if len(ng) == 0:
                return False
            s, e = int(ng[0]), int(ng[-1])
            ref_starts[ri] = s
            ref_ends[ri] = e
            ref_clean[ri, :s] = 0
            ref_clean[ri, e + 1:] = 0

        # Reject fast path if any reference is <85% identical to master over its aligned span
        for ri in range(len(ref_ids)):
            s, e = int(ref_starts[ri]), int(ref_ends[ri])
            r_sl = ref_clean[ri, s:e + 1]
            m_sl = M[s:e + 1]
            comp = (r_sl != 78) & (r_sl != 45) & (m_sl != 78) & (m_sl != 45)
            n_comp = int(comp.sum())
            if n_comp < max(100, int(0.5 * width)):
                return False
            ident = float((r_sl[comp] == m_sl[comp]).sum()) / n_comp
            if ident < 0.85:
                return False

        diff_mask = (ref_clean != M[None, :]) & (ref_clean != 0) & (ref_clean != 78)
        pref_diff = np.zeros((len(ref_ids), width + 1), dtype=np.int16)
        np.cumsum(diff_mask, axis=1, out=pref_diff[:, 1:])

        diff_refs_at_pos = [np.where(diff_mask[:, c])[0].astype(np.int32) for c in range(width)]
        pad_or_n_at_pos = [np.where((ref_clean[:, c] == 0) | (ref_clean[:, c] == 78))[0].astype(np.int32) for c in range(width)]
        same_allele_refs = {}
        for c in range(width):
            drefs = diff_refs_at_pos[c]
            if len(drefs) == 0:
                continue
            vals = ref_clean[drefs, c]
            for b in (65, 67, 71, 84, 45):
                if b == int(M[c]):
                    continue
                sub = drefs[vals == b]
                if len(sub):
                    same_allele_refs[(c, b)] = sub.astype(np.int32)

        ref_ins_raw = _load_insertions_csv(os.path.join(ref_out_dir, "refs.insertions.csv"))
        ref_id_to_idx = {rid: idx for idx, rid in enumerate(ref_ids)}
        ins_tmp = {}
        ref_ins_penalty = np.zeros(len(ref_ids), dtype=np.int32)
        for rid, tok_list in ref_ins_raw.items():
            r_idx = ref_id_to_idx.get(rid)
            if r_idx is None:
                continue
            ref_ins_penalty[r_idx] = 2 * len(tok_list)
            for tok, _w in tok_list:
                ins_tmp.setdefault(tok, []).append(r_idx)
        ins_refs_map = {tok: np.array(idxs, dtype=np.int32) for tok, idxs in ins_tmp.items()}

        _G_REF_IDS = tuple(ref_ids)
        _G_MASTER = M
        _G_REF_STARTS = ref_starts
        _G_REF_ENDS = ref_ends
        _G_PREF_DIFF = pref_diff
        _G_DIFF_REFS = diff_refs_at_pos
        _G_PAD_OR_N_REFS = pad_or_n_at_pos
        _G_SAME_ALLELE = same_allele_refs
        _G_INS_REFS = ins_refs_map
        _G_REF_INS_PENALTY = ref_ins_penalty

        print(f"[fast-ref-match] Indexed {len(ref_ids):,} references ({width:,} cols, "
              f"{len(ins_refs_map):,} insertion signatures) in {time.time() - t0:.2f}s")

        score_workers = min(32, max(1, int(threads)))
        total_seen = 0
        total_matched = 0
        total_fallback = 0
        chunk_idx = 0

        with open(tophits_out, "w", encoding="utf-8") as out_hits, \
             open(unaligned_fasta_out, "w", encoding="utf-8") as out_fb, \
             ProcessPoolExecutor(max_workers=score_workers) as pool:

            batch_records = []

            def _process_batch(records):
                nonlocal total_matched, total_fallback, chunk_idx
                chunk_idx += 1
                b_t0 = time.time()
                batch_dir = os.path.join(tmpdir, f"batch_{chunk_idx}")
                os.makedirs(batch_dir, exist_ok=True)
                in_fa = os.path.join(batch_dir, "in.fa")
                strand_by_id = {}
                orig_seq_by_id = {}

                with open(in_fa, "w", encoding="utf-8") as bfh:
                    for qid, seq in records:
                        orig_seq_by_id[qid] = seq
                        plus_hits = sum(1 for km in probes_plus if km in seq)
                        minus_hits = 0
                        if plus_hits < 2 and probes_minus:
                            minus_hits = sum(1 for km in probes_minus if km in seq)
                        if minus_hits >= 2 and minus_hits > plus_hits:
                            strand_by_id[qid] = "minus"
                            bfh.write(f">{qid}\n{reverse_complement(seq)}\n")
                        else:
                            strand_by_id[qid] = "plus"
                            bfh.write(f">{qid}\n{seq}\n")

                cmd_q = [
                    "nextalign", "run",
                    "--jobs", str(max(1, int(threads))),
                    "--input-ref", master_fa,
                    "--output-all", batch_dir,
                    "--output-basename", "q",
                    in_fa,
                ]
                q_res = subprocess.run(cmd_q, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                q_aln_fa = os.path.join(batch_dir, "q.aligned.fasta")

                aligned_ids = []
                if q_res.returncode == 0 and os.path.isfile(q_aln_fa):
                    q_ins_map = _load_insertions_csv(os.path.join(batch_dir, "q.insertions.csv"))
                    aligned_rows = []
                    for qid, aseq in iter_fasta_blocks(q_aln_fa):
                        ab = aseq.encode("ascii", errors="replace")
                        if len(ab) == width and aseq.count("-") < width:
                            aligned_ids.append(qid)
                            aligned_rows.append(ab)

                    if aligned_ids:
                        n_aln = len(aligned_ids)
                        sub_size = max(250, (n_aln + score_workers - 1) // score_workers)
                        tasks = []
                        for start in range(0, n_aln, sub_size):
                            end = min(n_aln, start + sub_size)
                            sub_ids = aligned_ids[start:end]
                            blob = b"".join(aligned_rows[start:end])
                            sub_strands = [strand_by_id.get(x, "plus") for x in sub_ids]
                            sub_ins = [q_ins_map.get(x, ()) for x in sub_ids]
                            tasks.append((sub_ids, blob, width, sub_strands, sub_ins))

                        written_this_batch = set()
                        for result_lines in pool.map(_score_worker_slice, tasks):
                            for line in result_lines:
                                out_hits.write(line)
                                written_this_batch.add(line.split("\t", 1)[0])
                        total_matched += len(written_this_batch)
                        aligned_set = written_this_batch
                    else:
                        aligned_set = set()
                else:
                    aligned_set = set()

                for qid, _ in records:
                    if qid not in aligned_set:
                        out_fb.write(f">{qid}\n{orig_seq_by_id[qid]}\n")
                        total_fallback += 1

                shutil.rmtree(batch_dir, ignore_errors=True)
                elapsed = time.time() - b_t0
                print(
                    f"[fast-ref-match] Batch {chunk_idx}: {len(aligned_set):,}/{len(records):,} matched "
                    f"in {elapsed:.1f}s ({len(records) / max(0.01, elapsed):,.0f} seq/s; "
                    f"cumulative matched={total_matched:,}, blast_fallback={total_fallback:,})",
                    flush=True,
                )

            for qid, seq in iter_fasta_blocks(query_fasta):
                total_seen += 1
                batch_records.append((qid, seq))
                if len(batch_records) >= batch_size:
                    _process_batch(batch_records)
                    batch_records = []
            if batch_records:
                _process_batch(batch_records)

        return True
