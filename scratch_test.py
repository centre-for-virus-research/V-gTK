import sys
sys.path.insert(0, "/home3/oml4h/RABV-gTK/scripts")
from PadAlignment import PadAlignment
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

processor = PadAlignment(
    reference_alignment=None,
    input_dir=None,
    base_dir=None,
    output_dir=None,
    keep_intermediate_files=True
)

reference_aligned = "ACGT--ACGT"
subalignment = [
    SeqRecord(Seq("AC--ACGT"), id="REF1"),
    SeqRecord(Seq("AC--ACGT"), id="Q1")
]

ref_aligned_str = str(reference_aligned)
nextalign_ref_rec = subalignment[0]
nextalign_ref_seq = str(nextalign_ref_rec.seq)

master_dq_raw = ref_aligned_str.replace('-', '')
query_dq_raw = nextalign_ref_seq.replace('-', '')

print("master_dq_raw:", master_dq_raw)
print("query_dq_raw:", query_dq_raw)

import difflib
matcher = difflib.SequenceMatcher(None, master_dq_raw, query_dq_raw)
mapping = {i: None for i in range(len(master_dq_raw))}
for a, b, size in matcher.get_matching_blocks():
    print(f"Match block: a={a}, b={b}, size={size}")
    for offset in range(size):
        mapping[a + offset] = b + offset

print("mapping:", mapping)

nextalign_raw_to_col = {}
raw_idx = 0
for col_idx, char in enumerate(nextalign_ref_seq):
    if char != '-':
        nextalign_raw_to_col[raw_idx] = col_idx
        raw_idx += 1

print("nextalign_raw_to_col:", nextalign_raw_to_col)

updated = processor.insert_gaps(reference_aligned, subalignment, "REF1")
print("Q1 projected sequence:", str(updated[1].seq))
