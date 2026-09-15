# Coordinate spaces in feature projection

How a master's annotated CDS becomes a row in `features`, which number lines are
in play, and the two defects that came from mixing them. Read this before
changing `CalcAlignmentCord.py`, `CalcGenomeCords.py`, or anything that consumes
the `features` table.

---

## The three spaces

| # | space | what a number means | who speaks it |
|---|---|---|---|
| 1 | **record coordinates** | position along the master's own submitted record, 1-based | the master's GFF; `region`, `CDS` lines |
| 2 | **alignment columns** | column number in the merged MSA, 1-based | `get_gap_ranges`, `count_gaps_before_position`, every row of one alignment file |
| 3 | **row coordinates** | position along one row's own residues, gaps removed | `cds_start_OG_seq` / `cds_end_OG_seq`, and `sequence_alignment.alignment` degapped |

Space 1 and space 2 coincide **only** when the master's alignment row is both
ungapped and flush with its record. That is true for RABV and HCV and false for
seven of influenza's eight segments, which is why every defect below was
invisible until a segmented build was examined.

## The fourth thing, which is not a space: the 5' trim

The master's row in the merged MSA comes from the guide alignment
(`ref_set_aligned`), and a curated reference there is frequently **trimmed** at
the 5' end. Measured on the influenza asset:

| segment | master | gap columns | 5' trim |
|---|---|---|---|
| 1 | KJ889313 | 3 | 27 |
| 2 | CY087822 | 3 | 10 |
| 3 | MW333949 | 6 | 24 |
| 4 | AB573800 | 132 | 69 |
| 5 | NC_007369 | 0 | 45 |
| 6 | AB472016 | 102 | 0 |
| 7 | NC_002016 | 9 | 25 |
| 8 | KY243296 | 6 | 26 |

The GFF numbers bases from the start of the **record**; counting residues down
the alignment row numbers them from the start of the **row**. Those differ by the
trim. `compute_trim_offset` locates the row inside the record and the coordinate
map is keyed by `residue_index + trim`; the covered span from
`CalculateGenomeCoordinates` is shifted by the same amount, because it counts
master residues down the row too and the span clamp must compare like with like.

Locating the trim needs the master's record, which arrives as `--master_seq_dir`
(`BLAST_ALIGNMENT`'s `master_seq`). Without it the trim is **reported** and not
corrected — it cannot be inferred from lengths alone, because a length mismatch
does not say whether the missing bases came off the 5' or the 3' end.

## What each column holds

| column | space |
|---|---|
| `aln_start` / `aln_end` | 2 — the row's covered span, first and last column carrying a residue |
| `cds_start` / `cds_end` | 2 — the feature's alignment columns |
| `cds_start_OG_seq` / `cds_end_OG_seq` | 3 — the feature in this row's own numbering |

`aln_*` and `cds_*` are deliberately the same space:
`ValidateDbTree.validate_feature_projection_integrity` clamps the master's
`cds_*` by a row's `aln_*` and compares, which is only meaningful if both are
columns.

`cds_*_OG_seq` is **not** the record coordinate whenever the trim is non-zero. It
indexes the alignment row as stored in `sequence_alignment`, degapped. Slice the
raw `sequences` row with it and you will be wrong by the trim — and, because
nextalign also strips insertions relative to the reference, the degapped
alignment is frequently not even a substring of the submitted record (430 of 496
influenza rows are a substring; 50 are not).

## The two defects this encodes

**Gap counts taken against a record coordinate.** `count_gaps_before_position`
measures gaps in columns, and was being handed space-1 coordinates. A no-op for
an ungapped master; for a gapped one it shifted every OG coordinate upstream by
the master's gap count. Those counts are all multiples of three, so the reading
frame survived and the only symptom was a CDS starting on a UTR codon — an
internal stop whenever that codon happened to be one. Influenza PB1 read 97.6%
internal stops; NP and PA, whose masters have no gaps, read 1.5% and 4.3%.

**Record coordinates counted down a trimmed row.** With the trim uncorrected,
`map[gff_cds_start]` landed on the start codon for segment 6 alone — the only
flush master. The rows then mixed two spaces: `og_start` was the record
coordinate, surviving as a pass-through because the conversion degenerates when
the master's gaps all sit downstream, while `og_end` was clamped to the row's
length. HA translated as a correct-looking protein missing its C-terminus,
truncated by exactly its 69 nt trim.

After both fixes, on the `segmented_test` build: the OG reading and the
alignment-column reading name the identical sequence for **756 of 756** feature
rows (previously 1 of 670), no row carries `og_start < 1` (previously 86), and
PB2, PB1, PA, NP, NA and M1 are 100% ATG-initiated with no internal stops.

## Clamping

A feature start clamps **forwards** (`count_gaps_strictly_before`) to the first
residue at or after its column; a feature end clamps **backwards**
(`count_gaps_before_position`) to the last residue at or before it. Sharing the
inclusive count for both put a start that landed inside a gap one or more
residues *before* the feature, and could return `0`, which is not a valid 1-based
coordinate.

## Known limitations

Each has a standing test in `tests/unit/test_exotic_virus_genomes.py`.

- **Spliced CDS.** Two GFF lines for one product become two independent feature
  rows; nothing records the joined coding length. Influenza M2, NEP and PA-X are
  already in this shape, which is why they sit at exactly 50% ATG — exon 1 starts
  with ATG and exon 2 does not. `AnnotateMutations.choose_feature_entry` breaks a
  same-priority tie by preferring the **shortest** span, so M2 resolves to its
  26 nt first exon.
- **Circular genomes.** A CDS wrapping the origin has GFF start > end and is
  dropped by the same branch that discards features outside a partial sequence's
  span — silently, and with no use made of `meta_data.topology`.
- **Ambisense genomes.** `GffToDictionary` parses the strand column and discards
  it, so a minus-strand gene is projected like a plus-strand one and translates
  to the wrong protein without any error.
- **A CDS starting before the guide row.** Influenza HA's curated row begins 69
  bases into its CDS, past the start codon, so no projection can recover the
  initiator. This is an asset property, not an arithmetic one.

## Where this is tested

| module | covers |
|---|---|
| `test_alignment_coordinate_integrity.py` | the projection arithmetic alone, including fuzz |
| `test_alignment_coordinate_contract.py` | the whole script over files, and the contract with `AnnotateMutations` |
| `test_segmented_worst_case.py` | eight segments, eight trims, eight widths, in one pass |
| `test_exotic_virus_genomes.py` | circular, spliced, ambisense, and degenerate annotations |
| `test_update_adversarial.py` | what a coordinate correction does to an existing database |
