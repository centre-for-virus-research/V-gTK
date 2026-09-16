"""Reading a master codon from the record, when the alignment misplaces it.

Every expected residue in the hand-written cases is read off the record by eye,
never from the module under test. The chaos tests at the end are the important
ones: frame_chaos builds records from a known history and damages their
alignments the way a frame-unaware aligner does, and a codon may come back
unknown, but never as a residue the record does not carry at that codon.
"""

import pytest

from AnnotateMutations import translate_codon
import frame_chaos as chaos
import reading_frame as rf


def frame_for(master, alignment, record, cds_start=1, cds_end=None, master_alignment=None):
    mapping = rf.RecordMapping(alignment, record)
    coord_map = {i + 1: i for i in range(len(master))}
    grid = rf.codon_grid(coord_map, cds_start, len(master) if cds_end is None else cds_end)
    frame = rf.ProteinFrame(mapping, grid, translate_codon, 'X',
                            master if master_alignment is None else master_alignment)
    return mapping, grid, frame


def residues(frame, grid, numbers=None):
    return {number: frame.read(grid[number]) for number in (numbers or sorted(grid))}


def test_placements_bound_every_valid_source_position():
    assert rf.placements('ACGT', 'ACGT') == ([0, 1, 2, 3], [0, 1, 2, 3])
    # The stripped G can only have come from one place.
    assert rf.placements('ACT', 'ACGT') == ([0, 1, 3], [0, 1, 3])
    # Inside a run, it could have been any of them.
    early, late = rf.placements('AA', 'AAA')
    assert (early, late) == ([0, 1], [1, 2])
    assert rf.placements('ACGT', 'ACG') is None


def test_record_stored_in_the_other_orientation_is_still_placed():
    mapping = rf.RecordMapping('ATGAAACCC', 'GGGTTTCAT')
    assert mapping.mapped and mapping.orientation == '-'
    assert mapping.record == 'ATGAAACCC'


def test_records_stored_as_rna_lower_case_or_wrapped_are_still_placed():
    for record in ['augaaaccc', 'ATGAAA\nCCC', 'ATG AAA CCC\r\n', 'GGGUUUCAU']:
        mapping = rf.RecordMapping('atgaaaccc', record)
        assert mapping.mapped, record
        assert mapping.record == 'ATGAAACCC'


def test_a_clean_record_reads_exactly_as_the_columns_do():
    master = 'ATGAAACCCGGGTTTTAA'
    _, grid, frame = frame_for(master, master, master)
    assert frame.status() == 'in_frame'
    assert not frame.has_issue
    assert residues(frame, grid) == {1: ('M', rf.IN_FRAME), 2: ('K', rf.IN_FRAME), 3: ('P', rf.IN_FRAME),
                                     4: ('G', rf.IN_FRAME), 5: ('F', rf.IN_FRAME), 6: ('*', rf.IN_FRAME)}


def test_compensated_artefact_is_reread_and_both_indels_are_named():
    """One base gapped at codon 5, one stripped at codon 2: codons 3-4 shift.

    The record reads M K S R V L * with no internal stop, so the shift is the
    aligner's and the right residues are still there to be read.
    """
    master = 'ATGAAACCCGGGTTTCTGTAA'
    alignment = 'ATGAAACCCGGGT-TCTGTAA'
    record = 'ATGAAATCCCGGGTTCTGTAA'
    _, grid, frame = frame_for(master, alignment, record)
    assert frame.status() == 'alignment_shift_corrected'
    assert (frame.frame, frame.off_frame, frame.split, frame.displaced) == (0, [3, 4], [], [3, 4])
    assert frame.internal_stops == 0
    read = residues(frame, grid)
    assert read[3] == ('S', rf.CORRECTED) and read[4] == ('R', rf.CORRECTED)
    assert read[1] == ('M', rf.IN_FRAME) and read[6] == ('L', rf.IN_FRAME)
    # 'T-T' in the columns, but the record has not lost a base between codons 2
    # and 6, so the gap is the artefact's and GTT is codon 5.
    assert read[5] == ('V', rf.CORRECTED)
    row = frame.scan_row()
    assert row['insertions'] == 'c2:1nt:stripped' and row['deletions'] == 'c5:1nt'
    assert (row['off_frame_codon_runs'], row['frameshift_indels'], row['unresolved_codons']) == ('3-4', 2, 0)


def test_a_two_base_shift_reads_the_codon_the_anchors_point_to():
    """The bug this module was rewritten for.

    Two gaps at codon 3 and two bases dropped at codon 5. Codon 4's columns hold
    record bases 7-9 ('CCG'); rounding their middle base to the frame picks
    record 6-8 (CCC, P). Codons 2 and 6 are in register, so codon 4 is record
    9-11: GGG, G.
    """
    master = 'ATGAAACCCGGGTTTCTGGCATAA'
    alignment = 'ATGAAA--CCCGGGTCTGGCATAA'
    _, grid, frame = frame_for(master, alignment, master)
    read = residues(frame, grid)
    assert read[3] == ('P', rf.CORRECTED)   # '--C' in the columns
    assert read[4] == ('G', rf.CORRECTED)   # 'CCG' in the columns
    assert read[5] == ('F', rf.CORRECTED)   # split by the dropped bases
    assert read[6] == ('L', rf.IN_FRAME)


def test_a_dropped_base_before_its_gap_shifts_the_other_way():
    master = 'ATGAAACATGGATTTCTGGCATAA'
    record = 'ATGAAACATGGATTTCTGGCATAA'
    alignment = 'ATGAAA' + 'ATGGA' + '-' + 'TTTCTGGCATAA'  # C at 6 dropped, gap at 11
    _, grid, frame = frame_for(master, alignment, record)
    read = residues(frame, grid)
    assert read[3] == ('H', rf.CORRECTED)   # columns say ATG, M
    assert read[4] == ('G', rf.CORRECTED)   # columns say GA-
    assert read[5] == ('F', rf.IN_FRAME)


def test_two_artefacts_summing_to_a_codon_leave_the_displaced_codons_unknown():
    """Three gaps at codon 3, three bases dropped at codon 8: codons 4-7 move a whole codon.

    They are in frame, so nothing about their phase looks wrong, and the columns
    read H G F L where the record says G F L A. A real codon deletion plus a
    real codon insertion would look identical, so the only safe answer is
    unknown - and the all-gap codon 3 is not a deletion either.
    """
    master = 'ATGAAACATGGATTTCTGGCAAGCTGGTAA'
    alignment = 'ATGAAA' + '---' + 'CATGGATTTCTG' + 'AGCTGGTAA'
    _, grid, frame = frame_for(master, alignment, master)
    # Codon 7's G could be record base 17 or 18 (CTG|GCA), and codon 8's A base
    # 20 or 21 (GCA|AGC), so neither is placed with certainty - but they sit
    # inside the displaced stretch all the same.
    assert frame.displaced == [4, 5, 6]
    read = residues(frame, grid)
    assert read[3] == ('X', rf.GUARD_AMBIGUOUS)
    assert all(read[n] == ('X', rf.GUARD_AMBIGUOUS) for n in (4, 5, 6, 7, 8))
    assert read[9] == ('W', rf.IN_FRAME)


def test_a_real_codon_deletion_is_still_left_to_the_alignment():
    master = 'ATGAAACATGGATTTCTGGCAAGCTGGTAA'
    record = 'ATGAAAGGATTTCTGGCAAGCTGGTAA'
    alignment = 'ATGAAA---GGATTTCTGGCAAGCTGGTAA'
    _, grid, frame = frame_for(master, alignment, record)
    assert frame.read(grid[3]) == (None, rf.DEFER)
    assert frame.read(grid[4]) == ('G', rf.IN_FRAME)


def test_a_shift_from_the_utr_across_most_of_a_short_protein_does_not_win_the_vote():
    """A gap in the 5' UTR and a base dropped in codon 6 shift codons 1-5 by one.

    Five of eight codons sit in the shifted phase, and the record reads cleanly
    in both phases, so a vote over the protein alone would call M H G F L as
    H A W I S. The flanks show the shift returning, so codons 1-6 are displaced
    and read from the anchors either side.
    """
    utr5, orf, utr3 = 'CGTACGGATCCGTAC', 'ATGCATGGATTTCTGGCAAGCTAA', 'GGTCCATGCAAGCTT'
    master = utr5 + orf + utr3
    alignment = master[:6] + '-' + master[6:31] + master[32:]
    _, grid, frame = frame_for(master, alignment, master, cds_start=16, cds_end=39)
    assert frame.frame == 0
    assert frame.displaced == [1, 2, 3, 4, 5]
    assert [frame.read(grid[n])[0] for n in range(1, 9)] == list('MHGFLAS*')


def test_an_even_split_is_settled_by_which_frame_translates_the_record():
    """Three codons each side of an uncompensated extra T. Frame 1 stops at once; frame 0 is clean.

    Codons 4-6 are off that frame with an anchor on one side only, so where the
    extra base really went - and whether codon 4 is TGG or GGG - is unknowable.
    """
    master = 'ATGAAACCCGGGTTTTAA'
    record = 'ATGAAACCCTGGGTTTTAA'
    _, grid, frame = frame_for(master, master, record)
    assert (frame.frame, frame.tie, frame.tie_broken_by_orf) == (0, False, True)
    assert frame.off_frame == [4, 5, 6]
    assert frame.read(grid[4]) == ('X', rf.GUARD_UNANCHORED)
    assert frame.read(grid[3]) == ('P', rf.IN_FRAME)


def test_an_unplaceable_base_inside_a_run_is_not_guessed_without_both_anchors():
    master = 'ATGCCCGGGAAAAAATTTTAA'
    record = 'ATGCCCGGGAAAAAAATTTTAA'
    _, grid, frame = frame_for(master, master, record)
    assert (frame.frame, frame.off_frame, frame.ambiguous_codons) == (0, [6, 7], 2)
    assert frame.read(grid[3]) == ('G', rf.IN_FRAME)
    assert frame.read(grid[4]) == ('X', rf.GUARD_UNANCHORED)
    assert frame.read(grid[6]) == ('X', rf.GUARD_UNANCHORED)


def test_a_record_whose_own_orf_is_broken_is_not_read_in_either_frame():
    """An uncompensated extra base: the record itself is frameshifted.

    Codons 3-6 outvote 1-2, but read in their frame the record stops at once,
    while the frame codons 1-2 sit in reads it cleanly. Nothing says which is
    the coding frame, so nothing is claimed.
    """
    master = 'ATGAAACCCGGGTTTTAA'
    record = 'ATGAAATCCCGGGTTTTAA'
    _, grid, frame = frame_for(master, master, record)
    assert (frame.frame, frame.off_frame, frame.frame_conflict) == (1, [1, 2], True)
    assert frame.status() == 'frame_undetermined'
    assert set(residues(frame, grid).values()) == {('X', rf.GUARD_UNDETERMINED)}


def test_without_a_record_every_codon_is_left_to_the_alignment():
    master = 'ATGAAACCCGGGTTTTAA'
    _, grid, frame = frame_for(master, master, '')
    assert frame.status() == 'record_unmapped'
    assert frame.read(grid[1]) == (None, rf.DEFER)


def test_columns_outside_the_protein_are_left_to_the_alignment():
    master = 'ATGAAACCCGGGTTTTAA'
    _, grid, frame = frame_for(master, master, master, cds_end=9)
    assert frame.read((9, 10, 11)) == (None, rf.DEFER)


def test_codon_grid_follows_the_masters_own_columns():
    master = 'AT-GAAACCC'
    coord_map = {position: column for position, column
                 in zip(range(1, 10), [i for i, base in enumerate(master) if base != '-'])}
    assert rf.codon_grid(coord_map, 1, 9) == {1: (0, 1, 3), 2: (4, 5, 6), 3: (7, 8, 9)}


def test_displaced_runs_take_the_shortest_artefact_out_first():
    # Two one-codon shifts (3, 6) around a two-codon correct stretch (4-5): the
    # shifts are shorter, so they go first and 4-5 rejoin the codons around them.
    starts = {1: 3, 2: 6, 3: 10, 4: 12, 5: 15, 6: 19, 7: 21, 8: 24}
    assert rf.displaced_numbers(starts) == {3, 6}
    # Swap the lengths: two-codon shifts (3-4, 6-7) around one correct codon (5).
    # Now codon 5 looks like the artefact. Once it is out, 3-7 is four codons
    # against two either side - not shorter than what it would be displaced
    # from - so it stays, and the frame is left to the vote to decide.
    starts = {1: 3, 2: 6, 3: 10, 4: 13, 5: 15, 6: 19, 7: 22, 8: 24, 9: 27}
    assert rf.displaced_numbers(starts) == {5}
    # Two short runs agreeing around a long one do not displace it.
    starts = {1: 5, 2: 8, **{n: 3 * n + 3 for n in range(3, 12)}, 12: 41, 13: 44}
    assert rf.displaced_numbers(starts) == set()


def test_one_artefact_and_two_real_indels_can_leave_identical_evidence():
    """The limit of what a row and its record can say, stated as a test.

    Truth A: the record is the master, and the aligner dropped its C at 6 and
    gapped column 11 - codons 3 and 4 are H and G. Truth B: a different master
    (T in column 11), a record with a real extra C and a real missing T, and a
    correct alignment - codon 3 is M. Row and record are identical; only the
    master's bases differ. The reader takes the explanation with fewer events,
    so under B it is wrong. Comparing against the master's bases could tell A
    from B, but would pull every call towards wild type, which is the one bias
    a resistance call cannot afford.
    """
    record = 'ATGAAACATGGATTTCTGGCATAA'
    row = 'ATGAAA' + 'ATGGA' + '-' + 'TTTCTGGCATAA'
    master_a, master_b = record, 'ATGAAAATGGATTTTCTGGCATAA'
    _, grid_a, frame_a = frame_for(master_a, row, record)
    _, grid_b, frame_b = frame_for(master_b, row, record)
    assert residues(frame_a, grid_a) == residues(frame_b, grid_b)
    assert frame_a.read(grid_a[3]) == ('H', rf.CORRECTED)  # right under A, wrong under B


# -- chaos ------------------------------------------------------------------
#
# frame_chaos builds each case from a known history, so truth never comes from
# the code under test. Two families of modes:
#
# DECIDABLE (and benign): every artefact has more clean sequence around it than
# it spans, away from real indels, record ends and other artefacts. Row and
# record then determine the answer, so any wrong residue at all is a bug.
#
# HOSTILE: artefacts overlap, run off record ends, sit on real indels. Some of
# that is undecidable in principle (see the test above): a real frameshift an
# artefact happens to cancel, two same-direction shifts around a short correct
# stretch, a whole-codon artefact at the edge of coverage. There the reader is
# held to never doing worse than the columns it replaces, by a wide margin.

SEEDS = 600


def read_case(case):
    mapping = rf.RecordMapping(case.row, case.record)
    assert mapping.mapped, case.describe()
    frame = rf.ProteinFrame(mapping, case.grid, translate_codon, 'X', case.master)
    for number, columns in case.grid.items():
        residue, verdict = frame.read(columns)
        column = chaos.column_reading(case.row, columns)
        yield number, column if verdict == rf.DEFER else residue, verdict, column


def wrong(case, residue, number):
    return residue != 'X' and residue not in case.truth[number]


@pytest.mark.parametrize('mode', ('benign',) + chaos.DECIDABLE)
def test_no_codon_is_ever_read_as_a_residue_the_record_does_not_carry(mode):
    failures = []
    for seed in range(SEEDS):
        case = chaos.build(seed, mode)
        bad = [(number, residue, verdict, sorted(case.truth[number]))
               for number, residue, verdict, _ in read_case(case) if wrong(case, residue, number)]
        if bad:
            failures.append(f'{case.describe()}\n  wrong (codon, read, verdict, truth): {bad}')
    assert not failures, f'{len(failures)} of {SEEDS} cases read a wrong residue:\n\n' + '\n\n'.join(failures[:5])


@pytest.mark.parametrize('mode', chaos.HOSTILE)
def test_under_undecidable_chaos_the_reader_still_never_does_worse_than_the_columns(mode):
    codons = introduced = reader_wrong = column_wrong = 0
    examples = []
    for seed in range(SEEDS):
        case = chaos.build(seed, mode)
        for number, residue, verdict, column in read_case(case):
            codons += 1
            column_wrong += wrong(case, column, number)
            if wrong(case, residue, number):
                reader_wrong += 1
                if residue != column:
                    introduced += 1
                    examples.append((seed, number, verdict, column, residue, sorted(case.truth[number])))
    # Measured when written: introduced 0-25 per ~25,000 codons, and the reader
    # left 10-30% of the columns' wrong residues (the rest became right or X).
    assert introduced <= 0.0015 * codons, (mode, introduced, codons, examples[:10])
    assert reader_wrong <= 0.35 * column_wrong, (mode, reader_wrong, column_wrong)


def test_isolated_artefacts_are_always_recovered_not_just_hidden():
    """Unknown is safe but useless if it is all the reader ever says.

    Benign cases carry only one- and two-base artefacts with clean codons
    around them, so every codon has agreeing anchors and must be read exactly.
    """
    failures, corrected = [], 0
    for seed in range(SEEDS):
        case = chaos.build(seed, 'benign')
        for number, residue, verdict, _ in read_case(case):
            (truth,) = case.truth[number]
            corrected += verdict == rf.CORRECTED
            if residue != truth:
                failures.append((seed, number, residue, verdict, truth))
    assert not failures, failures[:10]
    assert corrected > SEEDS  # the artefacts really were there to correct


@pytest.mark.parametrize('mode,floor', [
    ('decidable_chaos', 0.9), ('decidable_repeat_hell', 0.85), ('decidable_majority_shift', 0.85),
    # A real uncompensated frameshift leaves no frame the whole protein can be
    # read in, so it is mostly unknown by design.
    ('decidable_frameshift', 0.5),
])
def test_most_readable_codons_are_still_read_where_the_data_can_decide(mode, floor):
    readable = read = 0
    for seed in range(SEEDS):
        case = chaos.build(seed, mode)
        for number, residue, _, _ in read_case(case):
            if len(case.truth[number] - {'X'}) == 1:
                readable += 1
                read += residue in case.truth[number]
    assert read / readable > floor, (mode, read, readable)
