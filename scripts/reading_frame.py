"""Read a master codon in the submitted record's own reading frame.

Why this exists: nextalign aligns bases, not codons, and the pipeline runs it
with no gene map. It can place a gap and a stripped insertion a few bases apart,
both out of frame, and then the three alignment columns of one master codon hold
record bases that are not one codon. Translating those columns invents a residue
the record does not carry - measured on the shipped HCV build, 2,932 published
calls.

Coordinates do not move. A codon is still located by the master's alignment
columns, exactly as before; only the bases read out of it come from the record.

Terms used below:

record
    the submitted sequence (``sequences.sequence``), reverse-complemented if
    that is the orientation the alignment row was built from.
placement
    where each aligned base sits in the record. The leftmost and rightmost
    placements bound every valid one; where they differ, a base has more than
    one possible source and nothing is claimed from it.
offset
    record position minus three times the master codon number, for a codon whose
    three bases sit side by side at one certain place. Codons with the same
    offset are in register with each other: no base was gained or lost between
    them.
displaced run
    certainly placed codons in register with each other, with codons either
    side that agree with each other but not with them. That is the footprint
    of an aligner artefact - a gap and a dropped base a few codons apart - or,
    much more rarely, of a real insertion and deletion close together. Nothing
    tells those apart, so a displaced codon is never read where the alignment
    put it. Runs are judged over the protein and FLANK_CODONS of alignment
    either side, so a shift that starts outside the protein is still seen
    returning. The shortest displaced runs are taken out first and the rest
    re-judged, so a correct stretch between two artefacts is not mistaken for
    the artefact.
frame
    the phase (offset mod 3) of a protein's record codons: the *majority* phase
    of its certainly placed codons that are not displaced. A majority, because
    one local misalignment must not redefine the frame of every correct codon
    around it (anchoring on the first readable codon instead mis-flagged 20 HCV
    sites that GenBank confirms are in frame); displaced codons left out,
    because a shift covering most of a short protein would otherwise win the
    vote. A tie is broken by which candidate frame translates the record with no
    internal stops, and is otherwise undetermined.
anchor
    an in-frame, certainly placed, not displaced codon. Everything else is read
    *between two anchors*: if the left and right anchors share an offset, no
    base was gained or lost between them and the codon's record position is
    exact. If they differ (by a whole number of codons, since both are in frame)
    every in-frame position between the two offsets is tried, and the codon is
    only read if they all translate to the same residue.

    Why not just round an off-frame placement to the nearest in-frame codon:
    a stretch shifted by two bases has its middle base in the *neighbouring*
    codon, so rounding reads the wrong residue in the record and labels it
    corrected. Only the anchors know which way the shift went.
guard
    what happens when the record cannot answer: the residue is the caller's
    unknown token, never a guess. Guards fire for an undetermined frame, for a
    codon with no anchor on one side, for candidate positions that translate
    differently, and for a record whose own ORF is broken in the chosen frame -
    a real frameshift, where no in-frame codon is the homologue of the master's.

The module takes a ``translate`` callable and an ``unknown`` token from its
caller rather than importing them, so each caller keeps one residue vocabulary
and this file stays free of import cycles.
"""

from bisect import bisect_left
from collections import Counter
import re

#: Enough to reverse-complement a stored record, including the IUPAC codes
#: GenBank submissions carry. Anything else maps to N rather than raising.
_COMPLEMENT = str.maketrans('ACGTURYSWKMBDHVN-', 'TGCAAYRSWMKVHDBN-')

#: Master codons either side of a protein searched for displaced runs and anchors.
FLANK_CODONS = 30

#: Codons a rival phase needs before a clean read in it contradicts a clean
#: majority. One or two stray placements in a repeat are not a rival frame.
MIN_RIVAL_CODONS = 3


def reverse_complement(text):
    return text.translate(_COMPLEMENT)[::-1]


def normalise_bases(text):
    """Upper case, RNA read as DNA, whitespace dropped: one spelling to compare."""
    return re.sub(r'\s+', '', str(text or '')).upper().replace('U', 'T')


def placements(ungapped, record):
    """``(leftmost, rightmost)`` record index of each base of ``ungapped``.

    ``None`` when the aligned bases are not a subsequence of the record at all,
    which is the caller's signal that this row cannot be checked. Where the two
    lists agree at a position, that base has exactly one possible source.
    """
    if ungapped == record:
        shared = list(range(len(record)))
        return shared, shared
    early, offset = [], 0
    for base in ungapped:
        offset = record.find(base, offset)
        if offset == -1:
            return None
        early.append(offset)
        offset += 1
    late, offset = [], len(record)
    for base in reversed(ungapped):
        offset = record.rfind(base, 0, offset)
        late.append(offset)
    late.reverse()
    return early, late


def codon_grid(coord_map, cds_start, cds_end):
    """``{codon number: (column, column, column)}`` for one feature.

    ``coord_map`` is the master's ``{nucleotide position: alignment column}``
    map, so this is the same arithmetic the annotator already used to locate a
    codon - the grid simply covers the whole feature rather than one position.
    A codon the master's alignment cannot place is left out.
    """
    try:
        first, last = int(cds_start), int(cds_end)
    except (TypeError, ValueError):
        return {}
    grid = {}
    for number in range(1, (last - first + 2) // 3 + 1):
        start = first + (number - 1) * 3
        columns = tuple(coord_map.get(start + offset, -1) for offset in range(3))
        if min(columns) >= 0:
            grid[number] = columns
    return grid


class RecordMapping:
    """One alignment row placed back onto its own submitted record."""

    def __init__(self, alignment, record):
        self.alignment = normalise_bases(alignment)
        self.residue_columns = [i for i, base in enumerate(self.alignment) if base != '-']
        self.orientation = '+'
        self.record = normalise_bases(record)
        ungapped = self.alignment.replace('-', '')
        placed = placements(ungapped, self.record) if self.record and ungapped else None
        if placed is None and self.record and ungapped:
            flipped = reverse_complement(self.record)
            placed = placements(ungapped, flipped)
            if placed is not None:
                self.record, self.orientation = flipped, '-'
        self.early, self.late = placed if placed is not None else (None, None)

    @property
    def mapped(self):
        return self.early is not None

    def residue_index(self, column):
        index = bisect_left(self.residue_columns, column)
        if index < len(self.residue_columns) and self.residue_columns[index] == column:
            return index
        return None

    def positions(self, columns):
        """``(leftmost, rightmost)`` record positions of a codon's columns.

        ``None`` when any column is a gap in this row or past its end: that is a
        deletion or missing coverage, and the caller decides what it means.
        """
        if not self.mapped:
            return None
        indexes = []
        for column in columns:
            if column < 0 or column >= len(self.alignment) or self.alignment[column] == '-':
                return None
            indexes.append(self.residue_index(column))
        return [self.early[i] for i in indexes], [self.late[i] for i in indexes]

    def certain_start(self, columns):
        """Record start of a codon whose three bases sit side by side at one place."""
        placed = self.positions(columns)
        if placed is None or placed[0] != placed[1]:
            return None
        start = placed[0][0]
        return start if placed[0] == [start, start + 1, start + 2] else None


#: read() verdicts. 'alignment' means the record has nothing to say about this
#: codon and the caller's own alignment reading stands.
DEFER = 'alignment'
IN_FRAME = 'record_in_frame'
CORRECTED = 'record_corrected'
GUARD_UNDETERMINED = 'guard_frame_undetermined'
GUARD_AMBIGUOUS = 'guard_placement_ambiguous'
GUARD_UNANCHORED = 'guard_unanchored'
GUARD_FRAMESHIFT = 'guard_record_frameshift'
GUARDS = (GUARD_UNDETERMINED, GUARD_AMBIGUOUS, GUARD_UNANCHORED, GUARD_FRAMESHIFT)


class ProteinFrame:
    """One protein's reading frame in one record, and the indels around it."""

    def __init__(self, mapping, grid, translate, unknown='X', master_alignment=None,
                 flank_codons=FLANK_CODONS):
        self.mapping, self.grid = mapping, grid or {}
        self.translate, self.unknown = translate, unknown
        self.frame = None
        self.tie = False
        self.tie_broken_by_orf = False
        # The majority frame reads the record with internal stops while a frame
        # some codons were placed in reads it cleanly: no way to say which.
        self.frame_conflict = False
        self.split, self.off_frame, self.displaced = [], [], []
        self.ambiguous_codons = self.gapped_codons = 0
        self.internal_stops = 0
        self.record_frameshift = False
        self.insertions, self.deletions = [], []
        self._number_of = {tuple(columns): number for number, columns in self.grid.items()}
        self._answers = {}
        self._anchor_numbers, self._anchors, self._displaced_offsets = [], {}, {}
        self._displaced_all, self._unexplained = [], []
        self._ranges = {}
        self._master = str(master_alignment or '').upper()
        if not self.grid or mapping is None or not mapping.mapped:
            return

        certain = {}
        for number in sorted(self.grid):
            placed = mapping.positions(self.grid[number])
            if placed is None:
                self.gapped_codons += 1
                continue
            early, late = placed
            if early != late:
                self.ambiguous_codons += 1
            elif early != [early[0], early[0] + 1, early[0] + 2]:
                self.split.append(number)
            else:
                certain[number] = early[0]
        if not certain:
            return

        window = self._flanks(str(master_alignment or '').upper(), flank_codons)
        window.update((n, (start, self.grid[n][0])) for n, start in certain.items())
        displaced = displaced_numbers({n: start for n, (start, _) in window.items()})
        votes = Counter(start % 3 for n, start in certain.items() if n not in displaced)
        if not votes:
            # The whole protein is displaced; what it was displaced from decides.
            votes = Counter(start % 3 for n, (start, _) in window.items() if n not in displaced)
        if not votes:
            return

        phases = votes.most_common(2)
        self.frame = phases[0][0]
        self.tie = len(phases) == 2 and phases[0][1] == phases[1][1]
        if self.tie:
            # An even split leaves the majority meaningless, but a frame that
            # translates the record without internal stops is still the only
            # one that can be its coding frame.
            clean = [phase for phase, _ in phases if self._orf_stops(phase) == 0]
            if len(clean) == 1:
                self.frame, self.tie, self.tie_broken_by_orf = clean[0], False, True
        self.internal_stops = self._orf_stops(self.frame)
        if not self.tie:
            # A phase the majority disagrees with, held by codons nothing marks
            # as displaced, that reads the record cleanly: a shift running off
            # the end of the record reads its majority cleanly too, often
            # enough on a short protein to matter, so a clean majority does not
            # settle it either.
            others = Counter(start % 3 for n, start in certain.items()
                             if n not in displaced and start % 3 != self.frame)
            self.frame_conflict = any(
                self._orf_stops(phase) == 0 and (self.internal_stops or count >= MIN_RIVAL_CODONS)
                for phase, count in others.items())
        self.off_frame = sorted(n for n, position in certain.items() if position % 3 != self.frame)
        self.displaced = sorted(n for n in displaced if n in self.grid)
        # A record whose own ORF is broken here has no in-frame homologue for
        # the master's codon, so the off-frame ones are guarded rather than
        # re-read. An intact ORF means the shift is the aligner's, not the
        # submitter's, and re-reading recovers the right residue.
        self.record_frameshift = bool(self.off_frame or self.split) and self.internal_stops > 0
        self._anchors = {n: (start - 3 * n, column) for n, (start, column) in window.items()
                         if start % 3 == self.frame and n not in displaced}
        self._anchor_numbers = sorted(self._anchors)
        # Where displaced in-frame codons sit is one hypothesis for every codon
        # between the anchors around them, not only for themselves.
        self._displaced_offsets = {n: start - 3 * n for n, (start, _) in window.items()
                                   if n in displaced and start % 3 == self.frame}
        self._displaced_all = sorted(n for n in window if n in displaced)
        # Placed codons that are neither anchors nor displaced: a shift nothing
        # explains. Anchors either side of one do not vouch for what is between.
        self._unexplained = sorted(n for n in window if n not in displaced and n not in self._anchors)
        if master_alignment is not None:
            self._scan_indels(str(master_alignment or '').upper())

    # -- frame and anchors ---------------------------------------------------
    def _flanks(self, master_alignment, flank_codons):
        """Certain codon-width placements just outside the protein, both sides.

        ``{virtual codon number: (record start, first column)}``, numbered on
        from the protein's own grid (0, -1, ... upstream). They share the
        protein's register because they step through the master's residues in
        threes from its first and last codon.
        """
        placed = {}
        if not master_alignment or not flank_codons:
            return placed
        master_columns = [i for i, base in enumerate(master_alignment) if base != '-']
        numbers = sorted(self.grid)
        for edge, direction in ((numbers[0], -1), (numbers[-1], 1)):
            column = self.grid[edge][0]
            index = bisect_left(master_columns, column)
            if index >= len(master_columns) or master_columns[index] != column:
                continue
            for step in range(1, flank_codons + 1):
                first = index + direction * 3 * step
                if first < 0 or first + 3 > len(master_columns):
                    break
                columns = tuple(master_columns[first:first + 3])
                start = self.mapping.certain_start(columns)
                if start is not None:
                    placed[edge + direction * step] = (start, columns[0])
        return placed

    def _around(self, number):
        index = bisect_left(self._anchor_numbers, number)
        left = self._anchor_numbers[index - 1] if index > 0 else None
        if index < len(self._anchor_numbers) and self._anchor_numbers[index] == number:
            index += 1
        right = self._anchor_numbers[index] if index < len(self._anchor_numbers) else None
        return left, right

    def _mapped_span(self):
        """First and last certainly-placed record position across the feature."""
        if not hasattr(self, '_span'):
            positions = []
            for number in sorted(self.grid):
                placed = self.mapping.positions(self.grid[number])
                if placed is not None and placed[0] == placed[1]:
                    positions.extend(placed[0])
            self._span = (min(positions), max(positions)) if positions else None
        return self._span

    def _orf_stops(self, frame):
        """Internal stops when the record's own bases are read in ``frame``.

        The last codon of the span is excluded: a CDS ends on a stop, and a
        mature peptide's neighbour may begin with one.
        """
        span = self._mapped_span()
        if span is None:
            return 0
        start, end = span
        start += (frame - start) % 3
        residues = [self.translate(self.mapping.record[i:i + 3])
                    for i in range(start, end - 1, 3)]
        return sum(1 for residue in residues[:-1] if residue == '*')

    # -- indels ------------------------------------------------------------
    def _scan_indels(self, master_alignment):
        """Insertions and deletions inside the feature, in master codon terms.

        Three kinds, all expressed as "after master codon N": bases the aligner
        dropped from the row (stripped insertions, the ones that split codons),
        bases held in columns where the master is a gap (kept insertions, what a
        guide alignment preserves), and columns where the master has a base and
        this row does not (deletions). A length that is not a multiple of three
        is what moves the frame.
        """
        columns = sorted(column for triple in self.grid.values() for column in triple)
        if not columns:
            return
        first, last = columns[0], columns[-1]
        codon_of = {}
        for number, triple in self.grid.items():
            for column in triple:
                codon_of[column] = number

        def codon_at_or_before(column):
            # The nearest master codon column at or upstream of this one; a
            # multi-column master gap has no codon of its own.
            index = bisect_left(columns, column + 1) - 1
            return codon_of[columns[max(index, 0)]]

        alignment = self.mapping.alignment
        covered = self.mapping.residue_columns
        if not covered:
            return
        inside = range(first, min(last, len(alignment) - 1) + 1)
        for match in re.finditer(r'-+', alignment[first:last + 1]):
            start, end = first + match.start(), first + match.end() - 1
            if start <= covered[0] or end >= covered[-1]:
                continue  # padding, not a deletion
            length = sum(1 for column in range(start, end + 1)
                         if column < len(master_alignment) and master_alignment[column] != '-')
            if length:
                self.deletions.append((codon_at_or_before(start), length))
        for match in re.finditer(r'-+', master_alignment[first:last + 1] if master_alignment else ''):
            start, end = first + match.start(), first + match.end() - 1
            length = sum(1 for column in range(start, end + 1)
                         if column < len(alignment) and alignment[column] != '-')
            if length:
                self.insertions.append((codon_at_or_before(start), length, 'kept'))
        # Stripped insertions: record bases between two certainly-placed
        # residues that the row does not carry at all.
        anchors = [(index, self.mapping.early[index]) for index, column in enumerate(covered)
                   if column in inside and self.mapping.early[index] == self.mapping.late[index]]
        for (left_index, left_position), (right_index, right_position) in zip(anchors, anchors[1:]):
            extra = (right_position - left_position) - (right_index - left_index)
            if extra > 0:
                self.insertions.append((codon_at_or_before(covered[left_index]), extra, 'stripped'))

    # -- reading -----------------------------------------------------------
    def _codon_at(self, start):
        if start < 0 or start + 3 > len(self.mapping.record):
            return self.unknown
        return self.translate(self.mapping.record[start:start + 3])

    def _displaced_between(self, left, right):
        return [offset for n, offset in self._displaced_offsets.items() if left < n < right]

    @staticmethod
    def _any_between(numbers, left, right):
        index = bisect_left(numbers, left + 1)
        return index < len(numbers) and numbers[index] < right

    def _register_range(self, left, right):
        """Lowest and highest offset any codon between two anchors could have.

        Moving right from the left anchor, a codon's offset can only have risen
        by bases the record has there that the master does not (stripped from
        the row, or held in master-gap columns) and only fallen by master
        columns the row gapped. Moving left from the right anchor, the reverse.
        Whatever is a real indel and whatever the aligner's artefact, the truth
        is inside both bounds. For a one- or two-base artefact that leaves one
        in-frame offset; for a whole codon's worth of indels, several.
        """
        key = (left, right)
        if key not in self._ranges:
            (left_offset, left_column), (right_offset, right_column) = self._anchors[left], self._anchors[right]
            row, master = self.mapping.alignment, self._master
            gapped = kept = 0
            for column in range(left_column, min(right_column, len(row))):
                master_base = master[column] if column < len(master) else ''
                if master_base == '-':
                    kept += row[column] != '-'
                elif row[column] == '-':
                    gapped += 1
            record_bases = (3 * right + right_offset) - (3 * left + left_offset)
            row_bases = self.mapping.residue_index(right_column) - self.mapping.residue_index(left_column)
            gained = record_bases - row_bases + kept
            self._ranges[key] = (max(left_offset - gapped, right_offset - gained),
                                 min(left_offset + gained, right_offset + gapped), gapped)
        return self._ranges[key]

    def _between_anchors(self, number, own_offset=None):
        """Read a codon from its anchors, or guard. See the module docstring."""
        left, right = self._around(number)
        if left is None or right is None:
            return self.unknown, GUARD_UNANCHORED
        if self._any_between(self._unexplained, left, right):
            return self.unknown, GUARD_AMBIGUOUS
        (left_offset, _), (right_offset, _) = self._anchors[left], self._anchors[right]
        low, high, gapped = self._register_range(left, right)
        if gapped >= 3:
            # Enough master columns went missing here for a whole codon, and
            # nothing says which: this codon may be the one the record lost.
            return self.unknown, GUARD_AMBIGUOUS
        offsets = {offset for offset in range(low, high + 1) if offset % 3 == self.frame}
        offsets.update(self._displaced_between(left, right))
        offsets.update([left_offset, right_offset] + ([] if own_offset is None else [own_offset]))
        left_start, right_start = 3 * left + left_offset, 3 * right + right_offset
        starts = [3 * number + offset for offset in range(min(offsets), max(offsets) + 1, 3)]
        # A candidate overlapping an anchor codon cannot be this codon.
        starts = [start for start in starts if left_start + 3 <= start <= right_start - 3]
        residues = {self._codon_at(start) for start in starts}
        if len(residues) != 1:
            return self.unknown, GUARD_AMBIGUOUS
        if len(starts) > 1 and self.mapping.positions(self.grid[number]) is None:
            # Several registers fit, so this gapped codon may be the one the
            # record really lost: a residue, however unanimous, is not safe.
            return self.unknown, GUARD_AMBIGUOUS
        return residues.pop(), CORRECTED

    def _read_gapped(self, number):
        """A codon with a gap column: a deletion, missing coverage, or an artefact.

        Between two anchors, the record bases the row does not carry can be
        counted. None, with nothing displaced between, means every gap here is a
        real deletion in the record, and the alignment's deletion and coverage
        rules own the codon. Otherwise the row gapped bases the record still has
        nearby: a whole-codon gap there is not evidence of a deletion, and a
        part-gapped codon is only read when its anchors agree - if they do not,
        this codon may be the one really deleted.
        """
        columns = self.grid[number]
        covered = self.mapping.residue_columns
        if (not covered or columns[0] < covered[0] or columns[-1] > covered[-1]
                or self.frame is None):
            return None, DEFER
        left, right = self._around(number)
        if left is None or right is None:
            # Nothing on one side to count bases against: a gap here could be
            # half of an artefact whose other half is out of sight.
            return self.unknown, GUARD_UNANCHORED
        (left_offset, left_column), (right_offset, right_column) = self._anchors[left], self._anchors[right]
        record_bases = (3 * right + right_offset) - (3 * left + left_offset)
        row_bases = self.mapping.residue_index(right_column) - self.mapping.residue_index(left_column)
        # A displaced in-frame run between offers another register this codon
        # could be in; an unexplained shift means the anchors vouch for nothing.
        # A displaced off-frame run is just the artefact this gap belongs to.
        displaced = bool(self._displaced_between(left, right)) or self._any_between(self._unexplained, left, right)
        if record_bases == row_bases and not displaced:
            return None, DEFER
        if self.tie or self.frame_conflict:
            return self.unknown, GUARD_UNDETERMINED
        if self.internal_stops:
            return self.unknown, GUARD_FRAMESHIFT
        if (all(self.mapping.alignment[column] == '-' for column in columns)
                or left_offset != right_offset or displaced):
            return self.unknown, GUARD_AMBIGUOUS
        return self._between_anchors(number)

    def _resolve(self, number):
        placed = self.mapping.positions(self.grid[number])
        if placed is None:
            return self._read_gapped(number)
        if self.frame is None or self.tie or self.frame_conflict:
            return self.unknown, GUARD_UNDETERMINED
        if number in self._anchors:
            return self._codon_at(placed[0][0]), IN_FRAME
        if self.internal_stops:
            return self.unknown, GUARD_FRAMESHIFT
        own = None
        if number in self.displaced and placed[0] == placed[1] and placed[0][0] % 3 == self.frame:
            # In frame but out of register: where the alignment put it is one
            # of the candidates, not the answer.
            own = placed[0][0] - 3 * number
        return self._between_anchors(number, own)

    def read(self, columns):
        """``(residue, verdict)`` for one master codon of this protein.

        ``DEFER`` leaves the codon to the caller's own alignment reading, which
        owns the deletion and coverage rules. Every other verdict carries the
        residue to use. Columns that are not one of this protein's codons are
        deferred: there is no frame to judge them against.
        """
        if self.mapping is None or not self.mapping.mapped:
            return None, DEFER
        number = self._number_of.get(tuple(columns))
        if number is None:
            return None, DEFER
        if number not in self._answers:
            self._answers[number] = self._resolve(number)
        return self._answers[number]

    # -- reporting ---------------------------------------------------------
    @property
    def frameshift_indels(self):
        return ([event for event in self.insertions if event[1] % 3]
                + [event for event in self.deletions if event[1] % 3])

    @property
    def has_issue(self):
        return bool(self.split or self.off_frame or self.displaced or self.frameshift_indels
                    or self.record_frameshift or self.tie or self.frame_conflict
                    or (self.grid and self.frame is None))

    def status(self):
        if self.mapping is None or not self.mapping.mapped:
            return 'record_unmapped'
        if self.frame is None:
            return 'frame_unreadable'
        if self.tie or self.frame_conflict:
            return 'frame_undetermined'
        if self.record_frameshift:
            return 'record_frameshift'
        if self.split or self.off_frame or self.displaced:
            return 'alignment_shift_corrected'
        return 'in_frame'

    def scan_row(self, limit=25):
        """Flat summary of this protein in this record, for the scan table."""
        def events(items):
            # 'c31:1nt:stripped' = one base after master codon 31 the row dropped.
            return ';'.join(':'.join([f'c{event[0]}', f'{event[1]}nt', *event[2:]])
                            for event in items)[:2000]

        unresolved = sum(1 for columns in self.grid.values() if self.read(columns)[1] in GUARDS)
        return {
            'frame_status': self.status(),
            'record_orientation': self.mapping.orientation if self.mapping else '',
            'record_frame_phase': '' if self.frame is None else self.frame,
            'frame_tie_broken_by_orf': int(self.tie_broken_by_orf),
            'off_frame_codons': len(self.off_frame),
            'off_frame_codon_runs': ';'.join(f'{a}-{b}' for a, b in runs(self.off_frame))[:2000],
            'split_codons': len(self.split),
            'split_codon_list': ';'.join(str(number) for number in self.split[:limit]),
            'displaced_codons': len(self.displaced),
            'unresolved_codons': unresolved,
            'gapped_codons': self.gapped_codons,
            'ambiguous_codons': self.ambiguous_codons,
            'insertions': events(self.insertions[:limit]),
            'deletions': events(self.deletions[:limit]),
            'frameshift_indels': len(self.frameshift_indels),
            'record_internal_stops': self.internal_stops,
        }


def displaced_numbers(starts):
    """Codon numbers in displaced runs, from ``{codon number: certain record start}``.

    Runs are consecutive codons (in number order) sharing an offset. A run is
    displaced when its two neighbouring runs share an offset it does not have,
    or share a phase it does not have - the second catches a shift whose window
    also holds a real whole-codon indel, which leaves the neighbours a codon
    apart. The shortest displaced runs go first, their neighbours merge, and
    the rest are judged again.
    """
    runs_ = []
    for number in sorted(starts):
        offset = starts[number] - 3 * number
        if runs_ and runs_[-1][0] == offset:
            runs_[-1][1].append(number)
        else:
            runs_.append([offset, [number]])
    displaced = set()
    while True:
        # An artefact covers less sequence than the in-register sequence it
        # was displaced from; two short, themselves-shifted runs agreeing
        # around a long one are not evidence against it.
        candidates = [i for i in range(1, len(runs_) - 1)
                      if (runs_[i - 1][0] == runs_[i + 1][0] != runs_[i][0]
                          or runs_[i - 1][0] % 3 == runs_[i + 1][0] % 3 != runs_[i][0] % 3)
                      and len(runs_[i][1]) < len(runs_[i - 1][1]) + len(runs_[i + 1][1])]
        if not candidates:
            return displaced

        def weight(i):
            return len(runs_[i][1]) / (len(runs_[i - 1][1]) + len(runs_[i + 1][1]))

        # Least sequence against the most either side goes first. Two runs next
        # to each other are never taken in the same round: when two equal
        # artefacts sit around an equal correct stretch, taking one artefact
        # first lets the correct stretch rejoin the sequence around it.
        lightest = min(weight(i) for i in candidates)
        flagged = set()
        for i in sorted(candidates):
            if weight(i) == lightest and i - 1 not in flagged:
                flagged.add(i)
        merged = []
        for i, (offset, members) in enumerate(runs_):
            if i in flagged:
                displaced.update(members)
            elif merged and merged[-1][0] == offset:
                merged[-1][1].extend(members)
            else:
                merged.append([offset, list(members)])
        runs_ = merged


def runs(numbers):
    """[3, 4, 5, 9] -> [(3, 5), (9, 9)]"""
    grouped = []
    for number in sorted(numbers):
        if grouped and number == grouped[-1][1] + 1:
            grouped[-1][1] = number
        else:
            grouped.append([number, number])
    return [tuple(pair) for pair in grouped]
