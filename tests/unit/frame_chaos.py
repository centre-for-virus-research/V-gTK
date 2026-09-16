"""Ground truth for adversarial reading-frame tests. Imports nothing from scripts/.

A case is built the way the pipeline meets one in the wild:

1. a master ORF with UTRs, optionally full of the sequence that makes placement
   hard: homopolymers, dinucleotide repeats, repeated codons;
2. a submitted record descended from it by *known* events: substitutions,
   ambiguity codes, nonsense stops, whole-codon insertions and deletions,
   optionally a real uncompensated frameshift, optionally cut short at either
   end;
3. the true alignment row (each record base in its homologous master column,
   inserted bases stripped as nextalign strips them), then damaged the way a
   frame-unaware aligner damages it: gaps slid sideways, a gap in one place
   paired with dropped bases elsewhere in either order, overlapping, long
   range, and across the CDS boundary;
4. the stored record in whatever spelling a submission arrives in.

Truth for a master codon comes from the event history, never from an alignment:
the residue of its three homologous record bases when they sit side by side, a
deletion when all three were deleted, and nothing otherwise (any answer but
unknown is then wrong). An indel whose position is not unique - the same record
results from placing it elsewhere - contributes every placement's truth.
"""

from itertools import product
import random

AMINO = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
CODE = {''.join(c): AMINO[i] for i, c in enumerate(product('TCAG', repeat=3))}
SENSE = sorted(codon for codon, residue in CODE.items() if residue != '*')
IUPAC = 'RYKMSWN'
COMPLEMENT = str.maketrans('ACGTRYKMSWN', 'TGCAYRMKSWN')


def translate(codon):
    return CODE.get(codon.upper(), 'X') if len(codon) == 3 else 'X'


def column_reading(row, columns):
    """What the annotator reads from the alignment columns alone (its DEFER path)."""
    row = row.upper()
    codon = ''.join(row[c] for c in columns)
    if codon == '---':
        covered = [i for i, base in enumerate(row) if base != '-']
        inside = covered and covered[0] < columns[0] and columns[-1] < covered[-1]
        return '-' if inside else 'X'
    return 'X' if '-' in codon else translate(codon)


class Case:
    def __init__(self, seed, mode):
        self.seed, self.mode = seed, mode
        self.log = []

    def describe(self):
        return f'seed={self.seed} mode={self.mode}\n  ' + '\n  '.join(self.log) + (
            f'\n  master={self.master}\n  row   ={self.row}\n  record={self.record}')

    @property
    def grid(self):
        return {n: (self.cds_start + 3 * (n - 1), self.cds_start + 3 * (n - 1) + 1,
                    self.cds_start + 3 * (n - 1) + 2) for n in range(1, self.codons + 1)}


def _body(rng, count, repeats):
    codons = []
    while len(codons) < count:
        roll = rng.random()
        if roll >= repeats:
            codons.append(rng.choice(SENSE))
        elif roll < repeats * 0.35:
            codons.extend([rng.choice('ACGT') * 3] * rng.randint(2, 5))
        elif roll < repeats * 0.7:
            run = rng.choice(['AC', 'AG', 'CA', 'GT', 'TG', 'CT', 'GA', 'TC']) * rng.randint(3, 9)
            run = run[:len(run) // 3 * 3]
            codons.extend(run[i:i + 3] for i in range(0, len(run), 3))
        else:
            codons.extend([rng.choice(SENSE)] * rng.randint(2, 6))
    return [codon if CODE[codon] != '*' else rng.choice(SENSE) for codon in codons[:count]]


def _equivalent(labelled, event):
    """Every placement of an indel that leaves the same record string."""
    kind, position, payload = event
    bases = ''.join(base for base, _ in labelled)
    if kind == 'del':
        length = payload
        result = bases[:position] + bases[position + length:]
        window = range(max(0, position - 3 * length - 9), min(len(bases) - length, position + 3 * length + 9) + 1)
        return [p for p in window if bases[:p] + bases[p + length:] == result]
    result = bases[:position] + payload + bases[position:]
    length = len(payload)
    window = range(max(0, position - 3 * length - 9), min(len(bases), position + 3 * length + 9) + 1)
    return [p for p in window if bases[:p] + result[p:p + length] + bases[p:] == result]


def _apply(labelled, event, position=None):
    kind, canonical, payload = event
    position = canonical if position is None else position
    if kind == 'del':
        return labelled[:position] + labelled[position + payload:]
    bases = ''.join(base for base, _ in labelled)
    inserted = (bases[:canonical] + payload + bases[canonical:])[position:position + len(payload)]
    return labelled[:position] + [(base, None) for base in inserted] + labelled[position:]


def make_master(rng, mode):
    """``(master, cds_start (0-based), codons)`` with UTRs either side of an ORF."""
    hostile = mode != 'benign'
    repeats = {'benign': 0.0, 'repeat_hell': 0.9}.get(mode, 0.35)
    codons = rng.randint(14, 70) if mode != 'majority_shift' else rng.randint(12, 30)
    utr5 = ''.join(rng.choice('ACGT') for _ in range(rng.randint(0, 45) if hostile else 30))
    utr3 = ''.join(rng.choice('ACGT') for _ in range(rng.randint(0, 45) if hostile else 30))
    if mode == 'majority_shift':
        utr5 = utr5 or 'ACGTTGCA' * 4
        utr3 = utr3 or 'TTGACCAG' * 4
    orf = 'ATG' + ''.join(_body(rng, codons - 2, repeats)) + 'TAA'
    return utr5 + orf + utr3, len(utr5), codons


def build(seed, mode, master=None):
    """One adversarial case. ``mode`` picks how hostile it is; see MODES.

    ``master`` from make_master() shares one master between cases, as rows of
    one alignment do; by default each case gets its own.
    """
    rng = random.Random(f'{mode}:{seed}')
    case = Case(seed, mode)
    spaced, mode = mode.startswith('decidable_'), behaviour(mode)
    hostile = mode != 'benign'
    case.master, case.cds_start, case.codons = master or make_master(rng, mode)
    real_columns = []
    codons = case.codons
    cds_end = case.cds_start + 3 * codons

    # --- 2. the record's real history --------------------------------------
    bases = list(case.master)
    if hostile:
        for i in range(len(bases)):
            if rng.random() < 0.04:
                bases[i] = rng.choice([b for b in 'ACGT' if b != bases[i]])
            elif rng.random() < 0.006:
                bases[i] = rng.choice(IUPAC)
        if rng.random() < 0.15:
            k = rng.randint(2, codons - 1)
            at = case.cds_start + 3 * (k - 1)
            bases[at:at + 3] = 'TAA'
            case.log.append(f'nonsense stop at codon {k}')
    labelled = [(base, index) for index, base in enumerate(bases)]

    events = []
    if hostile:
        for _ in range(rng.choice([0, 0, 1, 1, 2, 3])):
            k = rng.randint(2, codons - 1)
            index = next((i for i, (_, label) in enumerate(labelled)
                          if label == case.cds_start + 3 * (k - 1)), None)
            if index is None:
                continue
            if rng.random() < 0.5:
                event = ('del', index, 3 * rng.randint(1, 2))
            else:
                event = ('ins', index, ''.join(rng.choice(SENSE) for _ in range(rng.randint(1, 2))))
            events.append(event)
            real_columns.append(case.cds_start + 3 * (k - 1))
            labelled = _apply(labelled, event)
            case.log.append(f'real {event[0]} before codon {k}: {event[2]}')
    if mode == 'frameshift':
        index = rng.randint(case.cds_start + 3, cds_end - 4)
        index = min(index, len(labelled) - 3)
        if rng.random() < 0.5:
            event = ('del', index, rng.randint(1, 2))
        else:
            event = ('ins', index, ''.join(rng.choice('ACGT') for _ in range(rng.randint(1, 2))))
        events.append(event)
        real_columns.append(next((label for _, label in labelled[index:] if label is not None), len(case.master) - 1))
        labelled = _apply(labelled, event)
        case.log.append(f'real frameshift {event[0]} at record {index}: {event[2]}')

    # Every equivalent placement of each event, holding the others fixed.
    histories = [None]
    replay = [(base, index) for index, base in enumerate(bases)]
    for number, event in enumerate(events):
        for alternative in _equivalent(replay, event):
            if alternative != event[1]:
                histories.append((number, alternative))
        replay = _apply(replay, event)
    variants = []
    for history in histories:
        current = [(base, index) for index, base in enumerate(bases)]
        for number, event in enumerate(events):
            position = history[1] if history and history[0] == number else None
            current = _apply(current, event, position)
        variants.append(current)

    trim5 = trim3 = 0
    if hostile and rng.random() < 0.3:
        trim5 = rng.randint(0, case.cds_start + 3 * codons // 2)
    if hostile and rng.random() < 0.3:
        trim3 = rng.randint(0, len(labelled) - case.cds_start - 3 * codons // 2)
    if trim5 or trim3:
        case.log.append(f'record cut: {trim5} bases off the 5\' end, {trim3} off the 3\' end')
    variants = [v[trim5:len(v) - trim3] for v in variants]
    labelled = variants[0]

    truth = {n: set() for n in range(1, codons + 1)}
    for variant in variants:
        where = {label: i for i, (_, label) in enumerate(variant) if label is not None}
        record = ''.join(base for base, _ in variant)
        for n in truth:
            members = [case.cds_start + 3 * (n - 1) + u for u in range(3)]
            present = [where.get(m) for m in members]
            if None not in present and present == list(range(present[0], present[0] + 3)):
                truth[n].add(translate(record[present[0]:present[0] + 3]))
            elif present == [None, None, None] and _deleted(members, variant, where, bases, trim5, trim3):
                truth[n].add('-')
    case.truth = truth

    # --- 3. the true alignment row, then the aligner's damage ---------------
    row = ['-'] * len(case.master)
    for base, label in labelled:
        if label is not None:
            row[label] = base
    covered = [i for i, base in enumerate(row) if base != '-']
    if covered:
        _damage(rng, case, row, covered, mode, real_columns if spaced else None)
    case.row = ''.join(row)
    record = ''.join(base for base, _ in labelled)

    # --- 4. how the record is stored ------------------------------------------
    if hostile:
        if rng.random() < 0.2:
            record = record.translate(COMPLEMENT)[::-1]
            case.log.append('record stored reverse-complemented')
        if rng.random() < 0.15:
            record = record.lower()
        if rng.random() < 0.1:
            record = record.replace('T', 'U')
            case.log.append('record stored as RNA')
        if rng.random() < 0.1 and len(record) > 10:
            cut = rng.randrange(1, len(record) - 1)
            record = record[:cut] + '\n' + record[cut:]
        if rng.random() < 0.1:
            case.row = case.row.lower()
    case.record = record
    return case


def _deleted(members, variant, where, bases, trim5, trim3):
    """True when a codon's bases were removed by an event, not cut off an end."""
    labels = [label for _, label in variant if label is not None]
    if not labels:
        return False
    return min(labels) < members[0] and members[-1] < max(labels)


def _damage(rng, case, row, covered, mode, occupied=None):
    first, last = covered[0], covered[-1]
    cds_first, cds_last = case.cds_start, case.cds_start + 3 * case.codons - 1
    spaced = occupied is not None
    real = list(occupied or [])
    occupied = []

    def room(x, y):
        # Clean sequence either side must outweigh the artefact and keep a few
        # codons, nothing else may sit inside that margin, and a real indel must
        # be further off still: an artefact a few codons from a real frameshift
        # or codon indel can cancel it, and then two equally short stories fit
        # the same row and record.
        width, left, right = y - x, x - first, last - y
        if not spaced:
            return True
        near = max(width, 12)
        far = max(3 * width, 30)
        return not (left < 12 or right < 12 or left + right <= width + 6
                    or any(x - near <= c < y + near for c in occupied)
                    or any(x - far <= c < y + far for c in real))

    def shuffle(x, y, length, gap_first):
        if not room(x, y):
            return
        occupied.extend(range(x, y))
        segment = row[x:y]
        if gap_first:
            kept = [b for b in segment[:len(segment) - length]]
            row[x:y] = ['-'] * length + kept
        else:
            row[x:y] = segment[length:] + ['-'] * length
        case.log.append(f"artefact: {'gap then drop' if gap_first else 'drop then gap'} "
                        f"{length}nt over columns {x}-{y - 1}")

    if mode == 'benign':
        # Isolated one- and two-base artefacts well inside the CDS, apart from
        # each other: every codon has an agreeing anchor on each side.
        start = cds_first + 12
        while start + 20 < cds_last - 12 and rng.random() < 0.8:
            length = rng.randint(1, 2)
            width = rng.randint(length + 3, 14)
            shuffle(start, start + width, length, rng.random() < 0.5)
            start += width + rng.randint(9, 20)
        return

    if mode == 'majority_shift':
        length = rng.randint(1, 2)
        span = int((cds_last - cds_first) * rng.uniform(0.55, 0.85))
        if rng.random() < 0.5:
            x = max(first, cds_first - rng.randint(3, 15))
            y = min(last + 1, cds_first + span)
        else:
            y = min(last + 1, cds_last + rng.randint(3, 15))
            x = max(first, cds_last - span)
        if y - x > length + 3:
            shuffle(x, y, length, rng.random() < 0.5)
        return

    count = {'repeat_hell': rng.randint(2, 6)}.get(mode, rng.randint(0, 4))
    for _ in range(count):
        roll = rng.random()
        length = rng.choice([1, 1, 2, 2, 3, 4, 5])
        if roll < 0.2:
            if spaced:
                continue  # a slide sits on a real indel by definition
            # Slide an existing gap run sideways: no base lost, just moved.
            gaps = [i for i in range(first + 1, last) if row[i] == '-']
            if not gaps:
                continue
            g = rng.choice(gaps)
            x, y = max(first, g - rng.randint(1, 6)), min(last + 1, g + rng.randint(1, 6))
            segment = row[x:y]
            moved = [b for b in segment if b != '-']
            holes = len(segment) - len(moved)
            row[x:y] = (['-'] * holes + moved) if rng.random() < 0.5 else (moved + ['-'] * holes)
            case.log.append(f'artefact: gaps slid within columns {x}-{y - 1}')
            continue
        if roll < 0.35:
            width = rng.randint(length + 3, max(length + 4, (last - first) // 2))
        elif roll < 0.5:
            # Across the CDS boundary.
            edge = rng.choice([cds_first, cds_last])
            x = max(first, edge - rng.randint(1, 20))
            width = rng.randint(length + 3, 40)
            y = min(last + 1, x + width)
            if y - x > length + 1:
                shuffle(x, y, length, rng.random() < 0.5)
            continue
        else:
            width = rng.randint(length + 2, 18)
        x = rng.randint(first, max(first, last - width))
        y = min(last + 1, x + width)
        if y - x > length + 1:
            shuffle(x, y, length, rng.random() < 0.5)


HOSTILE = ('chaos', 'repeat_hell', 'frameshift', 'majority_shift')
#: Same hostility, but every artefact keeps more clean sequence either side than
#: it spans, away from real indels, record ends and other artefacts: the
#: conditions under which the row and record alone can say what happened.
DECIDABLE = tuple(f'decidable_{mode}' for mode in HOSTILE)
MODES = ('benign',) + DECIDABLE + HOSTILE


def behaviour(mode):
    return mode[len('decidable_'):] if mode.startswith('decidable_') else mode
