"""Adversarial datasets through padding -> coordinates -> DB -> calls -> audit.

Expected calls are hand-authored in the fixture or come from frame_chaos's
event histories, never calculated with the annotator/verifier under test.
Corruption tests prove that the independent auditor can reject
plausible-looking but wrong outputs. All writes use tmp_path. No Nextflow,
network, real viral genomes or external aligner is needed.
"""

import csv
import hashlib
import json
from pathlib import Path
import random
import sqlite3
import sys

from Bio import SeqIO
from Bio.Seq import Seq
import pandas as pd
import pytest

import AnnotateMutations
import AuditAlignmentCoordinates as Audit
from CalcAlignmentCord import CalculateAlignmentCoordinates
import frame_chaos as chaos
from PadAlignment import PadAlignment


FIXTURE = Path(__file__).resolve().parents[2] / 'test_data/unit/alignment_red_team/cases.json'
CATALOG_FIELDS = ['mutation_id', 'protein_name', 'segment', 'aa_position', 'alt_residue',
                  'reference_accession', 'mutation_type', 'signature_id', 'signature_kind',
                  'combination_id', 'combination_size', 'phenotype', 'resistance_category', 'drug']


def write_fasta(path, rows):
    path.write_text(''.join(f'>{name}\n{seq}\n' for name, seq in rows.items()))


def run_pipeline(tmp_path, monkeypatch, master, cds, guide_rows, query_rows, raw_records, catalog):
    """Pad, project coordinates, load a DB and annotate it, as the workflow does.

    ``cds`` is the 1-based (start, end) of POL in ``master``; ``catalog`` is a
    list of (aa_position, alt_residue) pairs.
    """
    guide = tmp_path / 'MASTER.fasta'
    write_fasta(guide, guide_rows)
    inputs = tmp_path / 'inputs' / 'REF_INS'
    inputs.mkdir(parents=True)
    write_fasta(inputs / 'REF_INS.aligned.fasta', query_rows)
    padded = tmp_path / 'padded'
    padder = PadAlignment(str(guide), str(inputs.parent), str(tmp_path), str(padded), True)
    padder.process_master_alignment(str(guide), str(inputs.parent), str(tmp_path), str(padded), True)
    # Use the production dedup step: PadAlignment emits a guide reference and
    # its nextalign reference row; this is covered separately in its own suite.
    padder.remove_redundant_sequences()
    merged = padded / 'MASTER_merged_MSA.fasta'
    alignments = {r.id: str(r.seq) for r in SeqIO.parse(merged, 'fasta')}
    assert alignments == {**guide_rows, **query_rows}
    coordinate_input = tmp_path / 'coordinate_input'
    coordinate_input.mkdir()
    write_fasta(coordinate_input / 'MASTER.aligned_merged_MSA.fasta', alignments)
    master_seq = tmp_path / 'master_seq'
    master_seq.mkdir()
    write_fasta(master_seq / 'MASTER.fasta', {'MASTER': master})
    gff = tmp_path / 'MASTER.gff3'
    gff.write_text(
        f'##gff-version 3\nMASTER\ttest\tregion\t1\t{len(master)}\t.\t+\t.\tID=MASTER\n'
        f'MASTER\ttest\tCDS\t{cds[0]}\t{cds[1]}\t.\t+\t0\tID=pol;product=POL\n'
    )
    masters = tmp_path / 'masters.tsv'
    masters.write_text('MASTER\n')
    hits = tmp_path / 'hits.tsv'
    hits.write_text(''.join(f'{acc}\tREF_INS\t99.0\tplus\n' for acc in query_rows))
    CalculateAlignmentCoordinates(
        str(coordinate_input), [str(gff)], str(tmp_path), 'Tables', 'features.tsv',
        str(masters), str(hits), master_seq_dir=str(master_seq),
    ).find_gaps_in_fasta()
    features = pd.read_csv(tmp_path / 'Tables/features.tsv', sep='\t', dtype=str, keep_default_na=False)
    features['segment'] = '1'
    db = tmp_path / 'red_team.db'
    with sqlite3.connect(db) as conn:
        features.to_sql('features', conn, index=False)
        pd.DataFrame([
            {'primary_accession': acc, 'sequence_id': acc, 'alignment_name':
             'MASTER' if acc == 'MASTER' else 'REF_INS', 'alignment': sequence, 'segment': '1'}
            for acc, sequence in alignments.items()
        ]).to_sql('sequence_alignment', conn, index=False)
        pd.DataFrame([
            {'primary_accession': acc, 'accession_type': 'master' if acc == 'MASTER' else 'query'}
            for acc in alignments
        ]).to_sql('meta_data', conn, index=False)
        pd.DataFrame([
            {'header': acc, 'sequence': raw_records.get(acc, seq.replace('-', '')),
             'segment': '1'} for acc, seq in alignments.items()
        ]).to_sql('sequences', conn, index=False)

    catalog_path = tmp_path / 'catalog.tsv'
    with catalog_path.open('w') as handle:
        writer = csv.DictWriter(handle, fieldnames=CATALOG_FIELDS, delimiter='\t')
        writer.writeheader()
        for pos, alt in catalog:
            mutation = f"POL:{pos}{'del' if alt == '-' else alt}"
            writer.writerow(dict(mutation_id=mutation, protein_name='POL', segment='1',
                                 aa_position=pos, alt_residue=alt, reference_accession='MASTER',
                                 mutation_type='snp', signature_id=mutation, signature_kind='single'))
    monkeypatch.setattr(sys, 'argv', ['AnnotateMutations.py', '--db', str(db),
                                     '--mutation_catalog', str(catalog_path), '--virus', 'generic'])
    AnnotateMutations.main()
    return db


def build_dataset(tmp_path, monkeypatch, gaps=0, trim=False, shuffle=False):
    data = json.loads(FIXTURE.read_text())
    utr, orf = data['utr'], data['master']
    prefix = '' if trim else utr

    def project(sequence, insertion=''):
        # Deliberately split codon 2 across insertion columns. Column slicing
        # and residue slicing can no longer accidentally give the same codon.
        return sequence[:4] + insertion + sequence[4:]

    guide_rows = {
        'MASTER': prefix + project(orf, '-' * gaps),
        'REF_INS': prefix + project(orf, 'C' * gaps),
    }
    query_rows = {'REF_INS': guide_rows['REF_INS']}
    # Submitted records. 'raw' carries bases nextalign stripped as insertions;
    # they sit after codon 2, so project() places the guide insertion the same way.
    raw_records = {'MASTER': utr + orf}
    for row in data['rows']:
        partial = row['sequence'].startswith('-')
        row_prefix = '-' * len(prefix) if partial else prefix
        inserted = '-' * gaps if partial else 'C' * gaps
        query_rows[row['id']] = row_prefix + project(row['sequence'], inserted)
        raw_records[row['id']] = (row_prefix + project(row.get('raw', row['sequence']), inserted)).replace('-', '')
    if shuffle:
        items = list(query_rows.items())
        random.Random(17).shuffle(items)
        query_rows = dict(items)
    catalog = [(2, 'V'), (3, '-'), (5, '*'), (5, 'V'), (9, 'V')]
    db = run_pipeline(tmp_path, monkeypatch, utr + orf, (len(utr) + 1, len(utr + orf)),
                      guide_rows, query_rows, raw_records, catalog)
    return db, data


def expected_calls(data, gaps):
    key = 'calls_codon2_split' if gaps else 'calls'
    expected = {row['id']: set(row.get(key, row['calls'])) for row in data['rows']}
    return {acc: calls for acc, calls in expected.items() if calls}


@pytest.mark.parametrize('gaps,trim,shuffle', [
    (0, False, False), (1, False, False), (2, False, False),
    (3, False, False), (9, False, False), (3, True, False), (3, True, True), (2, True, True),
])
def test_fixed_red_team_dataset_across_coordinate_settings(tmp_path, monkeypatch, gaps, trim, shuffle):
    db, data = build_dataset(tmp_path, monkeypatch, gaps, trim, shuffle)
    with sqlite3.connect(db) as conn:
        actual = {acc: set(mutations.rstrip(';').split(';'))
                  for acc, mutations in conn.execute(
                      'SELECT primary_accession,relevant_mutations_present FROM sequence_relevant_mutation_summary')}
    expected = expected_calls(data, gaps)
    report = Audit.audit_database(db)
    assert report['ok'], report
    counts = report['counts']
    if gaps % 3 == 0:
        assert actual == expected  # Both false positives and false negatives.
        assert counts['mutation_call_rows'] == sum(len(calls) for calls in expected.values())
        assert counts['record_corrected_residues_verified'] == 1  # STRIPPED_INS_FRAMESHIFT's V
    else:
        # A one- or two-base guide insertion is a real frameshift inside codon 2
        # of every full-length row. Codons past it are still their own bases, so
        # the right answers do not change, but whether each can be read now
        # depends on how much frame evidence survives. Missing is allowed; wrong
        # is not.
        assert all(acc in expected and calls <= expected[acc] for acc, calls in actual.items()), (actual, expected)
    assert counts.get('mutation_residues_verified', 0) + counts.get('record_corrected_residues_verified', 0) \
        == counts['mutation_call_rows']
    assert counts['alignment_rows'] == 16
    assert not report['skipped']
    assert 'emitted_call_not_one_raw_codon' not in report['warnings'], report['warnings']
    assert 'record_corrected_residue_ambiguous' not in report['warnings']
    # A stop and a one-base deletion were planted deliberately. Geometry can
    # pass while coding diagnostics fail; these must not get conflated.
    assert report['warnings']['internal_stop_in_degapped_feature']['count'] >= 1
    assert report['warnings']['coding_length_not_multiple_of_three']['count'] >= 1


def test_calls_are_not_emitted_from_codons_absent_from_the_raw_record(tmp_path, monkeypatch):
    db, _ = build_dataset(tmp_path, monkeypatch)
    with sqlite3.connect(db) as conn:
        called = dict(conn.execute(
            'SELECT primary_accession, relevant_mutations_present FROM sequence_relevant_mutation_summary'))
        read = dict(conn.execute(
            "SELECT primary_accession, codon_read FROM sequence_mutation_calls WHERE aa_position = 5"))
    assert 'STRIPPED_INS_SPLIT5' not in called  # its TAA columns span a stripped G
    assert called['STRIPPED_INS_FRAMESHIFT'].rstrip(';') == 'POL:5V'  # not the TAA its columns hold
    assert read['STRIPPED_INS_FRAMESHIFT'] == 'record_corrected'
    assert read['STOP5'] == 'record'


@pytest.fixture
def red_team_db(tmp_path, monkeypatch):
    return build_dataset(tmp_path, monkeypatch, gaps=3, trim=True)[0]


@pytest.mark.parametrize('statement,error', [
    ("UPDATE features SET cds_start_OG_seq=cds_start_OG_seq+1 WHERE accession='ALT2'", 'feature_og_mismatch'),
    ("UPDATE features SET cds_end_OG_seq=cds_end_OG_seq-3 WHERE accession='ALT2'", 'feature_slice_mismatch'),
    ("UPDATE features SET aln_start='2' WHERE accession='ALT2'", 'covered_span_mismatch'),
    ("UPDATE features SET cds_start='1.5' WHERE accession='ALT2'", 'invalid_feature_coordinate'),
    ("UPDATE features SET cds_start_OG_seq='0' WHERE accession='ALT2'", 'invalid_feature_coordinate'),
    ("UPDATE features SET master_ref_accession='MISSING' WHERE accession='ALT2'", 'missing_master_feature'),
    ("UPDATE features SET segment='2' WHERE accession='MASTER'", 'missing_master_feature'),
    ("UPDATE sequence_alignment SET alignment=substr(alignment,1,length(alignment)-1) WHERE primary_accession='ALT2'", 'ragged_alignment'),
    ("UPDATE sequence_alignment SET alignment=replace(alignment,'A','?') WHERE primary_accession='ALT2'", 'invalid_alignment_alphabet'),
    ("INSERT INTO sequence_alignment SELECT * FROM sequence_alignment WHERE primary_accession='ALT2'", 'duplicate_alignment'),
    ("UPDATE sequence_mutation_calls SET observed_residue='K' WHERE primary_accession='STOP5'", 'mutation_residue_mismatch'),
    ("UPDATE sequence_mutation_calls SET aa_position=9 WHERE primary_accession='STOP5'", 'mutation_outside_feature'),
    ("UPDATE sequence_mutation_calls SET primary_accession='PARTIAL5' WHERE primary_accession='DELETION3'", 'mutation_residue_mismatch'),
    ("UPDATE sequences SET sequence='AAA' WHERE header='ALT2'", 'alignment_not_raw_subsequence'),
    # A record-corrected call is checked against the raw record, not waved through.
    ("UPDATE sequence_mutation_calls SET observed_residue='L', alt_residue='L' WHERE codon_read='record_corrected'",
     'record_corrected_residue_unverified'),
    ("UPDATE sequence_mutation_calls SET aa_position=4 WHERE codon_read='record_corrected'",
     'record_corrected_residue_unverified'),
    ("UPDATE sequence_mutation_calls SET codon_read='record_corrected' WHERE primary_accession='DELETION3'",
     'record_corrected_residue_unverified'),
    # STOP5's columns already read its own in-frame TAA: the residue is right,
    # but a claim that the record corrected it is false provenance.
    ("UPDATE sequence_mutation_calls SET codon_read='record_corrected' WHERE primary_accession='STOP5'",
     'record_corrected_residue_unverified'),
    ("UPDATE sequence_mutation_calls SET codon_read='alignment' WHERE codon_read='record_corrected'",
     'mutation_residue_mismatch'),
])
def test_auditor_rejects_deliberate_corruption(red_team_db, statement, error):
    with sqlite3.connect(red_team_db) as conn:
        assert conn.execute(statement).rowcount, 'the corruption must touch a row'
    report = Audit.audit_database(red_team_db)
    if error is None:
        assert report['ok'], report
        return
    assert not report['ok']
    assert error in report['errors'], report


def test_auditor_is_read_only_and_coding_warnings_can_fail_cli(red_team_db, tmp_path, capsys):
    before = hashlib.sha256(red_team_db.read_bytes()).hexdigest()
    output = tmp_path / 'audit.json'
    assert Audit.main(['--db', str(red_team_db), '--json', str(output)]) == 0
    assert Audit.main(['--db', str(red_team_db), '--strict-coding']) == 1
    assert json.loads(output.read_text())['ok']
    assert hashlib.sha256(red_team_db.read_bytes()).hexdigest() == before


def test_auditor_rejects_empty_schema(tmp_path):
    db = tmp_path / 'empty.db'
    sqlite3.connect(db).close()
    assert 'missing_schema' in Audit.audit_database(db)['errors']


def test_auditor_never_creates_a_missing_database(tmp_path):
    missing = tmp_path / 'missing.db'
    with pytest.raises(sqlite3.OperationalError):
        Audit.audit_database(missing)
    assert not missing.exists()


def test_reference_codons_can_pass_geometry_and_cross_raw_reading_frame(tmp_path):
    """An intact raw ORF acquires a stop in reference columns after projection.

    The aligned TAG uses raw bases 6..8; raw codons run 4..6 and 7..9.
    Checking stored mutation calls against these same columns cannot expose
    the frame error. The separate coding diagnostic must expose it.
    """
    master = 'ATGAAACCCGGGTTTTAA'
    raw = 'ATGCTTAGCCCTAGGTAA'  # M L S P R *, no internal stop.
    query = 'ATGTAGCCCTAGG--TAA'
    assert '*' not in str(Seq(raw).translate())[:-1]
    assert query[3:6] == raw[5:8] == 'TAG'
    report = Audit.audit_database(minimal_frame_db(tmp_path, master, query, raw, [(2, '*')]))
    assert report['ok'], report
    assert report['counts']['mutation_residues_verified'] == 1
    warning = report['warnings']['internal_stop_in_reference_codons']
    assert warning['count'] == 1
    assert warning['examples'][0]['accession'] == 'QUERY'
    assert warning['examples'][0]['positions'] == [2, 4]
    # The raw-frame check places the shift directly: codons 3-4 sit two raw
    # bases off codons 1 and 6, an even split, so the frame is a tie. Codon 2's
    # T has two possible raw sources (4 or 5), so the stored call there is
    # unassessed rather than guessed at.
    shifted = report['warnings']['reading_frame_shift_in_raw']['examples']
    assert [(e['accession'], e['codon_runs'], e['frame_tie']) for e in shifted] == [('QUERY', [[3, 4]], True)]
    assert report['counts']['call_codons_ambiguous'] == 1
    assert 'emitted_call_not_one_raw_codon' not in report['warnings']


def test_numeric_segment_spellings_and_independent_segment_widths(red_team_db):
    with sqlite3.connect(red_team_db) as conn:
        conn.execute("UPDATE sequence_mutation_calls SET segment='1.0'")
        for table in ['features', 'sequence_alignment', 'sequences']:
            frame = pd.read_sql_query(f'SELECT * FROM {table}', conn)
            frame['segment'] = '2'
            if table == 'sequence_alignment':
                frame['alignment'] = frame['alignment'] + '---'
            frame.to_sql(table, conn, index=False, if_exists='append')
    # Identical accession strings in two segments must remain separate. Extra
    # terminal padding changes only the width, not covered/feature coordinates.
    report = Audit.audit_database(red_team_db)
    assert report['ok'], report
    assert report['counts']['alignment_rows'] == 32


def minimal_frame_db(tmp_path, master, query, raw, calls):
    """One master and one query over a single POL feature spanning every column."""
    db = tmp_path / 'frame.db'
    rows = [('MASTER', master), ('QUERY', query)]
    with sqlite3.connect(db) as conn:
        pd.DataFrame([dict(primary_accession=acc, alignment_name='MASTER', alignment=seq)
                      for acc, seq in rows]).to_sql('sequence_alignment', conn, index=False)
        pd.DataFrame([
            dict(accession=acc, master_ref_accession='MASTER', product='POL',
                 aln_start=1, aln_end=len(master), cds_start=1, cds_end=len(master),
                 cds_start_OG_seq=1, cds_end_OG_seq=len(seq.replace('-', '')))
            for acc, seq in rows
        ]).to_sql('features', conn, index=False)
        pd.DataFrame([dict(header='MASTER', sequence=master), dict(header='QUERY', sequence=raw)]
                     ).to_sql('sequences', conn, index=False)
        pd.DataFrame([dict(primary_accession='QUERY', protein_name='POL', aa_position=pos,
                           alt_residue=alt, observed_residue=alt, mutation_id=f'POL:{pos}{alt}')
                      for pos, alt in calls]).to_sql('sequence_mutation_calls', conn, index=False)
    return db


def test_frame_warnings_name_the_codons_a_stripped_insertion_breaks(tmp_path, monkeypatch):
    db, _ = build_dataset(tmp_path, monkeypatch)
    warnings = Audit.audit_database(db)['warnings']
    split = {e['accession']: e['codons'] for e in warnings['codon_split_by_raw_insertion']['examples']}
    shifted = {e['accession']: e['codon_runs'] for e in warnings['reading_frame_shift_in_raw']['examples']}
    assert split == {'STRIPPED_INS_SPLIT5': [5]}
    # PARTIAL_CODON2 lost a base in codon 2, so its six later codons set the
    # frame and codon 1 is off it. STRIPPED_INS_INFRAME gained a whole codon.
    assert shifted == {'PARTIAL_CODON2': [[1, 1]], 'STRIPPED_INS_SPLIT5': [[6, 8]],
                       'STRIPPED_INS_FRAMESHIFT': [[4, 5]]}


@pytest.mark.parametrize('reverse', [False, True])
def test_compensated_frameshift_flags_only_the_codons_between(tmp_path, reverse):
    """+1 raw base stripped after codon 2, -1 base gapped in codon 5: only codons 3-4 shift."""
    master = 'ATGAAACCCGGGTTTCTGTAA'
    query = 'ATGAAACCCGGGT-TCTGTAA'
    raw = 'ATGAAATCCCGGGTTCTGTAA'
    if reverse:
        raw = str(Seq(raw).reverse_complement())
    report = Audit.audit_database(minimal_frame_db(tmp_path, master, query, raw, [(3, 'P'), (6, 'L')]))
    assert report['ok'], report
    assert report['counts']['raw_reverse_subsequences' if reverse else 'raw_forward_subsequences'] == 1
    shifted = report['warnings']['reading_frame_shift_in_raw']['examples']
    assert [(e['accession'], e['codon_runs']) for e in shifted] == [('QUERY', [[3, 4]])]
    calls = report['warnings']['emitted_call_not_one_raw_codon']['examples']
    assert [(e['aa_position'], e['problem']) for e in calls] == [(3, 'out_of_frame')]
    assert report['counts']['call_codons_in_frame'] == 1


def test_ambiguous_raw_positions_are_unassessed_not_flagged(tmp_path):
    """An extra A inside a homopolymer could have come from any of seven bases."""
    master = 'ATGCCCGGGAAAAAATTTTAA'
    raw = 'ATGCCCGGGAAAAAAATTTTAA'
    report = Audit.audit_database(minimal_frame_db(tmp_path, master, master, raw, [(4, 'K')]))
    assert report['ok'], report
    assert 'codon_split_by_raw_insertion' not in report['warnings']
    assert 'emitted_call_not_one_raw_codon' not in report['warnings']
    assert report['counts']['frame_codons_ambiguous'] == 2
    assert report['counts']['call_codons_ambiguous'] == 1
    # Unambiguous codons downstream still show the one-base shift.
    shifted = report['warnings']['reading_frame_shift_in_raw']['examples']
    assert [e['codon_runs'] for e in shifted] == [[[6, 7]]]


# -- chaos, end to end --------------------------------------------------------

EVERY_RESIDUE = 'ACDEFGHIKLMNPQRSTVWY*-'


@pytest.mark.parametrize('family', range(12))
def test_chaos_family_never_publishes_a_residue_its_record_does_not_carry(tmp_path, monkeypatch, family):
    """Fifteen hostile rows sharing one master, every residue catalogued at every codon.

    Rows come from the modes where row and record can decide what happened
    (frame_chaos.DECIDABLE, plus benign). Any call, emitted or suppressed,
    whose residue is not the truth from the row's own event history is a wrong
    published residue. The auditor must also pass the whole database.
    """
    rng = random.Random(f'family:{family}')
    master = chaos.make_master(rng, 'repeat_hell' if family % 3 == 0 else 'chaos')
    sequence, cds_start, codons = master
    modes = ('benign',) + chaos.DECIDABLE
    cases = {f'Q{i}': chaos.build(family * 100 + i, modes[i % len(modes)], master)
             for i in range(15)}
    guide_rows = {'MASTER': sequence, 'REF_INS': sequence}
    query_rows = {'REF_INS': sequence, **{acc: case.row for acc, case in cases.items()}}
    raw_records = {'MASTER': sequence, **{acc: case.record for acc, case in cases.items()}}
    catalog = [(pos, alt) for pos in range(1, codons + 1) for alt in EVERY_RESIDUE]
    db = run_pipeline(tmp_path, monkeypatch, sequence, (cds_start + 1, cds_start + 3 * codons),
                      guide_rows, query_rows, raw_records, catalog)
    with sqlite3.connect(db) as conn:
        calls = conn.execute('SELECT primary_accession, aa_position, observed_residue, codon_read, call_status '
                             'FROM sequence_mutation_calls').fetchall()
    wrong = {}
    for acc, position, observed, read, status in calls:
        if acc in cases and (observed == 'X' or observed not in cases[acc].truth[int(position)]):
            wrong.setdefault(acc, []).append((int(position), observed, read, status,
                                              sorted(cases[acc].truth[int(position)])))
    assert not wrong, '\n\n'.join(f'{acc}: {items}\n{cases[acc].describe()}' for acc, items in wrong.items())
    assert {acc for acc, *_ in calls} & set(cases), 'nothing was called at all'
    report = Audit.audit_database(db)
    assert report['ok'], report['errors']
