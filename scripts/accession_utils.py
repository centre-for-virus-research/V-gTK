"""Single source of truth for accession identity.

GenBank gives every record two spellings of its name, and this pipeline uses
**both**, deliberately, for different jobs:

==========================  ==========================================
``primary_accession``       BARE - ``PV547761``. The identity key.
``locus``                   BARE - the same value.
``accession_version``       VERSIONED - ``PV547761.1``. Metadata only.
==========================  ==========================================

The bare form is what ``meta_data``, ``sequence_alignment``, ``features``, the
FASTA headers, the newick tip labels and the on-disk filenames all join on. The
versioned form exists so that :mod:`GenBankFetcher` can notice on an update run
that a record has been *revised* - ``PV547761.1`` becoming ``PV547761.2`` - and
re-fetch it. Confusing the two breaks the pipeline in opposite directions:

* a versioned value used as an identity key joins against nothing, so the
  sequence silently drops out of the alignment, the tree, or the DB merge;
* a bare value used for revision detection makes every record look unrevised,
  so an update run never picks up corrected sequences.

Before this module the bare form was produced ad hoc by ``split('.')[0]``
scattered across a dozen scripts. That idiom truncates at the *first* dot, which
is only equivalent to "drop the version" when the input happens to have exactly
one dot. It is wrong for every other shape it meets:

======================================  ================  ================
input                                   ``split('.')[0]``  this module
======================================  ================  ================
``'PV547761.1'``                        ``'PV547761'``     ``'PV547761'``
``'PV547761'``                          ``'PV547761'``     ``'PV547761'``
``'NC_001542.1.fasta'`` (a filename)    ``'NC_001542'``    ``'NC_001542'``
``'NC_001542.aligned.fasta'``           ``'NC_001542'``    ``'NC_001542'``
``'.DS_Store'``                         ``''``             ``'.DS_Store'``
``'EPI_ISL_402124'`` (GISAID)           ``'EPI_ISL_402124'`` same
``'A/swine/Iowa/4.1/1976'`` (a strain)  ``'A/swine/Iowa/4'`` unchanged
======================================  ================  ================

Two of those rows are silent data corruption rather than mere inconsistency.
``'.DS_Store'`` collapsing to ``''`` produced an empty accession that was then
used as a directory name and as the value of a command-line flag. A strain name
containing a dot - influenza strain names routinely do - was truncated
mid-name.

The rules here are therefore deliberately conservative, in the same spirit as
:mod:`segment_utils`:

* A version suffix is only stripped from a string that is *entirely* an
  accession followed by ``.<digits>``. Anything else is a label and is returned
  untouched. We never scrape at a dot just because one is present.
* Filename suffixes are only stripped when they are *recognised* suffixes.
  ``.fasta`` is an extension; ``.1`` is not.
* Absent is always ``None``, never ``''``.
"""

import os
import re


#: Accession shapes GenBank/RefSeq actually issue: one to three letters, an
#: optional RefSeq underscore, then the digits. ``NC_001542``, ``KX148218``,
#: ``PV547761``, ``AF123456``. Anchored on purpose - a bare alphanumeric token
#: that does not match this is a *label* (a strain name, a cluster id, a
#: segment group), and labels must survive verbatim.
_ACCESSION_RE = re.compile(r'^[A-Z]{1,3}_?\d{4,9}$', re.IGNORECASE)

#: The same, carrying a version suffix. The version is bounded: no GenBank
#: record has a four-digit version, and without the bound a long digit run
#: would be mistaken for a version rather than part of a label.
_VERSIONED_RE = re.compile(r'^([A-Z]{1,3}_?\d{4,9})\.(\d{1,3})$', re.IGNORECASE)

#: Suffixes the pipeline actually writes. Stripped right to left, repeatedly,
#: so ``NC_001542.aligned.fasta`` and ``NC_001542.fasta.gz`` both reduce
#: correctly. ``.1`` is deliberately absent: a version is not a file extension
#: and is removed separately, only when the remainder is accession-shaped.
_KNOWN_SUFFIXES = frozenset({
    '.fasta', '.fa', '.fna', '.faa', '.fas', '.seq',
    '.aln', '.aligned', '.padded', '.dedup', '.unique',
    '.tsv', '.csv', '.txt', '.json', '.xml', '.gb', '.gbk', '.gff', '.gff3',
    '.vcf', '.nwk', '.newick', '.tree', '.treefile', '.log',
    '.gz', '.bz2', '.zip',
})


def normalise_accession(value):
    """Return the canonical (bare) identity key for ``value``, or ``None``.

    Only an accession-shaped string loses a version suffix. Everything else is
    a label and comes back with nothing but surrounding whitespace removed.

    >>> normalise_accession('PV547761.1')
    'PV547761'
    >>> normalise_accession('PV547761')
    'PV547761'
    >>> normalise_accession('  NC_001542.2  ')
    'NC_001542'
    >>> normalise_accession('EPI_ISL_402124')       # GISAID, not GenBank
    'EPI_ISL_402124'
    >>> normalise_accession('A/swine/Iowa/4.1/1976')  # a strain name
    'A/swine/Iowa/4.1/1976'
    >>> normalise_accession('cluster.1')            # a label, not an accession
    'cluster.1'
    >>> normalise_accession('  ') is None
    True
    """
    if value is None:
        return None
    text = str(value).strip()
    if not text:
        return None

    versioned = _VERSIONED_RE.match(text)
    if versioned:
        return versioned.group(1)
    return text


def split_accession_version(value):
    """Split ``value`` into ``(bare_accession, version_or_None)``.

    The inverse view of :func:`normalise_accession`, for the one caller that
    genuinely needs the version: detecting a revised record on an update run.

    >>> split_accession_version('PV547761.1')
    ('PV547761', 1)
    >>> split_accession_version('PV547761')
    ('PV547761', None)
    >>> split_accession_version('cluster.1')
    ('cluster.1', None)
    >>> split_accession_version(None)
    (None, None)
    """
    if value is None:
        return None, None
    text = str(value).strip()
    if not text:
        return None, None

    versioned = _VERSIONED_RE.match(text)
    if versioned:
        return versioned.group(1), int(versioned.group(2))
    return text, None


def is_accession(value):
    """Does ``value`` look like a GenBank/RefSeq accession, versioned or not?

    >>> is_accession('NC_001542'), is_accession('NC_001542.1')
    (True, True)
    >>> is_accession('EPI_ISL_402124'), is_accession('')
    (False, False)
    """
    if value is None:
        return False
    text = str(value).strip()
    return bool(_ACCESSION_RE.match(text) or _VERSIONED_RE.match(text))


def strip_known_suffixes(name):
    """Remove recognised file extensions from ``name``, right to left.

    Unrecognised suffixes are left alone, which is the whole point: ``.1`` is
    not an extension.

    >>> strip_known_suffixes('NC_001542.fasta')
    'NC_001542'
    >>> strip_known_suffixes('NC_001542.aligned.fasta')
    'NC_001542'
    >>> strip_known_suffixes('NC_001542.1.fasta')
    'NC_001542.1'
    >>> strip_known_suffixes('NC_001542.errors.csv')
    'NC_001542.errors'
    >>> strip_known_suffixes('.DS_Store')
    '.DS_Store'
    """
    text = str(name)
    while True:
        base, dot, suffix = text.rpartition('.')
        if not dot or not base:
            # No dot left, or the dot is leading (a dotfile) - stop.
            return text
        if ('.' + suffix).casefold() not in _KNOWN_SUFFIXES:
            return text
        text = base


def accession_from_filename(file_path):
    """The canonical accession a pipeline file is named after.

    Filenames in ``tmp/Blast/{grouped_fasta,ref_seqs,master_seq}`` and the
    nextalign output tree are ``<accession><ext>``. This is what turns one back
    into the identity key that the reference list, the update DB and the
    alignment headers all use.

    >>> accession_from_filename('/tmp/x/NC_001542.fasta')
    'NC_001542'
    >>> accession_from_filename('NC_001542.1.fasta')
    'NC_001542'
    >>> accession_from_filename('NC_001542.1.aligned.fasta')
    'NC_001542'
    >>> accession_from_filename('KX148218')
    'KX148218'
    >>> accession_from_filename('.DS_Store')      # not silently empty
    '.DS_Store'
    >>> accession_from_filename('segment_4_plus.fa')
    'segment_4_plus'
    """
    name = os.path.basename(str(file_path))
    if not name:
        return ''
    stripped = strip_known_suffixes(name)
    normalised = normalise_accession(stripped)
    return stripped if normalised is None else normalised


def normalise_accession_series(values):
    """:func:`normalise_accession` over an iterable, preserving order.

    ``None`` and blanks become ``''`` rather than being dropped, so the result
    lines up positionally with the input - callers are usually assigning back
    into a dataframe column.

    >>> normalise_accession_series(['PV547761.1', None, 'cluster.1'])
    ['PV547761', '', 'cluster.1']
    """
    out = []
    for value in values:
        normalised = normalise_accession(value)
        out.append('' if normalised is None else normalised)
    return out
