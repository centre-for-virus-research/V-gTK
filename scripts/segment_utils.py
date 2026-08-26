"""Single source of truth for segment identity.

Before this module six scripts each carried their own ``_normalize_segment``
and they did not agree. The divergences were not cosmetic:

======================  ========  =========  =========  ==========  ==========
input                   GenBank-  Blast-     Pad-       Nextalign   SegmentPivot
                        Parser    Alignment  Alignment
======================  ========  =========  =========  ==========  ==========
``'04'``                ``'04'``  ``'04'``   ``'4'``    ``'4'``     ``'4'``
``'4.0'``               ``'4.0'`` ``'40'``   ``'40'``   ``'40'``    ``'4'``
``'segment 7'``         verbatim  ``'7'``    ``'7'``    ``'7'``     verbatim
``'RNA 4'``             verbatim  ``'4'``    ``'4'``    ``'4'``     verbatim
``'PB2'`` (segment 1)   verbatim  ``'2'``    ``'2'``    ``'2'``     verbatim
``'PB1'`` (segment 2)   verbatim  ``'1'``    ``'1'``    ``'1'``     verbatim
missing                 ``''``    ``'0'``    ``None``   ``''``      ``None``
======================  ========  =========  =========  ==========  ==========

Three of those rows are silent data corruption rather than mere inconsistency:

* ``'4.0'`` became segment **40**, because the old rule scraped every digit out
  of the string instead of parsing it.
* ``'PB2'`` became segment **2** and ``'PB1'`` became segment **1** - exactly
  inverted, because digit-scraping picks up the digit inside a gene name.
  PB2 is segment 1 and PB1 is segment 2. ``M1``/``M2`` (both segment 7) and
  ``NS1`` (segment 8) were mangled the same way.
* missing became ``'0'`` in BlastAlignment, which reads downstream as a real
  segment rather than as an absent one.

The rules here are therefore deliberately conservative:

* Digits are only taken from a string that is *entirely* a number, or that is a
  recognised ``<word> <number>`` form such as ``'segment 7'`` / ``'RNA 4'``.
  A bare alphanumeric token like ``PB1`` is a **label**, never a number.
* Absent is always ``None``. Never ``''``, never ``'0'``.
* Whether a label counts as "no segment" is a question about the *virus*, not
  about the string, so it is answered by :func:`is_unavailable` against the
  reference list's valid label set - not by a hardcoded sentinel list. This
  matters because ``NA`` is a legitimate influenza segment name (neuraminidase,
  segment 6) while also being a conventional spelling of "missing".
"""

import csv
import os
import re


#: Spellings of "no value here" that are never a real segment label. Used only
#: when the caller cannot supply the virus's valid label set.
FALLBACK_MISSING_TOKENS = frozenset(
    {'', 'not found', 'na', 'n/a', 'none', 'null', 'nan', 'unknown', '-'}
)

#: Strings that are only ever the *stringified form of a null*, never a real
#: label. Distinct from :data:`FALLBACK_MISSING_TOKENS` because that set contains
#: ``na``, which is a genuine influenza segment name. Call sites that receive
#: values via pandas need to discard these without discarding neuraminidase.
PANDAS_NULL_TOKENS = frozenset({'nan', 'none', '<na>', 'null'})

#: ``<word><sep><number>`` forms seen in GenBank ``/segment=`` qualifiers.
#: Deliberately an explicit word list: scraping digits out of arbitrary text is
#: what inverted PB1/PB2. ``s``/``S`` is excluded because a bare ``S`` is a real
#: segment name for several bunyaviruses.
_PREFIXED_NUMBER_RE = re.compile(
    r'^(?:segment|segments|seg|sgt|rna|vrna|gene)[\s._:#-]*(\d+)$',
    re.IGNORECASE,
)

#: A string that is wholly a number, optionally with a redundant decimal part.
#: ``4``, ``04``, ``4.0`` all mean segment 4. ``4.5`` does not and is rejected.
#: Length-bounded: no genome has a 10-digit segment, and without the bound a long
#: digit string reaches ``int(float(text))``, where float() returns inf and int()
#: raises OverflowError - breaking the "never raises" contract these callers rely
#: on. Anything longer is treated as a label, which is what it is.
_WHOLE_NUMBER_RE = re.compile(r'^\d{1,9}(?:\.0*)?$')


def normalise_segment(value):
    """Return the canonical form of a segment identifier, or ``None``.

    Purely syntactic. It answers "what is this string, written canonically",
    not "does this virus have such a segment" - see :func:`is_unavailable`
    for the latter.

    >>> normalise_segment('04')
    '4'
    >>> normalise_segment('4.0')
    '4'
    >>> normalise_segment('Segment 7')
    '7'
    >>> normalise_segment('RNA 4')
    '4'
    >>> normalise_segment('PB1')          # a label, not the number 1
    'PB1'
    >>> normalise_segment('  ') is None
    True
    """
    if value is None:
        return None
    text = str(value).strip()
    if not text:
        return None

    if _WHOLE_NUMBER_RE.match(text):
        return str(int(float(text)))

    prefixed = _PREFIXED_NUMBER_RE.match(text)
    if prefixed:
        return str(int(prefixed.group(1)))

    # Anything else is a label. Preserve it verbatim (minus surrounding space)
    # so that L/M/S, HA, PB1 and friends survive intact.
    return text


def is_unavailable(value, valid_labels=None):
    """Is ``value`` "no segment was assigned"?

    ``valid_labels`` should be the set of labels the reference list declares for
    this virus (case-folded). When supplied it is authoritative, which is what
    lets influenza's ``NA`` segment be a real segment while ``NA`` in a virus
    that uses numeric labels still reads as missing.

    Without ``valid_labels`` we fall back to :data:`FALLBACK_MISSING_TOKENS`,
    preserving the historical behaviour of callers that have no reference list.
    """
    normalised = normalise_segment(value)
    if normalised is None:
        return True
    if valid_labels is not None:
        return normalised.casefold() not in valid_labels
    return normalised.casefold() in FALLBACK_MISSING_TOKENS


def normalise_or_none(value, valid_labels=None):
    """:func:`normalise_segment`, but collapsing "unavailable" to ``None``."""
    if is_unavailable(value, valid_labels):
        return None
    return normalise_segment(value)


def load_segment_names(path):
    """Load a segment-name mapping asset.

    The file is a TSV with a header and columns ``segment``, ``segment_name``
    and an optional comma-separated ``aliases``. Returns
    ``(by_segment, by_name)`` where ``by_segment`` maps canonical segment ->
    display name and ``by_name`` maps a case-folded name *or alias* -> canonical
    segment.

    ``keep_default_na=False`` semantics are inherent here: this is hand-parsed
    with :mod:`csv` precisely so that influenza's ``NA`` segment name survives.
    pandas would read it as a null.
    """
    by_segment = {}
    by_name = {}
    if not path or not os.path.isfile(path):
        return by_segment, by_name

    with open(path, 'r', encoding='utf-8', newline='') as handle:
        for row in csv.DictReader(handle, delimiter='\t'):
            segment = normalise_segment(row.get('segment'))
            name = (row.get('segment_name') or '').strip()
            if segment is None or not name:
                continue
            by_segment[segment] = name
            by_name[name.casefold()] = segment
            for alias in (row.get('aliases') or '').split(','):
                alias = alias.strip()
                if alias:
                    by_name.setdefault(alias.casefold(), segment)
    return by_segment, by_name


def segment_to_name(segment, by_segment):
    """Display name for a segment, or ``''`` when unknown.

    Non-numeric segment labels are already names (Lassa's ``L``/``S``), so they
    are returned unchanged rather than looked up and lost.
    """
    normalised = normalise_segment(segment)
    if normalised is None:
        return ''
    if normalised in by_segment:
        return by_segment[normalised]
    if not normalised.isdigit():
        return normalised
    return ''


def name_to_segment(name, by_name):
    """Canonical segment for a declared segment name, or ``None``.

    This is what lets a submitter-declared ``HA`` be reconciled with the
    pipeline's inferred ``4``.
    """
    if name is None:
        return None
    text = str(name).strip()
    if not text:
        return None
    return by_name.get(text.casefold())
