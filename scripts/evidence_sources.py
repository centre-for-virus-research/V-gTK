"""Evidence references in a mutation catalogue: what an id is and where it points.

A long-format catalogue carries one evidence reference per row. These rules are
shared by the catalogue builder (BuildHcvEvidenceCatalog.py) and the database
writer (AnnotateMutations.py), so the ids in ``mutation_catalog.evidence_id``
and the keys of the ``publications`` / ``clinical_trials`` lookup tables are
produced by the same code and always join.

Nothing here is HCV-specific. A virus with no evidence columns never calls it.
"""

import re

#: ``data_source`` vocabulary.
PUBMED = 'pubmed'
CONFERENCE_ABSTRACT = 'conference_abstract'
CLINICAL_TRIAL = 'clinical_trial'
PUBLICATION_SOURCES = (PUBMED, CONFERENCE_ABSTRACT)

#: Catalogue columns describing one evidence reference, in file order. They
#: vary between the rows of a single catalogue entry; every other column is
#: the entry itself.
EVIDENCE_COLUMNS = [
    'evidence_id', 'data_source', 'evidence_url', 'evidence_label',
    'evidence_type', 'regimens', 'linked_evidence_ids', 'finding_ids',
]

#: Separator for multi-valued cells, as everywhere else in the catalogue. A
#: comma would be ambiguous: 'Magellan-1, Part 1' is one trial name.
VALUE_SEP = ';'

_PMID = re.compile(r'^\d+$')
_NCT = re.compile(r'^NCT\d{8}$')


def publication_source(publication_id):
    """``27773808`` -> pubmed; ``AASLD_2017_Abs_1176`` -> conference_abstract."""
    return PUBMED if _PMID.match((publication_id or '').strip()) else CONFERENCE_ABSTRACT


def trial_registry_id(nct_id, curator_id):
    """The id a trial is keyed on: its NCT number, else the curator's own id.

    PHDR keys trials on a curator label, and several labels share one NCT number
    (trial arms, sponsor codes), so the registry number is the key. A trial with
    no NCT number (UMIN000015627, a Japanese UMIN-CTR registration) is still a
    trial and keeps the id it was registered under.
    """
    return (nct_id or '').strip() or (curator_id or '').strip()


def trial_url(registry_id):
    """ClinicalTrials.gov URL for an NCT id; '' for registries needing more than the id."""
    registry_id = (registry_id or '').strip()
    return f'https://clinicaltrials.gov/study/{registry_id}' if _NCT.match(registry_id) else ''
