Building Rabies(RABV) resource
==============================

Prerequisite Files for Building the RABV Resource
-------------------------------------------------

Building the RABV resource requires a set of prerequisite files that define how
sequences are selected, organised, mapped, and annotated during the V-gTK
database build process. These files provide the foundation for constructing a
consistent and reusable RABV genomic resource.

The main prerequisite files include:

* reference accession list
* reference tree
* clade assignment information
* sequence-to-clade or sequence-to-group mapping files
* master reference sequence information

Together, these files help V-gTK identify representative RABV sequences, assign
new or query sequences to appropriate reference groups, and support downstream
alignment, annotation, and phylogenetic workflows.

Reference Accession List
------------------------

The reference accession list contains a diverse set of GenBank accessions that
represent the known genetic diversity of rabies virus. These accessions are used
as the reference sequence set during the database build process.

The purpose of the reference set is to provide representative sequences against
which query sequences can be compared. Ideally, the reference list should include
sequences from different clades, lineages, geographic regions, and host
associations, depending on the scope of the database being built.

During the workflow, query sequences downloaded from GenBank are compared against
the reference set. If a query sequence aligns well to one of the existing
reference sequences, it can be associated with that reference group. This allows
V-gTK to organise related sequences together and maintain a structured resource.

However, if a query sequence does not align well to any sequence in the current
reference set, this may indicate that the sequence represents a divergent or
under-represented part of the dataset. In such cases, that sequence should be
considered for inclusion in the reference set. Adding it as a new reference helps
ensure that future sequences similar to it can be detected and grouped correctly.

The opposite case can also occur. A sequence may be included in the reference
set, but no query sequences align to it during the build process. In this
situation, the sequence may not currently act as an informative reference for
the available query dataset. Such sequences can be reviewed and, if appropriate,
treated as part of the query sequence set instead of the reference set.

In this way, the reference accession list is not necessarily fixed. It can be
updated and curated as new data become available, allowing the RABV resource to
remain flexible and representative over time.

Master Reference Sequence
-------------------------

In addition to the general reference accessions, the reference list should also
include a master reference sequence. The master sequence is used during the
post-alignment step for genome region annotation.

The master reference provides the coordinate framework required to map aligned
sequences to corresponding genomic regions, such as genes and coding sequences.
For RABV, this is important because the workflow needs to identify where each
sequence aligns relative to the standard genome organisation.

For example, the master reference can be used to annotate regions corresponding
to RABV genes such as:

* N gene
* P gene
* M gene
* G gene
* L gene

The master sequence should therefore be a complete and well-annotated reference
genome. For RABV, ``NC_001542`` is commonly used as the master reference
sequence.

Reference List Format
---------------------

The reference accession file is a tab-delimited file with two columns:

1. GenBank accession
2. sequence type

The sequence type should indicate whether the accession is a standard
``reference`` sequence or the ``master`` reference sequence.

Example format:

.. code-block:: text

   AB041966    reference
   AB247428    reference
   AB362483    reference
   AB383163    reference
   AB383164    reference
   AB383165    reference
   NC_001542   master

In this example, the accessions labelled ``reference`` are used as representative
RABV reference sequences. The accession labelled ``master`` is used as the main
coordinate reference for post-alignment annotation and genome-region mapping.

Important Notes
---------------

The reference accession list should be reviewed carefully before building the
database. A good reference set should be diverse enough to capture the known
variation in the dataset, but not so large that it becomes redundant or difficult
to maintain.

When curating the reference set, consider the following:

* include representative sequences from major RABV clades or lineages
* include geographically and genetically diverse sequences where possible
* ensure that the master reference is complete and correctly labelled as
  ``master``
* review query sequences that fail to align to the current reference set
* consider adding divergent unaligned query sequences as new references
* review reference sequences that do not map to any query sequences
* update the reference list as new RABV diversity becomes available

Maintaining a high-quality reference accession list is an important step in
building a reliable RABV genomic database. The quality of this file directly
affects sequence grouping, downstream annotation, clade assignment, and the
overall usability of the generated database.
