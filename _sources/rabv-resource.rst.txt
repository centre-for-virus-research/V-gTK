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

The reference accession list contains a diverse set of accessions that
represent the known genetic diversity of rabies virus. These accessions are used
as the reference sequence set during the database build process.

The purpose of the reference set is to provide representative sequences against
which query sequences can be compared. Ideally, the reference list should include
sequences from different clades, lineages, geographic regions, and host
associations, depending on the scope of the database being built.

During the workflow, if a query sequence does not align well to any sequence in the current
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

Reference Tree
--------------

The reference tree is one of the key prerequisite files required for building
the RABV resource. It contains a phylogenetic tree generated from all selected
reference sequences. These reference sequences should represent the known
genetic diversity of rabies virus and provide the framework used for downstream
clade assignment.

For RABV, the reference tree should cover the major and minor diversity present
within the virus. For example, RABV can be represented by several major clades
and many minor clades. The reference tree should therefore be built from a
carefully selected and diverse reference sequence set so that the major
evolutionary structure of the dataset is captured.

The reference tree acts as the main phylogenetic backbone of the database. New
or query sequences can later be compared against this reference structure to
support clade assignment, lineage interpretation, and downstream phylogenetic
analysis.

Importance of the Reference Tree
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The reference tree is used as the primary tree for assigning major and minor
clade information. Because of this, it is important that the tree is generated
from a representative set of reference sequences.

A well-curated reference tree should:

* include sequences covering the known RABV diversity
* represent all major clades where possible
* include representative sequences from minor clades
* avoid unnecessary redundancy
* be generated using an appropriate phylogenetic model
* be linked to clear metadata describing how the tree was created

If the reference tree does not adequately represent the diversity of the
available RABV sequences, newly added or query sequences may be incorrectly
placed or remain difficult to classify. Therefore, maintaining a diverse and
well-constructed reference tree is mandatory for reliable database construction
and clade assignment.

Reference Tree Metadata
~~~~~~~~~~~~~~~~~~~~~~~

Along with the tree file, a metadata file must also be provided. This metadata
describes the tree file and provides information such as the tree name, tree
type, model used, and the genome segment or chromosome to which the tree
belongs.

The metadata file allows V-gTK to identify which tree should be used as the main
reference tree and how it should be interpreted during the workflow.

An example tree metadata file is shown below:

.. code-block:: text

   chromosome    segment_number    tree_type    tree_name                 tree_model
   1             1                 reference    ref_tree_am3c_am5.treefile  GTR+F+R6

In this example:

* ``chromosome`` indicates the chromosome or genome component.
* ``segment_number`` indicates the segment number. For RABV, this is usually
  ``1`` because RABV has a non-segmented genome.
* ``tree_type`` defines the type of tree. The value ``reference`` indicates
  that this is the main reference tree.
* ``tree_name`` provides the filename of the tree file.
* ``tree_model`` records the phylogenetic model used to generate the tree, if
  available.

The ``tree_type`` field is particularly important because the tree labelled as
``reference`` is treated as the main tree for clade assignment.

Additional Trees
~~~~~~~~~~~~~~~~

In addition to the main reference tree, users can also include other
phylogenetic trees in the resource. These additional trees can be useful for
more focused analyses.

Examples include:

* a complete tree containing a larger set of available sequences
* a region-specific tree
* a country-specific tree
* a host-specific tree
* a clade-specific tree

For example, a tree generated using Latin American RABV sequences can be added
to the database. This type of regional tree can help users analyse newly
generated sequences from Latin America in greater detail. While the main
reference tree provides broad clade assignment, the regional tree can provide
additional resolution within a specific geographic or evolutionary context.

Similarly, users may add trees focused on particular clades, hosts, or outbreak
datasets depending on the purpose of the database.

Reference Tree Versus Additional Trees
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

It is important to distinguish the main reference tree from additional trees.

The ``reference`` tree is the primary tree used for assigning major and minor
clades. This tree must be diverse and carefully curated because it forms the
main phylogenetic framework for the RABV resource.

Additional trees are optional and are mainly used to support specialised
analysis. They can provide extra context, but they do not replace the main
reference tree unless they are explicitly defined as the reference tree in the
metadata.

For this reason, the main reference tree should always be built from a diverse
set of representative reference sequences. This ensures that the database has a
stable and reliable framework for classifying new RABV sequences.

Example Tree Metadata with Multiple Trees
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

An extended metadata file may include more than one tree:

.. code-block:: text

   chromosome    segment_number    tree_type       tree_name                         tree_model
   1             1                 reference       ref_tree_am3c_am5.treefile        GTR+F+R6
   1             1                 complete        complete_rabv_tree.treefile       GTR+F+R6
   1             1                 regional        latin_america_rabv.treefile       GTR+F+R6

In this example, ``ref_tree_am3c_am5.treefile`` is the main reference tree used
for clade assignment. The complete and regional trees provide additional
phylogenetic context for downstream analysis.

Important Notes
~~~~~~~~~~~~~~~

When preparing reference tree files, consider the following:

* the reference tree should be built from the curated reference accession set
* the tree should represent the full diversity expected in the database
* the metadata file must correctly identify the main tree as ``reference``
* additional trees can be included for specialised downstream analysis
* regional or clade-specific trees can improve interpretation within specific
  sequence groups
* the tree model should be recorded where possible
* the tree filename in the metadata must match the actual tree file name

A diverse and well-documented reference tree is essential for building a robust
RABV genomic database. Since major and minor clade assignments are based on the
reference tree, the quality and diversity of this tree directly affect the
accuracy and usefulness of the final database.

Mapping Files
-------------

Mapping files are used to standardise and correct metadata values during the
RABV database build process. They are especially useful when the same incorrect,
inconsistent, or non-standard metadata values occur many times across hundreds
or thousands of records.

For example, country names may appear in different formats such as ``USA``,
``United States of America``, or ``US``. Similarly, host names may be recorded
using common names, incomplete names, spelling variations, or extended
descriptions. Mapping files allow these values to be replaced in bulk with a
single standardised value.

This helps keep the database clean, consistent, and easier to query.

Purpose of Mapping Files
~~~~~~~~~~~~~~~~~~~~~~~~

Metadata downloaded from public databases such as GenBank can contain variation
in how fields are recorded. The same country, region, or host species may be
written in different ways by different submitters.

Mapping files provide a simple way to correct these inconsistencies before the
final database is generated.

They can be used to:

* correct country names
* standardise host names
* replace outdated geographic names
* convert abbreviations into full names
* fix spelling variations
* merge multiple equivalent terms into one accepted value
* improve consistency across large datasets

Country Mapping File
~~~~~~~~~~~~~~~~~~~~

The country mapping file is used to standardise country names. It contains two
columns:

1. the original value present in the metadata
2. the corrected or standardised value to replace it with

Example:

.. code-block:: text

   country          replaced_by
   USA              United States
   Czechoslovakia   Czechia
   Laos             Lao
   UK               United Kingdom
   CSRF             Czechia


Host Mapping File
~~~~~~~~~~~~~~~~~

The host mapping file is used to standardise host information. Like the country
mapping file, it contains two columns:

1. the original host value present in the metadata
2. the corrected or standardised host value to replace it with

Example:

.. code-block:: text

   host                replaced_by
   Indian Gaur         Bos gaurus
   bovine; male        Bos taurus
   bovine; female      Bos taurus
   Bos taurus taurus   Bos taurus


Clade Assignment Files
----------------------

Clade assignment files are used to assign major and minor clade information to
query sequences during the RABV database build and analysis workflow. These files
provide the reference framework required for classifying new or unassigned RABV
sequences based on their relationship to known reference sequences.

To perform clade assignment, a curated reference sequence set is required. This
can be the same reference sequence set used to build the reference tree. Each
reference sequence should have an associated major clade and, where available, a
minor clade. These assignments are then used by V-gTK to interpret where query
sequences belong within the RABV diversity framework.

For RABV, clade assignment is usually organised at two levels:

* **major clade**: broad evolutionary group
* **minor clade**: more specific sub-group within a major clade

Major Clade Assignment
~~~~~~~~~~~~~~~~~~~~~~

The major clade file defines the broad clade membership of each reference
sequence. Each sequence should be assigned to one major clade based on its
phylogenetic placement and curated classification.

The file should contain two columns:

1. ``sequence_id``
2. ``major_clade``

Example:

.. code-block:: text

   sequence_id    major_clade
   AB041966       Indian-Sub
   AB247428       Bats
   AB362483       Cosmopolitan
   AB383163       Bats
   AB383164       Bats


Minor Clade Assignment
~~~~~~~~~~~~~~~~~~~~~~

The minor clade file provides a more detailed classification within the major
clade structure. Minor clades represent finer-scale phylogenetic groups and can
be useful for more detailed epidemiological, geographic, or evolutionary
interpretation.

The file should contain two columns:

1. ``sequence_id``
2. ``minor_clade``

Example:

.. code-block:: text

   sequence_id    minor_clade
   AB041966       NULL
   AB247428       DR
   AB362483       AM3b
   AB383163       LC
   AB383164       LC

Important Notes
~~~~~~~~~~~~~~~

When preparing clade assignment files, consider the following:

* each reference sequence should have a major clade assignment
* minor clade assignment can be set to ``NULL`` if unavailable
* sequence identifiers must match across the reference files, tree, and clade
  files
* major and minor clades should be curated using the reference tree
* clade names should be written consistently across all files
* avoid spelling differences or duplicated labels for the same clade
* update clade files when new reference sequences are added
* review clade assignments after updating or rebuilding the reference tree

Clade assignment depends heavily on the quality of the reference sequence set and
reference tree. Therefore, the reference tree should include a diverse set of
sequences covering all known major and minor RABV clades where possible.
