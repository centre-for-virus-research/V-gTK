Functions
=========

GenBank XML Downloader
----------------------

This script downloads and updates GenBank XML files for a given species.  
It supports batch downloads, update mode, and configurable directories.

Usage
^^^^^

.. code-block:: bash

   python script.py --taxid <TAXID> [options]

Arguments
^^^^^^^^^

.. option:: -t, --taxid <TAXID>

   **(Required)** Taxonomic ID of the species.  
   Example: ``11292``

.. option:: -o, --tmp_dir <DIR>

   Temporary output directory where XML files are stored.  
   Default: ``tmp``

.. option:: -b, --batch_size <N>

   Maximum number of XML files to download and merge into a single file.  
   Default: ``100``

.. option:: --update

   Run in *update mode* to download only new sequences.  
   Expects a TSV file with a column named ``Accession Version``.

.. option:: -u, --base_url <URL>

   Base URL used to download the XML files.  
   Default: ``https://eutils.ncbi.nlm.nih.gov/entrez/eutils/``

.. option:: -e, --email <EMAIL>

   Email ID for NCBI requests.  
   Default: ``your_email@example.com``

.. option:: -s, --sleep_time <N>

   Delay (in seconds) between fetch requests.  
   Default: ``2``

.. option:: -d, --base_dir <DIR>

   Directory where all GenBank XML files are stored.  
   Default: ``GenBank-XML``

Example
^^^^^^^

.. code-block:: bash

   python script.py --taxid 11292 --email myname@domain.com --batch_size 200 --update


GenBank XML Parser
---------------------------

This command extracts information from GenBank XML files and writes a tab-separated table (TSV).

Usage
^^^^^

.. code-block:: bash

   python GenBankParser.py --ref_list <REFS_TSV> [options]

Arguments
^^^^^^^^^

.. option:: -d, --input_dir <DIR>

   Input directory containing GenBank XML files.  
   Default: ``tmp/GenBank-XML``

.. option:: -b, --base_dir <DIR>

   Base working directory used by the extractor.  
   Default: ``tmp``

.. option:: -o, --output_dir <DIR>

   Output directory where TSV tables are written.  
   Default: ``GenBank-matrix``

.. option:: -r, --ref_list <FILE>

   **(Required)** File containing a set of reference accessions.

.. option:: -e, --exclusion_list <FILE>

   Optional file with accessions to exclude from the output.

.. option:: -s, --is_segmented_virus {Y,N}

   Whether the virus is segmented (``Y`` or ``N``).  
   Default: ``N``

Examples
^^^^^^^^

Basic run:

.. code-block:: bash

   python GenBankParser.py \
     --ref_list refs.tsv

Custom directories and segmented virus:

.. code-block:: bash

   python GenBankParser.py \
     --input_dir tmp/GenBank-XML \
     --base_dir tmp \
     --output_dir GenBank-matrix \
     --ref_list refs.tsv \
     --exclusion_list exclude.tsv \
     --is_segmented_virus Y

Notes
^^^^^

- Make sure your ``--ref_list`` file contains valid accession identifiers (one per line, or as expected by your parser).
- When ``--is_segmented_virus Y`` is used, the output may include one row per segment depending on your parsing logic.


Download GFF3 from NCBI
-----------------------

The script ``DownloadGFF.py`` downloads **GFF3 files** from NCBI using *efetch*.  
It accepts one or more accession IDs and saves the results locally.

Usage
^^^^^

.. code-block:: bash

   python DownloadGFF.py --accession_ids <ACCESSION> [options]

Arguments
^^^^^^^^^

.. option:: -id, --accession_ids <ID>

   **(Required)** One or more NCBI accession IDs to download.  
   Example: ``NC_002645``

.. option:: -b, --base_dir <DIR>

   Base working directory.  
   Default: ``tmp``

.. option:: -o, --output_dir <DIR>

   Directory where the downloaded GFF3 files are saved.  
   Default: ``Gff``

Examples
^^^^^^^^

Download a single accession:

.. code-block:: bash

   python DownloadGFF.py --accession_ids NC_002645

Download into a custom directory:

.. code-block:: bash

   python DownloadGFF.py \
     --accession_ids NC_002645 \
     --base_dir data \
     --output_dir ncbi_gff

Validate GenBank Matrix
-----------------------

The script ``ValidateMatrix.py`` validates a GenBank matrix (``gB_matrix``) using
country, date, and host columns. It can download the NCBI taxonomy dump to
standardize names and produces a cleaned/validated TSV.

Usage
^^^^^

.. code-block:: bash

   python ValidateMatrix.py [options]

Arguments
^^^^^^^^^

.. option:: -u, --url <URL>

   URL to download the NCBI taxonomy dump (``taxdump.tar.gz``).  
   Default: ``https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz``

.. option:: -t, --taxa_path <DIR>

   Directory to store unpacked taxonomy files (names, nodes, etc.).  
   Default: ``Taxa``

.. option:: -b, --base_dir <DIR>

   Base working directory.  
   Default: ``tmp-test``

.. option:: -o, --output_dir <DIR>

   Output directory for validation results.  
   Default: ``Validate-matrix``

.. option:: -g, --gb_matrix <FILE>

   Input GenBank matrix TSV file.  
   Default: ``tmp-test/GenBank-matrix/gB_matrix_raw.tsv``

.. option:: -c, --country <FILE>

   M49 country reference CSV used for normalization/validation.  
   Default: ``assets/m49_country.csv``

.. option:: -a, --assets <DIR>

   Directory containing auxiliary assets (e.g., reference CSVs).  
   Default: ``assets/``

.. option:: --no-debug

   Disable verbose/fuzzy debug prints.  
   (Debug logging is **enabled by default**; this flag turns it off.)

Examples
^^^^^^^^

Basic validation with defaults:

.. code-block:: bash

   python ValidateMatrix.py

Specify custom locations for inputs/outputs:

.. code-block:: bash

   python ValidateMatrix.py \
     --gb_matrix data/GenBank-matrix/gB_matrix_raw.tsv \
     --country assets/m49_country.csv \
     --base_dir work \
     --output_dir Validate-matrix

Run quietly (no fuzzy debug prints):

.. code-block:: bash

   python ValidateMatrix.py --no-debug

Notes
^^^^^

- By default, the script enables debug messages; use ``--no-debug`` to suppress them.
- Ensure the input matrix contains the expected columns (country, collection date,
  host) so validation rules can be applied effectively.
- The taxonomy dump will be downloaded to ``--taxa_path`` if not present.


Add Missing Data to GenBank Matrix
----------------------------------

The script ``AddMissingData.py`` processes a GenBank matrix to handle missing
values and apply bulk replacements. It can use a *fillup file* to insert
additional information and a *bulk file* for systematic replacements.

Usage
^^^^^

.. code-block:: bash

   python AddMissingData.py [options]

Arguments
^^^^^^^^^

.. option:: -d, --tmp_dir <DIR>

   Directory where the validated GenBank matrix is saved.  
   Recommended: ``tmp/AddMissingData``  
   Default: ``tmp/AddMissingData``

.. option:: -g, --gb_matrix <FILE>

   Input GenBank matrix TSV generated by ``genbank_to_tsv.py``.  
   Default: ``tmp/GenBank-matrix/gB_matrix_raw.tsv``

.. option:: -f, --fillup_file <FILE>

   TSV file containing additional metadata to fill missing values.  
   The file **must** contain a primary accession column; other columns are optional.  
   Default: ``None``

.. option:: -b, --bulk_file <FILE>

   TSV file specifying bulk replacements with ``host`` and ``replaced_by`` columns.  
   Default: ``None``

Examples
^^^^^^^^

Run with default paths:

.. code-block:: bash

   python AddMissingData.py

Provide a fillup file to supplement missing metadata:

.. code-block:: bash

   python AddMissingData.py \
     --fillup_file data/fillup.tsv

Apply bulk replacements in addition to fillups:

.. code-block:: bash

   python AddMissingData.py \
     --gb_matrix data/gB_matrix_raw.tsv \
     --fillup_file data/fillup.tsv \
     --bulk_file data/bulk_replace.tsv \
     --tmp_dir tmp/AddMissingData

Notes
^^^^^

- The **fillup file** helps enrich the GenBank matrix with missing metadata.  
- The **bulk file** allows replacing multiple values systematically.  
- Output is written into the directory specified by ``--tmp_dir``.

Filter and Extract Sequences
----------------------------

The script ``FilterAndExtractSequence.py`` filters sequences from a GenBank matrix
and prepares them for downstream BLAST alignment. It supports filtering by sequence
length, ambiguity, GenBank divisions, and reference sets.

Usage
^^^^^

.. code-block:: bash

   python FilterAndExtractSequence.py --genbank_matrix <FILE> --ref_file <FILE> [options]

Arguments
^^^^^^^^^

.. option:: -g, --genbank_matrix <FILE>

   **(Required)** Input GenBank matrix file.

.. option:: -sf, --sequence_file <FILE>

   Input FASTA file containing sequences.  
   Default: ``tmp/GenBank-matrix/sequences.fa``

.. option:: -f, --genbank_matrix_filtered <DIR>

   Directory where the filtered GenBank matrix will be stored.  
   Default: ``tmp/GenBank-matrix``

.. option:: -r, --ref_file <FILE>

   **(Required)** Text file containing a list of reference sequence accessions.

.. option:: -b, --base_dir <DIR>

   Base working directory.  
   Default: ``tmp``

.. option:: -o, --output_dir <DIR>

   Directory where processed sequences will be written.  
   Default: ``Sequences``

.. option:: -l, --total_length <N>

   Minimum sequence length (including ambiguous characters).  
   Default: ``1``

.. option:: -n, --real_length <N>

   Minimum sequence length excluding ambiguous ``N`` bases.  
   Default: ``1``

.. option:: -a, --prop_ambigious_data <VALS>

   Proportion(s) of ambiguous data to exclude.  
   Accepts one or more values.  
   Default: ``None``

.. option:: -v, --segmented_virus <Y/N>

   Whether the dataset is from a segmented virus.  
   Default: ``N``

.. option:: -d, --genbank_division <LIST>

   GenBank divisions to exclude from analysis.  
   Default: ``["VRL", "PAT", "SYN", "ENV"]``

   Common divisions include:

   - ``VRL`` — Viral sequences  
   - ``ENV`` — Environmental sequences  
   - ``PAT`` — Patented sequences  
   - ``SYN`` — Synthetic/chimeric sequences  

.. option:: -vd, --valid_divisions <LIST>

   GenBank divisions to **include** for analysis.  
   Default: ``["VRL", "ENV"]``

.. option:: -s, --seq_type <TYPE>

   Sequence type (e.g., nucleotide, protein).  
   Default: ``None``

Examples
^^^^^^^^

Filter sequences longer than 5000 bases:

.. code-block:: bash

   python FilterAndExtractSequence.py \
     --genbank_matrix data/gB_matrix.tsv \
     --ref_file refs.txt \
     --total_length 5000

Exclude environmental and patented sequences:

.. code-block:: bash

   python FilterAndExtractSequence.py \
     --genbank_matrix data/gB_matrix.tsv \
     --ref_file refs.txt \
     --genbank_division VRL PAT SYN ENV

Filter by valid divisions and only keep nucleotide sequences:

.. code-block:: bash

   python FilterAndExtractSequence.py \
     --genbank_matrix data/gB_matrix.tsv \
     --ref_file refs.txt \
     --valid_divisions VRL \
     --seq_type nucleotide

Notes
^^^^^

- Both ``--genbank_division`` and ``--valid_divisions`` can be provided as lists
  (space-separated).  
- ``--total_length`` and ``--real_length`` help enforce sequence quality.  
- ``--prop_ambigious_data`` allows filtering sequences with excessive ambiguity.

BLAST Alignment
---------------

The script ``BlastAlignment.py`` performs BLAST alignments of query sequences
against a given set of reference sequences. It can create a BLAST database
from the provided reference FASTA file and produces a table of top hits.

Usage
^^^^^

.. code-block:: bash

   python BlastAlignment.py --master_acc <ACCESSION> [options]

Arguments
^^^^^^^^^

.. option:: -q, --query_fa <FILE>

   Query FASTA file containing sequences to align.  
   Default: ``tmp/Sequences/query_seq.fa``

.. option:: -r, --ref_fa <FILE>

   Reference FASTA file used to create a BLAST database.  
   (This file will be indexed for BLAST searches.)  
   Default: ``tmp/Sequences/ref_seq.fa``

.. option:: -b, --base_dir <DIR>

   Base working directory.  
   Default: ``tmp``

.. option:: -t, --output_dir <DIR>

   Output directory for BLAST results.  
   Default: ``Blast``

.. option:: -o, --output_file <FILE>

   Output file name for BLAST top hits.  
   Default: ``query_tophits.tsv``

.. option:: -s, --is_segmented_virus {Y,N}

   Whether the virus is segmented (``Y`` or ``N``).  
   Default: ``N``

.. option:: -f, --segment_file <FILE>

   File containing information about viral segments.  
   Required only when ``--is_segmented_virus Y``.

.. option:: -m, --master_acc <ACCESSION>

   **(Required)** Master reference accession.  
   Example: Rabies virus uses ``NC_001542``.

.. option:: -u, --is_update {Y,N}

   If set to ``Y``, skips BLAST on sequences that were already aligned,
   useful for incremental updates.  
   Default: ``N``

.. option:: -k, --keep_blast_tmp_dir {Y,N}

   Whether to retain the temporary BLAST working directory for debugging.  
   Default: ``N``

.. option:: -g, --gb_matrix <FILE>

   Input GenBank matrix file.  
   Default: ``tmp/GenBank-matrix/gB_matrix_raw.tsv``

Examples
^^^^^^^^

Run a basic BLAST alignment:

.. code-block:: bash

   python BlastAlignment.py \
     --master_acc NC_001542 \
     --query_fa data/query.fa \
     --ref_fa data/ref.fa

Run in update mode to skip previously aligned sequences:

.. code-block:: bash

   python BlastAlignment.py \
     --master_acc NC_001542 \
     --query_fa new_sequences.fa \
     --is_update Y

Debugging with retained BLAST temp directory:

.. code-block:: bash

   python BlastAlignment.py \
     --master_acc NC_001542 \
     --keep_blast_tmp_dir Y

Notes
^^^^^

- The **reference FASTA file** is automatically converted into a BLAST database.  
- The **master accession** ensures consistency across analyses (e.g., Rabies virus = ``NC_001542``).  
- Use ``--is_update Y`` for incremental analyses to save runtime.  
- Set ``--keep_blast_tmp_dir Y`` only when you need to inspect raw BLAST output for debugging.

Nextalign Alignment
-------------------

The script ``NextalignAlignment.py`` runs **Nextalign** for each query sequence,
aligning it to a master reference. It can optionally reuse a precomputed
reference alignment instead of generating one on the fly.

Usage
^^^^^

.. code-block:: bash

   python NextalignAlignment.py --master_ref <ACCESSION> [options]

Arguments
^^^^^^^^^

.. option:: -g, --gB_matrix <FILE>

   GenBank matrix (metadata) TSV file.  
   Default: ``tmp/GenBank-matrix/gB_matrix_raw.tsv``

.. option:: -q, --query_dir <DIR>

   Directory containing grouped query FASTA files.  
   Default: ``tmp/Blast/grouped_fasta``

.. option:: -r, --ref_dir <DIR>

   Directory containing reference FASTA files.  
   Default: ``tmp/Blast/ref_seqs``

.. option:: -f, --ref_fa_file <FILE>

   Combined reference FASTA file.  
   This file is used to perform Nextalign against the master reference sequence.  
   Default: ``tmp/Sequences/ref_seq.fa``

.. option:: -ms, --master_seq_dir <DIR>

   Directory where the master sequence (derived from the master accession) is stored.  
   Default: ``tmp/Blast/master_seq``

.. option:: -t, --tmp_dir <DIR>

   Temporary working directory.  
   Default: ``tmp``

.. option:: -m, --master_ref <ACCESSION>

   **(Required)** Master reference accession (e.g., RefSeq).  
   Example (Rabies virus): ``NC_001542``

.. option:: -n, --nextalign_dir <DIR>

   Output directory for Nextalign results.  
   Default: ``Nextalign``

.. option:: -ra, --ref_alignment_file <FILE>

   Optional **reference alignment** file to use instead of having Nextalign
   align the reference to the master on the fly.

Examples
^^^^^^^^

Run Nextalign using the default directories:

.. code-block:: bash

   python NextalignAlignment.py \
     --master_ref NC_001542

Specify a precomputed reference alignment file:

.. code-block:: bash

   python NextalignAlignment.py \
     --master_ref NC_001542 \
     --ref_alignment_file data/reference_alignment.fasta

Use custom locations for queries and references:

.. code-block:: bash

   python NextalignAlignment.py \
     --master_ref NC_001542 \
     --query_dir work/Blast/grouped_fasta \
     --ref_dir work/Blast/ref_seqs \
     --ref_fa_file work/Sequences/ref_seq.fa \
     --nextalign_dir results/Nextalign

Notes
^^^^^

- Ensure **Nextalign** is installed and available on your PATH before running the script.
- ``--ref_alignment_file`` is optional; if omitted, the script will let Nextalign
  align reference sequences against the master reference automatically.
- The **master reference accession** should correspond to the canonical reference
  genome used for your analyses (e.g., Rabies virus = ``NC_001542``).

Pad Alignment
-------------

The script ``PadAlignment.py`` inserts gaps from a **master alignment**
into corresponding subalignments (e.g., from Nextalign output).  
It ensures subalignments are padded consistently with the reference alignment
and merges them into a final alignment file.  

Usage
^^^^^

.. code-block:: bash

   python PadAlignment.py --reference_alignment <FILE> [options]

Arguments
^^^^^^^^^

.. option:: -r, --reference_alignment <FILE>

   **(Required)** Master alignment file in FASTA format.

.. option:: -i, --input_dir <DIR>

   Directory containing subalignment files (e.g., Nextalign output).  
   Default: ``tmp/Nextalign/query_aln``

.. option:: -d, --base_dir <DIR>

   Base working directory.  
   Default: ``tmp``

.. option:: -o, --output_dir <DIR>

   Directory to save padded subalignments and merged files.  
   Default: ``Pad-alignment``

.. option:: --keep_intermediate_files

   Retain intermediate padded subalignment files.  
   Default: disabled (files will be removed after processing).

Examples
^^^^^^^^

Run padding with defaults:

.. code-block:: bash

   python PadAlignment.py \
     --reference_alignment data/master_alignment.fasta

Keep intermediate padded subalignments:

.. code-block:: bash

   python PadAlignment.py \
     --reference_alignment data/master_alignment.fasta \
     --keep_intermediate_files

Use custom directories for input and output:

.. code-block:: bash

   python PadAlignment.py \
     --reference_alignment data/master_alignment.fasta \
     --input_dir work/Nextalign/query_aln \
     --output_dir results/Pad-alignment

Notes
^^^^^

- The **reference alignment** is used as a template: its gaps are inserted into
  each subalignment so all sequences share the same coordinate space.  
- By default, intermediate padded subalignment files are cleaned up automatically.  
  Use ``--keep_intermediate_files`` to retain them for debugging.  
- After padding, the script also removes redundant sequences to produce a clean
  merged alignment.

Calculate Alignment Coordinates
-------------------------------

The script ``CalcAlignmentCords.py`` calculates genome and CDS coordinates
for sequences based on a padded alignment and a master GFF3 annotation file.
It aligns features to the master reference and outputs a table of coordinates.

Usage
^^^^^

.. code-block:: bash

   python CalcAlignmentCords.py --paded_alignment <FASTA> --master_accession <ACC> --master_gff <FILE> [options]

Arguments
^^^^^^^^^

.. option:: -i, --paded_alignment <FILE(S)>

   **(Required)** Input FASTA alignment file(s).  
   Can be a single file or multiple sequence files.

.. option:: -b, --tmp_dir <DIR>

   Base working directory.  
   Default: ``tmp``

.. option:: -d, --output_dir <DIR>

   Directory to store processed data and results.  
   Default: ``Tables``

.. option:: -o, --output_file <FILE>

   Output TSV file containing calculated coordinates.  
   Default: ``features.tsv``

.. option:: -m, --master_accession <ACC>

   **(Required)** Master reference accession ID.  
   Example: ``NC_001542``

.. option:: -bh, --blast_uniq_hits <FILE>

   BLAST unique hits file used to refine coordinate mapping.  
   Default: ``tmp/Blast/query_uniq_tophits.tsv``

.. option:: -g, --master_gff <FILE>

   **(Required)** Master GFF3 annotation file for the reference sequence.

Examples
^^^^^^^^

Calculate coordinates with defaults:

.. code-block:: bash

   python CalcAlignmentCords.py \
     --paded_alignment data/padded_alignment.fa \
     --master_accession NC_001542 \
     --master_gff data/master.gff

Specify a custom output location:

.. code-block:: bash

   python CalcAlignmentCords.py \
     --paded_alignment data/aln.fa \
     --master_accession NC_001542 \
     --master_gff data/master.gff \
     --output_dir results/Tables \
     --output_file coords.tsv

Use a custom BLAST hits file:

.. code-block:: bash

   python CalcAlignmentCords.py \
     --paded_alignment data/aln.fa \
     --master_accession NC_001542 \
     --master_gff data/master.gff \
     --blast_uniq_hits work/Blast/query_uniq_tophits.tsv

Notes
^^^^^

- The **padded alignment** should be generated by previous steps (e.g., Nextalign + padding).  
- The **master GFF3** file provides the annotation used for feature coordinate mapping.  
- The **master accession** defines the coordinate reference system.  
- The script also leverages BLAST unique hits to refine CDS boundaries where applicable.

Software Version Checker
------------------------

The script ``SoftwareVersion.py`` collects information about available software
versions and generates a TSV summary table. This is useful for reproducibility
and workflow documentation.

Usage
^^^^^

.. code-block:: bash

   python SoftwareVersion.py [options]

Arguments
^^^^^^^^^

.. option:: -d, --tmp_dir <DIR>

   Temporary working directory.  
   Default: ``tmp``

.. option:: -o, --output_dir <DIR>

   Directory where the version summary table will be saved.  
   Default: ``Software_info``

.. option:: -f, --table_name <FILE>

   Name of the TSV file to store software version information.  
   Default: ``software_info.tsv``

Examples
^^^^^^^^

Generate a software version table with defaults:

.. code-block:: bash

   python SoftwareVersion.py

Save the version table with a custom name:

.. code-block:: bash

   python SoftwareVersion.py \
     --table_name versions.tsv

Write results to a custom directory:

.. code-block:: bash

   python SoftwareVersion.py \
     --output_dir results/Software_info

Notes
^^^^^

- The output is a TSV file listing the detected software and their versions.  
- Useful for documenting pipeline environments and ensuring reproducibility.  
- The script can be run at the end of a workflow to capture the environment state.

Generate Database-Ready Tables
------------------------------

The script ``GenerateTable.py`` produces a set of **TSV tables** ready for
loading into an SQLite database. It consolidates information from the GenBank
matrix, BLAST unique hits, and (padded) multiple-sequence alignments, and can
use Nextalign outputs as needed.

Usage
^^^^^

.. code-block:: bash

   python GenerateTable.py [options]

Arguments
^^^^^^^^^

.. option:: -g, --genbank_matrix <FILE>

   GenBank matrix TSV file.  
   Default: ``tmp/GenBank-matrix/gB_matrix_raw.tsv``

.. option:: -b, --base_dir <DIR>

   Base working directory.  
   Default: ``tmp``

.. option:: -o, --output_dir <DIR>

   Output directory where all DB-ready TSV files are written.  
   Default: ``Tables``

.. option:: -f, --host_taxa <FILE>

   Host taxa lookup table.  
   Default: ``host_taxa.tsv``

.. option:: -bh, --blast_hits <FILE>

   BLASTN unique hits table (top hits per query).  
   Default: ``tmp/Blast/query_uniq_tophits.tsv``

.. option:: -p, --paded_aln <FILE>

   **Padded** alignment FASTA file (reference-anchored, gap-padded).  
   Default: ``tmp/Pad-alignment/NC_001542.aligned_merged_MSA.fasta``

.. option:: -n, --nextalign_dir <DIR>

   Directory containing Nextalign alignment outputs.  
   Default: ``tmp/Nextalign/``

.. option:: -e, --email <EMAIL>

   Contact email (e.g., for NCBI-related requests/logging).  
   Default: ``your-email@example.com``

Examples
^^^^^^^^

Run with defaults (writes TSVs into ``Tables/``):

.. code-block:: bash

   python GenerateTable.py

Specify custom inputs and output location:

.. code-block:: bash

   python GenerateTable.py \
     --genbank_matrix data/gB_matrix_raw.tsv \
     --blast_hits results/Blast/query_uniq_tophits.tsv \
     --paded_aln results/Pad-alignment/aligned_merged_MSA.fasta \
     --nextalign_dir results/Nextalign \
     --output_dir results/Tables \
     --email me@lab.org

Notes
^^^^^

- Ensure the **padded alignment** corresponds to the same master reference used
  in previous steps, so feature coordinates remain consistent.
- The **BLAST unique hits** file should contain one best hit per query (or as
  expected by your pipeline).
- Output TSVs are organized in ``--output_dir`` and are intended for direct
  ingestion into your SQLite schema.

Create SQLite Database
----------------------

The script ``CreateSqliteTable.py`` creates an **SQLite database** from a
collection of input tables, metadata, alignments, and reference files.  
It consolidates processed results into a structured database ready for analysis.

Usage
^^^^^

.. code-block:: bash

   python CreateSqliteTable.py [options]

Arguments
^^^^^^^^^

.. option:: -m, --meta_data <FILE>

   Metadata table (GenBank matrix TSV).  
   Default: ``tmp/GenBank-matrix/gB_matrix_raw.tsv``

.. option:: -b, --base_dir <DIR>

   Base working directory.  
   Default: ``tmp``

.. option:: -o, --output_dir <DIR>

   Directory where the SQLite database will be stored.  
   Default: ``SqliteDB``

.. option:: -rf, --features <FILE>

   Features table TSV.  
   Default: ``tmp/Tables/features.tsv``

.. option:: -p, --pad_aln <FILE>

   Padded alignment table (TSV).  
   Default: ``tmp/Tables/sequence_alignment.tsv``

.. option:: -g, --gene_info <FILE>

   Gene information table (CSV).  
   Default: ``generic/rabv/Tables/gene_info.csv``

.. option:: -mc, --m49_countries <FILE>

   UN M49 countries CSV file.  
   Default: ``assets/m49_country.csv``

.. option:: -mir, --m49_interm_region <FILE>

   UN M49 intermediate regions CSV file.  
   Default: ``assets/m49_intermediate_region.csv``

.. option:: -mr, --m49_regions <FILE>

   UN M49 regions CSV file.  
   Default: ``assets/m49_region.csv``

.. option:: -msr, --m49_sub_regions <FILE>

   UN M49 sub-regions CSV file.  
   Default: ``assets/m49_sub_region.csv``

.. option:: -s, --proj_settings <FILE>

   Project/software settings TSV.  
   Default: ``tmp/Software_info/software_info.tsv``

.. option:: -fa, --fasta_sequences <FILE>

   FASTA sequence file used in the project.  
   Default: ``tmp/GenBank-matrix/sequences.fa``

.. option:: -i, --insertion_file <FILE>

   Nextalign insertions TSV file.  
   Default: ``tmp/Tables/insertions.tsv``

.. option:: -d, --db_name <STR>

   Name of the SQLite database (without ``.db`` extension).  
   Default: ``gdb``

Examples
^^^^^^^^

Run with default settings:

.. code-block:: bash

   python CreateSqliteTable.py

Specify a custom database name and output directory:

.. code-block:: bash

   python CreateSqliteTable.py \
     --db_name rabv_database \
     --output_dir results/SqliteDB

Use custom metadata and features tables:

.. code-block:: bash

   python CreateSqliteTable.py \
     --meta_data data/gB_matrix_raw.tsv \
     --features data/features.tsv \
     --gene_info data/gene_info.csv \
     --db_name custom_gdb

Notes
^^^^^

- The script integrates **metadata, features, alignments, gene info, country codes, and insertions** into a unified SQLite database.  
- Ensure that all required TSV/CSV inputs exist before running.  
- The output database will be stored in the directory specified by ``--output_dir`` with the name ``<db_name>.db``.

GenBank Sequence Submitter
--------------------------

The script ``GenBankSequenceSubmitter.py`` prepares sequence data and metadata
for **GenBank submission**. It uses a padded alignment, metadata, a master GFF3,
and a submission template to generate submission-ready files.

Usage
^^^^^

.. code-block:: bash

   python GenBankSequenceSubmitter.py --sequence_dir <DIR> --metadata <FILE> --ncbi_submission_template <FILE> --gff_file <FILE> --vgtk-db <DB> [options]

Arguments
^^^^^^^^^

.. option:: -q, --sequence_dir <DIR>

   **(Required)** Directory containing input FASTA sequence files  
   (single or multiple).

.. option:: -t, --tmp_dir <DIR>

   Temporary working directory.  
   Default: ``tmp``

.. option:: -o, --output_dir <DIR>

   Output directory where processed results and submission files are saved.  
   Default: ``Table2asn``

.. option:: -m, --metadata <FILE>

   **(Required)** Metadata file (tab-delimited).  
   Must include required fields for submission.

.. option:: -n, --ncbi_submission_template <FILE>

   **(Required)** NCBI submission template file (generated via:  
   `https://submit.ncbi.nlm.nih.gov/genbank/template/submission/`).

.. option:: -gff, --gff_file <FILE>

   **(Required)** Master reference GFF3 annotation file.

.. option:: -db, --vgtk-db <FILE>

   **(Required)** VGTK SQLite database file.

.. option:: -gp, --gaps_to_ignore <N>

   Length of alignment gaps to ignore during processing.  
   Default: ``30``

Examples
^^^^^^^^

Prepare sequences for submission using defaults:

.. code-block:: bash

   python GenBankSequenceSubmitter.py \
     --sequence_dir data/sequences \
     --metadata data/metadata.tsv \
     --ncbi_submission_template templates/template.sbt \
     --gff_file reference/master.gff3 \
     --vgtk-db db/vgtk.sqlite

Specify a custom gap threshold:

.. code-block:: bash

   python GenBankSequenceSubmitter.py \
     --sequence_dir data/sequences \
     --metadata data/metadata.tsv \
     --ncbi_submission_template templates/template.sbt \
     --gff_file reference/master.gff3 \
     --vgtk-db db/vgtk.sqlite \
     --gaps_to_ignore 50

Notes
^^^^^

- Ensure the **submission template** matches your organism and project type.  
- The **VGTK database** provides reference metadata for consistency.  
- The generated output in ``--output_dir`` will include files ready for upload to NCBI’s submission portal.  
- Adjust ``--gaps_to_ignore`` if small insertions/deletions should be tolerated.

