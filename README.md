# RABV-gTK (Virus Genome ToolKit)

RABV-gTK is a Nextflow-based bioinformatics pipeline designed to build structured database resources, coordinate mappings, phylogenetic trees, and mutation annotations for viral genomes.

---

## Supported Viruses

The pipeline has full out-of-the-box support for the following viruses, with predefined Nextflow profiles:

| Virus | Segmented? | Features Supported |
|---|---|---|
| **Rabies Virus (RABV)** | No | Initial build, database update, phylogenetic tree placement (UShER) |
| **Hepatitis C Virus (HCV)** | No | Database updates, mutation catalog annotation & verification (`VerifyMutations`) |
| **Influenza A Virus (IAV)** | Yes | Segment validation, strain validation, GISAID and GenBank XML merging |

---

## Environment Setup (Conda)

The pipeline requires a Conda environment containing all bioconda-native tools (Blast, Nextalign, VeryFastTree, UShER, MMseqs2, IQ-TREE, Nextflow, etc.).

### 1. Create the Environment
Initialize the environment using the provided `environment.yml` file:
```bash
conda env create -f environment.yml
```
*This creates a conda environment named `vgtk`.*

### 2. Activate the Environment
You must activate this environment before running the Nextflow pipeline:
```bash
conda activate vgtk
```

---

## Running the Pipeline

Ensure the `vgtk` conda environment is active. You can run the pipeline by specifying one of the predefined profiles listed below.

```bash
nextflow run vgtk-init.nf -profile <profile_name>
```

### Predefined Profiles by Run Type

The pipeline profiles in `nextflow.config` are grouped by their operational purposes:

| Run Type | Profile Name | Virus | Description / Inputs |
|---|---|---|---|
| **Test Runs**<br>*(Quick validation using small subsets)* | `test` | RABV | Basic run with a small set of sequences |
| | `HCV_test` | HCV | Small scale run validating mutation annotations |
| | `segmented_test` | IAV | Validation run for segmented influenza mapping |
| | `flu_gisaid_test` | IAV | Test run incorporating simulated GISAID metadata |
| **Fresh Runs**<br>*(Build full DB from scratch)* | `setup_rabv_full` | RABV | Complete database construction from scratch |
| | `HCV_full` | HCV | Complete database construction with mutation cataloging |
| **Update Runs**<br>*(Incremental updates to existing DB)* | `HCV_update` | HCV | Updates an existing HCV database with new sequences |
| | `setup_rabv_test_update` | RABV | Test run validating incremental rabies database updates |
| | `update_test` | RABV | Quick regression test checking database updates |

> [!IMPORTANT]
> **Custom XML Local Infrastructure Profiles:**
> Profiles utilizing local GenBank XML caches (such as `HCV_XML_full`, `flu_gisaid_xml_laura`, and `h10n8_xml_test`) are custom profiles configured to read from pre-downloaded directories on our local servers. These require adjustments to the `xml_dir` parameter in `nextflow.config` to run outside of our local infrastructure.

---

## Custom Configuration / Direct Runs
To run a new virus or bypass the profiles, override parameters directly on the command line:
```bash
nextflow run vgtk-init.nf \
  --tax_id <NCBI_TAX_ID> \
  --db_name <DATABASE_NAME> \
  --ref_list <PATH_TO_TSV> \
  --is_segmented <Y_OR_N> \
  --mmseqs_min_seq_id <CLUSTERING_THRESHOLD> \
  --publish_dir <OUTPUT_DIRECTORY>
```

---

## Mandatory Files to Run a New Virus

To construct a database for a new virus, you must provide the following input files:

### 1. Reference List TSV (`--ref_list`)
A tab-separated values file listing the references and masters for the virus. It can be headerless or have a header.

* **With Header (Auto-detected)**:
  Columns can use any of the following aliases:
  * **Accession ID**: `primary_accession`, `accession`, `accession_id`, `acc`, `id`
  * **Type / Status**: `accession_type`, `type`, `status`
  * **Segment**: `segment` (mandatory for segmented viruses)
  * **Genotype**: `genotype` (optional)
  * **Subtype**: `subtype` (optional)

* **Without Header**:
  Columns must be placed in this exact order:
  1. `primary_accession`
  2. `accession_type`
  3. `segment`
  4. `genotype`
  5. `subtype`

* **Segment and Accession Type Constraints**:
  * **`accession_type`** values must only be:
    * `master`: The principal reference sequence used to establish coordinate mapping space. **Exactly one** master is required per segment (unless a segment only contains exclusions).
    * `reference`: Valid representative reference sequences.
    * `exclusion_list`: Sequences to ignore/quarantine.
  * **`segment`** values must contain segment names or numbers if `--is_segmented Y` is set. Numeric digits are automatically normalized to their clean integer form (e.g., `"Segment 4"` or `"04"` maps to `"4"`). Non-numeric labels such as `L`, `M` and `S` are kept verbatim and matched case-insensitively.
  * For segmented builds the **`master`** rows also define the expected segment set for the per-isolate completeness table (see below), so adding a master raises the bar for calling an isolate complete.

### Per-isolate segment completeness

Every `--is_segmented Y` build publishes `gB_matrix_pivoted_segments.tsv`: one row
per isolate, one column per expected segment, each cell holding the accession(s)
contributing that segment, plus a `Complete_status` of `Complete`/`Incomplete`.
This used to be influenza-only; it now works for any segmented genome.

* The **expected segment set** comes from the `master` rows of the reference list
  (`exclusion_list` rows are ignored, so decoy references do not inflate the bar).
  Override with `--pivot_required_segments L,S`.
* The **isolate key** is elected once per run from the first sufficiently populated
  column of `Parsed_strain`, `isolate`, `strain`. `Parsed_strain` only exists on
  influenza runs, so non-flu builds normally group on the GenBank `isolate`
  qualifier. Override with `--pivot_isolate_key <col1,col2,...>`.
* A row whose isolate key is blank becomes its own single-accession isolate rather
  than being merged with every other unkeyed record.
* `gB_matrix_pivoted_segments.summary.tsv` records the segment set and its source,
  which key column was elected and the coverage of each candidate, and counts of
  excluded / unkeyed / segment-less / unexpected-segment rows. Check it first if
  the completeness numbers look wrong.

### 2. Gene Info TSV (`--gene_info`)
This file defines display metadata for the genes/proteins associated with the virus database. It must contain a header and the following columns:
* `name`: The canonical gene/protein identifier.
* `display_name`: Normalized/abbreviated name displayed in summary tables.
* `description`: The detailed description of the gene/protein.
* `parent_name`: Parent classification (typically `whole_genome`).

### 3. Aligned Reference Directory (`--ref_set_aligned`)
**Mandatory when `--is_segmented Y` is enabled.**
This directory must contain precomputed multiple sequence alignment (MSA) FASTA files of the reference sequences for each segment to guide padding and alignment.
* **Naming Convention**:
  Files must be named **`refset_<segment>_aln.fasta`** where `<segment>` is the normalized segment value (e.g. `refset_1_aln.fasta`, `refset_2_aln.fasta`).
  *Fallback:* If the specific pattern is not found, the pipeline scans the directory for a FASTA file whose name contains the segment number digits.
