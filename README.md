# V-gTK
V-gTK is a bioinformatics framework designed to support viral genomic analysis through a structured, reproducible, and modular workflow. The framework provides tools for managing viral genome datasets, preparing sequence alignments, integrating custom sequences with reference databases, and supporting downstream phylogenetic and clade-assignment analyses. V-gTK is developed with a focus on supporting both segnemted and non-segmented virus. By combining database-driven sequence handling, automated alignment preparation, metadata integration, and tree-based analysis support, V-gTK aims to simplify complex genomic workflows and improve consistency across analyses.

## Installation
Install Conda using one of the following options: Miniconda, Miniforge, or Anaconda.
After installation, create the V-gTK environment using the provided environment.yml file:

```shell
conda env create --file environment.yml
```

Once the environment is created, activate it using:

```shell
conda activate vgtk
```

## Working with a Real-Time RABV Example

V-gTK includes a real-time rabies virus example workflow for building and using a local RABV genomic database. This workflow is organised inside the `rabv-gdb_build` repo and provides the files required to construct a database from curated RABV reference sequences.

Once the RABV database has been created, the same database can be maintained, updated, curated, and modified as required. V-gTK provides a set of easy-to-use scripts that allow users to manage the database without needing to rebuild the entire workflow from the beginning each time.

These scripts support routine database maintenance tasks such as adding new RABV sequences, updating existing metadata, correcting sequence records, modifying clade or lineage annotations, and removing or replacing outdated entries. This makes the database suitable for real-time genomic surveillance, where new sequences and updated information may need to be incorporated regularly.

The framework is designed to keep the database flexible and reusable. Users can curate the reference dataset as new information becomes available, apply manual corrections where needed, and regenerate analysis-ready files for downstream workflows. By providing dedicated scripts for database updates and curation, V-gTK helps ensure that the RABV genomic database remains consistent, reproducible, and up to date.

In this way, V-gTK separates the initial database-building step from the ongoing database management process. After the first build, users can continue to refine and extend the same database through simple scripted operations, making it easier to support long-term virus specific genomic analysis and surveillance.

### Required input files

Before running the database build workflow, the required input files should be available inside `rabv-gdb_build`. These may include:

- curated RABV reference accessions
- reference tree with meta data
- clade, lineage, or designation files
- configuration files used by the database-building scripts

These files are used to link each sequence with its corresponding metadata and classification information. The database-building scripts then organise the input data into a structured local database that can be used by the rest of the V-gTK framework.

### Building the RABV database

To build the RABV database, navigate to the `rabv-gdb_build` directory and run the provided database construction scripts.

```shell
cd rabv-gdb_build
