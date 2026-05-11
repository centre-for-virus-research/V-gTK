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

V-gTK includes a real-time rabies virus example workflow for building and using a local RABV genomic database. This workflow is organised inside the `rabv-gdb_build` directory and provides the files and scripts required to construct a database from curated RABV reference sequences and their associated metadata.

The RABV database acts as the core reference resource for downstream V-gTK analyses. Once created, it can be used to compare newly generated or custom RABV sequences against existing reference data, prepare combined datasets, generate alignments, and support phylogenetic or clade-assignment workflows.

### Required input files

Before running the database build workflow, the required input files should be available inside `rabv-gdb_build`. These may include:

- curated RABV reference sequences in FASTA format
- sequence metadata files
- accession or GenBank record information
- clade, lineage, or designation files
- configuration files used by the database-building scripts

These files are used to link each sequence with its corresponding metadata and classification information. The database-building scripts then organise the input data into a structured local database that can be used by the rest of the V-gTK framework.

### Building the RABV database

To build the RABV database, navigate to the `rabv-gdb_build` directory and run the provided database construction scripts.

```shell
cd rabv-gdb_build
