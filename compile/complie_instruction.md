# V-gTK Conda Build

This directory contains the Conda recipe used to build V-gTK locally.

## Requirements

Create a Conda build environment:

```bash
conda create -n vgtk-build -c conda-forge conda-build
```

Activate it:

```bash
conda activate vgtk-build
```

## Build the package

From the root of the V-gTK repository:

```bash
conda build recipe \
    -c conda-forge \
    -c bioconda
```

## Get the package location

After a successful build:

```bash
conda build recipe \
    -c conda-forge \
    -c bioconda \
    --output
```

The generated package should be located under a directory similar to:

```text
~/miniconda3/envs/vgtk-build/conda-bld/linux-64/
```

For example:

```text
vgtk-1.0.0-1.conda
```

## Test the local package

Create a new test environment:

```bash
conda create -n vgtk-test \
    -c local \
    -c bioconda \
    -c conda-forge \
    vgtk
```

Activate it:

```bash
conda activate vgtk-test
```

Check that V-gTK is installed:

```bash
vgtk
```

Example:

```bash
vgtk DownloadGFF --help
```

## Clean previous builds

If required, remove previous Conda build files:

```bash
conda build purge
```

Then rebuild:

```bash
conda build recipe \
    -c conda-forge \
    -c bioconda
```

## Supported platforms

The Conda recipe is intended to support:

```text
linux-64
osx-arm64
```

`linux-64` is for Linux systems running on Intel/AMD x86-64 processors.

`osx-arm64` is for Apple Silicon Macs such as M1, M2, M3 and M4.

