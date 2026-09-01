# Pipeline Testing

This document describes the automated tests for the VGTK pipeline.

## Test Processes

The pipeline includes automated validation tests that run when `--test 1` is set in the parameters.

### TEST_SEGMENTED_OUTPUT

**When**: Runs for segmented viruses (e.g., influenza) when `is_segmented="Y"` and `test="1"`

**Validates**:
1. **Annotated BLAST file structure** - Verifies that `query_uniq_tophit_annotated.tsv` has 5 columns:
   - query accession
   - reference accession
   - alignment score
   - strand orientation
   - **segment** (the key column for segmented viruses; a number for influenza,
     but any label the reference list uses, e.g. `L`/`M`/`S`)

2. **Segment validation in matrix** - Checks that:
   - `segment_validated` column exists in the GenBank matrix
   - At least some records have valid segment assignments
   - Reports count of validated vs. total records

3. **Pivoted segments matrix** - runs for every `is_segmented = "Y"` build, not just
   influenza. Verifies that:
   - Pivoted matrix has `Complete_status` column
   - Segment columns are present. These are the expected segment set derived from
     the `master` rows of the reference list (1-8 for influenza, `L`/`S` for an
     arenavirus, and so on), so the count is virus-dependent.
   - Reports counts of Complete vs. Incomplete genomes

   The first column is the elected isolate key: `Parsed_strain` on flu runs,
   otherwise the GenBank `isolate` or `strain` qualifier. Which column was elected,
   and how well populated each candidate was, is recorded in
   `gB_matrix_pivoted_segments.summary.tsv` alongside the table.

**Output**: `test_segmented_results.txt` in the publish directory

### TEST_NON_SEGMENTED_OUTPUT

**When**: Runs for non-segmented viruses (e.g., RABV) when `is_segmented="N"` and `test="1"`

**Validates**:
1. **BLAST file structure** - Verifies that `query_uniq_tophits.tsv` has 4 columns only:
   - query accession
   - reference accession
   - alignment score
   - strand orientation
   - (NO segment column - this is correct for non-segmented)

2. **No segment column** - Confirms that `segment_validated` column is NOT present in the matrix (as expected)

3. **Sequence processing** - Verifies that:
   - BLAST hits were generated
   - Matrix records were created
   - Reports counts for both

**Output**: `test_non_segmented_results.txt` in the publish directory

## What test mode changes

`--test 1` does not just run extra checks. It substitutes the cheapest available
version of several steps, so **a test run's outputs are not comparable to a
production build**. Every substitution is gated on `params.test` and is off by
default.

| step | production | test mode | measured effect |
|---|---|---|---|
| clustering input | everything | all references + up to `test_max_cluster_seqs` queries | see below |
| MMseqs clustering | `mmseqs cluster` (sensitive) | `mmseqs linclust` (linear time, less sensitive) | 29.3s → 3.4s on the HCV input; identical 276 representatives there, but that is not guaranteed |
| IQ-TREE | full ML search | `--fast` ("search to resemble FastTree") | ~13 min → 14s on 277 HCV representatives, at a slightly worse likelihood |
| DB validation | strict | `--test-mode` relaxations | segment-tree coverage and strict consistency failures are downgraded |

**`test_max_cluster_seqs` caps QUERIES, not the total.** Every reference and the
master is kept unconditionally, because references seed the clusters, define the
tree topology, anchor UShER placement, and are what genotype calls are made
against. A plain `seqkit sample` over the whole alignment does not know that: it
used to drop 153 of 237 references on the HCV profile, producing a database whose
UShER tree held 85 of 238 references. `TEST_SUBSAMPLE_CLUSTER_INPUT` now fails
outright if any reference is missing from its own output rather than warning,
because nobody reads Nextflow warnings.

Consequences worth knowing:

- On most profiles the reference set is **larger** than the cap
  (`segmented_test` 418 vs 40, `HCV_test` 238 vs 120). That is expected, not a
  problem - the cap only ever bounded the query half.
- Keeping all references made the HCV test cluster into 277 representatives
  instead of 106, because HCV references span genotypes 1-8 at roughly 30%
  divergence so nearly every one is its own cluster. `--fast` and `linclust`
  are what keep the run quick despite that.
- `--test-mode` is passed to `ValidateDbTree` **only** when `params.test` is 1.
  It used to be hardcoded at both call sites, which made production validation
  weaker than test validation.

Never read a topology, a clustering, or a likelihood off a test run.

### Avoiding the NCBI fetch: `segmented_xml_test`

`segmented_xml_test` is a drop-in replacement for `segmented_test` that reads a
frozen GenBank XML snapshot instead of calling NCBI. Same taxid, same reference
list, same `publish_dir`, so a CI assertion on `test_out/IAV_DB/...` needs only
the profile name changed.

```bash
nextflow run vgtk-init.nf -profile conda,segmented_xml_test
```

Measured on a clean run (no `-resume`), total task time:

| | with fetch | with fixture |
|---|---|---|
| `FETCH_GENBANK` | 129.0s (37.8%) | **not run** |
| everything else | 212s | 212s |
| wall clock | 3m54s | **1m48s** |

Two reasons to prefer it in CI, and the second matters more:

1. The fetch is the single largest cost of the run.
2. Its output depends on what NCBI returns that day. A difference there has
   already produced a CI failure that could not be reproduced locally - a fresh
   local run produced a database identical to CI's on every statistic (518
   meta_data rows, 9 UShER trees, 406 centroids, 22 `exclusion_list` entries)
   yet differed in whether those 22 carried `exclusion_status=1`.

**A fixture hides that class of bug rather than fixing it.** `_classify_accession`
in `FilterAndExtractSequences.py` falls back to a reference-list lookup keyed on
`row['gi_number']`; if a fetch ever returns versioned accessions (`NC_002204.1`)
where the reference list holds bare ones, the lookup misses, the type defaults to
`query`, the exclusion flag is never set, and influenza B references are then
required to have alignments against influenza A data. That is the leading theory
for the CI failure and it is still unconfirmed. Freezing the input stops the test
flapping; it does not make the lookup robust.

### Refreshing the snapshot

```bash
scripts/fetch_segmented_fixture.sh --email you@example.com --clean
```

Writes `test_data/iav_11320/` (5.8 MB, 52 XML files, 518 records) with a
`manifest.txt` and `checksums.sha256`. Refresh deliberately - the point of a
fixture is that it does not move under you. `segmented_test` still exists and
still exercises the live fetch; run it when you want to know NCBI has not
changed under the pipeline.

The same pattern already exists for H10N8 - see the fixture section below.

## Running Tests

## Integration fixture test: H10N8 (taxid 286285)

This integration flow uses a frozen XML snapshot under `test_data/h10n8_286285/GenBank-XML` so upstream NCBI updates do not make tests non-reproducible.

### 1) Create or refresh fixture XML snapshot

```bash
scripts/fetch_h10n8_fixture.sh --email your_email@example.com --clean
```

### 2) Run end-to-end Nextflow integration checks

```bash
scripts/run_h10n8_xml_test.sh --email your_email@example.com
```

What this runner validates:
- pipeline completes using profile `h10n8_xml_test`
- SQLite DB is created and required tables are non-empty
- `ValidateDbTree.py` passes
- optional Robinson-Foulds regression check if `test_data/h10n8_286285/expected/iqtree.treefile` exists

To enable RF regression checks, place a baseline tree at:

```text
test_data/h10n8_286285/expected/iqtree.treefile
```

Optional RF threshold override:

```bash
MAX_NORMALIZED_RF=0.20 scripts/run_h10n8_xml_test.sh --email your_email@example.com
```

## Python Unit Tests (new)

Pytest-based unit tests now exist for the first two scripts in [scripts](scripts):
- [scripts/AddMissingData.py](scripts/AddMissingData.py)
- [scripts/BlastAlignment.py](scripts/BlastAlignment.py)

Additional unit tests now also cover:
- [scripts/CalcAlignmentCord.py](scripts/CalcAlignmentCord.py)
- [scripts/CalcGenomeCords.py](scripts/CalcGenomeCords.py)

Test files:
- [tests/unit/test_add_missing_data.py](tests/unit/test_add_missing_data.py)
- [tests/unit/test_blast_alignment.py](tests/unit/test_blast_alignment.py)
- [tests/unit/test_calc_alignment_cord.py](tests/unit/test_calc_alignment_cord.py)
- [tests/unit/test_calc_genome_cords.py](tests/unit/test_calc_genome_cords.py)

Deterministic input/expected fixtures:
- [test_data/unit/add_missing_data](test_data/unit/add_missing_data)
- [test_data/unit/blast_alignment](test_data/unit/blast_alignment)
- [test_data/unit/calc_alignment_cord](test_data/unit/calc_alignment_cord)
- [test_data/unit/calc_genome_cords](test_data/unit/calc_genome_cords)

Run locally:

```bash
pytest -q tests/unit
```

Coverage is enabled by default for pytest runs via [pytest.ini](pytest.ini) and [.coveragerc](.coveragerc).
Each run prints per-script coverage for [scripts](scripts) in the terminal and writes:
- [coverage.xml](coverage.xml) (Cobertura XML)
- [htmlcov/index.html](htmlcov/index.html) (interactive HTML report)

Run unit tests with coverage explicitly:

```bash
pytest -q tests/unit
```

Optional: enforce minimum script coverage (example 80%):

```bash
pytest -q tests/unit --cov-fail-under=80
```

CI execution on every push/PR:
- [.github/workflows/python-unit-tests.yml](.github/workflows/python-unit-tests.yml)

### For Segmented Virus (Influenza)
```bash
nextflow run vgtk-init.nf -profile segmented_test
```

Make sure your config includes:
```groovy
params.is_segmented = "Y"
params.is_flu = "Y"
params.test = "1"
```

### For Non-Segmented Virus (RABV)
```bash
nextflow run vgtk-init.nf -profile test
```

Make sure your config includes:
```groovy
params.is_segmented = "N"
params.test = "1"
```

## Test Output

Test results are published to `${params.publish_dir}/tests/` and include:
- ✓ PASS: Test passed successfully
- ✗ FAIL: Test failed (pipeline will exit with error)
- ⚠ WARNING: Potential issue detected (pipeline continues)

## Interpreting Results

### Segmented Virus Tests

**Expected output for successful flu run**:
```
=== Testing Segmented Virus Pipeline Output ===

Test 1: Checking annotated BLAST file structure...
✓ PASS: Annotated BLAST file has 5 columns (query, reference, score, strand, segment)

Test 2: Checking segment_validated column in matrix...
✓ PASS: segment_validated column exists
  - Found 150 records with valid segments out of 200 total
✓ PASS: At least some records have segment assignments

Test 3: Checking pivoted segments matrix...
✓ PASS: Pivoted matrix has Complete_status column
  - Found 8 segment columns          # virus-dependent: the ref_list master segments
  - Complete genomes: 15
  - Incomplete genomes: 10
✓ PASS: Pivoted matrix contains strain data

=== All segmented virus tests completed ===
```

### Non-Segmented Virus Tests

**Expected output for successful RABV run**:
```
=== Testing Non-Segmented Virus Pipeline Output ===

Test 1: Checking BLAST file structure...
✓ PASS: BLAST file has 4 columns (query, reference, score, strand)

Test 2: Verifying no segment column for non-segmented virus...
✓ PASS: No segment_validated column (correct for non-segmented virus)

Test 3: Checking sequence counts...
  - BLAST hits: 95
  - Matrix records: 95
✓ PASS: Pipeline processed sequences

=== All non-segmented virus tests completed ===
```

## Troubleshooting

### Common Issues

1. **"Annotated BLAST file has 4 columns, expected 5"**
   - The segmented virus pipeline isn't creating the annotated file
   - Check that `is_segmented="Y"` in params
   - Verify BlastAlignment.py has segment_file parameter

2. **"No records have valid segment assignments"**
   - BLAST couldn't match sequences to reference segments
   - Check reference list has segment information
   - Verify BLAST e-value threshold isn't too stringent

3. **"Pivoted matrix is empty"**
   - No sequences passed segment validation
   - Check exclusion criteria in ValidateSegment.py
   - Review BLAST alignment quality

4. **"segment_validated column found (unexpected for non-segmented)"**
   - Pipeline ran segmented steps for non-segmented virus
   - Verify `is_segmented="N"` in config
   - Check workflow conditional logic

## Adding New Tests

To add additional validation tests, create a new process following this template:

```groovy
process TEST_NEW_VALIDATION{
    publishDir "${params.publish_dir}/tests"
    when:
        // Your conditions here
    input:
        // Your inputs
    output:
        path "test_results.txt"
    shell:
    '''
    #!/bin/bash
    set -e
    
    echo "Test description" > test_results.txt
    
    # Your validation logic
    if [ condition ]; then
        echo "✓ PASS: Test passed" >> test_results.txt
    else
        echo "✗ FAIL: Test failed" >> test_results.txt
        exit 1
    fi
    
    cat test_results.txt
    '''
}
```

Then call it in the workflow section with appropriate inputs.
