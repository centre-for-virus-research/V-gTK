#!/bin/bash
# =============================================================================
# vgtk-rabv.sh  –  bash equivalent of the Nextflow vgtk-init.nf pipeline
#
# Default settings mirror the `test` Nextflow profile:
#   tax_id              = 11292  (RABV)
#   db_name             = rabv-jul0425
#   is_segmented        = N
#   test                = 1  (test-run mode)
#   max_threads         = 8
#   mmseqs_min_seq_id   = 0.98
#   test_max_cluster_seqs = 250
#   publish_dir         = test_out/basic_test
#   ref_list            = test_data/rabv_test_ref_list.txt
#
# Usage:
#   ./vgtk-rabv.sh [options]
#
# Options (all optional – defaults match the `test` profile):
#   --tax_id            <id>    NCBI taxonomy ID            (default: 11292)
#   --db_name           <name>  SQLite DB name              (default: rabv-jul0425)
#   --ref_list          <path>  Reference list file         (default: test profile)
#   --publish_dir       <path>  Output directory            (default: test_out/basic_test)
#   --is_segmented      Y|N     Segmented virus flag        (default: N)
#   --test              0|1     Test-run mode               (default: 1)
#   --update_db         <path>  Existing .db to update      (default: none)
#   --xml_dir           <path>  Pre-fetched GenBank XML dir (default: fetch from NCBI)
#   --tree_free                 Skip tree building entirely
#   --max_threads       <n>     CPU threads                 (default: 8)
#   --mmseqs_min_seq_id <f>     MMseqs2 min seq identity    (default: 0.98)
#   --test_max_cluster_seqs <n> Cap seqs for clustering     (default: 250)
#   --email             <addr>  Contact email for NCBI      (default: your_email@example.com)
#   --mutation_catalog  <path>  Mutation catalog TSV        (default: none)
#   --mutation_virus    <name>  Virus name for catalog      (default: none)
#   --max_aln_gap_proportion <f> Max gap fraction in aln   (default: 0.96)
#   --min_seq_length_ratio   <f> Min length fraction       (default: 0.05)
#   --start_step             <n> Resume from step N         (default: 1)
#                                Steps: 1=FETCH_GENBANK 2=DOWNLOAD_GFF 3=GENBANK_PARSER
#                                       4=FILTER_AND_EXTRACT 5=BLAST 6=NEXTALIGN
#                                       7=COLLECT_FILTERED 8=PAD_ALIGNMENT 9=DEDUP
#                                       10=MMSEQS/TREES 11=TREE_MANIFEST 12=CALC_ALN_CORD
#                                       13=SOFTWARE_VERSION 14=HOST_TAXA 15=GENERATE_TABLES
#                                       16=CREATE_DB 17=VALIDATE_DB 18=TEST_OUTPUT 19=VERIFY
# =============================================================================

set -euo pipefail

# ---------------------------------------------------------------------------
# Resolve script / project root
# ---------------------------------------------------------------------------
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SCRIPTS="${SCRIPT_DIR}/scripts"
ASSETS="${SCRIPT_DIR}/assets"

# ---------------------------------------------------------------------------
# Default parameters  (= `test` Nextflow profile)
# ---------------------------------------------------------------------------
TAX_ID="11292"
DB_NAME="rabv-jul0425"
REF_LIST="${SCRIPT_DIR}/test_data/rabv_test_ref_list.txt"
PUBLISH_DIR="${SCRIPT_DIR}/test_out/basic_test"
IS_SEGMENTED="N"
TEST_MODE="1"
UPDATE_DB=""
XML_DIR=""
TREE_FREE="false"
MAX_THREADS=8
MMSEQS_MIN_SEQ_ID="0.98"
TEST_MAX_CLUSTER_SEQS="250"
EMAIL="your_email@example.com"
MUTATION_CATALOG=""
MUTATION_VIRUS=""
MAX_ALN_GAP_PROPORTION="0.96"
MIN_SEQ_LENGTH_RATIO="0.05"

# ---------------------------------------------------------------------------
# Parse CLI arguments
# ---------------------------------------------------------------------------
while [[ $# -gt 0 ]]; do
    case "$1" in
        --tax_id)                 [[ $# -lt 2 ]] && { echo "[error] --tax_id requires a value" >&2; exit 1; }; TAX_ID="$2";               shift 2 ;;
        --db_name)                [[ $# -lt 2 ]] && { echo "[error] --db_name requires a value" >&2; exit 1; }; DB_NAME="$2";              shift 2 ;;
        --ref_list)               [[ $# -lt 2 ]] && { echo "[error] --ref_list requires a value" >&2; exit 1; }; REF_LIST="$2";             shift 2 ;;
        --publish_dir)            [[ $# -lt 2 ]] && { echo "[error] --publish_dir requires a value" >&2; exit 1; }; PUBLISH_DIR="$2";          shift 2 ;;
        --is_segmented)           [[ $# -lt 2 ]] && { echo "[error] --is_segmented requires a value" >&2; exit 1; }; IS_SEGMENTED="$2";         shift 2 ;;
        --test)                   [[ $# -lt 2 ]] && { echo "[error] --test requires a value" >&2; exit 1; }; TEST_MODE="$2";            shift 2 ;;
        --update_db)              [[ $# -lt 2 ]] && { echo "[error] --update_db requires a value" >&2; exit 1; }; UPDATE_DB="$2";            shift 2 ;;
        --xml_dir)                [[ $# -lt 2 ]] && { echo "[error] --xml_dir requires a value" >&2; exit 1; }; XML_DIR="$2";              shift 2 ;;
        --tree_free)              TREE_FREE="true";          shift   ;;
        --max_threads)            [[ $# -lt 2 ]] && { echo "[error] --max_threads requires a value" >&2; exit 1; }; MAX_THREADS="$2";          shift 2 ;;
        --mmseqs_min_seq_id)      [[ $# -lt 2 ]] && { echo "[error] --mmseqs_min_seq_id requires a value" >&2; exit 1; }; MMSEQS_MIN_SEQ_ID="$2";   shift 2 ;;
        --test_max_cluster_seqs)  [[ $# -lt 2 ]] && { echo "[error] --test_max_cluster_seqs requires a value" >&2; exit 1; }; TEST_MAX_CLUSTER_SEQS="$2"; shift 2 ;;
        --email)                  [[ $# -lt 2 ]] && { echo "[error] --email requires a value" >&2; exit 1; }; EMAIL="$2";                shift 2 ;;
        --mutation_catalog)       [[ $# -lt 2 ]] && { echo "[error] --mutation_catalog requires a value" >&2; exit 1; }; MUTATION_CATALOG="$2";     shift 2 ;;
        --mutation_virus)         [[ $# -lt 2 ]] && { echo "[error] --mutation_virus requires a value" >&2; exit 1; }; MUTATION_VIRUS="$2";       shift 2 ;;
        --max_aln_gap_proportion) [[ $# -lt 2 ]] && { echo "[error] --max_aln_gap_proportion requires a value" >&2; exit 1; }; MAX_ALN_GAP_PROPORTION="$2"; shift 2 ;;
        --min_seq_length_ratio)   [[ $# -lt 2 ]] && { echo "[error] --min_seq_length_ratio requires a value" >&2; exit 1; }; MIN_SEQ_LENGTH_RATIO="$2";  shift 2 ;;
        --start_step)             [[ $# -lt 2 ]] && { echo "[error] --start_step requires a value" >&2; exit 1; }; START_STEP="$2";             shift 2 ;;
        -h|--help)
            sed -n '3,/^# ===/p' "${BASH_SOURCE[0]}"
            exit 0
            ;;
        *) echo "[error] Unknown option: $1" >&2; exit 1 ;;
    esac
done

# ---------------------------------------------------------------------------
# Derived flags (mirrors Nextflow workflow variable logic)
# ---------------------------------------------------------------------------
UPDATE_MODE="false"
[[ -n "$UPDATE_DB" ]] && UPDATE_MODE="true"

if [[ "$UPDATE_MODE" == "true" ]]; then
    [[ -f "$UPDATE_DB" ]] || { echo "[error] --update_db must be an existing file: $UPDATE_DB" >&2; exit 1; }
    [[ "$UPDATE_DB" == *.db ]] || { echo "[error] --update_db must end with .db: $UPDATE_DB" >&2; exit 1; }
fi

if [[ -n "$XML_DIR" && ! -d "$XML_DIR" ]]; then
    echo "[error] --xml_dir must be an existing directory: $XML_DIR" >&2; exit 1
fi

[[ "$IS_SEGMENTED" == "Y" || "$IS_SEGMENTED" == "N" ]] || {
    echo "[error] --is_segmented must be Y or N" >&2; exit 1
}

# ---------------------------------------------------------------------------
# Helper: print a banner before each step
# ---------------------------------------------------------------------------
run_step() {
    local name="$1"; shift
    echo ""
    echo "================================================================================"
    echo "[step] ${name}"
    echo "================================================================================"
    "$@"
    echo "[done] ${name}"
}

# ---------------------------------------------------------------------------
# Dependency check – fail fast before wasting compute
# ---------------------------------------------------------------------------
_check_dep() { command -v "$1" >/dev/null 2>&1 || { echo "[error] Required tool not found in PATH: $1" >&2; return 1; }; }
_dep_ok=true
_check_dep python    || _dep_ok=false
_check_dep seqkit    || _dep_ok=false
_check_dep mmseqs    || _dep_ok=false
if ! command -v iqtree3 >/dev/null 2>&1; then
    echo "[error] Required tool not found in PATH: iqtree3" >&2; _dep_ok=false
fi
if ! command -v VeryFastTree >/dev/null 2>&1 && ! command -v FastTree >/dev/null 2>&1; then
    echo "[error] Neither VeryFastTree nor FastTree found in PATH" >&2; _dep_ok=false
fi
[[ "$_dep_ok" == true ]] || exit 1
unset _dep_ok _check_dep

# ---------------------------------------------------------------------------
# Working directory setup + START_STEP option
# ---------------------------------------------------------------------------
START_STEP="${START_STEP:-1}"

WORK_DIR="${PUBLISH_DIR}/work_sh"
mkdir -p "$PUBLISH_DIR" "$WORK_DIR"
cd "$WORK_DIR"

# ---------------------------------------------------------------------------
# PRE-DECLARE all intermediate path variables (static locations)
# These are always defined so that any --start_step value works without
# needing to comment out earlier blocks.
# ---------------------------------------------------------------------------
GENBANK_XML_DIR="${XML_DIR:-${WORK_DIR}/GenBank-XML}"
GFF_FILE=""          # discovered below
GB_MATRIX="${WORK_DIR}/gB_matrix_validated.tsv"
SEQUENCES_FA="${WORK_DIR}/sequences.fa"
QUERY_FA="${WORK_DIR}/query_seq.fa"
REF_FA="${WORK_DIR}/ref_seq.fa"
BLAST_TOPHITS="${WORK_DIR}/query_uniq_tophits.tsv"
GROUPED_FASTA_DIR="${WORK_DIR}/grouped_fasta"
REF_SEQS_DIR="${WORK_DIR}/ref_seqs"
REF_SEQ_FILTERED_FA="${WORK_DIR}/ref_seq_filtered.fa"
MASTER_SEQ_DIR="${WORK_DIR}/master_seq"
NEXTALIGN_DIR="${WORK_DIR}/Nextalign"
FILTERED_IDS="${WORK_DIR}/filtered_sequences_ids.txt"
FILTERED_TSV="${WORK_DIR}/filtered_sequences.tsv"
PADDED_MSA_FILES=()  # populated by step 8 or rediscovered below
DEDUP_FILES=()       # populated by step 9 or rediscovered below
MMSEQ_DIRS=()        # populated by step 10 or rediscovered below
IQTREE_DIRS=()       # populated by step 10 or rediscovered below
USHER_DIRS=()        # populated by step 10 or rediscovered below
TREE_MANIFEST="${WORK_DIR}/tree_manifest.tsv"
IQTREE_INPUT_DIR="${WORK_DIR}/iqtree_inputs"
USHER_INPUT_DIR="${WORK_DIR}/usher_inputs"
PADDED_ALN_STAGING="${WORK_DIR}/padded_alignments"
FEATURES_TSV="${WORK_DIR}/features.tsv"
SOFTWARE_INFO="${WORK_DIR}/Software_info/software_info.tsv"
HOST_TAXA="${WORK_DIR}/HostTaxa/Host_taxa.tsv"
SEQ_ALN="${WORK_DIR}/Tables/sequence_alignment.tsv"
INSERTIONS="${WORK_DIR}/Tables/insertions.tsv"
SQLITE_DB="${WORK_DIR}/${DB_NAME}.db"

# ---------------------------------------------------------------------------
# REDISCOVER DYNAMIC INTERMEDIATES when skipping early steps
# Scans existing work_sh outputs to populate variables that are normally set
# as side-effects of running the step loops.
# ---------------------------------------------------------------------------
if (( START_STEP > 2 )); then
    GFF_FILE="$(find "${WORK_DIR}" -maxdepth 1 -iname '*.gff*' | head -n1)"
    [[ -z "$GFF_FILE" ]] && echo "[warn] No .gff3 found in ${WORK_DIR} – step 14 will fail if reached" >&2
fi

if (( START_STEP > 8 )); then
    mapfile -t PADDED_MSA_FILES < <(find "${WORK_DIR}" -maxdepth 1 -iname '*_merged_msa.fasta' | sort)
    [[ ${#PADDED_MSA_FILES[@]} -eq 0 ]] && echo "[warn] No *_merged_MSA.fasta found – step 9 will fail if reached" >&2
fi

if (( START_STEP > 9 )); then
    mapfile -t DEDUP_FILES < <(
        find "${WORK_DIR}" -maxdepth 1 \( -name '*_dedup_cluster_input.fasta' -o -name '*_dedup.fasta' \) \
        | sort | awk '!seen[gensub(/_dedup(_.+)?\.fasta$/,"","1")]++' )
    [[ ${#DEDUP_FILES[@]} -eq 0 ]] && echo "[warn] No dedup FASTA found – steps 10-12 will fail if reached" >&2
fi

if (( START_STEP > 12 )); then
    mapfile -t MMSEQ_DIRS  < <(find "${WORK_DIR}" -maxdepth 1 -type d -name 'MMseqClusters_*'      | sort)
    mapfile -t IQTREE_DIRS < <(find "${WORK_DIR}" -maxdepth 1 -type d -name 'IQTree_MMseqClusters_*'  | sort)
    mapfile -t USHER_DIRS  < <(find "${WORK_DIR}" -maxdepth 1 -type d -name 'Usher*'                   | sort)
    mkdir -p "${IQTREE_INPUT_DIR}" "${USHER_INPUT_DIR}"
    for _d in "${IQTREE_DIRS[@]+"${IQTREE_DIRS[@]}"}"; do
        ln -sfn "$_d" "${IQTREE_INPUT_DIR}/$(basename "$_d")" 2>/dev/null || true
    done
    for _d in "${USHER_DIRS[@]+"${USHER_DIRS[@]}"}"; do
        ln -sfn "$_d" "${USHER_INPUT_DIR}/$(basename "$_d")" 2>/dev/null || true
    done
fi

if (( START_STEP > 13 )); then
    mapfile -t PADDED_MSA_FILES < <(find "${WORK_DIR}" -maxdepth 1 -iname '*_merged_msa.fasta' | sort)
fi

echo "[info] Starting from step ${START_STEP}"
STEP=0  # incremented at the start of each major step block

# ---------------------------------------------------------------------------
# STEP 1 – FETCH_GENBANK  (skipped when --xml_dir is provided)
# ---------------------------------------------------------------------------
STEP=$(( STEP + 1 ))
if (( STEP >= START_STEP )); then
    if [[ -n "$XML_DIR" ]]; then
        echo "[info] Using pre-fetched GenBank XML: ${XML_DIR}"
    else
        FETCH_EXTRA=()
        [[ "$TEST_MODE" == "1" ]] && FETCH_EXTRA+=( --test_run --ref_list "${REF_LIST}" )
        [[ "$UPDATE_MODE" == "true" ]] && FETCH_EXTRA+=( --update "${UPDATE_DB}" )
        run_step "FETCH_GENBANK" \
            python "${SCRIPTS}/GenBankFetcher.py" \
                --taxid "${TAX_ID}" -b 50 \
                -e "${EMAIL}" \
                -o . \
                "${FETCH_EXTRA[@]+${FETCH_EXTRA[@]}}"
        GENBANK_XML_DIR="${WORK_DIR}/GenBank-XML"
    fi
fi

# ---------------------------------------------------------------------------
# STEP 2 – DOWNLOAD_GFF
# ---------------------------------------------------------------------------
STEP=$(( STEP + 1 ))
if (( STEP >= START_STEP )); then
    DOWNLOAD_GFF_EXTRA=()
    [[ "$UPDATE_MODE" == "true" ]] && DOWNLOAD_GFF_EXTRA+=( --update_db "${UPDATE_DB}" )
    run_step "DOWNLOAD_GFF" \
        python "${SCRIPTS}/DownloadGFF.py" \
            --accession_ids "${REF_LIST}" \
            -o . -b . \
            "${DOWNLOAD_GFF_EXTRA[@]+${DOWNLOAD_GFF_EXTRA[@]}}"
fi
# Always (re)discover after step 2 in case it just ran or we skipped to here
GFF_FILE="$(find "${WORK_DIR}" -maxdepth 1 -iname '*.gff*' | head -n1)"

# ---------------------------------------------------------------------------
# STEP 3 – GENBANK_PARSER + VALIDATE_MATRIX
# ---------------------------------------------------------------------------
STEP=$(( STEP + 1 ))
if (( STEP >= START_STEP )); then
    PARSER_EXTRA=( --require_refs )
    if [[ -n "$XML_DIR" && "$TEST_MODE" == "1" ]]; then
        PARSER_EXTRA+=( --test_run )
    fi
    [[ "$UPDATE_MODE" == "true" ]] && PARSER_EXTRA+=( --update "${UPDATE_DB}" )

    run_step "GENBANK_PARSER" \
        python "${SCRIPTS}/GenBankParser.py" \
            -r "${REF_LIST}" \
            -d "${GENBANK_XML_DIR}" \
            -o . -b . \
            -s "${IS_SEGMENTED}" \
            --min_length_ratio "${MIN_SEQ_LENGTH_RATIO}" \
            "${PARSER_EXTRA[@]}"

    run_step "VALIDATE_MATRIX" \
        python "${SCRIPTS}/ValidateMatrix.py" \
            -o . -a "${ASSETS}" -b . \
            -g gB_matrix_raw.tsv \
            -m "${ASSETS}/host_mapping.tsv" \
            -n "${ASSETS}/country_mapping.tsv" \
            -c "${ASSETS}/m49_country.csv"
fi

# ---------------------------------------------------------------------------
# STEP 4 – FILTER_AND_EXTRACT_SEQUENCES
# ---------------------------------------------------------------------------
STEP=$(( STEP + 1 ))
if (( STEP >= START_STEP )); then
    FILTER_EXTRA=()
    [[ "$UPDATE_MODE" == "true" ]] && FILTER_EXTRA+=( --update_db "${UPDATE_DB}" )
    run_step "FILTER_AND_EXTRACT" \
        python "${SCRIPTS}/FilterAndExtractSequences.py" \
            -b . -o . \
            -r "${REF_LIST}" \
            -v "${IS_SEGMENTED}" \
            -g "${GB_MATRIX}" \
            -sf "${SEQUENCES_FA}" \
            "${FILTER_EXTRA[@]+${FILTER_EXTRA[@]}}"
fi

# ---------------------------------------------------------------------------
# STEP 5 – BLAST_ALIGNMENT
# ---------------------------------------------------------------------------
STEP=$(( STEP + 1 ))
if (( STEP >= START_STEP )); then
    BLAST_EXTRA=()
    [[ "$UPDATE_MODE" == "true" ]] && BLAST_EXTRA+=( --update_db "${UPDATE_DB}" )
    if [[ "$IS_SEGMENTED" == "Y" ]]; then
        run_step "BLAST_ALIGNMENT" \
            python "${SCRIPTS}/BlastAlignment.py" \
                -s Y \
                -f "${REF_LIST}" \
                -q "${QUERY_FA}" -r "${REF_FA}" \
                -t . -b . -m "${REF_LIST}" \
                -g "${GB_MATRIX}" \
                "${BLAST_EXTRA[@]+${BLAST_EXTRA[@]}}"
    else
        run_step "BLAST_ALIGNMENT" \
            python "${SCRIPTS}/BlastAlignment.py" \
                -f "${REF_LIST}" \
                -q "${QUERY_FA}" -r "${REF_FA}" \
                -b . -t . -m "${REF_LIST}" \
                -g "${GB_MATRIX}" \
                "${BLAST_EXTRA[@]+${BLAST_EXTRA[@]}}"
    fi
fi

# ---------------------------------------------------------------------------
# STEP 6 – NEXTALIGN_ALIGNMENT
# ---------------------------------------------------------------------------
STEP=$(( STEP + 1 ))
if (( STEP >= START_STEP )); then
    NEXTALIGN_EXTRA=()
    [[ "$UPDATE_MODE" == "true" ]] && NEXTALIGN_EXTRA+=( --update_db "${UPDATE_DB}" )
    run_step "NEXTALIGN_ALIGNMENT" \
        python "${SCRIPTS}/NextalignAlignment.py" \
            -r "${REF_SEQS_DIR}" \
            -q "${GROUPED_FASTA_DIR}" \
            -g "${GB_MATRIX}" \
            -t . \
            -f "${REF_SEQ_FILTERED_FA}" \
            -m "${REF_LIST}" \
            -ms "${MASTER_SEQ_DIR}" \
            "${NEXTALIGN_EXTRA[@]+${NEXTALIGN_EXTRA[@]}}"
fi

# ---------------------------------------------------------------------------
# STEP 7 – COLLECT_FILTERED_SEQUENCES
# ---------------------------------------------------------------------------
STEP=$(( STEP + 1 ))
if (( STEP >= START_STEP )); then
    run_step "COLLECT_FILTERED_SEQUENCES" \
        python "${SCRIPTS}/CollectFilteredSequences.py" \
            -n "${NEXTALIGN_DIR}" \
            -o filtered_sequences.tsv \
            -b . \
            --max_gap_proportion "${MAX_ALN_GAP_PROPORTION}"
fi
# Guarantee files exist (CollectFilteredSequences may not create them if nothing is filtered)
touch "${FILTERED_IDS}" "${FILTERED_TSV}"

# ---------------------------------------------------------------------------
# STEP 8 – PAD_ALIGNMENT
# ---------------------------------------------------------------------------
STEP=$(( STEP + 1 ))
if (( STEP >= START_STEP )); then
    run_step "PAD_ALIGNMENT" \
        python "${SCRIPTS}/PadAlignment.py" \
            -nd "${NEXTALIGN_DIR}" \
            -m "${REF_LIST}" \
            -o . -d . \
            -i "${NEXTALIGN_DIR}/query_aln" \
            --keep_intermediate_files \
            --update_db "${UPDATE_DB:-null}" \
            --segment_manifest_out pad_alignment_manifest.tsv \
            --skip_ids "${FILTERED_IDS}"
fi
# Always rediscover after step 8 in case it just ran
mapfile -t PADDED_MSA_FILES < <(find "${WORK_DIR}" -maxdepth 1 -iname '*_merged_msa.fasta' | sort)

# ---------------------------------------------------------------------------
# STEP 9 – DEDUP_ALIGNMENT (seqkit rmdup) + optional test-mode subsampling
# ---------------------------------------------------------------------------
STEP=$(( STEP + 1 ))
if (( STEP >= START_STEP )); then
    DEDUP_FILES=()
    for msa_file in "${PADDED_MSA_FILES[@]}"; do
        base="$(basename "${msa_file}" .fasta)"
        dedup_out="${WORK_DIR}/${base}_dedup.fasta"
        run_step "DEDUP_ALIGNMENT (${base})" \
            seqkit rmdup -n "${msa_file}" -o "${dedup_out}"

        # Mirror TEST_SUBSAMPLE_CLUSTER_INPUT: cap cluster input in test mode
        cluster_input="${dedup_out}"
        if [[ "$TEST_MODE" == "1" && -n "$TEST_MAX_CLUSTER_SEQS" && "$TEST_MAX_CLUSTER_SEQS" -gt 0 ]]; then
            TOTAL_SEQS=$(seqkit seq -n "${dedup_out}" | wc -l)
            cap_out="${WORK_DIR}/${base}_dedup_cluster_input.fasta"
            if [[ "$TOTAL_SEQS" -le "$TEST_MAX_CLUSTER_SEQS" ]]; then
                echo "[test-mode] ${TOTAL_SEQS} seqs <= ${TEST_MAX_CLUSTER_SEQS}; keeping all"
                cp "${dedup_out}" "${cap_out}"
            else
                echo "[test-mode] ${TOTAL_SEQS} seqs > ${TEST_MAX_CLUSTER_SEQS}; subsampling"
                seqkit head -n "${TEST_MAX_CLUSTER_SEQS}" "${dedup_out}" > "${cap_out}"
            fi
            cluster_input="${cap_out}"
        fi
        DEDUP_FILES+=("${cluster_input}")
    done
fi

# ---------------------------------------------------------------------------
# STEPS 10-12 – MMSEQS_CLUSTERING → VERY_FAST_TREE → IQ_TREE → USHER_PLACEMENT
#               Update mode: skip MMseqs/VFT/IQ-TREE; run USHER_PLACEMENT only
# ---------------------------------------------------------------------------
STEP=$(( STEP + 1 ))  # step 10
if (( STEP >= START_STEP )); then
    IQTREE_DIRS=()
    USHER_DIRS=()
    MMSEQ_DIRS=()
    mkdir -p "${WORK_DIR}/mmseqs_inputs"

    if [[ "$TREE_FREE" == "true" && "$UPDATE_MODE" == "true" ]]; then
        echo "[info] tree_free + update_db: retaining existing trees from update DB"
    elif [[ "$UPDATE_MODE" == "true" ]]; then
        # Update mode: USHER only (no MMseqs / VFT / IQ-TREE rebuild)
        for cluster_input in "${DEDUP_FILES[@]}"; do
            base="$(basename "${cluster_input}" .fasta)"
            usher_out_dir="${WORK_DIR}/UsherUpdate_${base}"
            run_step "USHER_PLACEMENT/update (${base})" \
                python "${SCRIPTS}/UsherPlacement.py" \
                    --padded_aln "${cluster_input}" \
                    --mmseq_cluster_dir "UNSET" \
                    --iqtree_dir "UNSET" \
                    --output_dir "${usher_out_dir}" \
                    --update_db "${UPDATE_DB}" \
                    --threads "${MAX_THREADS}" \
                    --test_mode "${TEST_MODE}"
            USHER_DIRS+=("${usher_out_dir}")
        done
    elif [[ "$TREE_FREE" == "true" ]]; then
        echo "[info] tree_free mode: skipping clustering and tree steps"
    else
        # Normal fresh build
        for cluster_input in "${DEDUP_FILES[@]}"; do
            base="$(basename "${cluster_input}" .fasta)"
            mmseq_dir="${WORK_DIR}/MMseqClusters_${base}"

            # Copy input into a staging dir (mirrors Nextflow process staging)
            mmseq_input_dir="${WORK_DIR}/mmseqs_input_${base}"
            mkdir -p "${mmseq_input_dir}"
            cp "${cluster_input}" "${mmseq_input_dir}/"

            # Remove any stale output from a previous run — MMseqs errors if output DBs already exist
            rm -rf "${mmseq_dir}"

            run_step "MMSEQS_CLUSTERING (${base})" \
                python "${SCRIPTS}/MMseqsClustering.py" \
                    -i "${mmseq_input_dir}" \
                    -o "${mmseq_dir}" \
                    --min-seq-id "${MMSEQS_MIN_SEQ_ID}" \
                    --threads "${MAX_THREADS}"
            MMSEQ_DIRS+=("${mmseq_dir}")

            # Symlink into mmseqs_inputs staging dir for BuildTreeManifest later
            ln -sfn "${mmseq_dir}" "${WORK_DIR}/mmseqs_inputs/$(basename "${mmseq_dir}")"

            # Find cluster representative FASTA
            CLUSTER_REP="$(find -L "${mmseq_dir}" -name '*_cluster_rep.fasta' -print -quit 2>/dev/null || true)"
            [[ -z "$CLUSTER_REP" ]] && CLUSTER_REP="$(find -L "${mmseq_dir}" -name '*.fasta' -print -quit)"

            # VERY_FAST_TREE: generate guide tree
            guide_tree_dir="${WORK_DIR}/VeryFastTree_MMseqClusters_${base}"
            mkdir -p "${guide_tree_dir}"
            guide_tree="${guide_tree_dir}/guide_tree.nwk"
            run_step "VERY_FAST_TREE (${base})" bash -c "
                if command -v VeryFastTree >/dev/null 2>&1; then
                    VeryFastTree -threads ${MAX_THREADS} -nt -gtr -double-precision '${CLUSTER_REP}' > '${guide_tree}'
                elif command -v FastTree >/dev/null 2>&1; then
                    FastTree -nt -gtr '${CLUSTER_REP}' > '${guide_tree}'
                else
                    echo '[error] Neither VeryFastTree nor FastTree found in PATH' >&2
                    exit 1
                fi"

            # IQ_TREE: maximum-likelihood tree
            iqtree_dir="${WORK_DIR}/IQTree_MMseqClusters_${base}"
            rm -rf "${iqtree_dir}"
            mkdir -p "${iqtree_dir}"
            if command -v iqtree3 >/dev/null 2>&1; then
                IQTREE_BIN="iqtree3"
            else
                echo "[error] iqtree3 not found in PATH" >&2; exit 1
            fi
            run_step "IQ_TREE (${base})" bash -c "
                ulimit -s unlimited
                '${IQTREE_BIN}' -s '${CLUSTER_REP}' -t '${guide_tree}' \
                    -T ${MAX_THREADS} -m GTR \
                    -pre '${iqtree_dir}/iqtree' -mem 300G"
            IQTREE_DIRS+=("${iqtree_dir}")

            # USHER_PLACEMENT
            usher_out_dir="${WORK_DIR}/Usher_MMseqClusters_${base}"
            rm -rf "${usher_out_dir}"
            run_step "USHER_PLACEMENT (${base})" \
                python "${SCRIPTS}/UsherPlacement.py" \
                    --padded_aln "${cluster_input}" \
                    --mmseq_cluster_dir "${mmseq_dir}" \
                    --iqtree_dir "${iqtree_dir}" \
                    --output_dir "${usher_out_dir}" \
                    --update_db "${UPDATE_DB:-null}" \
                    --threads "${MAX_THREADS}" \
                    --test_mode "${TEST_MODE}"
            USHER_DIRS+=("${usher_out_dir}")
        done
    fi
fi

# ---------------------------------------------------------------------------
# STEP 11 – BUILD_TREE_MANIFEST
# ---------------------------------------------------------------------------
STEP=$(( STEP + 1 ))
if (( STEP >= START_STEP )); then
    mkdir -p "${IQTREE_INPUT_DIR}" "${USHER_INPUT_DIR}"
    for d in "${IQTREE_DIRS[@]+"${IQTREE_DIRS[@]}"}"; do
        ln -sfn "$d" "${IQTREE_INPUT_DIR}/$(basename "$d")"
    done
    for d in "${USHER_DIRS[@]+"${USHER_DIRS[@]}"}"; do
        ln -sfn "$d" "${USHER_INPUT_DIR}/$(basename "$d")"
    done

    run_step "BUILD_TREE_MANIFEST" \
        python "${SCRIPTS}/BuildTreeManifest.py" \
            --output "${TREE_MANIFEST}" \
            --iqtree-dir "${IQTREE_INPUT_DIR}" \
            --usher-dir "${USHER_INPUT_DIR}"
fi

# ---------------------------------------------------------------------------
# STEP 12 – CALC_ALIGNMENT_CORD
# ---------------------------------------------------------------------------
STEP=$(( STEP + 1 ))
if (( STEP >= START_STEP )); then
    # Create a staging dir containing only the merged MSA fastas (mirrors Nextflow process)
    mkdir -p "${PADDED_ALN_STAGING}"
    for _msa in "${PADDED_MSA_FILES[@]}"; do
        cp "${_msa}" "${PADDED_ALN_STAGING}/"
    done

    # -m is the ref list file (Nextflow resolves master_acc_str → master_file_opt when file exists)
    CALC_EXTRA=()
    if [[ "$UPDATE_MODE" == "true" ]]; then
        CALC_EXTRA+=( --update_db "${UPDATE_DB}" --update_scope_tsv "${GB_MATRIX}" --segment_map_tsv "${GB_MATRIX}" )
    fi

    run_step "CALC_ALIGNMENT_CORD" \
        python "${SCRIPTS}/CalcAlignmentCord.py" \
            -i "${PADDED_ALN_STAGING}" \
            -m "${REF_LIST}" \
            -g "${GFF_FILE}" \
            -bh "${BLAST_TOPHITS}" \
            -b . -d . \
            -o "${FEATURES_TSV}" \
            "${CALC_EXTRA[@]+${CALC_EXTRA[@]}}"
fi

# ---------------------------------------------------------------------------
# STEP 13 – SOFTWARE_VERSION
# ---------------------------------------------------------------------------
STEP=$(( STEP + 1 ))
if (( STEP >= START_STEP )); then
    run_step "SOFTWARE_VERSION" \
        python "${SCRIPTS}/SoftwareVersion.py" \
            -d . -o Software_info -f software_info.tsv
fi

# ---------------------------------------------------------------------------
# STEP 14 – HOST_TAXA_TABLE
# ---------------------------------------------------------------------------
STEP=$(( STEP + 1 ))
if (( STEP >= START_STEP )); then
    # Warn if taxonomy dump files are missing — they are not generated by this pipeline;
    # they must be pre-downloaded (e.g. via NCBI taxdump) and placed in Taxa/.
    if [[ ! -f "${WORK_DIR}/Taxa/names.dmp" || ! -f "${WORK_DIR}/Taxa/nodes.dmp" ]]; then
        echo "[warn] Taxonomy dump files not found in ${WORK_DIR}/Taxa/ – HOST_TAXA_TABLE may fail." >&2
        echo "[warn] Download from: https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz" >&2
    fi
    run_step "HOST_TAXA_TABLE" \
        python "${SCRIPTS}/HostTaxaTable.py" \
            -g "${GB_MATRIX}" \
            -n "${WORK_DIR}/Taxa/names.dmp" \
            -s "${WORK_DIR}/Taxa/nodes.dmp" \
            -b . -o HostTaxa
fi

# ---------------------------------------------------------------------------
# STEP 15 – GENERATE_TABLES
# ---------------------------------------------------------------------------
STEP=$(( STEP + 1 ))
if (( STEP >= START_STEP )); then
    EXTRA_TABLE_ARGS=()
    [[ -f "$REF_LIST" ]] && EXTRA_TABLE_ARGS+=( -r "${REF_LIST}" )

    # Build -p flags: GenerateTables.py expects the padded MSA FASTA file(s),
    # not the work directory.  Mirrors PAD_ALIGNMENT.out.merged_msa.collect() in Nextflow.
    PADDED_ALN_ARGS=()
    for _f in "${PADDED_MSA_FILES[@]}"; do
        PADDED_ALN_ARGS+=( -p "$_f" )
    done

    run_step "GENERATE_TABLES" \
        python "${SCRIPTS}/GenerateTables.py" \
            -g "${GB_MATRIX}" \
            -bh "${BLAST_TOPHITS}" \
            "${PADDED_ALN_ARGS[@]}" \
            -n "${NEXTALIGN_DIR}" \
            -b . -o Tables \
            -e "${EMAIL}" \
            "${EXTRA_TABLE_ARGS[@]+${EXTRA_TABLE_ARGS[@]}}"
fi

# ---------------------------------------------------------------------------
# STEP 16 – CREATE_SQLITE_DB (+ optional ANNOTATE_MUTATIONS)
# ---------------------------------------------------------------------------
STEP=$(( STEP + 1 ))
if (( STEP >= START_STEP )); then
    IQTREE_ARG=()
    USHER_ARG=()
    CLUSTER_ARG=()
    FILTERED_ARG=()
    FILTERED_DETAILS_ARG=()
    REFERENCE_ARG=()
    UPDATE_ARGS=()
    TREE_MANIFEST_ARG=()

    IQTREE_TREEFILE="$(find -L "${IQTREE_INPUT_DIR}" -name '*.treefile' -print -quit 2>/dev/null || true)"
    [[ -n "$IQTREE_TREEFILE" ]] && IQTREE_ARG=( -it "${IQTREE_TREEFILE}" )

    USHER_TREEFILE=""
    for ud in "${USHER_DIRS[@]+"${USHER_DIRS[@]}"}"; do
        if [[ -f "${ud}/final-tree.nh" ]]; then
            USHER_TREEFILE="${ud}/final-tree.nh"; break
        elif [[ -f "${ud}/uncondensed-final-tree.nh" ]]; then
            USHER_TREEFILE="${ud}/uncondensed-final-tree.nh"; break
        fi
    done
    [[ -n "$USHER_TREEFILE" ]] && USHER_ARG=( -ut "${USHER_TREEFILE}" )

    MERGED_CLUSTER_TSV="${WORK_DIR}/merged_mmseqs_clusters.tsv"
    find -L "${WORK_DIR}/mmseqs_inputs" -type f -name '*_clusters.tsv' -print0 2>/dev/null \
        | sort -z | xargs -0 -r cat > "${MERGED_CLUSTER_TSV}" || true
    [[ -s "$MERGED_CLUSTER_TSV" ]] && CLUSTER_ARG=( -ct "${MERGED_CLUSTER_TSV}" -ci "${MMSEQS_MIN_SEQ_ID}" )

    [[ -f "$FILTERED_IDS" && -s "$FILTERED_IDS" ]] && FILTERED_ARG=( -fi "${FILTERED_IDS}" )
    [[ -f "$FILTERED_TSV" ]] && FILTERED_DETAILS_ARG=( -fd "${FILTERED_TSV}" )
    [[ -f "$REF_LIST" ]] && REFERENCE_ARG=( --reference_tsv "${REF_LIST}" )
    [[ "$UPDATE_MODE" == "true" ]] && \
        UPDATE_ARGS=( --update --update_db "${UPDATE_DB}" --batch_id "sh_run_$(date +%Y%m%d_%H%M%S)" )

    MANIFEST_ROWS="$(wc -l < "${TREE_MANIFEST}" 2>/dev/null || echo 0)"
    [[ "$MANIFEST_ROWS" -gt 1 ]] && TREE_MANIFEST_ARG=( --tree_manifest "${TREE_MANIFEST}" )

    # EPA-ng fallback: only when this run produced NO usable tree. With a tree the
    # tree neighbourhood is both cheaper and more consistent, so EPA-ng stays off.
    # tree_free is excluded deliberately: it opts out of phylogenetics for speed,
    # and EPA-ng (reference tree build + placement of every query) is the slowest
    # way to assign clades. Those runs fall through to the BLAST top hit.
    CLADE_ARG=()
    if [[ "$TREE_FREE" == "true" ]]; then
        echo "[info] tree_free mode: skipping EPA-ng clade assignment; using BLAST top-hit genotypes."
    elif [[ -z "$IQTREE_TREEFILE" && -z "$USHER_TREEFILE" && -f "$REF_LIST" \
          && ${#PADDED_MSA_FILES[@]} -gt 0 && -f "${PADDED_MSA_FILES[0]}" ]]; then
        echo "[info] No IQ-TREE/UShER tree in this run; running EPA-ng clade assignment fallback."
        # Clade labels come from the ref_list genotype/subtype columns for every
        # dataset - no separate per-organism taxon files.
        if python "${SCRIPTS}/CladeAssignment.py" \
                -p "${PADDED_MSA_FILES[0]}" -b "${WORK_DIR}" -o CladeAssignment \
                -r "${REF_LIST}" -a "${WORK_DIR}/clade_assignments.tsv" \
                -t "${THREADS}"; then
            CLADE_ARG=( --clade_assignments "${WORK_DIR}/clade_assignments.tsv" )
        else
            echo "[warn] EPA-ng clade assignment failed; falling back to BLAST top-hit genotypes." >&2
        fi
    fi

    run_step "CREATE_SQLITE_DB" \
        python "${SCRIPTS}/CreateSqliteDB.py" \
            -m "${GB_MATRIX}" \
            -rf "${FEATURES_TSV}" \
            -p "${SEQ_ALN}" \
            -i "${INSERTIONS}" \
            -ht "${HOST_TAXA}" \
            -s "${SOFTWARE_INFO}" \
            -fa "${SEQUENCES_FA}" \
            -mc "${ASSETS}/m49_country.csv" \
            -mir "${ASSETS}/m49_intermediate_region.csv" \
            -mr "${ASSETS}/m49_region.csv" \
            -msr "${ASSETS}/m49_sub_region.csv" \
            -d "${DB_NAME}" \
            -b . -o . \
            "${IQTREE_ARG[@]}" \
            "${USHER_ARG[@]}" \
            "${CLUSTER_ARG[@]}" \
            "${FILTERED_ARG[@]}" \
            "${FILTERED_DETAILS_ARG[@]}" \
            "${TREE_MANIFEST_ARG[@]}" \
            "${REFERENCE_ARG[@]}" \
            "${CLADE_ARG[@]+"${CLADE_ARG[@]}"}" \
            "${UPDATE_ARGS[@]}"

    # ANNOTATE_MUTATIONS (optional – only when mutation_catalog provided)
    if [[ -n "$MUTATION_CATALOG" && -f "$MUTATION_CATALOG" ]]; then
        VIRUS_ARG=()
        [[ -n "$MUTATION_VIRUS" ]] && VIRUS_ARG=( --virus "${MUTATION_VIRUS}" )
        CATALOG_PROFILE="${MUTATION_VIRUS:-all_columns}"
        run_step "ANNOTATE_MUTATIONS" \
            python "${SCRIPTS}/AnnotateMutations.py" \
                --db "${SQLITE_DB}" \
                --mutation_catalog "${MUTATION_CATALOG}" \
                --catalog_column_profile "${CATALOG_PROFILE}" \
                "${VIRUS_ARG[@]}"
    fi
fi

# ---------------------------------------------------------------------------
# STEP 17 – VALIDATE_DB_TREE  (test mode only)
# ---------------------------------------------------------------------------
STEP=$(( STEP + 1 ))
if (( STEP >= START_STEP && TEST_MODE == 1 )); then
    mkdir -p "${PUBLISH_DIR}/tests"
    VALIDATE_EXTRA=( --check-update-integrity )
    [[ "$IS_SEGMENTED" == "Y" ]] && VALIDATE_EXTRA+=( --expect-segment-trees )
    [[ "$TREE_FREE" == "true" ]] && VALIDATE_EXTRA+=( --allow-no-trees )
    run_step "VALIDATE_DB_TREE" \
        python "${SCRIPTS}/ValidateDbTree.py" \
            --db "${SQLITE_DB}" \
            --outdir "${PUBLISH_DIR}/tests" \
            --test-mode \
            "${VALIDATE_EXTRA[@]}"
fi

# ---------------------------------------------------------------------------
# STEP 18 – TEST_NON_SEGMENTED_OUTPUT (test mode, non-segmented only)
# ---------------------------------------------------------------------------
STEP=$(( STEP + 1 ))
if (( STEP >= START_STEP )) && [[ "$TEST_MODE" == "1" && "$IS_SEGMENTED" == "N" ]]; then
    run_step "TEST_NON_SEGMENTED_OUTPUT" \
        python "${SCRIPTS}/TestPipelineOutput.py" \
            --mode non_segmented \
            --blast_hits "${BLAST_TOPHITS}" \
            --gb_matrix "${GB_MATRIX}" \
            --sqlite_db "${SQLITE_DB}" \
            --output "${PUBLISH_DIR}/tests/test_non_segmented_results.txt"
fi

# ---------------------------------------------------------------------------
# STEP 19 – VERIFY_MUTATIONS  (only when mutation_catalog provided)
# ---------------------------------------------------------------------------
STEP=$(( STEP + 1 ))
if (( STEP >= START_STEP )) && [[ -n "$MUTATION_CATALOG" && -f "$MUTATION_CATALOG" ]]; then
    run_step "VERIFY_MUTATIONS" bash -c "
        python '${SCRIPTS}/VerifyMutations.py' \
            --db '${SQLITE_DB}' \
            --mutation_catalog '${MUTATION_CATALOG}' 2>&1 \
        | tee '${PUBLISH_DIR}/tests/mutation_verification.txt'"
fi

# ---------------------------------------------------------------------------
# Publish outputs to PUBLISH_DIR
# ---------------------------------------------------------------------------
cp -f "${SQLITE_DB}" "${PUBLISH_DIR}/${DB_NAME}.db"
[[ -f "${WORK_DIR}/db_summary.txt" ]] && cp -f "${WORK_DIR}/db_summary.txt" "${PUBLISH_DIR}/db_summary.txt"

# ---------------------------------------------------------------------------
# Done
# ---------------------------------------------------------------------------
echo ""
echo "================================================================================"
echo "Pipeline completed successfully."
echo "Output DB : ${PUBLISH_DIR}/${DB_NAME}.db"
[[ -f "${PUBLISH_DIR}/db_summary.txt" ]] && echo "" && cat "${PUBLISH_DIR}/db_summary.txt"
echo "================================================================================"

