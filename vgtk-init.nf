
//params defined in nextflow.config override with --params..

// profiles=conda,condaMamba,test,setup_rabv_full
//use conda if not running in conda env alraedy,
// use condaMamba for mamba if installed (check with mamba --version)

scripts_dir     = "${projectDir}/scripts"

def parsePositiveInt = { value, paramName ->
    try {
        def n = (value as Integer)
        if( n < 1 ){
            error("ERROR: params.${paramName} must be a positive integer, got: ${value}")
        }
        return n
    } catch(Exception e){
        error("ERROR: params.${paramName} must be a positive integer, got: ${value}")
    }
}

def maxThreadsRaw = params.max_threads ?: Math.min(8, Runtime.getRuntime().availableProcessors())
def MAX_THREADS = parsePositiveInt(maxThreadsRaw, 'max_threads')

def MMSEQS_THREADS_REQUESTED = parsePositiveInt(params.mmseqs_threads ?: MAX_THREADS, 'mmseqs_threads')
params.mmseqs_threads = Math.min(MMSEQS_THREADS_REQUESTED, MAX_THREADS)

def SEGMENT_PARALLEL_THREADS = Math.max(1, (int)Math.floor(MAX_THREADS / 8))
// 1. List your script's explicitly defined parameters (keep this in sync!)
def scriptDefinedParams = [
    'tax_id', 'db_name', 'is_segmented', 'extra_info_fill', 'test',
    "scripts_dir", "publish_dir", "email", "ref_list", "bulk_fillup_table", "is_flu", "gene_info",
    "xml_dir", "update", "update_file",
    "mmseqs_min_seq_id", "mmseqs_threads", "mmseqs_trim_cds_file",
    "gisaid_dir", "previous_db", "conda_path", "test_max_cluster_seqs", "max_threads", "ref_set_aligned"
    // Add all parameter names defined above
]

// 2. List common/allowed built-in Nextflow parameters (customize as needed)
//    You might need to add others depending on features you use (e.g., cloud options)
def knownNextflowParams = [
    'help', 'version', 'profile', 'config', 'workDir', 'resume',
    'entry', 'validateParams', 'cached', 'dumpHashes', 'listInputs',
    // Add cloud-specific ones if used: 'awsqueue', 'clusterOptions', etc.
    // Add tower related ones if used: 'tower', 'computeEnvId', etc.
    // Check Nextflow documentation for a more exhaustive list if required
]

// 3. Combine allowed parameters
def allowedParams = (scriptDefinedParams + knownNextflowParams).unique()

// Backward-compatible alias: update_file -> update
if( params.update_file && !params.update ){
    params.update = params.update_file
}

// Normalize optional params used as `val` process inputs
if( params.ref_set_aligned == null ){
    params.ref_set_aligned = 'UNSET'
}

// 4. Get all parameters provided via command line and config files
//    The 'params' map holds the merged view of all parameters
def providedParams = params.collect { it.key }

// 5. Find parameters that were provided but are not in the allowed list
def unexpectedParams = providedParams.findAll { ! (it in allowedParams) }

// 6. Error out if any unexpected parameters were found
if (unexpectedParams) {
    // Construct a helpful error message
    def errorMsg = """
    ERROR: Unknown command line parameter(s) provided: ${unexpectedParams.join(', ')}

    Please check your command. Allowed parameters are:
    Script parameters: ${scriptDefinedParams.join(', ')}
    Nextflow options : ${knownNextflowParams.join(', ')} (subset shown)
    """.stripIndent() // Use stripIndent for cleaner multi-line string formatting

    error(errorMsg) // Stop the pipeline
}

process TEST_DEPENDENCIES{
    output:
        path "dependency_test.txt"
    shell:
    '''
    set +e

    missing_count=0
    {
        echo "=== vgtk dependency preflight ==="
        echo "date: $(date -Iseconds)"
        echo "python executable: $(command -v python || echo 'NOT FOUND')"
        echo "python version: $(python --version 2>&1)"
        echo "CONDA_PREFIX=${CONDA_PREFIX:-<not set>}"
        echo
        echo "[required command-line tools]"
    } > dependency_test.txt

    check_cmd() {
        local cmd="$1"
        shift
        if command -v "$cmd" >/dev/null 2>&1; then
            local ver
            ver="$($cmd "$@" 2>&1 | head -n 1)"
            echo "OK   $cmd :: ${ver}" >> dependency_test.txt
        else
            echo "MISS $cmd :: not found in PATH" >> dependency_test.txt
            missing_count=$((missing_count + 1))
        fi
    }

    # Required by the current pipeline processes
    check_cmd python --version
    check_cmd blastn -version
    check_cmd makeblastdb -version
    check_cmd seqkit version
    check_cmd nextalign --version
    check_cmd mmseqs --version
    if command -v iqtree3 >/dev/null 2>&1; then
        iqv="$(iqtree3 --version 2>&1 | head -n 1)"
        echo "OK   iqtree3 :: ${iqv}" >> dependency_test.txt
        if echo "$iqv" | grep -qi 'COVID-edition'; then
            echo "MISS iqtree3 :: COVID-edition build detected (use standard iqtree3 build)" >> dependency_test.txt
            missing_count=$((missing_count + 1))
        fi
    else
        echo "MISS iqtree3 :: not found in PATH" >> dependency_test.txt
        missing_count=$((missing_count + 1))
    fi
    check_cmd usher --version
    check_cmd faToVcf 2>/dev/null
    check_cmd efetch -version

    {
        echo
        echo "[required python modules]"
    } >> dependency_test.txt

    check_py() {
        local module="$1"
        local label="$2"
        python - <<PY >> dependency_test.txt 2>&1
import importlib
module = "${module}"
label = "${label}"
try:
    m = importlib.import_module(module)
    ver = getattr(m, "__version__", "<no __version__>")
    print(f"OK   {label} ({module}) :: {ver}")
except Exception as e:
    print(f"MISS {label} ({module}) :: {e}")
    raise
PY
        if [ $? -ne 0 ]; then
            missing_count=$((missing_count + 1))
        fi
    }

    check_py Bio biopython
    check_py pandas pandas
    check_py numpy numpy
    check_py openpyxl openpyxl
    check_py requests requests
    check_py dateutil python-dateutil
    check_py matplotlib matplotlib
    check_py tqdm tqdm

    {
        echo
        echo "[optional tools used by utility/legacy scripts]"
    } >> dependency_test.txt

    if command -v VeryFastTree >/dev/null 2>&1; then
        echo "OK   VeryFastTree :: $(VeryFastTree -help 2>&1 | head -n 1)" >> dependency_test.txt
    elif command -v FastTree >/dev/null 2>&1; then
        echo "OK   FastTree :: $(FastTree -help 2>&1 | head -n 1)" >> dependency_test.txt
    else
        echo "WARN VeryFastTree/FastTree :: not found (only needed if VERY_FAST_TREE is enabled)" >> dependency_test.txt
    fi

    if command -v mafft >/dev/null 2>&1; then
        echo "OK   mafft :: $(mafft --version 2>&1 | head -n 1)" >> dependency_test.txt
    else
        echo "WARN mafft :: not found (used by utility scripts, not main flow)" >> dependency_test.txt
    fi

    if command -v ncbi-acc-download >/dev/null 2>&1; then
        echo "OK   ncbi-acc-download :: $(ncbi-acc-download --help 2>&1 | head -n 1)" >> dependency_test.txt
    else
        echo "WARN ncbi-acc-download :: not found (pip helper; not required for all runs)" >> dependency_test.txt
    fi

    echo >> dependency_test.txt
    echo "missing_required_dependencies=${missing_count}" >> dependency_test.txt

    if [ "$missing_count" -gt 0 ]; then
        echo "[ERROR] Missing required dependencies: ${missing_count}" >&2
        cat dependency_test.txt >&2
        exit 1
    fi
    '''
}


process VALIDATE_REF_LIST{
    input:
        val ref_list
        val is_segmented
    output:
        path "ref_list_validated.txt"
    shell:
    '''
    # check there are three columns, 
    #1st is accession, 2nd contains master or reference, 3rd is segment (if segmented)
        python !{scripts_dir}/ValidateRefList.py -r !{ref_list} -s !{is_segmented} -o ref_list_validated.txt

    '''
}
process FETCH_GENBANK{
    publishDir "${params.publish_dir}"
    input:
        val TAX_ID
        path ref_list
    output:
        path 'GenBank-XML', type: 'dir', emit: gen_bank_XML
    shell:
    '''
    extra=""
    if( [ "!{params.test}" -eq "1" ] )
    then
        extra="${extra} --test_run --ref_list !{ref_list}"
    fi

    if [ "!{params.update}" != "null" ] && [ -n "!{params.update}" ]; then
        extra="${extra} --update !{params.update}"
    fi

    python !{scripts_dir}/GenBankFetcher.py --taxid !{TAX_ID} -b 50 \
             ${extra} -e !{params.email} -o . 
    #--update tmp/GenBank-matrix/gB_matrix_raw.tsv is gonna be problematic for this!
    #what's update doing with a tmp dir?

    '''
}

process DOWNLOAD_GFF{
    input:
        val master_acc
        path master_file_opt
    output:
         path "*.gff3"
    shell:
    '''
    MASTER_ARG="!{master_acc}"
    if [ -f "!{master_file_opt}" ]; then
        MASTER_ARG="!{master_file_opt}"
    fi
    python "!{scripts_dir}/DownloadGFF.py" --accession_ids "$MASTER_ARG" -o . -b .

    '''
}


process GENBANK_PARSER{
    publishDir "${params.publish_dir}"
    input:
        path ref_list_path
        path gen_bank_XML
    output:
        path "gB_matrix_validated.tsv" , emit: gb_matrix
        path "sequences.fa", emit: sequences_out
    shell:
    '''
        echo "DEBUG: params.test is '!{params.test}'"
        echo "DEBUG: params.xml_dir is '!{params.xml_dir}'"
        
        # force re-run after update
        # Always require refs to ensure master sequences range are found (even in test mode)
        extra="--require_refs"
        
        # Logic: If xml_dir is provided (mimicking old XML_source=XML), apply test flags if needed
        if [ -n "!{params.xml_dir}" ] && [ "!{params.xml_dir}" != "null" ]; then
            if [ "!{params.test}" -eq "1" ]; then
                extra="${extra} --test_run"
            fi
        fi
        
        echo "DEBUG: Final extra args: ${extra}"
        
        python !{scripts_dir}/GenBankParser.py -r !{ref_list_path} -d !{gen_bank_XML} -o . -b . ${extra}
        python !{scripts_dir}/ValidateMatrix.py -o . -a !{projectDir}/assets -b . \
        -g gB_matrix_raw.tsv \
        -m !{projectDir}/assets/host_mapping.tsv -n !{projectDir}/assets/country_mapping.tsv  \
        -c !{projectDir}/assets/m49_country.csv
        
    '''
}

process TIDY_GISAID {
    input:
       path gisaid_dir
       path db_file
    output:
       path "metadata.tsv", emit: gisaid_meta
       path "all_nuc.fas", emit: gisaid_nuc
    shell:
       '''
       db_arg=""
       # Check if db_file is a valid file (not a directory or empty placeholder)
       if [ -f "!{db_file}" ]; then
           db_arg="--db_file !{db_file}"
       fi

       extra_args=""
       if [ "!{params.test}" = "1" ]; then
           extra_args="--test"
       fi

       # Default to tsv if no xls/xlsx found, logic handled largely by glob inside but the arg switches mode
       # To make it robust we could check, but for now assuming tsv as per GISAID standard downloads
       python !{scripts_dir}/gisaid_tidy.py --data_dir !{gisaid_dir} --output_dir . --filetype tsv $db_arg $extra_args
       '''
}

process MERGE_GISAID {
    input:
       path gb_matrix
       path gisaid_meta
       path gisaid_nuc
       path column_mapping
    output:
       path "gB_matrix_merged.tsv", emit: merged_matrix
    shell:
       '''
       python !{scripts_dir}/merge_into_gB_matrix.py -g !{gb_matrix} \
           -t !{gisaid_meta} -f !{gisaid_nuc} -o gB_matrix_merged.tsv \
           -k Segment_Id --dataset_source gisaid -m !{column_mapping}
       '''
}

process CAT_FASTA {
    input:
        path fa1
        path fa2
    output:
        path "combined_sequences.fa", emit: combined_fa
    shell:
        '''
        cat !{fa1} !{fa2} > combined_sequences.fa
        '''
}

process ADD_MISSING_DATA{
    input:
        path gen_bank_table
    output:
        path "*.tsv", emit: tsvs_out
    when:
        params.extra_info_fill.toBoolean()
    shell:
    '''
        python !{scripts_dir}/AddMissingData.py -b !{gen_bank_table} \
         -t !{params.bulk_fillup_table}  -d . -f !{params.bulk_fillup_table}
    '''
}


process FILTER_AND_EXTRACT{
    input:
        path table_in
        path seqs_in
    output:
        path "query_seq.fa", emit: query_seqs_out
        path "ref_seq.fa", emit: ref_seqs_out
    shell:
    '''

        python !{scripts_dir}/FilterAndExtractSequences.py -b . -o . -r !{params.ref_list} \
         -v !{params.is_segmented} -g !{table_in} -sf !{seqs_in}
    '''
}



process BLAST_ALIGNMENT{
    input:
        path query_seqs
        path ref_seqs
        path gb_matrix
    output:
        path "query_tophits.tsv",type: "file", emit: query_tophits
        path "query_uniq_tophits.tsv", type: "file", emit: query_uniq_tophits
        path "query_uniq_tophit_annotated.tsv", type: "file", optional: true, emit: query_uniq_tophit_annotated
        path "grouped_fasta", type: 'dir', emit: grouped_fasta
        path "ref_seqs", type: 'dir', emit: ref_seqs_dir
        path "ref_seq_filtered.fa", type: 'file', emit: ref_seqs_fasta
        path "master_seq", type: 'dir', emit: master_seq_dir
        
    shell:
    '''
        if [ "!{params.is_segmented}" = "Y" ]; then
            python "!{scripts_dir}/BlastAlignment.py" -s Y -f "!{params.ref_list}" -q !{query_seqs} -r !{ref_seqs} \
             -t . -b . -m !{params.ref_list}  -g !{gb_matrix}
        else
            python "!{scripts_dir}/BlastAlignment.py" -f "!{params.ref_list}" -q !{query_seqs} -r !{ref_seqs} \
             -b . -t . -m !{params.ref_list} -g !{gb_matrix}
        fi
    '''
}
//workdir files:
//DB                       grouped_fasta  merged_fasta  query_tophits.tsv       ref_seq.fa  sorted_all
//gB_matrix_validated.tsv  master_seq     query_seq.fa  query_uniq_tophits.tsv  ref_seqs    sorted_fasta


//python "${scripts_dir}/NextalignAlignment.py" -m $master_acc #-gff "tmp/Gff/NC_001542.gff3"
process NEXTALIGN_ALIGNMENT{
    input:
        path genbank_matrix
        path grouped_fasta_dir
        path ref_seqs
        path ref_seqs_fasta
        path master_seq_dir
        val master_acc_str
        path master_file_opt
    output:
        path "Nextalign", type: 'dir'
    shell:
    '''
        TARGET_M="!{master_acc_str}"
        if [ -f "!{master_file_opt}" ]; then
             TARGET_M="!{master_file_opt}"
        fi

        python !{scripts_dir}/NextalignAlignment.py  -r !{ref_seqs}  \
         -q !{grouped_fasta_dir} -g !{genbank_matrix} -t . \
         -f !{ref_seqs_fasta} -m "$TARGET_M" -ms !{master_seq_dir} 
    '''
}


//"${scripts_dir}/PadAlignment-1.py" -r "/home3/sk312p/task_dir/projects/VGTK/dev_version-jun-09/TING/alUnc509RefseqsMafftHandModified.fa
process PAD_ALIGNMENT{
    publishDir "${params.publish_dir}"
    input:
        path nextalign_dir 
        val master_acc_str
        path master_file_opt
        val ref_set_aligned_dir
    output:
        path "*_merged_MSA.fasta", emit: merged_msa
    shell:
    '''
        TARGET_M="!{master_acc_str}"
        if [ -f "!{master_file_opt}" ]; then
             TARGET_M="!{master_file_opt}"
        fi

        if [ "!{ref_set_aligned_dir}" != "UNSET" ] && [ "!{ref_set_aligned_dir}" != "null" ] && [ -n "!{ref_set_aligned_dir}" ] && [ ! -d "!{ref_set_aligned_dir}" ]; then
            echo "[error] params.ref_set_aligned points to a non-existent directory: !{ref_set_aligned_dir}" >&2
            exit 1
        fi

        PRECOMP_ARGS=""
        if [ -n "!{ref_set_aligned_dir}" ] && [ "!{ref_set_aligned_dir}" != "null" ] && [ "!{ref_set_aligned_dir}" != "UNSET" ] && [ -d "!{ref_set_aligned_dir}" ]; then
            PRECOMP_ARGS="--precomputed_ref_dir !{ref_set_aligned_dir}"
            echo "Using precomputed segment alignments from !{ref_set_aligned_dir}"
        fi

        python !{scripts_dir}/PadAlignment.py -nd !{nextalign_dir} \
        -m "$TARGET_M" \
        -o . -d . -i !{nextalign_dir}/query_aln --keep_intermediate_files ${PRECOMP_ARGS}
    '''
}

process COLLECT_FILTERED_SEQUENCES {
    publishDir "${params.publish_dir}"
    input:
        path nextalign_dir
    output:
        path "filtered_sequences.tsv", emit: filtered_tsv
        path "filtered_sequences_ids.txt", emit: filtered_ids
        path "filtered_sequences_summary.txt", emit: filtered_summary
    shell:
    '''
    python !{scripts_dir}/CollectFilteredSequences.py \
        -n !{nextalign_dir} \
        -o filtered_sequences.tsv \
        -b .
    '''
}

process DEDUP_ALIGNMENT{
    publishDir "${params.publish_dir}"
    input:
        path padded_aln
    output:
        path "${padded_aln.baseName}_dedup.fasta", emit: dedup_msa
    shell:
    '''
        seqkit rmdup -n !{padded_aln} -o !{padded_aln.baseName}_dedup.fasta
    '''
}

process TEST_SUBSAMPLE_CLUSTER_INPUT {
    publishDir "${params.publish_dir}"
    input:
        path dedup_msa
    output:
        path "${dedup_msa.baseName}_cluster_input.fasta", emit: dedup_for_cluster
    shell:
    '''
        MAX_SEQS="!{params.test_max_cluster_seqs}"
        OUT_FILE="!{dedup_msa.baseName}_cluster_input.fasta"
        TOTAL_SEQS=$(seqkit seq -n "!{dedup_msa}" | wc -l)

        if [ -n "$MAX_SEQS" ] && [ "$MAX_SEQS" != "null" ] && [ "$MAX_SEQS" -gt 0 ] 2>/dev/null; then
            if [ "$TOTAL_SEQS" -le "$MAX_SEQS" ]; then
                echo "[test-mode] Input has ${TOTAL_SEQS} sequences (<= ${MAX_SEQS}); keeping all sequences for clustering"
                cp "!{dedup_msa}" "$OUT_FILE"
            else
                echo "[test-mode] Subsampling !{dedup_msa} from ${TOTAL_SEQS} to ${MAX_SEQS} sequences using deterministic random sampling (seed=42)"
                seqkit sample -n "$MAX_SEQS" -s 42 "!{dedup_msa}" -o "$OUT_FILE"
            fi
        else
            cp "!{dedup_msa}" "$OUT_FILE"
        fi
    '''
}

process MMSEQS_CLUSTERING{
    publishDir "${params.publish_dir}"
    cpus { params.is_segmented == 'Y' ? SEGMENT_PARALLEL_THREADS : MAX_THREADS }
    input:
        path padded_aln
    output:
        path "MMseqClusters_${padded_aln.baseName}", emit: mmseq_clusters
    shell:
    '''
        mkdir -p mmseqs_input
        cp !{padded_aln} mmseqs_input/

        python !{scripts_dir}/MMseqsClustering.py \
            -i mmseqs_input \
            -o MMseqClusters_!{padded_aln.baseName} \
            --min-seq-id !{params.mmseqs_min_seq_id} \
            --threads !{task.cpus}
    '''
}

process IQ_TREE{
    publishDir "${params.publish_dir}"
    cpus { params.is_segmented == 'Y' ? SEGMENT_PARALLEL_THREADS : MAX_THREADS }
    input:
        path mmseq_cluster_dir
    output:
        path "IQTree_${mmseq_cluster_dir.baseName}", emit: iqtree_out
    shell:
    '''
        ulimit -s unlimited
        CLUSTER_REP=$(find -L !{mmseq_cluster_dir} -name "*_cluster_rep.fasta" -print -quit)
        if [ -z "$CLUSTER_REP" ]; then
            CLUSTER_REP=$(find -L !{mmseq_cluster_dir} -name "*.fasta" -print -quit)
            if [ -n "$CLUSTER_REP" ]; then
                echo "[warn] *_cluster_rep.fasta not found; falling back to $CLUSTER_REP" >&2
            else
                echo "[error] No FASTA found in !{mmseq_cluster_dir}" >&2
                exit 1
            fi
        fi

        LEN_STATS=$(seqkit fx2tab -n -s "$CLUSTER_REP" | awk '
            BEGIN{count=0; min=-1; max=0}
            {
                l=length($2)
                count++
                if(min==-1 || l<min) min=l
                if(l>max) max=l
            }
            END{
                if(count==0){
                    print "0\t0\t0"
                } else {
                    print count "\t" min "\t" max
                }
            }
        ')

        ALN_COUNT=$(echo "$LEN_STATS" | cut -f1)
        ALN_MIN_LEN=$(echo "$LEN_STATS" | cut -f2)
        ALN_MAX_LEN=$(echo "$LEN_STATS" | cut -f3)

        if [ "$ALN_COUNT" -eq 0 ]; then
            echo "[error] Cluster representative FASTA is empty: $CLUSTER_REP" >&2
            exit 1
        fi

        if [ "$ALN_MIN_LEN" -ne "$ALN_MAX_LEN" ]; then
            echo "[error] Alignment length mismatch in $CLUSTER_REP (min=${ALN_MIN_LEN}, max=${ALN_MAX_LEN})." >&2
            echo "[error] IQ-TREE requires aligned sequences of equal length." >&2
            exit 1
        fi

        IQTREE_BIN=""
        if command -v iqtree3 >/dev/null 2>&1; then
            IQTREE_BIN="iqtree3"
        else
            echo "[error] iqtree3 not found in PATH. Activate the vgtk conda environment / conda profile." >&2
            exit 1
        fi

        IQTREE_VER=$($IQTREE_BIN --version 2>&1 | head -n 1 || true)
        echo "[info] Using $IQTREE_BIN at $(command -v $IQTREE_BIN): $IQTREE_VER"
        if echo "$IQTREE_VER" | grep -qi 'COVID-edition'; then
            echo "[error] Detected COVID-edition IQ-TREE build, which is incompatible here. Install standard iqtree3 in vgtk env." >&2
            exit 1
        fi

        mkdir -p IQTree_!{mmseq_cluster_dir.baseName}
        "$IQTREE_BIN" -s "$CLUSTER_REP" -nt !{task.cpus} -m GTR -pre IQTree_!{mmseq_cluster_dir.baseName}/iqtree -mem 100G
    '''
}

process USHER_PLACEMENT{
    publishDir "${params.publish_dir}"
    cpus { params.is_segmented == 'Y' ? SEGMENT_PARALLEL_THREADS : MAX_THREADS }
    input:
        tuple path(mmseq_cluster_dir), path(iqtree_dir), path(padded_aln)
    output:
        path "Usher_${mmseq_cluster_dir.baseName}", emit: usher_out
    shell:
    '''
        CLUSTER_REP=$(find -L !{mmseq_cluster_dir} -name "*_cluster_rep.fasta" -print -quit)
        TREE_FILE=$(find -L !{iqtree_dir} -name "*.treefile" -print -quit)

        if [ -z "$CLUSTER_REP" ] || [ -z "$TREE_FILE" ]; then
            echo "[error] Missing MMseqs centroid FASTA or IQ-TREE treefile." >&2
            echo "[error] centroid: $CLUSTER_REP" >&2
            echo "[error] tree:     $TREE_FILE" >&2
            exit 1
        fi

        REF_ID=$(seqkit seq -n "$CLUSTER_REP" | head -n 1)
        if [ -z "$REF_ID" ]; then
            echo "[error] Could not determine reference ID from centroid FASTA." >&2
            exit 1
        fi

        # faToVcf may fail if the chosen reference ID is not present in the input MSA
        if ! seqkit seq -n "!{padded_aln}" | grep -Fxq "$REF_ID"; then
            echo "[warn] Reference ID '$REF_ID' is not present in !{padded_aln}; using first MSA ID instead." >&2
            REF_ID=$(seqkit seq -n "!{padded_aln}" | head -n 1)
            if [ -z "$REF_ID" ]; then
                echo "[error] Could not determine fallback reference ID from !{padded_aln}." >&2
                exit 1
            fi
        fi

        mkdir -p Usher_!{mmseq_cluster_dir.baseName}
        export OMP_NUM_THREADS=!{task.cpus}
        export OPENBLAS_NUM_THREADS=!{task.cpus}
        export MKL_NUM_THREADS=!{task.cpus}
        export NUMEXPR_NUM_THREADS=!{task.cpus}
        seqkit seq -n "$CLUSTER_REP" > Usher_!{mmseq_cluster_dir.baseName}/centroid_ids.txt
        seqkit seq -n "!{padded_aln}" > Usher_!{mmseq_cluster_dir.baseName}/aln_ids.txt
        awk -v ref="$REF_ID" '$0 != ref' Usher_!{mmseq_cluster_dir.baseName}/centroid_ids.txt > Usher_!{mmseq_cluster_dir.baseName}/exclude_ids.txt

        # If exclude list removes all non-reference sequences, faToVcf can crash.
        ALN_COUNT=$(wc -l < Usher_!{mmseq_cluster_dir.baseName}/aln_ids.txt)
        EXCLUDE_COUNT=$(wc -l < Usher_!{mmseq_cluster_dir.baseName}/exclude_ids.txt)

        if [ "$ALN_COUNT" -le 1 ]; then
            echo "[error] Need at least 2 sequences in !{padded_aln}, found ${ALN_COUNT}." >&2
            exit 1
        elif [ "$EXCLUDE_COUNT" -ge $((ALN_COUNT - 1)) ]; then
            if [ "!{params.test}" = "1" ]; then
                echo "[warn] Exclude list would remove all non-reference sequences (${EXCLUDE_COUNT}/${ALN_COUNT}) in test mode; likely no query sequences to place for this segment. Running faToVcf without -excludeFile." >&2
            else
                echo "[warn] Exclude list would remove all non-reference sequences (${EXCLUDE_COUNT}/${ALN_COUNT}); likely no query sequences to place for this segment. Running faToVcf without -excludeFile." >&2
            fi
            faToVcf -ref="$REF_ID" "!{padded_aln}" Usher_!{mmseq_cluster_dir.baseName}/all_samples.vcf
        elif ! faToVcf -ref="$REF_ID" -excludeFile=Usher_!{mmseq_cluster_dir.baseName}/exclude_ids.txt \
            "!{padded_aln}" Usher_!{mmseq_cluster_dir.baseName}/all_samples.vcf; then
            echo "[warn] faToVcf failed with -excludeFile; retrying without exclude filter." >&2
            faToVcf -ref="$REF_ID" "!{padded_aln}" Usher_!{mmseq_cluster_dir.baseName}/all_samples.vcf
        fi

        if usher --help 2>&1 | grep -q -- ' -T '; then
            usher \
                -v Usher_!{mmseq_cluster_dir.baseName}/all_samples.vcf \
                -t "$TREE_FILE" \
                -d Usher_!{mmseq_cluster_dir.baseName} \
                -o Usher_!{mmseq_cluster_dir.baseName}/usher.pb \
                -C -u -T !{task.cpus}
        else
            usher \
                -v Usher_!{mmseq_cluster_dir.baseName}/all_samples.vcf \
                -t "$TREE_FILE" \
                -d Usher_!{mmseq_cluster_dir.baseName} \
                -o Usher_!{mmseq_cluster_dir.baseName}/usher.pb \
                -C -u
        fi
    '''
}

process CALC_ALIGNMENT_CORD {
    input:
        path padded_fasta
        path gff_file
        path blast_hits
        val master_acc_str
        path master_file_opt
    output:
        path "features.tsv", emit: features
    shell:
    '''
    TARGET_M="!{master_acc_str}"
    if [ -f "!{master_file_opt}" ]; then
            TARGET_M="!{master_file_opt}"
    fi

    # CalcAlignmentCord expects a directory for -i
    mkdir padded_alignments
    cp !{padded_fasta} padded_alignments/
    
    python !{scripts_dir}/CalcAlignmentCord.py -i padded_alignments \
    -m "$TARGET_M" -g !{gff_file} -bh !{blast_hits} \
    -b . -d . -o features.tsv
    '''
}

process SOFTWARE_VERSION {
    output:
        path "Software_info/software_info.tsv", emit: software_info
    shell:
    '''
    python !{scripts_dir}/SoftwareVersion.py -d . -o Software_info \
     -f software_info.tsv
    '''
}

process VERY_FAST_TREE{
    publishDir "${params.publish_dir}"
    cpus { params.is_segmented == 'Y' ? SEGMENT_PARALLEL_THREADS : MAX_THREADS }
    when: 
        params.update == null
    input:
        path padded_aln
    output:
        path "tree.nwk"
    shell:
    '''
        seqkit rmdup !{padded_aln} -o !{padded_aln}_dedup.fa
        VeryFastTree -threads !{task.cpus} -nt -gtr -double-precision !{padded_aln}_dedup.fa > tree.nwk
    '''
}


process GENERATE_TABLES {
    input:
        path gb_matrix
        path blast_hits
        path padded_aln
        path nextalign_dir
    output:
        path "Tables/sequence_alignment.tsv", emit: sequence_alignment
        path "Tables/insertions.tsv", emit: insertions
        path "Tables/host_taxa.tsv", emit: host_taxa
        path "Tables", emit: tables_dir
    shell:
    '''
    python !{scripts_dir}/GenerateTables.py -g !{gb_matrix} \
    -bh !{blast_hits} -p !{padded_aln} -n !{nextalign_dir} \
    -b . -o Tables -e !{params.email}
    '''
}

process CREATE_SQLITE_DB {
    publishDir "${params.publish_dir}"
    input:
        path meta_data
        path features
        path sequence_alignment
        path insertions
        path host_taxa
        path software_info
        path fasta_sequences
        path iqtree_dirs          // Can be multiple directories for segmented viruses
        path mmseq_cluster_dirs   // Can be multiple directories for segmented viruses
        path usher_dirs           // Can be multiple directories for segmented viruses
        path filtered_tsv
        path filtered_ids
    output:
        path "${params.db_name}.db"
    shell:
    '''
    # For segmented viruses, there may be multiple directories - collect all trees with segment keys
    IQTREE_FILE=$(find -L . -name "*.treefile" -print -quit || true)
    CLUSTER_TSV=$(find -L . -name "*_clusters.tsv" -print -quit || true)
    USHER_FILE=$(find -L . -name "final-tree.nh" -print -quit || true)
    if [ -z "$USHER_FILE" ]; then
        USHER_FILE=$(find -L . -name "uncondensed-final-tree.nh" -print -quit || true)
    fi

    TREE_MANIFEST="tree_manifest.tsv"
    printf "source\tname\tsegment_key\tpath\n" > "$TREE_MANIFEST"

    for d in IQTree_MMseqClusters_*; do
        [ -d "$d" ] || continue
        tf=$(find -L "$d" -name "*.treefile" -print -quit || true)
        [ -n "$tf" ] || continue
        seg_key=$(basename "$d" | sed -E 's/^IQTree_MMseqClusters_//' | sed -E 's/_dedup.*$//' | cut -d'.' -f1)
        printf "iqtree\tiqtree_%s\t%s\t%s\n" "$seg_key" "$seg_key" "$tf" >> "$TREE_MANIFEST"
    done

    for d in Usher_MMseqClusters_*; do
        [ -d "$d" ] || continue
        tf=$(find -L "$d" -name "final-tree.nh" -print -quit || true)
        if [ -z "$tf" ]; then
            tf=$(find -L "$d" -name "uncondensed-final-tree.nh" -print -quit || true)
        fi
        [ -n "$tf" ] || continue
        seg_key=$(basename "$d" | sed -E 's/^Usher_MMseqClusters_//' | sed -E 's/_dedup.*$//' | cut -d'.' -f1)
        printf "usher\tusher_%s\t%s\t%s\n" "$seg_key" "$seg_key" "$tf" >> "$TREE_MANIFEST"
    done

    TREE_MANIFEST_ARG=""
    if [ "$(wc -l < "$TREE_MANIFEST")" -gt 1 ]; then
        TREE_MANIFEST_ARG="--tree_manifest $TREE_MANIFEST"
    fi

    echo "Using IQ-TREE file: $IQTREE_FILE ; MMseqs cluster TSV: $CLUSTER_TSV ; USHER file: $USHER_FILE ; tree manifest rows: $(($(wc -l < "$TREE_MANIFEST")-1))"
    IQTREE_ARG=""
    CLUSTER_ARG=""
    USHER_ARG=""
    FILTERED_ARG=""
    FILTERED_DETAILS_ARG=""
    if [ -n "$IQTREE_FILE" ] && [ -f "$IQTREE_FILE" ]; then
        IQTREE_ARG="-it $IQTREE_FILE"
    fi
    if [ -n "$CLUSTER_TSV" ] && [ -f "$CLUSTER_TSV" ]; then
        CLUSTER_ARG="-ct $CLUSTER_TSV -ci !{params.mmseqs_min_seq_id}"
    fi
    if [ -n "$USHER_FILE" ] && [ -f "$USHER_FILE" ]; then
        USHER_ARG="-ut $USHER_FILE"
    fi
    if [ -f "!{filtered_ids}" ] && [ -s "!{filtered_ids}" ]; then
        FILTERED_ARG="-fi !{filtered_ids}"
        echo "Excluding $(wc -l < !{filtered_ids}) filtered sequences from DB"
    fi
    if [ -f "!{filtered_tsv}" ]; then
        FILTERED_DETAILS_ARG="-fd !{filtered_tsv}"
    fi

    python !{scripts_dir}/CreateSqliteDB.py -m !{meta_data} \
    -rf !{features} -p !{sequence_alignment} \
    -i !{insertions} -ht !{host_taxa} \
    -s !{software_info} -fa !{fasta_sequences} \
    -g !{params.gene_info} \
    -mc !{projectDir}/assets/m49_country.csv \
    -mir !{projectDir}/assets/m49_intermediate_region.csv \
    -mr !{projectDir}/assets/m49_region.csv \
    -msr !{projectDir}/assets/m49_sub_region.csv \
    -d !{params.db_name} -b . -o . ${IQTREE_ARG} ${USHER_ARG} ${CLUSTER_ARG} ${FILTERED_ARG} ${FILTERED_DETAILS_ARG} ${TREE_MANIFEST_ARG}
    # need to make gene info come from gff file rather than hardcoded rabv one
    '''
}

process TEST_DB_VALIDATION {
    publishDir "${params.publish_dir}/tests"
    when:
        params.test == "1"
    input:
        path sqlite_db
    output:
        path "db_tree_validation.txt"
        path "db_tree.png"
    shell:
    '''
    EXTRA_ARGS=""
    if [ "!{params.is_segmented}" = "Y" ]; then
        EXTRA_ARGS="--expect-segment-trees"
    fi

    python !{scripts_dir}/ValidateDbTree.py \
        --db !{sqlite_db} \
        --outdir . \
        --test-mode \
        ${EXTRA_ARGS}
    '''
}

process VALIDATE_SEGMENT{

    when:
        params.is_segmented.toBoolean()
    input:
        path gb_matrix
        path blast_hits
    output:
        path "gB_matrix_validated_segment.tsv", emit: validated_matrix
    shell:
    '''
    python !{scripts_dir}/ValidateSegment.py \
        -g !{gb_matrix} \
        -s !{blast_hits} \
        -o gB_matrix_validated_segment.tsv
    '''
}

process VALIDATE_STRAIN{

    when:
        params.is_flu == "Y"
    input:
        path gb_matrix
    output:
        path "gB_matrix_validated_strain.tsv", emit: validated_matrix
    shell:
    '''
    python !{scripts_dir}/ValidateStrain.py \
        -g !{gb_matrix} \
        -o gB_matrix_validated_strain.tsv \
        -m !{projectDir}/generic/influenza/serotype_mapping.tsv
    '''
}

process PIVOT_TABLE_SEGMENTS{
    publishDir "${params.publish_dir}"
    when:
        params.is_segmented =="Y" &&
        params.is_flu == "Y"
    input:
        path gb_matrix
    output:
        path "gB_matrix_pivoted_segments.tsv", emit: pivoted_matrix
    shell:
    '''
    python !{scripts_dir}/FluPivotTable.py \
        -g !{gb_matrix} \
        -o gB_matrix_pivoted_segments.tsv
    '''
}

process TEST_SEGMENTED_OUTPUT{
    publishDir "${params.publish_dir}/tests"
    when:
        params.is_segmented == "Y" && params.test == "1"
    input:
        path annotated_blast
        path validated_matrix
        path pivoted_matrix
        path sqlite_db
    output:
        path "test_segmented_results.txt"
    shell:
    '''
    python !{scripts_dir}/TestPipelineOutput.py \
        --mode segmented \
        --annotated_blast !{annotated_blast} \
        --validated_matrix !{validated_matrix} \
        --pivoted_matrix !{pivoted_matrix} \
        --sqlite_db !{sqlite_db} \
        --output test_segmented_results.txt

    cat test_segmented_results.txt
    '''
}

process TEST_NON_SEGMENTED_OUTPUT{
    publishDir "${params.publish_dir}/tests"
    when:
        params.is_segmented == "N" && params.test == "1"
    input:
        path blast_hits
        path gb_matrix
        path sqlite_db
    output:
        path "test_non_segmented_results.txt"
    shell:
    '''
    python !{scripts_dir}/TestPipelineOutput.py \
        --mode non_segmented \
        --blast_hits !{blast_hits} \
        --gb_matrix !{gb_matrix} \
        --sqlite_db !{sqlite_db} \
        --output test_non_segmented_results.txt

    cat test_non_segmented_results.txt
    '''
}

workflow {

    // check some params are in right form
    // params.is_segmented should be either Y or N
    TEST_DEPENDENCIES()

    def missingParams = []
    if( !params.tax_id ) missingParams << 'tax_id'
    if( !params.db_name ) missingParams << 'db_name'
    if( !params.ref_list ) missingParams << 'ref_list'
    if( !params.gene_info ) missingParams << 'gene_info'
    if( params.mmseqs_min_seq_id == null ) missingParams << 'mmseqs_min_seq_id'
    if( missingParams ){
        error("ERROR: Missing required parameter(s): ${missingParams.join(', ')}")
    }

    // check input params
    if( !(params.is_segmented in ['Y','N']) ){
        error("ERROR: params.is_segmented should be either Y or N")
    }

    if( params.xml_dir ){
        def xmlDirFile = file(params.xml_dir)
        if( !xmlDirFile.exists() || !xmlDirFile.isDirectory() ){
            error("ERROR: params.xml_dir must be an existing directory: ${params.xml_dir}")
        }
    }
    
    if( params.update ){
        def updateFile = file(params.update)
        if( !updateFile.exists() ){
            error("ERROR: params.update file not found: ${params.update}")
        }
        def headerLine = updateFile.text.readLines().find { it?.trim() }
        if( !headerLine ){
            error("ERROR: params.update file is empty: ${params.update}")
        }
        def headerCols = headerLine.split('\t')
        if( !headerCols.contains('primary_accession') ){
            error("ERROR: params.update must be a TSV with header containing 'primary_accession'")
        }
    }
    if( !params.gene_info ){
        error("ERROR: params.gene_info is required. Provide a gene_info TSV with columns: description, display_name, name, parent_name")
    }

    // decide on run mode, update or initial run

    VALIDATE_REF_LIST(params.ref_list, params.is_segmented)
    def ref_list_file = file(params.ref_list)
    def genbank_xml_dir
    
    // Logic: If xml_dir is provided, use it. Otherwise fetch from GenBank.
    if( params.xml_dir ){
        genbank_xml_dir = file(params.xml_dir)
    } else {
        FETCH_GENBANK(params.tax_id, ref_list_file)
        genbank_xml_dir = FETCH_GENBANK.out.gen_bank_XML
    }

    DOWNLOAD_GFF(params.ref_list, ref_list_file)

    GENBANK_PARSER(params.ref_list, genbank_xml_dir)

    def gb_matrix_ch = GENBANK_PARSER.out.gb_matrix
    def gb_seqs_ch = GENBANK_PARSER.out.sequences_out

    if (params.gisaid_dir) {
            def fallbackDbPath = params.scripts_dir ?: "${projectDir}/scripts"
            def db_in = params.previous_db ? file(params.previous_db) : file(fallbackDbPath)
         TIDY_GISAID(params.gisaid_dir, db_in)
         
         // Define column mapping path relative to project
         col_map = file("${projectDir}/generic/influenza/column_mapping.tsv")
         MERGE_GISAID(gb_matrix_ch, TIDY_GISAID.out.gisaid_meta, TIDY_GISAID.out.gisaid_nuc, col_map)
         gb_matrix_ch = MERGE_GISAID.out.merged_matrix
         
         CAT_FASTA(gb_seqs_ch, TIDY_GISAID.out.gisaid_nuc)
         gb_seqs_ch = CAT_FASTA.out.combined_fa
    }

    if(params.extra_info_fill){
        data=ADD_MISSING_DATA(gb_matrix_ch)
    }else{
        data=gb_matrix_ch
    }

    FILTER_AND_EXTRACT(data, 
                        gb_seqs_ch)

    BLAST_ALIGNMENT(FILTER_AND_EXTRACT.out.query_seqs_out,
                    FILTER_AND_EXTRACT.out.ref_seqs_out,
                    data)

    // Add VALIDATE_SEGMENT here
    if (params.is_segmented == 'Y') {
        VALIDATE_SEGMENT(data, BLAST_ALIGNMENT.out.query_uniq_tophit_annotated)
        // Update 'data' to point to the new validated matrix for downstream steps
        data = VALIDATE_SEGMENT.out.validated_matrix
        
        if (params.is_flu == "Y") {
            // Flu pivoting requires Parsed_strain, created by VALIDATE_STRAIN
            VALIDATE_STRAIN(data)
            data = VALIDATE_STRAIN.out.validated_matrix

            PIVOT_TABLE_SEGMENTS(data)
        }
    }


    NEXTALIGN_ALIGNMENT(data,
                        BLAST_ALIGNMENT.out.grouped_fasta,
                        BLAST_ALIGNMENT.out.ref_seqs_dir,
                        BLAST_ALIGNMENT.out.ref_seqs_fasta,
                        BLAST_ALIGNMENT.out.master_seq_dir,
                        params.ref_list,
                        ref_list_file)
    PAD_ALIGNMENT(NEXTALIGN_ALIGNMENT.out,
                  params.ref_list,
                  ref_list_file,
                  (params.ref_set_aligned ?: 'UNSET'))
    
    // Collect sequences that were filtered during nextalign alignment
    COLLECT_FILTERED_SEQUENCES(NEXTALIGN_ALIGNMENT.out)

    // For segmented viruses, PAD_ALIGNMENT emits multiple fasta files (one per segment).
    // Use flatten() to create a channel where each file is processed independently in parallel.
    // For non-segmented viruses, this just passes through the single file.
    padded_msa_ch = PAD_ALIGNMENT.out.merged_msa.flatten()

    DEDUP_ALIGNMENT(padded_msa_ch)

    // Keep clustered input small in test mode to speed up CI and avoid long MMseqs runs
    cluster_input_ch = DEDUP_ALIGNMENT.out.dedup_msa
    if (params.test == "1") {
        TEST_SUBSAMPLE_CLUSTER_INPUT(cluster_input_ch)
        cluster_input_ch = TEST_SUBSAMPLE_CLUSTER_INPUT.out.dedup_for_cluster
    }

    MMSEQS_CLUSTERING(cluster_input_ch)
    IQ_TREE(MMSEQS_CLUSTERING.out.mmseq_clusters)
    
    // Join the channels by segment name for USHER_PLACEMENT
    // Create tuples of (basename, file) for proper matching
    mmseq_with_key = MMSEQS_CLUSTERING.out.mmseq_clusters
        .map { dir ->
            def key = dir.name
                .replaceFirst(/^MMseqClusters_/, '')
                .replaceFirst(/_dedup$/, '')
                .tokenize('.')[0]
            [key, dir]
        }
    iqtree_with_key = IQ_TREE.out.iqtree_out
        .map { dir ->
            def key = dir.name
                .replaceFirst(/^IQTree_MMseqClusters_/, '')
                .replaceFirst(/_dedup$/, '')
                .tokenize('.')[0]
            [key, dir]
        }
    dedup_with_key = cluster_input_ch
        .map { fasta ->
            def key = fasta.name
                .replaceFirst(/_dedup\.fasta$/, '')
                .tokenize('.')[0]
            [key, fasta]
        }
    
    // Join the three channels by their segment key
    usher_input_ch = mmseq_with_key
        .join(iqtree_with_key)
        .join(dedup_with_key)
        .map { key, mmseq_dir, iqtree_dir, dedup_fasta -> 
            tuple(mmseq_dir, iqtree_dir, dedup_fasta)
        }
    
    USHER_PLACEMENT(usher_input_ch)
    
    // VERY_FAST_TREE(PAD_ALIGNMENT.out.merged_msa)
    
    // For CALC_ALIGNMENT_CORD, collect all MSA files back together
    CALC_ALIGNMENT_CORD(PAD_ALIGNMENT.out.merged_msa.collect(), 
                        DOWNLOAD_GFF.out, 
                        BLAST_ALIGNMENT.out.query_uniq_tophits,
                        params.ref_list,
                        ref_list_file)
                        
    SOFTWARE_VERSION()
    
    GENERATE_TABLES(data, 
                    BLAST_ALIGNMENT.out.query_uniq_tophits, 
                    PAD_ALIGNMENT.out.merged_msa.collect(), 
                    NEXTALIGN_ALIGNMENT.out)

    // Collect all per-segment outputs for the database
    // For non-segmented viruses, these will have single items
    iqtree_collected = IQ_TREE.out.iqtree_out.collect()
    mmseq_collected = MMSEQS_CLUSTERING.out.mmseq_clusters.collect()
    usher_collected = USHER_PLACEMENT.out.usher_out.collect()
                    
    CREATE_SQLITE_DB(data, 
                     CALC_ALIGNMENT_CORD.out.features, 
                     GENERATE_TABLES.out.sequence_alignment, 
                     GENERATE_TABLES.out.insertions, 
                     GENERATE_TABLES.out.host_taxa, 
                     SOFTWARE_VERSION.out.software_info, 
                     GENBANK_PARSER.out.sequences_out,
                     iqtree_collected,
                     mmseq_collected,
                     usher_collected,
                     COLLECT_FILTERED_SEQUENCES.out.filtered_tsv,
                     COLLECT_FILTERED_SEQUENCES.out.filtered_ids
                     )

    if (params.test == "1") {
        if (params.is_segmented == "Y" && params.is_flu == "Y") {
            TEST_SEGMENTED_OUTPUT(
                BLAST_ALIGNMENT.out.query_uniq_tophit_annotated,
                VALIDATE_SEGMENT.out.validated_matrix,
                PIVOT_TABLE_SEGMENTS.out.pivoted_matrix,
                CREATE_SQLITE_DB.out
            )
        } else if (params.is_segmented == "N") {
            TEST_NON_SEGMENTED_OUTPUT(
                BLAST_ALIGNMENT.out.query_uniq_tophits,
                data,
                CREATE_SQLITE_DB.out
            )
        }
    }

    TEST_DB_VALIDATION(CREATE_SQLITE_DB.out)
}


// if you wanted it to do an update run, would have to swap "."s for all the directories for a pre-made one

// notes, the python scripts arguments change a lot -d vs -b vs -o etc
// there's too much directory structure, I'd really strip that all out. 
// It could be base=tmp then every function has a subdir in tmp to keep things clear.
