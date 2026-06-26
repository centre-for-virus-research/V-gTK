
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

def parseCpuSetCount = { String raw ->
    if( !raw ) return null
    int total = 0
    raw.split(',').each { token ->
        def t = token.trim()
        if( !t ) return
        if( t.contains('-') ){
            def parts = t.split('-')
            if( parts.size() == 2 ){
                def start = parts[0] as Integer
                def end = parts[1] as Integer
                if( end >= start ) total += (end - start + 1)
            }
        } else {
            total += 1
        }
    }
    return total > 0 ? total : null
}

def detectEffectiveCpuLimit = {
    def candidates = []

    def runtimeCpus = Runtime.getRuntime().availableProcessors()
    if( runtimeCpus > 0 ) candidates << runtimeCpus

    def cgroupV2 = new File('/sys/fs/cgroup/cpu.max')
    if( cgroupV2.exists() ){
        def txt = cgroupV2.text.trim()
        def parts = txt.split(/\s+/)
        if( parts.size() >= 2 && parts[0] != 'max' ){
            def quota = parts[0] as Long
            def period = parts[1] as Long
            if( quota > 0 && period > 0 ){
                candidates << Math.max(1, (int)Math.floor((double)quota / (double)period))
            }
        }
    }

    def cgroupV1Quota = new File('/sys/fs/cgroup/cpu/cpu.cfs_quota_us')
    def cgroupV1Period = new File('/sys/fs/cgroup/cpu/cpu.cfs_period_us')
    if( cgroupV1Quota.exists() && cgroupV1Period.exists() ){
        def quota = cgroupV1Quota.text.trim() as Long
        def period = cgroupV1Period.text.trim() as Long
        if( quota > 0 && period > 0 ){
            candidates << Math.max(1, (int)Math.floor((double)quota / (double)period))
        }
    }

    def cpusetCandidates = ['/sys/fs/cgroup/cpuset.cpus.effective', '/sys/fs/cgroup/cpuset/cpuset.cpus']
    cpusetCandidates.each { p ->
        def f = new File(p)
        if( f.exists() ){
            def count = parseCpuSetCount(f.text.trim())
            if( count != null && count > 0 ) candidates << count
        }
    }

    return candidates ? candidates.min() : 1
}

def EFFECTIVE_CPU_LIMIT = detectEffectiveCpuLimit()
def maxThreadsRequestedRaw = params.max_threads ?: EFFECTIVE_CPU_LIMIT
def maxThreadsRequested = parsePositiveInt(maxThreadsRequestedRaw, 'max_threads')
def MAX_THREADS = Math.min(maxThreadsRequested, EFFECTIVE_CPU_LIMIT)
if( maxThreadsRequested > EFFECTIVE_CPU_LIMIT ){
    log.warn("params.max_threads=${maxThreadsRequested} exceeds available CPUs (${EFFECTIVE_CPU_LIMIT}); clamping to ${MAX_THREADS}")
}
params.max_threads = MAX_THREADS

def SEGMENT_PARALLEL_THREADS = Math.max(1, (int)Math.floor(MAX_THREADS / 8))
def HEAVY_TASK_CPUS = params.is_segmented == 'Y' ? SEGMENT_PARALLEL_THREADS : MAX_THREADS
def HEAVY_TASK_MAX_FORKS = Math.max(1, (int)Math.floor((double)MAX_THREADS / (double)Math.max(1, HEAVY_TASK_CPUS)))
// 1. List your script's explicitly defined parameters (keep this in sync!)
def scriptDefinedParams = [
    'tax_id', 'db_name', 'is_segmented', 'extra_info_fill', 'test',
    "scripts_dir", "publish_dir", "email", "ref_list", "bulk_fillup_table", "is_flu", "gene_info",
    "xml_dir", "update", "update_file", "update_db",
    "mmseqs_min_seq_id", "mmseqs_trim_cds_file","mutation_catalog", "mutation_virus",
    "gisaid_dir", "previous_db", "conda_path", "test_max_cluster_seqs", "max_threads", "ref_set_aligned",
    "min_seq_length_ratio", "max_aln_gap_proportion", "tree_free"
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

// Backward-compatible aliases: update_file/update -> update_db
if( params.update_file && !params.update_db ){
    params.update_db = params.update_file
}
if( params.update && !params.update_db ){
    params.update_db = params.update
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

process VALIDATE_REF_LIST_DB {
    when:
        params.update_db
    input:
        path ref_list
        path update_db
    output:
        path "ref_list_db_validated.txt"
    shell:
    '''
    python !{scripts_dir}/ValidateRefListAgainstDb.py --ref_list !{ref_list} --db !{update_db}
    echo "ok" > ref_list_db_validated.txt
    '''
}

process TEST_INPUT_DB_VALIDATION {
    publishDir "${params.publish_dir}/tests"
    when:
        params.test == "1" && params.update_db
    input:
        path sqlite_db
    output:
        path "input_db_tree_validation.txt"
        path "input_db_tree.png"
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

    mv db_tree_validation.txt input_db_tree_validation.txt
    mv db_tree.png input_db_tree.png
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

    if [ "!{params.update_db}" != "null" ] && [ -n "!{params.update_db}" ]; then
        extra="${extra} --update !{params.update_db}"
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
    EXTRA_ARGS=""
    if [ -f "!{master_file_opt}" ]; then
        MASTER_ARG="!{master_file_opt}"
    fi
    if [ "!{params.update_db}" != "null" ] && [ -n "!{params.update_db}" ]; then
        EXTRA_ARGS="--update_db !{params.update_db}"
    fi
    python "!{scripts_dir}/DownloadGFF.py" --accession_ids "$MASTER_ARG" -o . -b . ${EXTRA_ARGS}

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
        path "Taxa/names.dmp", emit: taxa_names
        path "Taxa/nodes.dmp", emit: taxa_nodes
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

        if [ "!{params.update_db}" != "null" ] && [ -n "!{params.update_db}" ]; then
            extra="${extra} --update !{params.update_db}"
        fi
        
        echo "DEBUG: Final extra args: ${extra}"
        
            python !{scripts_dir}/GenBankParser.py -r !{ref_list_path} -d !{gen_bank_XML} -o . -b . -s !{params.is_segmented} ${extra} \
            --min_length_ratio !{params.min_seq_length_ratio}
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
        path ref_list_path
    output:
        path "query_seq.fa", emit: query_seqs_out
        path "ref_seq.fa", emit: ref_seqs_out
    shell:
    '''
        EXTRA_ARGS=""
        if [ "!{params.update_db}" != "null" ] && [ -n "!{params.update_db}" ]; then
            EXTRA_ARGS="--update_db !{params.update_db}"
        fi

        python !{scripts_dir}/FilterAndExtractSequences.py -b . -o . -r !{ref_list_path} \
         -v !{params.is_segmented} -g !{table_in} -sf !{seqs_in} ${EXTRA_ARGS}
    '''
}



process BLAST_ALIGNMENT{
    input:
        path query_seqs
        path ref_seqs
        path gb_matrix
        path ref_list_path
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
        UPDATE_BLAST_ARGS=""
        if [ "!{params.update_db}" != "null" ] && [ -n "!{params.update_db}" ]; then
            UPDATE_BLAST_ARGS="--update_db !{params.update_db}"
        fi

        if [ "!{params.is_segmented}" = "Y" ]; then
            python "!{scripts_dir}/BlastAlignment.py" -s Y -f "!{ref_list_path}" -q !{query_seqs} -r !{ref_seqs} \
             -t . -b . -m !{ref_list_path}  -g !{gb_matrix} ${UPDATE_BLAST_ARGS}
        else
            python "!{scripts_dir}/BlastAlignment.py" -f "!{ref_list_path}" -q !{query_seqs} -r !{ref_seqs} \
             -b . -t . -m !{ref_list_path} -g !{gb_matrix} ${UPDATE_BLAST_ARGS}
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
        EXTRA_ARGS=""
        if [ -f "!{master_file_opt}" ]; then
             TARGET_M="!{master_file_opt}"
        fi

        if [ "!{params.update_db}" != "null" ] && [ -n "!{params.update_db}" ]; then
            EXTRA_ARGS="--update_db !{params.update_db}"
        fi

        python !{scripts_dir}/NextalignAlignment.py  -r !{ref_seqs}  \
         -q !{grouped_fasta_dir} -g !{genbank_matrix} -t . \
         -f !{ref_seqs_fasta} -m "$TARGET_M" -ms !{master_seq_dir} ${EXTRA_ARGS}
    '''
}


//"${scripts_dir}/PadAlignment-1.py" -r "/home3/sk312p/task_dir/projects/VGTK/dev_version-jun-09/TING/alUnc509RefseqsMafftHandModified.fa
process PAD_ALIGNMENT{
    publishDir "${params.publish_dir}", mode: 'copy'
    input:
        path nextalign_dir 
        val master_acc_str
        path master_file_opt
        val ref_set_aligned_dir
        path filtered_ids
    output:
        path "*_merged_MSA.fasta", emit: merged_msa
    shell:
    '''
        python !{scripts_dir}/PadAlignment.py -nd !{nextalign_dir} \
        -m !{master_file_opt} \
        -o . -d . -i !{nextalign_dir}/query_aln --keep_intermediate_files \
        --precomputed_ref_dir "!{ref_set_aligned_dir}" --update_db "!{params.update_db}" \
        --segment_manifest_out pad_alignment_manifest.tsv \
        --skip_ids !{filtered_ids}
    '''
}

process COLLECT_FILTERED_SEQUENCES {
    publishDir "${params.publish_dir}" , mode: 'copy'
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
        -b . \
        --max_gap_proportion !{params.max_aln_gap_proportion}
    '''
}

process DEDUP_ALIGNMENT{
    publishDir "${params.publish_dir}" , mode: 'copy'
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
    publishDir "${params.publish_dir}" , mode: 'copy'
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
    publishDir "${params.publish_dir}" , mode: 'copy'
    cpus { params.is_segmented == 'Y' ? SEGMENT_PARALLEL_THREADS : MAX_THREADS }
    maxForks HEAVY_TASK_MAX_FORKS
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
    publishDir "${params.publish_dir}" , mode: 'copy'
    cpus { params.is_segmented == 'Y' ? SEGMENT_PARALLEL_THREADS : MAX_THREADS }
    maxForks HEAVY_TASK_MAX_FORKS
    input:
        tuple path(mmseq_cluster_dir), path(guide_tree_dir)
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

        GUIDE_TREE=$(find -L !{guide_tree_dir} -name "*.nwk" -print -quit)
        if [ -z "$GUIDE_TREE" ]; then
            echo "[error] No guide tree found in !{guide_tree_dir}" >&2
            exit 1
        fi

        mkdir -p IQTree_!{mmseq_cluster_dir.baseName}
        "$IQTREE_BIN" -s "$CLUSTER_REP" -t "$GUIDE_TREE" -T !{task.cpus} -m GTR -pre IQTree_!{mmseq_cluster_dir.baseName}/iqtree -mem 300G
    '''
}

process USHER_PLACEMENT{
    publishDir "${params.publish_dir}" , mode: 'copy'
    cpus { params.is_segmented == 'Y' ? SEGMENT_PARALLEL_THREADS : MAX_THREADS }
    maxForks HEAVY_TASK_MAX_FORKS
    input:
        tuple val(mmseq_cluster_dir), val(iqtree_dir), path(padded_aln), val(usher_output_dir)
    output:
        path "${usher_output_dir}", emit: usher_out
    shell:
    '''
        python !{scripts_dir}/UsherPlacement.py --padded_aln !{padded_aln} \
        --mmseq_cluster_dir "!{mmseq_cluster_dir}" --iqtree_dir "!{iqtree_dir}" \
        --output_dir "!{usher_output_dir}" --update_db "!{params.update_db}" \
        --threads !{task.cpus} --test_mode !{params.test}
    '''
}

process CALC_ALIGNMENT_CORD {
    input:
        path padded_fasta
        path gff_file
        path blast_hits
        path update_scope_tsv
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

    EXTRA_CALC_ARGS=""
    if [ "!{params.update_db}" != "null" ] && [ -n "!{params.update_db}" ]; then
        EXTRA_CALC_ARGS="--update_db !{params.update_db} --update_scope_tsv !{update_scope_tsv} --segment_map_tsv !{update_scope_tsv}"
    fi
    
    python !{scripts_dir}/CalcAlignmentCord.py -i padded_alignments \
    -m "$TARGET_M" -g !{gff_file} -bh !{blast_hits} \
    -b . -d . -o features.tsv ${EXTRA_CALC_ARGS}
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
    maxForks HEAVY_TASK_MAX_FORKS
    when: 
        params.update_db == null
    input:
        path mmseq_cluster_dir
    output:
        path "VeryFastTree_${mmseq_cluster_dir.baseName}", emit: guide_tree_out
    shell:
    '''
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

        mkdir -p VeryFastTree_!{mmseq_cluster_dir.baseName}
        GUIDE_TREE=VeryFastTree_!{mmseq_cluster_dir.baseName}/guide_tree.nwk
        if command -v VeryFastTree >/dev/null 2>&1; then
            VeryFastTree -threads !{task.cpus} -nt -gtr -double-precision "$CLUSTER_REP" > "$GUIDE_TREE"
        elif command -v FastTree >/dev/null 2>&1; then
            FastTree -nt -gtr "$CLUSTER_REP" > "$GUIDE_TREE"
        else
            echo "[error] Neither VeryFastTree nor FastTree was found in PATH. Activate the vgtk conda environment / conda profile." >&2
            exit 1
        fi
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
        path "Tables", emit: tables_dir
    shell:
    '''
    python !{scripts_dir}/GenerateTables.py -g !{gb_matrix} \
    -bh !{blast_hits} -p !{padded_aln} -n !{nextalign_dir} \
    -b . -o Tables -e !{params.email}
    '''
}

process HOST_TAXA_TABLE {
    input:
        path gb_matrix
        path names_dmp
        path nodes_dmp
    output:
        path "HostTaxa/Host_taxa.tsv", emit: host_taxa
    shell:
    '''
    python !{scripts_dir}/HostTaxaTable.py \
        -g !{gb_matrix} \
        -n !{names_dmp} \
        -s !{nodes_dmp} \
        -b . \
        -o HostTaxa
    '''
}

process CREATE_SQLITE_DB {
    publishDir "${params.publish_dir}" , mode: 'copy'
    input:
        path meta_data
        path features
        path sequence_alignment
        path insertions
        path host_taxa
        path software_info
        path fasta_sequences
        path iqtree_dirs, stageAs: 'iqtree_inputs/*'          // Can be multiple directories for segmented viruses
        path mmseq_cluster_dirs, stageAs: 'mmseq_inputs/*'   // Can be multiple directories for segmented viruses
        path usher_dirs, stageAs: 'usher_inputs/*'           // Can be multiple directories for segmented viruses
        path filtered_tsv
        path filtered_ids
    output:
        path "${params.db_name}.db", emit: sqlite_db
        path "db_summary.txt", emit: db_summary_out, optional: true
    shell:
    '''
    # For segmented viruses, there may be multiple directories - collect all trees with segment keys
    IQTREE_FILE=""
    if [ -d iqtree_inputs ]; then
        IQTREE_FILE=$(find -L iqtree_inputs -name "*.treefile" -print -quit || true)
    fi
    CLUSTER_TSV=$(find -L . -name "*_clusters.tsv" -print -quit || true)
    USHER_FILE=""
    if [ -d usher_inputs ]; then
        for USHER_DIR in usher_inputs/*; do
            if [ ! -d "$USHER_DIR" ]; then
                continue
            fi
            if [ -f "$USHER_DIR/final-tree.nh" ]; then
                USHER_FILE="$USHER_DIR/final-tree.nh"
                break
            fi
            if [ -f "$USHER_DIR/uncondensed-final-tree.nh" ]; then
                USHER_FILE="$USHER_DIR/uncondensed-final-tree.nh"
                break
            fi
        done
    fi

    TREE_MANIFEST="tree_manifest.tsv"
    python !{scripts_dir}/BuildTreeManifest.py \
        --output "$TREE_MANIFEST" \
        --iqtree-dir iqtree_inputs \
        --usher-dir usher_inputs

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
    REFERENCE_ARG=""
    if [ "!{params.ref_list}" != "null" ] && [ -n "!{params.ref_list}" ]; then
        if [ ! -f "!{params.ref_list}" ]; then
            echo "[error] reference TSV not found: !{params.ref_list}" >&2
            exit 1
        fi
        REFERENCE_ARG="--reference_tsv !{params.ref_list}"
    fi

    UPDATE_ARGS=""
    if [ "!{params.update_db}" != "null" ] && [ -n "!{params.update_db}" ]; then
        UPDATE_ARGS="--update --update_db !{params.update_db} --batch_id nf_!{workflow.runName}"
    fi

    set -o pipefail

    python !{scripts_dir}/CreateSqliteDB.py -m !{meta_data} \
    -rf !{features} -p !{sequence_alignment} \
    -i !{insertions} -ht !{host_taxa} \
    -s !{software_info} -fa !{fasta_sequences} \
    -g !{params.gene_info} \
    -mc !{projectDir}/assets/m49_country.csv \
    -mir !{projectDir}/assets/m49_intermediate_region.csv \
    -mr !{projectDir}/assets/m49_region.csv \
    -msr !{projectDir}/assets/m49_sub_region.csv \
    -d !{params.db_name} -b . -o . ${IQTREE_ARG} ${USHER_ARG} ${CLUSTER_ARG} ${FILTERED_ARG} ${FILTERED_DETAILS_ARG} ${TREE_MANIFEST_ARG} ${REFERENCE_ARG} ${UPDATE_ARGS} | tee db_summary.txt
    if [ "!{params.mutation_catalog}" != "null" ] && [ -n "!{params.mutation_catalog}" ]; then
        if [ ! -f "!{params.mutation_catalog}" ]; then
            echo "[error] mutation catalog not found: !{params.mutation_catalog}" >&2
            exit 1
        fi

        CATALOG_PROFILE="all_columns"
        VIRUS_ARG=""
        if [ "!{params.mutation_virus}" != "null" ] && [ -n "!{params.mutation_virus}" ]; then
            VIRUS_ARG="--virus !{params.mutation_virus}"
            CATALOG_PROFILE="!{params.mutation_virus}"
        fi

        python !{scripts_dir}/AnnotateMutations.py \
            --db !{params.db_name}.db \
            --mutation_catalog !{params.mutation_catalog} \
            --catalog_column_profile ${CATALOG_PROFILE} \
            ${VIRUS_ARG}
    fi
    # need to make gene info come from gff file rather than hardcoded rabv one

    '''
}

process VERIFY_MUTATIONS {
    publishDir "${params.publish_dir}/tests", mode: 'copy'

    input:
        path sqlite_db
        path mutation_catalog

    output:
        path "mutation_verification.txt", emit: verification_txt

    shell:
    '''
    set -o pipefail

    python !{scripts_dir}/VerifyMutations.py \
        --db !{sqlite_db} \
        --mutation_catalog !{mutation_catalog} 2>&1 | tee mutation_verification.txt
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
    EXTRA_ARGS="${EXTRA_ARGS} --check-update-integrity"
    if [ "!{params.tree_free}" = "true" ]; then
        EXTRA_ARGS="${EXTRA_ARGS} --allow-no-trees"
    fi

    # Merged DB validation: table consistency + tree/update integrity checks
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
    publishDir "${params.publish_dir}" , mode: 'copy'
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
    
    def UPDATE_MODE = params.update_db && params.update_db != 'null'
    def TREE_FREE_MODE = params.tree_free.toString().toBoolean()

    if( UPDATE_MODE ){
        def updateDbFile = file(params.update_db)
        if( !updateDbFile.exists() || !updateDbFile.isFile() ){
            error("ERROR: params.update_db must be an existing SQLite DB file: ${params.update_db}")
        }
        if( !updateDbFile.name.toLowerCase().endsWith('.db') ){
            error("ERROR: params.update_db must point to a .db file: ${params.update_db}")
        }

        def outDbFile = file("${params.publish_dir}/${params.db_name}.db")
        if( updateDbFile.toAbsolutePath().normalize().toString() == outDbFile.toAbsolutePath().normalize().toString() ){
            error("ERROR: update DB path matches output DB path. Use a different --publish_dir or --db_name for update runs.")
        }
    }
    if( !params.gene_info ){
        error("ERROR: params.gene_info is required. Provide a gene_info TSV with columns: description, display_name, name, parent_name")
    }

    // decide on run mode, update or initial run

    if( UPDATE_MODE ){
        TEST_INPUT_DB_VALIDATION(file(params.update_db))
    }

    VALIDATE_REF_LIST(params.ref_list, params.is_segmented)
    def ref_list_file = file(params.ref_list)
    def genbank_xml_dir
    def ref_backbone_dir = (params.ref_set_aligned ?: 'UNSET')
    def effective_ref_list = ref_list_file

    if( UPDATE_MODE ){
        VALIDATE_REF_LIST_DB(ref_list_file, file(params.update_db))
    }
    
    // Logic: If xml_dir is provided, use it. Otherwise fetch from GenBank.
    if( params.xml_dir ){
        genbank_xml_dir = file(params.xml_dir)
    } else {
        FETCH_GENBANK(params.tax_id, effective_ref_list)
        genbank_xml_dir = FETCH_GENBANK.out.gen_bank_XML
    }

    DOWNLOAD_GFF(params.ref_list, effective_ref_list)

    GENBANK_PARSER(effective_ref_list, genbank_xml_dir)

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
                        gb_seqs_ch,
                        effective_ref_list)

    BLAST_ALIGNMENT(FILTER_AND_EXTRACT.out.query_seqs_out,
                    FILTER_AND_EXTRACT.out.ref_seqs_out,
                    data,
                    effective_ref_list)

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
                    effective_ref_list)
    // Collect sequences that were filtered during nextalign alignment (runs first; feeds skip_ids into PAD)
    COLLECT_FILTERED_SEQUENCES(NEXTALIGN_ALIGNMENT.out)

    PAD_ALIGNMENT(NEXTALIGN_ALIGNMENT.out,
                  params.ref_list,
                effective_ref_list,
                  ref_backbone_dir,
                  COLLECT_FILTERED_SEQUENCES.out.filtered_ids)

    // For segmented viruses, PAD_ALIGNMENT emits multiple fasta files (one per segment).
    // Use flatten() to create a channel where each file is processed independently in parallel.
    // For non-segmented viruses, this just passes through the single file.
    padded_msa_ch = PAD_ALIGNMENT.out.merged_msa.flatten()

    DEDUP_ALIGNMENT(padded_msa_ch)

    // Keep clustered input small in test mode to speed up CI and avoid long MMseqs runs
    cluster_input_ch = DEDUP_ALIGNMENT.out.dedup_msa


    def iqtree_collected
    def mmseq_collected
    def usher_collected

    def usher_input_ch

    if( UPDATE_MODE ){
        if( TREE_FREE_MODE ){
            log.warn("params.tree_free=true with params.update_db: skipping tree placement and retaining any existing trees already stored in the update DB")
            iqtree_collected = Channel.value([])
            mmseq_collected = Channel.value([])
            usher_collected = Channel.value([])
        } else {
            usher_input_ch = cluster_input_ch
                .map { fasta ->
                    tuple('UNSET', 'UNSET', fasta, "UsherUpdate_${fasta.baseName}")
                }
        }
    } else {
        MMSEQS_CLUSTERING(cluster_input_ch)
        mmseq_collected = MMSEQS_CLUSTERING.out.mmseq_clusters.collect()

        if( TREE_FREE_MODE ){
            iqtree_collected = Channel.value([])
            usher_collected = Channel.value([])
        } else {
            VERY_FAST_TREE(MMSEQS_CLUSTERING.out.mmseq_clusters)
            mmseq_tree_input_ch = MMSEQS_CLUSTERING.out.mmseq_clusters
                .map { dir ->
                    def key = dir.name
                        .replaceFirst(/^MMseqClusters_/, '')
                        .replaceFirst(/_dedup$/, '')
                        .tokenize('.')[0]
                    [key, dir]
                }
                .join(
                    VERY_FAST_TREE.out.guide_tree_out
                        .map { dir ->
                            def key = dir.name
                                .replaceFirst(/^VeryFastTree_MMseqClusters_/, '')
                                .replaceFirst(/_dedup$/, '')
                                .tokenize('.')[0]
                            [key, dir]
                        }
                )
                .map { key, mmseq_dir, guide_tree_dir ->
                    tuple(mmseq_dir, guide_tree_dir)
                }

            IQ_TREE(mmseq_tree_input_ch)
            
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
                    tuple(mmseq_dir.toString(), iqtree_dir.toString(), dedup_fasta, "Usher_${mmseq_dir.name}")
                }

            iqtree_collected = IQ_TREE.out.iqtree_out.collect()
            USHER_PLACEMENT(usher_input_ch)
            usher_collected = USHER_PLACEMENT.out.usher_out.collect()
        }
    }

    if( UPDATE_MODE && !TREE_FREE_MODE ){
        USHER_PLACEMENT(usher_input_ch)
        usher_collected = USHER_PLACEMENT.out.usher_out.collect()
        iqtree_collected = usher_collected
        mmseq_collected = usher_collected
    }
    
    // For CALC_ALIGNMENT_CORD, collect all MSA files back together
    CALC_ALIGNMENT_CORD(PAD_ALIGNMENT.out.merged_msa.collect(), 
                        DOWNLOAD_GFF.out, 
                        BLAST_ALIGNMENT.out.query_uniq_tophits,
                        data,
                        params.ref_list,
                        effective_ref_list)
                        
    SOFTWARE_VERSION()
    
    GENERATE_TABLES(data, 
                    BLAST_ALIGNMENT.out.query_uniq_tophits, 
                    PAD_ALIGNMENT.out.merged_msa.collect(), 
                    NEXTALIGN_ALIGNMENT.out)

    HOST_TAXA_TABLE(data,
                    GENBANK_PARSER.out.taxa_names,
                    GENBANK_PARSER.out.taxa_nodes)

    CREATE_SQLITE_DB(data, 
                     CALC_ALIGNMENT_CORD.out.features, 
                     GENERATE_TABLES.out.sequence_alignment, 
                     GENERATE_TABLES.out.insertions, 
                     HOST_TAXA_TABLE.out.host_taxa, 
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
                CREATE_SQLITE_DB.out.sqlite_db
            )
        } else if (params.is_segmented == "N") {
            TEST_NON_SEGMENTED_OUTPUT(
                BLAST_ALIGNMENT.out.query_uniq_tophits,
                data,
                CREATE_SQLITE_DB.out.sqlite_db
            )
        }
    }

    TEST_DB_VALIDATION(CREATE_SQLITE_DB.out.sqlite_db)

    if (params.mutation_catalog != null && params.mutation_catalog != "" && params.mutation_catalog != "null") {
        VERIFY_MUTATIONS(CREATE_SQLITE_DB.out.sqlite_db, file(params.mutation_catalog))
    }
}

workflow.onComplete {
    def reportFile = file("${params.publish_dir}/db_summary.txt")
    if (reportFile.exists()) {
        println ""
        println reportFile.text
    } else {
        println "\n[info] No DB summary report found at ${params.publish_dir}/db_summary.txt"
    }
    println "================================================================================"
    println "Pipeline completed at: $workflow.complete"
    println "Execution status     : ${workflow.success ? 'OK' : 'failed'}"
    println "================================================================================\n"
}

// if you wanted it to do an update run, would have to swap "."s for all the directories for a pre-made one

// notes, the python scripts arguments change a lot -d vs -b vs -o etc
// there's too much directory structure, I'd really strip that all out. 
// It could be base=tmp then every function has a subdir in tmp to keep things clear.
