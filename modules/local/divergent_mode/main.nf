// Processes for the opt-in divergent run mode.
//
// Everything here is reachable only when params.run_divergent_mode is true.
// Nothing in this file is imported by the legacy path, and no process here
// shares a name with one in vgtk-init.nf.
//
// Scripts live in scripts/divergent_mode/. They are invoked by path, never
// staged as `path` inputs: staging moves sys.path[0] and breaks sibling
// imports, the same trap documented for COLLECT_FILTERED_SEQUENCES.
//
// Steps up to the codon alignment are family-wide and run once. Everything
// from gene selection onwards is per segment and carries a `segment` key:
// segments have their own accessions, their own dense clusters and - in a
// family that reassorts - their own topology, so one backbone cannot stand
// for all of them. Per-segment outputs are published under
// `divergent_mode/by_segment/<segment>/`.

process DM_SPLIT_SEGMENTS {
    publishDir "${params.publish_dir}/divergent_mode", mode: 'copy'
    input:
        path genbank_flatfile
    output:
        path "segments", type: 'dir', emit: segments
    shell:
    '''
    python !{params.scripts_dir}/divergent_mode/split_segments.py \
        -i !{genbank_flatfile} \
        -o segments \
        --segments "!{params.divergent_segments}" \
        --floors "!{params.divergent_segment_floors}"
    '''
}

process DM_BUILD_PROTEIN_PROFILES {
    publishDir "${params.publish_dir}/divergent_mode", mode: 'copy'
    cpus { Math.min(params.max_threads as int, 8) }
    input:
        path genbank_flatfile
        path segments_dir
    output:
        path "orthogroups", type: 'dir', emit: orthogroups
    shell:
    '''
    # Reference accessions are the RefSeq entries that survived the
    # completeness floor. Deriving them here rather than hard-coding keeps the
    # step usable on a family whose RefSeq set differs.
    ACC=$(awk -F'\\t' 'NR>1 && $7=="Y" && $1 ~ /^NC_/ {printf "%s,",$1}' \
          !{segments_dir}/segment_manifest.tsv | sed 's/,$//')
    if [ -z "$ACC" ]; then
        echo "[error] no RefSeq reference segments found in !{segments_dir}/segment_manifest.tsv" >&2
        exit 1
    fi
    SEED="!{params.divergent_seed_accessions}"
    if [ "$SEED" = "null" ] || [ -z "$SEED" ]; then SEED="$ACC"; fi

    python !{params.scripts_dir}/divergent_mode/BuildProteinProfiles.py \
        -g !{genbank_flatfile} \
        -a "$ACC" \
        -s "$SEED" \
        -o orthogroups \
        --threads !{task.cpus} \
        --evalue !{params.divergent_cluster_evalue} \
        --min-coverage !{params.divergent_cluster_coverage} \
        --min-identity !{params.divergent_cluster_identity}
    '''
}

process DM_EXTRACT_CDS {
    publishDir "${params.publish_dir}/divergent_mode", mode: 'copy'
    cpus { Math.min(params.max_threads as int, 12) }
    input:
        path segments_dir
        path orthogroup_dir
    output:
        path "cds", type: 'dir', emit: cds
    shell:
    '''
    python !{params.scripts_dir}/divergent_mode/ExtractDivergentCds.py \
        --segment-dir !{segments_dir} \
        --orthogroup-dir !{orthogroup_dir} \
        -o cds \
        --segments "!{params.divergent_segments}" \
        --threads !{task.cpus} \
        --iterations !{params.divergent_extract_iterations}
    '''
}

process DM_CODON_ALIGN {
    publishDir "${params.publish_dir}/divergent_mode", mode: 'copy'
    cpus { Math.min(params.max_threads as int, 12) }
    input:
        path cds_dir
    output:
        path "codon_aln", type: 'dir', emit: codon_aln
    shell:
    '''
    python !{params.scripts_dir}/divergent_mode/CodonAwareAlign.py \
        --cds-dir !{cds_dir} \
        -o codon_aln \
        --threads !{task.cpus} \
        --min-taxa !{params.divergent_min_taxa}

    # Fail loud on any frameshift. The whole point of back-projection is that
    # this count is structurally zero; a non-zero value means the invariant
    # broke, and a tree built on a frameshifted codon alignment is worse than
    # no tree at all.
    BAD=$(awk -F'\\t' 'NR>1 {s+=$5} END {print s+0}' codon_aln/codon_alignment_report.tsv)
    if [ "$BAD" -ne 0 ]; then
        echo "[error] ${BAD} frameshift gaps in codon alignments; back-projection invariant violated" >&2
        exit 1
    fi
    echo "[ok] codon alignments carry 0 frameshift gaps"
    '''
}

process DM_SEED_GENE {
    tag "${segment}"
    publishDir "${params.publish_dir}/divergent_mode/by_segment/${segment}", mode: 'copy'
    input:
        path cds_dir
        path codon_aln_dir
        path segments_dir
        val segment
    output:
        tuple val(segment), path("seed_alignment.fasta"), emit: seed_aln, optional: true
        tuple val(segment), path("${segment}.seed_status.tsv"), emit: status
    shell:
    '''
    python !{params.scripts_dir}/divergent_mode/SelectBackboneGenes.py \
        --codon-dir !{codon_aln_dir} \
        --report !{cds_dir}/cds_extraction_report.tsv \
        -o seed_gene.txt \
        --pick-seed \
        --allow-empty \
        --status !{segment}.seed_status.tsv \
        --manifest !{segments_dir}/segment_manifest.tsv \
        --segments "!{segment}"

    GENE=$(tr -d '\\n' < seed_gene.txt)
    if [ -z "$GENE" ]; then
        # No usable gene on this segment. Not fatal: other segments may still
        # carry one, and DM_SEGMENT_SUMMARY reports the gap at the end.
        echo "[skip] segment !{segment}: no gene to seed a partition" >&2
        exit 0
    fi
    cp "!{codon_aln_dir}/${GENE}.codon.fasta" seed_alignment.fasta
    echo "[ok] segment !{segment} seed gene ${GENE}"
    '''
}

process DM_PARTITION {
    tag "${segment}/${outdir}"
    publishDir "${params.publish_dir}/divergent_mode/by_segment/${segment}/${outdir}", mode: 'copy'
    input:
        tuple val(segment), path(alignment)
        val outdir
    output:
        tuple val(segment), path("partition", type: 'dir'), emit: partition
    shell:
    '''
    python !{params.scripts_dir}/divergent_mode/PartitionByDivergence.py \
        -a !{alignment} \
        -o partition \
        --cutoff "!{params.divergent_cutoff}" \
        --auto-ceiling !{params.divergent_auto_ceiling} \
        --auto-min-band-frac !{params.divergent_auto_min_band_frac} \
        --prefix partition
    '''
}

process DM_SELECT_BACKBONE_GENES {
    tag "${segment}"
    publishDir "${params.publish_dir}/divergent_mode/by_segment/${segment}", mode: 'copy'
    input:
        path cds_dir
        path codon_aln_dir
        path segments_dir
        tuple val(segment), path(partition_dir)
    output:
        tuple val(segment), path("backbone_genes.txt"), emit: genes, optional: true
        tuple val(segment), path("${segment}.backbone_gene_selection.tsv"), emit: summary, optional: true
        tuple val(segment), path("${segment}.selection_status.tsv"), emit: status
    shell:
    '''
    python !{params.scripts_dir}/divergent_mode/SelectBackboneGenes.py \
        --codon-dir !{codon_aln_dir} \
        --report !{cds_dir}/cds_extraction_report.tsv \
        -o backbone_genes.txt \
        --summary !{segment}.backbone_gene_selection.tsv \
        --status !{segment}.selection_status.tsv \
        --allow-empty \
        --segments "!{segment}" \
        --manifest !{segments_dir}/segment_manifest.tsv \
        --assignment !{partition_dir}/partition_assignment.tsv \
        --min-coverage !{params.divergent_min_gene_coverage}

    # An empty list means no gene on this segment clears the coverage bar.
    # Drop the file so the segment leaves the channel instead of emitting a
    # supermatrix of nothing; the summary step reports it.
    if [ ! -s backbone_genes.txt ]; then
        rm -f backbone_genes.txt
        echo "[skip] segment !{segment}: no gene above --divergent_min_gene_coverage" >&2
    fi
    '''
}

process DM_CONCATENATE {
    tag "${segment}"
    publishDir "${params.publish_dir}/divergent_mode/by_segment/${segment}", mode: 'copy'
    input:
        path codon_aln_dir
        tuple val(segment), path(genes_file)
    output:
        tuple val(segment), path("backbone", type: 'dir'), emit: backbone
    shell:
    '''
    GENES=$(tr -d '\\n' < !{genes_file})
    if [ -z "$GENES" ]; then
        echo "[error] segment !{segment}: empty backbone gene list" >&2
        exit 1
    fi
    python !{params.scripts_dir}/divergent_mode/ConcatenateCodonAlignments.py \
        --codon-dir !{codon_aln_dir} \
        --genes "$GENES" \
        -o backbone \
        --prefix backbone
    '''
}

process DM_BACKBONE_ML {
    tag "${segment}"
    publishDir "${params.publish_dir}/divergent_mode/by_segment/${segment}", mode: 'copy'
    cpus { Math.min(params.max_threads as int, 12) }
    input:
        tuple val(segment), path(backbone_dir), path(partition_dir)
    output:
        tuple val(segment), path("backbone_ml", type: 'dir'), emit: backbone_ml
    shell:
    '''
    python !{params.scripts_dir}/divergent_mode/BuildDivergentBackbone.py \
        -a !{backbone_dir}/backbone.codon.fasta \
        --assignment !{partition_dir}/partition_assignment.tsv \
        --coords !{backbone_dir}/backbone.gene_coordinates.tsv \
        -o backbone_ml \
        --model "!{params.divergent_model}" \
        --threads !{task.cpus} \
        --mem "!{params.iqtree_mem}" \
        --bootstrap !{params.divergent_bootstrap} \
        --prefix backbone
    '''
}

process DM_GRAFT {
    tag "${segment}"
    publishDir "${params.publish_dir}/divergent_mode/by_segment/${segment}", mode: 'copy'
    cpus { Math.min(params.max_threads as int, 8) }
    input:
        tuple val(segment), path(backbone_dir), path(partition_dir), path(backbone_ml_dir)
    output:
        tuple val(segment), path("grafted", type: 'dir'), emit: grafted
    shell:
    '''
    python !{params.scripts_dir}/divergent_mode/GraftDenseClusters.py \
        --backbone-tree !{backbone_ml_dir}/backbone.treefile \
        -a !{backbone_dir}/backbone.codon.fasta \
        --assignment !{partition_dir}/partition_assignment.tsv \
        -o grafted \
        --threads !{task.cpus} \
        --prefix grafted
    '''
}

process DM_SEGMENT_SUMMARY {
    publishDir "${params.publish_dir}/divergent_mode", mode: 'copy'
    input:
        path status_files
        val built_segments
    output:
        path "backbone_summary.tsv", emit: summary
    shell:
    '''
    python !{params.scripts_dir}/divergent_mode/SummariseSegmentBackbones.py \
        --status !{status_files} \
        --segments "!{params.divergent_segments}" \
        --built "!{built_segments.join(',')}" \
        -o backbone_summary.tsv
    '''
}
