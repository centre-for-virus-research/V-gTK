#!/usr/bin/env nextflow
/*
 * vgtk-reference-update.nf - surgically add new references to a COPY of a database.
 *
 * ******************************************************************************
 * *  WORK IN PROGRESS - NOT EXPECTED TO WORK YET.                              *
 * *  PLAN, EXPORT_BACKBONE and VERIFY_UNCHANGED are implemented. FETCH, ALIGN  *
 * *  and INSERT are sketches (see scripts/AddReferencesToDb.py); INSERT stops  *
 * *  the run with NotImplementedError. Do not use its output.                  *
 * ******************************************************************************
 *
 * Purpose: when a new, hyper-divergent lineage needs a reference, add it with the
 * smallest possible change - no new trees, no clade reassignment, no re-annotation,
 * existing alignments keep their columns - and prove nothing that already existed
 * changed. For new trees and clades, run a full build (vgtk-init.nf) instead.
 *
 * This is deliberately NOT part of vgtk-init.nf: the main pipeline's update mode
 * refuses a reference list that adds references.
 *
 *   nextflow run vgtk-reference-update.nf \
 *       --source_db existing.db --ref_list new_ref_list.txt \
 *       --db_name existing_plus_refs --publish_dir out/ --email you@example.org
 */

nextflow.enable.dsl = 2

params.source_db   = null
params.ref_list    = null
params.db_name     = null
params.publish_dir = "${launchDir}/reference_update_out"
params.email       = null
params.max_threads = 2

def scripts_dir = "${projectDir}/scripts"

log.warn "vgtk-reference-update.nf is WORK IN PROGRESS and not expected to work yet"
[ 'source_db', 'ref_list', 'db_name', 'email' ].each { name ->
    if( !params[name] ) error("ERROR: --${name} is required")
}
if( !file(params.source_db).exists() ) error("ERROR: --source_db not found: ${params.source_db}")
if( !file(params.ref_list).exists() )  error("ERROR: --ref_list not found: ${params.ref_list}")
if( file(params.source_db).name == "${params.db_name}.db" ) {
    error("ERROR: --db_name must differ from the source database; the original is never modified")
}

// Which references are new. Fails unless the list is a pure addition to the
// database: no removed references, no changed master, segment or type, no new master.
process PLAN {
    publishDir "${params.publish_dir}", mode: 'copy'
    input:
        path source_db
        path ref_list
    output:
        path "new_references.tsv"
    shell:
    '''
    python !{scripts_dir}/AddReferencesToDb.py plan --db !{source_db} --ref_list !{ref_list} \
        --output new_references.tsv
    '''
}

// The database's own per-segment backbone, exactly as update mode derives it.
process EXPORT_BACKBONE {
    input:
        path source_db
    output:
        path "db_backbone", type: 'dir'
    shell:
    '''
    python !{scripts_dir}/AddReferencesToDb.py export-backbone --db !{source_db} --output_dir db_backbone
    '''
}

// [SKETCH] Sequences for the new references only.
process FETCH {
    input:
        path plan_tsv
    output:
        path "new_references.fasta"
    shell:
    '''
    python !{scripts_dir}/AddReferencesToDb.py fetch --plan !{plan_tsv} --output new_references.fasta \
        --email "!{params.email}"
    '''
}

// [SKETCH] New references added to the backbone with its columns fixed.
// TODO: one alignment per segment, using the plan's segment column.
process ALIGN {
    cpus params.max_threads
    input:
        path backbone_dir
        path new_fasta
    output:
        path "aligned_new_references.fasta"
    shell:
    '''
    python !{scripts_dir}/AddReferencesToDb.py align --backbone !{backbone_dir}/refset_*_aln.fasta \
        --new_fasta !{new_fasta} --output aligned_new_references.fasta --threads !{task.cpus}
    '''
}

// [SKETCH] Copy of the source database with the new references' rows. Not written yet.
process INSERT {
    input:
        path source_db
        path plan_tsv
        path aligned
    output:
        path "${params.db_name}.db"
    shell:
    '''
    python !{scripts_dir}/AddReferencesToDb.py insert --db !{source_db} --new_db !{params.db_name}.db \
        --plan !{plan_tsv} --aligned !{aligned}
    '''
}

// Every row of every original table must be byte-identical in the copy.
process VERIFY_UNCHANGED {
    publishDir "${params.publish_dir}", mode: 'copy'
    input:
        path source_db, stageAs: 'original/*'
        path new_db
    output:
        path new_db
        path "verify_unchanged.txt"
    shell:
    '''
    set -o pipefail
    python !{scripts_dir}/AddReferencesToDb.py verify-unchanged --original_db !{source_db} \
        --new_db !{new_db} 2>&1 | tee verify_unchanged.txt
    '''
}

workflow {
    source_db = file(params.source_db)
    plan_tsv  = PLAN(source_db, file(params.ref_list))
    backbone  = EXPORT_BACKBONE(source_db)
    new_fasta = FETCH(plan_tsv)
    aligned   = ALIGN(backbone, new_fasta)
    new_db    = INSERT(source_db, plan_tsv, aligned)
    VERIFY_UNCHANGED(source_db, new_db)
}
