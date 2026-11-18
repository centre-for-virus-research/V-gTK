// Opt-in divergent run mode for structurally divergent, multi-segment
// families (first testbed Cystoviridae, NCBI TaxID 10877; the per-segment
// path is exercised by Hantaviridae, TaxID 1980413).
//
// Reachable only from `if (params.run_divergent_mode)` in vgtk-init.nf. When
// the flag is false this file is never imported, no process here is
// instantiated, and the legacy path is untouched.
//
// Why this is a separate path rather than a set of swapped-in steps: the
// legacy path resolves every coordinate against one master accession, and the
// whole point here is that no such master exists for this family. Bolting the
// divergent steps into the middle of that path would mean changing the
// coordinate contract that CreateSqliteDB, AnnotateMutations and the
// insertions table all depend on. This mode therefore produces alignments, a
// tree and parent-child pairs as published artefacts, and stops short of the
// database build. Wiring it into CreateSqliteDB is deliberately left as
// separate work with its own schema review.
//
// Shape: the extraction and alignment steps are family-wide and run once.
// From gene selection onwards the workflow fans out over segments, because a
// segment is its own analysis - its own accessions, its own dense clusters,
// and in a reassorting family its own topology. Each segment that carries a
// gene recoverable across the divergent taxa gets its own backbone; each that
// does not is named in `backbone_summary.tsv` with the reason.

include {
    DM_SPLIT_SEGMENTS;
    DM_BUILD_PROTEIN_PROFILES;
    DM_EXTRACT_CDS;
    DM_CODON_ALIGN;
    DM_SEED_GENE;
    DM_PARTITION as DM_PARTITION_SEED;
    DM_PARTITION as DM_PARTITION_FINAL;
    DM_SELECT_BACKBONE_GENES;
    DM_CONCATENATE;
    DM_BACKBONE_ML;
    DM_GRAFT;
    DM_SEGMENT_SUMMARY;
} from '../../../modules/local/divergent_mode/main.nf'

workflow DIVERGENT_MODE {
    take:
        genbank_flatfile

    main:
        // 1. Segment split, with per-segment completeness floors.
        DM_SPLIT_SEGMENTS(genbank_flatfile)

        // 2. Orthologue groups from the annotated references, used as the
        //    tblastn query set. Product strings are not reliable gene keys,
        //    so groups come from sequence similarity and only their names
        //    come from annotation. Groups never span segments.
        DM_BUILD_PROTEIN_PROFILES(genbank_flatfile, DM_SPLIT_SEGMENTS.out.segments)

        // 3. CDS extraction, resolved to whole codons by construction.
        DM_EXTRACT_CDS(DM_SPLIT_SEGMENTS.out.segments,
                       DM_BUILD_PROTEIN_PROFILES.out.orthogroups)

        // 4. Protein-guided alignment and codon back-projection. This step
        //    aborts if any frameshift gap appears in its own output.
        DM_CODON_ALIGN(DM_EXTRACT_CDS.out.cds)

        // The two directories below are consumed once per segment, so they
        // have to be value channels; a queue channel would be drained by the
        // first segment and the rest would never run.
        cds_dir   = DM_EXTRACT_CDS.out.cds.first()
        codon_dir = DM_CODON_ALIGN.out.codon_aln.first()
        seg_dir   = DM_SPLIT_SEGMENTS.out.segments.first()

        segments_ch = Channel.fromList(
            params.divergent_segments.toString()
                  .split(',').collect{ it.trim() }.findAll{ it } as List
        )

        // 5. Two-stage gene choice, per segment. A gene looks universal if
        //    scored over all taxa when one dense cluster dominates the
        //    dataset, so each segment is first partitioned on its own most
        //    complete gene, and real coverage is then measured over that
        //    segment's cluster representatives. Scoping one segment's
        //    selection to another segment's partition scores every gene at
        //    zero; SelectBackboneGenes.py refuses that outright.
        DM_SEED_GENE(cds_dir, codon_dir, seg_dir, segments_ch)
        DM_PARTITION_SEED(DM_SEED_GENE.out.seed_aln, 'seed')

        DM_SELECT_BACKBONE_GENES(cds_dir, codon_dir, seg_dir,
                                 DM_PARTITION_SEED.out.partition)

        // 6. Partitioned supermatrix over the genes chosen for that segment.
        //    A segment with no qualifying gene emits nothing here and leaves
        //    the channel; it is reported in step 8 rather than dropped.
        DM_CONCATENATE(codon_dir, DM_SELECT_BACKBONE_GENES.out.genes)

        // 7. Final partition on the supermatrix, then ML backbone + ASR over
        //    representatives only, then parsimony placement inside each dense
        //    cluster and grafting back onto that segment's backbone.
        DM_PARTITION_FINAL(
            DM_CONCATENATE.out.backbone.map { seg, dir ->
                tuple(seg, "${dir}/backbone.codon.fasta")
            },
            'final')

        DM_BACKBONE_ML(
            DM_CONCATENATE.out.backbone.join(DM_PARTITION_FINAL.out.partition))

        DM_GRAFT(
            DM_CONCATENATE.out.backbone
                .join(DM_PARTITION_FINAL.out.partition)
                .join(DM_BACKBONE_ML.out.backbone_ml))

        // 8. One row per requested segment: backbone built, or the reason it
        //    was not. Fails the run if no segment produced a backbone.
        status_files = DM_SEED_GENE.out.status.map { seg, f -> f }
            .mix(DM_SELECT_BACKBONE_GENES.out.status.map { seg, f -> f })
            .collect()
        built = DM_GRAFT.out.grafted.map { seg, dir -> seg }
            .collect().ifEmpty([])

        DM_SEGMENT_SUMMARY(status_files, built)

    emit:
        segments    = DM_SPLIT_SEGMENTS.out.segments
        codon_aln   = DM_CODON_ALIGN.out.codon_aln
        backbone    = DM_CONCATENATE.out.backbone
        partition   = DM_PARTITION_FINAL.out.partition
        backbone_ml = DM_BACKBONE_ML.out.backbone_ml
        grafted     = DM_GRAFT.out.grafted
        summary     = DM_SEGMENT_SUMMARY.out.summary
}
