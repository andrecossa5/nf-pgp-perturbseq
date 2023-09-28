// sc subworkflow

nextflow.enable.dsl = 2

// Include here
include { MERGE_TENX } from "./modules/merge_tenx.nf"
include { MERGE_GBC } from "./modules/merge_gbc.nf"
include { FASTA_FROM_REF } from "./modules/fasta_from_ref.nf"
include { BOWTIE_INDEX_REF } from "./modules/create_bowtie_index_ref.nf"
include { SOLO } from "./modules/Solo.nf"
include { GET_GBC_ELEMENTS } from "./modules/filter_and_extract_from_GBC.nf"
include { GBC_TO_FASTA } from "./modules/gbc_to_fasta.nf"
include { ALIGN_GBC } from "./modules/align_GBC_to_ref.nf"
include { CELL_ASSIGNMENT } from "./modules/cell_assignment.nf"
include { generate_run_summary_sc } from "./modules/run_summary.nf"
include { publish_sc } from "./modules/publish.nf"

 
//


workflow sc {
    
    take:
        ch_tenx
        ch_gbc

    main:
 
        // Merge reads
        MERGE_TENX(ch_tenx)
        MERGE_GBC(ch_gbc)

        // STARSolo
        SOLO(MERGE_TENX.out.reads)
 
        // Assign cells to clones
        FASTA_FROM_REF(ch_tenx)
        BOWTIE_INDEX_REF(FASTA_FROM_REF.out.fasta)
        // GET_GBC_ELEMENTS(MERGE_GBC.out.reads.combine(SOLO.out.filtered, by:0))
        // GBC_TO_FASTA(GET_GBC_ELEMENTS.out.elements.map{ it -> tuple(it[0], it[3]) })
        // ALIGN_GBC(BOWTIE_INDEX_REF.out.index.combine(GBC_TO_FASTA.out.fasta, by:0))
        // CELL_ASSIGNMENT(GET_GBC_ELEMENTS.out.elements.combine(ALIGN_GBC.out.names, by:0))

        // Summary
        // summary_input = MERGE_TENX.out.reads.map{ it -> tuple(it[0], it[1]) }
        //     .combine(MERGE_GBC.out.reads.map{ it -> tuple(it[0], it[1]) }, by:0)
        //     .combine(GET_GBC_ELEMENTS.out.elements.map{ it -> tuple(it[0], it[3]) }, by:0)
        //     .combine(SOLO.out.filtered, by:0)
        //     .combine(CELL_ASSIGNMENT.out.cells_summary, by:0)
        //     .combine(CELL_ASSIGNMENT.out.clones_summary, by:0)
        // generate_run_summary_sc(summary_input)

        // Publishing
        // publish_input = CELL_ASSIGNMENT.out.CBC_GBC_combos
        //     .combine(CELL_ASSIGNMENT.out.plot, by:0)
        //     .combine(CELL_ASSIGNMENT.out.cells_summary, by:0)
        //     .combine(CELL_ASSIGNMENT.out.clones_summary, by:0)
        //     .combine(SOLO.out.bam, by:0)
        //     .combine(SOLO.out.stats, by:0)
        //     .combine(SOLO.out.summary, by:0)
        //     .combine(SOLO.out.filtered, by:0)
        //     .combine(SOLO.out.raw, by:0)
        //     .combine(generate_run_summary_sc.out.summary, by:0)
        // publish_sc(publish_input)

    emit:
        // summary = generate_run_summary_sc.out.summary
        summary = FASTA_FROM_REF.out.fasta

}