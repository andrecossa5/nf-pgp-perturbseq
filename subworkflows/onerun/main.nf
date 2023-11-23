// Onerun workflow

// Include here
nextflow.enable.dsl = 2

// BULK
include { SEARCH_PATTERNS } from "../bulk/modules/generate_search_patterns.nf"
include { EXTRACT_GBC } from "../bulk/modules/extract_GBC.nf"
include { CORRECT_AND_COUNT } from "../bulk/modules/correct_and_count.nf"
include { INFER_PREVALENCES } from "../bulk/modules/infer_clone_prevalences.nf"
include { generate_run_summary_bulk } from "../bulk/modules/run_summary.nf"
include { publish_bulk } from "../bulk/modules/publish.nf"

// SC
include { MERGE_TENX } from "../sc/modules/merge_tenx.nf"
include { MERGE_GBC } from "../sc/modules/merge_gbc.nf"
include { BOWTIE_INDEX_REF } from "../sc/modules/create_bowtie_index_ref.nf"
include { SOLO } from "../sc/modules/Solo.nf"
include { GET_GBC_ELEMENTS } from "../sc/modules/filter_and_extract_from_GBC.nf"
include { GBC_TO_FASTA } from "../sc/modules/gbc_to_fasta.nf"
include { ALIGN_GBC } from "../sc/modules/align_GBC_to_ref.nf"
include { CELL_ASSIGNMENT } from "../sc/modules/cell_assignment.nf"
include { publish_sc } from "../sc/modules/publish.nf"

// ONERUN
include { FASTA_FROM_REF } from "./modules/fasta_from_ref.nf" // NB
include { generate_run_summary_sc } from "./modules/summary.nf"


//


//----------------------------------------------------------------------------//
// Onerun preprocessing subworkflow
//----------------------------------------------------------------------------//

workflow onerun {

    take:
        ch_input_bulk
        ch_tenx
        ch_gbc
        map_bulk_sc

    main:

        // BULK
        SEARCH_PATTERNS()
        EXTRACT_GBC(ch_input, SEARCH_PATTERNS.out.search_patterns)
        CORRECT_AND_COUNT(EXTRACT_GBC.out.GBC)
        INFER_PREVALENCES(CORRECT_AND_COUNT.out.counts)

        generate_run_summary_bulk(
          EXTRACT_READS.out.reads, 
          CORRECT_AND_COUNT.out.read_counts,
          CORRECT_AND_COUNT.out.correction_df,
          INFER_PREVALENCES.out.stats_table
        )
        publish_bulk(
          CORRECT_AND_COUNT.out.counts,
          CORRECT_AND_COUNT.out.correction_df,
          CORRECT_AND_COUNT.out.whitelist,
          INFER_PREVALENCES.out.stats_table,
          INFER_PREVALENCES.out.prevalences_plot,
          INFER_PREVALENCES.out.spikeins_plot,
          INFER_PREVALENCES.out.df_spikeins,
          generate_run_summary_bulk.out.summary
        )
        
        // SC
        MERGE_TENX(ch_tenx)
        MERGE_GBC(ch_gbc)
        SOLO(MERGE_TENX.out.reads)

        // MERGE BULK OUTPUTs 
        combined_read_counts = map_bulk_sc
            .combine(COUNT.out.read_counts, by:0)
            .map{ it -> tuple(it[1], it[2]) }
            .concat(COUNT.out.read_counts)

        // FINISH UP SC
        FASTA_FROM_REF(combined_read_counts)
        BOWTIE_INDEX_REF(FASTA_FROM_REF.out.fasta)
        GET_GBC_ELEMENTS(MERGE_GBC.out.reads.combine(SOLO.out.filtered, by:0))
        GBC_TO_FASTA(GET_GBC_ELEMENTS.out.elements.map{ it -> tuple(it[0], it[3]) })
        ALIGN_GBC(BOWTIE_INDEX_REF.out.index.combine(GBC_TO_FASTA.out.fasta, by:0))
        CELL_ASSIGNMENT(GET_GBC_ELEMENTS.out.elements.combine(ALIGN_GBC.out.names, by:0))

        summary_input = MERGE_TENX.out.reads.map{ it -> tuple(it[0], it[1]) }
            .combine(MERGE_GBC.out.reads.map{ it -> tuple(it[0], it[1]) }, by:0)
            .combine(combined_read_counts, by:0)
            .combine(GET_GBC_ELEMENTS.out.elements.map{ it -> tuple(it[0], it[3]) }, by:0)
            .combine(SOLO.out.filtered, by:0)
            .combine(CELL_ASSIGNMENT.out.cells_summary, by:0)
            .combine(CELL_ASSIGNMENT.out.clones_summary, by:0)
        generate_run_summary_sc(summary_input)

        publish_input = CELL_ASSIGNMENT.out.CBC_GBC_combos
            .combine(CELL_ASSIGNMENT.out.plot, by:0)
            .combine(CELL_ASSIGNMENT.out.cells_summary, by:0)
            .combine(CELL_ASSIGNMENT.out.clones_summary, by:0)
            .combine(SOLO.out.bam, by:0)
            .combine(SOLO.out.stats, by:0)
            .combine(SOLO.out.summary, by:0)
            .combine(SOLO.out.filtered, by:0)
            .combine(SOLO.out.raw, by:0)
            .combine(generate_run_summary_sc.out.summary, by:0)
        publish_sc(publish_input)

    emit:
        summary_bulk = generate_run_summary_bulk.out.summary
        summary_sc = generate_run_summary_sc.out.summary

}
