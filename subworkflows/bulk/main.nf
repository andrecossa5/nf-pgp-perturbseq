// Bulk workflow

// Include here
nextflow.enable.dsl = 2
include { SEARCH_PATTERNS } from "./modules/generate_search_patterns.nf"
include { EXTRACT_READS } from "./modules/extract_reads.nf"
include { FIND_GBC } from "./modules/find_GBC.nf"
include { CORRECT_AND_COUNT } from "./modules/correct_GBC.nf"
include { INFER_PREVALENCES } from "./modules/infer_clone_prevalences.nf"
include { generate_run_summary_bulk } from "./modules/run_summary.nf"
include { publish_bulk } from "./modules/publish.nf"
include { collapse_output } from "./modules/collapse_out.nf"

 
//


workflow bulk {

  take:
      ch_input

  main:
      SEARCH_PATTERNS()
      EXTRACT_READS(ch_input)
      FIND_GBC(SEARCH_PATTERNS.out.search_patterns, EXTRACT_READS.out.reads)
      CORRECT_AND_COUNT(FIND_GBC.out.GBC)
      INFER_PREVALENCES(CORRECT_AND_COUNT.out.counts)

      // Summary and cleanup
      generate_run_summary_bulk(
        EXTRACT_READS.out.reads, 
        CORRECT_AND_COUNT.out.counts,
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
      collapse_output(
        publish_bulk.out.finish_flag.collect().last() // Need to fire only at the end
      )

  emit:
      flags = publish_bulk.out.finish_flag.collect()

}