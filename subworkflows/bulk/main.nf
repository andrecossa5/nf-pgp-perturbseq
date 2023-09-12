// Bulk workflow

// Include here
nextflow.enable.dsl = 2
include { SEARCH_PATTERNS } from "./modules/generate_search_patterns.nf"
include { EXTRACT_READS } from "./modules/extract_reads.nf"
include { FIND_GBC } from "./modules/find_GBC.nf"
include { REMOVE_SHORT } from "./modules/remove_GBC_shorter_than_18bp.nf"
include { WHITELIST } from "./modules/generate_GBC_whitelist.nf"
include { FORMAT_WHITELIST } from "./modules/reformat_GBC_whitelist_for_error_correction.nf"
include { CORRECT } from "./modules/correct_GBC.nf"
include { COUNT } from "./modules/count_reads_for_each_GBC.nf"
include { INFER_PREVALENCES } from "./modules/infer_clone_prevalences.nf"
include { generate_run_summary_bulk } from "./modules/run_summary.nf"
include { publish_bulk } from "./modules/publish.nf"

 
//


workflow bulk {

  take:
      ch_input

  main:
      SEARCH_PATTERNS()
      EXTRACT_READS(ch_input)
      FIND_GBC(SEARCH_PATTERNS.out.search_patterns, EXTRACT_READS.out.reads)
      REMOVE_SHORT(FIND_GBC.out.GBC)
      WHITELIST(REMOVE_SHORT.out.GBC)
      FORMAT_WHITELIST(WHITELIST.out.whitelist)
      CORRECT(REMOVE_SHORT.out.GBC, FORMAT_WHITELIST.out.formatted_whitelist)
      COUNT(CORRECT.out.GBC, WHITELIST.out.whitelist, REMOVE_SHORT.out.GBC)
      INFER_PREVALENCES(COUNT.out.read_counts)

      // Summary and cleanup
      generate_run_summary_bulk(
        EXTRACT_READS.out.reads, 
        FIND_GBC.out.GBC, 
        REMOVE_SHORT.out.GBC, 
        COUNT.out.read_counts,
        INFER_PREVALENCES.out.stats_table
      )
      publish_bulk(
        CORRECT.out.GBC, 
        COUNT.out.read_counts, 
        WHITELIST.out.whitelist,
        INFER_PREVALENCES.out.stats_table,
        INFER_PREVALENCES.out.prevalences_plot,
        INFER_PREVALENCES.out.spikeins_plot,
        INFER_PREVALENCES.out.df_spikeins,
        generate_run_summary_bulk.out.summary
      )

  emit:
      read_counts = COUNT.out.read_counts
      summary = generate_run_summary_bulk.out.summary

}