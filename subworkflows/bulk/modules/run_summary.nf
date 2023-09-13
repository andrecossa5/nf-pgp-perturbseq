// generate_run_summary_bulk module

nextflow.enable.dsl = 2

//

// Process
process generate_run_summary_bulk {

  tag "${sample_name}"

  input:
  tuple val(sample_name), path(reads)
  tuple val(sample_name), path(counts)
  tuple val(sample_name), path(GBC_not_corrected_18bp)
  tuple val(sample_name), path(read_count_by_GBC_corrected)
  tuple val(sample_name), path(stats_table)

        FIND_GBC.out.GBC, 
        REMOVE_SHORT.out.GBC, 
        COUNT.out.read_counts,
        INFER_PREVALENCES.out.stats_table

  output:
  tuple val(sample_name), path('run_summary.txt'), emit: summary

  script:
  """
  echo "Summary Step 1, sample ${sample_name}" > run_summary.txt
  echo "-------------------------------------" >> run_summary.txt
  echo "" >> run_summary.txt
  echo "Overview" >> log.txt
  echo "- Date of analysis:  \$(date)" >> run_summary.txt
  echo "- User:              ${USER}" >> run_summary.txt
  echo "- Working directory: ${PWD}" >> run_summary.txt
  echo "" >> run_summary.txt
  echo "Parameters" >> run_summary.txt
  echo "--indir:        ${params.bulk_indir}" >> run_summary.txt
  echo "--outdir:        ${params.bulk_outdir}" >> run_summary.txt
  echo "--anchor_sequence: ${params.anchor_sequence}" >> run_summary.txt
  echo "" >> run_summary.txt
  echo "Numbers" >> run_summary.txt
  echo "- Reads in input:             \$(zcat ${reads} | wc -l | LC_ALL=en_US.UTF-8 awk '{ printf("%'"'"'d", \$0) }')" >> run_summary.txt
  echo "- Reads with GBC:             \$(zcat ${GBC_not_corrected} | wc -l | LC_ALL=en_US.UTF-8 awk '{ printf("%'"'"'d", \$0) }')" >> run_summary.txt
  echo "- Reads with full length GBC: \$(zcat ${GBC_not_corrected_18bp} | wc -l | LC_ALL=en_US.UTF-8 awk '{ printf("%'"'"'d", \$0) }')" >> run_summary.txt
  echo "- Unique GBCs:" >> run_summary.txt
  echo "  - before error-correction:  \$(zcat ${GBC_not_corrected_18bp} | sort | uniq | wc -l | LC_ALL=en_US.UTF-8 awk '{ printf("%'"'"'d", \$0) }')" >> run_summary.txt
  echo "  - after error-correction (i.e., the used as reference for single-cell): \$(tail -n +2 ${read_count_by_GBC_corrected} | wc -l | LC_ALL=en_US.UTF-8 awk '{ printf("%'"'"'d", \$0) }')" >> run_summary.txt
  echo "- Max number of retained GBC: \$(cat ${stats_table} | wc -l | LC_ALL=en_US.UTF-8 awk '{ printf("%'"'"'d", \$0) }')" >> run_summary.txt
  """

}