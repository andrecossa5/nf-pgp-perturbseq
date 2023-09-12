// generate_run_summary_sc

nextflow.enable.dsl = 2

//

process generate_run_summary_sc {

    tag "${sample_name}"
 
    input:
    tuple val(sample_name), 
        path(R1_tenx), 
        path(R1_gbc),
        path(ref_GBCs),
        path(GBCs), 
        path(filtered), 
        path(cells_summary),
        path(clones_summary)
  
    output:
    tuple val(sample_name), path("run_summary.txt"), emit: summary

    script:
    """
    echo "Summary Step 2, sample ${sample_name}" > run_summary.txt
    echo "-------------------------------------" >> run_summary.txt
    echo "" >> run_summary.txt
    echo "Overview" >> run_summary.txt
    echo "- Date of analysis:  \$(date)" >> run_summary.txt
    echo "- User:              ${USER}" >> run_summary.txt
    echo "- Working directory: ${PWD}" >> run_summary.txt
    echo "" >> run_summary.txt
    echo "Parameters" >> run_summary.txt
    echo "--sc_tenx:                ${params.sc_tenx}" >> run_summary.txt
    echo "--sc_gbc:                ${params.sc_gbc}" >> run_summary.txt
    echo "--sc_outdir:               ${params.sc_outdir}" >> run_summary.txt
    echo "--step_1_out:           ${params.sc_outdir}" >> run_summary.txt
    echo "--pattern:              ${params.sc_pattern}" >> run_summary.txt
    echo "--ref:                  ${params.ref}" >> run_summary.txt
    echo "Numbers" >> run_summary.txt
    echo "- Transcriptomic reads:            \$(zcat ${R1_tenx} | awk 'END{print NR/4}' | LC_ALL=en_US.UTF-8 awk '{ printf("%'"'"'d", \$0) }')" >> run_summary.txt
    echo "- GBC reads:                       \$(zcat ${R1_gbc} | awk 'END{print NR/4}' | LC_ALL=en_US.UTF-8 awk '{ printf("%'"'"'d", \$0) }')" >> run_summary.txt
    echo "- Unique GBC in bulk reference:         \$(cat ${ref_GBCs} | wc -l | LC_ALL=en_US.UTF-8 awk '{ printf("%'"'"'d", \$0) }')" >> run_summary.txt
    echo "- Unique GBC found in sc: \$(cat ${GBCs} | wc -l | LC_ALL=en_US.UTF-8 awk '{ printf("%'"'"'d", \$0) }')" >> run_summary.txt
    echo "- Putative cell n (Solo cell-calling): \$(zcat ${filtered}/barcodes.tsv.gz | wc -l | LC_ALL=en_US.UTF-8 awk '{ printf("%'"'"'d", \$0) }')" >> run_summary.txt
    echo "- Total number of transcripts:     \$(zcat ${filtered}/features.tsv.gz | wc -l | LC_ALL=en_US.UTF-8 awk '{ printf("%'"'"'d", \$0) }')" >> run_summary.txt
    echo "- n clones:     \$(cat ${clones_summary} | wc -l | LC_ALL=en_US.UTF-8 awk '{ printf("%'"'"'d", \$0) }')" >> run_summary.txt
    echo "- n cells confidently assigned to GBC clones: \$(cat ${cells_summary} | wc -l | LC_ALL=en_US.UTF-8 awk '{ printf("%'"'"'d", \$0) }')" >> run_summary.txt
    """

}