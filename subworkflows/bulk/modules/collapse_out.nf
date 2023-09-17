// collapse_output module

nextflow.enable.dsl = 2

//

process collapse_output {

  // Publish
  publishDir "${params.bulk_outdir}/summary", mode: 'copy'

  input:
  val last

  script:
  """
  python3 \
  ${baseDir}/bin/bulk/collapse_outputs.py -i ${params.bulk_outdir}
  """

  stub:
  """
  echo "Collapsing output in ${params.bulk_outdir}..."
  """

}
