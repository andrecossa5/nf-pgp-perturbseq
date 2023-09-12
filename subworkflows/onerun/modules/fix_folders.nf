// fix_bulk_dir process

nextflow.enable.dsl = 2

//

process fix_bulk_dir {

    input:
    tuple val(bulk), val(sc)

    output:
    stdout

    script:
    """
    python ${baseDir}/bin/fix_bulk_dir.py ${params.bulk_outdir} ${bulk} ${sc}
    """

    stub:
    """
    echo ${bulk} ${sc}
    """

}