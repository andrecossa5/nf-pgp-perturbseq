// nf-perturbseq-pipeline

nextflow.enable.dsl = 2
include { bulk } from "./subworkflows/bulk/main"
include { sc } from "./subworkflows/sc/main"
include { onerun } from "./subworkflows/onerun/main"

//

// Input bulk
ch_input_bulk = Channel
    .fromPath("${params.bulk_indir}/*", type:'dir') 
    .map{ tuple(it.getName(), it) }

// Input sc_tenx
ch_input_sc_tenx = Channel
    .fromPath("${params.sc_tenx}/*", type:'dir') 
    .map{ tuple(it.getName(), it) }

// Input sc_gbc
ch_input_sc_gbc = Channel
    .fromPath("${params.sc_gbc}/*", type:'dir')
    .map{ tuple(it.getName(), it) }

// Input map bulk-sc
map_bulk_sc = Channel
    .fromPath("${params.bulk_sc_map}")
    .splitCsv(skip:1)
    .map{ tuple(it[0], it[1]) }

//

//----------------------------------------------------------------------------//
// Perturb-seq pipeline entry points and main workflow
//----------------------------------------------------------------------------//

workflow bulk_only {

    bulk(ch_input_bulk)
    bulk.out.flags.view()

}

//

workflow sc_only {

    sc(ch_input_sc_tenx, ch_input_sc_gbc)
    sc.out.summary.view()

}

//

workflow {

    onerun(
        ch_input_bulk, 
        ch_input_sc_tenx, 
        ch_input_sc_gbc, 
        map_bulk_sc
    )

}

//----------------------------------------------------------------------------//
