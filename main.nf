// nf-perturbseq

nextflow.enable.dsl = 2
include { bulk } from "./subworkflows/bulk/main"
include { sc } from "./subworkflows/sc/main"

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

//

//----------------------------------------------------------------------------//
// Perturb-seq pipeline entry points and main workflow: bulk_is_the_key
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

    Channel.of(1,2,3,4) | view

}

//----------------------------------------------------------------------------//
