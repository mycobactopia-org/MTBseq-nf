include { PARALLEL_ANALYSIS } from "./workflows/parallel_analysis.nf"
include { NORMAL_ANALYSIS } from "./workflows/normal_analysis.nf"


workflow MTBSEQ_NF {

    //--------------------------------------------
    // Analysis mode
    //--------------------------------------------

    if( params.parallel ) {

        PARALLEL_ANALYSIS(reads_ch,
                          [params.resilist,
                           params.intregions,
                           params.categories,
                           params.basecalib])

    } else {

        //NOTE: Defaults to the normal analysis as implemented in MTBseq
        NORMAL_ANALYSIS(reads_ch,
                       [params.resilist,
                        params.intregions,
                        params.categories,
                        params.basecalib])

    }
}
