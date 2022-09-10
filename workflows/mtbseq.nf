include { FASTQC } from "../modules/qc/fastqc.nf" addParams (params.FASTQC)
include { MULTIQC } from "../modules/qc/multiqc.nf" addParams (params.MULTIQC)

include { PARALLEL_ANALYSIS } from "./subworkflows/parallel_analysis.nf"
include { NORMAL_ANALYSIS } from "./subworkflows/normal_analysis.nf"



workflow MTBSEQ_NF {
    take: 
        reads_ch


    main:

        //--------------------------------------------
        // Optional QC analysis
        //--------------------------------------------

        if( params.only_qc ) {

             FASTQC(reads_ch)
             MULTIQC(FASTQC.out.html_zip_ch.collect())

        }

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
