nextflow.enable.dsl = 2
// NOTE: To properly setup the gatk inside the docker image
// - Download the gatk-3.8.0 tar file from here https://console.cloud.google.com/storage/browser/gatk-software/package-archive/gatk;tab=objects?prefix=&forceOnObjectsSortingFiltering=false
// - tar -xvf GATK_TAR_FILE
// - gatk-register gatk_folder/gatk_jar


include { PARALLEL_ANALYSIS } from "./workflows/parallel_analysis.nf"
include { BATCH_ANALYSIS } from "./workflows/batch_analysis.nf"

workflow {

    //--------------------------------------------
    // Source the input reads
    //
    //TODO: Refactor this section later to only rely upon a samplesheet.
    // Determine dynamically the samples origin i.e. SRA/LOCAL etc
    //--------------------------------------------

    if( params.run_type == "sra" ) {

        reads_ch = Channel.fromSRA(params.genomeIds, cache: true, apiKey: params.ncbi_api_key)

    } else if( params.run_type == "local" ) {

        reads_ch = Channel.fromFilePairs("${params.local_location}/*{R1,R2}*gz")

    } else {

        // Default to reading from a samplesheet

        reads_ch = Channel.fromPath(params.input_samplesheet)
                          .splitCsv(header: false, skip: 1)

    }


    //--------------------------------------------
    // Analysis mode
    //--------------------------------------------

    if( params.analysis_mode == "parallel" ) {

        PARALLEL_ANALYSIS(reads_ch,
                          [params.resilist,
                           params.intregions,
                           params.categories,
                           params.basecalib])

    } else {

        //NOTE: Use the traditional Batch analysis, similar to MTBseq
        BATCH_ANALYSIS(reads_ch,
                       [params.resilist,
                        params.intregions,
                        params.categories,
                        params.basecalib])

    }
}

//=======================================
// TESTING
//=======================================

workflow TEST {


}
