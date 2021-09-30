nextflow.enable.dsl = 2
// NOTE: To properly setup the gatk inside the docker image
// - Download the gatk-3.8.0 tar file from here https://console.cloud.google.com/storage/browser/gatk-software/package-archive/gatk;tab=objects?prefix=&forceOnObjectsSortingFiltering=false
// - tar -xvf GATK_TAR_FILE
// - gatk-register gatk_folder/gatk_jar


include { PARALLEL_ANALYSIS } from "./workflows/parallel_analysis/parallel_analysis.nf"
include { BATCH_ANALYSIS } from "./workflows/batch_analysis/batch_analysis.nf"
include { TRIMMOMATIC } from "./modules/trimmomatic/trimmomatic.nf"

workflow {
    if( params.run_type == "sra" ) {
        reads_ch = Channel.fromSRA(params.genomeIds, cache: true, apiKey: params.ncbi_api_key)
    } else if( params.run_type == "local" ) {
        reads_ch = Channel.fromFilePairs(params.reads)
    }

    if( params.use_trimmomatic == "false" &&  params.analysis_mode == "parallel") {

        PARALLEL_ANALYSIS(reads_ch)

    } else if ( params.use_trimmomatic == "false" &&  params.analysis_mode == "batch" ) {

        BATCH_ANALYSIS(reads_ch)

    } else if ( params.use_trimmomatic == "true" &&  params.analysis_mode == "parallel" ) {

        TRIMMOMATIC(reads_ch)
        PARALLEL_ANALYSIS(TRIMMOMATIC.out)

    } else if (params.use_trimmomatic == "true" &&  params.analysis_mode == "batch") {

        TRIMMOMATIC(reads_ch)
        BATCH_ANALYSIS(TRIMMOMATIC.out)}
}

//=======================================
// TESTING
//=======================================

workflow test {
    reads_ch = Channel.fromFilePairs("${params.local_location}/*{R1,R2}*gz")

    if( params.analysis_mode == "parallel" ) {
        PARALLEL_ANALYSIS(reads_ch)
    } else if( params.analysis_mode == "batch" ) {
        BATCH_ANALYSIS(reads_ch)
    }

}
