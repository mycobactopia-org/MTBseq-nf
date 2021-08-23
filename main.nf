nextflow.enable.dsl = 2
// NOTE: To properly setup the gatk inside the docker image
// - Download the gatk-3.8.0 tar file from here https://console.cloud.google.com/storage/browser/gatk-software/package-archive/gatk;tab=objects?prefix=&forceOnObjectsSortingFiltering=false
// - tar -xvf GATK_TAR_FILE
// - gatk-register gatk_folder/gatk_jar


include { PER_SAMPLE_ANALYSIS } from "./workflows/per_sample_analysis/per_sample_analysis.nf"

include { COHORT_ANALYSIS } from "./workflows/cohort_analysis/cohort_analysis.nf"

include { TBFULL_ANALYSIS } from "./workflows/tbfull_analysis/tbfull_analysis.nf"

workflow {

    COHORT_ANALYSIS()

}

//=======================================
// TESTING
//=======================================

include { TBFULL } from './modules/mtbseq/tbfull/tbfull.nf' addParams(params.TBFULL)

include { TBJOIN } from './modules/mtbseq/tbjoin/tbjoin.nf' addParams(params.TBJOIN)
include { TBAMEND } from './modules/mtbseq/tbamend/tbamend.nf' addParams(params.TBAMEND)
include { TBGROUPS } from './modules/mtbseq/tbgroups/tbgroups.nf' addParams(params.TBGROUPS)


workflow test {
    reads_ch = Channel.fromFilePairs("${params.local_location}/*{R1,R2}*gz")
    // reads_ch = Channel.fromSRA(params.genomeIds, cache: true, apiKey: params.apiKey)


    TBFULL(reads_ch,
           params.gatk38_jar,
           params.user)

    // PER_SAMPLE_ANALYSIS(reads_ch) // DONE

}
