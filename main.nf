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

include { TBBWA } from './modules/mtbseq/tbbwa/tbbwa.nf' addParams(params.TBBWA)
include { TBREFINE } from './modules/mtbseq/tbrefine/tbrefine.nf' addParams(params.TBREFINE)
include { TBPILE } from './modules/mtbseq/tbpile/tbpile.nf' addParams(params.TBPILE)
include { TBLIST } from './modules/mtbseq/tblist/tblist.nf' addParams(params.TBLIST)
include { TBVARIANTS } from './modules/mtbseq/tbvariants/tbvariants.nf' addParams(params.TBVARIANTS)
include { TBSTATS } from './modules/mtbseq/tbstats/tbstats.nf' addParams(params.TBSTATS)
include { TBSTRAINS } from './modules/mtbseq/tbstrains/tbstrains.nf' addParams(params.TBSTRAINS)
include { TBJOIN } from './modules/mtbseq/tbjoin/tbjoin.nf' addParams(params.TBJOIN)
include { TBAMEND } from './modules/mtbseq/tbamend/tbamend.nf' addParams(params.TBAMEND)
include { TBGROUPS } from './modules/mtbseq/tbgroups/tbgroups.nf' addParams(params.TBGROUPS)


workflow test {
//    reads_ch = Channel.fromFilePairs(params.reads)
//    reads_ch = Channel.fromSRA(params.genomeIds, cache: true, apiKey: params.apiKey)


    reads_ch = Channel.of([
        ["5765", ["${baseDir}/../5765_R1.p.fastq.gz", "${baseDir}/../5765_R2.p.fastq.gz"]],
        ["5800", ["${baseDir}/../5800_R1.p.fastq.gz", "${baseDir}/../5800_R2.p.fastq.gz"]]
    ])

    env_user_ch = Channel.value("root")

    // TBBWA(reads_ch,
    //       params.gatk38_jar,
    //       env_user_ch)


    reads_ch.view()

}
