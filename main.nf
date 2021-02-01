nextflow.enable.dsl = 2
// NOTE: To properly setup the gatk inside the docker image
// - Download the gatk-3.8.0 tar file from here https://console.cloud.google.com/storage/browser/gatk-software/package-archive/gatk;tab=objects?prefix=&forceOnObjectsSortingFiltering=false
// - tar -xvf GATK_TAR_FILE
// - gatk-register gatk_folder/gatk_jar


include { MTBSEQ } from "./modules/mtbseq/mtbseq.nf"

workflow mtbseq {
    reads_ch = Channel.fromFilePairs(params.reads)
    gatk38_jar_ch = Channel.value(params.gatk38_jar)
    env_user_ch = Channel.value("root")

    TRIMMOMATIC(reads_ch)
    MTBSEQ(TRIMMOMATIC.out,
            gatk38_jar_ch,
            env_user_ch)

}


