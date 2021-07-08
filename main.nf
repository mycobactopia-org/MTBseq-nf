nextflow.enable.dsl = 2
// NOTE: To properly setup the gatk inside the docker image
// - Download the gatk-3.8.0 tar file from here https://console.cloud.google.com/storage/browser/gatk-software/package-archive/gatk;tab=objects?prefix=&forceOnObjectsSortingFiltering=false
// - tar -xvf GATK_TAR_FILE
// - gatk-register gatk_folder/gatk_jar


include { MTBSEQ } from "./modules/mtbseq/mtbseq.nf"

include { TBBWA } from './modules/mtbseq/tbbwa/tbbwa.nf'
include { TBREFINE } from './modules/mtbseq/tbrefine/tbrefine.nf'
include { TBPILE } from './modules/mtbseq/tbpile/tbpile.nf'
include { TBLIST } from './modules/mtbseq/tblist/tblist.nf'
include { TBVARIANTS } from './modules/mtbseq/tbvariants/tbvariants.nf'
include { TBSTATS } from './modules/mtbseq/tbstats/tbstats.nf'
include { TBSTRAINS } from './modules/mtbseq/tbstrains/tbstrains.nf'

workflow mtbseq {
    reads_ch = Channel.fromFilePairs(params.reads)
    gatk38_jar_ch = Channel.value(params.gatk38_jar)
    env_user_ch = Channel.value("root")

    TRIMMOMATIC(reads_ch)
    MTBSEQ(TRIMMOMATIC.out,
            gatk38_jar_ch,
            env_user_ch)

}

workflow per_sample {
    reads_ch = Channel.fromFilePairs(params.reads)

    TBBWA(reads_ch, params.gatk38_jar, params.user)
    TBREFINE(TBBWA.out.next_step, params.gatk38_jar, params.user)
    TBPILE(TBREFINE.out.next_step, params.gatk38_jar, params.user)
    TBLIST(TBPILE.out.next_step, params.gatk38_jar, params.user)
    TBVARIANTS(TBLIST.out.next_step, params.gatk38_jar, params.user)
    TBSTATS(TBBWA.out.next_step,TBLIST.out.next_step, params.gatk38_jar, params.user)
    TBSTRAINS(TBLIST.out.next_step, params.gatk38_jar, params.user)


}
