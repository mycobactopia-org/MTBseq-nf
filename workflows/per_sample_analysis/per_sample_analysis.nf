nextflow.enable.dsl = 2

include { TBBWA } from '../../modules/mtbseq/tbbwa/tbbwa.nf' addParams (params.TBBWA)
include { TBREFINE } from '../../modules/mtbseq/tbrefine/tbrefine.nf' addParams (params.TBREFINE)
include { TBPILE } from '../../modules/mtbseq/tbpile/tbpile.nf' addParams (params.TBPILE)
include { TBLIST } from '../../modules/mtbseq/tblist/tblist.nf' addParams (params.TBLIST)
include { TBVARIANTS } from '../../modules/mtbseq/tbvariants/tbvariants.nf' addParams (params.TBVARIANTS)
include { TBSTATS } from '../../modules/mtbseq/tbstats/tbstats.nf' addParams (params.TBSTATS)
include { TBSTRAINS } from '../../modules/mtbseq/tbstrains/tbstrains.nf' addParams (params.TBSTRAINS)

workflow PER_SAMPLE_ANALYSIS {
    reads_ch = Channel.fromFilePairs(params.reads)

    TBBWA(reads_ch, params.gatk38_jar, params.user)
    TBREFINE(TBBWA.out.bam, params.gatk38_jar, params.user)
    TBPILE(TBREFINE.out.gatk_bam, params.gatk38_jar, params.user)
    TBLIST(TBPILE.out.mpileup, params.gatk38_jar, params.user)
    TBVARIANTS(TBLIST.out.position_table, params.gatk38_jar, params.user)
    TBSTATS(TBBWA.out.bam.join(TBLIST.out.position_table), params.gatk38_jar, params.user)
    TBSTRAINS(TBLIST.out.position_table, params.gatk38_jar, params.user)

}

