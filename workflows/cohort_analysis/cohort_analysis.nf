nextflow.enable.dsl = 2

include { TBBWA } from '../../modules/mtbseq/tbbwa/tbbwa.nf'
include { TBREFINE } from '../../modules/mtbseq/tbrefine/tbrefine.nf'
include { TBPILE } from '../../modules/mtbseq/tbpile/tbpile.nf'
include { TBLIST } from '../../modules/mtbseq/tblist/tblist.nf'
include { TBVARIANTS } from '../../modules/mtbseq/tbvariants/tbvariants.nf'
include { TBSTATS } from '../../modules/mtbseq/tbstats/tbstats.nf'
include { TBSTRAINS } from '../../modules/mtbseq/tbstrains/tbstrains.nf'
include { TBJOIN } from '../../modules/mtbseq/tbjoin/tbjoin.nf'
include { TBAMEND } from '../../modules/mtbseq/tbamend/tbamend.nf'
include { TBGROUPS } from '../../modules/mtbseq/tbgroups/tbgroups.nf'

workflow COHORT_ANALYSIS {
    reads_ch = Channel.fromFilePairs(params.reads)

    TBBWA(reads_ch, params.gatk38_jar, params.user)
    TBREFINE(TBBWA.out.bam, params.gatk38_jar, params.user)
    TBPILE(TBREFINE.out.gatk_bam, params.gatk38_jar, params.user)
    TBLIST(TBPILE.out.mpileup, params.gatk38_jar, params.user)
    TBVARIANTS(TBLIST.out.position_table, params.gatk38_jar, params.user)
    TBSTATS(TBBWA.out.bam.join(TBLIST.out.position_table), params.gatk38_jar, params.user)
    TBSTRAINS(TBLIST.out.position_table, params.gatk38_jar, params.user)

    samples_tsv_file = TBBWA.out.genomes_names
            .collect()
            .flatten().map { n -> "$n" + "\t" + "${params.library_name}" + "\n" }
            .collectFile(name: 'samples.tsv', newLine: false, storeDir: "${params.outdir}")

    TBJOIN(samples_tsv_file,TBVARIANTS.out.tbjoin_input.collect(), TBLIST.out.tbjoin_input.collect(), params.gatk38_jar, params.user)
    TBAMEND(TBJOIN.out.joint_samples, params.gatk38_jar, params.user)
    TBGROUPS(TBAMEND.out.samples_amended, params.gatk38_jar, params.user)

}
