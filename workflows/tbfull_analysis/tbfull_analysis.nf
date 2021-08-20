nextflow.enable.dsl = 2

include { TBFULL } from '../../modules/mtbseq/tbfull/tbfull.nf'
include { TBJOIN } from '../../modules/mtbseq/tbjoin/tbjoin.nf'
include { TBAMEND } from '../../modules/mtbseq/tbamend/tbamend.nf'
include { TBGROUPS } from '../../modules/mtbseq/tbgroups/tbgroups.nf'

workflow TBFULL_ANALYSIS {
    reads_ch = Channel.fromFilePairs(params.reads)

    TBFULL(reads_ch, params.gatk38_jar, params.user)
    samples_tsv_file = TBFULL.out.genomes_names
            .collect()
            .flatten().map { n -> "$n" + "\t" + "${params.library_name}" + "\n" }
            .collectFile(name: 'samples.tsv', newLine: false, storeDir: "${params.outdir}")

    TBJOIN(samples_tsv_file,TBFULL.out.variants_table.collect(), TBLIST.out.position_table.collect(), params.gatk38_jar, params.user)
    TBAMEND(TBJOIN.out.joint_samples, params.gatk38_jar, params.user)
    TBGROUPS(TBAMEND.out.samples_amended, params.gatk38_jar, params.user)

}
