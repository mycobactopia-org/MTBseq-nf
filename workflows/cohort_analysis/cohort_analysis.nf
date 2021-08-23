nextflow.enable.dsl = 2

include { TBBWA } from '../../modules/mtbseq/tbbwa/tbbwa.nf' addParams (params.TBBWA)
include { TBREFINE } from '../../modules/mtbseq/tbrefine/tbrefine.nf' addParams (params.TBREFINE)
include { TBPILE } from '../../modules/mtbseq/tbpile/tbpile.nf' addParams (params.TBPILE)
include { TBLIST } from '../../modules/mtbseq/tblist/tblist.nf' addParams (params.TBLIST)
include { TBVARIANTS } from '../../modules/mtbseq/tbvariants/tbvariants.nf' addParams (params.TBVARIANTS)
include { TBSTATS } from '../../modules/mtbseq/tbstats/tbstats.nf' addParams (params.TBSTATS)
include { TBSTRAINS } from '../../modules/mtbseq/tbstrains/tbstrains.nf' addParams (params.TBSTRAINS)
include { TBJOIN } from '../../modules/mtbseq/tbjoin/tbjoin.nf' addParams (params.TBJOIN)
include { TBAMEND } from '../../modules/mtbseq/tbamend/tbamend.nf' addParams (params.TBAMEND)
include { TBGROUPS } from '../../modules/mtbseq/tbgroups/tbgroups.nf' addParams (params.TBGROUPS)

workflow COHORT_ANALYSIS {
    take:
        reads_ch

    main:

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
                .collectFile(name: params.samplesheet_name, newLine: false, storeDir: "${params.outdir}")

    //TODO: Consume the emitted channels from per_sample_analysis workflow
        TBJOIN(samples_tsv_file,
               TBVARIANTS.out.tbjoin_input.collect(),
               TBLIST.out.tbjoin_input.collect(),
               params.gatk38_jar,
               params.user)

        TBAMEND(TBJOIN.out.joint_samples, params.gatk38_jar, params.user)
        TBGROUPS(TBAMEND.out.samples_amended, params.gatk38_jar, params.user)

}
