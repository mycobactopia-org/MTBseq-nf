nextflow.enable.dsl = 2

include { TBFULL } from '../../modules/mtbseq/tbfull/tbfull.nf' addParams (params.TBFULL)
include { TBJOIN } from '../../modules/mtbseq/tbjoin/tbjoin.nf' addParams (params.TBJOIN)
include { TBAMEND } from '../../modules/mtbseq/tbamend/tbamend.nf' addParams (params.TBAMEND)
include { TBGROUPS } from '../../modules/mtbseq/tbgroups/tbgroups.nf' addParams (params.TBGROUPS)

workflow BATCH_ANALYSIS {

    take:
        reads_ch

    main:
        TBFULL(reads_ch,
               params.gatk38_jar,
               params.user)

        samples_tsv_file = TBFULL.out.genome_names
                .collect()
                .flatten().map { n -> "$n" + "\t" + "${params.library_name}" + "\n" }
                .collectFile(name: params.samplesheet_name, newLine: false, storeDir: "${params.outdir}")

        TBJOIN(TBFULL.out.position_variants.collect(),
               TBFULL.out.position_tables.collect(),
               samples_tsv_file,
               params.gatk38_jar,
               params.user)

        TBAMEND(TBJOIN.out.joint_samples,
                samples_tsv_file,
                params.gatk38_jar,
                params.user)

        TBGROUPS(TBAMEND.out.samples_amended,
                 samples_tsv_file,
                 params.gatk38_jar,
                 params.user)



}
