nextflow.enable.dsl = 2

include { TBJOIN } from '../modules/mtbseq/tbjoin.nf' addParams (params.TBJOIN)
include { TBAMEND } from '../modules/mtbseq/tbamend.nf' addParams (params.TBAMEND)
include { TBGROUPS } from '../modules/mtbseq/tbgroups.nf' addParams (params.TBGROUPS)

workflow COHORT_ANALYSIS {
    take:
        genome_names
        position_variants
        position_tables
        references_ch

    main:
        samples_tsv_file = genome_names
                .collect()
                .flatten().map { n -> "$n" + "\t" + "${params.library_name}" + "\n" }
                .collectFile(name: params.cohort_tsv, newLine: false, storeDir: "${params.outdir}", cache: false)

        TBJOIN(position_variants.collect(),
               position_tables.collect(),
               samples_tsv_file,
               params.gatk38_jar,
               params.user,
               references_ch)

        TBAMEND(TBJOIN.out.joint_samples,
                samples_tsv_file,
                params.gatk38_jar,
                params.user,
                references_ch)

        TBGROUPS(TBAMEND.out.samples_amended,
                 samples_tsv_file,
                 params.gatk38_jar,
                 params.user,
                 references_ch)

}
