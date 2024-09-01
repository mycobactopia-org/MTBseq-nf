include { TBFULL } from '../../../modules/mtbseq/tbfull.nf'
include { TBJOIN } from '../../../modules/mtbseq/tbjoin.nf'
include { TBAMEND } from '../../../modules/mtbseq/tbamend.nf'
include { TBGROUPS } from '../../../modules/mtbseq/tbgroups.nf'


workflow NORMAL_MODE {

    take:
        reads_ch
        derived_cohort_tsv
        references_ch

    main:
        ch_versions = Channel.empty()
        ch_multiqc_files = Channel.empty()


        TBFULL(reads_ch.collect(),
               params.user,
               references_ch)


    // COHORT STEPS


        TBJOIN(TBFULL.out.position_variants.collect(sort:true),
               TBFULL.out.position_tables.collect(sort:true),
               derived_cohort_tsv,
               params.user,
               references_ch)

        TBAMEND(TBJOIN.out.joint_samples,
                derived_cohort_tsv,
                params.user,
                references_ch)

        TBGROUPS(TBAMEND.out.samples_amended,
                 derived_cohort_tsv,
                 params.user,
                 references_ch)

        ch_versions = ch_versions.mix(TBFULL.out.versions)
        ch_multiqc_files = ch_multiqc_files.mix(TBFULL.out.statistics)
                                .mix(TBFULL.out.statistics)
                                .mix(TBGROUPS.out.distance_matrix.first())
                                .mix(TBGROUPS.out.groups.first())




    emit:
        versions       = ch_versions
        multiqc_files  = ch_multiqc_files.collect()


}
