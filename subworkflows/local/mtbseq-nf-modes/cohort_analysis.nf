include { TBJOIN } from '../../../modules/mtbseq/tbjoin.nf'
include { TBAMEND } from '../../../modules/mtbseq/tbamend.nf'
include { TBGROUPS } from '../../../modules/mtbseq/tbgroups.nf'

workflow COHORT {
    take:
        samples_tsv_file
        position_variants
        position_tables
        references_ch

    main:
        ch_versions = Channel.empty()

        TBJOIN(position_variants.collect(sort:true),
               position_tables.collect(sort:true),
               samples_tsv_file,
               params.user,
               references_ch)

        TBAMEND(TBJOIN.out.joint_samples,
                samples_tsv_file,
                params.user,
                references_ch)

        TBGROUPS(TBAMEND.out.samples_amended,
                 samples_tsv_file,
                 params.user,
                 references_ch)
        ch_versions = ch_versions.mix(TBGROUPS.out.versions.first())



    emit:
        versions         = ch_versions
        groups           = TBGROUPS.out.groups.first()
        distance_matrix  = TBGROUPS.out.distance_matrix.first()
}
