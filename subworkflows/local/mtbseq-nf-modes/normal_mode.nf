include { TBFULL } from '../../../modules/mtbseq/tbfull.nf' addParams (params.TBFULL)

include { COHORT } from "./cohort_analysis.nf"

workflow NORMAL_MODE {

    take:
        reads_ch
        samples_tsv_file
        references_ch

    main:
        ch_versions = Channel.empty()
        ch_multiqc_files = Channel.empty()


        TBFULL(RENAME_FILES.out.collect(),
               params.user,
               references_ch)


        COHORT(samples_tsv_file,
               TBFULL.out.position_variants,
               TBFULL.out.position_tables,
               references_ch)


        ch_versions = ch_versions.mix(TBFULL.out.versions)
        ch_multiqc_files = ch_multiqc_files.mix(TBFULL.out.statistics)
                                .mix(TBFULL.out.classification)
                                .mix(COHORT.out.distance_matrix)
                                .mix(COHORT.out.groups)


    emit:
        versions       = ch_versions
        multiqc_files  = ch_multiqc_files.collect()


}
