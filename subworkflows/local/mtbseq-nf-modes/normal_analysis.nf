include { RENAME_FILES } from '../../../modules/utils/rename_files.nf' addParams (params.RENAME_FILES)
include { TBFULL } from '../../../modules/mtbseq/tbfull.nf' addParams (params.TBFULL)
include { TBJOIN } from '../../../modules/mtbseq/tbjoin.nf' addParams (params.TBJOIN)
include { TBAMEND } from '../../../modules/mtbseq/tbamend.nf' addParams (params.TBAMEND)
include { TBGROUPS } from '../../../modules/mtbseq/tbgroups.nf' addParams (params.TBGROUPS)

workflow NORMAL_ANALYSIS {

    take:
        reads_ch
        references_ch

    main:
        ch_versions = Channel.empty()
        ch_multiqc_files = Channel.empty()

        samples_tsv_file = reads_ch
                .map {it -> it[0].id}
                .collect()
                .flatten()
                .map { n -> "$n" + "\t" + "${params.library_name}" + "\n" }
                .collectFile(name: params.cohort_tsv, newLine: false, storeDir: "${params.outdir}", cache: false)

        RENAME_FILES(reads_ch)

        //NOTE: Requires atleast 5_CPU/16_MEM
        TBFULL(RENAME_FILES.out.collect(),
               params.user,
               references_ch)

        TBJOIN(TBFULL.out.position_variants.collect(),
               TBFULL.out.position_tables.collect(),
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

        ch_versions = ch_versions.mix(TBGROUPS.out.versions)
        ch_multiqc_files.mix(
            TBFULL.out.statistics,
            TBFULL.out.classification,
            TBGROUPS.out.distance_matrix,
            TBGROUPS.out.groups
            )

    emit:
        versions       = ch_versions
        multiqc_files  = ch_multiqc_files


}
