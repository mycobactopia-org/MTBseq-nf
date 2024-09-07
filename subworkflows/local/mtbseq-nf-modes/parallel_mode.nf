include { TBBWA } from '../../../modules/mtbseq/tbbwa.nf'
include { TBREFINE } from '../../../modules/mtbseq/tbrefine.nf'
include { TBPILE } from '../../../modules/mtbseq/tbpile.nf'
include { TBLIST } from '../../../modules/mtbseq/tblist.nf'
include { TBVARIANTS } from '../../../modules/mtbseq/tbvariants.nf'
include { TBSTATS } from '../../../modules/mtbseq/tbstats.nf'
include { TBSTRAINS } from '../../../modules/mtbseq/tbstrains.nf'
include { TBJOIN } from '../../../modules/mtbseq/tbjoin.nf'
include { TBAMEND } from '../../../modules/mtbseq/tbamend.nf'
include { TBGROUPS } from '../../../modules/mtbseq/tbgroups.nf'



workflow PARALLEL_MODE {
    take:
        reads_ch
        derived_cohort_tsv
        references_ch

    main:

        ch_versions = Channel.empty()
        ch_multiqc_files = Channel.empty()


        reads_ch.dump(tag:'PARALLEL_MODE.reads_ch')

        TBBWA(reads_ch.collate(2), params.user, references_ch)
        TBREFINE(TBBWA.out.bam_tuple, params.user, references_ch)
        TBPILE(TBREFINE.out.gatk_bam, params.user, references_ch)
        TBLIST(TBPILE.out.mpileup, params.user, references_ch)
        TBVARIANTS(TBLIST.out.position_table_tuple, params.user, references_ch)

    // NOTE: These are part of per-sample-analysis but they need the output from
    // all other processes to compute the metrics for the cohort.
        TBSTATS(TBBWA.out.bam.collect(),
                TBLIST.out.position_table.collect(),
                params.user,
                references_ch)

        TBSTRAINS(TBLIST.out.position_table.collect(),
                  params.user,
                  references_ch)

    // COHORT STEPS


        TBJOIN(TBVARIANTS.out.tbjoin_input.collect(sort:true),
               TBLIST.out.position_table.collect(sort:true),
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

        ch_versions = ch_versions.mix(TBGROUPS.out.versions.first())
        ch_multiqc_files = ch_multiqc_files.mix(TBSTRAINS.out.classification)
                                .mix(TBSTATS.out.statistics)
                                .mix(TBGROUPS.out.distance_matrix.first())
                                .mix(TBGROUPS.out.groups.first())

    emit:
        versions       = ch_versions
        multiqc_files  = ch_multiqc_files.collect()
}
