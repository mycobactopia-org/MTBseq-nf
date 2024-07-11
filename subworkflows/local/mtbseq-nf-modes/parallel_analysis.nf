include { COHORT_ANALYSIS } from "./cohort_analysis.nf"
include { PER_SAMPLE_ANALYSIS } from "./per_sample_analysis.nf"

workflow PARALLEL_ANALYSIS {
    take:
        reads_ch
        references_ch

    main:
        PER_SAMPLE_ANALYSIS(reads_ch, references_ch)

        COHORT_ANALYSIS(PER_SAMPLE_ANALYSIS.out.genome_names,
                        PER_SAMPLE_ANALYSIS.out.position_variants,
                        PER_SAMPLE_ANALYSIS.out.position_tables,
                        references_ch)
        ch_versions = ch_versions.mix(COHORT_ANALYSIS.out.versions)

    emit:
        versions       = ch_versions
}
