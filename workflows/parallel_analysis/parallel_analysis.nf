nextflow.enable.dsl = 2

include { COHORT_ANALYSIS } from "../cohort_analysis/cohort_analysis.nf"
include { PER_SAMPLE_ANALYSIS } from "../per_sample_analysis/per_sample_analysis.nf"

workflow PARALLEL_ANALYSIS {
    take:
        reads_ch

    main:
        refereces_ch = Channel.of([params.ref,params.resilist,params.intregions,params.categories,params.basecalib])
        PER_SAMPLE_ANALYSIS(reads_ch)
        COHORT_ANALYSIS(PER_SAMPLE_ANALYSIS.out.genome_names,
                    PER_SAMPLE_ANALYSIS.out.position_variants,
                    PER_SAMPLE_ANALYSIS.out.position_tables)
}
