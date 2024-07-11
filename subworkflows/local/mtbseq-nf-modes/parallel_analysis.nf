include { COHORT_ANALYSIS } from "./cohort_analysis.nf"
include { PER_SAMPLE_ANALYSIS } from "./per_sample_analysis.nf"

workflow PARALLEL_ANALYSIS {
    take:
        reads_ch
        references_ch

    main:
        ch_versions = Channel.empty()
        ch_multiqc_files = Channel.empty()

        PER_SAMPLE_ANALYSIS(reads_ch, references_ch)

        COHORT_ANALYSIS(PER_SAMPLE_ANALYSIS.out.genome_names,
                        PER_SAMPLE_ANALYSIS.out.position_variants,
                        PER_SAMPLE_ANALYSIS.out.position_tables,
                        references_ch)
        ch_versions = ch_versions.mix(COHORT_ANALYSIS.out.versions)
   
    fn: "Strain_Classification.tsv"

    fn: "Mapping_and_Variant_Statistics.tsv"

    fn: "ClusterGroups.tsv"
 
    fn: "distance_matrix.txt"

    emit:
        versions       = ch_versions
        multiqc_files  = ch_multiqc_files
}
