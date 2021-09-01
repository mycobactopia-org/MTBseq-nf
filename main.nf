nextflow.enable.dsl = 2
// NOTE: To properly setup the gatk inside the docker image
// - Download the gatk-3.8.0 tar file from here https://console.cloud.google.com/storage/browser/gatk-software/package-archive/gatk;tab=objects?prefix=&forceOnObjectsSortingFiltering=false
// - tar -xvf GATK_TAR_FILE
// - gatk-register gatk_folder/gatk_jar


include { PER_SAMPLE_ANALYSIS } from "./workflows/per_sample_analysis/per_sample_analysis.nf"
include { COHORT_ANALYSIS } from "./workflows/cohort_analysis/cohort_analysis.nf"
include { BATCH_ANALYSIS } from "./workflows/batch_analysis/batch_analysis.nf"

workflow {
    reads_ch = Channel.fromSRA(params.genomeIds, cache: true, apiKey: params.apiKey)
    // reads_ch = Channel.fromFilePairs("${params.local_location}/*{R1,R2}*gz")


    //NOTE: Parallel Analysis
    PER_SAMPLE_ANALYSIS(reads_ch)
    COHORT_ANALYSIS(PER_SAMPLE_ANALYSIS.out.genome_names,
                    PER_SAMPLE_ANALYSIS.out.position_variants,
                    PER_SAMPLE_ANALYSIS.out.position_tables)
}

//=======================================
// TESTING
//=======================================

workflow test {
    reads_ch = Channel.fromFilePairs("${params.local_location}/*{R1,R2}*gz")

    reads_ch
        .map { it -> it[1] }
        .collect()
        .flatten()
        .view()

    // BATCH_ANALYSIS(reads_ch)

    // PER_SAMPLE_ANALYSIS(reads_ch)

    // COHORT_ANALYSIS(PER_SAMPLE_ANALYSIS.out.genome_names,
    //                 PER_SAMPLE_ANALYSIS.out.position_variants,
    //                 PER_SAMPLE_ANALYSIS.out.position_tables)
}
