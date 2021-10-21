nextflow.enable.dsl = 2
// NOTE: To properly setup the gatk inside the docker image
// - Download the gatk-3.8.0 tar file from here https://console.cloud.google.com/storage/browser/gatk-software/package-archive/gatk;tab=objects?prefix=&forceOnObjectsSortingFiltering=false
// - tar -xvf GATK_TAR_FILE
// - gatk-register gatk_folder/gatk_jar


include { PARALLEL_ANALYSIS } from "./workflows/parallel_analysis.nf"
include { BATCH_ANALYSIS } from "./workflows/batch_analysis.nf"

workflow {

    //NOTE: Create a channel for all reference files
    references_ch = Channel.of[params.global_mtb_ref,
                               params.global_resilist,
                               params.global_intregions,
                               params.global_categories,
                               params.global_basecalib]

    annotations_ch = Channel.of[params.global_resilist,
                                params.global_intregions,
                                params.global_categories,
                                params.global_basecalib]



    //TODO: Refactor this section later to only rely upon a samplesheet.
    // Determine dynamically the samples origin i.e. SRA/

    if( params.run_type == "sra" ) {

        reads_ch = Channel.fromSRA(params.genomeIds, cache: true, apiKey: params.ncbi_api_key)

    } else if( params.run_type == "local" ) {

        reads_ch = Channel.fromPath(params.input_samplesheet)
        .splitCsv(header: false, skip: 1)
        .map { row -> {
                    sampleName           = row[0]
                    read1                = row[1]
                    read2                = row[2]

                    return tuple(sampleName, tuple(file(read1), file(read2)))
                }
            }
    }

    if( params.analysis_mode == "parallel" ) {

        PARALLEL_ANALYSIS(reads_ch,references_ch)

    } else if( params.analysis_mode == "batch" ) {

        BATCH_ANALYSIS(reads_ch,references_ch)

    }
}

//=======================================
// TESTING
//=======================================

workflow TEST {

    reads_ch = Channel.fromPath("${projectDir}/data/mock_data/input_samplesheet.csv")
        .splitCsv(header: false, skip: 1)

    reads_ch.view()

}
