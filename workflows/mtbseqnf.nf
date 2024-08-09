/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { QUALITY_CONTROL        } from '../subworkflows/local/quality_control'
include { REPORT                 } from '../subworkflows/local/report'
include { PARALLEL_ANALYSIS      } from "../subworkflows/local/mtbseq-nf-modes/parallel_analysis.nf"
include { NORMAL_ANALYSIS        } from "../subworkflows/local/mtbseq-nf-modes/normal_analysis.nf"
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow MTBSEQ_NF {

    take:
    ch_samplesheet // channel: samplesheet read in from --input

    main:

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()


    QUALITY_CONTROL(ch_samplesheet)

    ch_versions = ch_versions.mix(QUALITY_CONTROL.out.versions)
    ch_multiqc_files = ch_multiqc_files.mix(QUALITY_CONTROL.out.multiqc_files)


    if(!params.only_qc) {
        // MTBSEQ run modes
        if( params.parallel ) {

                PARALLEL_ANALYSIS(ch_samplesheet,
                                  [params.resilist,
                                   params.intregions,
                                   params.categories,
                                   params.basecalib])


                ch_versions.mix(PARALLEL_ANALYSIS.out.versions)
                ch_multiqc_files.mix(PARALLEL_ANALYSIS.out.multiqc_files)

            } else {

                //NOTE: Defaults to the normal analysis as implemented in MTBseq
                NORMAL_ANALYSIS(ch_samplesheet,
                               [params.resilist,
                                params.intregions,
                                params.categories,
                                params.basecalib])

                ch_versions.mix(NORMAL_ANALYSIS.out.versions)
                ch_multiqc_files.mix(NORMAL_ANALYSIS.out.multiqc_files)

        }
    }

    REPORT (ch_multiqc_files.collect(), ch_versions.collect())

    emit:
    multiqc_report = REPORT.out.multiqc_report
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
