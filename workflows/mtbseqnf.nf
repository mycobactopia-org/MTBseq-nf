/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { QC                     } from '../subworkflows/local/qc'
include { REPORT                 } from '../subworkflows/local/report'
include { PARALLEL_ANALYSIS } from "../subworkflows/local/mtbseq-nf-modes/parallel_analysis.nf"
include { NORMAL_ANALYSIS } from "../subworkflows/local/mtbseq-nf-modes/normal_analysis.nf"
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

    QC(ch_samplesheet)
    ch_versions.mix(QC.out.ch_versions)

    // MTBSEQ run modes
    if( params.parallel && !params.only_qc ) {

                PARALLEL_ANALYSIS(ch_samplesheet,
                                  [params.resilist,
                                   params.intregions,
                                   params.categories,
                                   params.basecalib])
           //     ch_versions = ch_versions.mix(PARALLEL_ANALYSIS.out.versions)

            } else {

                //NOTE: Defaults to the normal analysis as implemented in MTBseq
                NORMAL_ANALYSIS(ch_samplesheet,
                               [params.resilist,
                                params.intregions,
                                params.categories,
                                params.basecalib])
           //     ch_versions = ch_versions.mix(NORMAL_ANALYSIS.out.versions)

    }
    /// END MTBSEQ ANALYSIS

    REPORT (QC.out.ch_multiqc_files, ch_versions)
    multiqc_report = REPORT.out.multiqc_report

    emit:
    multiqc_report
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
