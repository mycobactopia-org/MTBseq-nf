/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { QC                     } from '../subworkflows/local/qc'
include { REPORT                 } from '../subworkflows/local/report'
include { PREPARE_MULTIQC } from '../modules/utils/prepare_multiqc'
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
    ch_multiqc_files = Channel.empty()
    QC(ch_samplesheet)
    ch_versions.mix(QC.out.ch_versions)
    ch_multiqc_files.mix(QC.out.ch_multiqc_files)

    // MTBSEQ run modes
    if( params.parallel && !params.only_qc ) {

                PARALLEL_ANALYSIS(ch_samplesheet,
                                  [params.resilist,
                                   params.intregions,
                                   params.categories,
                                   params.basecalib])

               ch_versions = ch_versions.mix(PARALLEL_ANALYSIS.out.versions)
               ch_mtbseq_multiqc_files = PARALLEL_ANALYSIS.out.multiqc_files
            } else {

                //NOTE: Defaults to the normal analysis as implemented in MTBseq
                NORMAL_ANALYSIS(ch_samplesheet,
                               [params.resilist,
                                params.intregions,
                                params.categories,
                                params.basecalib])

                ch_versions = ch_versions.mix(NORMAL_ANALYSIS.out.versions)
                ch_mtbseq_multiqc_files = NORMAL_ANALYSIS.out.multiqc_files

    }
    /// END MTBSEQ ANALYSIS
    PREPARE_MULTIQC(ch_mtbseq_multiqc_files)
    ch_multiqc_files.mix(PREPARE_MULTIQC.out.multiqc_files)

    REPORT (ch_multiqc_files, ch_versions)
    multiqc_report = REPORT.out.multiqc_report

    emit:
    multiqc_report
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
