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

    QC(ch_samplesheet)

    // MTBSEQ run modes
    if( params.parallel && !params.only_qc ) {

                PARALLEL_ANALYSIS(ch_samplesheet,
                                  [params.resilist,
                                   params.intregions,
                                   params.categories,
                                   params.basecalib])


                //ch_versions = Channel.empty()
                //ch_multiqc_files = Channel.empty()

                //ch_versions.mix(QC.out.ch_versions)
                //ch_multiqc_files.mix(QC.out.ch_multiqc_files)

                ch_versions = Channel.empty().mix(PARALLEL_ANALYSIS.out.versions)
                ch_multiqc_files = Channel.empty().mix(PARALLEL_ANALYSIS.out.multiqc_files)

                ch_multiqc_files.view()

                REPORT (ch_multiqc_files.collect(), ch_versions.collect())
            } else {

                //NOTE: Defaults to the normal analysis as implemented in MTBseq
                NORMAL_ANALYSIS(ch_samplesheet,
                               [params.resilist,
                                params.intregions,
                                params.categories,
                                params.basecalib])

                ch_versions = Channel.empty()
                ch_multiqc_files = Channel.empty()

                ch_versions.mix(QC.out.ch_versions)
                ch_multiqc_files.mix(QC.out.ch_multiqc_files)


                ch_versions.mix(NORMAL_ANALYSIS.out.versions)
                ch_multiqc_files.mix(NORMAL_ANALYSIS.out.multiqc_files)

                REPORT (ch_multiqc_files.collect(), ch_versions.collect())
    }



    emit:
    multiqc_report = REPORT.out.multiqc_report
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
