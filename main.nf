nextflow.enable.dsl = 2

include { FASTQC } from "./modules/qc/fastqc.nf" addParams (params.FASTQC)
include { MULTIQC } from "./modules/qc/multiqc.nf" addParams (params.MULTIQC)
include { MTBSEQ_NF } from "./workflows/mtbseq.nf"

//--------------------------------------------
// Source the input reads
//--------------------------------------------


if( params.run_type == "folder" ) {

    // Optional QC workflow
    reads_ch = Channel.fromFilePairs(params.input_folder)

} else {

    // Default to reading from a samplesheet
    reads_ch = Channel.fromPath(params.input_samplesheet)
                        .splitCsv(header: false, skip: 1)
}


workflow QC_REPORTS {
    take:
        reads_ch

    main:
         FASTQC(reads_ch)
         MULTIQC(FASTQC.out.html_zip_ch.collect())

}


workflow {

    MTBSEQ_NF(reads_ch)

}
