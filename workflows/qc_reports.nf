include { FASTQC } from "../modules/utils/fastqc.nf" addParams (params.FASTQC)
include { MULTIQC } from "../modules/utils/multiqc.nf" addParams (params.MULTIQC)


workflow QC_REPORTS {
    take:
        reads_ch

    main:
         FASTQC(reads_ch)
         MULTIQC(FASTQC.out.collect())

}
