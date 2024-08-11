include { FASTQC                 } from '../../../modules/nf-core/fastqc/main'




workflow QUALITY_CONTROL {
    take:
    ch_samplesheet

    main:

    FASTQC (ch_samplesheet)

    emit:

    multiqc_files = FASTQC.out.zip.collect{it[1]}
    versions = FASTQC.out.versions.first()


}
