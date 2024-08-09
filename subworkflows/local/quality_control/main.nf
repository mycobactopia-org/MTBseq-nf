include { FASTQC                 } from '../../../modules/nf-core/fastqc/main'




workflow QUALITY_CONTROL {
    take:
    ch_samplesheet

    main:

    FASTQC (ch_samplesheet)

    emit:

    multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]})
    versions = ch_versions.mix(FASTQC.out.versions.first())


}
