include { FASTQC                 } from '../../modules/nf-core/fastqc/main'




workflow QC {
    take:
    ch_samplesheet

    main:
    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    FASTQC (ch_samplesheet)

    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]})
    ch_versions = ch_versions.mix(FASTQC.out.versions.first())

    emit:
    ch_multiqc_files
    ch_versions
}
