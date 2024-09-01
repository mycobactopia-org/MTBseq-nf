include { FASTQC                 } from '../../../modules/nf-core/fastqc/main'
include { RENAME_FILES           } from '../../../modules/utils/rename_files.nf' addParams (params.RENAME_FILES)



workflow QUALITY_CONTROL {
    take:
        ch_samplesheet

    main:


        samples_tsv_file = ch_samplesheet
                    .map {it -> it[0].id+"\t"+it[0].library}
                    .collectFile(name: params.cohort_tsv, newLine: true, storeDir: "${params.outdir}/misc", cache: false)


        RENAME_FILES (ch_samplesheet)

        RENAME_FILES.out.files.collect().dump(tag: 'RENAME_FILES.out.files')
        RENAME_FILES.out.meta_and_files.collect().dump(tag: 'RENAME_FILES.out.meta_and_files')

        FASTQC (RENAME_FILES.out.meta_and_files)

    emit:
    reads_and_meta_ch      = RENAME_FILES.out.meta_and_files.collect()
    reads_ch               = RENAME_FILES.out.files.collect()
    multiqc_files          = FASTQC.out.zip.collect{it[1]}
    versions               = FASTQC.out.versions.first()
    derived_cohort_tsv     = samples_tsv_file

}
