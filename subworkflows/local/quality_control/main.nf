include { FASTQC                 } from '../../../modules/nf-core/fastqc/main'
include { RENAME_FILES           } from '../../../modules/utils/rename_files.nf'



workflow QUALITY_CONTROL {
    take:
        ch_samplesheet

    main:

        samples_tsv_file = Channel.empty()

        RENAME_FILES (ch_samplesheet)

        FASTQC (RENAME_FILES.out.meta_and_files)


        if (!params.cohort_tsv) {

            samples_tsv_file = ch_samplesheet
                        .map {it -> it[0].id+"\t"+it[0].library}
                        .collectFile(name: "derived_cohort.tsv", newLine: true, sort: true, storeDir: "${params.outdir}/misc", cache: true)


        } else {
            samples_tsv_file = Channel.fromPath( params.cohort_tsv )
        }


    emit:
    reads_and_meta_ch      = RENAME_FILES.out.meta_and_files.collect()
    reads_ch               = RENAME_FILES.out.files.collect()
    multiqc_files          = FASTQC.out.zip.collect{it[1]}
    versions               = FASTQC.out.versions.first()
    derived_cohort_tsv     = samples_tsv_file.first()

}
