nextflow.enable.dsl = 2

params.results_dir = "${params.outdir}/utils/rename_files"
params.save_mode = 'copy'
params.should_publish = true

process RENAME_FILES {
    tag "${genomeFileName}"
    publishDir params.results_dir, mode: params.save_mode, enabled: params.should_publish

    input:
    tuple val(genomeFileName), path("${genomeFileName}_${params.library_name}_R?.fastq.gz")

    output:
    path("*fastq.gz")

    script:
    """
    echo "Renaming files as per MTBseq requirements using Nextflow file staging techniques"
    """
}
