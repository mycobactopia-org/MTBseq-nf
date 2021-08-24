nextflow.enable.dsl = 2

params.results_dir = "${params.outdir}/utils/rename_files"
params.save_mode = 'copy'
params.should_publish = true

process RENAME_FILES {
    tag "${genomeFileName}"
    publishDir params.results_dir, mode: params.save_mode, enabled: params.should_publish

    input:
    tuple val(genomeFileName), path(reads)

    output:
    path("${genomeFileName}_${params.library_name}_R?.fastq.gz")

    script:
    """
    echo "Renaming files as per MTBseq requirements."

    cp ${reads[0]} ${genomeFileName}_${params.library_name}_R1.fastq.gz
    cp ${reads[1]} ${genomeFileName}_${params.library_name}_R2.fastq.gz
    """
}
