process FASTQC {
    tag "${genomeName}"
    publishDir params.results_dir, mode: params.save_mode, enabled: params.should_publish

    input:
    tuple val(genomeName), path(genomeReads)

    output:
    tuple path('*.html'), path('*.zip')
    path("*.{html,zip}"), emit: html_zip_ch

    script:

    """
    fastqc *fastq* -t ${task.cpus}
    """

    stub:

    """
    echo "fastqc *fastq*"
    touch ${genomeName}.html
    touch ${genomeName}.zip
    """
}
