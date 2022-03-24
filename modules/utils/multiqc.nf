nextflow.enable.dsl = 2

process MULTIQC {
    publishDir params.results_dir, mode: params.save_mode, enabled: params.should_publish

    input:
    path("*")

    output:
    tuple path("multiqc_data"), path("multiqc_report.html")

    script:

    """
    multiqc .
    """

    stub:
    """
    mkdir multiqc_data
    touch multiqc_report.html
    """
}
