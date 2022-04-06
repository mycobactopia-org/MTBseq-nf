nextflow.enable.dsl = 2



params.gatk_jar_link = "https://storage.googleapis.com/gatk-software/package-archive/gatk/GenomeAnalysisTK-3.8-0-ge9d806836.tar.bz2"

process DOWNLOAD_GATK_JAR {
    publishDir params.results_dir, mode: params.save_mode, enabled: params.should_publish

    output:
    path('**/GenomeAnalysisTK.jar')


    script:

    """
    wget ${params.gatk_jar_link}
    tar -xfv *.tar.bz2
    """

    stub:

    """
    echo "wget ${params.gatk_jar_link}"

    touch GenomeAnalysisTK.jar
    """
}
