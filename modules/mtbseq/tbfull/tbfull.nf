nextflow.enable.dsl = 2
// NOTE: To properly setup the gatk inside the docker image
// - Download the gatk-3.8.0 tar file from here https://console.cloud.google.com/storage/browser/gatk-software/package-archive/gatk;tab=objects?prefix=&forceOnObjectsSortingFiltering=false
// - tar -xvf GATK_TAR_FILE
// - gatk-register gatk_folder/gatk_jar

params.resultsDir = "${params.outdir}/tbfull"
params.saveMode = 'copy'
params.shouldPublish = true
params.minbqual = 13
params.mincovf = 4
params.mincovr = 4
params.minphred = 4
params.minfreq = 75

process TBFULL {
    tag "${genomeName}"
    publishDir params.resultsDir, mode: params.saveMode, enabled: params.shouldPublish

    input:
    tuple val(genomeFileName), path("${genomeFileName}_${params.mtbseq_library_name}_R?.fastq.gz")
    path(gatk_jar)
    env USER

    output:
    val("${genomeFileName}")
    path("${genomeFileName}")
    path("${genomeFileName}/Called/*tab")
    path("${genomeFileName}/Position_Tables/*tab")

    script:

    """

    gatk-register ${gatk_jar}

    MTBseq --step TBfull --thread ${task.cpus} \
        --minbqual ${params.minbqual} \
        --mincovf ${params.mincovf} \
        --mincovr ${params.mincovr} \
        --minphred ${params.minphred} \
        --minfreq ${params.minfreq}
    """

    stub:
    """
    sleep \$[ ( \$RANDOM % 10 )  + 1 ]s

    mkdir ${genomeFileName}
    mkdir ${genomeFileName}/GATK_Bam
    touch ${genomeFileName}/GATK_Bam/${genomeFileName}_${params.library_name}.gatk.bam
    touch ${genomeFileName}/GATK_Bam/${genomeFileName}_${params.library_name}.gatk.bai
    touch ${genomeFileName}/GATK_Bam/${genomeFileName}_${params.library_name}.gatk.bamlog
    touch ${genomeFileName}/GATK_Bam/${genomeFileName}_${params.library_name}.gatk.grp
    touch ${genomeFileName}/GATK_Bam/${genomeFileName}_${params.library_name}.gatk.intervals
    mkdir ${genomeFileName}/Bam
    mkdir ${genomeFileName}/Bam/${genomeFileName}
    touch ${genomeFileName}/Bam/${genomeFileName}_${params.library_name}.bam
    touch ${genomeFileName}/Bam/${genomeFileName}_${params.library_name}.bai
    touch ${genomeFileName}/Bam/${genomeFileName}_${params.library_name}.bamlog
    mkdir ${genomeFileName}/Called
    touch ${genomeFileName}/Called/${genomeFileName}_${params.library_name}.gatk_position_uncovered_cf${params.mincovf}_cr${params.mincovr}_fr${params.minfreq}_ph${params.minphred}_outmode000.tab
    touch ${genomeFileName}/Called/${genomeFileName}_${params.library_name}.gatk_position_variants_cf${params.mincovf}_cr${params.mincovr}_fr${params.minfreq}_ph${params.minphred}_outmode000.tab
    mkdir ${genomeFileName}/Mpileup
    touch ${genomeFileName}/Mpileup/${genomeFileName}_${params.library_name}.gatk.mpileup
    touch ${genomeFileName}/Mpileup/${genomeFileName}_${params.library_name}.gatk.mpileuplog
    mkdir ${genomeFileName}/Classification
    touch ${genomeFileName}/Classification/Strain_Classification.tab
    mkdir ${genomeFileName}/Position_Tables
    touch ${genomeFileName}/Position_Tables/${genomeFileName}_${params.library_name}.gatk_position_table.tab
    mkdir ${genomeFileName}/Statistics
    touch ${genomeFileName}/Statistics/Mapping_and_Variant_Statistics.tab

    echo " MTBseq --step TBfull --threads ${task.cpus} \
        --minbqual ${params.minbqual} \
        --mincovf ${params.mincovf} \
        --mincovr ${params.mincovr} \
        --minphred ${params.minphred} \
        --minfreq ${params.minfreq} "
    """
}
