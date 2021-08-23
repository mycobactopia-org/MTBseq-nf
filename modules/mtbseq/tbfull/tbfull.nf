nextflow.enable.dsl = 2

params.results_dir = "${params.outdir}/tbfull"
params.save_mode = 'copy'
params.should_publish = true
params.minbqual = 13
params.mincovf = 4
params.mincovr = 4
params.minphred = 4
params.minfreq = 75

process TBFULL {
    tag "${genomeFileName}"
    publishDir params.results_dir, mode: params.save_mode, enabled: params.should_publish

    input:
    tuple val(genomeFileName), path("${genomeFileName}_${params.library_name}_R?.fastq.gz")
    path(gatk_jar)
    env(USER)

    output:
    tuple val("${genomeFileName}"), path("${genomeFileName}")
    path("${genomeFileName}/Called/*tab"), emit: variants_table
    path("${genomeFileName}/Position_Tables/*tab"), emit: position_table

    script:

    """

    gatk-register ${gatk_jar}

    MTBseq --step TBfull \
    --thread ${task.cpus} \
    --minbqual ${params.minbqual} \
    --mincovf ${params.mincovf} \
    --mincovr ${params.mincovr} \
    --minphred ${params.minphred} \
    --minfreq ${params.minfreq} \
    1>>.command.out \
    2>>.command.err \
    || true               # NOTE This is a hack to overcome the exit status 1 thrown by mtbseq


    mkdir ${genomeFileName}
    mv Amend ${genomeFileName}
    mv Bam ${genomeFileName}
    mv Called ${genomeFileName}
    mv Classification ${genomeFileName}
    mv GATK_Bam ${genomeFileName}
    mv Groups ${genomeFileName}
    mv Joint ${genomeFileName}
    mv Mpileup ${genomeFileName}
    mv Position_Tables ${genomeFileName}
    mv Statistics ${genomeFileName}

    """

    stub:
    """
    echo " MTBseq --step TBfull --threads ${task.cpus} \
    --minbqual ${params.minbqual} \
    --mincovf ${params.mincovf} \
    --mincovr ${params.mincovr} \
    --minphred ${params.minphred} \
    --minfreq ${params.minfreq} "

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

    """

}
