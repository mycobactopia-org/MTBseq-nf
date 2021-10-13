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
    tag "${params.project_name}"
    publishDir params.results_dir, mode: params.save_mode, enabled: params.should_publish

    input:
    path("*")
    path(gatk_jar)
    tuple path("${params.mtb_ref_name}.*"), path(resilist), path(intregions), path(categories), path(basecalib)
    env(USER)

    output:
    path("Called")
    path("Position_Tables")
    path("Classification")
    path("Statistics")
    path("Called/*_${params.library_name}*gatk_position_variants*.tab"), emit: position_variants
    path("Position_Tables/*_${params.library_name}*.gatk_position_table.tab"), emit: position_tables

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

    mkdir GATK_Bam
    touch GATK_Bam/${genomeFileName}_${params.library_name}.gatk.bam
    touch GATK_Bam/${genomeFileName}_${params.library_name}.gatk.bai
    touch GATK_Bam/${genomeFileName}_${params.library_name}.gatk.bamlog
    touch GATK_Bam/${genomeFileName}_${params.library_name}.gatk.grp
    touch GATK_Bam/${genomeFileName}_${params.library_name}.gatk.intervals
    mkdir Bam
    mkdir Bam/${genomeFileName}
    touch Bam/${genomeFileName}_${params.library_name}.bam
    touch Bam/${genomeFileName}_${params.library_name}.bai
    touch Bam/${genomeFileName}_${params.library_name}.bamlog
    mkdir Called
    touch Called/${genomeFileName}_${params.library_name}.gatk_position_uncovered_cf${params.mincovf}_cr${params.mincovr}_fr${params.minfreq}_ph${params.minphred}_outmode000.tab
    touch Called/${genomeFileName}_${params.library_name}.gatk_position_variants_cf${params.mincovf}_cr${params.mincovr}_fr${params.minfreq}_ph${params.minphred}_outmode000.tab
    mkdir Mpileup
    touch Mpileup/${genomeFileName}_${params.library_name}.gatk.mpileup
    touch Mpileup/${genomeFileName}_${params.library_name}.gatk.mpileuplog
    mkdir Classification
    touch Classification/Strain_Classification.tab
    mkdir Position_Tables
    touch Position_Tables/${genomeFileName}_${params.library_name}.gatk_position_table.tab
    mkdir Statistics
    touch Statistics/Mapping_and_Variant_Statistics.tab

    """

}
