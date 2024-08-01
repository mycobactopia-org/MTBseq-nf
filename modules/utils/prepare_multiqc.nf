include { MULTIQC                 } from '../nf-core/multiqc/main'

process PREPARE_MULTIQC {

    label 'process_medium'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/fastqc:0.12.1--hdfd78af_0' :
        'biocontainers/fastqc:0.12.1--hdfd78af_0' }"

 input:
        path(mapping_variants_statistics)
        path(Strain_Classification)
        path(groups)
        path(distance_matrix)

    output:
        path  "*.tsv"           , emit: multiqc_files

    shell:
        '''
        awk '{sub(/[^\t]*/,"");sub(/\t/,"")} 1' Mapping_and_Variant_Statistics.tab | awk '{gsub(/'\''/,"")} 1' >> Mapping_and_Variant_Statistics.tsv
        awk '{sub(/[^\t]*/,"");sub(/\t/,"")} 1' Strain_Classification.tab | awk '{gsub(/'\''/,"")} 1'>> Strain_Classification.tsv
        awk '/### Output as lists:/,0{if (!/### Output as lists:/) print $0}' *.groups | sed '1i Sample\tClustering Group' >> ClusterGroups.tsv
        make_symmetric_matrix.py !{distance_matrix}
        '''

}
