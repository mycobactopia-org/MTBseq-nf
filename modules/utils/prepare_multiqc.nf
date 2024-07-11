process PREPARE_MULTIQC {

    input:
        path(mapping_variants_statistics)
        path(Strain_Classification)
        path(groups)
        path(distance_matrix)

    output:
        path  "*.tsv"           , emit: multiqc_files


    script:
        """
        awk '{sub(/[^\t]*/,"");sub(/\t/,"")} 1' Mapping_and_Variant_Statistics.tab | awk '{gsub(/'\''/,"")} 1' >> Mapping_and_Variant_Statistics.tsv
        awk '{sub(/[^\t]*/,"");sub(/\t/,"")} 1' Strain_Classification.tab | awk '{gsub(/'\''/,"")} 1'>> Strain_Classification.tsv
        awk '/### Output as lists:/,0{if (!/### Output as lists:/) print $0}' *.groups | sed '1i Sample\tClustering Group' >> ClusterGroups.tsv
        make_symmetric_matrix.py ${distance_matrix}
        """

}
