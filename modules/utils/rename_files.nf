process RENAME_FILES {
    tag "${meta.id}"
    conda "bioconda::mtbseq=1.1.0"
    container "${'quay.io/biocontainers/mtbseq:1.1.0--hdfd78af_0'}"

    input:
        tuple val(meta), path(reads)

    output:
        path("${meta.id}_${params.library_name}_R?.fastq.gz")

    script:
        """
        echo "Renaming ${reads} files as per MTBseq requirements."

        cp ${reads[0]} ${meta.id}_${params.library_name}_R1.fastq.gz
        cp ${reads[1]} ${meta.id}_${params.library_name}_R2.fastq.gz
        """

}
