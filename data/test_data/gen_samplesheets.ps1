$csv = Import-Csv ./original_dataset.ena.deduplicated.csv

$DATASET_FOLDER="/home/ubuntu/DATASET/publication/fastq"

foreach ($r in $csv) {

    $libName = $r.run_alias.split("_")[1]
    $sampleName = $r.run_alias.split("_")[0]

    Write-Output "$sampleName,$libName,$DATASET_FOLDER/$($r.experiment_accession)_$($r.run_accession)_1.fastq.gz,$DATASET_FOLDER/$($r.experiment_accession)_$($r.run_accession)_2.fastq.gz" >> ./publication.90samples.csv

}
