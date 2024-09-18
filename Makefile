# https://makefiletutorial.com/

run_stub:
	bash ./data/mock_data/generate_mock_files.sh && nextflow run main.nf -stub-run

run_dev:
	nextflow run main.nf -profile dev -resume -with-tower
