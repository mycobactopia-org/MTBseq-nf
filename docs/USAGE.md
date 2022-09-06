# USAGE

The simplest use case is to analyze a few genomes on a local environment (laptop/server). Almost all aspects are customizable but for the sake of brevity, a bare bones guide for any beginner user is as shown below

## Minimal guide for local environment

- [ ] 1. Make sure that you have installed MTBseq-nf following our [INSTALL]("./INSTALL") docs

- [ ] 2. Create a `fastqs` folder and add your `fastq` files inside it

They should follow the pattern `SAMPLE_R1.fastq.gz` and `SAMPLE_R2.fastq.gz`

> **Note**
> Please don't use the `_` character on sample names i.e. don't use `sample_001_R1.fastq.gz`. Using `sample001_R1.fastq.gz` would be better

> **Note**
> We suggest the usage of [fetchngs](https://nf-co.re/fetchngs) for fetching fastq
> data from private or public repositories

> **Note**
> You may specify another folder containing the genomes by the `--folder` parameter

- [ ] 3.1 If you're using conda:

  ```terminal
  nextflow run main.nf -profile standard,conda
  ```

- [ ] 3.2 If you're using Docker

  ```terminal
  nextflow run main.nf -profile standard,docker
  ```

>**Note**
> If you want to use `parallel` mode, you may add the `--parallel` parameter to the end of the command

```terminal
nextflow run main.nf -profile standard,docker --parallel
```
