# USAGE

The simplest use case is to analyze a few genomes on a local environment. Almost all aspects are customizable but for the sake of brevity, a bare bones guide for any beginner user is as shown below

- [ ] 1. Make sure that you have installed MTBseq-nf following our [INSTALL]("./INSTALL") docs


- [ ] 2. Create a `fastqs` folder and add your genomes inside it

They should follow the pattern `SAMPLE_R1.fastq.gz` and `SAMPLE_R2.fastq.gz`

> **Note**
> Please don't use the `_` character on sample names. (ie.: Don't use `sample_001_R1.fastq.gz`, using `sample001_R1.fastq.gz` would be better)

> **Note**
> We suggest the usage of [fetchngs](https://nf-co.re/fetchngs) for fetching fastq
> data from private or public respositories

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

>**Warn**
> If you want to use `parallel` mode, you may add the `--parallel` parameter to the end of the comand 

## Normal and Parallel workflows

This pipeline has two execution types: normal and parallel and here is a dag example for them!

The execution type is determined by adding or not the `parallel` parameter

![](./resources/mtbseq_nf_workflow.png)

### What are the differences between `Normal` and `Parallel` analysis modes?

A normal MTBseq run would use `MTBseq --step full` on each sample, not allowing parallelization of secondary steps like `TB BWA` and `TB Variants` and leaving the resource control to MTBseq.

Using `parallel` run we enforce the parallelization of each step. The main advantage of it is the precise resource usage as the steps are controled by Nextflow, and some steps require less CPUs and RAM than other, this allow us to optimize the run time and resource costs.

