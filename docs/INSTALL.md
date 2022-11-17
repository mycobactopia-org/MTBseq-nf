# Installation instructions

## Prerequisites

For a successful Installation of this pipeline, you'll need the following programs

1. [Java version 11 or above](https://www.nextflow.io/docs/latest/getstarted.html#requirements)
1. [Nextflow](https://www.nextflow.io/docs/latest/getstarted.html#installation)
1. [Docker](https://docs.docker.com/engine/install/) or [Conda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html)

## Clone the project from Github

```terminal
git clone https://github.com/mtb-bioinformatics/MTBseq-nf.git
```

> **Note**
> If you're using Docker the pipeline is ready to use and
> it's not necessary to setup the conda environment

## Conda environment creation

> **Note**
> These instructions will generate a conda environment in the project folder

If you're planing to run the pipeline using conda, you'll need to setup the
conda enviroment prior running it.

- [ ] At the project root folder run the following code:

```terminal
cd conda_envs
bash setup_conda_envs.sh
```

## Docker based execution

For environments where docker is available, the pipeline will download a custom container and run the analysis.

This step will generate the environment at this folder, nextflow will take care of activating this environment when executing the pipeline

You can confirm the setup by activating that environment and using the `nextflow info`  command

```terminal
conda activate ./mtbseq-nf-env

nextflow info 
```
