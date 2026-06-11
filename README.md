## Introduction

**MTBseq-nf**: Enabling Scalable Tuberculosis Genomics “Big Data” Analysis Through a User-Friendly Nextflow Wrapper for MTBseq Pipeline

> [!IMPORTANT]
> **For reproducing the publication-specific analysis, use the [`v1.0.0`](https://github.com/mycobactopia-org/MTBseq-nf/releases/tag/v1.0.0) tag of the pipeline:**
>
> ```bash
> nextflow run mycobactopia-org/MTBseq-nf -r v1.0.0 -profile <docker/singularity/.../institute> --input samplesheet.csv --outdir <OUTDIR>
> ```

> [!NOTE]
> This pipeline is being actively maintained and updated to follow the [nf-core](https://nf-co.re) community standards. The `master` branch (i.e. everything after the `v1.0.0` tag) reflects the latest changes and may differ from the published version.

## Usage

> [!NOTE]
> If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/usage/installation) on how to set-up Nextflow. Make sure to [test your setup](https://nf-co.re/docs/usage/introduction#how-to-run-a-pipeline) with `-profile test` before running the workflow on actual data.

<!-- TODO nf-core: Describe the minimum required steps to execute the pipeline, e.g. how to prepare samplesheets.
     Explain what rows and columns represent. For instance (please edit as appropriate):

First, prepare a samplesheet with your input data that looks as follows:

`samplesheet.csv`:

```csv
sample,fastq_1,fastq_2
CONTROL_REP1,AEG588A1_S1_L002_R1_001.fastq.gz,AEG588A1_S1_L002_R2_001.fastq.gz
```

Each row represents a fastq file (single-end) or a pair of fastq files (paired end).

-->

Now, you can run the pipeline using:

<!-- TODO nf-core: update the following command to include all required parameters for a minimal example -->

```bash
nextflow run mycobactopia-org/MTBseq-nf \
   -profile <docker/singularity/.../institute> \
   --input samplesheet.csv \
   --outdir <OUTDIR>
```

> [!WARNING]
> Please provide pipeline parameters via the CLI or Nextflow `-params-file` option. Custom config files including those provided by the `-c` Nextflow option can be used to provide any configuration _**except for parameters**_;
> see [docs](https://nf-co.re/usage/configuration#custom-configuration-files).

## Credits

mycobactopia-org/MTBseq-nf was originally written by Abhinav Sharma (@abhi18av) and Davi Marcon (@mxrcon).

We thank the following people for their extensive assistance in the development of this pipeline:

<!-- TODO nf-core: If applicable, make list of people who have also contributed -->

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

## Documentation

For detailed instructions on installation, usage, and output interpretation, please refer to the [documentation](docs/README.md).

## Citations

If you use mycobactopia-org/MTBseq-nf for your analysis, please cite it using the following publication:

> Sharma, A.; Marcon, D.J.; Loubser, J.; Lima, K.V.B.; van der Spuy, G.; Conceição, E.C. MTBseq-nf: Enabling Scalable Tuberculosis Genomics "Big Data" Analysis Through a User-Friendly Nextflow Wrapper for MTBseq Pipeline. *Microorganisms* 2025, 13, 2685. [https://doi.org/10.3390/microorganisms13122685](https://doi.org/10.3390/microorganisms13122685)

BibTeX entry:

```bibtex
@article{Sharma2025,
   author = {Sharma, Abhinav and Marcon, Davi Josué and Loubser, Johannes and Lima, Karla Valéria Batista and van der Spuy, Gian and Conceição, Emilyn Costa},
   title = {MTBseq-nf: Enabling Scalable Tuberculosis Genomics "Big Data" Analysis Through a User-Friendly Nextflow Wrapper for MTBseq Pipeline},
   journal = {Microorganisms},
   volume = {13},
   number = {12},
   ISSN = {2076-2607},
   DOI = {10.3390/microorganisms13122685},
   year = {2025},
   type = {Journal Article}
}
```

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

This pipeline uses code and infrastructure developed and maintained by the [nf-core](https://nf-co.re) community, reused here under the [MIT license](https://github.com/nf-core/tools/blob/master/LICENSE).

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
