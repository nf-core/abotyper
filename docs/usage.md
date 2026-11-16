# nf-core/abotyper: Usage

## :warning: Please read this documentation on the nf-core website: [https://nf-co.re/abotyper/usage](https://nf-co.re/abotyper/usage)

> _Documentation of pipeline parameters is generated automatically from the pipeline schema and can no longer be found in markdown files._

## Introduction

nf-core/abotyper types the ABO blood group from Oxford Nanopore amplicon reads covering ABO exons 2–7. Reads are aligned to the ABO RefSeqGene (NG_006669.2), variants are called with Clair3 and checked against an independent BAM pileup, and the calls are scored against a diagnostic variant panel to give the ABO phenotype, genotype and ISBT allele names, including A1/A2 subtypes. Results are for research use only; see [intended use](intended_use.md).

## Samplesheet input

You will need to create a samplesheet with information about the samples you would like to analyse before running the pipeline. Use this parameter to specify its location. It has to be a comma-separated file with a header row and the columns shown below. The pipeline types Oxford Nanopore reads (see [Sequencing platform compatibility](../README.md#sequencing-platform-compatibility)).

```bash
--input '[path to samplesheet file]'
```

```csv title="samplesheet.csv"
sample,fastq_1,control,barcode
SAMPLE1_barcode01,/path/to/SAMPLE1_barcode01.fastq.gz,,
SAMPLE2,/path/to/SAMPLE2.fastq.gz,,02
NTC,/path/to/NTC.fastq.gz,negative,95
POSCTRL,/path/to/POSCTRL.fastq.gz,positive,96
```

| Column    | Description                                                                                                                                                                                                                          |
| --------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| `sample`  | Sample name, used in every result table: letters, digits, `.`, `_` and `-`. Without a `barcode` value, a trailing `_barcode<NN>` is reported as the barcode. See below.                                                              |
| `fastq_1` | Path to the sample's demultiplexed Nanopore reads, with the extension `.fastq`, `.fq`, `.fastq.gz` or `.fq.gz`.                                                                                                                      |
| `control` | Optional. `negative` for a no-template control, `positive` for a positive control, empty for a sample. Controls are only recognised from this column, never from the sample name.                                                    |
| `barcode` | Optional. The barcode to report for the sample (`7`, `07` or `barcode07`). When it is given, the sample is reported under its full `sample` name; when it is empty, the barcode is taken from a trailing `_barcode<NN>` in the name. |

An [example samplesheet](../assets/samplesheet.csv) has been provided with the pipeline.

### Sample names

Sample names may contain letters, digits, `.`, `_` and `-`, for example `ABO_AA_barcode01`, `Sample-45-44874_barcode25`, `patient-01` or `S1.r1`. The sample name and barcode in the result tables are set when the samplesheet is read: with a `barcode` value, the sample is reported under its full name with that barcode; otherwise a name ending in `_barcode<digits>` is reported as that sample and barcode (`S1_barcode07` is sample `S1`, barcode 07), and any other name is reported whole, with barcode `barcode00`. Libraries that share a sample name but have different barcodes are compared as replicates.

### Controls

A sample marked `negative` is expected to have few or no ABO reads: low depth is reported as the expected result, and reads above the depth threshold are flagged as possible contamination. A sample marked `positive` must reach normal depth, or it is flagged as a failed control. Every other sample is typed and checked in the same way, whatever its name.

### Multiple runs of the same sample

Rows with the same `sample` name are treated as one sample: their reads are concatenated before alignment. The rows must have the same `control` and `barcode` values.

```csv title="samplesheet.csv"
sample,fastq_1,control,barcode
SAMPLE1_barcode01,/path/to/run1/SAMPLE1_barcode01.fastq.gz,,
SAMPLE1_barcode01,/path/to/run2/SAMPLE1_barcode01.fastq.gz,,
```

To type the same sample twice as separate replicates, give each library its own name (for example `SAMPLE1_barcode01` and `SAMPLE1_barcode13`); replicate libraries of the same sample are compared in the result tables.

## Running the pipeline

The typical command for running the pipeline is as follows:

```bash
nextflow run nf-core/abotyper --input ./samplesheet.csv --outdir ./results -profile docker
```

This will launch the pipeline with the `docker` configuration profile. See below for more information about profiles.

Note that the pipeline will create the following files in your working directory:

```bash
work                # Directory containing the nextflow working files
<OUTDIR>            # Finished results in specified location (defined with --outdir)
.nextflow_log       # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs.
```

If you wish to repeatedly use the same parameters for multiple runs, rather than specifying each flag in the command, you can specify these in a params file.

Pipeline settings can be provided in a `yaml` or `json` file via `-params-file <file>`.

> [!WARNING]
> Do not use `-c <file>` to specify parameters as this will result in errors. Custom config files specified with `-c` must only be used for [tuning process resource specifications](https://nf-co.re/docs/running/run-pipelines#configuring-pipelines), other infrastructural tweaks (such as output directories), or module arguments (args).

The above pipeline run specified with a params file in yaml format:

```bash
nextflow run nf-core/abotyper -profile docker -params-file params.yaml
```

with:

```yaml title="params.yaml"
input: './samplesheet.csv'
outdir: './results/'
<...>
```

You can also generate such `YAML`/`JSON` files via [nf-core/launch](https://nf-co.re/launch).

### Clair3 model

Clair3 needs a model that matches the flow cell and basecaller. The pipeline picks one in this order:

1. `--clair3_model_url <url>`: download the model folder (`pileup.pt`, `full_alignment.pt`) from this URL.
2. `--clair3_model_path <folder>`: use your own model folder.
3. Otherwise the model packaged in the Clair3 container, named by `--clair3_model` (default `r1041_e82_400bps_sup_v500`, for R10.4.1 flow cells basecalled with a super-accuracy model).

Give a URL or a folder, not both. The packaged models are only in the Docker, Singularity and Podman containers; with other profiles (conda, apptainer) the `--clair3_model` model is downloaded by name instead, and the log says so.

```bash
# packaged model (default)
nextflow run nf-core/abotyper -profile docker --input samplesheet.csv --outdir results
# your own model folder
nextflow run nf-core/abotyper -profile docker --input samplesheet.csv --outdir results --clair3_model_path /path/to/r1041_e82_400bps_sup_v500
```

### Updating the pipeline

When you run the above command, Nextflow automatically pulls the pipeline code from GitHub and stores it as a cached version. When running the pipeline after this, it will always use the cached version if available - even if the pipeline has been updated since. To make sure that you're running the latest version of the pipeline, make sure that you regularly update the cached version of the pipeline:

```bash
nextflow pull nf-core/abotyper
```

### Reproducibility

It is a good idea to specify the pipeline version when running the pipeline on your data. This ensures that a specific version of the pipeline code and software are used when you run your pipeline. If you keep using the same tag, you'll be running the same version of the pipeline, even if there have been changes to the code since.

First, go to the [nf-core/abotyper releases page](https://github.com/nf-core/abotyper/releases) and find the latest pipeline version - numeric only (eg. `1.3.1`). Then specify this when running the pipeline with `-r` (one hyphen) - eg. `-r 1.3.1`. Of course, you can switch to another version by changing the number after the `-r` flag.

This version number will be logged in reports when you run the pipeline, so that you'll know what you used when you look back in the future. For example, at the bottom of the MultiQC reports.

To further assist in reproducibility, you can use share and reuse [parameter files](#running-the-pipeline) to repeat pipeline runs with the same settings without having to write out a command with every single parameter.

> [!TIP]
> If you wish to share such profile (such as upload as supplementary material for academic publications), make sure to NOT include cluster specific paths to files, nor institutional specific profiles.

## Core Nextflow arguments

> [!NOTE]
> These options are part of Nextflow and use a _single_ hyphen (pipeline parameters use a double-hyphen)

### `-profile`

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments.

Several generic profiles are bundled with the pipeline which instruct the pipeline to use software packaged using different methods (Docker, Singularity, Podman, Shifter, Charliecloud, Apptainer, Conda) - see below.

> [!IMPORTANT]
> We highly recommend the use of Docker or Singularity containers for full pipeline reproducibility, however when this is not possible, Conda is also supported.

The pipeline also dynamically loads configurations from [https://github.com/nf-core/configs](https://github.com/nf-core/configs) when it runs, making multiple config profiles for various institutional clusters available at run time. For more information and to check if your system is supported, please see the [nf-core/configs documentation](https://github.com/nf-core/configs#documentation).

Note that multiple profiles can be loaded, for example: `-profile test,docker` - the order of arguments is important!
They are loaded in sequence, so later profiles can overwrite earlier profiles.

If `-profile` is not specified, the pipeline will run locally and expect all software to be installed and available on the `PATH`. This is _not_ recommended, since it can lead to different results on different machines dependent on the computer environment.

- `test`
  - A profile with a complete configuration for automated testing
  - Includes links to test data so needs no other parameters
- `docker`
  - A generic configuration profile to be used with [Docker](https://docker.com/)
- `singularity`
  - A generic configuration profile to be used with [Singularity](https://sylabs.io/docs/)
- `podman`
  - A generic configuration profile to be used with [Podman](https://podman.io/)
- `shifter`
  - A generic configuration profile to be used with [Shifter](https://nersc.gitlab.io/development/shifter/how-to-use/)
- `charliecloud`
  - A generic configuration profile to be used with [Charliecloud](https://charliecloud.io/)
- `apptainer`
  - A generic configuration profile to be used with [Apptainer](https://apptainer.org/)
- `wave`
  - A generic configuration profile to enable [Wave](https://seqera.io/wave/) containers. Use together with one of the above (requires Nextflow `24.03.0-edge` or later).
- `conda`
  - A generic configuration profile to be used with [Conda](https://conda.io/docs/). Please only use Conda as a last resort i.e. when it's not possible to run the pipeline with Docker, Singularity, Podman, Shifter, Charliecloud, or Apptainer.
- `conda_lock`
  - Use with `conda` or `mamba` to install the exact package builds the nf-core modules were tested with. See [Conda](#conda).

### Conda

With `-profile conda` (or `mamba`), each tool's environment is built from its module's `environment.yml` when the pipeline runs. That fixes the tool version, but its dependencies can change as new builds are published to conda-forge and bioconda. To install the exact builds the modules were tested with, add `conda_lock`:

```bash
nextflow run nf-core/abotyper -profile conda,conda_lock --input samplesheet.csv --outdir <OUTDIR>
```

`conda_lock` covers the nf-core modules that ship lock files (bcftools, samtools and MultiQC) on linux-64 and linux-aarch64, chosen from the machine Nextflow runs on. The other tools still come from their `environment.yml`.

> [!WARNING]
> Do not use `conf/containers_conda_lock_files_amd64.config` or `conf/containers_conda_lock_files_arm64.config`. nf-core tools generates them, but they point at lock files that conda cannot read, so environment creation fails. `conda_lock` uses converted copies of the same lock files.

### `-resume`

Specify this when restarting a pipeline. Nextflow will use cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously. For input to be considered the same, not only the names must be identical but the files' contents as well. For more info about this parameter, see [this blog post](https://www.nextflow.io/blog/2019/demystifying-nextflow-resume.html).

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

### `-c`

Specify the path to a specific config file (this is a core Nextflow command). See the [nf-core website documentation](https://nf-co.re/usage/configuration) for more information.

## Custom configuration

### Resource requests

Whilst the default requirements set within the pipeline will hopefully work for most people and with most input data, you may find that you want to customise the compute resources that the pipeline requests. Each step in the pipeline has a default set of requirements for number of CPUs, memory and time. For most of the pipeline steps, if the job exits with any of the error codes specified [here](https://github.com/nf-core/rnaseq/blob/4c27ef5610c87db00c3c5a3eed10b1d161abf575/conf/base.config#L18) it will automatically be resubmitted with higher resources request (2 x original, then 3 x original). If it still fails after the third attempt then the pipeline execution is stopped.

To change the resource requests, please see the [max resources](https://nf-co.re/docs/running/configuration/nextflow-for-your-system#set-max-resources) and [customise process resources](https://nf-co.re/docs/running/configuration/nextflow-for-your-system#customize-process-resources) section of the nf-core website.

### Custom Containers

In some cases, you may wish to change the container or conda environment used by a pipeline steps for a particular tool. By default, nf-core pipelines use containers and software from the [biocontainers](https://biocontainers.pro/) or [bioconda](https://bioconda.github.io/) projects. However, in some cases the pipeline specified version maybe out of date.

To use a different container from the default container or conda environment specified in a pipeline, please see the [updating tool versions](https://nf-co.re/docs/running/configuration/nextflow-for-your-system#update-tool-versions) section of the nf-core website.

### Custom Tool Arguments

A pipeline might not always support every possible argument or option of a particular tool used in pipeline. Fortunately, nf-core pipelines provide some freedom to users to insert additional parameters that the pipeline does not include by default.

To learn how to provide additional arguments to a particular tool of the pipeline, please see the [customising tool arguments](https://nf-co.re/docs/running/configuration/nextflow-for-your-system#modifying-tool-arguments) section of the nf-core website.

### nf-core/configs

In most cases, you will only need to create a custom config as a one-off but if you and others within your organisation are likely to be running nf-core pipelines regularly and need to use the same settings regularly it may be a good idea to request that your custom config file is uploaded to the `nf-core/configs` git repository. Before you do this please can you test that the config file works with your pipeline of choice using the `-c` parameter. You can then create a pull request to the `nf-core/configs` repository with the addition of your config file, associated documentation file (see examples in [`nf-core/configs/docs`](https://github.com/nf-core/configs/tree/master/docs)), and amending [`nfcore_custom.config`](https://github.com/nf-core/configs/blob/master/nfcore_custom.config) to include your custom profile.

See the main [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html) for more information about creating your own configuration files.

If you have any questions or issues please send us a message on [Slack](https://nf-co.re/join/slack) on the [`#configs` channel](https://nfcore.slack.com/channels/configs).

## Running in the background

Nextflow handles job submissions and supervises the running jobs. The Nextflow process must run until the pipeline is finished.

The Nextflow `-bg` flag launches Nextflow in the background, detached from your terminal so that the workflow does not stop if you log out of your session. The logs are saved to a file.

Alternatively, you can use `screen` / `tmux` or similar tool to create a detached session which you can log back into at a later time.
Some HPC setups also allow you to run nextflow within a cluster job submitted your job scheduler (from where it submits more jobs).

## Nextflow memory requirements

In some cases, the Nextflow Java virtual machines can start to request a large amount of memory.
We recommend adding the following line to your environment to limit this (typically in `~/.bashrc` or `~./bash_profile`):

```bash
NXF_OPTS='-Xms1g -Xmx4g'
```
