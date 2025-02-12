
<h1>
  <picture>
    <source media="(prefers-color-scheme: light)" srcset="docs/images/nf-core-abotyper_logo_dark.png">
    <img alt="nf-core/abotyper" src="docs/images/nf-core-abotyper_logo_light.png">
  </picture>
</h1>


# ABO blood typing using Oxford Nanopore MinION sequencing
![nf-core/abotyper metro map](docs/images/nf-core-abotyper-metro-map.jpg)

ABO sequences were acquired from the NCBI dbRBC database:

[https://www.ncbi.nlm.nih.gov/projects/gv/mhc/xslcgi.cgi?cmd=bgmut/home](https://www.ncbi.nlm.nih.gov/projects/gv/mhc/xslcgi.cgi?cmd=bgmut/home)

See [https://ftp.ncbi.nlm.nih.gov/pub/mhc/rbc/Final%20Archive/Excel_and_PowerPoint/](https://ftp.ncbi.nlm.nih.gov/pub/mhc/rbc/Final%20Archive/Excel_and_PowerPoint/) for some literature.

# Required tools

The pipeline makes use of the following core dependencies:

```yaml
- bioconda::fastqc=0.12.1
- bioconda::bwa=0.7.17
- conda-forge::ncurses
- bioconda::samtools=1.19.2
- bioconda::minimap2=2.26
- conda-forge::biopython=1.83
- python=3.10
- pip
- pip:
    - numpy>=1.26.0
    - Bio>=1.6.0
    - biopython>=1.8o
    - openpyxl>=3.1.0
    - pandas>=2.2.0
    - pysam>=0.22.0
    - matplotlib>=3.8.0
    - XlsxWriter>=3.2.0
    - multiqc>=1.18
```

A complete list of dependencies is found in the assets folder `assets/conda.yml`.

# Required input files structure

Ensure that all input fastq files have a naming convention that matches this regular expression (`regex`)

```python
## python regex for matching samples
pattern = r"^(IMM|INGS|NGS|[A-Z0-9]+)(-[0-9]+-[0-9]+)?_barcode\d+$"
```

The regex does the following:

- `^(IMM|INGS|NGS|[A-Z0-9]+)` allows for files strating with the prefixes IMM, INGS, NGS, or any combination of letters `A-to-Z` and digits `0-to-9`.
- `(-[0-9]+-[0-9]+)?` handles optional segments of digits separated by a dash(-).
- `_barcode\d+$` ensures the filename ends with_barcode followed by digits to denote barcode numbers.

There is a filename handling logic in the code `filename.split("_")` assumes that the barcode is always the last part of the filename. The names are split into basename and barcode which are then used in later reporting. Please Adjust this if necessary based on actual filename structure in your assays.

Here are a few examples of acceptable input file names:

```txt
NGSPOS_barcode13.fastq
NGSNEG_barcode12.fastq
INGSPOS_barcode01.fastq
INGSNEG_barcode96.fastq
IBTGSPOS_barcode19.fastq
2025705_barcode14.fastq
IMM-24-44988_barcode51.fastq
Sample1-2024-12345_barcode22.fastq
```

# Testing without `nextflow`

The pipeline can be tested of single input file by cloning this repo and installing all dependencies above, then running the following commands:

```python
python bin/AnalyzeAbo_Main.py  \
  --reference="assets/A1_01_01_1_reference_Exon6.fasta" \
  --alleles="assets/ABO_Database.fasta" \
  --output="SampleName/exon6" \
  --analysis-type="READS" \
  --reads="SampleName.fastq" \

 python bin/AnalyzeAbo_Main.py  \
  --reference="assets/reads_bc51/A1_01_01_1_reference_Exon7.fasta" \
  --alleles="assets/input/ABO_Database.fasta" \
  --output="SampleName/exon7" \
  --analysis-type="READS" \
  --reads="SampleName.fastq"
```

Looping through a couple of samples with the above command will generate the following outputs:

## Output data structure from the testing

```yaml
OutputDirectoryName/
├── Sample1
│   ├── exon6
│   │   └── alignment
│   └── exon7
│       └── alignment
├── Sample2
│   ├── exon6
│   │   └── alignment
│   └── exon7
│       └── alignment
├── Sample3
│   ├── exon6
│   │   └── alignment
│   └── exon7
│       └── alignment
```

With individual files named as follows:

```yaml
OutputDirectoryName/
├── Sample1
│   ├── exon6
│   │   ├── ABOPhenotype.txt
│   │   ├── ABOReadPolymorphisms.txt
│   │   ├── alignment
│   │   │   ├── alignment.bam
│   │   │   ├── alignment.bam.bai
│   │   │   ├── AlignmentReference.fasta
│   │   │   ├── AlignmentReference.fasta.amb
│   │   │   ├── AlignmentReference.fasta.ann
│   │   │   ├── AlignmentReference.fasta.bwt
│   │   │   ├── AlignmentReference.fasta.pac
│   │   │   └── AlignmentReference.fasta.sa
│   │   └── ReadAlignmentSpreadsheet.csv
│   ├── exon7
│   │   ├── ABOPhenotype.txt
│   │   ├── ABOReadPolymorphisms.txt
│   │   ├── alignment
│   │   │   ├── alignment.bam
│   │   │   ├── alignment.bam.bai
│   │   │   ├── AlignmentReference.fasta
│   │   │   ├── AlignmentReference.fasta.amb
│   │   │   ├── AlignmentReference.fasta.ann
│   │   │   ├── AlignmentReference.fasta.bwt
│   │   │   ├── AlignmentReference.fasta.pac
│   │   │   └── AlignmentReference.fasta.sa
│   │   └── ReadAlignmentSpreadsheet.csv
│   ├── Sample1_exon6.log.txt
│   └── Sample1_exon7.log.txt
```

The `ABOPhenotype.txt` files from each sampe can then be collated using:

`python bin/Aggregate_ABO_reports.py OutputDirectoryName`

# The `nextflow` workflow

The steps above are simplified in a `NextFlow; https://www.nextflow.io/` pipeline that does all the above steps and streamlines installation of requisite software and tools with a single command.

Besides reproducability, nextflow offers several advantages over conventional `for loops`, including scalability, portability, and debugging/resumption of failed tasks.

Input files and output directory can be defined in the config files or provided directly in the commandline.

To analyse files with config, run:

- `nextflow run main.nf -resume` (user can override inputs and output using `--reads '*.fastq' --outdir 'ABO_results'` on the commandline).

We have also added the ability for the pipeline to automatically set-up a conda or docker based environment with all required tools and libraries.

Users may also opt for a workload manager such as `-profile slurm,docker|-profile slurm,conda`, is which case, all required modules docker/conda must be installed and loaded. The config slurm parameters must also be defined to ensure tasks are submitted to the correct resource queue/account.

For `conda` environment, it is advisable to prepare the working computer using mamba for easy resolution of environments.
Follow these steps to achieve better results.

```bash
mamba create -y -n abo-analysis-env
conda activate abo-analysis-env
mamba env update --file abo-analysis/assets/conda.yml --prune
conda deactivate

# If conda cativate fails, run:
source {path_to_anaconda}/anaconda3/etc/profile.d/conda.sh
```

To run without the workload manager but with a specific containerization, use:

- `nextflow run abo-analysis/main.nf -resume --outdir "$PWD/230128R_ABO_results" -with-conda abo-analysis-env` or
  `nextflow run abo-analysis/main.nf -resume --outdir "$PWD/230128R_ABO_results" -profile conda`
- `nextflow run abo-analysis/main.nf -resume --outdir "$PWD/230128R_ABO_results" -with-docker fmobegi/abo-analysis` or
  `nextflow run abo-analysis/main.nf -resume --outdir "$PWD/230128R_ABO_results" -profile docker`

# Renaming samples

The code by default renames samples using a tab file with `sequencingID` and `sampleName` (see `nextflow.config` file under `$params.renaming_file`).
This option is controlled by the parameter `$params.skip_renaming` and can be overridden via the commandline using option `--skip_renaming true` to skip the process.

# Results from the `Nextflow` pipeline will look something like this

```yaml
230128R_ABO_results/
├── ABO_result.txt
├── ABO_result.xlsx
├── execution_report.html
├── execution_timeline.html
├── execution_trace.txt
├── SampleName
│   ├── exon6
│   │   ├── ABOPhenotype.txt
│   │   ├── ABOReadPolymorphisms.txt
│   │   ├── alignment
│   │   │   ├── alignment.bam
│   │   │   ├── alignment.bam.bai
│   │   │   ├── AlignmentReference.fasta
│   │   │   ├── AlignmentReference.fasta.amb
│   │   │   ├── AlignmentReference.fasta.ann
│   │   │   ├── AlignmentReference.fasta.bwt
│   │   │   ├── AlignmentReference.fasta.pac
│   │   │   └── AlignmentReference.fasta.sa
│   │   └── ReadAlignmentSpreadsheet.csv
│   ├── exon7
│   │   ├── ABOPhenotype.txt
│   │   ├── ABOReadPolymorphisms.txt
│   │   ├── alignment
│   │   │   ├── alignment.bam
│   │   │   ├── alignment.bam.bai
│   │   │   ├── AlignmentReference.fasta
│   │   │   ├── AlignmentReference.fasta.amb
│   │   │   ├── AlignmentReference.fasta.ann
│   │   │   ├── AlignmentReference.fasta.bwt
│   │   │   ├── AlignmentReference.fasta.pac
│   │   │   └── AlignmentReference.fasta.sa
│   │   └── ReadAlignmentSpreadsheet.csv
│   ├── SampleName_exon6.log.txt
│   └── SampleName_exon7.log.txt
├── software_versions.txt
└── workflow.oncomplete.txt
```

Feel free to raise an issue or reach out if you need any support getting this tool running, or with suggestions for improvement.

## Credits

nf-core/abotyper was originally written by @fmobegi.

We thank the following people for their extensive assistance in the development of this pipeline:

<!-- TODO nf-core: If applicable, make list of people who have also contributed -->

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

For further information or help, don't hesitate to get in touch on the [Slack `#abotyper` channel](https://nfcore.slack.com/channels/abotyper) (you can join with [this invite](https://nf-co.re/join/slack)).

## Citations

<!-- TODO nf-core: Add citation for pipeline after first release. Uncomment lines below and update Zenodo doi and badge at the top of this file. -->
<!-- If you use nf-core/abotyper for your analysis, please cite it using the following doi: [10.5281/zenodo.XXXXXX](https://doi.org/10.5281/zenodo.XXXXXX) -->

<!-- TODO nf-core: Add bibliography of tools and data used in your pipeline -->

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

You can cite the `nf-core` publication as follows:

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
