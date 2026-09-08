<h1>
  <picture>
    <source media="(prefers-color-scheme: light)" srcset="docs/images/nf-core-abotyper_logo_dark.png">
    <img alt="nf-core/abotyper" src="docs/images/nf-core-abotyper_logo_light.png">
  </picture>
</h1>

[![Open in GitHub Codespaces](https://img.shields.io/badge/Open_In_GitHub_Codespaces-black?labelColor=grey&logo=github)](https://github.com/codespaces/new/nf-core/abotyper)
[![GitHub Actions CI Status](https://github.com/nf-core/abotyper/actions/workflows/nf-test.yml/badge.svg)](https://github.com/nf-core/abotyper/actions/workflows/nf-test.yml)
[![GitHub Actions Linting Status](https://github.com/nf-core/abotyper/actions/workflows/linting.yml/badge.svg)](https://github.com/nf-core/abotyper/actions/workflows/linting.yml)[![AWS CI](https://img.shields.io/badge/CI%20tests-full%20size-FF9900?labelColor=000000&logo=Amazon%20AWS)](https://nf-co.re/abotyper/results)[![Cite with Zenodo](http://img.shields.io/badge/DOI-10.5281/zenodo.XXXXXXX-1073c8?labelColor=000000)](https://doi.org/10.5281/zenodo.XXXXXXX)
[![nf-test](https://img.shields.io/badge/unit_tests-nf--test-337ab7.svg)](https://www.nf-test.com)
[![Open in GitHub Codespaces](https://img.shields.io/badge/Open_In_GitHub_Codespaces-black?labelColor=grey&logo=github)](https://github.com/codespaces/new/nf-core/abotyper)

[![Nextflow](https://img.shields.io/badge/version-%E2%89%A525.10.4-green?style=flat&logo=nextflow&logoColor=white&color=%230DC09D&link=https%3A%2F%2Fnextflow.io)](https://www.nextflow.io/)
[![nf-core template version](https://img.shields.io/badge/nf--core_template-4.1.0-green?style=flat&logo=nfcore&logoColor=white&color=%2324B064&link=https%3A%2F%2Fnf-co.re)](https://github.com/nf-core/tools/releases/tag/4.1.0)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
[![Launch on Seqera Platform](https://img.shields.io/badge/Launch%20%F0%9F%9A%80-Seqera%20Platform-%234256e7)](https://cloud.seqera.io/launch?pipeline=https://github.com/nf-core/abotyper)

[![Get help on Slack](http://img.shields.io/badge/slack-nf--core%20%23abotyper-4A154B?labelColor=000000&logo=slack)](https://nfcore.slack.com/channels/abotyper)[![Follow on Bluesky](https://img.shields.io/badge/bluesky-%40nf__core-1185fe?labelColor=000000&logo=bluesky)](https://bsky.app/profile/nf-co.re)[![Follow on Mastodon](https://img.shields.io/badge/mastodon-nf__core-6364ff?labelColor=FFFFFF&logo=mastodon)](https://mstdn.science/@nf_core)[![Watch on YouTube](http://img.shields.io/badge/youtube-nf--core-FF0000?labelColor=000000&logo=youtube)](https://www.youtube.com/c/nf-core)

# ABO blood typing using Oxford Nanopore MinION sequencing

**nf-core/abotyper** is a bioinformatics pipeline that analyses data obtained from Third Generation Sequencing of the `Homo sapiens ABO, alpha 1-3-N-acetylgalactosaminyltransferase and alpha 1-3-galactosyltransferase` (ABO) gene to deduce the ABO blood type.<br/>
It takes a samplesheet and FASTQ files as input, performs quality control (QC), mapping to the reference sequences, variant characterisation, and finally deduce the Blood Group Statistics based on known ABO-related Single nucleotide variants (SVNs).

![nf-core/abotyper metro map](docs/images/nf-core-abotyper-metro-map.jpg)

ABO sequences were acquired from the NCBI RefSeq and dbRBC databases:

- [ABO RefSeqGene (NG_006669.2)](https://www.ncbi.nlm.nih.gov/nuccore/NG_006669.2) - single combined reference used for alignment
- [dbMHC and IHWG data](https://ftp.ncbi.nlm.nih.gov/pub/mhc/mhc/Final%20Archive/)

## Pipeline steps

> [!NOTE]
> If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/get_started/environment_setup/overview) on how to set-up Nextflow. Make sure to [test your setup](https://nf-co.re/docs/get_started/run-your-first-pipeline) with `-profile test` before running the workflow on actual data.

1. **Read quality control** - Quality assessment of input FASTQ files using FastQC ([`FASTQC`](modules/nf-core/fastqc/))
2. **Read alignment** - Align reads to a single combined ABO reference (NG_006669.2) using Minimap2 ([`MINIMAP2_ALIGN`](modules/nf-core/minimap2/align/))
3. **Alignment statistics** - Generate comprehensive alignment metrics including coverage, flagstat, and detailed statistics:
   - Coverage analysis ([`SAMTOOLS_COVERAGE`](modules/nf-core/samtools/coverage/))
   - Flagstat metrics ([`SAMTOOLS_FLAGSTAT`](modules/nf-core/samtools/flagstat/))
   - Detailed statistics ([`SAMTOOLS_STATS`](modules/nf-core/samtools/stats/))
4. **Variant quantification** - Compute per-position allele frequencies and per-read haplotypes directly from the BAM file in a single pass ([`HAPLOSCAN`](modules/local/haploscan/))
5. **SNP extraction** - Extract and analyze ABO-relevant single nucleotide variants using the panel-driven coordinate system ([`ABO_GETABOSNPS`](modules/local/abo/getabosnps/))
6. **Phenotype prediction** - Predict ABO blood group phenotype from combined SNP and haplotype patterns ([`ABO_SNPS2PHENO`](modules/local/abo/snps2pheno/))
7. **Quality control reporting** - Compile comprehensive QC report with alignment and variant metrics ([`MULTIQC`](modules/nf-core/multiqc/))

The pipeline uses NG_006669.2 (RefSeqGene, LRG_792) as a single combined ABO reference. Variant marker positions and coordinate calibration are defined in `assets/refs/abo_variant_panel.yaml` (see [`bin/README.md`](bin/README.md) for details).

## Summary of tools and version used in the pipeline

The pipeline makes use of the following core dependencies:

| Dependency             | Minimum version |
| ---------------------- | --------------- |
| **Core Tools**         |                 |
| fastqc                 | 0.12.1          |
| minimap2               | 2.29-r1283      |
| samtools               | 1.21            |
| multiqc                | 1.30            |
| **Python Environment** |                 |
| python                 | 3.13.5          |
| pandas                 | 2.3.1           |
| numpy                  | 2.3.2           |
| xlsxwriter             | 3.2.5           |
| json5                  | 0.12.0          |
| openpyxl               | 3.1.2           |
| re                     | 2.2.1           |

## Required input files structure

Ensure that all input fastq files have a naming convention that matches this regular expression (`regex`)

```python
## python regex for matching samples
pattern = r"^(IMM|INGS|NGS|[A-Z0-9]+)(-[0-9]+-[0-9]+)?_barcode\d+$"
```

The regex does the following:

- `^(IMM|INGS|NGS|[A-Z0-9]+)` allows for files strating with the prefixes IMM, INGS, NGS, or any combination of letters `A-to-Z` and digits `0-to-9`.
- `(-[0-9]+-[0-9]+)?` handles optional segments of digits separated by a dash(-).
- `_barcode\d+$` ensures the filename ends with_barcode followed by digits to denote barcode numbers.

There is a file handling logic in the code `filename.split("_")` that assumes the barcode is always the last part of the filename.
The names are split into `basename` and `barcode` which are then used in later reporting.<br/>Please Adjust this if necessary based on actual filename structure in your assays.

Here are a few examples of acceptable input file names:

```txt
NGSPOS_barcode13.fastq
NGSNEG_barcode12.fastq
INGSPOS_barcode01.fastq
INGSNEG_barcode96.fastq
BTGSPOS_barcode19.fastq
2025705_barcode14.fastq
IMM-45-44874_barcode25.fastq
Sample1-2024-12345_barcode22.fastq
```

## Sequencing platform compatibility

This pipeline was originally developed to process amplicon sequencing data from Oxford Nanopore Technologies platforms where multiple samples are expected to be barcoded. The pipeline has been extensively validated using Oxford Nanopore MinION data targeting ABO exons 6 and 7, which are the primary regions containing clinically relevant polymorphisms for ABO blood group determination.

While we recommend using the above naming convention for optimal compatibility, the pipeline can also handle FASTQ files from other sequencing platforms including:

- **PacBio** (currently undergoing testing)
- **Ion Torrent** (currently undergoing testing)
- **Illumina** (currently undergoing testing)

The pipeline will attempt to extract the sample name and barcode from the filenames using standard genomic sequence naming conventions, but will fall back to a default barcode00 if filenames lack the expected barcode format. For non-Nanopore platforms, ensure your FASTQ files contain reads spanning the ABO exon 6 and exon 7 regions (within the combined NG_006669.2 reference) for accurate genotyping.

## Running `nf-core/abotyper`

This pipeline has been extensively tested using conda, docker, and singularity profiles. Other containerisation methods are being improved,tested and documented. It is adisable to run a minimal test run to check your environment before submitting big jobs.

To run this pipeline, use:

```bash
# Quick test
nextflow run nf-core-abotyper/main.nf --outdir test_working -profile test,conda -resume

# Complete run
nextflow nf-core/abotyper \
  -resume \
  -profile <docker/singularity/.../institute> \
  --input samplesheet.csv \
  --outdir <OUTDIR>
```

> [!WARNING]
> Please provide pipeline parameters via the CLI or Nextflow `-params-file` option. Custom config files including those provided by the `-c` Nextflow option can be used to provide any configuration _**except for parameters**_; see [docs](https://nf-co.re/docs/running/run-pipelines#using-parameter-files).

The code by permits renaming of samples using a tab-delimited file with `sequencingID` and `sampleName` (see `nextflow.config` file under `$params.renaming_file`).
This option is controlled by the parameter `$params.skip_renaming` and can be overridden via the commandline using option `--skip_renaming true` to skip the process.

## Output

For each sample, the pipeline aligns reads once against a single combined ABO reference (NG_006669.2, spanning exons 2-7) and generates `BAM` files, `BAM metrics`, and per-position variant statistics.

The output directory generated by this `Nextflow` pipeline will look something like this:

```
OUTDIR/
├── ABO_results.log
├── ABO_result.txt
├── ABO_result.xlsx
├── final_export.csv
├── per_sample_processing
│   └── SAMPLE1_barcode01
│       └── combined
│           ├── ABOReadPolymorphisms.txt
│           ├── alignment
│           │   ├── SAMPLE1_barcode01.bam
│           │   ├── SAMPLE1_barcode01.bam.bai
│           │   ├── SAMPLE1_barcode01.coverage.txt
│           │   ├── SAMPLE1_barcode01.flagstat
│           │   └── SAMPLE1_barcode01.stats
│           ├── SAMPLE1_barcode01.ABOPhenotype.txt
│           ├── SAMPLE1_barcode01.AlignmentStatistics.tsv
│           ├── SAMPLE1_barcode01.Haplotypes.tsv
│           └── SAMPLE1_barcode01.log.txt
├── pipeline_info
│   ├── execution_report_DATETIME.html
│   ├── execution_timeline_DATETIME.html
│   ├── execution_trace_DATETIME.txt
│   ├── nf_core_pipeline_software_mqc_versions.yml
│   ├── params_DATETIME.json
│   └── pipeline_dag_DATETIME.html
└── qc-reports
    ├── fastqc
    │   ├── SAMPLE1_barcode01_fastqc.html
    │   ├── SAMPLE1_barcode01_fastqc.zip
    └── multiqc
        ├── multiqc_data
        ├── multiqc_plots
        │   ├── pdf
        │   ├── png
        │   └── svg
        └── multiqc_report.html
```

The `ABO_result.xlsx` Excel worksheet contains details of all SNVs and metrics used to deduce the ABO phenotype for each sample.

A summary of the ABO typing results is provided in `final_export.csv`

Feel free to raise an issue or reach out if you need any support getting this tool running, or with suggestions for improvement.

## Credits

nf-core/abotyper was originally written by Fredrick M. Mobegi: [@fmobegi](https://github.com/fmobegi) at the Department of Clinical Immunology, [PathWest Laboratory Medicine WA](https://pathwest.health.wa.gov.au/).

We thank the following people for their extensive assistance in the development and testing of this pipeline:

- [Benedict Matern](https://github.com/bmatern)
- [Mathijs Groeneweg](https://orcid.org/0000-0002-6615-9239)
- [Filipe Ayora](https://github.com/fayora)

Maintenance and future developements will be led by Fredrick Mobegi.

## Acknowledgements

<p float="center">
  <img src = "docs/images/pathwest_logo.png", width="400", height="90"/>
  <img src = "docs/images/uwa_logo.png", width="400", height="90">
</p>

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](docs/CONTRIBUTING.md).

For further information or help, don't hesitate to get in touch on the [Slack `#abotyper` channel](https://nfcore.slack.com/channels/abotyper) (you can join with [this invite](https://nf-co.re/join/slack)).

## Further reading

Results generated from this pipeline should be interpreted together with the corresponding publication and literature on ABO genotyping.

To get up to speed with ABO genotyping, there is detailed reading material [here](https://ftp.ncbi.nlm.nih.gov/pub/mhc/rbc/Final%20Archive/Excel_and_PowerPoint/).

Verified SNVs relevant to ABO blood group genotyping have also been documented extensively [here](https://bloodgroupdatabase.org/groups/details/?group_name=ABO)

## Citations

If you use nf-core/abotyper for your analysis, please cite it using the following publication:

> **Characterisation of the ABO Blood Group Phenotypes Using Third-Generation Sequencing.**
>
> Fredrick M. Mobegi, Samuel Bruce, Naser El-Lagta, Felipe Ayora, Benedict M. Matern, Mathijs Groeneweg, Lloyd J. D'Orsogna & Dianne De Santis.
>
> _Int. J. Mol. Sci._ 2025 Jun 06. doi: [10.3390/ijms26125443](https://doi.org/10.3390/ijms26125443).

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

You can cite the `nf-core` publication as follows:

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
>
> _Int. J. Mol. Sci._ 2025 Jun 06. doi: [10.3390/ijms26125443](https://doi.org/10.3390/ijms26125443).

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

You can cite the `nf-core` publication as follows:

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
