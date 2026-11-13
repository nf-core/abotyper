# Conda lock files

These files pin the conda environments of the nf-core modules that ship lock files (bcftools, samtools and MultiQC). Use them with the `conda_lock` profile:

```bash
nextflow run nf-core/abotyper -profile conda,conda_lock ...
```

`conda_lock` loads `conf/conda_lock_explicit_amd64.config` or `conf/conda_lock_explicit_arm64.config`, depending on the machine Nextflow runs on, and those point here.

## Why they exist

nf-core tools generates `conf/containers_conda_lock_files_{amd64,arm64}.config`, which point at the lock files in each module (`modules/nf-core/<module>/.conda-lock/*.txt`). Those are pixi locks (`version: 6`), which conda cannot read: Nextflow 25.10 passes them to `conda create --file` and environment creation fails, and Nextflow 26.10 does not accept them either. The MultiQC lock is also invalid YAML in nf-core/modules. The generated configs are left as nf-core tools makes them, because `nf-core pipelines lint` checks that they are unchanged.

Each file here is a conda `@EXPLICIT` list (package URL and md5) with the same packages as the module lock. Nothing is solved again.

## After a module update

Regenerate from the pipeline root:

```bash
python assets/conda-lock/make_explicit_locks.py
```

The script reads the module lock paths from the generated configs, writes the files here, rewrites `conf/conda_lock_explicit_{amd64,arm64}.config` and removes files no longer used.

Once nf-core tools ships lock files that conda can read, delete this folder, the two `conf/conda_lock_explicit_*` configs and the `conda_lock` profile.
