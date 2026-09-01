# HLA association analysis using Regenie (Nextflow pipeline)

[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A524.04.2-brightgreen.svg)](https://www.nextflow.io/)
[![Regenie](https://img.shields.io/badge/regenie-4.0-blue.svg)](https://rgcgithub.github.io/regenie/)

## Introduction

The pipeline runs association analysis of imputed HLA alleles, amino acid residues and SNPs in the MHC region against binary phenotypes using Regenie.
See the [Paper](https://www.nature.com/articles/s41588-021-00870-7) and the [GitHub](https://github.com/rgcgithub/regenie) repo.

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. Each phenotype is analysed in parallel, and the results are plotted, cleaned and merged into a single table at the end.

The pipeline takes in PLINK binary files (`.bed/.bim/.fam`) of imputed HLA variants as input, together with a phenotype file and a covariate file in Regenie format, and a plain-text list of the phenotype columns to analyse.

### Pipeline steps

| Step | Process | Tool | Description |
| ---- | ------- | ---- | ----------- |
| 1 | `step_1` | Regenie | Whole-genome ridge regression (`--step 1`, `--bt`, `--cv 3`, `--lowmem`) to fit the null model and produce the `*_pred.list` predictions for each phenotype |
| 2 | `step_2` | Regenie | Association testing of each HLA variant (`--step 2`, `--bt`, `--minMAC 1`) with Firth correction, conditioned on the step 1 predictions |
| 3 | `plot_manhattan` | R (ggplot2) | Manhattan plot of the MHC region, coloured by classical HLA gene and shaped by variant type (HLA allele, amino acid, SNP, rsID) |
| 4 | `top_alleles` | R | Cleans each step 2 output (`CHROM, GENPOS, ID, BETA, SE, LOG10P, A1FREQ, ALLELE0, ALLELE1`) and tags it with the phenotype name |
| 5 | `combineFiles` | R | Row-binds the cleaned results of all phenotypes into one `combined_file.tsv` |

Processes for a two-round conditional analysis (extracting the top allele, adding its dosage as a covariate and re-running Regenie) are drafted in `templates/conditional_analysis.R` and commented out in `main.nf`; they are not yet wired into the workflow.

## Installation

1. Nextflow

```
wget -qO- https://get.nextflow.io | bash
```

2. Regenie

```
conda install -c bioconda regenie
```

or load it as a module on your HPC, e.g.

```
module load Regenie/4.0-GCC-12.3.0
```

3. R (≥ 4.2) with the packages used by the plotting and cleaning templates

```
install.packages(c("pacman", "tidyverse", "ggplot2", "dplyr", "plyr", "vroom", "readr", "qqman", "gridExtra", "RColorBrewer"))
```

The `pacman` call at the top of each template installs any missing package automatically on first run.

## Input files

### Genotypes

PLINK binary fileset of imputed HLA variants for the MHC region. Variant IDs in the `.bim` file should follow the usual HLA imputation convention (`HLA_*`, `AA_*`, `SNP_*`, `rs*`), as the Manhattan plot uses the ID prefix to distinguish variant types.

### Phenotype and covariate files

Both files follow the Regenie format: whitespace-delimited, with `FID` and `IID` as the first two columns, one column per phenotype/covariate and `NA` for missing values.

```
FID IID T2D CAD Asthma
1   1   0   1   NA
2   2   1   0   0
```

Binary phenotypes must be coded 0/1.

### Phenotype list

A plain-text file with one phenotype column name per line. Each line becomes an independent Regenie run.

```
T2D
CAD
Asthma
```

## Running the pipeline

The pipeline does not require installation as `Nextflow` will automatically fetch it from `GitHub`.

### Configure your data

Copy `conf/test.config` and edit the paths to point to your files:

```groovy
def cov_pheno_file = "/path/to/phenotype_list.txt"
def cov_pheno_list = []

new File(cov_pheno_file).eachLine { line ->
    cov_pheno_list.add([line.trim(), "/path/to/phenotypes.tsv", "/path/to/covariates.tsv"])
}

params {
    cov_pheno = cov_pheno_list

    plink_files = [
        ['MyCohort',
         '/path/to/hla_imputed.bed',
         '/path/to/hla_imputed.bim',
         '/path/to/hla_imputed.fam'],
    ]

    outdir = "./output"
}
```

`cov_pheno` is a list of `[phenotype_name, phenotype_file, covariate_file]` and `plink_files` is a list of `[dataset_name, bed, bim, fam]`. Several datasets can be listed; each will be run against every phenotype.

### Run on a Slurm cluster

```
nextflow run nanjalaruth/HLA-association-analysis-with-Regenie -profile slurm -c <path to your edited config file> -resume
```

A submission script is provided in `run.sh`, which loads the required modules and launches the pipeline:

```
sbatch run.sh
```

### Run locally

```
nextflow run nanjalaruth/HLA-association-analysis-with-Regenie -profile standard -c <path to your edited config file> -resume
```

## To run the updated version of this pipeline, run:

```
nextflow pull nanjalaruth/HLA-association-analysis-with-Regenie
```

## Arguments

### Required Arguments

| Argument       | Usage                                              | Description |
| -------------- | -------------------------------------------------- | ----------- |
| -profile       | <standard,slurm,singularity,conda,debug>           | Configuration profile to use. `standard` runs locally; `slurm` submits each process to the scheduler |
| -c             | <path to config>                                   | Config file defining `cov_pheno` and `plink_files` (see above) |
| --cov_pheno    | [[pheno, pheno_file, covar_file], ...]             | List of phenotypes with their phenotype and covariate files |
| --plink_files  | [[dataset, bed, bim, fam], ...]                    | PLINK binary fileset(s) of imputed HLA variants |

### Optional Arguments

| Argument       | Default                        | Description |
| -------------- | ------------------------------ | ----------- |
| --outdir       | `./output`                     | Directory for published results |
| --tracedir     | `${outdir}/pipeline_info`      | Directory for the Nextflow timeline, report, trace and DAG |
| --max_memory   | `270.GB`                       | Upper bound for process memory |
| --max_cpus     | `16`                           | Upper bound for process CPUs |
| --max_time     | `240.h`                        | Upper bound for process wall time |

Regenie step 2 is labelled `bigmem` (120 GB, 9 CPUs, 24 h, doubled on retry). Adjust the `withLabel` blocks in `nextflow.config` if your cluster limits differ.

## Output

```
output/
├── Regression_results/
│   ├── <dataset>_<pheno>_regenie-step1_pred.list
│   ├── <dataset>_<pheno>_regenie-step1_*.loco
│   └── <dataset>_<pheno>_regenie-step2-bin-out-firth_<pheno>.regenie
├── plots/
│   └── <dataset>_<pheno>_phewas.pdf
├── top_alleles/
│   └── clean_<pheno>.tsv
├── merged/
│   └── combined_file.tsv
└── pipeline_info/
    ├── execution_timeline.html
    ├── execution_report.html
    ├── execution_trace.txt
    └── pipeline_dag.png
```

`combined_file.tsv` contains the association statistics for every variant across all phenotypes and is the main table for downstream work (e.g. filtering by `LOG10P` or extracting classical HLA alleles with `grep '^HLA'`).

## Notes

- The gene boundaries used to colour the Manhattan plot (`templates/manhattan.R`) are GRCh37 coordinates. Update them if your genotypes are on GRCh38.
- The R templates assume `Rscript` is on the `PATH`; edit the shebang line if your environment differs.

## Support

I track open tasks using github's [issues](https://github.com/nanjalaruth/HLA-association-analysis-with-Regenie/issues)
