#!/bin/bash

#SBATCH --job-name='regression analysis'
#SBATCH --output=reg-%j-stdout.log
#SBATCH --error=reg-%j-stderr.log

#Run command

module purge
module load Nextflow/24.04.2
module load BCFtools/1.18-GCC-12.3.0 
module load Regenie/4.0-GCC-12.3.0

nextflow run main.nf -resume -profile slurm -c conf/test.config

