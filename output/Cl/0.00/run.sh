#!/bin/bash

#SBATCH --job-name=Cl-mag-fields-gto
#SBATCH -e out_%j
#SBATCH -o out_%j
#SBATCH --mem-per-cpu=300
#SBATCH -t 01:00:00
#SBATCH -n 8
#SBATCH -p normal

erkale_complex_orbs_omp HGBSP1-5_8.run
