#!/bin/bash
#SBATCH -A p31628               # Allocation
#SBATCH -p long                # Queue
#SBATCH -t 70:00:00             # Walltime/duration of the job
#SBATCH -N 5                   # Number of Nodes
#SBATCH --mem=150G               # Memory per node in GB needed for a job. Also see --mem-per-cpu
#SBATCH --ntasks-per-node=4     # Number of Cores (Processors)
#SBATCH --mail-user=zidanwang2025@u.northwestern.edu  # Designate email address for job communications
#SBATCH --mail-type=END         # Events options are job BEGIN, END, NONE, FAIL, REQUEUE
#SBATCH --job-name="zidanwang_tca"  # Name of job

# unload any modules that carried over from your command line session
module purge
module spider R
module load R/4.1.1

# load modules you need to use
# A command you actually want to execute:
Rscript TCA_generate.R
