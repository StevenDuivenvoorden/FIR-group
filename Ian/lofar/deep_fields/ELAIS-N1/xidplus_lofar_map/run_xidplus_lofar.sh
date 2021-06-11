#!/bin/bash

# -- The shell used to interpret this script
#$ -S /bin/bash
# -- Execute this job from the current working directory.
#$ -cwd
#$ -q smp.q@@lfs212
# -- Allocate more memory
#$ -l m_mem_free=5G
# -- Job output to stderr will be merged into standard out. Remove this line if
# -- you want to have separate stderr and stdout log files
#$ -j y
#$ -o output/
#$ -cwd


#sumbit with the following command
# qsub -t 1-10 -l m_mem_free=20G run_xidplus_lofar.sh

. /etc/profile.d/modules.sh
module load sge
module load easybuild/software
module load Anaconda3/2020.02

source activate /its/home/im281/.conda/envs/herschelhelp_v2

cd /mnt/lustre/projects/astro/general/im281/FIR-group/Ian/lofar/deep_fields/ELAIS-N1/xidplus_lofar_map
which python
echo running_python

~/.conda/envs/herschelhelp_v2/bin/python run_xidplus_lofar_random_positions.py
echo finished

