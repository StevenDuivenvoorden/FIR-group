#!/bin/bash


#$ -cwd

#$ -o /lustre/scratch/astro/im281/FIR-group/Ian/lofar/deep_fields/ELAIS-N1/log/out
#$ -e /lustre/scratch/astro/im281/FIR-group/Ian/lofar/deep_fields/ELAIS-N1/log/err

#sumbit with the following command
# qsub -t 1-10 -l m_mem_free=10G run_xidplus_PACS.sh

. /etc/profile.d/modules.sh
module load sge

module load easybuild/software
module load python/intelpython3/3.5.3
#module load Anaconda3/4.0.0

#source activate /its/home/im281/.conda/envs/herschelhelp
source activate herschelhelp

cd /lustre/scratch/astro/im281/FIR-group/Ian/lofar/deep_fields/ELAIS-N1

echo running_python

#/its/home/im281/.conda/envs/herschelhelp/bin/python3.6 XID+_run_PACS.py
python XID+_run_PACS.py
echo finished
