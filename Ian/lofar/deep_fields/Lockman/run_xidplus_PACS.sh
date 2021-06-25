#!/bin/bash
#$ -S /bin/bash

#$ -cwd

#$ -o /lustre/scratch/astro/im281/FIR-group/Ian/lofar/deep_fields/Lockman/log/
#$ -e /lustre/scratch/astro/im281/FIR-group/Ian/lofar/deep_fields/Lockman/log/

# sumbit with the following command
# qsub -t 1-10 -l m_mem_free=20G run_xidplus_elais_n1.sh
echo starting script
. /etc/profile.d/modules.sh
module load sge

module load easybuild/software
module load python/intelpython3/3.5.3
# module load Anaconda3/4.0.0
echo conda_activate
source activate /its/home/im281/.conda/envs/herschelhelp_v2

cd /lustre/scratch/astro/im281/FIR-group/Ian/lofar/deep_fields/Lockman
# export PATH="/its/home/im281/.conda/envs/herschelhelp/bin/python":$PATH
echo $SGE_TASK_ID
echo running_python

#export PATH=$PATH':/research/astro/fir/HELP/help_python/miniconda3/bin'

#python3 XID+_run_PACS.py
/its/home/im281/.conda/envs/herschelhelp/bin/python3.6 XID+_run_PACS.py
# /its/home/im281/.conda/envs/herschelhelp/bin/python3.6 XID+_run_ELAIS-N1_PACS.py
echo finished

