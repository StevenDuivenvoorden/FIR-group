#!/bin/bash


#$ -cwd

#$ -o /lustre/scratch/astro/im281/FIR-group/Ian/lofar/deep_fields/Lockman/log/SPIRE/
#$ -e /lustre/scratch/astro/im281/FIR-group/Ian/lofar/deep_fields/Lockman/log/SPIRE/

# sumbit with the following command
# qsub -t 1-10 -l m_mem_free=20G run_xidplus_elais_n1.sh
echo starting script
. /etc/profile.d/modules.sh
module load sge

#module load easybuild/software
#module load python/intelpython3/3.5.3
# module load Anaconda3/4.0.0
echo conda_activate
#source activate /its/home/im281/.conda/envs/herschelhelp

cd /lustre/scratch/astro/im281/FIR-group/Ian/lofar/deep_fields/Lockman
# export PATH="/its/home/im281/.conda/envs/herschelhelp/bin/python":$PATH

export PATH=$PATH':/research/astro/fir/HELP/help_python/miniconda3/bin'

echo running_python
python3 XID+_run_SPIRE.py
#/its/home/im281/.conda/envs/herschelhelp/bin/python3.6 XID+_run_SPIRE.py
# /its/home/im281/.conda/envs/herschelhelp/bin/python3.6 XID+_run_ELAIS-N1_PACS.py
echo finished

