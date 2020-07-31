#!/bin/bash


#$ -cwd

#$ -o /lustre/scratch/astro/im281/FIR-group/Ian/lofar/deep_fields/ELAIS-N1/log/
#$ -e /lustre/scratch/astro/im281/FIR-group/Ian/lofar/deep_fields/ELAIS-N1/log/

# sumbit with the following command
# qsub -t 1-10 -l m_mem_free=10G run_xidplus_SPIRE.sh
#echo starting script
#. /etc/profile.d/modules.sh
#module load sge

#module load easybuild/software
#module load python/intelpython3/3.5.3
# module load Anaconda3/4.0.0

source ~/.bashrc

source activate /its/home/im281/.conda/envs/herschelhelp
#source activate hhelp

cd /lustre/scratch/astro/im281/FIR-group/Ian/lofar/deep_fields/ELAIS-N1


echo running_python
python XID+_run_SPIRE.py
#/its/home/im281/.conda/envs/herschelhelp/bin/python3.6 XID+_run_SPIRE.py

echo finished
