#!/bin/bash


#$ -cwd

#$ -o /lustre/scratch/astro/im281/FIR-group/Ian/lofar/deep_fields/ELAIS-N1/log/out
#$ -e /lustre/scratch/astro/im281/FIR-group/Ian/lofar/deep_fields/ELAIS-N1/log/err

#sumbit with the following command
# qsub -t 1-10 -l m_mem_free=20G run_xidplus_MIPS.sh

. /etc/profile.d/modules.sh
module load sge

module load easybuild/software
module load python/intelpython3/3.5.3
#module load Anaconda3/4.0.0
#echo conda_activate
source activate /its/home/im281/.conda/envs/herschelhelp
#conda activate herschelhelp
#conda list | grep webencoding >> file
cd /lustre/scratch/astro/im281/FIR-group/Ian/lofar/deep_fields/ELAIS-N1
#export PATH="/its/home/im281/.conda/envs/herschelhelp/bin/python":$PATH
echo running_python
#/its/home/im281/.conda/envs/herschelhelp/bin/python3.6 XID+_run_ELAIS-N1_apollo.py
/its/home/im281/.conda/envs/herschelhelp/bin/python3.6 XID+_run_MIPS.py
echo finished
