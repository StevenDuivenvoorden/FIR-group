#!/bin/bash


#$ -cwd

#$ -o /lustre/scratch/astro/im281/FIR-group/Ian/lofar/deep_fields/Bootes/log/
#$ -e /lustre/scratch/astro/im281/FIR-group/Ian/lofar/deep_fields/Bootes/log/

# sumbit with the following command
# qsub -t 1-10 -l m_mem_free=20G run_xidplus_elais_n1.sh
#echo starting script
#. /etc/profile.d/modules.sh
#module load sge

#module load easybuild/software
#module load python/intelpython3/3.5.3
# module load Anaconda3/4.0.0
#echo conda_activate
#source activate herschelhelp_v2

#cd /lustre/scratch/astro/im281/FIR-group/Ian/lofar/deep_fields/Bootes
#echo running_python
#/its/home/im281/.conda/envs/herschelhelp/bin/python3.6 XID+_run_PACS.py
#echo finished



. /etc/profile.d/modules.sh
module load sge

module load easybuild/software
module load python/intelpython3/3.5.3
#module load Anaconda3/4.0.0

source activate herschelhelp_v2

cd /lustre/scratch/astro/im281/FIR-group/Ian/lofar/deep_fields/Bootes
echo running_python

python XID+_run_PACS.py
echo finished