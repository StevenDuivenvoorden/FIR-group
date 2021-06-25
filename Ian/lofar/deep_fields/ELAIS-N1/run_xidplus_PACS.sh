#!/bin/bash


#$ -cwd

#$ -o /lustre/scratch/astro/im281/FIR-group/Ian/lofar/deep_fields/ELAIS-N1/log/
#$ -e /lustre/scratch/astro/im281/FIR-group/Ian/lofar/deep_fields/ELAIS-N1/log/

#sumbit with the following command
# qsub -t 1-10 -l m_mem_free=10G run_xidplus_PACS.sh

. /etc/profile.d/modules.sh
module load sge

module load easybuild/software
module load python/intelpython3/3.5.3
#module load Anaconda3/4.0.0

source activate herschelhelp_v2

cd /lustre/scratch/astro/im281/FIR-group/Ian/lofar/deep_fields/ELAIS-N1
echo running_python

python XID+_run_PACS.py
echo finished


