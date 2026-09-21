#!/bin/sh 
#SBATCH -o daemon_log
#SBATCH -J USPEX_CALCULATION_NAME_DAEMON
#SBATCH -t 1-00:00:00
#SBATCH -N 1
#SBATCH -n 1

##module load gnu-openmpi

while [ ! -f ./USPEX_IS_DONE ]; do
   date >> log 
   USPEX -r >> log
   sleep 10
done

 
