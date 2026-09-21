#!/bin/bash
#SBATCH -J CALCULATION_NAME
#SBATCH -o job-%J.o
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=8
#SBATCH -A general
#SBATCH -N 1
#SBATCH -t 4:00:00
#SBATCH --mem-per-cpu=5G

export JOB=CALCULATION_NAME
export DIR=$SLURM_SUBMIT_DIR
export scratch=$SCRATCH/crys17/

echo "submit directory: "
echo $SLURM_SUBMIT_DIR

ml -* CRYSTAL/17

mkdir  -p $scratch/$JOB
cp $DIR/$JOB.d12  $scratch/$JOB/INPUT
cd $scratch/$JOB

srun Pcrystal 2>&1 >& $DIR/${JOB}.out
cp fort.9 ${DIR}/${JOB}.f9
