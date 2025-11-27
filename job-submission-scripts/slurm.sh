#!/usr/bin/env bash

## timelimit in min or min:sec or h:min:sec or d-h or d-h:min or d-h:min:sec
#SBATCH -t 30:00:00

## multiple executions of the same command, match total number of tasks accordingly, all number commands below will apply to each job in the array
## !!! use this to set the number of events !!!
#SBATCH --array=1-128
## number of nodes
#SBATCH --nodes=1
## number of tasks, usually number of MPI ranks, number of processes to be run
# #SBATCH --ntasks=1
## number of cpus per task, actually number of cores per task
#SBATCH --cpus-per-task=8
## amount of memory per cpu (core)
#SBATCH --mem-per-cpu 2G
## number of tasks per node
#SBATCH --ntasks-per-node=1


## quality of service, basically priority level, test express fpgasynthesis
# #SBATCH -q

## partition to use, normal or gpu
#SBATCH -p normal

## mail of user
## SBATCH --mail-user example@email.com
## mail type: NONE BEGIN END FAIL REQUEUE ALL
#SBATCH --mail-type END

## do not kill job entire job if one task fails
#SBATCH --no-kill

## do not requeue if job fails
#SBATCH --no-requeue

## output location of .out files
#SBATCH --output=slurm/%A_%a.out
## output location of .err files
# #SBATCH --error=slurm/%A_%a.err

## -the slum system might need certain includes to run
## -compiling statically might make these unnecessary
# ml compiler/GCC/11.3.0
# ml numlib/GSL/2.7-GCC-11.3.0

## command to be run (remove everything besides the eic command to actually run)
cat <<- EOF
	Would run:
	"./eic-sde -A 1 -H 3 -rH2 0.7 -p --add-to-seed ${SLURM_ARRAY_TASK_ID}" -t $SLURM_CPUS_PER_TASK"
EOF
exit 1
