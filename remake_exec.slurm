#!/bin/bash
# Created by the University of Melbourne job script generator for SLURM
# Thu Oct 22 2020 00:07:31 GMT+1100 (Australian Eastern Daylight Time)

# Partition for the job:
#SBATCH --partition=long

# Multithreaded (SMP) job: must run on one node 
#SBATCH --nodes=1

# The name of the job:
#SBATCH --job-name="WD_elec_cap"

# The project ID which this job should run under:
#SBATCH --account="punim1360"

# Maximum number of tasks/CPU cores used by the job:
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1

#SBATCH -o ./spart_output/make_out.out # STDOUT 

# The maximum running time of the job in days-hours:mins:sec
#SBATCH --time=0-0:5:00

# check that the script is launched with sbatch
if [ "x$SLURM_JOB_ID" == "x" ]; then
   echo "You need to submit your job to the queuing system with sbatch"
   exit 1
fi

# Run the job from the directory where it was launched (default)

# The modules to load:
module purge
module load gcc/8.3.0
module load openmpi/3.1.4
module load gsl/2.5


# The job command(s):
make clean
make
