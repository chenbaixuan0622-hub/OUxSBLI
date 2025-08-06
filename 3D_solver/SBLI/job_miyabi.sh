#!/bin/bash
#PBS -q tutorial1-g
#PBS -l select=1:mpiprocs=4
#PBS -l walltime=00:60:00
#PBS -W group_list=gt01
#PBS -j oe

module purge
module load nvidia/24.9 nv-hpcx/24.9

cd ${PBS_O_WORKDIR}
mpiexec ./a.out

exit 0

