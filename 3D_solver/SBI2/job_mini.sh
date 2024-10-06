#!/bin/bash
#PJM -L rscgrp=tutorial1-a
#PJM -L node=1
#PJM --mpi proc=4
#PJM -L elapse=00:05:00
#PJM -g gt01

module load nvidia/22.7 cuda/11.4 ompi-cuda

mpiexec -machinefile $PJM_O_NODEINF -n \
$PJM_MPI_PROC -npernode 4 ./a.out

