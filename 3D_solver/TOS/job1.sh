#!/bin/bash
#PJM -L rscgrp=tutorial-share
#PJM -L gpu=1
#PJM --mpi proc=2
#PJM -L elapse=00:60:00
#PJM -g gt01

module load nvidia/24.1 nvmpi
mpiexec -machinefile $PJM_O_NODEINF -n \
$PJM_MPI_PROC -npernode 2 ./a.out

