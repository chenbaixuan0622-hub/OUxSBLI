#!/bin/bash
#PBS -q tutorial1-g
#PBS -l select=1
#PBS -l walltime=00:60:00
#PBS -W group_list=gt01
#PBS -j oe

module purge
module load nvidia/24.9 nv-hpcx/24.9 python/3.10.16

export CUPY_CUDA_PATH=/work/opt/local/aarch64/cores/nvidia/24.9/Linux_aarch64/24.9
export PATH=$CUPY_CUDA_PATH/compilers/bin:$PATH
export CUPY_NVCC_FLAGS="-I$CUPY_CUDA_PATH/include"

source venv/bin/activate

cd ${PBS_O_WORKDIR}
python3 main.py

exit 0

