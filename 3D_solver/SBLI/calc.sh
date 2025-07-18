export UCX_MEMTYPE_CACHE=n
export UCX_IB_GPU_DIRECT_RDMA=y
export OMPI_MCA_opal_cuda_support=true

nohup mpiexec -n 4 ./a.out &
cp mod_globals.f90 ./data
cp set.f90 ./data

