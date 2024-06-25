nohup nsys profile -t cuda,nvtx,openacc,osrt -f true -o my_report_shared mpirun -np 4 ./a.out &

