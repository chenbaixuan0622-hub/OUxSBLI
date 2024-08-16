nohup mpirun -n 4 a.out &

rm *.mod *.o

cp mod_globals.f90 ./data
cp set.f90 ./data

