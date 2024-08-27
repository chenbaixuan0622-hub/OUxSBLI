rm ./data/*.d
rm ./data/*.vtr
rm ./data/1d/*.d
nohup mpiexec -n 2 ./a.out &
cp mod_globals.f90 ./data
cp set.f90 ./data

