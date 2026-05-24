#!/bin/bash
mkdir -p data/1 data/3 data
nohup mpiexec -n 4 ./a.out &
cp mod_globals.f90 ./data
cp set.f90 ./data
