mpicc -o task2 ./mpitsp.c -lm
mpirun -np 7 ./task2 3500 2 chain best 20000