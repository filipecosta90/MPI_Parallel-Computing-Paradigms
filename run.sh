#!/bin/sh#
mpirun -np 1 valgrind --track-origins=yes --leak-check=full bin/pcp_tp2_mpi 4 128 128
