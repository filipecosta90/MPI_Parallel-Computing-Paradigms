#!/bin/sh#
mpiexec -n 1 -mca btl self,sm,tcp bin/pcp_tp2_mpi 8192 1 1 
