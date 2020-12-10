#!/bin/bash

############################################################################
# Definition of variables
############################################################################
EXE=lmp
Partitions=2
coresPerPartition=2
totalCores=`echo ${Partitions}*${coresPerPartition} | bc -l`
############################################################################

mpirun -np ${totalCores} ${EXE} -p ${Partitions}x${coresPerPartition} -in start.lmp > out.lmp
