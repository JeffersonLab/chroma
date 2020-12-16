#! /bin/bash
#####################################################
# Trivial Runner Script 
# Suitable for running on scalar machines or places
# where no special invocation is needed for running
# parallel code
####################################################
export OMP_NUM_THREADS=1
mpirun -genvall -n 1 $*

