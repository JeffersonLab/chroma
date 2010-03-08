#!/bin/bash

##########################################################
# Sanity Check - want to potentially run this from PBS
# Interactive Queue: 
#
# Assumptions - I am in the working directory
#               PBS_NODEFILE is defined
##########################################################

# Hardwire MPI for now
MPI_HOME=/usr/local/mvapich-0.9.9
NODES="qcd7n1411 qcd7n1411 qcd7n0122 qcd7n0122"
# Nex is the name of the program
PROG=$1

# Now get the args
shift
ARGS=$*

echo Want to run program: $PROG with args $ARGS on nodes $NODES
# -----------------------------------------------
# Internal configuration
# -----------------------------------------------

# Count nodes
MPIRUN=${MPI_HOME}/bin/mpirun_rsh
MPIRUN_ARGS="-rsh -np 4 ${NODES} ${PROG} ${ARGS}"

# Run program
${MPIRUN} ${MPIRUN_ARGS}



