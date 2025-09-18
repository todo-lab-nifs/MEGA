#!/bin/bash
ID=${OMPI_COMM_WORLD_RANK}

export ROCR_VISIBLE_DEVICES="${OMPI_COMM_WORLD_LOCAL_RANK}"
# ./mega2025_open3amd-4mpi.go
./mega2025_open3amd-4mpi.go 1> stdout.$ID 2> stderr.$ID
