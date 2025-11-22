#!/bin/bash
ID=${OMPI_COMM_WORLD_RANK}

export ROCR_VISIBLE_DEVICES="${OMPI_COMM_WORLD_LOCAL_RANK}"
./mega2025_open4amd_kti-8mpi.go 1> stdout.$ID 2> stderr.$ID
