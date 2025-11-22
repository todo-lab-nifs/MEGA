#!/bin/bash

#PBS -N opn020
#PBS -P NIFS24KISC014
#PBS -q A_S
#PBS -l select=4:ncpus=256:mpiprocs=256:mem=752gb:ompthreads=1
#PBS -l walltime=10:00:00

# export F_UFMTENDIAN=60
export FORT1='./opn020_001.in'
export FORT7='/home/todo/MEGA_OPEN/opn020/opn020_001.out'
export FORT30='./opn020_001.moments'
export FORT31='./opn020_001.harmonics'
export FORT40='/home/todo/MEGA_OPEN/opn020/opn020_001.energy_phys.txt'
export FORT41='/home/todo/MEGA_OPEN/opn020/opn020_001.energy_n.txt'
export FORT60='/data/todo/241216.eql062.lendian.d'

  mkdir -p /data/todo/opn020
  cd /data/todo/opn020
  cp /home/todo/MEGA_OPEN/mega2025_open4_kti-1024mpi.go .
  cp /home/todo/MEGA_OPEN/opn020/opn020_001.in .

  module load intel/2025.1
  export I_MPI_DEBUG=4
  export I_MPI_OFI_PROVIDER=psm3

  mpirun -np 1024 ./mega2025_open4_kti-1024mpi.go
