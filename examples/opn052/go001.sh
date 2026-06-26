#!/bin/bash

#PBS -N opn052
#PBS -P NIFS26KISC031
#PBS -q A_S
#PBS -l select=4:ncpus=256:mpiprocs=256:mem=752gb:ompthreads=1
#PBS -l walltime=10:00:00

# export F_UFMTENDIAN=60
export FORT1='./opn052_001.in'
export FORT7='/home/todo/MEGA_OPEN/opn052/opn052_001.out'
export FORT30='./opn052_001.moments'
export FORT31='./opn052_001.harmonics'
export FORT40='/home/todo/MEGA_OPEN/opn052/opn052_001.energy_phys.txt'
export FORT41='/home/todo/MEGA_OPEN/opn052/opn052_001.energy_n.txt'
export FORT60='/data/todo/241216.eql062.lendian.d'

  mkdir -p /data/todo/opn052
  cd /data/todo/opn052
  cp /home/todo/MEGA_OPEN/mega2026_open6_kti-1024mpi.go .
  cp /home/todo/MEGA_OPEN/opn052/opn052_001.in .

  module load intel/2025.1
  export I_MPI_DEBUG=4
  export I_MPI_OFI_PROVIDER=psm3

  mpirun -np 1024 ./mega2026_open6_kti-1024mpi.go
