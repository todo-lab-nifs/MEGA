#!/bin/bash

#PBS -N opn015
#PBS -P project_name
#PBS -q A_S
#PBS -l select=4:ncpus=256:mpiprocs=256:ompthreads=1
#PBS -l walltime=05:00:00

# export F_UFMTENDIAN=60
export FORT1='./opn015_001.in'
export FORT7='/home/user_name/MEGA_OPEN/opn015/opn015_001.out'
export FORT30='./opn015_001.moments'
export FORT31='./opn015_001.harmonics'
export FORT40='/home/user_name/MEGA_OPEN/opn015/opn015_001.energy_phys.txt'
export FORT41='/home/user_name/MEGA_OPEN/opn015/opn015_001.energy_n.txt'
export FORT60='/data/user_name/241216.eql062.lendian.d'

  mkdir -p /data/user_name/opn015
  cd /data/user_name/opn015
  cp /home/user_name/MEGA_OPEN/mega2025_open3-1024mpi.go .
  cp /home/user_name/MEGA_OPEN/opn015/opn015_001.in .

  module load intel/2025.1
  export I_MPI_DEBUG=4
  export I_MPI_OFI_PROVIDER=psm3

  mpirun -np 1024 ./mega2025_open3-1024mpi.go
