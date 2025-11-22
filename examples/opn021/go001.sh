#!/bin/bash

#PBS -N opn021
#PBS -P NIFS24KISC014
#PBS -q B_S
#PBS -l select=8:ncpus=24:ngpus=1:mem=124gb
#PBS -l walltime=24:00:00
export OMP_STACKSIZE=4G
export OMP_NUM_THREADS=16
export HSA_XNACK=1
export UCX_LOG_LEVEL=error

# export F_UFMTENDIAN=60
# export FORT1='./opn021_001.in'
# export FORT7='/home/todo/MEGA_OPEN/opn021/opn021_001.out'
# export FORT30='./opn021_001.moments'
# export FORT31='./opn021_001.harmonics'
# export FORT40='/home/todo/MEGA_OPEN/opn021/opn021_001.energy_phys.txt'
# export FORT41='/home/todo/MEGA_OPEN/opn021/opn021_001.energy_n.txt'
# export FORT60='/data/todo/241216.eql062.lendian.d'

  mkdir -p /data/todo/opn021
  cd /data/todo/opn021
  cp /home/todo/MEGA_OPEN/mega2025_open4amd_kti-8mpi.go .
  cp /home/todo/MEGA_OPEN/opn021/opn021_001.in ./fort.1
  cp /data/todo/241216.eql062.lendian.d ./fort.60
  cp /home/todo/MEGA_OPEN/opn021/kick.sh .
#  cp /home/todo/MEGA_OPEN/opn021/myrankfile .

  module load openmpi/5.0.7/rocm6.3.3_amdflang_afar

  mpirun -np 8 --cpus-per-rank 24 ./kick.sh 1> stdout 2> stderr

  mv ./fort.7 /home/todo/MEGA_OPEN/opn021/opn021_001.out
  mv ./fort.30 ./opn021_001.moments
  mv ./fort.31 ./opn021_001.harmonics
  mv ./fort.40 /home/todo/MEGA_OPEN/opn021/opn021_001.energy_phys.txt
  mv ./fort.41 /home/todo/MEGA_OPEN/opn021/opn021_001.energy_n.txt
