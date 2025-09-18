#!/bin/bash

#PBS -N opn016
#PBS -P project_name
#PBS -q B_S
#PBS -l select=4:ncpus=24:ngpus=1:mem=124gb
#PBS -l walltime=10:00:00
export OMP_STACKSIZE=4G
export OMP_NUM_THREADS=16
export HSA_XNACK=1
export UCX_LOG_LEVEL=error

# export F_UFMTENDIAN=60
# export FORT1='./opn016_001.in'
# export FORT7='/home/user_name/MEGA_OPEN/opn016/opn016_001.out'
# export FORT30='./opn016_001.moments'
# export FORT31='./opn016_001.harmonics'
# export FORT40='/home/user_name/MEGA_OPEN/opn016/opn016_001.energy_phys.txt'
# export FORT41='/home/user_name/MEGA_OPEN/opn016/opn016_001.energy_n.txt'
# export FORT60='/data/user_name/241216.eql062.lendian.d'

  mkdir -p /data/user_name/opn016
  cd /data/user_name/opn016
  cp /home/user_name/MEGA_OPEN/mega2025_open3amd-4mpi.go .
  cp /home/user_name/MEGA_OPEN/opn016/opn016_001.in ./fort.1
  cp /data/user_name/241216.eql062.lendian.d ./fort.60
  cp /home/user_name/MEGA_OPEN/opn016/kick.sh .
#  cp /home/user_name/MEGA_OPEN/opn016/myrankfile .

  module load openmpi/5.0.7/rocm6.3.3_amdflang_afar

  mpirun -np 4 --cpus-per-rank 24 ./kick.sh 1> stdout 2> stderr

  mv ./fort.7 /home/user_name/MEGA_OPEN/opn016/opn016_001.out
  mv ./fort.30 ./opn016_001.moments
  mv ./fort.31 ./opn016_001.harmonics
  mv ./fort.40 /home/user_name/MEGA_OPEN/opn016/opn016_001.energy_phys.txt
  mv ./fort.41 /home/user_name/MEGA_OPEN/opn016/opn016_001.energy_n.txt
