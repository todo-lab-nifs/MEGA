# MEGA

MEGA is a simulation code for magnetohydrodynamic (MHD) phenomena in magnetically confined plasmas, focusing on energetic particle driven instabilities, energetic particle transport, and kinetic effects of thermal ions and energetic particles.

---

## Citation

If you use this code, please cite the following work:

Todo, Yasushi et al. *Energetic particle driven Alfvén eigenmodes and associated energetic particle redistribution in a tokamak burning plasma.* Nuclear Fusion 65 (2025) 102003.  
https://doi.org/10.1088/1741-4326/ae059f  

Additionally, please cite the software itself:

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.17170100.svg)](https://doi.org/10.5281/zenodo.17170100)


---

## License

This code is distributed under the MIT License.  
See the [LICENSE](LICENSE) file for details.

---

## Compilation

1. Compile mega2025_open4.f90 using Intel® compiler for Subsystem A of Plasma Simulator (Intel® Xeon® 6980P)
+ % cd src
+ % module load intel/2025.1
+ % make
+ Makefile is for flat MPI parallelization on Subsystem A. The Fortran compiler is the Intel® compiler (mpiifx). The preprocessor 
‘-fpp –DPRM1024MPI –DPTCL1MPI -DXEON6900P’ means that PRM1024MPI, PTCL1MPI, and XEON6900P are designated in the code. PRM1024MPI means 1024 MPI processes for domain decomposition. Particle decomposition is also available with PTCL*MPI, but PTCL1MPI is usually used. If XEON6900P is designated, branches and modules optimized for Xeon® 6980P are used.
+ -fpp -DKTI: for simulations with kinetic thermal ions
+ -fpp -DEXP_EQ: for the equilibrium data constructed using experimental plasma profiles and EFT equilibrium

2. Compile mega2025_open4.f90 using AMD® compiler for Subsystem B of Plasma Simulator (AMD® MI300A)
+ % cd src
+ % module load openmpi/5.0.7/rocm6.3.3_amdflang_afar
+ % make -f Makefile_amdgpu
+ Makefile_amdgpu is for OpenMP and MPI hybrid parallelization on Subsystem B. The Fortran compiler is the AMD® compiler (amdflang). The preprocessor  
‘-cpp –DPRM4MPI–DPTCL1MPI -DAMDGPU -DOPEN_MP -DSMP16’ means that PRM4MPI, PTCL1MPI, AMDGPU, OPEN_MP, and SMP16 are designated in the code. PRM4MPI means 4 MPI processes for domain decomposition. Particle decomposition is also available with PTCL*MPI, but PTCL1MPI is usually used. MI300A is an APU consisting of GPU and CPU with a shared memory. If AMDGPU is designated, branches and modules optimized for MI300A are used. OPEN_MP and SMP16 means that the OpenMP branches are used on the CPU part of MI300A with 16 threads on each MPI process. One MPI process runs on one MI300A. 
+ -cpp -DKTI: for simulations with kinetic thermal ions
+ -cpp -DEXP_EQ: for the equilibrium data constructed using experimental plasma profiles and EFIT equilibrium
+ The built-in mathematical function ‘erfc’ is not available with the AMD® compiler. One may find the function on ‘https://www.netlib.org/specfun/erf’. In addition, a pseudo random number generator based on hip (hipRAND) is used using a wrapper written in C language (hiprandwrapper.cu). 

3. In addition to Plasma Simulator, optimized modules for the supercomputer Fugaku and NEC SX-Aurora TSUBASA, which are designated by the preprocessors with -DFUGAKU and -DAURORA, respectively, are available. The default module is ‘fx100’ which is optimized for Fujitsu Supercomputer PRIMEHPC FX100. 

---

## Usage

1. examples: example jobs for an energetic particle driven Alfvén eigenmode (AE) on Subsystems A  (Intel® Xeon® 6980P) and B (AMD® Instinct MI300A) on Plasma Simulator
+ opn015 without kinetic thermal ions using Subsystem A 
+ opn016 without kinetic thermal ions using Subsystem B
+ opn020 with kinetic thermal ions using Subsystem A
+ opn021 with kinetic thermal ions using Subsystem B

2. equilibrium: equilibrium data for MEGA simulation (241216.eql062.lendian.d) which was used in the simulations presented in Figures 1 and 2 in 
[Y. Todo et al., Plasma Physics and Controlled Fusion 63 (2021) 075018]

3. opn015, opn020
+ Example jobs on Subsystem A for 1024 MPI processes using 4 nodes (=8 CPUs). The job script is go001.sh. The job script should be modified to fit the user’s computer environment. The input parameters are given in “opn015_001.in”. This job focuses on an n=4 AE with an input parameter “PHIMODE=4.0d0” which defines the toroidal angle range 0 &leq; &phi; < 2&pi; / PHIMODE. 
+ The equilibrium data ‘equilibrium/241216.eql062.lendian.d’ is used in this run. 
+ In these runs, the source codes “mega2025_open4.f90” and “optimized_mod25xeon+amd.f90” are used. 

4. opn016, opn021
+ Example jobs on Subsystem B for 4 MPI processes using 1 node (=4 APUs) or for 8 MPI processes using 2 nodes (=8 APUs).

5. diagnostics
+ time evolution of each energy component is output to “opn*_001.energy_phys.txt”.
+ 2D figures in a poloidal plane: opn*_001.ks00100000.pdf 
(input data is opn*_001.moments)
+ radial MHD  velocity profile analysis: profile.f90  
input data: opn*_001.harmonics  
output data: opn*_001_kstep=0100000_n=+04-vrad.txt  
figure are contained in directories opn015 and opn017
+ time evolution of the radial MHD velocity: evolve.f90  
input data: opn*_001.harmonics  
output data: opn*_evol_m=06_n=+04_l=097-vrad.txt  
figures are contained in directories opn015 and opn017

---

## Repository Structure

MEGA/  
├─ src/ # Source code, Makefile, and diagnostics tool  
├─ equilibrium/ # Equilibrium data
├─ examples/ # Example input files and test cases  
├─ README.md # This file  
├─ LICENSE # License information  
├─ CITATION.cff # Citation metadata  
└─ .zenodo.json # Zenodo configuration  


---

## Contact

Author: **Yasushi Todo**  
Affiliation: *National Institute for Fusion Science*  
ORCID: [0000-0001-9323-8285](https://orcid.org/0000-0001-9323-8285)  

For questions or collaboration inquiries, please open an Issue or contact directly.

---

## A Brief Overview of MEGA

1. MEGA = hybrid simulation code: MHD (fluid) + Energetic Particles (PIC)  
[+ Kinetic Thermal Ions (PIC)]

+ MHD: 4th order finite difference + 4th Runge-Kutta for time integration
+ Energetic particles [and kinetic thermal ions]: delta-f PIC  (particle-in-cell)
+ References  
[1]	Y. Todo et al., Plasma Physics and Controlled Fusion 63 (2021) 075018.  
[2]	Y. Todo et al., Nuclear Fusion 65 (2025) 102003. https://doi.org/10.1088/1741-4326/ae059f

2. Computational language and parallelization  
+ Fortran90 + MPI (Message-Passing-Interface) + OpenMP
+ 3-dimensional domain decomposition

3. Coordinates and grid points
+ cylindrical coordinates, grid points are assigned to (R, &phi;, z)
+ magnetic flux coordinates (r, &phi;, &theta;) are also used for data analysis 

4. Important parameters to be specified in the code
+ mpi_proc_r (z, &phi;): number of decomposition in R (z, &phi;) direction
+ lr(z,phi)net: total number of grid points in R (z,phi) direction
+ marker_a0_total: total number of PIC particles for the energetic particles
+ marker_i0_total: total number of PIC particles for the thermal ions
+ minor_r: minor radius (normalized by rho_a=va/omega_a)
+ major_r: major radius (defined in subroutine sgrid)
+ simulation domain:  
major_r - minor_r &leq; R &leq; major_r + minor_r  
0 &leq; z &leq; 2*minor_r  
phimode: 0 &leq; &phi; < 2&pi;/phimode  
+ lpara: lpara=1 for flat MPI parallelization  
lpara=number of threads for MPI + OpenMP parallelization
+ normalization:  
Velocity, time, and length are normalized by Alfven velocity at the 
magnetic axis (va), energetic particle cyclotron frequency (omega_a), 
and gyro radius of energetic particle with Alfven velocity (rho_a=va/omega_a),
respectively.  
In the simulation, mass and charge of energetic particle are unity (m_a=1, e_a=1), 
and the magnetic field and Alfven velocity at the magnetic axis are unity (b0=1, 
va0=1).  
The magnetic permeability is set to be unity. Then, the mass density at 
the magnetic axis is unity. 

5. Input parameters using namelist  
+ job_seq: when starting the simulation, job_seq=1. 
if it is succeeded, job_seq=2 for the next job, and so on.  
+ kstep: initial time step of the job, when starting kstep=0
+ ksmax: final time step of the job  
+ etlim: elapse time limit for the main loop (sec)  
+ kwchk: step interval for subroutine write_check  
+ kwout: step interval for write_data (for job succession)  
+ ksnapsh: step interval for write_snapsh (2D data on a poloidal plane)  
+ kmovie: step interval for write_movie (3D data for movie), not used  
+ dt: time step width (with e_a*b0/m_a=1)  
+ phimode: explained above for the simulation domain  
+ nu0: viscosity normalized by va0*major_r  
+ eta0: resistivity normalized by va0*major_r*mu_0 
        (with magnetic permeability=1)  
+ nu_n0: diffusivity for MHD density  
+ chi0: diffusivity for MHD pressure
+ flag_FLR: FLR of EP and thermal ion, .TRUE. or .FALSE.
+ flag_HM: extended MHD (thermal ion diamag.), .TRUE. or .FALSE
+ flag_BDP: beam injection, .TRUE. or .FALSE.
+ type_a: EP distribution function, 0: Maxwell, 1: isotropic slowing down, 2: slowing down with mu=0, 3: anisotropic slowing down, -5: beam injection
+ clambda0: for type_a=3
+ dclambda: for type_a=3
+ hpower: beam deposition power [MW] for flag_BDP=.TRUE.
+ sd_accl: acceleration of slowing down time, not used
+ tmax0: total simulation time [ms], for flag_BDP=.TRUE.
+ t_classic: classical simulation period [ms], for flag_BDP=.TRUE.
+ t_mhd: hybrid simulation period [ms], for flag_BDP=.TRUE.
+ beta_a0: EP beta at the plasma center
+ scale_psi: scale length of EP spatial profile
+ valpha: maximum velocity for type_a=1, 2, 3, defined in the code
+ temp_a: energetic particle temperature for type_a=0, defined in the code

6. Subroutines called from the main program
+ start_mpi: initialization of MPI processes with domain decomposition 
+ wall_clock: measure elapse time with mpi_wtime
+ equi_solution: read the equilibrium data 
+ sgrid: construct simulation grid points
+ particle_parameters: set up important parameters
+ equilibrium: set up initial MHD equilibrium
+ flux_coordinates: construct magnetic flux coordinates for data analysis
+ beam_deposit: beam deposition profile for flag_BDP=.TRUE.
+ initial_particle: distribute PIC particles
+ initial_mhd_balance: neglect numerical error in the equilibrium data
+ injection: beam injection for flag_BDP=.TRUE.
+ scattering: pitch-angle scattering and energy diffusion of EPs for 
flag_BDP=.TRUE
+ perturb_fluid(_random): give initial perturbation to the MHD fluid
+ read_data: read the data of the previous job for succession
+ density_particle: calculate energetic-particle pressures from the PIC particles
+ e_field: give electric field from the Ohm's law
+ wirte_check or dissipation_analysis: output step number, energy evolution
+ write_snapsh: output 2D data on a poloidal plane
+ write_movie: output 3D data for movie (currently commented out)
+ harmonics: decompose data to Fourier modes on magnetic flux coordinates and output
+ write_data: write all the data necessary for job succession
+ order_particle: order particles for fast computation
+ satellite: finite Larmor radius effects for flag_FLR=.TRUE.
+ push_particle: drift kinetic equation of motion of energetic particles
+ bcptcl: boundary condition of energetic particles
+ mhd: MHD equations and time integration
+ com_particle: PIC particle data is transferred to the next MPI domain
+ switch_classic: for type_a=-5 with beam injection
+ mhd_lowpass: filter of toroidal mode number within |n|<=lphi_n for MHD
+ mhd_smoothing: numerical filter of the MHD data variation in each ksmth steps
+ end_mpi: finalization of the MPI processes




