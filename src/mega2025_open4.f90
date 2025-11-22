module parameters
#ifdef PRM65536MPI
  integer, parameter::mpi_proc_r=64,mpi_proc_z=64,mpi_proc_phi=16
#elif PRM32768MPI
  integer, parameter::mpi_proc_r=32,mpi_proc_z=64,mpi_proc_phi=16
#elif PRM16384MPI
  integer, parameter::mpi_proc_r=32,mpi_proc_z=64,mpi_proc_phi=8
#elif PRM8192MPI
  integer, parameter::mpi_proc_r=32,mpi_proc_z=32,mpi_proc_phi=8
#elif PRM4096MPI
  integer, parameter::mpi_proc_r=16,mpi_proc_z=32,mpi_proc_phi=8
#elif PRM2048MPI
  integer, parameter::mpi_proc_r=16,mpi_proc_z=16,mpi_proc_phi=8
#elif PRM1024MPI
  integer, parameter::mpi_proc_r=16,mpi_proc_z=16,mpi_proc_phi=4
#elif PRM512MPI
  integer, parameter::mpi_proc_r=8,mpi_proc_z=16,mpi_proc_phi=4
#elif PRM256MPI
  integer, parameter::mpi_proc_r=8,mpi_proc_z=16,mpi_proc_phi=2
#elif PRM128MPI
  integer, parameter::mpi_proc_r=8,mpi_proc_z=8,mpi_proc_phi=2
#elif PRM64MPI
  integer, parameter::mpi_proc_r=4,mpi_proc_z=8,mpi_proc_phi=2
#elif PRM32MPI
  integer, parameter::mpi_proc_r=2,mpi_proc_z=8,mpi_proc_phi=2
#elif PRM16MPI
  integer, parameter::mpi_proc_r=2,mpi_proc_z=8,mpi_proc_phi=1
#elif PRM8MPI
  integer, parameter::mpi_proc_r=2,mpi_proc_z=4,mpi_proc_phi=1
#elif PRM4MPI
! if mpi_proc_r=1, lr=lrnet,  if mpi_proc_z=1, lz=lznet
!  integer, parameter::mpi_proc_r=1,mpi_proc_z=1,mpi_proc_phi=4
  integer, parameter::mpi_proc_r=2,mpi_proc_z=2,mpi_proc_phi=1
#elif PRM1MPI
  integer, parameter::mpi_proc_r=1,mpi_proc_z=1,mpi_proc_phi=1
#else
  integer, parameter::mpi_proc_r=8,mpi_proc_z=8,mpi_proc_phi=4
#endif

!---------- modified for particle parallelization -------------------
!                     date:   2012-06-17
!--------------------------------------------------------------------
#ifdef PTCL16MPI
  integer, parameter::mpi_proc_ptcl=16
#elif PTCL8MPI
  integer, parameter::mpi_proc_ptcl=8
#elif PTCL4MPI
  integer, parameter::mpi_proc_ptcl=4
#elif PTCL2MPI
  integer, parameter::mpi_proc_ptcl=2
#elif PTCL1MPI
  integer, parameter::mpi_proc_ptcl=1
#else
  integer, parameter::mpi_proc_ptcl=1
#endif

  integer, parameter::mpi_proc_mhd=mpi_proc_r*mpi_proc_z*mpi_proc_phi
  integer, parameter::mpi_proc=mpi_proc_mhd*mpi_proc_ptcl
  integer, parameter::mpi_proc_pol=mpi_proc_r*mpi_proc_z

  integer, parameter::mpi_proc_mhd1=mpi_proc_mhd-1
  integer, parameter::mpi_proc1=mpi_proc-1
  integer, parameter::mpi_proc_pol1=mpi_proc_pol-1

!---- end: modified for particle parallelization 2012-06-17 ----------

  integer, parameter::lrnet=128,lznet=128,lphinet=32,lphinet4=lphinet+4
!  integer, parameter::lr_shd=2,lz_shd=2,lphi_shd=2 ! w/o FLR, 2017-09-16
  integer, parameter::lr_shd=6,lz_shd=6,lphi_shd=2 ! shadow width
#ifdef PRM1MPI
  integer, parameter::lr=lrnet
  integer, parameter::lz=lznet
!#elif PRM4MPI !if mpi_proc_r=1 and mpi_proc_z=1
!  integer, parameter::lr=lrnet
!  integer, parameter::lz=lznet
!  integer, parameter::lz=lznet/mpi_proc_z + 2*lz_shd
! #elif PRM4MPI !if mpi_proc_z=1
!  integer, parameter::lr=lrnet/mpi_proc_r + 2*lr_shd
!  integer, parameter::lz=lznet
#else
  integer, parameter::lr=lrnet/mpi_proc_r + 2*lr_shd
  integer, parameter::lz=lznet/mpi_proc_z + 2*lz_shd
#endif
  integer, parameter::lphi=lphinet/mpi_proc_phi + 2*lphi_shd
  integer, parameter::lrz=lr*lz,lrphi=lr*lphi,lzphi=lz*lphi,lrzphi=lr*lz*lphi
  integer, parameter::lr_loc=lr-2*lr_shd,lz_loc=lz-2*lz_shd,lphi_loc=lphi-2*lphi_shd !2016-01-08
  integer, parameter::lrzphi_loc=lr_loc*lz_loc*lphi_loc !2016-01-08
  integer, parameter::nr=lrnet,nz=lznet/2,nz2=lznet
  integer, parameter::lpsi=201,ltheta=256,lpsi2=51,lcfpphi=9 !2019-06-17

#ifdef SMP32
  integer, parameter::lpara=32
#elif SMP20
  integer, parameter::lpara=20
#elif SMP18
  integer, parameter::lpara=18
#elif SMP16
  integer, parameter::lpara=16
#elif SMP12
  integer, parameter::lpara=12
#elif SMP10
  integer, parameter::lpara=10
#elif SMP8
  integer, parameter::lpara=8
#elif SMP4
  integer, parameter::lpara=4
#elif SMP2
  integer, parameter::lpara=2
#else
  integer, parameter::lpara=1
#endif

  integer, parameter::marker_a0=(lrnet/mpi_proc_r)*(lznet/mpi_proc_z) &
                               *(lphinet/mpi_proc_phi) * 8
!                               *(lphinet/mpi_proc_phi) * 64
!  integer, parameter::marker_i0=marker_a0 !2021-04-29
#ifdef KTI
  integer, parameter::marker_i0=marker_a0
#else
  integer, parameter::marker_i0=2**0
#endif
  integer, parameter::marker_e0=2**0

  integer, parameter::marker_e=marker_e0
  integer, parameter::marker_i=marker_i0+marker_i0 !2024-03-30
  integer, parameter::marker_a=marker_a0+marker_a0 !2024-03-30
  integer, parameter::ncomm=max(marker_e0, marker_i0, marker_a0) !2024-03-30
  integer, parameter::marker_each=max(marker_e, marker_i, marker_a)
  integer, parameter::marker_local_max=128

!size for gc_e, gc_i, gc_a, 2015-06-23
! gc_x(:,n) = 1:r, 2:z, 3:phi, 4:p, 5:mu, 6:weight, 7:active, 8:p_phi, 9:fff, 10:fnrml
!  integer, parameter::ngc1=8,ngc2=10
!  integer, parameter::ngc1=9,ngc2=11
  integer, parameter::ngc1=11,ngc2=11 !dfnrml is calcualted

!size for flp and flp_gyro
  integer, parameter::nflp=33,nflp_gyro=7 !2025-04-04, add grad B0
!  integer, parameter::nflp=30,nflp_gyro=7 !2016-08-05
 
! _a: deuterium beam, _i: deuterium, _e: electron not used
! _ti: thermal ion
! normalization is for deuterium, then m_a=1, e_a=1
  real(8), parameter::m_e=0.50d0/1.836d3,e_e=-1.0d0
  real(8), parameter::m_ti=1.0d0,e_ti=1.0d0
  real(8), parameter::m_i=1.0d0,e_i=1.0d0
  real(8), parameter::m_a=1.0d0,e_a=1.0d0
  real(8), parameter::pi=3.14159265358979323846d0,twopi=2.0d0*pi
!  real(8), parameter::pi=3.14159265358979323846d0,twopi=2.0d0*pi,sqr_pi1&
!       &=1.0d0/sqrt(pi)
!  real(8), parameter::minor_r=2.0d1,phimode=1.0d0

! !Fredricson(2009), 2025-02-05
!  real(8), parameter::minor_r=1.22d1
!  real(8), parameter::aspect=1.344d0, dkappa=2.049d0

!Todo(2021), 2025-04-04
!  real(8), parameter::minor_r=1.6d1
  real(8), parameter::aspect=3.2d0, dkappa=1.0d0
  
  real(8), parameter::m_u=1.66053906660d-27
  real(8), parameter::m_H=1.00782503207d0*m_u
  real(8), parameter::m_D=2.0141017778*m_u
  real(8), parameter::m_T=3.0160492777d0*m_u
  real(8), parameter::m_He=4.00260325415d0*m_u
  real(8), parameter::m_elc=5.48579909065d-4*m_u
  real(8), parameter::e_u=1.602176634d-19
  real(8), parameter::munit=m_D,eunit=e_u !2025-02-05
!  real(8), parameter::munit=m_He,eunit=2.0d0*e_u !2021-07-15

  real(8), parameter::b0=1.0d0,eps_b=5.0d-1
  real(8), parameter::va0=1.0d0 !2012-02-06 Alfven vlc. normalized to unity on axis
  real(8)::va,omega_a,rho_a !2012-02-06 input parameters in equi_solution
  real(8)::bunit !2013-05-03 [T]
end module parameters

module job_cntrol
  integer::job_seq
end module

module mpiset
  use mpi
  use parameters
!!!HITACHI
  integer::new_d7i2_type
!!!HITACHI
  integer, parameter::lrequest=10  !2012-06-29
  integer::my_rank,nprocess,my_rank_r,my_rank_z,my_rank_phi,my_rank_pol !2015-07-16
  integer::my_rank_mhd,my_rank_ptcl !2012-06-17
  integer::lrstart,lrend,lzstart,lzend,lphistart,lphiend
  integer::my_status(mpi_status_size),mpi_err
!2012-06-29
  integer::my_request(lrequest),my_status_all(mpi_status_size,lrequest)
!2012-06-29 end
  integer::kr_offset(0:mpi_proc1),kz_offset(0:mpi_proc1)&
       &,kphi_offset(0:mpi_proc1)
  integer::kr_s(0:mpi_proc1),kr_e(0:mpi_proc1)
  integer::kz_s(0:mpi_proc1),kz_e(0:mpi_proc1)
  integer::isize(3),isubsize(3),isend_up(3),irecv_up(3) &
       &,isend_down(3),irecv_down(3)
  integer::ir_sendup,ir_senddown,ir_recvup,ir_recvdown
  integer::jz_sendup,jz_senddown,jz_recvup,jz_recvdown

! modified on 2011-01-29
  integer::r_world,z_world,phi_world,poloidal_world

!2012-06-17
  integer::ptcl_allreduce
  integer::mhd_world,ptcl_world
!2012-06-17 end
!2012-06-29
!  integer::mpi_satellite,mpi_ptcl
end module

module equi_sol
  use parameters
#ifdef EXPEQ
  real(8)::psisim(nr,nz2),brsim(nr,nz2),bzsim(nr,nz2),bhsim(nr,nz2)
  real(8)::ppsim(nr,nz2),pfasim(nr,nz2),pnbsim(nr,nz2)
  real(8)::rhosim(nr,nz2),dns_asim(nr,nz2),dns_isim(nr,nz2),vcritsim(nr,nz2)
  real(8)::nusdsim(nr,nz2),tisim(nr,nz2),tesim(nr,nz2),rotsim(nr,nz2) !2013-12-18
  real(8)::raxis,zaxis
  real(8)::gpsi(lpsi),rpsi(lpsi),qpsi(lpsi)
  real(8)::gpsi_nrm(lpsi)
  real(8)::pp_psi(lpsi),pfa_psi(lpsi),pnb_psi(lpsi),rho_psi(lpsi)
#else
  real(8)::psi1(nr,nz),curden(nr,nz),presur(nr,nz),di2(nr,nz)
  real(8)::presurh_para(nr,nz),presurh_perp(nr,nz),presurh_perp2(nr,nz)
  real(8)::gsbr(nr,nz),gsbphi(nr,nz),gsbz(nr,nz),gsrho(nr,nz),gste(nr,nz)
  real(8)::gsparacurr(nr,nz),gsparacurphi(nr,nz),gsparacurz(nr,nz)
  real(8)::gsbr2(nr,nz2),gsbphi2(nr,nz2),gsbz2(nr,nz2),gsrho2(nr,nz2),gste2(nr,nz2)
  real(8)::prsr2(nr,nz2),prsr2h_para(nr,nz2),prsr2h_perp(nr,nz2),prsr2h_perp2(nr,nz2)
  real(8)::psi2(nr,nz2),gscurr2(nr,nz2),gscurphi2(nr,nz2),gscurz2(nr,nz2)
  integer::iaxis,jaxis
  real(8)::draxis,brat,raxis
  real(8)::gpsi(lpsi),rpsi(lpsi),qpsi_exp(lpsi),qpsi(lpsi)
  real(8)::p_b_psi(lpsi),p_h_psi(lpsi),rho_psi(lpsi)
  real(8)::gpsi_nrm(lpsi)
  real(8)::ddi2psi(lpsi2)
  integer::iter_eq
#endif
end module equi_sol

module flux_c
  use parameters
  integer, parameter::n_flux_size=lpsi*ltheta/max(1,min(mpi_proc_r, mpi_proc_z)/2)
  integer, parameter::mpol=64,ntor=2 !ntor <= lphi_n, !2025-09-03
  integer, parameter::ntor_phi=ntor/mpi_proc_phi+1,rank_phi_max=ntor/ntor_phi

  real(8)::crdmagr(lpsi,ltheta),crdmagz(lpsi,ltheta),crdmagphi(lpsi,ltheta)
  real(8)::dtheta
  integer::n_flux_grid
  real(8)::r_flux_grid(n_flux_size),z_flux_grid(n_flux_size)
  integer::i_flux_grid(n_flux_size),l_flux_grid(n_flux_size)
  real(8)::dt_r(n_flux_size),dt_z(n_flux_size),dt_phi(n_flux_size)&
       &,dt_abs(n_flux_size)
  real(8)::dr_r(n_flux_size),dr_z(n_flux_size),dr_phi(n_flux_size)&
       &,dr_abs(n_flux_size)
  real(8)::dp_r(n_flux_size),dp_z(n_flux_size),dp_phi(n_flux_size)
!2012-07-11
!  real(8)::cos_flux_grid(0:mpol,-ntor:ntor,ltheta,lphi)
!  real(8)::sin_flux_grid(0:mpol,-ntor:ntor,ltheta,lphi)
  real(8)::cos_theta_flux(0:mpol,ltheta),sin_theta_flux(0:mpol,ltheta)
  real(8)::cos_phi_flux(-ntor:ntor,lphi),sin_phi_flux(-ntor:ntor,lphi)
!2012-07-11 end

end module

module field
   use parameters
! modified on 2011-01-29
   real(8)::fld(nflp,lr,lz,lphi)
   real(8)::epara(lr,lz,lphi),epara0(lr,lz,lphi)
   real(8)::psi(lr,lz,lphi),qvalue(lr,lz,lphi)
   real(8)::psi_r(lr,lz,lphi),psi_z(lr,lz,lphi),psi_phi(lr,lz,lphi)
   real(8)::psimax !2012-02-06
   real(8)::psivac !2012-08-24
   real(8)::sqrtpsi_chi_enh,dsqrtpsi_chi_enh !2013-02-17
   real(8)::sqrtpsi_pe,dsqrtpsi_pe !2021-07-25
   real(8)::psi_edge !2023-03-13
   real(8)::rho(lr,lz,lphi),prs(lr,lz,lphi),enrg(lr,lz,lphi) !2021-05-05
   real(8)::prs_r(lr,lz,lphi),prs_z(lr,lz,lphi),prs_phi(lr,lz,lphi) !2017-04-04
   real(8)::vr(lr,lz,lphi),vz(lr,lz,lphi),vphi(lr,lz,lphi)
   real(8)::vdrftr(lr,lz,lphi),vdrftz(lr,lz,lphi),vdrftphi(lr,lz,lphi) !2012-08-16
   real(8)::ur(lr,lz,lphi),uz(lr,lz,lphi),uphi(lr,lz,lphi),div_u(lr,lz,lphi) !2012-04-09
   real(8)::omegar(lr,lz,lphi),omegaz(lr,lz,lphi),omegaphi(lr,lz,lphi) !2012-04-16
   real(8)::er(lr,lz,lphi),ez(lr,lz,lphi),ephi(lr,lz,lphi)
   real(8)::er0(lr,lz,lphi),ez0(lr,lz,lphi),ephi0(lr,lz,lphi)

!2019-06-04s
   real(8)::er_res(lr,lz,lphi),ez_res(lr,lz,lphi),ephi_res(lr,lz,lphi)
   real(8)::er0_res(lr,lz,lphi),ez0_res(lr,lz,lphi),ephi0_res(lr,lz,lphi)
!2019-06-04e

   real(8)::br(lr,lz,lphi),bz(lr,lz,lphi),bphi(lr,lz,lphi),babs(lr,lz,lphi)
   real(8)::br0(lr,lz,lphi),bz0(lr,lz,lphi),bphi0(lr,lz,lphi),babs0(lr,lz,lphi)
   real(8)::rho0(lr,lz,lphi),prs0(lr,lz,lphi),enrg0(lr,lz,lphi) !2021-05-05
   real(8)::vtotphi(lr,lz,lphi),utotphi(lr,lz,lphi) !2013-12-18
   real(8)::vpir(lr,lz,lphi),vpiz(lr,lz,lphi),vpiphi(lr,lz,lphi) !2014-04-29

   real(8)::gradbr(lr,lz,lphi),gradbz(lr,lz,lphi),gradbphi(lr,lz,lphi)
   real(8)::curvbr(lr,lz,lphi),curvbz(lr,lz,lphi),curvbphi(lr,lz,lphi)
   real(8)::cr(lr,lz,lphi),cz(lr,lz,lphi),cphi(lr,lz,lphi)
   real(8)::cr0(lr,lz,lphi),cz0(lr,lz,lphi),cphi0(lr,lz,lphi)
!2012-01-04   real(8)::charge_density(lr,lz,lphi),current_density(lr,lz,lphi)
!2012-01-04   real(8)::c_dot_b(lr,lz,lphi),e_dot_b(lr,lz,lphi)
   real(8)::drho(lr,lz,lphi),rho1(lr,lz,lphi),rho2(lr,lz,lphi),drho0(lr,lz,lphi),rho_ref(lr,lz,lphi)
   real(8)::dprs(lr,lz,lphi),prs1(lr,lz,lphi),prs2(lr,lz,lphi),dprs0(lr,lz,lphi),prs_ref(lr,lz,lphi)
   real(8)::denrg(lr,lz,lphi),enrg1(lr,lz,lphi),enrg2(lr,lz,lphi),denrg0(lr,lz,lphi),enrg_ref(lr,lz,lphi) !2021-05-05
   real(8)::dvr(lr,lz,lphi),vr1(lr,lz,lphi),vr2(lr,lz,lphi),dvr0(lr,lz,lphi),vr_ref(lr,lz,lphi)
   real(8)::dvz(lr,lz,lphi),vz1(lr,lz,lphi),vz2(lr,lz,lphi),dvz0(lr,lz,lphi),vz_ref(lr,lz,lphi)
   real(8)::dvphi(lr,lz,lphi),vphi1(lr,lz,lphi),vphi2(lr,lz,lphi) &
           ,dvphi0(lr,lz,lphi),vphi_ref(lr,lz,lphi)
   real(8)::dbr(lr,lz,lphi),br1(lr,lz,lphi),br2(lr,lz,lphi),dbr0(lr,lz,lphi),br_ref(lr,lz,lphi)
   real(8)::dbz(lr,lz,lphi),bz1(lr,lz,lphi),bz2(lr,lz,lphi),dbz0(lr,lz,lphi),bz_ref(lr,lz,lphi)
   real(8)::dbphi(lr,lz,lphi),bphi1(lr,lz,lphi),bphi2(lr,lz,lphi) &
           ,dbphi0(lr,lz,lphi),bphi_ref(lr,lz,lphi)
   real(8)::nu0,eta0,nu_n0,chi0,nu,nu_vac,eta,nu_n,chi,gamma,chi_enhanced,nun_enhanced !2013-02-17
   real(8)::vac(lr,lz,lphi),nu_spt(lr,lz,lphi),eta_spt(lr,lz,lphi),chi_spt(lr,lz,lphi) !2013-02-17
   real(8)::vac_pe(lr,lz,lphi) !2021-07-25
   real(8)::nun_spt(lr,lz,lphi) !2013-02-27
   real(8)::bmin

!2016-08-05s
   real(8)::dns_e(lr,lz,lphi),mom_e(lr,lz,lphi) &
           ,ppara_e(lr,lz,lphi),pperp_e(lr,lz,lphi) &
           ,dns_e0(lr,lz,lphi),mom_e0(lr,lz,lphi),temp_e0(lr,lz,lphi) &
           ,ppara_e0(lr,lz,lphi),pperp_e0(lr,lz,lphi) &
           ,dns_e0_r(lr,lz,lphi),dns_e0_z(lr,lz,lphi),dns_e0_phi(lr,lz,lphi) &
           ,temp_e0_r(lr,lz,lphi),temp_e0_z(lr,lz,lphi),temp_e0_phi(lr,lz,lphi)
!           ,ppara_e_r(lr,lz,lphi),ppara_e_z(lr,lz,lphi),ppara_e_phi(lr,lz,lphi) &
!           ,ppara_e0_r(lr,lz,lphi),ppara_e0_z(lr,lz,lphi),ppara_e0_phi(lr,lz,lphi)
   real(8)::dns_i(lr,lz,lphi),mom_i(lr,lz,lphi) &
           ,ppara_i(lr,lz,lphi),pperp_i(lr,lz,lphi) &
           ,dns_i0(lr,lz,lphi),mom_i0(lr,lz,lphi),temp_i0(lr,lz,lphi) &
           ,ppara_i0(lr,lz,lphi),pperp_i0(lr,lz,lphi) &
           ,dns_i0_r(lr,lz,lphi),dns_i0_z(lr,lz,lphi),dns_i0_phi(lr,lz,lphi) &
           ,temp_i0_r(lr,lz,lphi),temp_i0_z(lr,lz,lphi),temp_i0_phi(lr,lz,lphi)
!           ,ppara_i_r(lr,lz,lphi),ppara_i_z(lr,lz,lphi),ppara_i_phi(lr,lz,lphi) &
!           ,ppara_i0_r(lr,lz,lphi),ppara_i0_z(lr,lz,lphi),ppara_i0_phi(lr,lz,lphi)
   real(8)::dns_a(lr,lz,lphi),mom_a(lr,lz,lphi) &
           ,ppara_a(lr,lz,lphi),pperp_a(lr,lz,lphi) &
           ,dns_a0(lr,lz,lphi),mom_a0(lr,lz,lphi),temp_a0(lr,lz,lphi) &
           ,ppara_a0(lr,lz,lphi),pperp_a0(lr,lz,lphi) &
           ,dns_a0_r(lr,lz,lphi),dns_a0_z(lr,lz,lphi),dns_a0_phi(lr,lz,lphi) &
           ,temp_a0_r(lr,lz,lphi),temp_a0_z(lr,lz,lphi),temp_a0_phi(lr,lz,lphi)
!           ,ppara_a_r(lr,lz,lphi),ppara_a_z(lr,lz,lphi),ppara_a_phi(lr,lz,lphi) &
!           ,ppara_a0_r(lr,lz,lphi),ppara_a0_z(lr,lz,lphi),ppara_a0_phi(lr,lz,lphi)

   real(8)::dns_e_min=1.0d-2,temp_e_min=1.0d-4
   real(8)::dns_i_min=1.0d-2,temp_i_min=1.0d-4
   real(8)::dns_a_min=1.0d-5,temp_a_min=1.0d-4

   real(8)::qpara_a(lr,lz,lphi),qperp_a(lr,lz,lphi) &
           ,qpara_a0(lr,lz,lphi),qperp_a0(lr,lz,lphi)
   real(8)::qpara_i(lr,lz,lphi),qperp_i(lr,lz,lphi) &
           ,qpara_i0(lr,lz,lphi),qperp_i0(lr,lz,lphi)

   real(8)::vcrit(lr,lz,lphi)
   real(8)::nusd(lr,lz,lphi),ti(lr,lz,lphi),te(lr,lz,lphi) !2013-05-02
   real(8)::zeff !2013-05-02
   real(8)::utor(lr,lz,lphi) !2013-12-18

   real(4)::movie_work(lr,lz,lphi,0:mpi_proc_mhd1) !2013-07-17
end module

module grid
   use parameters
#ifdef EXPEQ
   real(8)::major_r,minor_r,phimode
#else
   real(8)::major_r,phimode
   real(8)::minor_r=1.6d1
#endif
   real(8)::rleng,zleng,phileng,dr,dz,dphi
   real(8)::grr(lr,lz,lphi),gzz(lr,lz,lphi),gphi(lr,lz,lphi)
   real(8)::gtht(lr,lz,lphi),gaa(lr,lz,lphi)
   real(8)::dt,t
   integer::kstep,ksmax,kwout,kwchk,ksnapsh,ksmth,kmovie,kharmonics
   integer::kperturb,kscatter !2013-05-02
   logical::flag_classic,flag_stored !2013-05-02
end module

module particle
   use parameters
! setting parameters
   integer::type_e,type_d,type_t,type_i,type_a
   real(8)::temp_e,temp_i,temp_a,beta0,beta_a0,n_a0
   real(8)::temp_d,temp_t !2021-06-12
   real(8)::vmin_e,vmax_e,vmin_i,vmax_i,vmin_a,vmax_a
   real(8)::vmin_d,vmax_d,vmin_t,vmax_t !2021-06-12
   real(8)::valpha,vbeam,deltav  ! 2025-02-05
   real(8)::vbfactor !2024-03-19, control vbeam 
   real(8)::clambda0,dclambda !2012-06-17
   real(8)::vcrit0 !2025-02-05
   real(8)::scale_psi !2016-10-15
!2013-10-18 for beam injection
   real(8)::hpower,sd_accl,sd_time,nu_sd,tmax,tmax0,pas0,t_classic,t_mhd
   real(8)::nu_krook !2024-04-23
!      real(4), allocatable::r(:),z(:),lambda(:),e(:),zeta(:),v(:)
   real(4), allocatable::bdp_data(:,:),t_bd(:)
   integer, allocatable::injection_bd(:)
   integer::n_bd,n_inj
   real(8)::co_stored,cntr_stored,co_stored_total,cntr_stored_total
   real(8)::rmin,rmax,zmin,zmax,phimin,phimax

! electrons
   real(8)::gc_e(ngc2, marker_e)
   real(8)::dgc_e(ngc1,marker_e),gc1_e(ngc1,marker_e),gc2_e(ngc1,marker_e)
   real(8)::v_e(marker_e),lambda_e(marker_e) &
           ,dwpsi_e(marker_e),dwenrc_e(marker_e) &
           ,cf_pphi_e(0:lcfpphi),pphi_min_e,pphi_max_e,total_energy_e
   real(8)::flp_e(nflp,marker_e)
   integer::node_e(marker_e),node_now_e(marker_e)

! ions
   real(8)::gc_i(ngc2, marker_i)
   real(8)::dgc_i(ngc1,marker_i),gc1_i(ngc1,marker_i),gc2_i(ngc1,marker_i)
   real(8)::v_i(marker_i),lambda_i(marker_i) &
           ,dwpsi_i(marker_i),dwenrc_i(marker_i) &
           ,cf_pphi_i(0:lcfpphi),pphi_min_i,pphi_max_i,total_energy_i
   real(8)::flp_i(nflp,marker_i)
   integer::node_i(marker_i),node_now_i(marker_i)

! alpha particles
   real(8)::gc_a(ngc2, marker_a)
   real(8)::dgc_a(ngc1,marker_a),gc1_a(ngc1,marker_a),gc2_a(ngc1,marker_a)
   real(8)::v_a(marker_a),lambda_a(marker_a) &
           ,dwpsi_a(marker_a),dwenrc_a(marker_a) &
           ,cf_pphi_a(0:lcfpphi),pphi_min_a,pphi_max_a,total_energy_a
   real(8)::flp_a(nflp,marker_a)
   integer::node_a(marker_a),node_now_a(marker_a)
end module

! for FLR effects of alpha particles and D beam (ion)
module gyro
   use parameters
   integer, parameter::ngyro=4
   integer, parameter::marker_i_gyro=marker_i*ngyro
   integer, parameter::marker_a_gyro=marker_a*ngyro
   integer, parameter::marker_each_gyro=max(marker_i_gyro,marker_a_gyro)
   integer, parameter::ncomm_gyro=marker_each_gyro/ngyro
! ions
!   real(8)::r_i_gyro(marker_i_gyro),z_i_gyro(marker_i_gyro),phi_i_gyro(marker_i_gyro) &
!           ,p_i_gyro(marker_i_gyro),mu_i_gyro(marker_i_gyro) &
!           ,weight_i_gyro(marker_i_gyro),active_i_gyro(marker_i_gyro)
   real(8)::gyro_i(2,marker_i_gyro)
   real(8)::flp_i_gyro(nflp_gyro,marker_i_gyro)
!   integer::node_i_gyro(marker_i_gyro),node_now_i_gyro(marker_i_gyro)
!   integer::no_i_gyro(2,marker_i_gyro)

! alpha particles
!   real(8)::r_a_gyro(marker_a_gyro),z_a_gyro(marker_a_gyro),phi_a_gyro(marker_a_gyro) &
!           ,p_a_gyro(marker_a_gyro),mu_a_gyro(marker_a_gyro) &
!           ,weight_a_gyro(marker_a_gyro),active_a_gyro(marker_a_gyro)
   real(8)::gyro_a(2,marker_a_gyro)
   real(8)::flp_a_gyro(nflp_gyro,marker_a_gyro)
!   integer::node_a_gyro(marker_a_gyro),node_now_a_gyro(marker_a_gyro)
!   integer::no_a_gyro(2,marker_a_gyro)


!!!HITACHI
   type d7i2
     real(8) :: d(7)
     integer :: i(2)
   end type
   type(d7i2) ::tsend_up(ncomm_gyro),tsend_down(ncomm_gyro) &
               ,trecv_up(ncomm_gyro),trecv_down(ncomm_gyro)
!!!HITACHI

!2012-06-29
!   integer::nsend_up(ncomm_gyro),nsend_down(ncomm_gyro)
!   type sr_buf
!     integer::j_buf(3)
!     real(8)::r_buf(7)
!   end type sr_buf
!   type(sr_buf)::send_up(ncomm_gyro),send_down(ncomm_gyro) &
!                ,recv_up(ncomm_gyro),recv_down(ncomm_gyro)
!2012-06-29 end
end module

#ifdef MKL
   include 'mkl_vsl.f90'
#endif

module random_num !2012-01-28 with random_number of fortran95
   use parameters
#ifdef HITACHI
   integer(8)::mtbl(4224),mtbl2(4224)
   integer::mtblp,mtblp2
   integer::ix,iwksize,is,iopt1,iopt2(2)
   integer,allocatable::lis(:),lis2(:)
   integer,allocatable::lstpu(:),lstpu2(:)
   integer(kind=4) :: iflag=0
#elif MKL
   use mkl_vsl_type !2020-02-14
   use mkl_vsl
   integer(kind=4),parameter :: iseed = 7777777
   type(vsl_stream_state),save :: stream
   integer(kind=4) :: iflag=0, my_rank_seed
!2020-02-15s
   integer(kind=1), allocatable::memptr(:) !2020-02-16
   integer(kind=4) :: memsize !2020-02-16
   integer(kind=4) :: brng
   integer(kind=4) :: method=VSL_RNG_METHOD_UNIFORM_STD
!2020-02-15e
   real(kind=8) :: lb,rb
#elif MSMT
   use mt_stream
   type(mt_state),save  :: mts, mts_new
   integer(kind=4) :: iseed=1
   integer(kind=4) :: iflag=0, id
#elif DCMT
   integer(kind=4) :: w,p,id
   integer(kind=4) :: seed0,seed1
   integer(kind=4) :: iflag=0
#elif FUJITSU
! Fujitsu START
   integer(4),external :: omp_get_thread_num
!   integer(kind=4),parameter :: iseed = 7777777, NRNDMAX = lrzphi, NMAX = 10000, NUMT = 1
   integer(kind=4),parameter :: iseed = 7777777, NRNDMAX = max(lrzphi,marker_each), NMAX = 10000, NUMT = lpara
   integer(kind=4) :: IX, NBUF, NWORK, ICON
   integer(kind=4) :: iflag=0, mythread
   double precision DA(NRNDMAX,NUMT)
   double precision DWORK(NMAX,NUMT)
! Fujitsu END
#elif AURORA
   use asl_unified
   integer :: rng
   integer(kind=4) :: iflag=0 !2020-07-28
#else
   integer::nsize
   integer,allocatable::seed(:)
   real(8),allocatable::rand_no_f95(:)
   integer::iflag=0 !2017-01-23
#endif
end module

module check
   use parameters
   integer, parameter::lphi_n=2,lphi_n_size=(lphi_n+1)*2 !2025-09-03
   integer, parameter::l_kin_species=2
   integer, parameter::lkinetic=5
   integer::loss1=0,loss2=0
   integer::loss1_total,loss2_total
#ifdef KTI
   real(8)::e_trans(lkinetic),de_trans(lkinetic),e_trans1(lkinetic),e_trans2(lkinetic)
#else
   real(8)::e_trans(0:l_kin_species),de_trans(0:l_kin_species),e_trans1(0:l_kin_species),e_trans2(0:l_kin_species)
#endif
   real(8)::e_dissp,de_dissp,e_dissp1,e_dissp2
   real(8)::cos_phi(0:lphi_n,lphi),sin_phi(0:lphi_n,lphi)
end module

#ifdef TIMER
include 'mod_timer.f90' !2025-09-06
#endif

! #ifdef AMDGPU
! include 'erfc.f90'
! #endif
include 'optimized_mod25xeon+amd.f90' !2025-09-06

!--------------------------------------------------------------------
program mega2025
! MHD + Energetic Alpha Particles Simulation Code

! modified on 2025-09-07
!   giga2025_open2.f90 with kinetic thermal ions is unified in this version
!   preprecessor KTI for kinetic thermal ions  

! modified on 2025-09-06
!   use optimized_mod25XEON+AMD.f90, module amdgpu
!   optimized for Subsystem B of new Plasma Simulator with MI-300A

! modified on 2025-09-01
!   preprocessor EXPEQ for experimental equilibrium
  
! modified on 2025-07-11
!   use optimized_mod25XEON.f90, module xeon6900p
!   optimized for Subsystem A of new Plasma Simulator with XEON6900P
  
! modified on 2025-04-27
!   use optimized_mod25MAY15.f90
!   mu*grad (B_eq - B) is considered
  
! modified on 2025-04-06
!   MHD pressure is solved to avoid a numerical instability
!   from mega2025_v3dftest.f90

! modified on 2025-04-04
!   for benchmark case investigated in [Todo, PPCF 63 (2021) 075018]
!   from mega2025_v2.f90

! modified on 2025-02-05
!   for spherical tokamaks
!   from iter_mhd2vb_512PS.f90

! copyright: Y. Todo (National Institute for Fusion Science, Japan)

!--------------------------------------------------------------------
! grid points:
!      phi-direction: periodic
!      r,z-directions: 
!        my_rank_r=0:  1 2       lr-5 lr-4 lr-3 lr-2 lr-1 lr
!                      ^ (net start)
!        my_rank_r=1:             1    2    3    4
!                                           ^ (net start)
!
!        my_rank_r=i:                      lr-3 lr-2 lr-1 lr
!        my_rank_r=i+1:                     1    2    3    4
!                                                     ^ (net start)
!        my_rank_r=mpi_proc_r-2: lr-5 lr-4 lr-3 lr-2 lr-1 lr
!        my_rank_r=mpi_proc_r-1:  1    2    3    4    5    6
!                                                     ^ (net start)
!--------------------------------------------------------------------
      use mpiset
      use job_cntrol
      use grid
      use field
      use particle
      use gyro
      use random_num
      use check

#ifdef XEON6900P
      use xeon6900p
#elif AMDGPU
      use amdgpu
#elif AURORA
      use aurora
#elif FUGAKU
      use fugaku
#else
      use fx100
#endif
#ifdef TIMER
      use mod_timer
#endif

      implicit none
      real(8)::etlim !2012-02-24
      integer::kstep0,istep
      logical::flag_FLR,flag_HM,flag_BDP !2013-05-02

      namelist /param0/job_seq
      namelist /param1/kstep,ksmax,etlim !2012-02-24
      namelist /param2/kwchk,kwout,ksnapsh,kmovie
      namelist /param3/dt,phimode
      namelist /param4/nu0,eta0,nu_n0,chi0
      namelist /param5/flag_FLR,flag_HM,flag_BDP
      namelist /param6/type_a,clambda0,dclambda !2025-02-05
      namelist /param7/hpower,sd_accl,tmax0,t_classic,t_mhd !2013-10-18
      namelist /param8/beta_a0,scale_psi,valpha,temp_a !2025-04-04
      real(8)::dprep,dstart,dloop,dend,tstart,tloop,tend
      real(8)::t0,t1,t2 !2025-07-11

! start mpi

      call start_mpi

      call wall_clock(dprep)

      read(1,param0)
      read(1,param1)
      read(1,param2)
      read(1,param3)
      read(1,param4)
      read(1,param5)
      read(1,param6)
      read(1,param7)
      read(1,param8)

! --- input parameters ---
!      JOB_SEQ=1
!      KSTEP=0
!      KSMAX=100
!      KWCHK=10
!      KWOUT=10000
!      KSNAPSH=10000
!      KMOVIE=1000

!      DT=5.0D-2
!      BETA_A0=1.0D-2
!      NU0=1.0d-6
!      ETA0=1.0d-6
!      CHI0=1.0d-6
!      FLAG_FLR=.TRUE.
!      FLAG_FLR=.FALSE.
      flag_classic = .false.
      flag_stored  = .false.
! --- input parameters end ---

!      if(flag_HM)then
        psivac = 4.0d-2
!        psivac = 1.0d-1 !2025-02-05 for edge stability
!      else
!        psivac = 0.0d0
!      end if


      sqrtpsi_chi_enh = 0.95d0
      dsqrtpsi_chi_enh = 5.0d-2
      
      sqrtpsi_pe = 0.90d0
      psi_edge = 0.90d0**2 !psi, not sqrtpsi

      if(flag_FLR.and.my_rank.eq.0)then
        write(7,'(a)')'GK FLR effects both for alpha and D beam'
      end if

!      kstep = 0
!      ksmax = 100
!      kwchk = 1
!      kwout = ksmax
!      ksnapsh = ksmax
      ksmth = 1000
      kperturb = 0
      kharmonics = 1000
      kscatter = 1

      t = 0.0d0
!      dt= 1.0d2

!      call start_ransuu

! read equilibrium data created with EFIT and ASTRA
! the data should be converted for MEGA in advance
      call equi_solution
      call sgrid
      call particle_parameters
      call equilibrium
      call flux_coordinates

! electrons
!      call initial_particle(marker_e0,marker_e,m_e,e_e &
!        ,type_e,temp_e,vbeam,deltav,clambda0,dclambda &
!        ,vmin_e,vmax_e &
!        ,gc_e,v_e,lambda_e,node_e &
!        ,cf_pphi_e,pphi_min_e,pphi_max_e,total_energy_e,0)

#ifdef KTI
! ions
      if(my_rank.eq.0)then
        write(7,'(/a)')'load thermal ion particles'
      end if

      call initial_particle(marker_i0,marker_i,m_i,e_i &
        ,type_i,temp_i,vbeam,deltav,clambda0,dclambda &
        ,vmin_i,vmax_i &
        ,gc_i,v_i,lambda_i,node_i &
        ,cf_pphi_i,pphi_min_i,pphi_max_i,total_energy_i,1)
#endif
      
! alpha particles
      if(flag_BDP)then
        type_a =-5
        if(my_rank.eq.0)then
          write(7,'(/a)')'read beam deposition profile'
        end if

        call beam_deposit

      else

        if(my_rank.eq.0)then
          write(7,'(/a)')'load alpha particles'
        end if

        call initial_particle(marker_a0,marker_a,m_a,e_a &
          ,type_a,temp_a,valpha,deltav,clambda0,dclambda &
          ,vmin_a,vmax_a &
          ,gc_a,v_a,lambda_a,node_a &
          ,cf_pphi_a,pphi_min_a,pphi_max_a,total_energy_a,2)
      end if

! remove initial discrapancy for fluid part

      if(type_a.lt.0)then
        dt = omega_a/1.2d9 !dt*1.1e6/omega_a = 1ms
      end if

      call initial_mhd_balance(flag_HM)

      if(type_a.eq.-5)then
        call injection
      end if

! read data
      if(kstep.eq.0)then

        if(type_a.ge.0)then
!          if(perturb_type.eq.'random')then
!            call perturb_fluid_random(flag_HM)
!          else
            call perturb_fluid(flag_HM)
!          end if
        end if

      else
        call read_data
      end if


      call write_snapshot
!       call write_movie

! density and field
      flag_stored  = .true.
      call density_particle(flag_FLR)

#ifdef KTI
      call write_check_kti
#else      
      call write_check
#endif
      
      call harmonics

#ifdef KTI
          call e_field_kti(1)
#else
          call e_field(1)
#endif

!--------------- loop start ------------------ 
      call wall_clock(dstart)

      if(type_a.lt.0)then
        call switch_classic(flag_FLR)
      end if

      if(.not.flag_classic)then
        flag_stored  = .false.
        call density_particle(flag_FLR)
      end if

      kstep0 = kstep + 1

step: do kstep = kstep0, ksmax

      if(mod(kstep-1,10).eq.0)then
#ifdef KTI
        call order_particle(marker_i,gc_i,node_i)
#endif
        call order_particle(marker_a,gc_a,node_a)
      end if


! satellite should be called after order_particle
        if(flag_FLR)then
          call satellite
        end if

rkg:    do istep = 1, 4

          call push_particle(flag_FLR,istep)

          if(.not.flag_classic)then
#ifdef KTI
             call mhd_kti(istep,flag_HM) !subroutines mhd and hm are unified
#else             
             call mhd(istep,flag_HM) !subroutines mhd and hm are unified
#endif
          end if

          if(.not.flag_classic)then
            call density_particle(flag_FLR)

#ifdef KTI
            call e_field_kti(1)
#else             
            call e_field(1)
#endif
          end if

        end do rkg

!collision and injection
!this subroutine should be called before com_particle to refer to v_a(n)
        if(type_a.eq.-5)then
          call injection
          call scattering
        end if
! communication

!         call com_particle(marker_e,gc_e,node_e,node_now_e)
#ifdef KTI
         call com_particle(marker_i,gc_i,node_i,node_now_i)
#endif
         call com_particle(marker_a,gc_a,node_a,node_now_a)

!        t = dt*dble(kstep)
        t = t + dt

        if(type_a.lt.0)then
          call switch_classic(flag_FLR)
        end if

        if(.not.flag_classic)then
        if(mod(kstep,ksmth).eq.0)then
          call mhd_lowpass(flag_HM)
          call mhd_smoothing(flag_HM)

          flag_stored = .false.
          call density_particle(flag_FLR)

#ifdef KTI
          call e_field_kti(1)
#else
          call e_field(1)
#endif
        end if
        end if

        if(mod(kstep,kwout).eq.0)then
          call write_data
        end if 

        if(mod(kstep,ksnapsh).eq.0)then
          call write_snapshot
        end if 

!        if(mod(kstep,kmovie).eq.0)then
!          call write_movie
!        end if 


        if(mod(kstep,kwchk)*mod(kstep,kharmonics).eq.0)then
          flag_stored  = .true.
          call density_particle(flag_FLR)

          if(mod(kstep,kwchk).eq.0)then
#ifdef KTI
            call write_check_kti
#else     
            call write_check
#endif
          end if

          if(mod(kstep,kharmonics).eq.0)then
            call harmonics
          end if 

          if(.not.flag_classic)then
            flag_stored  = .false.
            call density_particle(flag_FLR)
          end if

        end if

        if(mod(kstep,1000).eq.0)then
          call wall_clock(dloop)
          if((etlim-(dloop-dprep)).lt.9.0d2)exit
        end if 

      end do step

9999  continue

      call wall_clock(dloop)
!--------------- loop end ------------------

        if(kstep.gt.ksmax)then
          kstep = ksmax
        end if

! output the last data

        if(mod(kstep,kwout).ne.0)then
          call write_data
        end if 

        if(mod(kstep,ksnapsh).ne.0)then
          call write_snapshot
        end if 

!        if(mod(kstep,kmovie).ne.0)then
!          call write_movie
!        end if 

        if(.not.flag_classic)then
          flag_stored  = .true.
          call density_particle(flag_FLR)
        end if

        if(mod(kstep,kwchk).ne.0)then
#ifdef KTI
            call write_check_kti
#else     
            call write_check
#endif
        end if

        if(mod(kstep,kharmonics).ne.0)then
          call harmonics
        end if 

      call wall_clock(dend)

      if(my_rank.eq.0)then
        tstart= dstart - dprep
        tloop = dloop - dstart
        tend  = dend  - dloop

        write(7,'()')
        write(7,'(a,i9)')'number of injected particles=',n_inj
        write(7,'(a,1pe14.5)')'tstart=',tstart
        write(7,'(a,1pe14.5)')'tloop =',tloop
        write(7,'(a,1pe14.5)')'tend  =',tend
      end if

#ifdef TIMER
      call print_timer()
#endif

! end random number generator
      call end_ransuu

! end mpi
      call end_mpi

end
!--------------------------------------------------------------------
subroutine switch_classic(flag_FLR)
!--------------------------------------------------------------------
      use parameters
      use mpiset
      use field
      use particle, only:t_classic,t_mhd
      use grid, only:t,dt,kstep,flag_classic,flag_stored
      implicit none
      logical::flag_FLR

      if(mod(t/omega_a*1.0d3, (t_classic + t_mhd) ).lt.t_classic)then
        if(.not.flag_classic)then
          if(my_rank.eq.0)then
            write(7,'(a,i8,a,f7.3,2a)')'kstep=',kstep,'  t=',t/omega_a*1.0d3,' ms ' &
                                      ,'  switch to classical simulation'
          end if
          rho = rho0
          prs = prs0
          vr = 0.0d0
          vz = 0.0d0
          vphi = 0.0d0
          ur = 0.0d0
          uz = 0.0d0
          uphi = 0.0d0
          er = 0.0d0 !2013-06-12
          ez = 0.0d0 !2013-06-12
          ephi = 0.0d0 !2013-06-12
          epara = 0.0d0 !2013-06-12
          br = br0
          bz = bz0
          bphi = bphi0
          enrg = enrg0 !2021-05-05

          vr_ref = vr
          vz_ref = vz
          vphi_ref = vphi
          br_ref = br
          bz_ref = bz
          bphi_ref = bphi
          rho_ref = rho
          prs_ref = prs
          enrg_ref = enrg !2021-05-05
        end if

        flag_classic = .true.
        flag_stored  = .true.
        dt = omega_a/3.0d7 !dt*3e5/omega_a = 10ms

      else
        if(flag_classic)then
          if(my_rank.eq.0)then
            write(7,'(a,i8,a,f7.3,2a)')'kstep=',kstep,'  t=',t/omega_a*1.0d3,' ms ' &

                                      ,'  switch to full simulation'
          end if

        flag_stored  = .false.
        call density_particle(flag_FLR)

#ifdef KTI
          call e_field_kti(1)
#else
          call e_field(1)
#endif

        end if

        flag_classic = .false.
        dt = omega_a/1.2d9 !dt*1.1e6/omega_a = 1ms !2014-09-17
!        dt = omega_a/1.1d9 !dt*1.1e6/omega_a = 1ms !2014-03-29
!2014-03-29        dt = omega_a/9.0d8 !dt*9e5/omega_a = 1ms
!2013-11-19        dt = omega_a/9.0d8*2.0d0 !dt*9e5/omega_a = 2ms !2013-10-18 for d3d4_64v3.f90

      end if

end
!--------------------------------------------------------------------
subroutine density_particle(flag_FLR)
!--------------------------------------------------------------------
      use parameters
      use mpiset
      use field
      use particle
      use gyro
      use grid

#ifdef XEON6900P
      use xeon6900p
#elif AMDGPU
      use amdgpu
#elif AURORA
      use aurora
#elif FUGAKU
      use fugaku
#else
      use fx100
#endif
#ifdef TIMER
      use mod_timer
#endif
      implicit none
      logical::flag_FLR


      if(flag_FLR)then
        call satellite

        if(flag_stored)then
#ifdef KTI
          call moments_gyro(marker_i,marker_i_gyro,m_i,gc_i,gyro_i & !2016-02-04
                    ,dns_i,mom_i,ppara_i,pperp_i,qpara_i,qperp_i &
                    ,dns_i0,mom_i0,ppara_i0,pperp_i0,qpara_i0,qperp_i0)
#endif
          call moments_gyro(marker_a,marker_a_gyro,m_a,gc_a,gyro_a & !2016-02-04
                    ,dns_a,mom_a,ppara_a,pperp_a,qpara_a,qperp_a &
                    ,dns_a0,mom_a0,ppara_a0,pperp_a0,qpara_a0,qperp_a0)

        else
#ifdef KTI
           call density_gyro(marker_i,marker_i_gyro,m_i,gc_i,gyro_i & !2022-01-11
                    ,dns_i,mom_i,ppara_i,pperp_i &
                    ,dns_i0,mom_i0,ppara_i0,pperp_i0)
#endif
           call density_gyro(marker_a,marker_a_gyro,m_a,gc_a,gyro_a & !2022-01-11
                    ,dns_a,mom_a,ppara_a,pperp_a &
                    ,dns_a0,mom_a0,ppara_a0,pperp_a0)

        end if

      else
        if(flag_stored)then
#ifdef KTI
          call moments(marker_i,m_i,gc_i &
                      ,dns_i,mom_i,ppara_i,pperp_i,qpara_i,qperp_i &
                      ,dns_i0,mom_i0,ppara_i0,pperp_i0,qpara_i0,qperp_i0)
#endif
          call moments(marker_a,m_a,gc_a &
                      ,dns_a,mom_a,ppara_a,pperp_a,qpara_a,qperp_a &
                      ,dns_a0,mom_a0,ppara_a0,pperp_a0,qpara_a0,qperp_a0)
        else

#ifdef KTI
          call density(marker_i,m_i,gc_i &
                      ,dns_i,mom_i,ppara_i,pperp_i &
                      ,dns_i0,mom_i0,ppara_i0,pperp_i0)
#endif
          call density(marker_a,m_a,gc_a &
                      ,dns_a,mom_a,ppara_a,pperp_a &
                      ,dns_a0,mom_a0,ppara_a0,pperp_a0)
        end if
      end if

! for electrons, 2022-04-09
!        if(flag_stored)then
!          call moments(marker_e,m_e,gc_e &
!                      ,dns_e,mom_e,ppara_e,pperp_e,qpara_e,qperp_e &
!                      ,dns_e0,mom_e0,ppara_e0,pperp_e0,qpara_e0,qperp_e0)

!        else
!          call density(marker_e,m_e,gc_e &
!                      ,dns_e,mom_e,ppara_e,pperp_e &
!                      ,dns_e0,mom_e0,ppara_e0,pperp_e0)
!        end if
      
end
!--------------------------------------------------------------------
subroutine push_particle(flag_FLR,istep)
!--------------------------------------------------------------------
      use parameters
      use mpiset
      use field
      use particle
      use gyro
      use grid

!2025-07-11
#ifdef XEON6900P
      use xeon6900p
#elif AMDGPU
      use amdgpu
#elif AURORA
      use aurora
#elif FUGAKU
      use fugaku
#else
      use fx100
#endif
      implicit none
      logical::flag_FLR
      integer::istep
      integer::i_FLR
      integer::e_FLR=0 !2022-04-09

      if(flag_FLR)then

!ion
#ifdef KTI
         call emf_gyro(marker_i,marker_i_gyro,gc_i,gyro_i &
                      ,flp_i,flp_i_gyro)
#endif
         
!alpha        
         call emf_gyro(marker_a,marker_a_gyro,gc_a,gyro_a &
                     ,flp_a,flp_a_gyro)
      end if

      if(flag_FLR)then
        i_FLR = 1
      else
        i_FLR = 0
      end if

!electron
!      call push(marker_e,m_e,e_e &
!               ,type_e,temp_e,vbeam,deltav,clambda0,dclambda &
!               ,gc_e,dgc_e,v_e &
!               ,cf_pphi_e,pphi_min_e,pphi_max_e &
!               ,flp_e,0,e_FLR)

!ion
#ifdef KTI
      call push(marker_i,m_i,e_i &
               ,type_i,temp_i,vbeam,deltav,clambda0,dclambda &
               ,gc_i,dgc_i,v_i &
               ,cf_pphi_i,pphi_min_i,pphi_max_i &
               ,flp_i,1,i_FLR)
#endif
      
!alpha
      call push(marker_a,m_a,e_a &
               ,type_a,temp_a,valpha,deltav,clambda0,dclambda &
               ,gc_a,dgc_a,v_a &
               ,cf_pphi_a,pphi_min_a,pphi_max_a &
               ,flp_a,2,i_FLR)


!      call t_integration(istep,type_e,marker_e &
!               ,gc_e,dgc_e,gc1_e,gc2_e)

#ifdef KTI
      call t_integration(istep,type_i,marker_i &
               ,gc_i,dgc_i,gc1_i,gc2_i)
#endif
      call t_integration(istep,type_a,marker_a &
                ,gc_a,dgc_a,gc1_a,gc2_a)

!      call bcptcl(marker_e,gc_e,node_e)
#ifdef KTI
      call bcptcl(marker_i,gc_i,node_i)
#endif
      call bcptcl(marker_a,gc_a,node_a)

end
!--------------------------------------------------------------------
subroutine read_data
!--------------------------------------------------------------------
      use parameters
      use job_cntrol
      use mpiset
      use field
      use particle
      use grid
      use check, only:e_trans,e_dissp
      use random_num !2013-05-02

      implicit none
      integer::kst,i
      character(len=64)::cfile
      integer, save::i_open=0
      integer::job_prev
      real(8)::dr1,dz1,dphi1
      integer::n,ia,ia1,ja,ja1,ka,ka1
      real(8)::ar,ar1,az,az1,aphi,aphi1
      real(8)::aaa1,aaa2,aaa3,aaa4,aaa5,aaa6,aaa7,aaa8
      real(8)::bre,bze,bphie,babse
      integer::ierr,i_allocate=0 !2020-02-16, MKL

!      save i_open
!       data i_open / 0 /

      cfile = 'job_num-xxxxx.pd'
      job_prev = job_seq - 1
      write(cfile(5:7),'(i3.3)') job_prev
      write(cfile(9:13),'(i5.5)') my_rank

      if(i_open.eq.0)then
        open(10,file=cfile,form='unformatted')
        i_open = 1
      end if

1000  continue
      read(10,end=9999)kst,t,e_trans,e_dissp,n_inj,co_stored,cntr_stored & !2017-04-04
           ,rho,prs,vr,vz,vphi,br,bz,bphi & 
           ,dns_e,mom_e,ppara_e,pperp_e &
           ,dns_i,mom_i,ppara_i,pperp_i &
           ,dns_a,mom_a,ppara_a,pperp_a &
           ,gc_e,node_e &
           ,gc_i,node_i &
           ,gc_a,node_a
!2013-05-02
#ifdef MSMT
      call read(mts,10)

#elif FUJITSU
      read(10)DWORK
      IX = 0
      iflag = 1

!2020-02-16s
#elif MKL
      read(10)memsize
      if(i_allocate.eq.0)then
        allocate(memptr(memsize) )
        i_allocate = 1
      end if
      read(10)memptr
!2020-02-16s

#endif
!2013-05-02 end

      if(kst.lt.kstep)then
        go to 1000
      end if

9999  continue

!2020-02-16s
#ifdef MKL
      ierr = vslloadstreamm( stream, memptr )
      iflag = 1
      deallocate(memptr)
#endif
!2020-02-16e

      kstep = kst

      babs = max(eps_b, sqrt(br**2 + bz**2 + bphi**2) )
      enrg = 0.50d0*(vr**2 + vz**2 + vphi**2)/rho &
           + 0.50d0*(br**2 + bz**2 + bphi**2) &
           + prs/(gamma - 1.0d0)

      call gradmg

#ifdef KTI
          call e_field_kti(1)
#else
          call e_field(1)
#endif

      do i = 1, lrzphi
        br_ref(i,1,1) = br(i,1,1)
        bz_ref(i,1,1) = bz(i,1,1)
        bphi_ref(i,1,1) = bphi(i,1,1)

        vr_ref(i,1,1) = vr(i,1,1)
        vz_ref(i,1,1) = vz(i,1,1)
        vphi_ref(i,1,1) = vphi(i,1,1)

        rho_ref(i,1,1) = rho(i,1,1)
        prs_ref(i,1,1) = prs(i,1,1)
        enrg_ref(i,1,1) = enrg(i,1,1)
      end do


! set v_i(n), v_a(n)
      call total_v(marker_i,m_i,gc_i,v_i)
      call total_v(marker_a,m_a,gc_a,v_a)

end
!--------------------------------------------------------------------
subroutine total_v(marker_num,amassp,gc,v)
!--------------------------------------------------------------------
      use parameters
      use mpiset
      use field
      use grid

      implicit none
      integer::marker_num
      real(8)::amassp
      real(8)::gc(ngc2,marker_num)
      real(8)::v(marker_num)

      real(8)::dr1,dz1,dphi1
      integer::n,ia,ia1,ja,ja1,ka,ka1
      real(8)::ar,ar1,az,az1,aphi,aphi1
      real(8)::aaa1,aaa2,aaa3,aaa4,aaa5,aaa6,aaa7,aaa8
      real(8)::bre,bze,bphie,babse

      dr1 = 1.0d0/dr
      dz1 = 1.0d0/dz
      dphi1 = 1.0d0/dphi

!$omp parallel do private(ia,ia1,ja,ja1,ka,ka1,ar,ar1,az,az1,aphi,aphi1 &
!$omp&      ,aaa1,aaa2,aaa3,aaa4,aaa5,aaa6,aaa7,aaa8,bre,bze,bphie,babse)
      do n = 1, marker_num
        ia = max(1, min(lr-1, &
                 int( ( gc(1,n)-(major_r - minor_r) )*dr1) + 1 &
                      - kr_offset(my_rank) ) )
        ia1= ia + 1
        ja = max(1, min(lz-1, &
                 int( gc(2,n)*dz1 ) + 1 - kz_offset(my_rank) ) )
        ja1= ja + 1
        ka = max(1, min(lphi-1, &
                 int( gc(3,n)*dphi1 ) + 3 - kphi_offset(my_rank) ) )
        ka1= ka + 1

        ar1 = (gc(1,n)-grr(ia,ja,ka) ) *dr1
        ar  = 1.0d0 - ar1
        az1 = (gc(2,n)-gzz(ia,ja,ka) ) *dz1
        az  = 1.0d0 - az1
        aphi1 = (gc(3,n)-gphi(ia,ja,ka) ) *dphi1
        aphi  = 1.0d0 - aphi1

        aaa1 = ar *az *aphi
        aaa2 = ar1*az *aphi
        aaa3 = ar *az1*aphi
        aaa4 = ar1*az1*aphi
        aaa5 = ar *az *aphi1
        aaa6 = ar1*az *aphi1
        aaa7 = ar *az1*aphi1
        aaa8 = ar1*az1*aphi1

! fields at each particle position

        bre = br(ia, ja, ka )*aaa1 + br(ia1,ja, ka )*aaa2 &
            + br(ia, ja1,ka )*aaa3 + br(ia1,ja1,ka )*aaa4 &
            + br(ia, ja, ka1)*aaa5 + br(ia1,ja, ka1)*aaa6 &
            + br(ia, ja1,ka1)*aaa7 + br(ia1,ja1,ka1)*aaa8
        bze = bz(ia, ja, ka )*aaa1 + bz(ia1,ja, ka )*aaa2 &
            + bz(ia, ja1,ka )*aaa3 + bz(ia1,ja1,ka )*aaa4 &
            + bz(ia, ja, ka1)*aaa5 + bz(ia1,ja, ka1)*aaa6 &
            + bz(ia, ja1,ka1)*aaa7 + bz(ia1,ja1,ka1)*aaa8
        bphie = bphi(ia, ja, ka )*aaa1 + bphi(ia1,ja, ka )*aaa2 &
            + bphi(ia, ja1,ka )*aaa3 + bphi(ia1,ja1,ka )*aaa4 &
            + bphi(ia, ja, ka1)*aaa5 + bphi(ia1,ja, ka1)*aaa6 &
            + bphi(ia, ja1,ka1)*aaa7 + bphi(ia1,ja1,ka1)*aaa8
        babse = sqrt(bre**2 + bze**2 + bphie**2)

!   gc_x(:,n) = 1:r, 2:z, 3:phi, 4:p, 5:mu, 6:weight, 7:active, 8:p_phi, 9:fff, 10:fnrml
        v(n) = sqrt(2.0d0*gc(5,n)*babse/amassp + (gc(4,n)/amassp)**2)
      end do


end
!--------------------------------------------------------------------
subroutine write_data
!--------------------------------------------------------------------
      use job_cntrol
      use mpiset
      use field
      use particle
      use grid
      use check, only:e_trans,e_dissp
      use random_num !2013-05-02
      implicit none

      character(len=64)::cfile
      integer, save::i_open=0
      integer::ierr !2020-02-16, MKL

!      save i_open
!       data i_open / 0 /

      cfile = 'job_num-xxxxx.pd'
      write(cfile(5:7),'(i3.3)') job_seq
      write(cfile(9:13),'(i5.5)') my_rank

      if(i_open.eq.0)then
        open(20,file=cfile,form='unformatted')
        i_open = 1
      end if

      write(20)kstep,t,e_trans,e_dissp,n_inj,co_stored,cntr_stored & !2017-04-04
           ,rho,prs,vr,vz,vphi,br,bz,bphi & 
           ,dns_e,mom_e,ppara_e,pperp_e &
           ,dns_i,mom_i,ppara_i,pperp_i &
           ,dns_a,mom_a,ppara_a,pperp_a &
           ,gc_e,node_e &
           ,gc_i,node_i &
           ,gc_a,node_a
!2013-05-02
#ifdef MSMT
      call save(mts,20)

#elif FUJITSU
      write(20)DWORK

!2020-02-16s
#elif MKL
      memsize = vslgetstreamsize( stream ) 
      allocate(memptr(memsize) )
      ierr = vslsavestreamm( stream, memptr )
      write(20)memsize
      write(20)memptr
      deallocate(memptr)
!2020-02-16s

#endif
!2013-05-02 end

      flush(20)

end
!--------------------------------------------------------------------
subroutine write_snapshot
! write out data on a poloidal plane for snapshot
!--------------------------------------------------------------------
      use mpiset
      use field
      use grid
      implicit none

      integer::kp
      real(8)::br_pol(lrnet,lznet)
      real(8)::bz_pol(lrnet,lznet)
      real(8)::bphi_pol(lrnet,lznet)
      real(8)::vr_pol(lrnet,lznet)
      real(8)::vz_pol(lrnet,lznet)
      real(8)::vphi_pol(lrnet,lznet)
      real(8)::er_pol(lrnet,lznet)
      real(8)::ez_pol(lrnet,lznet)
      real(8)::ephi_pol(lrnet,lznet)
      real(8)::epara_pol(lrnet,lznet)
      real(8)::rho_pol(lrnet,lznet)
      real(8)::prs_pol(lrnet,lznet)
      real(8)::ppara_i_pol(lrnet,lznet)
      real(8)::pperp_i_pol(lrnet,lznet)
      real(8)::ppara_a_pol(lrnet,lznet)
      real(8)::pperp_a_pol(lrnet,lznet)

! global toroidal grid number to plot
      kp = 3

      call data_snapshot(br,br_pol,kp)
      call data_snapshot(bz,bz_pol,kp)
      call data_snapshot(bphi,bphi_pol,kp)
      call data_snapshot(vr,vr_pol,kp)
      call data_snapshot(vz,vz_pol,kp)
      call data_snapshot(vphi,vphi_pol,kp)
      call data_snapshot(er,er_pol,kp)
      call data_snapshot(ez,ez_pol,kp)
      call data_snapshot(ephi,ephi_pol,kp)
      call data_snapshot(epara,epara_pol,kp)
      call data_snapshot(rho,rho_pol,kp)
      call data_snapshot(prs,prs_pol,kp)
      call data_snapshot(ppara_i,ppara_i_pol,kp)
      call data_snapshot(pperp_i,pperp_i_pol,kp)
      call data_snapshot(ppara_a,ppara_a_pol,kp)
      call data_snapshot(pperp_a,pperp_a_pol,kp)

      if(my_rank.eq.0)then
        write(30)kstep,t,br_pol,bz_pol,bphi_pol ,vr_pol,vz_pol&
             &,vphi_pol ,er_pol,ez_pol,ephi_pol,epara_pol ,rho_pol&
             &,prs_pol ,ppara_i_pol,pperp_i_pol,ppara_a_pol,pperp_a_pol

        flush(30)
      end if

      end
!--------------------------------------------------------------------
subroutine particle_parameters
!--------------------------------------------------------------------

      use particle
      use field
      use grid
      use mpiset
      
      implicit none
      integer::i

!Fredricson(2009), 2025-02-05
!      bunit=0.45d0 ![T]
!      va=1.07d6 ![m/s]
!      omega_a=eunit*bunit/munit
!      rho_a=va/omega_a

!      valpha=1.76d6/va
!      vcrit0=1.00d6/va
#ifdef EXPEQ
#else
      valpha = 1.2d0 !2025-04-024, Todo (2021)
      vcrit0 = 0.50d0 !2025-04-024, Todo (2021)
#endif
      deltav = valpha*1.0d-1
      
! type=0: maxwellian, type=1: slowing down, type=2: beam with mu=0
! type=3: anisotropic beam specified with clambda0 & dclambda 
! type=-5: full-f and beam injection with flag_BDP, 2013-05-2

      type_e = 0
      type_i = 0
!      type_a = 1

! temperature

      temp_e = 1.0d-2
      temp_i = 1.0d-2 !2025-02-05
      temp_a = 1.0d0
!      beta0  = 1.0d-3
!      beta_a0= 6.0d-4
!      beta_a0= 6.0d-2
!      n_a0 = n0*beta_a0/beta0*(temp_e + temp_i*abs(e_e)/e_i )/temp_a

! alfven velocity: va0=1 set in module parameter

!      va0 = b0/sqrt(m_ti*n0*abs(e_e)/e_i)

! conductivity, viscosity, resistivity, specific heat ratio

      nu_n = nu_n0*va0*major_r
      nun_enhanced = 1.0d-5*va0*major_r
      chi = chi0*va0*major_r
      chi_enhanced = 1.0d-5*va0*major_r
      nu  = nu0*va0*major_r
      nu_vac = nu
!      nu_vac = 1.0d-5*va0*major_r
      eta = eta0*va0*major_r
      gamma = 5.0d0/3.0d0

! beam injection parameters: 2013-05-02
! normalize by va,omega_a,rho_a in SI
      tmax = tmax0*1.0d-3*omega_a !tmax0 [ms]
!      sd_time = sd_time0*1.0d-3*omega_a !sd_time0 [ms]
!      nu_sd = 1.0d0/sd_time !nu_sd is not used 2013-05-02, arary nusd is employed

!      nu_krook = 1.0d4/omega_a ! [1/s]
!      nu_krook = 3.0d-2/major_r ! 3% of Alfven frequency (vA/R0), 2025-02-05
      nu_krook = 0.0d-2/major_r ! =0, 2025-02-05

!      bunit = omega_a*munit/eunit !defined in subroutine equi_solution
      zeff = 1.4d0 !for DIII-D 142111.525

! marker particle loading: minimum and maximum velocity

      vmin_e = 0.0d0
! for numerical stability
      vmax_e = sqrt(temp_e/m_e)
!      vmax_e = 3.0d0*sqrt(temp_e/m_e)

      vmin_i = 0.0d0
      vmax_i = 3.0d0*sqrt(3.0d0*temp_i/m_i)

! vbfactor controls vbeam
!      vbeam = vbeam * vbfactor
      vbeam = valpha !vbeam is not used

      if(type_a.eq.0)then
        vmin_a = 0.0d0
        vmax_a = 3.0d0*sqrt(3.0d0*temp_a/m_a)
      else
        vmin_a = valpha*1.0d-1
        vmax_a = valpha + 3.0d0*deltav
      end if

! The following cf_pphi and pphi_min are not used

      cf_pphi_e(0) = 1.0d0
      cf_pphi_i(0) = 1.0d0
      cf_pphi_a(0) = 1.0d0
      do i = 1, lcfpphi
        cf_pphi_e(i) = 0.0d0
        cf_pphi_i(i) = 0.0d0
        cf_pphi_a(i) = 0.0d0
      end do

      pphi_min_e = 0.0d0
      pphi_max_e = 1.0d0
!      pphi_min_i = 0.0d0
!      pphi_max_i = 1.0d0
!      pphi_min_a = 0.0d0
!      pphi_max_a = 1.0d0


! pitch-angle sattering rate
!=4.46d-3*z_eff*z_alpha**2/sqrt(a_alpha)*ln_lambda*(E_alpha/MeV)**(-3/2)*(ne/1d20 m**(-3)) sec**-1
! major_r=2.53m, va=1.27d7 m/s, z_eff=2.44, z_alpha=2, a_alpha=4,valpha(3.5MeV)=1.02,ne0=4d19
!ln_lamdba is multiplied in subroutine equilibrium
! pas = pas0*ln_lambda*(ne/ne0)/v**3

!      ln_lambda = 3.91d1 - 1.15d0*log10(4.0d19) + 2.3d0*log10(6.0d0)
!      pas0 = 4.46d-3*(2.53d0/1.27d7)*2.44d0*2.0d0*1.02d0**3*4.0d-1

end
!--------------------------------------------------------------------
subroutine sgrid
!--------------------------------------------------------------------
      use grid
      use mpiset
      use check
      implicit none
      real(8)::x,y
      integer::i,j,k,n

#ifdef EXPEQ
      minor_r = rleng*0.50d0
#else
      major_r = minor_r*aspect
      rleng = 2.0d0*minor_r
      zleng = rleng*dkappa
#endif      
      phileng = twopi/phimode

      dr = rleng/dble(lrnet-1)
      dz = zleng/dble(lznet-1)
      dphi = phileng/dble(lphinet)

!$omp parallel do private(y,x)
      do k = 1, lphi
      do j = 1, lz
      do i = 1, lr
        grr(i,j,k) = major_r - minor_r + dr*dble(i +&
             & kr_offset(my_rank) - 1)
        gzz(i,j,k) = dz*dble(j + kz_offset(my_rank) - 1)
        gphi(i,j,k) = dphi*dble(k + kphi_offset(my_rank) - lphi_shd - 1)
!        gphi(i,j,k) = dphi*dble(k + kphi_offset(my_rank) - 3)
        y = gzz(i,j,k) - 0.50d0*zleng
        x = grr(i,j,k) - major_r
        gtht(i,j,k) = atan2(y+1.0d-20,x)
        gaa(i,j,k) = sqrt(x**2 + y**2)
      end do
      end do
      end do

!$omp parallel do
      do n = 0, lphi_n
      do k = 1, lphi
        cos_phi(n,k) = cos(dble(n)*phimode*gphi(1,1,k) ) !2012-04-17
        sin_phi(n,k) = sin(dble(n)*phimode*gphi(1,1,k) ) !2012-04-17
      end do
      end do

end
!--------------------------------------------------------------------
subroutine equi_solution
!--------------------------------------------------------------------
      use mpiset
      use equi_sol
      use grid
      use flux_c
#ifdef EXPEQ
      use particle, only:valpha,vbeam
#else
      use particle, only:cf_pphi_a,pphi_min_a,pphi_max_a
#endif
      use field, only:psimax
      implicit none
      integer::i,ir,iz,irad,itheta

#ifdef EXPEQ
      read(60)ir,iz,psimax,rleng,zleng,major_r,raxis,zaxis &
                   ,va,omega_a,rho_a,valpha &
                   ,psisim,brsim,bzsim,bhsim,ppsim,pfasim &
                   ,rhosim,dns_asim,vcritsim,nusdsim,tisim,tesim,rotsim &
                   ,irad,itheta,rpsi,gpsi,qpsi &
                   ,crdmagr,crdmagz,crdmagphi

! va, omega_a, rho_a: given in subroutine equi_solution for realistic applications
      bunit = omega_a*munit/eunit !2025-08-31

      if(ir.ne.lrnet.or.iz.ne.lznet.or.irad.ne.lpsi.or.itheta.ne.ltheta)then
        write(7,*)'equilibrium data size is different'
        write(7,*)'equil. lr, lz, lpsi, ltheta=',ir,iz,irad,itheta
        write(7,*)'simul. lrnet, lznet, lpsi, ltheta=',lrnet,lznet,lpsi,ltheta
      end if
#else      
      read(60)psi1,curden,presur,presurh_para,di2 &
             ,gsbr,gsbphi,gsbz,gsrho,gste &
             ,iaxis,draxis &
             ,gpsi,rpsi,qpsi_exp,qpsi,p_b_psi,p_h_psi,rho_psi &
             ,crdmagr,crdmagz,crdmagphi &
             ,psimax,cf_pphi_a,iter_eq,ddi2psi

      raxis = minor_r*draxis
#endif
      

!$omp parallel do
      do i = 1, lpsi
        gpsi_nrm(i) = 1.0d0 - gpsi(i)/psimax
      end do

end
!--------------------------------------------------------------------
subroutine equilibrium
!--------------------------------------------------------------------
      use mpiset
      use field
      use particle
      use grid
      use equi_sol
      use check

      implicit none

      integer::i,j,k,i1,j1
      integer::k0,k1
      real(8)::energy_e,energy_i,energy_a,beta_a_max
      real(8)::bmin_global,b2
      real(8)::sqrtpsi
      real(8)::pmax,pmin
      real(8)::y

#ifdef EXPEQ
      pphi_max_i = e_i*psimax
      pphi_min_i = 0.0d0
      pphi_max_a = e_a*psimax
      pphi_min_a = 0.0d0

!$omp parallel do private(i1,j1)
      do k = 1, lphi
      do j = 1, lz
      do i = 1, lr
       i1 = i + kr_offset(my_rank)
       j1 = j + kz_offset(my_rank)
       psi(i,j,k) = psisim(i1,j1)
       br(i,j,k) = brsim(i1,j1)
       bz(i,j,k) = bzsim(i1,j1)
       bphi(i,j,k) = bhsim(i1,j1)
!       ppara_i(i,j,k) = pnbsim(i1,j1)
!       pperp_i(i,j,k) = pnbsim(i1,j1)

       if(type_a.lt.0)then
         ppara_a(i,j,k) = 0.0d0
         pperp_a(i,j,k) = 0.0d0
         dns_a(i,j,k) = 0.0d0
       else
!         ppara_a(i,j,k) = pfasim(i1,j1)*n_a0
!         pperp_a(i,j,k) = pfasim(i1,j1)*n_a0
!         dns_a(i,j,k) = dns_asim(i1,j1)*n_a0

! scale_psi is an input parameter
        y = (1.0d0 - psi(i,j,k)/psimax)/scale_psi
        ppara_a(i,j,k) = 0.50d0*b0**2*beta_a0*exp(-y)
        pperp_a(i,j,k) = 0.50d0*b0**2*beta_a0*exp(-y)
        dns_a(i,j,k) =(ppara_a(i,j,k) + 2.0d0*pperp_a(i,j,k))/(3.0d0*temp_a) 
       end if

       prs(i,j,k) = ppsim(i1,j1)
!       prs(i,j,k) = 5.0d-3 !2013-02-26
       rho(i,j,k) = rhosim(i1,j1)
!       dns_i(i,j,k) = dns_isim(i1,j1)
#ifdef KTI
       ppara_i(i,j,k) = 0.50d0*prs(i,j,k)
       pperp_i(i,j,k) = 0.50d0*prs(i,j,k)
       dns_i(i,j,k) = rho(i,j,k)/m_i
       mom_i(i,j,k) = 0.0d0

       ppara_e(i,j,k) = ppara_i(i,j,k)
       pperp_e(i,j,k) = ppara_i(i,j,k)
       dns_e(i,j,k) =-e_i*dns_i(i,j,k)/e_e
       mom_e(i,j,k) = 0.0d0

! prs is electron pressure
       prs(i,j,k) =(ppara_e(i,j,k) + 2.0d0*pperp_e(i,j,k) )/3.0d0
#else
       ppara_i(i,j,k) = 0.0d0
       pperp_i(i,j,k) = 0.0d0
       dns_i(i,j,k) = 0.0d0
       mom_i(i,j,k) = 0.0d0

       ppara_e(i,j,k) = 0.50d0*prs(i,j,k)
       pperp_e(i,j,k) = 0.50d0*prs(i,j,k)
       dns_e(i,j,k) = rho(i,j,k)/m_i
       mom_e(i,j,k) = 0.0d0

!       prs(i,j,k) is fluid pressure
#endif
       vcrit(i,j,k) = vcritsim(i1,j1)
       nusd(i,j,k) = nusdsim(i1,j1)*sd_accl
       ti(i,j,k) = tisim(i1,j1)
       te(i,j,k) = tesim(i1,j1)
       utor(i,j,k)=rotsim(i1,j1)
      end do
      end do
      end do
#else      
! brat is corrected
      brat = 1.0d0/gsbphi(iaxis,1)
      jaxis = lznet/2 + 1

!$omp parallel do private(k0,k1)
      do k = 1, nz
        k0 =-k + nz + 1
        k1 = k + nz
      do i = 1, nr
        GSBR2(I,K0) =-GSBR(I,K)*BRAT
        GSBR2(I,K1) = GSBR(I,K)*BRAT
        GSBphi2(I,K0) = GSBphi(I,K)*BRAT
        GSBphi2(I,K1) = GSBphi(I,K)*BRAT
        GSBZ2(I,K0) = GSBZ(I,K)*BRAT
        GSBZ2(I,K1) = GSBZ(I,K)*BRAT
        PRSR2(I,K0) = PRESUR(I,K)*BRAT**2
        PRSR2(I,K1) = PRESUR(I,K)*BRAT**2
        PRSR2H_PARA(I,K0) = PRESURH_PARA(I,K)*BRAT**2
        PRSR2H_PARA(I,K1) = PRESURH_PARA(I,K)*BRAT**2
        PRSR2H_PERP(I,K0) = PRESURH_PARA(I,K)*BRAT**2
        PRSR2H_PERP(I,K1) = PRESURH_PARA(I,K)*BRAT**2
        PRSR2H_PERP2(I,K0) = PRESURH_PARA(I,K)*BRAT**2
        PRSR2H_PERP2(I,K1) = PRESURH_PARA(I,K)*BRAT**2
!        PRSR2H_PERP(I,K0) = PRESURH_PERP(I,K)*BRAT**2
!        PRSR2H_PERP(I,K1) = PRESURH_PERP(I,K)*BRAT**2
!        PRSR2H_PERP2(I,K0) = PRESURH_PERP2(I,K)*BRAT**2
!        PRSR2H_PERP2(I,K1) = PRESURH_PERP2(I,K)*BRAT**2
        GSRHO2(I,K0) = GSRHO(I,K)
        GSRHO2(I,K1) = GSRHO(I,K)
        GSTE2(I,K0) = GSTE(I,K)
        GSTE2(I,K1) = GSTE(I,K)

        PSI2(I,K0) = PSI1(I,K)*BRAT*MINOR_R**2
        PSI2(I,K1) = PSI1(I,K)*BRAT*MINOR_R**2
!        GSCURR2(I,K0) =-GSPARACURR(I,K)*BRAT/MINOR_R
!        GSCURR2(I,K1) = GSPARACURR(I,K)*BRAT/MINOR_R
!        GSCURphi2(I,K0) = GSPARACURphi(I,K)*BRAT/MINOR_R
!        GSCURphi2(I,K1) = GSPARACURphi(I,K)*BRAT/MINOR_R
!        GSCURZ2(I,K0) = GSPARACURZ(I,K)*BRAT/MINOR_R
!        GSCURZ2(I,K1) = GSPARACURZ(I,K)*BRAT/MINOR_R
      end do
      end do

!      pphi_max_a = pphi_max_a*BRAT*MINOR_R**2
!      pphi_min_a = pphi_min_a*BRAT*MINOR_R**2
      psimax = psimax*BRAT*MINOR_R**2

      pphi_max_a = e_a*psimax
      pphi_min_a = 0.0d0

!2013-07-24s
      PMAX = PRSR2(IAXIS,JAXIS)
      PMIN = 1.0d-1*PMAX

! search beta_a_max
      beta_a_max = 0.0d0
      do i = 1, nr
        beta_a_max = max(beta_a_max, (prsr2h_para(i,nz) + 2.0d0*prsr2h_perp(i,nz) )*2.0d0/3.0d0 )
      end do
!2013-07-24e

!$omp parallel do
      do j = 1, nz2
      do i = 1, nr
        PRSR2H_PARA(i,j) = PRSR2H_PARA(i,j)*beta_a0/beta_a_max
        PRSR2H_PERP(i,j) = PRSR2H_PERP(i,j)*beta_a0/beta_a_max
      end do
      end do

!$omp parallel do private(i1,j1)
      do k = 1, lphi
      do j = 1, lz
      do i = 1, lr
        i1 = i + kr_offset(my_rank)
        j1 = j + kz_offset(my_rank)
        br(i,j,k) = GSBR2(i1,j1)
        bphi(i,j,k) = GSBphi2(i1,j1)
        bz(i,j,k) = GSBZ2(i1,j1)
        psi(i,j,k) = psi2(i1,j1)
        prs(i,j,k) = max(PMIN, PRSR2(i1,j1) )
        rho(i,j,k) = 1.0d0

        ppara_a(i,j,k) = PRSR2H_PARA(i1,j1)
        pperp_a(i,j,k) = PRSR2H_PERP(i1,j1)
        dns_a(i,j,k)   = (ppara_a(i,j,k) + 2.0d0*pperp_a(i,j,k) ) &
                        /(3.0d0*temp_a)
        mom_a(i,j,k) = 0.0d0

#ifdef KTI
        ppara_i(i,j,k) = 0.50d0*prs(i,j,k)
        pperp_i(i,j,k) = 0.50d0*prs(i,j,k)
        dns_i(i,j,k) = rho(i,j,k)/m_i
        mom_i(i,j,k) = 0.0d0

        ppara_e(i,j,k) = ppara_i(i,j,k)
        pperp_e(i,j,k) = ppara_i(i,j,k)
        dns_e(i,j,k) =-e_i*dns_i(i,j,k)/e_e
        mom_e(i,j,k) = 0.0d0

! prs is electron pressure
        prs(i,j,k) =(ppara_e(i,j,k) + 2.0d0*pperp_e(i,j,k) )/3.0d0
#else
        ppara_i(i,j,k) = 0.0d0
        pperp_i(i,j,k) = 0.0d0
        dns_i(i,j,k) = 0.0d0
        mom_i(i,j,k) = 0.0d0

        ppara_e(i,j,k) = 0.50d0*prs(i,j,k)
        pperp_e(i,j,k) = 0.50d0*prs(i,j,k)
        dns_e(i,j,k) = rho(i,j,k)/m_i
        mom_e(i,j,k) = 0.0d0
!       prs(i,j,k) is fluid pressure
#endif
 
        vcrit(i,j,k) = vcrit0 !vcrit is not used for type_a=0
        nusd(i,j,k) = 0.0d0
        ti(i,j,k)   = 0.0d0
        te(i,j,k)   = 0.0d0
        utor(i,j,k) = 0.0d0

!       vcrit(i,j,k) = vcritsim(i1,j1)
!       nusd(i,j,k) = nusdsim(i1,j1)*sd_accl !2013-10-18
!       ti(i,j,k) = tisim(i1,j1)
!       te(i,j,k) = tesim(i1,j1)
!       utor(i,j,k)=rotsim(i1,j1) !2013-12-18
      end do
      end do
      end do
#endif      

      call gradient(1,psi,psi_r,psi_z,psi_phi)

!$omp parallel do
      do k = 1, lphi
      do j = 1, lz
      do i = 1, lr
        bphi0(i,j,k) = bphi(i,j,k)
        br0(i,j,k) = br(i,j,k)
        bz0(i,j,k) = bz(i,j,k)
        babs0(i,j,k) = sqrt(br0(i,j,k)**2 + bz0(i,j,k)**2 &
                           +bphi0(i,j,k)**2 )
        babs(i,j,k) = babs0(i,j,k)
        prs0(i,j,k) = prs(i,j,k)
        rho0(i,j,k) = rho(i,j,k)

        ppara_e0(i,j,k) = ppara_e(i,j,k)
        pperp_e0(i,j,k) = pperp_e(i,j,k)
        ppara_i0(i,j,k) = ppara_i(i,j,k)
        pperp_i0(i,j,k) = pperp_i(i,j,k)
        ppara_a0(i,j,k) = ppara_a(i,j,k)
        pperp_a0(i,j,k) = pperp_a(i,j,k)

        qpara_i(i,j,k) = 0.0d0
        qperp_i(i,j,k) = 0.0d0
        qpara_i0(i,j,k) = 0.0d0
        qperp_i0(i,j,k) = 0.0d0

        qpara_a(i,j,k) = 0.0d0
        qperp_a(i,j,k) = 0.0d0
        qpara_a0(i,j,k) = 0.0d0
        qperp_a0(i,j,k) = 0.0d0

        dns_e0(i,j,k) = max(dns_e_min, dns_e(i,j,k) )
        dns_i0(i,j,k) = max(dns_i_min, dns_i(i,j,k) )
        dns_a0(i,j,k) = max(dns_a_min, dns_a(i,j,k) )

        mom_e0(i,j,k) = mom_e(i,j,k)
        mom_i0(i,j,k) = mom_i(i,j,k)
        mom_a0(i,j,k) = mom_a(i,j,k)

        temp_e0(i,j,k) = (ppara_e0(i,j,k) + 2.0d0*pperp_e0(i,j,k) ) &
                        /(3.0d0*dns_e0(i,j,k) )
        temp_i0(i,j,k) = (ppara_i0(i,j,k) + 2.0d0*pperp_i0(i,j,k) ) &
                        /(3.0d0*dns_i0(i,j,k) )
!        temp_a0(i,j,k) = (ppara_a0(i,j,k) + 2.0d0*pperp_a0(i,j,k) ) &
!                        /(3.0d0*dns_a0(i,j,k) )
        temp_a0(i,j,k) = temp_a !2015-09-28

        ephi(i,j,k) = 0.0d0
        er(i,j,k) = 0.0d0
        ez(i,j,k) = 0.0d0
        epara(i,j,k) = 0.0d0

        vphi(i,j,k) = 0.0d0
        vr(i,j,k) = 0.0d0
        vz(i,j,k) = 0.0d0

        br_ref(i,j,k)   = br(i,j,k)
        bz_ref(i,j,k)   = bz(i,j,k)
        bphi_ref(i,j,k) = bphi(i,j,k)

        vr_ref(i,j,k)   = vr(i,j,k)
        vz_ref(i,j,k)   = vz(i,j,k)
        vphi_ref(i,j,k) = vphi(i,j,k)

        rho_ref(i,j,k)  = rho(i,j,k)
        prs_ref(i,j,k)  = prs(i,j,k)
      end do
      end do
      end do

      bmin = sqrt(br(1,1,1)**2 + bz(1,1,1)**2 + bphi(1,1,1)**2)

!$omp parallel do reduction(min:bmin) private(b2)
      do k = 1, lphi
      do j = 1, lz
      do i = 1, lr
        b2 = br(i,j,k)**2 + bz(i,j,k)**2 + bphi(i,j,k)**2
        bmin = min(bmin, sqrt(b2) )
#ifndef KTI                             
        enrg(i,j,k) = 0.50d0*(vr(i,j,k)**2 + vz(i,j,k)**2 & 
                             +vphi(i,j,k)**2)/rho(i,j,k) &
                   + 0.50d0*b2 + prs(i,j,k)/(gamma - 1.0d0)
        enrg0(i,j,k) = enrg(i,j,k)
        enrg_ref(i,j,k) = enrg(i,j,k)
#endif
      end do
      end do
      end do

!2012-06-17
      call mpi_allreduce(bmin,bmin_global,1,mpi_real8 ,mpi_min&
           &,mhd_world,mpi_err)
!2012-06-17 end

      bmin = bmin_global


! initial electron momentum

!      call rotation(1,br,bz,bphi,cr,cz,cphi)
!      mom_e0 =(cr*br + cz*bz + cphi*bphi)/babs/e_e*m_e


! uniform temperature
!!$omp parallel do
!      do k = 1, lphi
!      do j = 1, lz
!      do i = 1, lr
!        temp_i0(i,j,k)  = temp_i
!        ppara_i0(i,j,k) = dns_i0(i,j,k)*temp_i
!        pperp_i0(i,j,k) = dns_i0(i,j,k)*temp_i

!        ppara_i(i,j,k) = ppara_i0(i,j,k)
!        pperp_i(i,j,k) = pperp_i0(i,j,k)
!      end do
!      end do
!      end do
!2020-04-20e

!      call gradient(1,dns_e0,dns_e0_r,dns_e0_z,dns_e0_phi)
      call gradient(1,dns_i0,dns_i0_r,dns_i0_z,dns_i0_phi)
      call gradient(1,dns_a0,dns_a0_r,dns_a0_z,dns_a0_phi)

!      call gradient(1,temp_e0,temp_e0_r,temp_e0_z,temp_e0_phi)
      call gradient(1,temp_i0,temp_i0_r,temp_i0_z,temp_i0_phi)
      call gradient(1,temp_a0,temp_a0_r,temp_a0_z,temp_a0_phi)

!      call gradient(1,ppara_e0,ppara_e0_r,ppara_e0_z,ppara_e0_phi)
!      call gradient(1,ppara_i0,ppara_i0_r,ppara_i0_z,ppara_i0_phi)
!      call gradient(1,ppara_a0,ppara_a0_r,ppara_a0_z,ppara_a0_phi)

! vaccuum

!$omp parallel do private(sqrtpsi)
      do k = 1, lphi
      do j = 1, lz
      do i = 1, lr
        eta_spt(i,j,k)= eta

        if(psi(i,j,k)/psimax.gt.psivac)then
          vac(i,j,k) = 1.0d0
          nu_spt(i,j,k) = nu
        else
          vac(i,j,k) = 0.0d0
          nu_spt(i,j,k) = nu_vac
          utor(i,j,k) = 0.0d0
        end if

        sqrtpsi = sqrt(max(0.0d0, min(1.0d0, 1.0d0 - psi(i,j,k)/psimax)))
        if(sqrtpsi.gt.sqrtpsi_chi_enh)then
          chi_spt(i,j,k) = chi_enhanced
        else
          chi_spt(i,j,k) = chi
        end if

        if(sqrtpsi.gt.sqrtpsi_chi_enh)then
          nun_spt(i,j,k) = nun_enhanced
        else
          nun_spt(i,j,k) = nu_n
        end if

!        chi_spt(i,j,k) = chi + (chi_enhanced - chi) &
!           *(0.50d0 + 0.50d0*erf((sqrtpsi-sqrtpsi_chi_enh)/dsqrtpsi_chi_enh) ) 
!        chi_spt(i,j,k) = chi + (chi_enhanced - chi) &
!           *exp(-((sqrtpsi-sqrtpsi_chi_enh)/dsqrtpsi_chi_enh)**2) 

!for pressure gradient term in Ohm's law        
        if(sqrtpsi.gt.sqrtpsi_pe)then
          vac_pe(i,j,k) = 0.0d0
        else
          vac_pe(i,j,k) = 1.0d0
        end if

      end do
      end do
      end do

      if(my_rank_r.eq.0)then 
!$omp parallel do
        do k = 1, lphi
        do j = 1, lz
          eta_spt(1,j,k)= 0.0d0
        end do
        end do
      end if

      if(my_rank_r.eq.mpi_proc_r-1)then 
!$omp parallel do
        do k = 1, lphi
        do j = 1, lz
          eta_spt(lr,j,k)= 0.0d0
        end do
        end do
      end if

      if(my_rank_z.eq.0)then 
!$omp parallel do
        do k = 1, lphi
        do i = 1, lr
          eta_spt(i,1,k)= 0.0d0
        end do
        end do
      end if

      if(my_rank_z.eq.mpi_proc_z-1)then 
!$omp parallel do
        do k = 1, lphi
        do i = 1, lr
          eta_spt(i,lz,k)= 0.0d0
        end do
        end do
      end if


! total energy for normalization

      energy_e = 0.0d0
      energy_i = 0.0d0
      energy_a = 0.0d0

!$omp parallel do reduction(+:energy_e,energy_i,energy_a)
      do k = lphistart, lphiend
      do j = lzstart, lzend
      do i = lrstart, lrend
        energy_e = energy_e + 0.50d0*(ppara_e(i,j,k) + 2.0d0&
             &*pperp_e(i,j,k) ) *grr(i,j,k)*dr*dz*dphi*vac(i,j,k)
        energy_i = energy_i + 0.50d0*(ppara_i(i,j,k) + 2.0d0&
             &*pperp_i(i,j,k) ) *grr(i,j,k)*dr*dz*dphi*vac(i,j,k)
        energy_a = energy_a + 0.50d0*(ppara_a(i,j,k) + 2.0d0&
             &*pperp_a(i,j,k) ) *grr(i,j,k)*dr*dz*dphi*vac(i,j,k)
      end do
      end do
      end do

      call mpi_allreduce(energy_e,total_energy_e,1,mpi_real8 ,mpi_sum&
           &,mhd_world,mpi_err)
      call mpi_allreduce(energy_i,total_energy_i,1,mpi_real8 ,mpi_sum&
           &,mhd_world,mpi_err)
      call mpi_allreduce(energy_a,total_energy_a,1,mpi_real8 ,mpi_sum&
           &,mhd_world,mpi_err)

! for MHD NL study

      if(my_rank.eq.0)then
        write(7,*)'beta_a0=',beta_a0
        write(7,*)'total_energy_e,i,a=' &
                  ,total_energy_e,total_energy_i,total_energy_a
      end if


! injected energy, 2013-05-03
      co_stored = 0.0d0
      cntr_stored = 0.0d0

! initial value of transferred energy from field to particle
      e_trans = 0.0d0
      e_dissp = 0.0d0

! magnetic field gradient and curvature

      call gradmg

#ifdef KTI
          call e_field_kti(0)
#else
          call e_field(0)
#endif

! for subroutine push
! copy invariable arrays

!$omp parallel
!$omp workshare
      fld(1,:,:,:) = er(:,:,:)
      fld(2,:,:,:) = ez(:,:,:)
      fld(3,:,:,:) = ephi(:,:,:)
      fld(4,:,:,:) = epara(:,:,:)
      fld(5,:,:,:) = br(:,:,:) - br0(:,:,:)
      fld(6,:,:,:) = bz(:,:,:) - bz0(:,:,:)
      fld(7,:,:,:) = bphi(:,:,:) - bphi0(:,:,:)
      fld(8,:,:,:) = br0(:,:,:)
      fld(9,:,:,:) = bz0(:,:,:)
      fld(10,:,:,:)= bphi0(:,:,:)
      fld(11,:,:,:)= gradbr(:,:,:)
      fld(12,:,:,:)= gradbz(:,:,:)
      fld(13,:,:,:)= gradbphi(:,:,:)
      fld(14,:,:,:)= curvbr(:,:,:)
      fld(15,:,:,:)= curvbz(:,:,:)
      fld(16,:,:,:)= curvbphi(:,:,:)
      fld(17,:,:,:)= psi(:,:,:)
      fld(18,:,:,:)= psi_r(:,:,:)
      fld(19,:,:,:)= psi_z(:,:,:)
      fld(20,:,:,:)= psi_phi(:,:,:)
      fld(21,:,:,:)= vcrit(:,:,:)
      fld(22,:,:,:)= nusd(:,:,:)
!grad B0 vector
      fld(31,:,:,:)= gradbr(:,:,:)
      fld(32,:,:,:)= gradbz(:,:,:)
      fld(33,:,:,:)= gradbphi(:,:,:)
!$omp end workshare
!$omp end parallel

end
!--------------------------------------------------------------------
subroutine beam_deposit
! read beam deposition profile on DIII-D; 2013-05-02
!--------------------------------------------------------------------
      use mpiset
      use field
      use particle
      use grid
      use random_num
      implicit none

      character(len=5)::char
      integer::npr,npr5,n,nrandom,i
      real(8)::e_total,factor

      real(8)::rand_no(marker_each) !2016-01-09
      real(8), allocatable::rand_bd(:)
      real(4)::bottom_sim=-1.01454 !in [m]


!      open(61,file='/mlng/todo/142111M34_combined_birth_t500_530ms.out' &
!             ,form='formatted')

      if(my_rank.eq.0)then
        read(61,*)char,npr
        write(7,*)'npr=',npr
!        if(npr.gt.marker_a0_total/2)then !2013-05-26
!          npr = marker_a0_total/2 !2013-05-26
        if(npr.gt.marker_a0*mpi_proc/2)then !2013-05-26
          npr = marker_a0*mpi_proc/2 !2013-05-26
          write(7,*)'npr is set to be equal to marker_a0_total/2=',npr
        end if
        write(7,*)'total normalized time to be simulated =',tmax
      end if

      call mpi_bcast(npr,1,mpi_integer,0,mpi_comm_world,mpi_err)

      allocate(bdp_data(5,npr))
      allocate(t_bd(npr))
      allocate(injection_bd(npr))
      allocate(rand_bd(npr))

      if(my_rank.eq.0)then
        read(61,*)char
        read(61,*)char

        do n = 1, npr
          read(61,*)(bdp_data(i,n),i=1,5)
! bdp_data(i,:) contains the following data
!          read(61,*)r(n),z(n),lambda(n),e(n),zeta(n)
! r(n) and z(n) are in [cm], convert to [m], and normalize by rho_a
          bdp_data(1,n) = bdp_data(1,n)*1.0e-2/rho_a
          bdp_data(2,n) =(bdp_data(2,n)*1.0e-2 - bottom_sim)/rho_a
! e(n) is in [eV], convert to v[m/s], and normalize by va
          bdp_data(4,n) = 0.9786e7*sqrt(bdp_data(4,n)*1.0d-6)/va  ! 0.9786e7 m/s for 1MeV
        end do
      end if

! for major radius
      call ransuu(npr,rand_bd)
      bdp_data(5,:) = rand_bd(:)

      npr5 = npr*5
      call mpi_bcast(bdp_data,npr5,mpi_real,0,mpi_comm_world,mpi_err)

      e_total = 0.0d0
      do n = 1, npr
        t_bd(n) = tmax/real(npr)*real(n) !injection time
        injection_bd(n) = 0 !not injected:0, injected:1
        e_total = e_total + bdp_data(4,n)**2
      end do
        e_total = 0.50d0*m_a*e_total/tmax*phimode
        factor = hpower*1.0d6/(bunit**2/(4.0d-7*pi)*rho_a**3*omega_a)/e_total
        n_bd = npr 
        n_inj = 0 !=number of particles already injected

!      if(my_rank.eq.0)then
!        write(7,*)'t_bd(1),t_bd(2)=',t_bd(1),t_bd(2)
!      end if
!--------------------------------------------------------------------
! hereafer, for buffer particles
!--------------------------------------------------------------------

! particle disribution region

      rmin = major_r - minor_r + rleng*dble(my_rank_r)/dble(mpi_proc_r)
      rmax = major_r - minor_r + rleng*dble(my_rank_r+1)/dble(mpi_proc_r)
      zmin = zleng*dble(my_rank_z)/dble(mpi_proc_z)
      zmax = zleng*dble(my_rank_z+1)/dble(mpi_proc_z)
      phimin = phileng*dble(my_rank_phi)/dble(mpi_proc_phi)
      phimax = phileng*dble(my_rank_phi+1)/dble(mpi_proc_phi)

      nrandom = marker_a

! gc_x(:,n) = 1:r, 2:z, 3:phi, 4:p, 5:mu, 6:weight, 7:active, 8:p_phi, 9:fff, 10:fnrml

! for major radius
      call ransuu(nrandom,rand_no)

      do n = 1, nrandom
        gc_a(1,n) = rmin + (rmax - rmin)*rand_no(n)
      end do

! for vertical position
      call ransuu(nrandom,rand_no)

      do n = 1, nrandom
        gc_a(2,n) = zmin + (zmax - zmin)*rand_no(n)
      end do

! for toroidal angle
      call ransuu(nrandom,rand_no)

      do n = 1, nrandom
        gc_a(3,n) = phimin + (phimax - phimin)*rand_no(n)
      end do


      do n = 1, marker_a
        v_a(n) = 1.0d-2
        lambda_a(n) = 1.0d0

        gc_a(4,n) = m_a*v_a(n)*lambda_a(n)
        gc_a(5,n)= 0.0d0
        gc_a(7,n) = 0.0d0
        gc_a(8,n) = 0.0d0
        gc_a(9,n) = 1.0d0
        gc_a(10,n) = factor

        gc_a(6,n) = gc_a(9,n)*gc_a(10,n) 
        node_a(n) = mod(my_rank + 1 + nprocess, nprocess)
      end do


end
!--------------------------------------------------------------------
subroutine flux_coordinates
! create flux coordinate grid points for data analysis 
!--------------------------------------------------------------------
      use mpiset
      use grid
      use flux_c
      use field, only:br0,bz0,bphi0
      implicit none

      integer::i,j,l,m,n,node_grid,i0,i1
      integer::nb,ia,ia1,ja,ja1
      real(8)::theta,phi
      real(8)::r_flux(ltheta,lpsi),z_flux(ltheta,lpsi)
      real(8)::ar,ar1,az,az1
      real(8)::aaa1,aaa2,aaa3,aaa4
      real(8)::value_r,value_z,value_phi


! ************************************************************
!
! create data analysis grid points
!
! ************************************************************


      dtheta = 2.0d0*pi/dble(ltheta)

      do l = 1, lpsi
      do i = 1, ltheta
#ifdef EXPEQ
! crdmagr, z are already set in to-MEGA512.f, 2021-05-02
        r_flux(i,l) = crdmagr(l,i)
        z_flux(i,l) = crdmagz(l,i)
#else         
        r_flux(i,l) = crdmagr(l,i)*minor_r
! on midplane z=0 in boozer while z=zleng/2 in the simulation
        z_flux(i,l) = crdmagz(l,i)*minor_r + 0.50d0*zleng
#endif
      end do
      end do

      nb = 0

      do l = 1, lpsi
      do i = 1, ltheta
!2012-06-17
        node_grid = max(0, min(mpi_proc_r-1, &
                    int(dble(mpi_proc_r)*(r_flux(i,l)-(major_r-minor_r) )/rleng ) &
                           ) ) &
                  + max(0, min(mpi_proc_z-1, &
                    int(dble(mpi_proc_z)*z_flux(i,l)/zleng ) &
                           ) )*mpi_proc_r &
                  + my_rank_phi*mpi_proc_pol &
                  + my_rank_ptcl*mpi_proc_mhd
!2012-06-17 end

        if(node_grid.eq.my_rank)then
          nb = min(nb + 1,n_flux_size)
          r_flux_grid(nb) = r_flux(i,l)
          z_flux_grid(nb) = z_flux(i,l)
          i_flux_grid(nb) = i
          l_flux_grid(nb) = l
        end if

      end do
      end do


      if(nb.eq.n_flux_size)then
        write(7,*)'overflow in n_flux_size for rank ',my_rank
      end if

      n_flux_grid = nb

      do nb = 1, n_flux_grid
! i: theta direction, l: radial direction
        i = i_flux_grid(nb)
        l = l_flux_grid(nb)
        i0 = max(1,i-1)
        i1 = min(ltheta,i+1)

        theta = dtheta*dble(i_flux_grid(nb)-1)

        ia = max(1, min(lr-1, &
                 int( ( r_flux_grid(nb)-(major_r - minor_r) )/dr) + 1 &
                      - kr_offset(my_rank) ) )
        ia1= ia + 1
        ja = max(1, min(lz-1, &
                 int( z_flux_grid(nb)/dz ) + 1 - kz_offset(my_rank) ) )
        ja1= ja + 1

        ar1 = (r_flux_grid(nb)-grr(ia,ja,1) ) /dr
        ar  = 1.0d0 - ar1
        az1 = (z_flux_grid(nb)-gzz(ia,ja,1) ) /dz
        az  = 1.0d0 - az1

        aaa1 = ar *az
        aaa2 = ar1*az
        aaa3 = ar *az1
        aaa4 = ar1*az1

! interpolated value

        value_r = br0(ia, ja, 1 )*aaa1  + br0(ia1,ja, 1 )*aaa2 &
                + br0(ia, ja1,1 )*aaa3  + br0(ia1,ja1,1 )*aaa4
        value_z = bz0(ia, ja, 1 )*aaa1  + bz0(ia1,ja, 1 )*aaa2 &
                + bz0(ia, ja1,1 )*aaa3  + bz0(ia1,ja1,1 )*aaa4
        value_phi = bphi0(ia, ja, 1 )*aaa1  + bphi0(ia1,ja, 1 )*aaa2 &
                  + bphi0(ia, ja1,1 )*aaa3  + bphi0(ia1,ja1,1 )*aaa4

        value_r = value_r*sign(1.0d0,value_phi)
        value_z = value_z*sign(1.0d0,value_phi)
        value_phi = value_phi*sign(1.0d0,value_phi)

        dt_r(nb) = r_flux(i1,l) - r_flux(i0,l)
        dt_phi(nb) = 0.0d0
        dt_z(nb) = z_flux(i1,l) - z_flux(i0,l)

        dr_r(nb) = dt_z(nb)*value_phi
        dr_phi(nb) =-dt_z(nb)*value_r + dt_r(nb)*value_z
        dr_z(nb) =-dt_r(nb)*value_phi

        dp_r(nb) = 0.0d0
        dp_phi(nb) = 1.0d0
        dp_z(nb) = 0.0d0

        dt_abs(nb) = max(1.0d-10, sqrt(dt_r(nb)**2 + dt_z(nb)**2 + dt_phi(nb)**2) )
        dr_abs(nb) = max(1.0d-10, sqrt(dr_r(nb)**2 + dr_z(nb)**2 + dr_phi(nb)**2) )

        dt_r(nb)   = dt_r(nb)/dt_abs(nb)
        dt_z(nb)   = dt_z(nb)/dt_abs(nb)
        dt_phi(nb) = dt_phi(nb)/dt_abs(nb)

        dr_r(nb)   = dr_r(nb)/dr_abs(nb)
        dr_z(nb)   = dr_z(nb)/dr_abs(nb)
        dr_phi(nb) = dr_phi(nb)/dr_abs(nb)
      end do

!2012-07-11
!$omp parallel do private(theta)
      do i = 1, ltheta
        theta= dtheta*dble(i-1)
        do m = 0, mpol
          cos_theta_flux(m,i) = cos(dble(m)*theta)
          sin_theta_flux(m,i) = sin(dble(m)*theta)
        end do
      end do

!$omp parallel do private(phi)
      do j = 1, lphi
        phi = gphi(1,1,j)
        do n =-ntor, ntor
          cos_phi_flux(n,j) = cos(dble(n)*phi*phimode)
          sin_phi_flux(n,j) = sin(dble(n)*phi*phimode)
        end do
      end do
!2012-07-11 end

end
!--------------------------------------------------------------------
subroutine harmonics
! analysis in the flux coordinates
! 2019-06-04
! dns_i,mom_i,ppara_i,pperp_i,dns_a,mom_a,ppara_a,pperp_a,er,ez,ephi are analyzed
!--------------------------------------------------------------------
      use mpiset
      use grid
      use flux_c
      use equi_sol
      use field, only:br,bz,bphi,ur,uz,uphi,prs,rho &
                     ,er,ez,ephi &
                     ,dns_i,mom_i,ppara_i,pperp_i &
                     ,dns_a,mom_a,ppara_a,pperp_a
      implicit none
!      integer, parameter::lfluid=17
      integer, parameter::lfluid=19 !2019-06-04
      real(8)::brad_harmonics(0:mpol,-ntor:ntor,lpsi,2)
      real(8)::btheta_harmonics(0:mpol,-ntor:ntor,lpsi,2)
      real(8)::bphi_harmonics(0:mpol,-ntor:ntor,lpsi,2)
      real(8)::vrad_harmonics(0:mpol,-ntor:ntor,lpsi,2)
      real(8)::vtheta_harmonics(0:mpol,-ntor:ntor,lpsi,2)
      real(8)::vphi_harmonics(0:mpol,-ntor:ntor,lpsi,2)
      real(8)::erad_harmonics(0:mpol,-ntor:ntor,lpsi,2)
      real(8)::etheta_harmonics(0:mpol,-ntor:ntor,lpsi,2)
      real(8)::ephi_harmonics(0:mpol,-ntor:ntor,lpsi,2)
      real(8)::prs_harmonics(0:mpol,-ntor:ntor,lpsi,2)
      real(8)::rho_harmonics(0:mpol,-ntor:ntor,lpsi,2)
      real(8)::dns_i_harmonics(0:mpol,-ntor:ntor,lpsi,2)
      real(8)::mom_i_harmonics(0:mpol,-ntor:ntor,lpsi,2)
      real(8)::ppara_i_harmonics(0:mpol,-ntor:ntor,lpsi,2)
      real(8)::pperp_i_harmonics(0:mpol,-ntor:ntor,lpsi,2)
      real(8)::dns_a_harmonics(0:mpol,-ntor:ntor,lpsi,2)
      real(8)::mom_a_harmonics(0:mpol,-ntor:ntor,lpsi,2)
      real(8)::ppara_a_harmonics(0:mpol,-ntor:ntor,lpsi,2)
      real(8)::pperp_a_harmonics(0:mpol,-ntor:ntor,lpsi,2)
!      real(8)::qpara_a_harmonics(0:mpol,-ntor:ntor,lpsi,2)
!      real(8)::qperp_a_harmonics(0:mpol,-ntor:ntor,lpsi,2)

      real(8)::fluid_harmonics(0:mpol,-ntor:ntor,lpsi,2,lfluid)
      real(8)::fluid_total(0:mpol,-ntor:ntor,lpsi,2,lfluid)
      integer::lsize
#ifdef AMDGPU
      integer::i,j,k,l
#endif

! vr, vz, vphi are momentum 2012-03-26
      call analyz_v(ur,uz,uphi,vrad_harmonics,vtheta_harmonics,vphi_harmonics)
! end 2012-03-26

      call analyz_v(br,bz,bphi,brad_harmonics,btheta_harmonics,bphi_harmonics)
      call analyz_v(er,ez,ephi,erad_harmonics,etheta_harmonics,ephi_harmonics)
      call analyz_s(prs,prs_harmonics)
      call analyz_s(rho,rho_harmonics)
      call analyz_s(dns_i,dns_i_harmonics)
      call analyz_s(mom_i,mom_i_harmonics)
      call analyz_s(ppara_i,ppara_i_harmonics)
      call analyz_s(pperp_i,pperp_i_harmonics)
      call analyz_s(dns_a,dns_a_harmonics)
      call analyz_s(mom_a,mom_a_harmonics)
      call analyz_s(ppara_a,ppara_a_harmonics)
      call analyz_s(pperp_a,pperp_a_harmonics)
!      call analyz_s(qpara_a,qpara_a_harmonics)
!      call analyz_s(qperp_a,qperp_a_harmonics)

#ifndef AMDGPU
      fluid_harmonics(:,:,:,:,1) = vrad_harmonics(:,:,:,:)
      fluid_harmonics(:,:,:,:,2) = vtheta_harmonics(:,:,:,:)
      fluid_harmonics(:,:,:,:,3) = vphi_harmonics(:,:,:,:)
      fluid_harmonics(:,:,:,:,4) = brad_harmonics(:,:,:,:)
      fluid_harmonics(:,:,:,:,5) = btheta_harmonics(:,:,:,:)
      fluid_harmonics(:,:,:,:,6) = bphi_harmonics(:,:,:,:)
      fluid_harmonics(:,:,:,:,7) = erad_harmonics(:,:,:,:)
      fluid_harmonics(:,:,:,:,8) = etheta_harmonics(:,:,:,:)
      fluid_harmonics(:,:,:,:,9) = ephi_harmonics(:,:,:,:)
      fluid_harmonics(:,:,:,:,10)= prs_harmonics(:,:,:,:)
      fluid_harmonics(:,:,:,:,11)= rho_harmonics(:,:,:,:)
      fluid_harmonics(:,:,:,:,12)= dns_i_harmonics(:,:,:,:)
      fluid_harmonics(:,:,:,:,13)= mom_i_harmonics(:,:,:,:)
      fluid_harmonics(:,:,:,:,14)= ppara_i_harmonics(:,:,:,:)
      fluid_harmonics(:,:,:,:,15)= pperp_i_harmonics(:,:,:,:)
      fluid_harmonics(:,:,:,:,16)= dns_a_harmonics(:,:,:,:)
      fluid_harmonics(:,:,:,:,17)= mom_a_harmonics(:,:,:,:)
      fluid_harmonics(:,:,:,:,18)= ppara_a_harmonics(:,:,:,:)
      fluid_harmonics(:,:,:,:,19)= pperp_a_harmonics(:,:,:,:)
#else
!$omp target teams distribute parallel do collapse(4) private(i,j,k,l)
      do l=1,2
      do k=1,lpsi
      do j=-ntor,ntor
      do i=0,mpol
      fluid_harmonics(i,j,k,l,1) = vrad_harmonics(i,j,k,l)
      fluid_harmonics(i,j,k,l,2) = vtheta_harmonics(i,j,k,l)
      fluid_harmonics(i,j,k,l,3) = vphi_harmonics(i,j,k,l)
      fluid_harmonics(i,j,k,l,4) = brad_harmonics(i,j,k,l)
      fluid_harmonics(i,j,k,l,5) = btheta_harmonics(i,j,k,l)
      fluid_harmonics(i,j,k,l,6) = bphi_harmonics(i,j,k,l)
      fluid_harmonics(i,j,k,l,7) = erad_harmonics(i,j,k,l)
      fluid_harmonics(i,j,k,l,8) = etheta_harmonics(i,j,k,l)
      fluid_harmonics(i,j,k,l,9) = ephi_harmonics(i,j,k,l)
      fluid_harmonics(i,j,k,l,10)= prs_harmonics(i,j,k,l)
      fluid_harmonics(i,j,k,l,11)= rho_harmonics(i,j,k,l)
      fluid_harmonics(i,j,k,l,12)= dns_i_harmonics(i,j,k,l)
      fluid_harmonics(i,j,k,l,13)= mom_i_harmonics(i,j,k,l)
      fluid_harmonics(i,j,k,l,14)= ppara_i_harmonics(i,j,k,l)
      fluid_harmonics(i,j,k,l,15)= pperp_i_harmonics(i,j,k,l)
      fluid_harmonics(i,j,k,l,16)= dns_a_harmonics(i,j,k,l)
      fluid_harmonics(i,j,k,l,17)= mom_a_harmonics(i,j,k,l)
      fluid_harmonics(i,j,k,l,18)= ppara_a_harmonics(i,j,k,l)
      fluid_harmonics(i,j,k,l,19)= pperp_a_harmonics(i,j,k,l)
      end do
      end do
      end do
      end do
#endif
!      fluid_harmonics(:,:,:,:,16)= qpara_a_harmonics(:,:,:,:)
!      fluid_harmonics(:,:,:,:,17)= qperp_a_harmonics(:,:,:,:)

! mpi summation

      lsize = (mpol+1)*(2*ntor+1)*lpsi*2*lfluid

      call mpi_reduce(fluid_harmonics(0,-ntor,1,1,1) &
                     ,fluid_total(0,-ntor,1,1,1),lsize,mpi_real8 &
                     ,mpi_sum,0,mhd_world,mpi_err)

      if(my_rank.eq.0)then
        write(31)kstep,t,rpsi,gpsi_nrm,qpsi,fluid_total

        flush(31)
      end if

end
!--------------------------------------------------------------------
subroutine gradmg
! magnetic field gradient and curvature
!--------------------------------------------------------------------
      use parameters
      use field, only:br,bz,bphi,babs,gradbr,gradbz,gradbphi &
                     ,curvbr,curvbz,curvbphi
      use grid
      implicit none
      integer::i
      real(8)::hmr(lr,lz,lphi),hmz(lr,lz,lphi),hmphi(lr,lz,lphi)
!      real(8)::dr1,dz1,dphi1,rg1

!      dr1= 0.50d0/dr
!      dphi1= 0.50d0/dphi
!      dz1= 0.50d0/dz

#ifdef AMDGPU
!$omp target teams distribute parallel do private(i)
#else        
!$omp parallel do
#endif
!      do k = 1, lphi
!      do j = 1, lz
      do i = 1, lrzphi
        babs(i,1,1) = max(eps_b,sqrt(br(i,1,1)**2 + bz(i,1,1)**2 + bphi(i,1,1)**2))
        hmr(i,1,1) = br(i,1,1)/babs(i,1,1)
        hmz(i,1,1) = bz(i,1,1)/babs(i,1,1)
        hmphi(i,1,1) = bphi(i,1,1)/babs(i,1,1)
      end do
!      end do
!      end do

      call gradient(1,babs,gradbr,gradbz,gradbphi)
      call rotation(1,hmr,hmz,hmphi,curvbr,curvbz,curvbphi)

end
!--------------------------------------------------------------------
subroutine initial_particle(marker_num0,marker_num,m,e &
        ,type,temp,vbirth,deltav,clambda0,dclambda,vmin,vmax & !2012-06-17
        ,gc,v,lambda,node &
        ,cf_pphi,pphi_min,pphi_max,total_energy,ispecies)
! type=0: maxwellian, type=1: slowing down, type=2: beam
! type=3: beam with finite pitch angle width, 2012-06-17
! ispecies=0: electron, 1: ion, 2: alpha
! modified on 2015-06-23
! gc_x(:,n) = 1:r, 2:z, 3:phi, 4:p, 5:mu, 6:weight, 7:active, 8:p_phi, 9:fff, 10:fnrml
!--------------------------------------------------------------------
      use mpiset
      use field
      use grid
      use particle, only:rmin,rmax,zmin,zmax,phimin,phimax
      use equi_sol, only:raxis
      use random_num
      implicit none

      integer::node_up,node_down
      integer::marker_num0,marker_num,type
      integer::ispecies
      integer::nrandom
      real(8)::m,e
      real(8)::gc(ngc2,marker_num)
      real(8)::v(marker_num),lambda(marker_num)
      integer::node(marker_num)
      real(8)::cf_pphi(0:lcfpphi),pphi_min,pphi_max,total_energy !2016-02-04
      real(8)::dpartr,vmin,vmax,dv
      real(8)::temp,vbirth,vcritp,deltav
      real(8)::clambda,pitchfact,clambda0,dclambda !2012-06-17
      integer::n,ia,ia1,ja,ja1,ka,ka1,n_buf
      real(8)::ar,ar1,az,az1,aphi,aphi1
      real(8)::aaa1,aaa2,aaa3,aaa4,aaa5,aaa6,aaa7,aaa8
      real(8)::bre,bze,bphie,babse,psip,psinrm,sptfact
      real(8)::total_e,energy_p
      real(8)::total_n,temp_p,number_p
      real(8)::r_p,z_p,phi_p,mom_p,total_r,total_z,total_phi,total_mom
!2017-09-13s
      real(8)::enpara_p,total_epara,ptotal
      integer::i,j,k
!2017-09-13e
      real(8)::vnrm,fraction
!      real(8)::rmin,rmax,zmin,zmax,phimin,phimax
      real(8)::pphi_n,rminor,bmax,ekin,vpara
      real(8)::factor,temp_local !2016-08-05
      real(8)::rand_no(marker_each)
      real(8)::ftotal

      node_up = mod(my_rank + 1 + nprocess, nprocess)
      node_down = mod(my_rank - 1 + nprocess, nprocess)

! particle disribution region

      rmin = major_r - minor_r + rleng*dble(my_rank_r)/dble(mpi_proc_r)
      rmax = major_r - minor_r + rleng*dble(my_rank_r+1)/dble(mpi_proc_r)
      zmin = zleng*dble(my_rank_z)/dble(mpi_proc_z)
      zmax = zleng*dble(my_rank_z+1)/dble(mpi_proc_z)
      phimin = phileng*dble(my_rank_phi)/dble(mpi_proc_phi)
      phimax = phileng*dble(my_rank_phi+1)/dble(mpi_proc_phi)

! for paticle loading: uniform in phase space

      nrandom = marker_num0
!      nrandom_total = marker_num0*mpi_proc

! for absolute velocity

      dv = (vmax**3-vmin**3)/dble(marker_num0)


      do n = 1, marker_num0
        v(n) = (vmin**3 + dv*(dble(n) - 0.50d0) )**(1.0d0/3.0d0)
      end do


! for major radius

      call ransuu(nrandom,rand_no)

!      dpartr = rmax**2 - rmin**2
      dpartr = rmax - rmin


      do n = 1, marker_num0
        gc(1,n) = rmin + rand_no(n)*dpartr
!        gc(1,n) = sqrt(rmin**2 + rand_no(n)*dpartr)
      end do

! for vertical position

      call ransuu(nrandom,rand_no)


      do n = 1, marker_num0
        gc(2,n) = zmin + (zmax - zmin)*rand_no(n)
      end do

! for toroidal angle

      call ransuu(nrandom,rand_no)


      do n = 1, marker_num0
        gc(3,n) = phimin + (phimax - phimin)*rand_no(n)
      end do

! for pitch angle

      call ransuu(nrandom,rand_no)

      if(type.eq.0.or.type.eq.1.or.type.eq.3)then
!      if(type.eq.0.or.type.eq.1)then !2025-05-01 for AUG

        do n = 1, marker_num0
          lambda(n) =-1.0d0 + 2.0d0*rand_no(n)
        end do
      else if(type.eq.2)then

        do n = 1, marker_num0
          lambda(n) = sign(1.0d0, rand_no(n) - 0.50d0)
        end do

!2025-05-01s
! for AUG, co-going particles with anisotropic distribution
! psi>0, B_phi<0, co-going particles are counter-going to magnetic field
!                 lambda < 0
!      else if(type.eq.3)then

!        do n = 1, marker_num0
!          lambda(n) =-1.0d0 + rand_no(n)
!        end do
!2025-05-01e
      end if

      do n = 1, marker_num0
        gc(7,n) = 1.0d0
        node(n) = my_rank
      end do



!ts!$omp parallel private(ia,ia1,ja,ja1,ka,ka1,ar1,ar,az1,az,aphi1,aphi &
!ts!$omp& ,aaa1,aaa2,aaa3,aaa4,aaa5,aaa6,aaa7,aaa8,bre,bze,bphie,psip,vcritp &
!ts!$omp& ,rminor,bmax,vpara,pphi_n,sptfact,babse,ekin,clambda,temp_local)
!ts!$omp do
      do n = 1, marker_num0
        ia = max(1, min(lr-1, &
             int( ( gc(1,n)-(major_r - minor_r) )/dr) + 1 - kr_offset(my_rank) ))
        ia1= ia + 1
        ja = max(1, min(lz-1, &
             int( gc(2,n)/dz ) + 1 - kz_offset(my_rank) ) )
        ja1= ja + 1
        ka = max(1, min(lphi-1, &
             int( gc(3,n)/dphi ) + 1 + lphi_shd - kphi_offset(my_rank) ) )
!             int( gc(3,n)/dphi ) + 3 - kphi_offset(my_rank) ) )
        ka1= ka + 1

        ar1 = (gc(1,n)-grr(ia,ja,ka) ) /dr
        ar  = 1.0d0 - ar1
        az1 = (gc(2,n)-gzz(ia,ja,ka) ) /dz
        az  = 1.0d0 - az1
        aphi1 = (gc(3,n)-gphi(ia,ja,ka) ) /dphi
        aphi  = 1.0d0 - aphi1

        aaa1 = ar *az *aphi
        aaa2 = ar1*az *aphi
        aaa3 = ar *az1*aphi
        aaa4 = ar1*az1*aphi
        aaa5 = ar *az *aphi1
        aaa6 = ar1*az *aphi1
        aaa7 = ar *az1*aphi1
        aaa8 = ar1*az1*aphi1

        bre = br(ia, ja, ka )*aaa1  + br(ia1,ja, ka )*aaa2 &
            + br(ia, ja1,ka )*aaa3  + br(ia1,ja1,ka )*aaa4 &
            + br(ia, ja, ka1)*aaa5  + br(ia1,ja, ka1)*aaa6 &
            + br(ia, ja1,ka1)*aaa7  + br(ia1,ja1,ka1)*aaa8
        bze = bz(ia, ja, ka )*aaa1  + bz(ia1,ja, ka )*aaa2 &
            + bz(ia, ja1,ka )*aaa3  + bz(ia1,ja1,ka )*aaa4 &
            + bz(ia, ja, ka1)*aaa5  + bz(ia1,ja, ka1)*aaa6 &
            + bz(ia, ja1,ka1)*aaa7  + bz(ia1,ja1,ka1)*aaa8
        bphie = bphi(ia, ja, ka )*aaa1  + bphi(ia1,ja, ka )*aaa2 &
              + bphi(ia, ja1,ka )*aaa3  + bphi(ia1,ja1,ka )*aaa4 &
              + bphi(ia, ja, ka1)*aaa5  + bphi(ia1,ja, ka1)*aaa6 &
              + bphi(ia, ja1,ka1)*aaa7  + bphi(ia1,ja1,ka1)*aaa8
        psip  = psi(ia, ja, ka )*aaa1  + psi(ia1,ja, ka )*aaa2 &
              + psi(ia, ja1,ka )*aaa3  + psi(ia1,ja1,ka )*aaa4 &
              + psi(ia, ja, ka1)*aaa5  + psi(ia1,ja, ka1)*aaa6 &
              + psi(ia, ja1,ka1)*aaa7  + psi(ia1,ja1,ka1)*aaa8
        vcritp  = vcrit(ia, ja, ka )*aaa1  + vcrit(ia1,ja, ka )*aaa2 &
                + vcrit(ia, ja1,ka )*aaa3  + vcrit(ia1,ja1,ka )*aaa4 &
                + vcrit(ia, ja, ka1)*aaa5  + vcrit(ia1,ja, ka1)*aaa6 &
                + vcrit(ia, ja1,ka1)*aaa7  + vcrit(ia1,ja1,ka1)*aaa8

!2016-08-06s
      if(ispecies.eq.0)then
        sptfact = dns_e0(ia, ja, ka )*aaa1  + dns_e0(ia1,ja, ka )*aaa2 &
                + dns_e0(ia, ja1,ka )*aaa3  + dns_e0(ia1,ja1,ka )*aaa4 &
                + dns_e0(ia, ja, ka1)*aaa5  + dns_e0(ia1,ja, ka1)*aaa6 &
                + dns_e0(ia, ja1,ka1)*aaa7  + dns_e0(ia1,ja1,ka1)*aaa8
        temp_local = temp_e0(ia, ja, ka )*aaa1  + temp_e0(ia1,ja, ka )*aaa2 &
                + temp_e0(ia, ja1,ka )*aaa3  + temp_e0(ia1,ja1,ka )*aaa4 &
                + temp_e0(ia, ja, ka1)*aaa5  + temp_e0(ia1,ja, ka1)*aaa6 &
                + temp_e0(ia, ja1,ka1)*aaa7  + temp_e0(ia1,ja1,ka1)*aaa8
        v(n) = v(n)*sqrt(temp_local/temp)

      else if(ispecies.eq.1)then
        sptfact = dns_i0(ia, ja, ka )*aaa1  + dns_i0(ia1,ja, ka )*aaa2 &
                + dns_i0(ia, ja1,ka )*aaa3  + dns_i0(ia1,ja1,ka )*aaa4 &
                + dns_i0(ia, ja, ka1)*aaa5  + dns_i0(ia1,ja, ka1)*aaa6 &
                + dns_i0(ia, ja1,ka1)*aaa7  + dns_i0(ia1,ja1,ka1)*aaa8
        temp_local = temp_i0(ia, ja, ka )*aaa1  + temp_i0(ia1,ja, ka )*aaa2 &
                + temp_i0(ia, ja1,ka )*aaa3  + temp_i0(ia1,ja1,ka )*aaa4 &
                + temp_i0(ia, ja, ka1)*aaa5  + temp_i0(ia1,ja, ka1)*aaa6 &
                + temp_i0(ia, ja1,ka1)*aaa7  + temp_i0(ia1,ja1,ka1)*aaa8
        if(type.eq.0)then
          v(n) = v(n)*sqrt(temp_local/temp)
        end if
      else if(ispecies.eq.2)then
        sptfact = dns_a0(ia, ja, ka )*aaa1  + dns_a0(ia1,ja, ka )*aaa2 &
                + dns_a0(ia, ja1,ka )*aaa3  + dns_a0(ia1,ja1,ka )*aaa4 &
                + dns_a0(ia, ja, ka1)*aaa5  + dns_a0(ia1,ja, ka1)*aaa6 &
                + dns_a0(ia, ja1,ka1)*aaa7  + dns_a0(ia1,ja1,ka1)*aaa8
        temp_local = temp_a0(ia, ja, ka )*aaa1  + temp_a0(ia1,ja, ka )*aaa2 &
                + temp_a0(ia, ja1,ka )*aaa3  + temp_a0(ia1,ja1,ka )*aaa4 &
                + temp_a0(ia, ja, ka1)*aaa5  + temp_a0(ia1,ja, ka1)*aaa6 &
                + temp_a0(ia, ja1,ka1)*aaa7  + temp_a0(ia1,ja1,ka1)*aaa8
        if(type.eq.0)then
          v(n) = v(n)*sqrt(temp_local/temp)
        end if
      end if
!2016-08-06e

        babse = sqrt(bre**2 + bze**2 + bphie**2)
        gc(4,n) = m*v(n)*lambda(n)
        ekin = 0.50d0*m*v(n)**2
        gc(5,n)= ekin*(1.0d0 - lambda(n)**2)/babse
        clambda = gc(5,n)*b0/ekin !2012-06-17
        gc(8,n) = e*psip + gc(4,n)*gc(1,n)*bphie/babse

! comment for itpa benchmark
! at t=0, abs(b) = abs(b0)

!        rminor = sqrt( (gc(1,n)-raxis)**2 &
!                     + ((gc(2,n)-0.50d0*zleng) )**2 &
!                     )
!        bmax = b0*raxis/(raxis-rminor)
!        vpara = sqrt(2.0d0*(ekin-gc(5,n)*bmin)/m) &
!               *0.50d0*(1.0d0 + sign(1.0d0, ekin-gc(5,n)*bmax) ) &
!                      *sign(1.0d0,lambda(n) )
!        pphi_n = gc(8,n) - m*raxis*vpara
!        sptfact = exp(pphi_n/(e*psimax*0.37d0) )

! distribution function

      if(type.eq.0)then
        gc(9,n) = sptfact*exp(-0.50d0*m*v(n)**2/temp_local)*temp_local**(-1.5d0)
        gc(10,n) = gc(1,n)*(temp_local/temp)**1.50d0
      else if(type.eq.1.or.type.eq.2)then
        gc(9,n) = sptfact/(v(n)**3 + vcritp**3) &
                *0.50d0*erfc((v(n)-vbirth)/deltav)
        gc(10,n) = gc(1,n)
      else if(type.eq.3)then !2012-06-17
        gc(9,n) = sptfact/(v(n)**3 + vcritp**3) &
                *0.50d0*erfc((v(n)-vbirth)/deltav) &
                *exp(-(clambda-clambda0)**2/dclambda**2)
!                *(sign(0.50d0, gc(4,n)) + 0.50d0) !2013-03-02
        gc(10,n) = gc(1,n)
      end if

      end do

! check
!      if(my_rank.eq.0)then

!        pphi_total = 0.0d0
!        do n = 1, marker_num0
!          pphi_total = pphi_total + gc(8,n)
!        end do
!        pphi_total = pphi_total/dble(marker_num0)

!        write(7,*)'pphi_min,pphi_max,pphi_average=', &
!                   pphi_min,pphi_max,pphi_total
!      end if

! --------- use randum number for buffer particles --------- !
      nrandom = marker_num - marker_num0
!      nrandom_total = nrandom*mpi_proc

! for major radius
      call ransuu(nrandom,rand_no)


      do n = 1, nrandom
        n_buf = n + marker_num0
        gc(1,n_buf) = rmin + rand_no(n)*dpartr
      end do

! for vertical position
      call ransuu(nrandom,rand_no)


      do n = 1, nrandom
        n_buf = n + marker_num0
        gc(2,n_buf) = zmin + (zmax - zmin)*rand_no(n)
      end do

! for toroidal angle
      call ransuu(nrandom,rand_no)

      do n = 1, nrandom
        n_buf = n + marker_num0
        gc(3,n_buf) = phimin + (phimax - phimin)*rand_no(n)
      end do

      do n = marker_num0 + 1, marker_num
        v(n) = 1.0d-2
        lambda(n) = 1.0d0
        gc(4,n) = m*v(n)*lambda(n)
        gc(5,n)= 0.0d0
        gc(7,n) = 0.0d0
        gc(8,n) = 0.0d0
        gc(9,n) = 0.0d0
        gc(10,n) = gc(1,n)
        node(n) = node_up
      end do

! moved after the definition of unused particles, 2016-01-06
      call bcptcl(marker_num,gc,node)

! normalize ffp

      number_p = 0.0d0
      energy_p = 0.0d0
      enpara_p = 0.0d0 !2017-09-13
      r_p = 0.0d0
      z_p = 0.0d0
      phi_p = 0.0d0
      mom_p = 0.0d0

      do n = 1, marker_num0
        ftotal = gc(7,n)*gc(9,n)*gc(10,n)
        number_p = number_p + ftotal
        energy_p = energy_p + m*v(n)**2*0.50d0*ftotal
        enpara_p = enpara_p + gc(4,n)**2*0.50d0/m*ftotal !2017-09-13
        r_p = r_p + gc(1,n)*ftotal
        z_p = z_p + gc(2,n)*ftotal
        phi_p = phi_p + gc(3,n)*ftotal
        mom_p = mom_p + gc(4,n)*ftotal
      end do

      call mpi_allreduce(number_p,total_n,1,mpi_real8 &
                        ,mpi_sum,mpi_comm_world,mpi_err)
      call mpi_allreduce(energy_p,total_e,1,mpi_real8 &
                        ,mpi_sum,mpi_comm_world,mpi_err)
!2017-09-13s
      call mpi_allreduce(enpara_p,total_epara,1,mpi_real8 &
                        ,mpi_sum,mpi_comm_world,mpi_err)
!2017-09-13e

      call mpi_allreduce(r_p,total_r,1,mpi_real8 &
                        ,mpi_sum,mpi_comm_world,mpi_err)
      call mpi_allreduce(z_p,total_z,1,mpi_real8 &
                        ,mpi_sum,mpi_comm_world,mpi_err)
      call mpi_allreduce(phi_p,total_phi,1,mpi_real8 &
                        ,mpi_sum,mpi_comm_world,mpi_err)
      call mpi_allreduce(mom_p,total_mom,1,mpi_real8 &
                        ,mpi_sum,mpi_comm_world,mpi_err)

!2017-09-13s
!      if(type.eq.0)then
!        vnrm = vmax/sqrt(2.0d0*temp/m)
!        fraction = erf(vnrm) - 2.0d0/sqrt(pi)*vnrm*exp(-vnrm**2) &
!                             - 4.0d0/(3.0d0*sqrt(pi))*vnrm**3*exp(-vnrm**2)
!      else 
        fraction = 1.0d0
!      end if

      factor = fraction*total_energy/total_e
      temp_p = 2.0d0*total_e/total_n/3.0d0

      if(my_rank.eq.0)then
        write(7,'(a,1pe14.5)')'total_e=',total_e
        write(7,'(a,1pe14.5)')'total_epara=',total_epara
        write(7,'(a,1pe14.5)')'temp_p=',temp_p
        write(7,'(a,1pe14.5)')'average r=',total_r/total_n
        write(7,'(a,1pe14.5)')'average z=',total_z/total_n
        write(7,'(a,1pe14.5)')'average phi=',total_phi/total_n
        write(7,'(a,1pe14.5)')'average mom=',total_mom/total_n
        write(7,'(a,1pe14.5)')'factor = ',factor
      end if

! for e, i, a, D, T with type=0, 2021-06-19
      if(type.eq.0)then
        do n = 1, marker_num
          gc(6,n) = 0.0d0
          gc(10,n) = gc(10,n)*factor
        end do
      else

! for ion, 2021-05-06
      if(ispecies.eq.1.and.type.gt.0)then
        do k = 1, lphi
        do j = 1, lz
        do i = 1, lr

          dns_i(i,j,k) = dns_i(i,j,k)*temp/temp_p
          dns_i(i,j,k) = max(dns_i_min, dns_i(i,j,k) )
          dns_i0(i,j,k) = dns_i(i,j,k)

          ptotal = ppara_i0(i,j,k) + 2.0d0*pperp_i0(i,j,k)
          ppara_i0(i,j,k) = ptotal*total_epara/total_e
          pperp_i0(i,j,k) =(ptotal - ppara_i0(i,j,k) )*0.50d0
          ppara_i(i,j,k) = ppara_i0(i,j,k)
          pperp_i(i,j,k) = pperp_i0(i,j,k)

          temp_i0(i,j,k) = max(temp_i_min, &
                          (ppara_i0(i,j,k) + 2.0d0*pperp_i0(i,j,k) ) &
                         /(3.0d0*dns_i0(i,j,k) ) &
                              )
        end do
        end do
        end do

        call gradient(1,dns_i0,dns_i0_r,dns_i0_z,dns_i0_phi)
        call gradient(1,temp_i0,temp_i0_r,temp_i0_z,temp_i0_phi)

        do n = 1, marker_num
          gc(6,n) = 0.0d0
          gc(9,n) = gc(9,n)*temp/temp_p
          gc(10,n) = gc(10,n)*factor*temp_p/temp
        end do
      end if

! for alpha, 2021-05-06
      if(ispecies.eq.2.and.type.gt.0)then
        do k = 1, lphi
        do j = 1, lz
        do i = 1, lr
          dns_a(i,j,k) = dns_a(i,j,k)*temp/temp_p
          dns_a(i,j,k) = max(dns_a_min, dns_a(i,j,k) )
          dns_a0(i,j,k) = dns_a(i,j,k)

          ptotal = ppara_a0(i,j,k) + 2.0d0*pperp_a0(i,j,k)
          ppara_a0(i,j,k) = ptotal*total_epara/total_e
          pperp_a0(i,j,k) =(ptotal - ppara_a0(i,j,k) )*0.50d0
          ppara_a(i,j,k) = ppara_a0(i,j,k)
          pperp_a(i,j,k) = pperp_a0(i,j,k)

          temp_a0(i,j,k) = max(temp_a_min, &
                          (ppara_a0(i,j,k) + 2.0d0*pperp_a0(i,j,k) ) &
                         /(3.0d0*dns_a0(i,j,k) ) &
                              )
        end do
        end do
        end do

        call gradient(1,dns_a0,dns_a0_r,dns_a0_z,dns_a0_phi)
        call gradient(1,temp_a0,temp_a0_r,temp_a0_z,temp_a0_phi)

        do n = 1, marker_num
          gc(6,n) = 0.0d0
          gc(9,n) = gc(9,n)*temp/temp_p
          gc(10,n) = gc(10,n)*factor*temp_p/temp
        end do
      end if

      end if

end
!--------------------------------------------------------------------
subroutine scattering
! pitch-angle scattering and energy diffusion of energetic particles
!--------------------------------------------------------------------
      use mpiset
      use field
      use grid
      use particle
      use random_num
      implicit none

      real(8)::dr1,dz1,dphi1
      integer::n,ia,ia1,ja,ja1,ka,ka1,nrandom
      real(8)::ar,ar1,az,az1,aphi,aphi1
      real(8)::aaa1,aaa2,aaa3,aaa4,aaa5,aaa6,aaa7,aaa8
      real(8)::kdt,nudt,lambda_old,lambda_new
      real(8)::vcrit3,ti_a,te_a,dl_pa,dv_drag,dv_dif
!      real(8)::pa_min,pa_max,mu_min,mu_max
!      real(8)::lambda_old_min,lambda_old_max
!      real(8)::lambda_new_min,lambda_new_max
      real(8)::bre,bze,bphie,babse
      real(8)::rand_no(marker_each)

      nrandom = marker_a

      call ransuu(nrandom,rand_no)

!$omp parallel do
      do n = 1, marker_a
         lambda_a(n) = sign(1.0d0, rand_no(n) - 0.50d0)
      end do

      call ransuu(nrandom,rand_no)


      dr1 = 1.0d0/dr
      dz1 = 1.0d0/dz
      dphi1 = 1.0d0/dphi
      kdt = dble(kscatter)*dt

!$omp parallel private(n,ia,ia1,ja,ja1,ka,ka1,ar1,ar,az1,az,aphi1,aphi &
!$omp& ,aaa1,aaa2,aaa3,aaa4,aaa5,aaa6,aaa7,aaa8,bre,bze,bphie,babse &
!$omp& ,nudt,vcrit3,ti_a,te_a,dl_pa,dv_drag,dv_dif,lambda_old,lambda_new)
!$omp do
      do n = 1, marker_a

        if(node_a(n).eq.my_rank)then

! gc_x(:,n) = 1:r, 2:z, 3:phi, 4:p, 5:mu, 6:weight, 7:active, 8:p_phi, 9:fff, 10:fnrml
        ia = max(1, min(lr-1, &
                 int( ( gc_a(1,n)-(major_r - minor_r) )*dr1) + 1 &
                      - kr_offset(my_rank) ) )
        ia1= ia + 1
        ja = max(1, min(lz-1, &
                 int( gc_a(2,n)*dz1 ) + 1 - kz_offset(my_rank) ) )
        ja1= ja + 1
        ka = max(1, min(lphi-1, &
                 int( gc_a(3,n)*dphi1 ) + 3 - kphi_offset(my_rank) ) )
        ka1= ka + 1

        ar1 = max(0.0d0, min(1.0d0, (gc_a(1,n)-grr(ia,ja,ka) ) *dr1))
        ar  = 1.0d0 - ar1
        az1 = max(0.0d0, min(1.0d0, (gc_a(2,n)-gzz(ia,ja,ka) ) *dz1))
        az  = 1.0d0 - az1
        aphi1 = max(0.0d0, min(1.0d0, (gc_a(3,n)-gphi(ia,ja,ka) ) *dphi1))
        aphi  = 1.0d0 - aphi1

        aaa1 = ar *az *aphi
        aaa2 = ar1*az *aphi
        aaa3 = ar *az1*aphi
        aaa4 = ar1*az1*aphi
        aaa5 = ar *az *aphi1
        aaa6 = ar1*az *aphi1
        aaa7 = ar *az1*aphi1
        aaa8 = ar1*az1*aphi1

! fields at each particle position

        bre = br(ia, ja, ka )*aaa1 + br(ia1,ja, ka )*aaa2 &
            + br(ia, ja1,ka )*aaa3 + br(ia1,ja1,ka )*aaa4 &
            + br(ia, ja, ka1)*aaa5 + br(ia1,ja, ka1)*aaa6 &
            + br(ia, ja1,ka1)*aaa7 + br(ia1,ja1,ka1)*aaa8
        bze = bz(ia, ja, ka )*aaa1 + bz(ia1,ja, ka )*aaa2 &
            + bz(ia, ja1,ka )*aaa3 + bz(ia1,ja1,ka )*aaa4 &
            + bz(ia, ja, ka1)*aaa5 + bz(ia1,ja, ka1)*aaa6 &
            + bz(ia, ja1,ka1)*aaa7 + bz(ia1,ja1,ka1)*aaa8
        bphie = bphi(ia, ja, ka )*aaa1 + bphi(ia1,ja, ka )*aaa2 &
            + bphi(ia, ja1,ka )*aaa3 + bphi(ia1,ja1,ka )*aaa4 &
            + bphi(ia, ja, ka1)*aaa5 + bphi(ia1,ja, ka1)*aaa6 &
            + bphi(ia, ja1,ka1)*aaa7 + bphi(ia1,ja1,ka1)*aaa8
        babse = sqrt(bre**2 + bze**2 + bphie**2)

        v_a(n) = sqrt(2.0d0*gc_a(5,n)*babse/m_a + (gc_a(4,n)/m_a)**2)


        nudt =(nusd(ia, ja, ka )*aaa1 + nusd(ia1,ja, ka )*aaa2 &
             + nusd(ia, ja1,ka )*aaa3 + nusd(ia1,ja1,ka )*aaa4 &
             + nusd(ia, ja, ka1)*aaa5 + nusd(ia1,ja, ka1)*aaa6 &
             + nusd(ia, ja1,ka1)*aaa7 + nusd(ia1,ja1,ka1)*aaa8 &
              )*kdt

        vcrit3=(vcrit(ia, ja, ka )*aaa1 + vcrit(ia1,ja, ka )*aaa2 &
              + vcrit(ia, ja1,ka )*aaa3 + vcrit(ia1,ja1,ka )*aaa4 &
              + vcrit(ia, ja, ka1)*aaa5 + vcrit(ia1,ja, ka1)*aaa6 &
              + vcrit(ia, ja1,ka1)*aaa7 + vcrit(ia1,ja1,ka1)*aaa8 &
               )**3/v_a(n)**3

        ti_a =(ti(ia, ja, ka )*aaa1 + ti(ia1,ja, ka )*aaa2 &
             + ti(ia, ja1,ka )*aaa3 + ti(ia1,ja1,ka )*aaa4 &
             + ti(ia, ja, ka1)*aaa5 + ti(ia1,ja, ka1)*aaa6 &
             + ti(ia, ja1,ka1)*aaa7 + ti(ia1,ja1,ka1)*aaa8 &
              )
        te_a =(te(ia, ja, ka )*aaa1 + te(ia1,ja, ka )*aaa2 &
             + te(ia, ja1,ka )*aaa3 + te(ia1,ja1,ka )*aaa4 &
             + te(ia, ja, ka1)*aaa5 + te(ia1,ja, ka1)*aaa6 &
             + te(ia, ja1,ka1)*aaa7 + te(ia1,ja1,ka1)*aaa8 &
              )

        dl_pa = nudt*vcrit3*gc_a(7,n)*zeff  ! factor 2*M_i/(2*M_fast)=1
        dv_drag = nudt/m_a*(te_a - 0.50d0*ti_a*vcrit3 )/v_a(n)*gc_a(7,n)
        dv_dif  = nudt/m_a*(te_a + ti_a*vcrit3 )*gc_a(7,n)

        lambda_old = max(-1.0d0, min(1.0d0, gc_a(4,n)/(v_a(n)*m_a) ) )
        lambda_new = lambda_old*(1.0d0 - dl_pa) + lambda_a(n)*sqrt((1.0d0-lambda_old**2)*dl_pa)

! be careful for the sign of dv_drag (simlar to the drag in pitch angle)
        v_a(n) = v_a(n) + dv_drag + sign(1.0d0, rand_no(n) - 0.50d0)*sqrt(dv_dif)

        gc_a(4,n) = v_a(n)*lambda_new*m_a
        gc_a(5,n) = 0.50d0*m_a*(1.0d0 - lambda_new**2)*v_a(n)**2/babse

        end if
      end do
!$omp end parallel


end
!--------------------------------------------------------------------
subroutine injection
!--------------------------------------------------------------------
      use mpiset
      use field
      use particle
      use grid
      use random_num
      implicit none

      integer::n,n_inj0
      integer::my_rank_deposit,nbeam,nrank_r,nrank_z
      integer::nstart

      integer::ia,ia1,ja,ja1,ka
      real(8)::ar,ar1,az,az1,aaa1,aaa2,aaa3,aaa4
      real(8)::bre,bze,bphie,psip,babse,ekin

      n_inj0 = n_inj
      do n = n_inj0 + 1, n_bd
        if(t_bd(n).lt.t.and.injection_bd(n).eq.0)then
!          if(my_rank.eq.0) write(7,*)'t_bd(',n,'),t=',t_bd(n),t
          injection_bd(n) = 1
        else
          exit
        end if
      end do
      n_inj = n - 1

!--------------------------------------------------------------------
! hereafer, substitute into the simulation arrays
!--------------------------------------------------------------------
!   real(4), allocatable::r(:),z(:),lambda(:),e(:),zeta(:),v(:)
!   real(4), allocatable::bdp_data(:,:),t_bd(:)
!   integer, allocatable::injection_bd(:)

      my_rank_deposit = my_rank_phi + mpi_proc_phi*my_rank_ptcl
      nstart = 1

      do nbeam = n_inj0 + 1, n_inj
        nrank_r = max(0, min(mpi_proc_r-1, &
               int(dble(mpi_proc_r)*(bdp_data(1,nbeam)-(major_r-minor_r) )/rleng ) &
                           ) ) 
        nrank_z = max(0, min(mpi_proc_z-1, &
               int(dble(mpi_proc_z)*bdp_data(2,nbeam)/zleng ) &
                           ) )

        if(nrank_r.eq.my_rank_r.and.nrank_z.eq.my_rank_z)then
          if(mod(nbeam,mpi_proc_phi*mpi_proc_ptcl).eq.my_rank_deposit)then
            do n = nstart, marker_a
              if(node_a(n).ne.my_rank) exit
            end do
            nstart = n
            if(n.eq.marker_a)then
              write(7,*)'number of particles reached marker_a at rank',my_rank
              exit
            end if

! gc_x(:,n) = 1:r, 2:z, 3:phi, 4:p, 5:mu, 6:weight, 7:active, 8:p_phi, 9:fff, 10:fnrml
            gc_a(1,n) = bdp_data(1,nbeam)
            gc_a(2,n) = bdp_data(2,nbeam)
!            lambda_a(n) = bdp_data(3,nbeam)
            lambda_a(n) =-bdp_data(3,nbeam) !2013-11-22 bdp_data is with respect to current direction in DIII-D
            v_a(n) = bdp_data(4,nbeam)
            gc_a(3,n) = phimin + (phimax - phimin)*bdp_data(5,nbeam)
            gc_a(4,n) = m_a*v_a(n)*lambda_a(n)
            gc_a(7,n) = 1.0d0
            node_a(n) = my_rank

            ia = max(1, min(lr-1, &
                 int( ( gc_a(1,n)-(major_r - minor_r) )/dr) + 1 - kr_offset(my_rank) ))
            ia1= ia + 1
            ja = max(1, min(lz-1, &
                 int( gc_a(2,n)/dz ) + 1 - kz_offset(my_rank) ) )
            ja1= ja + 1
            ka = 3

            ar1 = (gc_a(1,n)-grr(ia,ja,ka) ) /dr
            ar  = 1.0d0 - ar1
            az1 = (gc_a(2,n)-gzz(ia,ja,ka) ) /dz
            az  = 1.0d0 - az1

            aaa1 = ar *az
            aaa2 = ar1*az
            aaa3 = ar *az1
            aaa4 = ar1*az1

            bre = br0(ia, ja, ka )*aaa1  + br0(ia1,ja, ka )*aaa2 &
                + br0(ia, ja1,ka )*aaa3  + br0(ia1,ja1,ka )*aaa4
            bze = bz0(ia, ja, ka )*aaa1  + bz0(ia1,ja, ka )*aaa2 &
                + bz0(ia, ja1,ka )*aaa3  + bz0(ia1,ja1,ka )*aaa4
            bphie = bphi0(ia, ja, ka )*aaa1  + bphi0(ia1,ja, ka )*aaa2 &
                  + bphi0(ia, ja1,ka )*aaa3  + bphi0(ia1,ja1,ka )*aaa4
            psip  = psi(ia, ja, ka )*aaa1  + psi(ia1,ja, ka )*aaa2 &
                  + psi(ia, ja1,ka )*aaa3  + psi(ia1,ja1,ka )*aaa4

            babse = sqrt(bre**2 + bze**2 + bphie**2)
            ekin = 0.50d0*m_a*v_a(n)**2
! gc_x(:,n) = 1:r, 2:z, 3:phi, 4:p, 5:mu, 6:weight, 7:active, 8:p_phi, 9:fff, 10:fnrml
            gc_a(5,n)= ekin*(1.0d0 - lambda_a(n)**2)/babse
            gc_a(8,n) = e_a*psip + gc_a(4,n)*gc_a(1,n)*bphie/babse
          end if
        end if
      end do

! remove ash particles
      do n = 1, marker_a
        if(v_a(n).lt.1.0d-1*valpha)then
          gc_a(7,n) = 0.0d0
          node_a(n) = mod(my_rank + 1 + nprocess, nprocess)
        end if
      end do

end
!--------------------------------------------------------------------
subroutine t_integration(istep,type,marker_num &
                        ,gc,dgc,gc1,gc2  )
! modified on 2015-06-23
! gc_x(:,n) = 1:r, 2:z, 3:phi, 4:p, 5:mu, 6:weight, 7:active, 8:p_phi, 9:fff, 10:fnrml
!--------------------------------------------------------------------
      use parameters
      implicit none

      integer::marker_num,istep,type !2016-01-09
      real(8)::gc(ngc2,marker_num)
      real(8)::dgc(ngc1,marker_num)
      real(8)::gc1(ngc1,marker_num)
      real(8)::gc2(ngc1,marker_num)
      real(8)::c0,c1
      integer::n

      if(istep.eq.1)then
        c0 = 1.0d0/6.0d0
        c1 = 0.50d0

#ifdef AMDGPU
!$omp target teams distribute parallel do private(n)
#else        
!$omp parallel do private(n)
#endif
        do n = 1, marker_num
          gc1(1,n) = gc(1,n)
          gc1(2,n) = gc(2,n)
          gc1(3,n) = gc(3,n)
          gc1(4,n) = gc(4,n)
          gc1(6,n) = gc(6,n)
          gc1(8,n) = gc(8,n)

          gc2(1,n) = gc(1,n)
          gc2(2,n) = gc(2,n)
          gc2(3,n) = gc(3,n)
          gc2(4,n) = gc(4,n)
          gc2(6,n) = gc(6,n)
          gc2(8,n) = gc(8,n)
        end do

      else if(istep.eq.2)then
        c0 = 1.0d0/3.0d0
        c1 = 0.50d0
      else if(istep.eq.3)then
        c0 = 1.0d0/3.0d0
        c1 = 1.0d0
      else if(istep.eq.4)then
        c0 = 1.0d0/6.0d0
      end if

      if(istep.ne.4)then

#ifdef AMDGPU
!$omp target teams distribute parallel do private(n)
#else        
!$omp parallel do private(n)
#endif
        do n = 1, marker_num
          gc1(1,n) = gc1(1,n) + c0*dgc(1,n)
          gc1(2,n) = gc1(2,n) + c0*dgc(2,n)
          gc1(3,n) = gc1(3,n) + c0*dgc(3,n)
          gc1(4,n) = gc1(4,n) + c0*dgc(4,n)
          gc1(6,n) = gc1(6,n) + c0*dgc(6,n)
          gc1(8,n) = gc1(8,n) + c0*dgc(8,n)

          gc(1,n) = gc2(1,n) + c1*dgc(1,n)
          gc(2,n) = gc2(2,n) + c1*dgc(2,n)
          gc(3,n) = gc2(3,n) + c1*dgc(3,n)
          gc(4,n) = gc2(4,n) + c1*dgc(4,n)
          gc(6,n) = gc2(6,n) + c1*dgc(6,n)
          gc(8,n) = gc2(8,n) + c1*dgc(8,n)
        end do

      else


#ifdef AMDGPU
!$omp target teams distribute parallel do private(n)
#else        
!$omp parallel do private(n)
#endif
        do n = 1, marker_num
          gc(1,n) = gc1(1,n) + c0*dgc(1,n)
          gc(2,n) = gc1(2,n) + c0*dgc(2,n)
          gc(3,n) = gc1(3,n) + c0*dgc(3,n)
          gc(4,n) = gc1(4,n) + c0*dgc(4,n)
          gc(6,n) = gc1(6,n) + c0*dgc(6,n)
          gc(8,n) = gc1(8,n) + c0*dgc(8,n)
        end do

      end if

! for beam injection
      if(type.eq.-5)then
        if(istep.eq.1)then
          c0 = 1.0d0/6.0d0
          c1 = 0.50d0

#ifdef AMDGPU
!$omp target teams distribute parallel do private(n)
#else        
!$omp parallel do private(n)
#endif
          do n = 1, marker_num
            gc1(5,n)  = gc(5,n)
            gc2(5,n)  = gc(5,n)
            gc1(10,n) = gc(10,n)
            gc2(10,n) = gc(10,n)
          end do
        else if(istep.eq.2)then
          c0 = 1.0d0/3.0d0
          c1 = 0.50d0
        else if(istep.eq.3)then
          c0 = 1.0d0/3.0d0
          c1 = 1.0d0
        else if(istep.eq.4)then
          c0 = 1.0d0/6.0d0
        end if

        if(istep.ne.4)then

#ifdef AMDGPU
!$omp target teams distribute parallel do private(n)
#else        
!$omp parallel do private(n)
#endif
          do n = 1, marker_num
            gc1(5,n) = gc1(5,n)  + c0*dgc(5,n)
            gc(5,n)  = gc2(5,n)  + c1*dgc(5,n)
            gc1(10,n)= gc1(10,n) + c0*dgc(10,n)
            gc(10,n) = gc2(10,n) + c1*dgc(10,n)
          end do

        else

#ifdef AMDGPU
!$omp target teams distribute parallel do private(n)
#else        
!$omp parallel do private(n)
#endif
          do n = 1, marker_num
            gc(5,n)  = gc1(5,n)  + c0*dgc(5,n)
            gc(10,n) = gc1(10,n) + c0*dgc(10,n)
          end do
        end if

      end if

end
!--------------------------------------------------------------------
subroutine com_particle(marker_num,gc,node,node_now)
!     communication of particle data
! modified on 2015-06-23
! gc_x(:,n) = 1:r, 2:z, 3:phi, 4:p, 5:mu, 6:weight, 7:active, 8:p_phi, 9:fff, 10:fnrml
!--------------------------------------------------------------------
      use mpiset
      use grid
      implicit none

      integer::marker_num
      real(8)::gc(ngc2,marker_num)
      integer::node(marker_num),node_now(marker_num)

      integer::n,idirection

!$omp parallel do
      do n = 1, marker_num
        gc(3,n) = mod(gc(3,n)+phileng,phileng)
      end do

      if(mpi_proc_r.ge.2)then
        call com_particle3(0, &
             marker_num,gc,node,node_now)
      end if

      if(mpi_proc_z.ge.2)then
        call com_particle3(1, &
             marker_num,gc,node,node_now)
      end if

      if(mpi_proc_phi.ge.2)then
        call com_particle3(2, &
             marker_num,gc,node,node_now)
      end if

end
!--------------------------------------------------------------------
subroutine bcptcl(marker_num,gc,node)
! boundary condition for particles
! gc_x(:,n) = 1:r, 2:z, 3:phi, 4:p, 5:mu, 6:weight, 7:active, 8:p_phi, 9:fff, 10:fnrml
! new version, 2016-01-08
!--------------------------------------------------------------------
      use mpiset
      use grid
      use check
      implicit none

      integer::marker_num
      real(8)::gc(ngc2,marker_num)
      integer::node(marker_num)
      integer::n,node_up
      real(8)::rin,rout,zin,zout,phiin,phiout
      real(8)::rspan,zspan,phispan,rat
#ifdef AMDGPU
      real(8)::t0,t1
#endif      
      node_up = mod(my_rank + 1 + nprocess, nprocess)

! periodic in phi-direction
!
!      do n = 1, marker_num
!        gc(3,n) = mod(gc(3,n)+phileng,phileng)
!      end do
!
! wall condition at both the r-edge and z-edge
!
      rin  = grr(1,1,1)
      rout = grr(lr,1,1)
      zin  = gzz(1,1,1)
      zout = gzz(1,lz,1)
      phiin= gphi(1,1,1)
      phiout= gphi(1,1,lphi)

      rspan = rout - rin
      zspan = zout - zin
      phispan = phiout - phiin
 
#ifdef AMDGPU
!$omp target teams distribute parallel do private(n,rat,t0,t1)
#else
!$omp parallel do private(rat)
#endif
      do n = 1, marker_num
        if((gc(1,n).gt.rout.or.gc(1,n).lt.rin.or.gc(2,n).gt.zout.or.gc(2,n).lt.zin &
       .or.gc(3,n).gt.phiout.or.gc(3,n).lt.phiin).and.node(n).eq.my_rank)then !2015-08-28

          rat = (dble(n)-0.50d0)/dble(marker_num)
#ifdef AMDGPU
          ! mod for fp is not supported!
          t0 = rspan*rat*dble(lzphi)
          t1 = rspan
          gc(1,n) = rin + (t0 - int(t0/t1)*t1)
          t0 = zspan*rat*dble(lphi)
          t1 = zspan
          gc(2,n) = zin + (t0 - int(t0/t1)*t1)
#else
          gc(1,n) = rin + mod(rspan*rat*dble(lzphi), rspan)
          gc(2,n) = zin + mod(zspan*rat*dble(lphi), zspan)
#endif
          gc(3,n) = phiin + phispan*rat
          gc(7,n) = 0.0d0
          node(n) = node_up

        end if

      end do


end
!--------------------------------------------------------------------
subroutine satellite
!--------------------------------------------------------------------
   use parameters
   use particle, only:gc_i,node_i,flp_i &
                     ,gc_a,node_a,flp_a
   use gyro
   use mpiset

   implicit none

#ifdef KTI   
!ion
        call create_satellite(marker_i,marker_i_gyro,m_i,e_i &
                             ,gc_i,gyro_i)
#endif
        
!alpha
        call create_satellite(marker_a,marker_a_gyro,m_a,e_a &
                             ,gc_a,gyro_a)


end
!--------------------------------------------------------------------
subroutine create_satellite(marker_num,marker_num_gyro,mass,charge &
                           ,gc,gyro_phys)
!2015-06-26 version
!   gc_x(:,n) = 1:r, 2:z, 3:phi, 4:p, 5:mu, 6:weight, 7:active, 8:p_phi, 9:fff, 10:fnrml
!--------------------------------------------------------------------
      use parameters
      use field, only:br,bz,bphi,babs
      use grid
      use mpiset
      implicit none

      integer::marker_num,marker_num_gyro
      real(8)::mass,charge,charge1
      real(8)::gc(ngc2,marker_num)
      real(8)::gyro_phys(2,marker_num_gyro)

      real(8)::dr1,dz1,dphi1,babse,rho_g
      integer::n,ia,ja,ka
      real(8)::ar,ar1,az,az1,aphi,aphi1
      real(8)::aaa1,aaa2,aaa3,aaa4,aaa5,aaa6,aaa7,aaa8
      integer::kr,kz,kphi,ip
      real(8)::ma_mi_r

      dr1 = 1.0d0/dr
      dz1 = 1.0d0/dz
      dphi1 = 1.0d0/dphi
      charge1 = 1.0d0/charge

      kr   = 1 - kr_offset(my_rank)
      kz   = 1 - kz_offset(my_rank)
      kphi = 1 + lphi_shd - kphi_offset(my_rank)
!      kphi = 3 - kphi_offset(my_rank)
      ma_mi_r = major_r-minor_r

#ifdef AMDGPU
!$omp  target teams distribute parallel do
!$omp& private(n,ia,ja,ka,ar1,ar,az1,az,aphi1,aphi,aaa1,aaa2,aaa3,aaa4,aaa5,aaa6,aaa7,aaa8)
!$omp& private(babse,rho_g,ip)
#else      
!$omp parallel private(n,ia,ja,ka,ar1,ar,az1,az,aphi1,aphi &
!$omp& ,aaa1,aaa2,aaa3,aaa4,aaa5,aaa6,aaa7,aaa8,babse,rho_g,ip)
!$omp do schedule(dynamic,1000)
#endif
      do n = 1, marker_num
        ia=max(1,min(lr  -1,int((gc(1,n)-ma_mi_r)*dr1  ) + kr  ))
        ja=max(1,min(lz  -1,int(gc(2,n)          *dz1  ) + kz  ))
        ka=max(1,min(lphi-1,int(gc(3,n)        *dphi1) + kphi))

        ar1  = max(0.0d0, min(1.0d0, (gc(1,n) - ma_mi_r)*dr1 -  dble(ia - kr)  ) )
        ar   = 1.0d0 - ar1
        az1  = max(0.0d0, min(1.0d0, gc(2,n)*dz1 - dble(ja - kz) ) )
        az   = 1.0d0 - az1
        aphi1 = max(0.0d0, min(1.0d0, gc(3,n)*dphi1 - dble(ka - kphi) ) )
        aphi = 1.0d0 - aphi1

        aaa1 = ar *az *aphi
        aaa2 = ar1*az *aphi
        aaa3 = ar *az1*aphi
        aaa4 = ar1*az1*aphi
        aaa5 = ar *az *aphi1
        aaa6 = ar1*az *aphi1
        aaa7 = ar *az1*aphi1
        aaa8 = ar1*az1*aphi1

        babse = babs(ia, ja,  ka  )*aaa1 + babs(ia+1,ja,  ka  )*aaa2 &
              + babs(ia, ja+1,ka  )*aaa3 + babs(ia+1,ja+1,ka  )*aaa4 &
              + babs(ia, ja,  ka+1)*aaa5 + babs(ia+1,ja,  ka+1)*aaa6 &
              + babs(ia, ja+1,ka+1)*aaa7 + babs(ia+1,ja+1,ka+1)*aaa8
!        rho_g = sqrt(2.0d0*gc(5,n)*mass/babse)*charge1
        rho_g = sqrt(gc(5,n)*mass/babse)*charge1 !2015-07-10, satellites are located at 45 degree

        do ip = 1, 4
          gyro_phys(1, 4*(n-1)+ip) = gc(1,n) + rho_g*dble(1 - 2*mod(ip,2) ) ! -1, 1, -1, 1
          gyro_phys(2, 4*(n-1)+ip) = gc(2,n) + rho_g*dble(-1 +2*int( (ip-1)/2 ) ) ! -1, -1, 1, 1
        end do

      end do
#ifndef AMDGPU
!$omp end do
!$omp end parallel
#endif
      
end

!--------------------------------------------------------------------
#ifndef KTI
!subroutines mhd, e_field, and write_check
!--------------------------------------------------------------------


!--------------------------------------------------------------------
subroutine mhd(istep,flag_HM)
!modified on 2025-09-06 for optimization to AMDGPU
!--------------------------------------------------------------------
      use parameters
      use mpiset
      use field, only:br,bz,bphi,babs,er,ez,ephi,cr,cz,cphi,cr0,cz0,cphi0 &
                     ,rho,prs,vr,vz,vphi &
                     ,ur,uz,uphi,div_u,omegar,omegaz,omegaphi &
                     ,utotphi,vtotphi &
                     ,prs_r,prs_z,prs_phi &
                     ,ppara_a,pperp_a,ppara_i,pperp_i &
                     ,gradbr,gradbz,gradbphi,curvbr,curvbz,curvbphi &
                     ,br0,bz0,bphi0,babs0,rho0,prs0 &
                     ,dbr,dbz,dbphi,dvr,dvz,dvphi,drho,dprs &
                     ,br1,bz1,bphi1,vr1,vz1,vphi1,rho1,prs1 &
                     ,br2,bz2,bphi2,vr2,vz2,vphi2,rho2,prs2 &
                     ,dvr0,dvz0,dvphi0,drho0,dprs0,vac &
!                     ,enrg,enrg0,denrg,enrg1,enrg2,denrg0 &
                     ,nu_spt,eta_spt,nu_n,chi,gamma,chi_spt,nun_spt &
                     ,vpir,vpiz,vpiphi &
                     ,vdrftr,vdrftz,vdrftphi & !2025-09-06
                     ,er_res,ez_res,ephi_res
      use grid, only:dt,grr,dr,dz,dphi
      use check
      implicit none

      logical::flag_HM !2013-01-12
      integer::istep
      integer::i,j,k,m !2025-03-05
      real(8)::de_trans_local(0:l_kin_species),de_dissp_local !2025-03-05
      real(8)::rot2_ur(lr,lz,lphi),rot2_uz(lr,lz,lphi),rot2_uphi(lr,lz,lphi)
      real(8)::div_ur(lr,lz,lphi),div_uz(lr,lz,lphi),div_uphi(lr,lz,lphi)
      real(8)::lpl3_rho(lr,lz,lphi),lpl3_prs(lr,lz,lphi) !2013-01-12
      real(8)::cnetr(lr,lz,lphi),cnetz(lr,lz,lphi),cnetphi(lr,lz,lphi)
!      real(8)::magc_e_r(lr,lz,lphi),magc_e_z(lr,lz,lphi),magc_e_phi(lr,lz,lphi)
      real(8)::magc_i_r(lr,lz,lphi),magc_i_z(lr,lz,lphi),magc_i_phi(lr,lz,lphi)
      real(8)::magc_a_r(lr,lz,lphi),magc_a_z(lr,lz,lphi),magc_a_phi(lr,lz,lphi)

      real(8)::c_ir(lr,lz,lphi),c_iz(lr,lz,lphi),c_iphi(lr,lz,lphi)
      real(8)::c_ar(lr,lz,lphi),c_az(lr,lz,lphi),c_aphi(lr,lz,lphi)

!      real(8)::eflux_r(lr,lz,lphi),eflux_z(lr,lz,lphi),eflux_phi(lr,lz,lphi) !2021-05-05
      real(8)::etransfer(lr,lz,lphi,0:l_kin_species) !2025-03-05
      real(8)::vparar(lr,lz,lphi),vparaz(lr,lz,lphi),vparaphi(lr,lz,lphi)
      real(8)::vparar_r(lr,lz,lphi),vparar_z(lr,lz,lphi),vparar_phi(lr,lz,lphi)
      real(8)::vparaz_r(lr,lz,lphi),vparaz_z(lr,lz,lphi),vparaz_phi(lr,lz,lphi)
      real(8)::vparaphi_r(lr,lz,lphi),vparaphi_z(lr,lz,lphi),vparaphi_phi(lr,lz,lphi)
      real(8)::rho_vr(lr,lz,lphi),rho_vz(lr,lz,lphi),rho_vphi(lr,lz,lphi)
      real(8)::drho1,dvr1,dvz1,dvphi1,dprs1,dbr1,dbz1,dbphi1,denrg1 !2021-05-05
      real(8)::b21,ppara_i1,ppara_a1 !2025-04-29
      real(8)::eidealr,eidealz,eidealphi,ekin !2021-05-05
      real(8)::c0,c1
      real(8)::c43=4.0d0/3.0d0
      real(8)::cmratio_ti=e_ti/m_ti,mcratio_ti=m_ti/e_ti !2025-04-29, for thermal ion
      real(8)::ubr,ubz,ubphi,vdotb,rhoinv
#ifdef AMDGPU
      real(8)::t1,t2,t3,t4,t5
#endif

! time integration
      if(istep.eq.1)then
        c0 = 1.0d0/6.0d0
        c1 = 0.50d0
      else if(istep.eq.2)then
        c0 = 1.0d0/3.0d0
        c1 = 0.50d0
      else if(istep.eq.3)then
        c0 = 1.0d0/3.0d0
        c1 = 1.0d0
      else if(istep.eq.4)then
        c0 = 1.0d0/6.0d0
      end if

! ------------- preparation ---------------

! magnetization current
      call gradient(1,pperp_i,c_ir,c_iz,c_iphi)
      call gradient(1,pperp_a,c_ar,c_az,c_aphi)

#ifdef AMDGPU      
!$omp target teams distribute parallel do private(i,b21,ppara_i1,ppara_a1)
#else
!$omp parallel do private(i,b21,ppara_i1,ppara_a1)
#endif
      do i = 1, lrzphi
        b21 = 1.0d0/babs(i,1,1)**2
        magc_i_r(i,1,1) =-(c_iphi(i,1,1)*bz(i,1,1) - c_iz(i,1,1)*bphi(i,1,1) )*b21
        magc_i_z(i,1,1) =-(c_ir(i,1,1)*bphi(i,1,1) - c_iphi(i,1,1)*br(i,1,1) )*b21
        magc_i_phi(i,1,1) =-(c_iz(i,1,1)*br(i,1,1) - c_ir(i,1,1)*bz(i,1,1) )*b21
        ppara_i1 =(ppara_i(i,1,1) - pperp_i(i,1,1) )/babs(i,1,1)

        c_ir(i,1,1) = ppara_i1*curvbr(i,1,1) &
                    + magc_i_r(i,1,1)
        c_iz(i,1,1) = ppara_i1*curvbz(i,1,1) &
                    + magc_i_z(i,1,1)
        c_iphi(i,1,1) = ppara_i1*curvbphi(i,1,1) &
                      + magc_i_phi(i,1,1)

        magc_a_r(i,1,1) =-(c_aphi(i,1,1)*bz(i,1,1) - c_az(i,1,1)*bphi(i,1,1) )*b21
        magc_a_z(i,1,1) =-(c_ar(i,1,1)*bphi(i,1,1) - c_aphi(i,1,1)*br(i,1,1) )*b21
        magc_a_phi(i,1,1) =-(c_az(i,1,1)*br(i,1,1) - c_ar(i,1,1)*bz(i,1,1) )*b21
        ppara_a1 =(ppara_a(i,1,1) - pperp_a(i,1,1) )/babs(i,1,1)
        c_ar(i,1,1) = ppara_a1*curvbr(i,1,1) &
                    + magc_a_r(i,1,1)
        c_az(i,1,1) = ppara_a1*curvbz(i,1,1) &
                    + magc_a_z(i,1,1)
        c_aphi(i,1,1) = ppara_a1*curvbphi(i,1,1) &
                      + magc_a_phi(i,1,1)

        cnetr(i,1,1) = cr(i,1,1) - c_ir(i,1,1) - c_ar(i,1,1)
        cnetz(i,1,1) = cz(i,1,1) - c_iz(i,1,1) - c_az(i,1,1)
        cnetphi(i,1,1) = cphi(i,1,1) - c_iphi(i,1,1) - c_aphi(i,1,1)
      end do

! ------------- preparation end ---------------

! energy evolution
#ifdef AMDGPU      
!$omp target teams distribute parallel do private(i,eidealr,eidealz,eidealphi)
#else
!$omp parallel do private(i,eidealr,eidealz,eidealphi)
#endif
      do i = 1, lrzphi
        eidealr =-(uphi(i,1,1)*bz(i,1,1) - uz(i,1,1)*bphi(i,1,1))
        eidealz =-(ur(i,1,1)*bphi(i,1,1) - uphi(i,1,1)*br(i,1,1))
        eidealphi =-(uz(i,1,1)*br(i,1,1) - ur(i,1,1)*bz(i,1,1))

!        ekin = 0.50d0*(ur(i,1,1)**2 + uz(i,1,1)**2 + uphi(i,1,1)**2)*rho(i,1,1) &
!             + gamma*prs(i,1,1)/(gamma - 1.0d0)
!        eflux_r(i,1,1) =(ur(i,1,1)*ekin &
!                       + eidealphi*bz(i,1,1) - eidealz*bphi(i,1,1) &
!                        )
!        eflux_z(i,1,1) =(uz(i,1,1)*ekin &
!                       + eidealr*bphi(i,1,1) - eidealphi*br(i,1,1) &
!                        )
!        eflux_phi(i,1,1) =(uphi(i,1,1)*ekin &
!                       + eidealz*br(i,1,1) - eidealr*bz(i,1,1) &
!                        )

        etransfer(i,1,1,1) =-( &
                               c_ir(i,1,1)*eidealr &
                             + c_iz(i,1,1)*eidealz &
                             + c_iphi(i,1,1)*eidealphi &
                             )
        etransfer(i,1,1,2) =-( &
                               c_ar(i,1,1)*eidealr &
                             + c_az(i,1,1)*eidealz &
                             + c_aphi(i,1,1)*eidealphi &
                             )
     end do
! 2025-03-05e

!      call divergence(0,eflux_r,eflux_z,eflux_phi,denrg)

! energy transfer from field to particles
#ifdef AMDGPU
      t1 = 0.0d0; t2 = 0.0d0; t3 = 0.0d0; t4 = 0.0d0; t5 = 0.0d0
#else
      de_trans_local = 0.0d0
#endif

#ifdef AMDGPU
!$omp target teams distribute parallel do private(k,j,i) collapse(3) reduction(+:t1,t2)
      do k = lphistart, lphiend
      do j = lzstart, lzend
      do i = lrstart, lrend
        t1 = t1 + grr(i,j,k)*etransfer(i,j,k,1)
        t2 = t2 + grr(i,j,k)*etransfer(i,j,k,2)
      end do
      end do
      end do
#else
!$omp parallel do collapse(3) reduction(+:de_trans_local)
!      do m = 1, l_kin_species
      do k = lphistart, lphiend
      do j = lzstart, lzend
      do i = lrstart, lrend
        de_trans_local(1) = de_trans_local(1) + grr(i,j,k)*etransfer(i,j,k,1)
        de_trans_local(2) = de_trans_local(2) + grr(i,j,k)*etransfer(i,j,k,2)
      end do
      end do
      end do
!      end do
#endif

#ifdef AMDGPU      
      de_trans(1) =-t1*dr*dz*dphi*dt
      de_trans(2) =-t2*dr*dz*dphi*dt
#else
      de_trans =-de_trans_local*dr*dz*dphi*dt
#endif

! ------------- loop of preparation for HM simulation ---------------
      if(flag_HM)then
         
#ifdef AMDGPU
!$omp target teams distribute parallel do private(i,b21,ubr,ubz,ubphi,vdotb,rhoinv)
#else
!$omp parallel do private(b21,ubr,ubz,ubphi,vdotb,rhoinv)
#endif
      do i = 1, lrzphi
        b21 = 1.0d0/babs(i,1,1)**2
        ubr = br(i,1,1)/babs(i,1,1)
        ubz = bz(i,1,1)/babs(i,1,1)
        ubphi = bphi(i,1,1)/babs(i,1,1)
! vpi{r,z,phi}: diamag current
        vpir(i,1,1) =-(prs_phi(i,1,1)*bz(i,1,1) - prs_z(i,1,1)*bphi(i,1,1) )*b21
        vpiz(i,1,1) =-(prs_r(i,1,1)*bphi(i,1,1) - prs_phi(i,1,1)*br(i,1,1) )*b21
        vpiphi(i,1,1) =-(prs_z(i,1,1)*br(i,1,1) - prs_r(i,1,1)*bz(i,1,1) )*b21

        rho_vr(i,1,1) =(0.50d0*vpir(i,1,1)*mcratio_ti + rho(i,1,1)*ur(i,1,1) )*vac(i,1,1)
        rho_vz(i,1,1) =(0.50d0*vpiz(i,1,1)*mcratio_ti + rho(i,1,1)*uz(i,1,1) )*vac(i,1,1)
        rho_vphi(i,1,1)=(0.50d0*vpiphi(i,1,1)*mcratio_ti + rho(i,1,1)*utotphi(i,1,1) )*vac(i,1,1)

!        rhoinv = 1.0d0/(cmratio_ti*rho(i,1,1) + e_a*dns_a(i,1,1) )
        rhoinv = 1.0d0/(cmratio_ti*rho(i,1,1) )
        vpir(i,1,1) =(0.50d0*vpir(i,1,1) + c_ar(i,1,1) )*rhoinv
        vpiz(i,1,1) =(0.50d0*vpiz(i,1,1) + c_az(i,1,1) )*rhoinv
        vpiphi(i,1,1) =(0.50d0*vpiphi(i,1,1) + c_aphi(i,1,1) )*rhoinv

        vdrftr(i,1,1) =(ur(i,1,1) + vpir(i,1,1) )*vac(i,1,1)
        vdrftz(i,1,1) =(uz(i,1,1) + vpiz(i,1,1) )*vac(i,1,1)
        vdrftphi(i,1,1) =(utotphi(i,1,1) + vpiphi(i,1,1) )*vac(i,1,1)

! parallel velocity
!        vdotb = ur(i,1,1)*ubr + uz(i,1,1)*ubz + utotphi(i,1,1)*ubphi
        vdotb = ur(i,1,1)*ubr + uz(i,1,1)*ubz + uphi(i,1,1)*ubphi
        vparar(i,1,1) = vdotb*ubr
        vparaz(i,1,1) = vdotb*ubz
        vparaphi(i,1,1) = vdotb*ubphi
      end do

      end if
! ------------- loop end for HM simulation---------------

      
      if(flag_HM)then
! rho_vr, rho_vz, rho_vphi gives density evolution
        call divergence(0,rho_vr,rho_vz,rho_vphi,drho)
! for dvr, dvz, dvphi, dprs: pressure gradient and compression term included
        call divergence4hm
        call gradient(0,vparar,vparar_r,vparar_z,vparar_phi)
        call gradient(0,vparaz,vparaz_r,vparaz_z,vparaz_phi)
        call gradient(0,vparaphi,vparaphi_r,vparaphi_z,vparaphi_phi)

      else
! vr, vz, vphi are momentum
        call divergence(0,vr,vz,vtotphi,drho)
! for dvr, dvz, dvphi, dprs: pressure gradient and compression term included
        call divergence4
      end if

! second derivatives of velocity for viscosity
      call rotation5(0,nu_spt,rho,omegar,omegaz,omegaphi,rot2_ur,rot2_uz,rot2_uphi)
      call gradient5(0,nu_spt,rho,div_u,div_ur,div_uz,div_uphi)

! laplacian of rho and prs
      call laplacian5(0,rho,rho0,lpl3_rho,nun_spt) !2013-02-27
      call laplacian5(0,prs,prs0,lpl3_prs,chi_spt) !2013-02-17


! for induction equation
!      call rotation(0,er,ez,ephi,dbr,dbz,dbphi)
      call rotation(0,er_res,ez_res,ephi_res,dbr,dbz,dbphi) !2019-06-04

      de_dissp = 0.0d0

! ------------- loop of mhd equations ---------------

#ifdef AMDGPU      
!$omp  target teams distribute parallel do
!$omp& private(k,i,dvr1,dvz1,dvphi1,drho1,dprs1,dbr1,dbz1,dbphi1) collapse(2)
#else
!$omp  parallel do
!$omp& private(k,i,dvr1,dvz1,dvphi1,drho1,dprs1,dbr1,dbz1,dbphi1) collapse(2)
#endif
      do k = lphistart, lphiend
      do i = 1, lrz
! momentum equation
! @v/@t = -dvi(v*u)- grad(prs) + (j x B) - rot(nu rho rot(u)) + grad (nu rho div(u))

      if(flag_HM)then
      dvr1 =(- dvr(i,1,k) + vphi(i,1,k)*vdrftphi(i,1,k)/grr(i,1,k) & !2014-07-18
                   + (vpir(i,1,k)*vparar_r(i,1,k) + vpiz(i,1,k)*vparar_z(i,1,k) &
                     +vpiphi(i,1,k)*vparar_phi(i,1,k) - vpiphi(i,1,k)*vparaphi(i,1,k)/grr(i,1,k) &
                     )*rho(i,1,k) &
                   + (cnetphi(i,1,k)*bz(i,1,k) - cnetz(i,1,k)*bphi(i,1,k)) &
                   + (div_ur(i,1,k)*c43 - rot2_ur(i,1,k) ) &
                  )*dt - dvr0(i,1,k)

      dvz1 =(- dvz(i,1,k) &
                   + (vpir(i,1,k)*vparaz_r(i,1,k) + vpiz(i,1,k)*vparaz_z(i,1,k) &
                     +vpiphi(i,1,k)*vparaz_phi(i,1,k) &
                     )*rho(i,1,k) &
                   + (cnetr(i,1,k)*bphi(i,1,k) - cnetphi(i,1,k)*br(i,1,k)) &
                   + (div_uz(i,1,k)*c43 - rot2_uz(i,1,k) ) &
                  )*dt - dvz0(i,1,k)

      dvphi1 =(- dvphi(i,1,k) - vr(i,1,k)*vdrftphi(i,1,k)/grr(i,1,k) &
                     + (vpir(i,1,k)*vparaphi_r(i,1,k) + vpiz(i,1,k)*vparaphi_z(i,1,k) &
                       +vpiphi(i,1,k)*vparaphi_phi(i,1,k) + vpiphi(i,1,k)*vparar(i,1,k)/grr(i,1,k) &
                       )*rho(i,1,k) &
                     + (cnetz(i,1,k)*br(i,1,k) - cnetr(i,1,k)*bz(i,1,k)) &
                     + (div_uphi(i,1,k)*c43 - rot2_uphi(i,1,k) ) &
                    )*dt - dvphi0(i,1,k)

      else      
!      dvr1 =(- dvr(i,1,k) + vtotphi(i,1,k)*utotphi(i,1,k)/grr(i,1,k) &
      dvr1 =(- dvr(i,1,k) + vphi(i,1,k)*utotphi(i,1,k)/grr(i,1,k) & !2025-11-09
                   + (cnetphi(i,1,k)*bz(i,1,k) - cnetz(i,1,k)*bphi(i,1,k)) &
                   + (div_ur(i,1,k)*c43 - rot2_ur(i,1,k) ) &
                  )*dt - dvr0(i,1,k)

      dvz1 =(- dvz(i,1,k) &
                   + (cnetr(i,1,k)*bphi(i,1,k) - cnetphi(i,1,k)*br(i,1,k)) &
                   + (div_uz(i,1,k)*c43 - rot2_uz(i,1,k) ) &
                  )*dt - dvz0(i,1,k)

      dvphi1 =(- dvphi(i,1,k) - vr(i,1,k)*utotphi(i,1,k)/grr(i,1,k) &
                     + (cnetz(i,1,k)*br(i,1,k) - cnetr(i,1,k)*bz(i,1,k)) &
                     + (div_uphi(i,1,k)*c43 - rot2_uphi(i,1,k) ) &
                    )*dt - dvphi0(i,1,k)
      end if
! continuity equation
! @rho/@t = - div (rho v)
      drho1 = (-drho(i,1,k) + lpl3_rho(i,1,k) )*dt - drho0(i,1,k)

! @p/@t = - div(p v) - (gamma -1)*p div(v)
      dprs1 =( &
                  - dprs(i,1,k) &
                  + (gamma-1.0d0)*nu_spt(i,1,k)*rho(i,1,k) &
                   *(omegar(i,1,k)**2 + omegaz(i,1,k)**2 + omegaphi(i,1,k)**2 &
                    +div_u(i,1,k)**2*c43 &
                    ) &
                  + (gamma-1.0d0)*eta_spt(i,1,k) &
                   *(cr(i,1,k)*(cr(i,1,k) - cr0(i,1,k) ) &
                    +cz(i,1,k)*(cz(i,1,k) - cz0(i,1,k) ) &
                    +cphi(i,1,k)*(cphi(i,1,k) - cphi0(i,1,k) ) &
                    ) &
                  + lpl3_prs(i,1,k) &
                   )*dt &
                  - dprs0(i,1,k)
!      denrg1 = (-denrg(i,1,k) + lpl3_prs(i,1,k)/(gamma-1.0d0) + etransfer(i,1,k,1) + etransfer(i,1,k,2) )*dt &
!             - denrg0(i,1,k) !2021-05-05

! induction equation
! @b/@t = -rot e, dbr0 etc are considered by er = er - er0 etc
      dbr1 = -dbr(i,1,k)*dt
      dbz1 = -dbz(i,1,k)*dt
      dbphi1 = -dbphi(i,1,k)*dt

! time integration
      if(istep.eq.1)then
          br2(i,1,k) = br(i,1,k)
          bz2(i,1,k) = bz(i,1,k)
          bphi2(i,1,k) = bphi(i,1,k)
          vr2(i,1,k) = vr(i,1,k)
          vz2(i,1,k) = vz(i,1,k)
          vphi2(i,1,k) = vphi(i,1,k)
          rho2(i,1,k) = rho(i,1,k)
          prs2(i,1,k) = prs(i,1,k)
!          enrg2(i,1,k) = enrg(i,1,k)

          vr1(i,1,k) = vr(i,1,k) + c0*dvr1
          vr(i,1,k) =(vr(i,1,k) + c1*dvr1)*vac(i,1,k)

          vz1(i,1,k) = vz(i,1,k) + c0*dvz1
          vz(i,1,k) =(vz(i,1,k) + c1*dvz1)*vac(i,1,k)

          vphi1(i,1,k) = vphi(i,1,k) + c0*dvphi1
          vphi(i,1,k) =(vphi(i,1,k) + c1*dvphi1)*vac(i,1,k)

          rho1(i,1,k) = rho(i,1,k) + c0*drho1
          rho(i,1,k) = rho(i,1,k) + c1*drho1

          prs1(i,1,k) = prs(i,1,k) + c0*dprs1
          prs(i,1,k) = prs(i,1,k) + c1*dprs1
!          enrg1(i,1,k) = enrg(i,1,k) + c0*denrg1
!          enrg(i,1,k)  = enrg(i,1,k) + c1*denrg1

          br1(i,1,k) = br(i,1,k) + c0*dbr1
          br(i,1,k) = br(i,1,k) + c1*dbr1

          bz1(i,1,k) = bz(i,1,k) + c0*dbz1
          bz(i,1,k) = bz(i,1,k) + c1*dbz1

          bphi1(i,1,k) = bphi(i,1,k) + c0*dbphi1
          bphi(i,1,k) = bphi(i,1,k) + c1*dbphi1

      else if(istep.eq.2.or.istep.eq.3)then

          vr1(i,1,k) = vr1(i,1,k) + c0*dvr1
          vr(i,1,k) =(vr2(i,1,k) + c1*dvr1)*vac(i,1,k)

          vz1(i,1,k) = vz1(i,1,k) + c0*dvz1
          vz(i,1,k) =(vz2(i,1,k) + c1*dvz1)*vac(i,1,k)

          vphi1(i,1,k) = vphi1(i,1,k) + c0*dvphi1
          vphi(i,1,k) =(vphi2(i,1,k) + c1*dvphi1)*vac(i,1,k)

          rho1(i,1,k) = rho1(i,1,k) + c0*drho1
          rho(i,1,k) = rho2(i,1,k) + c1*drho1

          prs1(i,1,k) = prs1(i,1,k) + c0*dprs1
          prs(i,1,k) = prs2(i,1,k) + c1*dprs1
!          enrg1(i,1,k) = enrg1(i,1,k) + c0*denrg1
!          enrg(i,1,k)  = enrg2(i,1,k) + c1*denrg1

          br1(i,1,k) = br1(i,1,k) + c0*dbr1
          br(i,1,k) = br2(i,1,k) + c1*dbr1

          bz1(i,1,k) = bz1(i,1,k) + c0*dbz1
          bz(i,1,k) = bz2(i,1,k) + c1*dbz1

          bphi1(i,1,k) = bphi1(i,1,k) + c0*dbphi1
          bphi(i,1,k) = bphi2(i,1,k) + c1*dbphi1

      else if(istep.eq.4)then

          vr(i,1,k) =(vr1(i,1,k) + c0*dvr1)*vac(i,1,k)

          vz(i,1,k) =(vz1(i,1,k) + c0*dvz1)*vac(i,1,k)

          vphi(i,1,k) =(vphi1(i,1,k) + c0*dvphi1)*vac(i,1,k)

          rho(i,1,k) = rho1(i,1,k) + c0*drho1

          prs(i,1,k) = prs1(i,1,k) + c0*dprs1
!          enrg(i,1,k) = enrg1(i,1,k) + c0*denrg1

          br(i,1,k) = br1(i,1,k) + c0*dbr1

          bz(i,1,k) = bz1(i,1,k) + c0*dbz1

          bphi(i,1,k) = bphi1(i,1,k) + c0*dbphi1

      else if(istep.eq.0)then

          dvr0(i,1,k) = dvr1

          dvz0(i,1,k) = dvz1

          dvphi0(i,1,k) = dvphi1

          drho0(i,1,k) = drho1

          dprs0(i,1,k) = dprs1
!          denrg0(i,1,k) = denrg1

      end if

      end do
      end do

! ------------- loop end of mhd equations ---------------


      if(istep.ne.0)then
        e_trans  = e_trans + c0*de_trans
        e_dissp  = e_dissp + c0*de_dissp
      end if

      call mhd_boundary(flag_HM)

      call gradmg

end
!--------------------------------------------------------------------
subroutine e_field(isw)
! simple ohm's law
! modified on 2013-12-18, utotphi = uphi + utor
!--------------------------------------------------------------------
      use parameters
      use field
      use grid, only:major_r
      implicit none

      real(8)::ppara,ubr,ubz,ubphi,echarge2
      integer::i,isw

! current density

      call cc_density

! pressure gradient
        call gradient(1,prs,prs_r,prs_z,prs_phi)

! e = - v x b + eta j
      if(isw.eq.0)then

! store the equlibrium electric field and current density
#ifdef AMDGPU
!$omp target teams distribute parallel do private(i,ppara,ubr,ubz,ubphi,echarge2)
#else
!$omp parallel do private(i,ppara,ubr,ubz,ubphi,echarge2)
#endif
         do i = 1, lrzphi
          ubr = br(i,1,1)/babs(i,1,1)
          ubz = bz(i,1,1)/babs(i,1,1)
          ubphi = bphi(i,1,1)/babs(i,1,1)

!          ppara =(prs_r(i,1,1)*ubr + prs_z(i,1,1)*ubz + prs_phi(i,1,1)*ubphi)*vac(i,1,1)*0.50d0
          ppara = 0.0d0
          echarge2 = m_ti/(e_e*rho(i,1,1) )

! ur, uz, uphi are velocity
          ur(i,1,1) = vr(i,1,1)/rho(i,1,1)
          uz(i,1,1) = vz(i,1,1)/rho(i,1,1)
          uphi(i,1,1) = vphi(i,1,1)/rho(i,1,1)

          utotphi(i,1,1) = uphi(i,1,1) + utor(i,1,1)
          vtotphi(i,1,1) = utotphi(i,1,1)*rho(i,1,1)

          epara0(i,1,1) = ppara*echarge2

          er0(i,1,1) =-(uphi(i,1,1)*bz(i,1,1) - uz(i,1,1)*bphi(i,1,1)) &
                     + ppara*ubr*echarge2

          ez0(i,1,1) =-(ur(i,1,1)*bphi(i,1,1) - uphi(i,1,1)*br(i,1,1)) &
                     + ppara*ubz*echarge2

          ephi0(i,1,1) =-(uz(i,1,1)*br(i,1,1) - ur(i,1,1)*bz(i,1,1)) &
                       + ppara*ubphi*echarge2

          er0_res(i,1,1) = eta_spt(i,1,1)*cr(i,1,1)
          ez0_res(i,1,1) = eta_spt(i,1,1)*cz(i,1,1)
          ephi0_res(i,1,1) = eta_spt(i,1,1)*cphi(i,1,1)

          er(i,1,1) = 0.0d0
          ez(i,1,1) = 0.0d0
          ephi(i,1,1)= 0.0d0
          epara(i,1,1)= 0.0d0

          er_res(i,1,1) = 0.0d0
          ez_res(i,1,1) = 0.0d0
          ephi_res(i,1,1)= 0.0d0

          cr0(i,1,1) = cr(i,1,1)
          cz0(i,1,1) = cz(i,1,1)
          cphi0(i,1,1)= cphi(i,1,1)

        end do

      else

#ifdef AMDGPU
!$omp target teams distribute parallel do private(i,ppara,ubr,ubz,ubphi,echarge2)
#else
!$omp parallel do private(i,ppara,ubr,ubz,ubphi,echarge2)
#endif
         do i = 1, lrzphi
          ubr = br(i,1,1)/babs(i,1,1)
          ubz = bz(i,1,1)/babs(i,1,1)
          ubphi = bphi(i,1,1)/babs(i,1,1)

!          ppara =(prs_r(i,1,1)*ubr + prs_z(i,1,1)*ubz + prs_phi(i,1,1)*ubphi)*vac(i,1,1)*0.50d0
          ppara = 0.0d0
          echarge2 = m_ti/(e_e*rho(i,1,1) )

! ur, uz, uphi are velocity
          ur(i,1,1) = vr(i,1,1)/rho(i,1,1)
          uz(i,1,1) = vz(i,1,1)/rho(i,1,1)
          uphi(i,1,1) = vphi(i,1,1)/rho(i,1,1)

          utotphi(i,1,1) = uphi(i,1,1) + utor(i,1,1)
          vtotphi(i,1,1) = utotphi(i,1,1)*rho(i,1,1)

          epara(i,1,1) = ppara*echarge2 - epara0(i,1,1)

          er(i,1,1) =-(uphi(i,1,1)*bz(i,1,1) - uz(i,1,1)*bphi(i,1,1)) &
                    + ppara*ubr*echarge2 - er0(i,1,1)

          ez(i,1,1) =-(ur(i,1,1)*bphi(i,1,1) - uphi(i,1,1)*br(i,1,1)) &
                    + ppara*ubz*echarge2 - ez0(i,1,1)

          ephi(i,1,1) =-(uz(i,1,1)*br(i,1,1) - ur(i,1,1)*bz(i,1,1)) &
                      + ppara*ubphi*echarge2 - ephi0(i,1,1)

          er_res(i,1,1) = eta_spt(i,1,1)*cr(i,1,1) - er0_res(i,1,1) + er(i,1,1)
          ez_res(i,1,1) = eta_spt(i,1,1)*cz(i,1,1) - ez0_res(i,1,1) + ez(i,1,1)
          ephi_res(i,1,1) = eta_spt(i,1,1)*cphi(i,1,1) - ephi0_res(i,1,1) + ephi(i,1,1)
      end do

      end if


! velocity ur, uz, uphi and vorticity
      call rotation(0,ur,uz,utotphi,omegar,omegaz,omegaphi)
      call divergence(0,ur,uz,utotphi,div_u)
      call periodic_field_mlt7(er_res,ez_res,ephi_res,omegar,omegaz,omegaphi,div_u) !2025-10-05

end
!--------------------------------------------------------------------
subroutine write_check
!--------------------------------------------------------------------
      use mpiset
      use grid
      use equi_sol, only:raxis
      use field
      use particle !2013-05-02
      use check
      implicit none
!2025-03-05s
      real(8)::e_kin,e_mag,e_thr,e_ep(0:l_kin_species)
      real(8)::e_kin_total,e_mag_total,e_thr_total,e_ep_total(0:l_kin_species),e_total
!2025-03-05e
      real(8)::e_trans_total(0:l_kin_species),e_dissp_total !2025-03-05
      real(8)::wat
      real(8)::srho_vr(lr,lz,lphi),srho_vz(lr,lz,lphi),srho_vphi(lr&
           &,lz,lphi)
      real(8)::delta_br(lr,lz,lphi),delta_bz(lr,lz,lphi)&
           &,delta_bphi(lr,lz,lphi)
      real(8)::delta_prs(lr,lz,lphi)
      real(8)::srho_vr_n(lphi_n_size,lr,lz),srho_vz_n(lphi_n_size,lr&
           &,lz)
      real(8)::srho_vphi_n(lphi_n_size,lr,lz)
      real(8)::delta_br_n(lphi_n_size,lr,lz),delta_bz_n(lphi_n_size&
           &,lr,lz)
      real(8)::delta_bphi_n(lphi_n_size,lr,lz)&
           &,delta_prs_n(lphi_n_size,lr,lz)
      real(8)::energy_n(0:lphi_n),energy_n_total(0:lphi_n)
      integer::i,j,k,n,n1,n2
      integer::l_kin_mpi=l_kin_species+1 !2025-03-05
      integer::ititle=0
      
!      call mpi_reduce(loss1,loss1_total,1,mpi_integer !             
      !        ,mpi_sum,0,mpi_comm_world,mpi_err)
!      call mpi_reduce(loss2,loss2_total,1,mpi_integer !             
      !        ,mpi_sum,0,mpi_comm_world,mpi_err)


      if(my_rank.eq.0)then
        if(ititle.eq.0)then
          write(7,'(//)')
        end if
        write(7,600)'kstep=',kstep,'  t=',t !2017-09-16
      end if
600   format(a,i8,a,1pe14.5,a,2i8)


!--------------------------------------------------------------------
! check of energy evolution

      do i = 1, lrzphi
        delta_br(i,1,1) = br(i,1,1) - br0(i,1,1)
        delta_bz(i,1,1) = bz(i,1,1) - bz0(i,1,1)
        delta_bphi(i,1,1)=bphi(i,1,1)-bphi0(i,1,1)
        delta_prs(i,1,1)= prs(i,1,1)- prs0(i,1,1)
      end do


      e_kin = 0.0d0
      e_mag = 0.0d0
      e_thr = 0.0d0
      e_ep  = 0.0d0

      do k = lphistart, lphiend
      do j = lzstart, lzend
      do i = lrstart, lrend
        e_kin  = e_kin + 0.50d0*rho(i,j,k)*grr(i,j,k) *(ur(i,j,k)**2 &
             &+ uz(i,j,k)**2 + uphi(i,j,k)**2 )
        e_mag  = e_mag + 0.50d0*grr(i,j,k) *( delta_br(i,j,k)*(br(i,j&
             &,k)+br0(i,j,k)) +delta_bz(i,j,k)*(bz(i,j,k)+bz0(i,j,k))&
             & +delta_bphi(i,j,k)*(bphi(i,j,k)+bphi0(i,j,k)) )
        e_thr = e_thr + delta_prs(i,j,k)/(gamma-1.0d0)*grr(i,j,k)
        e_ep(1) = e_ep(1) + 0.50d0*(      (ppara_i(i,j,k)-ppara_i0(i,j,k) )&
             & +2.0d0*(pperp_i(i,j,k)-pperp_i0(i,j,k) ) )*grr(i,j,k)
        e_ep(2) = e_ep(2) + 0.50d0*(      (ppara_a(i,j,k)-ppara_a0(i,j,k) )&
             & +2.0d0*(pperp_a(i,j,k)-pperp_a0(i,j,k) ) )*grr(i,j,k)
      end do
      end do
      end do

      e_kin = e_kin*dr*dz*dphi
      e_mag = e_mag*dr*dz*dphi
      e_thr = e_thr*dr*dz*dphi
      e_ep  = e_ep *dr*dz*dphi

      call mpi_reduce(e_kin,e_kin_total,1,mpi_real8 ,mpi_sum,0&
           &,mhd_world,mpi_err)
      call mpi_reduce(e_mag,e_mag_total,1,mpi_real8 ,mpi_sum,0&
           &,mhd_world,mpi_err)
      call mpi_reduce(e_thr,e_thr_total,1,mpi_real8 ,mpi_sum,0&
           &,mhd_world,mpi_err)
      call mpi_reduce(e_ep,e_ep_total,l_kin_mpi,mpi_real8 ,mpi_sum,0& !2025-03-05
           &,mhd_world,mpi_err)

      if(ititle.eq.0)then
        ititle = 1
        if(my_rank.eq.0)then
#ifdef EXPEQ
          write(40,'(11(x,a))')'kstep','t[ms]','kin_energy','mag_energy'&
#else
               write(40,'(11(x,a))')'kstep','t[R/vA]','kin_energy','mag_energy'&
#endif
               ,'thr_energy' ,'beam_energy','beam_energy_trns','ep_energy','ep_energy_trns' &
               ,'dssp_energy','total_energy' 
#ifdef EXPEQ
          write(41,'(43(x,a))')'kstep','t[ms]','e_n=0','e_n=1','e_n=2'&
#else
          write(41,'(43(x,a))')'kstep','t[R/vA]','e_n=0','e_n=1','e_n=2'&
#endif
                ,'e_n=3','e_n=4','e_n=5','e_n=6','e_n=7','e_n=8' &
                ,'e_n=9','e_n=10','e_n=11','e_n=12','e_n=13','e_n=14' &
                ,'e_n=15','e_n=16','e_n=17','e_n=18','e_n=19','e_n=20' &
                ,'e_n=21','e_n=22','e_n=23','e_n=24','e_n=25','e_n=26' &
                ,'e_n=27','e_n=28','e_n=29','e_n=30','e_n=31','e_n=32' &
                ,'e_n=33','e_n=34','e_n=35','e_n=36','e_n=37','e_n=38' &
                ,'e_n=39','e_n=40'
        end if
      end if

      call mpi_reduce(e_trans,e_trans_total,l_kin_mpi,mpi_real8 ,mpi_sum,0 & !2025-03-05
           &,mhd_world,mpi_err)
      call mpi_reduce(e_dissp,e_dissp_total,1,mpi_real8 ,mpi_sum,0&
           &,mhd_world,mpi_err)

      e_total = e_kin_total + e_mag_total + e_thr_total + e_trans_total(1) + e_trans_total(2) &
              + e_dissp_total
      if(my_rank.eq.0)then
#ifdef EXPEQ
        wat = t/omega_a*1.0d3
#else
        wat = t/raxis
#endif
        write(40,'(i10,10(1pe14.5))')kstep,wat,e_kin_total,e_mag_total,e_thr_total &
            ,e_ep_total(1),e_trans_total(1),e_ep_total(2),e_trans_total(2) &
            ,e_dissp_total,e_total
        flush(40)
      end if

!--------------------------------------------------------------------
! check energy of each toroidal mode number

      do i = 1, lrzphi
        srho_vr(i,1,1) = sqrt(rho(i,1,1))*ur(i,1,1)
        srho_vz(i,1,1) = sqrt(rho(i,1,1))*uz(i,1,1)
        srho_vphi(i,1,1) = sqrt(rho(i,1,1))*uphi(i,1,1)
      end do

      call toroidal_mode_d(srho_vr,srho_vr_n)
      call toroidal_mode_d(srho_vz,srho_vz_n)
      call toroidal_mode_d(srho_vphi,srho_vphi_n)
      call toroidal_mode_d(delta_br,delta_br_n)
      call toroidal_mode_d(delta_bz,delta_bz_n)
      call toroidal_mode_d(delta_bphi,delta_bphi_n)
      call toroidal_mode_d(delta_prs,delta_prs_n)

      do n = 0, lphi_n
        energy_n(n) = 0.0d0
      end do

      n = 0
      n1 = 2*n + 1
      do j = lzstart, lzend
      do i = lrstart, lrend
        energy_n(n) = energy_n(n) + 0.50d0*grr(i,j,1) *(srho_vr_n(n1&
             &,i,j)**2 + srho_vz_n(n1,i,j)**2 +srho_vphi_n(n1,i,j)**2&
             & ) + 0.50d0*grr(i,j,1) *( delta_br_n(n1,i,j)&
             &*(delta_br_n(n1,i,j)+2.0d0*br0(i,j,1)) +delta_bz_n(n1,i&
             &,j)*(delta_bz_n(n1,i,j)+2.0d0*bz0(i,j,1)) &
             &+delta_bphi_n(n1,i,j)*(delta_bphi_n(n1,i,j)+2.0d0&
             &*bphi0(i,j,1)) ) + delta_prs_n(n1,i,j)/(gamma-1.0d0)&
             &*grr(i,j,1)
      end do
      end do

      do n = 1, lphi_n
        n1 = 2*n + 1
        n2 = 2*n + 2
      do j = lzstart, lzend
      do i = lrstart, lrend
        energy_n(n) = energy_n(n) + 0.50d0*grr(i,j,1) *( srho_vr_n(n1&
             &,i,j)**2 +srho_vz_n(n1,i,j)**2 +srho_vphi_n(n1,i,j)**2&
             & ) + 0.50d0*grr(i,j,1) *( delta_br_n(n1,i,j)**2 &
             &+delta_bz_n(n1,i,j)**2 +delta_bphi_n(n1,i,j)**2 ) +&
             & 0.50d0*grr(i,j,1) *( srho_vr_n(n2,i,j)**2 &
             &+srho_vz_n(n2,i,j)**2 +srho_vphi_n(n2,i,j)**2 ) +&
             & 0.50d0*grr(i,j,1) *( delta_br_n(n2,i,j)**2 &
             &+delta_bz_n(n2,i,j)**2 +delta_bphi_n(n2,i,j)**2 ) 
      end do
      end do
      end do

      energy_n(0) = energy_n(0)*dr*dz*phileng/dble(mpi_proc_phi)
      do n = 1, lphi_n
        energy_n(n) = energy_n(n)*dr*dz*phileng/dble(mpi_proc_phi)&
             &*0.50d0
      end do

      n = lphi_n + 1
!2012-06-17
      call mpi_reduce(energy_n(0),energy_n_total(0),n,mpi_real8 &
           &,mpi_sum,0,mhd_world,mpi_err) 
!2012-06-17 end

      if(my_rank.eq.0)then
#ifdef EXPEQ
        wat = t/omega_a*1.0d3
#else
        wat = t/raxis !2017-09-16
#endif
        write(41,'(i10,42(1pe14.5))')kstep,wat,energy_n_total
        flush(41)
      end if

end
!--------------------------------------------------------------------
#else 
!subroutines mhd_kti, e_field_kti, and write_check_kti
!--------------------------------------------------------------------
subroutine mhd_kti(istep,flag_HM)
! 2023-09-09: subroutines mhd and hm are unified
! 2025-03-23: electron pressure evolution is solved
!--------------------------------------------------------------------
      use parameters
      use mpiset
      use field, only:br,bz,bphi,babs,er,ez,ephi,cr,cz,cphi,cr0,cz0,cphi0 &
                     ,rho,prs,vr,vz,vphi,vdrftr,vdrftz,vdrftphi &
                     ,prs_r,prs_z,prs_phi,dns_e & !2015-10-11
                     ,ur,uz,uphi,div_u,omegar,omegaz,omegaphi &
                     ,utotphi,vtotphi,utor & !2017-05-05
                     ,temp_i0 & !2016-04-04
                     ,dns_a,mom_a,dns_i,mom_i & !2017-05-04
                     ,ppara_a,pperp_a,ppara_i,pperp_i &
                     ,gradbr,gradbz,gradbphi,curvbr,curvbz,curvbphi &
                     ,br0,bz0,bphi0,babs0,rho0,prs0 &
                     ,dbr,dbz,dbphi,dvr,dvz,dvphi,drho,dprs &
                     ,br1,bz1,bphi1,vr1,vz1,vphi1,rho1,prs1 &
                     ,br2,bz2,bphi2,vr2,vz2,vphi2,rho2,prs2 &
                     ,dvr0,dvz0,dvphi0,drho0,dprs0,vac &
                     ,nu_spt,eta_spt,nu_n,chi,gamma,chi_spt,nun_spt &
                     ,er_res,ez_res,ephi_res !2019-06-04
      use grid, only:dt,grr,dr,dz,dphi
      use check
      implicit none

      logical::flag_HM !2013-01-12
      integer::istep
      integer::i,j,k
      real(8)::de_trans_local(lkinetic),de_dissp_local
      real(8)::rot2_ur(lr,lz,lphi),rot2_uz(lr,lz,lphi),rot2_uphi(lr,lz,lphi)
      real(8)::div_ur(lr,lz,lphi),div_uz(lr,lz,lphi),div_uphi(lr,lz,lphi)
      real(8)::pur(lr,lz,lphi),puz(lr,lz,lphi),puphi(lr,lz,lphi) !2025-03-24
      real(8)::lpl3_prs(lr,lz,lphi) !2025-03-26
      real(8)::cnetr(lr,lz,lphi),cnetz(lr,lz,lphi),cnetphi(lr,lz,lphi)
      real(8)::c_ir(lr,lz,lphi),c_iz(lr,lz,lphi),c_iphi(lr,lz,lphi)
      real(8)::c_ar(lr,lz,lphi),c_az(lr,lz,lphi),c_aphi(lr,lz,lphi)
!      real(8)::magc_e_r(lr,lz,lphi),magc_e_z(lr,lz,lphi),magc_e_phi(lr,lz,lphi)
      real(8)::magc_i_r(lr,lz,lphi),magc_i_z(lr,lz,lphi),magc_i_phi(lr,lz,lphi)
      real(8)::magc_a_r(lr,lz,lphi),magc_a_z(lr,lz,lphi),magc_a_phi(lr,lz,lphi)
      real(8)::dvr1,dvz1,dvphi1,dbr1,dbz1,dbphi1,dprs1 !2025-03-23
      real(8)::b21,ppara_i1,ppara_a1 !2016-03-29
!      real(8)::eidealr,eidealz,eidealphi
      real(8)::eid_r(lr,lz,lphi),eid_z(lr,lz,lphi),eid_phi(lr,lz,lphi)
      real(8)::dbid_r(lr,lz,lphi),dbid_z(lr,lz,lphi),dbid_phi(lr,lz,lphi)
      real(8)::c0,c1
      real(8)::c43=4.0d0/3.0d0
      real(8)::mcratio_i=m_i/e_i, mcratio_a=m_a/e_a !2017-05-04
      real(8)::ubr,ubz,ubphi,vpir,vpiz,vpiphi,mom_sum
#ifdef AMDGPU
      real(8)::t1,t2,t3,t4,t5
#endif

! time integration

      if(istep.eq.1)then
        c0 = 1.0d0/6.0d0
        c1 = 0.50d0
      else if(istep.eq.2)then
        c0 = 1.0d0/3.0d0
        c1 = 0.50d0
      else if(istep.eq.3)then
        c0 = 1.0d0/3.0d0
        c1 = 1.0d0
      else if(istep.eq.4)then
        c0 = 1.0d0/6.0d0
      end if

! ------------- preparation ---------------

! magnetization current
      call gradient(1,pperp_i,c_ir,c_iz,c_iphi)
      call gradient(1,pperp_a,c_ar,c_az,c_aphi)

#ifdef AMDGPU      
!$omp target teams distribute parallel do private(i,b21,ppara_i1,ppara_a1)
#else
!$omp parallel do private(i,b21,ppara_i1,ppara_a1)
#endif
      do i = 1, lrzphi
        b21 = 1.0d0/babs(i,1,1)**2
        magc_i_r(i,1,1) =-(c_iphi(i,1,1)*bz(i,1,1) - c_iz(i,1,1)*bphi(i,1,1) )*b21
        magc_i_z(i,1,1) =-(c_ir(i,1,1)*bphi(i,1,1) - c_iphi(i,1,1)*br(i,1,1) )*b21
        magc_i_phi(i,1,1) =-(c_iz(i,1,1)*br(i,1,1) - c_ir(i,1,1)*bz(i,1,1) )*b21
        ppara_i1 =(ppara_i(i,1,1) - pperp_i(i,1,1) )/babs(i,1,1)

        c_ir(i,1,1) = ppara_i1*curvbr(i,1,1) &
                    + magc_i_r(i,1,1)

        c_iz(i,1,1) = ppara_i1*curvbz(i,1,1) &
                    + magc_i_z(i,1,1)

        c_iphi(i,1,1) = ppara_i1*curvbphi(i,1,1) &
                      + magc_i_phi(i,1,1)

        magc_a_r(i,1,1) =-(c_aphi(i,1,1)*bz(i,1,1) - c_az(i,1,1)*bphi(i,1,1) )*b21
        magc_a_z(i,1,1) =-(c_ar(i,1,1)*bphi(i,1,1) - c_aphi(i,1,1)*br(i,1,1) )*b21
        magc_a_phi(i,1,1) =-(c_az(i,1,1)*br(i,1,1) - c_ar(i,1,1)*bz(i,1,1) )*b21
        ppara_a1 =(ppara_a(i,1,1) - pperp_a(i,1,1) )/babs(i,1,1)
        c_ar(i,1,1) = ppara_a1*curvbr(i,1,1) &
                    + magc_a_r(i,1,1)
        c_az(i,1,1) = ppara_a1*curvbz(i,1,1) &
                    + magc_a_z(i,1,1)
        c_aphi(i,1,1) = ppara_a1*curvbphi(i,1,1) &
                      + magc_a_phi(i,1,1)

        cnetr(i,1,1) = cr(i,1,1) - c_ir(i,1,1) - c_ar(i,1,1)
        cnetz(i,1,1) = cz(i,1,1) - c_iz(i,1,1) - c_az(i,1,1)
        cnetphi(i,1,1) = cphi(i,1,1) - c_iphi(i,1,1) - c_aphi(i,1,1)

! for induction equation
        eid_r(i,1,1)   =-(uphi(i,1,1)*bz(i,1,1) - uz(i,1,1)*bphi(i,1,1))
        eid_z(i,1,1)   =-(ur(i,1,1)*bphi(i,1,1) - uphi(i,1,1)*br(i,1,1))
        eid_phi(i,1,1) =-(uz(i,1,1)*br(i,1,1)   - ur(i,1,1)*bz(i,1,1))
      end do

! laplacian of rho and prs
!      call laplacian5(0,rho,rho0,lpl3_rho,nun_spt) !2025-03-26
      call laplacian5(0,prs,prs0,lpl3_prs,chi_spt) !2025-03-26

      call rotation(0,eid_r,eid_z,eid_phi,dbid_r,dbid_z,dbid_phi)
!      call rotation(0,er,ez,ephi,dbr,dbz,dbphi)
      call rotation(0,er_res,ez_res,ephi_res,dbr,dbz,dbphi) !2019-06-04

! energy transfer from field to particles
#ifdef AMDGPU
      t1 = 0.0d0; t2 = 0.0d0; t3 = 0.0d0; t4 = 0.0d0; t5 = 0.0d0
#else
      de_trans_local = 0.0d0
#endif

#ifdef AMDGPU
!$omp target teams distribute parallel do private(k,j,i) collapse(3) reduction(+:t1,t2,t3,t4,t5)
#else
!$omp parallel do private(k,j,i) collapse(3) reduction(+:de_trans_local)
#endif
      do k = lphistart, lphiend
      do j = lzstart, lzend
      do i = lrstart, lrend
#ifdef AMDGPU
         t1 = t1 &
       + grr(i,j,k)*( &
            c_ar(i,j,k)*eid_r(i,j,k) &
          + c_az(i,j,k)*eid_z(i,j,k) &
          + c_aphi(i,j,k)*eid_phi(i,j,k) &
                    )

        t2 = t2 &
       + grr(i,j,k)*( &
            c_ir(i,j,k)*eid_r(i,j,k) &
          + c_iz(i,j,k)*eid_z(i,j,k) &
          + c_iphi(i,j,k)*eid_phi(i,j,k) &
                    )

! for electron pressure
        t3 = t3 &
       + grr(i,j,k)*( &
            ur(i,j,k)*prs_r(i,j,k) &
          + uz(i,j,k)*prs_z(i,j,k) &
          + uphi(i,j,k)*prs_phi(i,j,k) &
                    )

! for magnetic energy
        t4 = t4 &
       - grr(i,j,k)*( & !sign is corrected, the array dbr is actually -dbr
            br(i,j,k)*dbr(i,j,k) &
          + bz(i,j,k)*dbz(i,j,k) &
          + bphi(i,j,k)*dbphi(i,j,k) &
                    )
        t5 = t5 &
       - grr(i,j,k)*( & !sign is corrected, the array dbid_r is actually -dbid_r
            br(i,j,k)*dbid_r(i,j,k) &
          + bz(i,j,k)*dbid_z(i,j,k) &
          + bphi(i,j,k)*dbid_phi(i,j,k) &
                    )
#else
         de_trans_local(1) = de_trans_local(1) &
       + grr(i,j,k)*( &
            c_ar(i,j,k)*eid_r(i,j,k) &
          + c_az(i,j,k)*eid_z(i,j,k) &
          + c_aphi(i,j,k)*eid_phi(i,j,k) &
                    )

        de_trans_local(2) = de_trans_local(2) &
       + grr(i,j,k)*( &
            c_ir(i,j,k)*eid_r(i,j,k) &
          + c_iz(i,j,k)*eid_z(i,j,k) &
          + c_iphi(i,j,k)*eid_phi(i,j,k) &
                    )

! for electron pressure
        de_trans_local(3) = de_trans_local(3) &
       + grr(i,j,k)*( &
            ur(i,j,k)*prs_r(i,j,k) &
          + uz(i,j,k)*prs_z(i,j,k) &
          + uphi(i,j,k)*prs_phi(i,j,k) &
                    )

! for magnetic energy
        de_trans_local(4) = de_trans_local(4) &
       - grr(i,j,k)*( & !sign is corrected, the array dbr is actually -dbr
            br(i,j,k)*dbr(i,j,k) &
          + bz(i,j,k)*dbz(i,j,k) &
          + bphi(i,j,k)*dbphi(i,j,k) &
                    )
        de_trans_local(5) = de_trans_local(5) &
       - grr(i,j,k)*( & !sign is corrected, the array dbid_r is actually -dbid_r
            br(i,j,k)*dbid_r(i,j,k) &
          + bz(i,j,k)*dbid_z(i,j,k) &
          + bphi(i,j,k)*dbid_phi(i,j,k) &
                    )
#endif
      end do
      end do
      end do

#ifdef AMDGPU      
      de_trans(1) = t1*dr*dz*dphi*dt
      de_trans(2) = t2*dr*dz*dphi*dt
      de_trans(3) = t3*dr*dz*dphi*dt
      de_trans(4) = t4*dr*dz*dphi*dt
      de_trans(5) = t5*dr*dz*dphi*dt
#else
      de_trans = de_trans_local*dr*dz*dphi*dt
#endif

! ------------- preparation end ---------------

! pressure gradient is calculated in subroutine e_field_kti
!      call gradient(1,prs,prs_r,prs_z,prs_phi)

!2025-10-09
! second derivatives of velocity for viscosity, should be done before mom_sum is added to u
      call rotation5(0,nu_spt,rho,omegar,omegaz,omegaphi,rot2_ur,rot2_uz,rot2_uphi)
      call gradient5(0,nu_spt,rho,div_u,div_ur,div_uz,div_uphi)

! ------------- loop of preparation ---------------
      if(flag_HM)then
#ifdef AMDGPU
!$omp target teams distribute parallel do private(i,ubr,ubz,ubphi,vpir,vpiz,vpiphi,mom_sum)
#else
!$omp parallel do private(i,ubr,ubz,ubphi,vpir,vpiz,vpiphi,mom_sum)
#endif
      do i = 1, lrzphi
        ubr = br(i,1,1)/babs(i,1,1)
        ubz = bz(i,1,1)/babs(i,1,1)
        ubphi = bphi(i,1,1)/babs(i,1,1)

        vpir =(c_ir(i,1,1)*mcratio_i + c_ar(i,1,1)*mcratio_a)
        vpiz =(c_iz(i,1,1)*mcratio_i + c_az(i,1,1)*mcratio_a)
        vpiphi =(c_iphi(i,1,1)*mcratio_i + c_aphi(i,1,1)*mcratio_a)

        mom_sum =  mom_i(i,1,1) + mom_a(i,1,1)

        ur(i,1,1) =(ur(i,1,1) + mom_sum*ubr/rho(i,1,1) )*vac(i,1,1)
        uz(i,1,1) =(uz(i,1,1) + mom_sum*ubz/rho(i,1,1) )*vac(i,1,1)
        uphi(i,1,1) =(uphi(i,1,1) + mom_sum*ubphi/rho(i,1,1) )*vac(i,1,1)
        utotphi(i,1,1) =(utotphi(i,1,1) + mom_sum*ubphi/rho(i,1,1) )*vac(i,1,1)
        vtotphi(i,1,1) = utotphi(i,1,1)*rho(i,1,1)

        vdrftr(i,1,1) = ur(i,1,1) + vpir*vac(i,1,1)/rho(i,1,1)
        vdrftz(i,1,1) = uz(i,1,1) + vpiz*vac(i,1,1)/rho(i,1,1)
        vdrftphi(i,1,1)=utotphi(i,1,1) + vpiphi*vac(i,1,1)/rho(i,1,1)
      end do

      else

#ifdef AMDGPU
!$omp target teams distribute parallel do private(i,ubr,ubz,ubphi,mom_sum)
#else
!$omp parallel do private(i,ubr,ubz,ubphi,mom_sum)
#endif
      do i = 1, lrzphi
        ubr = br(i,1,1)/babs(i,1,1)
        ubz = bz(i,1,1)/babs(i,1,1)
        ubphi = bphi(i,1,1)/babs(i,1,1)

        mom_sum =  mom_i(i,1,1) + mom_a(i,1,1)

        ur(i,1,1) =(ur(i,1,1) + mom_sum*ubr/rho(i,1,1) )*vac(i,1,1)
        uz(i,1,1) =(uz(i,1,1) + mom_sum*ubz/rho(i,1,1) )*vac(i,1,1)
        uphi(i,1,1) =(uphi(i,1,1) + mom_sum*ubphi/rho(i,1,1) )*vac(i,1,1)
        utotphi(i,1,1) =(utotphi(i,1,1) + mom_sum*ubphi/rho(i,1,1) )*vac(i,1,1)
        vtotphi(i,1,1) = utotphi(i,1,1)*rho(i,1,1)

        vdrftr(i,1,1) = ur(i,1,1)
        vdrftz(i,1,1) = uz(i,1,1)
        vdrftphi(i,1,1)= utotphi(i,1,1)
      end do

      end if
! ------------- loop end ---------------

! for dvr, dvz, dvphi, dprs: pressure gradient and compression term included
      call divergence(0,ur,uz,utotphi,div_u) !used in divergence4hm
      call divergence4hm

! dissipated energy
      de_dissp_local = 0.0d0

#ifdef AMDGPU
!$omp target teams distribute parallel do private(k,j,i) collapse(3) reduction(+:de_dissp_local)
#else
!$omp parallel do private(k,j,i) collapse(3) reduction(+:de_dissp_local)
#endif
      do k = lphistart, lphiend
      do j = lzstart, lzend
      do i = lrstart, lrend
        de_dissp_local = de_dissp_local &
       + grr(i,j,k)*( &
                      nu_spt(i,j,k)*rho(i,j,k) &
                     *(omegar(i,j,k)**2 + omegaz(i,j,k)**2 + omegaphi(i,j,k)**2 &
                      +div_u(i,j,k)**2*4.0d0/3.0d0 &
                      ) &
!2019-06-22 commented out
!                     + (cr(i,j,k)*(er_res(i,j,k) - eidealr) &
!                       +cz(i,j,k)*(ez_res(i,j,k) - eidealz) &
!                       +cphi(i,j,k)*(ephi_res(i,j,k) - eidealphi) &
!                       ) &
!2019-06-22 commented out
!                    + eta_spt(i,j,k) &
!                     *(cr(i,j,k)*(cr(i,j,k) - cr0(i,j,k) ) &
!                      +cz(i,j,k)*(cz(i,j,k) - cz0(i,j,k) ) &
!                      +cphi(i,j,k)*(cphi(i,j,k) - cphi0(i,j,k) ) &
!                      ) &
!2015-10-11, include electric field in balance with pressure gradient
!                    +(cr(i,j,k)*prs_r(i,j,k) + cz(i,j,k)*prs_z(i,j,k) &
!                     +cphi(i,j,k)*prs_phi(i,j,k) )/(e_e*dns_e(i,j,k) ) &
                    )
      end do
      end do
      end do

      de_dissp = de_dissp_local*dr*dz*dphi*dt

! de_trans(4) = increment of real magnetic energy
! de_trans(5) = increment of ideal magnetic energy
! de_trans(4) - de_trans(5) < 0, de_trans(5) - de_trans(4) = dissipation
      de_dissp = de_dissp + de_trans(5) - de_trans(4)


! ------------- loop of mhd equations ---------------

#ifdef AMDGPU      
!$omp  target teams distribute parallel do
!$omp& private(k,i,dvr1,dvz1,dvphi1,dprs1,dbr1,dbz1,dbphi1) collapse(2)
#else
!$omp  parallel do
!$omp& private(k,i,dvr1,dvz1,dvphi1,dprs1,dbr1,dbz1,dbphi1) collapse(2)
#endif
      do k = lphistart, lphiend
      do i = 1, lrz
! momentum equation
! @v/@t = -dvi(v*u)- grad(prs) + (j x B) - rot(nu rho rot(u)) + grad (nu rho div(u))

      dvr1 =(- dvr(i,1,k) + vphi(i,1,k)*vdrftphi(i,1,k)/grr(i,1,k) &
                   + (cnetphi(i,1,k)*bz(i,1,k) - cnetz(i,1,k)*bphi(i,1,k)) &
                   + (div_ur(i,1,k)*c43 - rot2_ur(i,1,k) ) &
                  )*dt - dvr0(i,1,k)
      dvz1 =(- dvz(i,1,k) &
                   + (cnetr(i,1,k)*bphi(i,1,k) - cnetphi(i,1,k)*br(i,1,k)) &
                   + (div_uz(i,1,k)*c43 - rot2_uz(i,1,k) ) &
                  )*dt - dvz0(i,1,k)
      dvphi1 =(- dvphi(i,1,k) - vr(i,1,k)*vdrftphi(i,1,k)/grr(i,1,k) &
                     + (cnetz(i,1,k)*br(i,1,k) - cnetr(i,1,k)*bz(i,1,k)) &
                     + (div_uphi(i,1,k)*c43 - rot2_uphi(i,1,k) ) &
                    )*dt - dvphi0(i,1,k)

! @p/@t = - div(p v) - (gamma -1)*p div(v), calculated in subroutine divergence4hm
      dprs1 =( &
                  - dprs(i,1,k) &
!                  + (gamma-1.0d0)*nu_spt(i,1,k)*rho(i,1,k) &
!                   *(omegar(i,1,k)**2 + omegaz(i,1,k)**2 + omegaphi(i,1,k)**2 &
!                    +div_u(i,1,k)**2*c43 &
!                    ) &
!                  + (gamma-1.0d0)*eta_spt(i,1,k) &
!                   *(cnetr(i,1,k)*(cr(i,1,k) - cr0(i,1,k) ) &
!                    +cnetz(i,1,k)*(cz(i,1,k) - cz0(i,1,k) ) &
!                    +cnetphi(i,1,k)*(cphi(i,1,k) - cphi0(i,1,k) ) &
!                    ) &
                  + lpl3_prs(i,1,k) & !2025-03-26
                   )*dt &
                  - dprs0(i,1,k)

! induction equation
! @b/@t = -rot e, dbr0 etc are considered by er = er - er0 etc
      dbr1 = -dbr(i,1,k)*dt
      dbz1 = -dbz(i,1,k)*dt
      dbphi1 = -dbphi(i,1,k)*dt

! time integration
      if(istep.eq.1)then
          br2(i,1,k) = br(i,1,k)
          bz2(i,1,k) = bz(i,1,k)
          bphi2(i,1,k) = bphi(i,1,k)
          vr2(i,1,k) = vr(i,1,k)
          vz2(i,1,k) = vz(i,1,k)
          vphi2(i,1,k) = vphi(i,1,k)
          prs2(i,1,k) = prs(i,1,k) !2025-03-23

          vr1(i,1,k) = vr(i,1,k) + c0*dvr1
          vr(i,1,k) =(vr(i,1,k) + c1*dvr1)*vac(i,1,k)

          vz1(i,1,k) = vz(i,1,k) + c0*dvz1
          vz(i,1,k) =(vz(i,1,k) + c1*dvz1)*vac(i,1,k)

          vphi1(i,1,k) = vphi(i,1,k) + c0*dvphi1
          vphi(i,1,k) =(vphi(i,1,k) + c1*dvphi1)*vac(i,1,k)

          prs1(i,1,k) = prs(i,1,k) + c0*dprs1
          prs(i,1,k)  = prs(i,1,k) + c1*dprs1

          br1(i,1,k) = br(i,1,k) + c0*dbr1
          br(i,1,k) = br(i,1,k) + c1*dbr1

          bz1(i,1,k) = bz(i,1,k) + c0*dbz1
          bz(i,1,k) = bz(i,1,k) + c1*dbz1

          bphi1(i,1,k) = bphi(i,1,k) + c0*dbphi1
          bphi(i,1,k) = bphi(i,1,k) + c1*dbphi1

      else if(istep.eq.2.or.istep.eq.3)then
          vr1(i,1,k) = vr1(i,1,k) + c0*dvr1
          vr(i,1,k) =(vr2(i,1,k) + c1*dvr1)*vac(i,1,k)

          vz1(i,1,k) = vz1(i,1,k) + c0*dvz1
          vz(i,1,k) =(vz2(i,1,k) + c1*dvz1)*vac(i,1,k)

          vphi1(i,1,k) = vphi1(i,1,k) + c0*dvphi1
          vphi(i,1,k) =(vphi2(i,1,k) + c1*dvphi1)*vac(i,1,k)

          prs1(i,1,k) = prs1(i,1,k) + c0*dprs1
          prs(i,1,k)  = prs2(i,1,k) + c1*dprs1

          br1(i,1,k) = br1(i,1,k) + c0*dbr1
          br(i,1,k) = br2(i,1,k) + c1*dbr1

          bz1(i,1,k) = bz1(i,1,k) + c0*dbz1
          bz(i,1,k) = bz2(i,1,k) + c1*dbz1

          bphi1(i,1,k) = bphi1(i,1,k) + c0*dbphi1
          bphi(i,1,k) = bphi2(i,1,k) + c1*dbphi1

      else if(istep.eq.4)then
          vr(i,1,k) =(vr1(i,1,k) + c0*dvr1)*vac(i,1,k)
          vz(i,1,k) =(vz1(i,1,k) + c0*dvz1)*vac(i,1,k)
          vphi(i,1,k) =(vphi1(i,1,k) + c0*dvphi1)*vac(i,1,k)
          prs(i,1,k) = prs1(i,1,k) + c0*dprs1
          br(i,1,k) = br1(i,1,k) + c0*dbr1
          bz(i,1,k) = bz1(i,1,k) + c0*dbz1
          bphi(i,1,k) = bphi1(i,1,k) + c0*dbphi1

       else if(istep.eq.0)then
          dvr0(i,1,k) = dvr1
          dvz0(i,1,k) = dvz1
          dvphi0(i,1,k) = dvphi1
          dprs0(i,1,k) = dprs1

      end if

      end do
      end do

! ------------- loop end of mhd equations ---------------

      if(istep.ne.0)then
        e_trans  = e_trans + c0*de_trans
        e_dissp  = e_dissp + c0*de_dissp
      end if

      call mhd_boundary(flag_HM)

      call gradmg

end
!--------------------------------------------------------------------
subroutine e_field_kti(isw)
! simple ohm's law
! modified on 2021-07-25 for vac_pe for pressure gradient term
!--------------------------------------------------------------------
      use parameters
      use field
      use grid, only:major_r
      implicit none

      real(8)::vpara,ppara,ubr,ubz,ubphi,echarge1,echarge2
      real(8)::ppara_res !2021-08-30
      integer::i,isw

      echarge1 = 1.0/e_e

! current density

      call cc_density

! e = - v x b + eta j

! store the equlibrium electric field and current density
! rho, prs, and parallel momentum are given by ion partilces

#ifdef AMDGPU
!$omp target teams distribute parallel do private(i,ubr,ubz,ubphi,vpara)
#else
!$omp parallel do private(i,ubr,ubz,ubphi,vpara)
#endif
        do i = 1, lrzphi
          ubr = br(i,1,1)/babs(i,1,1)
          ubz = bz(i,1,1)/babs(i,1,1)
          ubphi = bphi(i,1,1)/babs(i,1,1)

          rho(i,1,1) = m_i*dns_i(i,1,1) + m_a*dns_a(i,1,1)
          vpara = vr(i,1,1)*ubr + vz(i,1,1)*ubz + vphi(i,1,1)*ubphi

! ur, uz, uphi are velocity, do not inculde parallel component
          vr(i,1,1) = vr(i,1,1) - vpara*ubr
          vz(i,1,1) = vz(i,1,1) - vpara*ubz
          vphi(i,1,1) = vphi(i,1,1) - vpara*ubphi

          ur(i,1,1) = vr(i,1,1)/rho(i,1,1)
          uz(i,1,1) = vz(i,1,1)/rho(i,1,1)
          uphi(i,1,1) = vphi(i,1,1)/rho(i,1,1)

          dns_e(i,1,1) = -(e_i*dns_i(i,1,1) + e_a*dns_a(i,1,1) )*echarge1
! prs is electron pressure
!          prs(i,1,1) = dns_e(i,1,1)*temp_e0(i,1,1)
        end do

! pressure gradient, used in subroutine mhd and hm
        call gradient(1,prs,prs_r,prs_z,prs_phi)

      if(isw.eq.0)then

! store the equlibrium electric field and current density
! rho, prs, and parallel momentum are given by ion partilces

#ifdef AMDGPU
!$omp target teams distribute parallel do private(i,ppara,ubr,ubz,ubphi,echarge2,ppara_res)
#else
!$omp parallel do private(i,ppara,ubr,ubz,ubphi,echarge2,ppara_res)
#endif
        do i = 1, lrzphi
          ubr = br(i,1,1)/babs(i,1,1)
          ubz = bz(i,1,1)/babs(i,1,1)
          ubphi = bphi(i,1,1)/babs(i,1,1)

          ppara =(prs_r(i,1,1)*ubr + prs_z(i,1,1)*ubz + prs_phi(i,1,1)*ubphi)*vac(i,1,1)
          ppara_res =(prs_r(i,1,1)*ubr + prs_z(i,1,1)*ubz + prs_phi(i,1,1)*ubphi)*vac_pe(i,1,1)
          
          echarge2 = 1.0d0/(e_e*dns_e(i,1,1) )

          utotphi(i,1,1) = uphi(i,1,1) + utor(i,1,1)
          vtotphi(i,1,1) = utotphi(i,1,1)*rho(i,1,1)

          epara0(i,1,1) = ppara*echarge2

          er0(i,1,1) =-(uphi(i,1,1)*bz(i,1,1) - uz(i,1,1)*bphi(i,1,1)) &
                     + ppara*ubr*echarge2

          ez0(i,1,1) =-(ur(i,1,1)*bphi(i,1,1) - uphi(i,1,1)*br(i,1,1)) &
                     + ppara*ubz*echarge2

          ephi0(i,1,1) =-(uz(i,1,1)*br(i,1,1) - ur(i,1,1)*bz(i,1,1)) &
                       + ppara*ubphi*echarge2

          er0_res(i,1,1) = eta_spt(i,1,1)*cr(i,1,1) &
                         -(uphi(i,1,1)*bz(i,1,1) - uz(i,1,1)*bphi(i,1,1)) &
                         + ppara_res*ubr*echarge2

          ez0_res(i,1,1) = eta_spt(i,1,1)*cz(i,1,1) &
                         -(ur(i,1,1)*bphi(i,1,1) - uphi(i,1,1)*br(i,1,1)) &
                         + ppara_res*ubz*echarge2

          ephi0_res(i,1,1) = eta_spt(i,1,1)*cphi(i,1,1) &
                           -(uz(i,1,1)*br(i,1,1) - ur(i,1,1)*bz(i,1,1)) &
                           + ppara_res*ubphi*echarge2

          er(i,1,1) = 0.0d0
          ez(i,1,1) = 0.0d0
          ephi(i,1,1)= 0.0d0
          epara(i,1,1)= 0.0d0

          er_res(i,1,1) = 0.0d0
          ez_res(i,1,1) = 0.0d0
          ephi_res(i,1,1)= 0.0d0

          cr0(i,1,1) = cr(i,1,1)
          cz0(i,1,1) = cz(i,1,1)
          cphi0(i,1,1)= cphi(i,1,1)

        end do

      else

! rho, prs, and parallel momentum are given by ion partilces
#ifdef AMDGPU
!$omp target teams distribute parallel do private(i,ppara,ubr,ubz,ubphi,echarge2,ppara_res)
#else
!$omp parallel do private(i,ppara,ubr,ubz,ubphi,echarge2,ppara_res)
#endif
      do i = 1, lrzphi
          ubr = br(i,1,1)/babs(i,1,1)
          ubz = bz(i,1,1)/babs(i,1,1)
          ubphi = bphi(i,1,1)/babs(i,1,1)

          ppara =(prs_r(i,1,1)*ubr + prs_z(i,1,1)*ubz + prs_phi(i,1,1)*ubphi)*vac(i,1,1)
          ppara_res =(prs_r(i,1,1)*ubr + prs_z(i,1,1)*ubz + prs_phi(i,1,1)*ubphi)*vac_pe(i,1,1)

          echarge2 = 1.0d0/(e_e*dns_e(i,1,1) )

          utotphi(i,1,1) = uphi(i,1,1) + utor(i,1,1)
          vtotphi(i,1,1) = utotphi(i,1,1)*rho(i,1,1)

          epara(i,1,1) = ppara*echarge2 - epara0(i,1,1)

          er(i,1,1) =-(uphi(i,1,1)*bz(i,1,1) - uz(i,1,1)*bphi(i,1,1)) &
                    + ppara*ubr*echarge2 - er0(i,1,1)

          ez(i,1,1) =-(ur(i,1,1)*bphi(i,1,1) - uphi(i,1,1)*br(i,1,1)) &
                    + ppara*ubz*echarge2 - ez0(i,1,1)

          ephi(i,1,1) =-(uz(i,1,1)*br(i,1,1) - ur(i,1,1)*bz(i,1,1)) &
                      + ppara*ubphi*echarge2 - ephi0(i,1,1)

          er_res(i,1,1) = eta_spt(i,1,1)*cr(i,1,1) - er0_res(i,1,1) &
                        -(uphi(i,1,1)*bz(i,1,1) - uz(i,1,1)*bphi(i,1,1)) &
                        + ppara_res*ubr*echarge2
          
          ez_res(i,1,1) = eta_spt(i,1,1)*cz(i,1,1) - ez0_res(i,1,1) &
                        -(ur(i,1,1)*bphi(i,1,1) - uphi(i,1,1)*br(i,1,1)) &
                        + ppara_res*ubz*echarge2
          
          ephi_res(i,1,1) = eta_spt(i,1,1)*cphi(i,1,1) - ephi0_res(i,1,1) &
                          -(uz(i,1,1)*br(i,1,1) - ur(i,1,1)*bz(i,1,1)) &
                          + ppara_res*ubphi*echarge2
      end do

      end if

! velocity ur, uz, uphi and vorticity
      call rotation(0,ur,uz,utotphi,omegar,omegaz,omegaphi)
      call divergence(0,ur,uz,utotphi,div_u)
      call periodic_field_mlt7(er_res,ez_res,ephi_res,omegar,omegaz,omegaphi,div_u) !2025-10-05

end
!--------------------------------------------------------------------
subroutine write_check_kti
!--------------------------------------------------------------------
      use mpiset
      use grid
      use equi_sol, only:raxis
      use field
      use check
      implicit none
      real(8)::e_kin,e_mag,e_thr,e_ep
      real(8)::e_kin_total,e_mag_total,e_thr_total,e_ep_total,e_total
      real(8)::e_trans_total(lkinetic),e_dissp_total !2014-10-24
      real(8)::e_0_free
      real(8)::wat
      real(8)::delta_br(lr,lz,lphi),delta_bz(lr,lz,lphi),delta_bphi(lr,lz,lphi)
      real(8)::delta_prs(lr,lz,lphi)
      real(8)::delta_vr(lr,lz,lphi),delta_vz(lr,lz,lphi),delta_vphi(lr,lz,lphi)
      real(8)::delta_rho(lr,lz,lphi)
      real(8)::delta_er(lr,lz,lphi),delta_ez(lr,lz,lphi),delta_ephi(lr,lz,lphi) !2015-10-09

      real(8)::br_n(lphi_n_size,lr,lz),bz_n(lphi_n_size,lr,lz)
      real(8)::bphi_n(lphi_n_size,lr,lz),prs_n(lphi_n_size,lr,lz)
      real(8)::vr_n(lphi_n_size,lr,lz),vz_n(lphi_n_size,lr,lz)
      real(8)::vphi_n(lphi_n_size,lr,lz),rho_n(lphi_n_size,lr,lz)
!2015-10-09s
      real(8)::er_n(lphi_n_size,lr,lz),ez_n(lphi_n_size,lr,lz)
      real(8)::ephi_n(lphi_n_size,lr,lz)
!2015-10-09e

      real(8)::energy_n(0:lphi_n),energy_n_total(0:lphi_n)
      real(8)::energy_0(3),energy_0_total(3),factor_nrm
      real(8)::dssp(4,0:lphi_n),dssp_total(4,0:lphi_n)
      real(8)::delta_cr(lr,lz,lphi),delta_cz(lr,lz,lphi),delta_cphi(lr,lz,lphi)
      real(8)::domegar(lr,lz,lphi),domegaz(lr,lz,lphi),domegaphi(lr,lz,lphi)
      real(8)::ddiv_v(lr,lz,lphi)
      integer::i,j,k,n,n1,n2
      integer::ititle=0


      if(my_rank.eq.0)then
        if(ititle.eq.0)then
          write(7,'(//)')
        end if
        write(7,600)'kstep=',kstep,'  t=',t
      end if
600   format(a,i8,a,1pe14.5,a,2i8)

!--------------------------------------------------------------------
! check of energy evolution

!$omp parallel do
      do i = 1, lrzphi
        delta_br(i,1,1) = br(i,1,1) - br0(i,1,1)
        delta_bz(i,1,1) = bz(i,1,1) - bz0(i,1,1)
        delta_bphi(i,1,1)=bphi(i,1,1)-bphi0(i,1,1)
        delta_prs(i,1,1)= prs(i,1,1)- prs0(i,1,1)
      end do

      e_kin = 0.0d0
      e_mag = 0.0d0
      e_thr = 0.0d0
      e_ep  = 0.0d0

!$omp parallel do reduction(+:e_kin,e_mag,e_thr,e_ep) 
      do k = lphistart, lphiend
      do j = lzstart, lzend
      do i = lrstart, lrend
        e_kin  = e_kin &
               + 0.50d0*rho(i,j,k)*grr(i,j,k) &
                       *(ur(i,j,k)**2 + uz(i,j,k)**2 + uphi(i,j,k)**2 )
        e_mag  = e_mag &
               + 0.50d0*grr(i,j,k) &
                *( delta_br(i,j,k)*(br(i,j,k)+br0(i,j,k)) &
                  +delta_bz(i,j,k)*(bz(i,j,k)+bz0(i,j,k)) &
                  +delta_bphi(i,j,k)*(bphi(i,j,k)+bphi0(i,j,k)) &
                 )

        e_thr = e_thr &
             + 0.50d0*(      (ppara_i(i,j,k)-ppara_i0(i,j,k) ) &
                      +2.0d0*(pperp_i(i,j,k)-pperp_i0(i,j,k) ) &
                      )*grr(i,j,k) &
             + delta_prs(i,j,k)/(gamma-1.0d0)*grr(i,j,k) !2025-04-27

        e_ep = e_ep &
             + 0.50d0*(      (ppara_a(i,j,k)-ppara_a0(i,j,k) ) &
                      +2.0d0*(pperp_a(i,j,k)-pperp_a0(i,j,k) ) &
                      )*grr(i,j,k)
      end do
      end do
      end do

      e_kin = e_kin*dr*dz*dphi
      e_mag = e_mag*dr*dz*dphi
      e_thr = e_thr*dr*dz*dphi
      e_ep  = e_ep *dr*dz*dphi

      call mpi_reduce(e_kin,e_kin_total,1,mpi_real8 &
                     ,mpi_sum,0,mhd_world,mpi_err)
      call mpi_reduce(e_mag,e_mag_total,1,mpi_real8 &
                     ,mpi_sum,0,mhd_world,mpi_err)
      call mpi_reduce(e_thr,e_thr_total,1,mpi_real8 &
                     ,mpi_sum,0,mhd_world,mpi_err)
      call mpi_reduce(e_ep,e_ep_total,1,mpi_real8 &
                     ,mpi_sum,0,mhd_world,mpi_err)

      if(ititle.eq.0)then
        ititle = 1
        if(my_rank.eq.0)then
#ifdef EXPEQ
          write(40,'(13(x,a))')'kstep','t[ms]','kin_energy','mag_energy'&
#else
          write(40,'(13(x,a))')'kstep','wat','kin_energy','mag_energy' &
#endif
                              ,'ep_energy','ion_energy' &
                              ,'ep_energy_trns','ion_energy_trns','elcl_energy_trns' &
                              ,'mag_int','mag_ideal' &
                              ,'dssp_energy','total_energy'

#ifdef EXPEQ
          write(41,'(32(x,a))')'kstep','t[ms]','e_n=0','e_n=00','e_n=0th' &
#else
          write(41,'(32(x,a))') 'kstep','wat','e_n=0','e_n=00','e_n=0th' &
#endif
                             ,'e_dissp','e_n=0-e_dissp' &
                             ,'e_n=0_free','dssp_n=0','dssp_o_n=0','dssp_c_n=0','dssp_j_n=0' &
                             ,'e_n=1','dssp_n=1','dssp_o_n=1','dssp_c_n=1','dssp_j_n=1' &
                             ,'e_n=2','dssp_n=2','dssp_o_n=2','dssp_c_n=2','dssp_j_n=2' !&
!                             ,'e_n=3','dssp_n=3','dssp_o_n=3','dssp_c_n=3','dssp_j_n=3' &
!                             ,'e_n=4','dssp_n=4','dssp_o_n=4','dssp_c_n=4','dssp_j_n=4'
        end if
      end if

      call mpi_reduce(e_trans,e_trans_total,lkinetic,mpi_real8 ,mpi_sum,0&
           &,mhd_world,mpi_err)
      call mpi_reduce(e_dissp,e_dissp_total,1,mpi_real8 ,mpi_sum,0&
           &,mhd_world,mpi_err)

!      e_total = e_kin_total + e_mag_total + e_thr_total
      e_total = e_kin_total + e_mag_total + e_dissp_total !2015-10-09
!      do i = 1, lkinetic
      do i = 1, 3 !2017-05-06
        e_total = e_total + e_trans_total(i)
      end do

      if(my_rank.eq.0)then
#ifdef EXPEQ
        wat = t/omega_a*1.0d3
#else
        wat = t/raxis
#endif
        write(40,'(i10,12(1pe14.5))')kstep,wat,e_kin_total,e_mag_total,e_ep_total &
                 ,e_thr_total,(e_trans_total(i),i=1,lkinetic),e_dissp_total,e_total

        flush(40)
      end if

!--------------------------------------------------------------------
! check energy of each toroidal mode number

      call toroidal_mode_d(ur,vr_n)
      call toroidal_mode_d(uz,vz_n)
      call toroidal_mode_d(uphi,vphi_n)
      call toroidal_mode_d(delta_br,br_n)
      call toroidal_mode_d(delta_bz,bz_n)
      call toroidal_mode_d(delta_bphi,bphi_n)
      call toroidal_mode_d(delta_prs,prs_n)
      call toroidal_mode_d(rho,rho_n)
      call toroidal_mode_d(er,er_n)
      call toroidal_mode_d(ez,ez_n)
      call toroidal_mode_d(ephi,ephi_n)


!$omp parallel do
      do n = 0, lphi_n
        energy_n(n) = 0.0d0
      end do

      do n = 1, 3
        energy_0(n) = 0.0d0
      end do

      n = 0
      n1 = 2*n + 1
!$omp parallel do
      do j = lzstart, lzend
      do i = lrstart, lrend
        energy_0(1) = energy_0(1) &
       + 0.50d0*grr(i,j,1) &
               *(vr_n(n1,i,j)**2 + vz_n(n1,i,j)**2 &
                +vphi_n(n1,i,j)**2 )*rho_n(n1,i,j) &
       + 0.50d0*grr(i,j,1) &
               *( br_n(n1,i,j)**2 &
                 +bz_n(n1,i,j)**2 &
                 +bphi_n(n1,i,j)**2 &
                )

        energy_0(2) = energy_0(2) &
       + 0.50d0*grr(i,j,1) &
               *( br_n(n1,i,j)*2.0d0*br0(i,j,1) &
                 +bz_n(n1,i,j)*2.0d0*bz0(i,j,1) &
                 +bphi_n(n1,i,j)*2.0d0*bphi0(i,j,1) &
                )

        energy_0(3) = energy_0(3) &
       + prs_n(n1,i,j)/(gamma-1.0d0)*grr(i,j,1)
      end do
      end do

      energy_n(0) = energy_0(1) + energy_0(2) + energy_0(3)

!$omp parallel do private(n1,n2)
      do n = 1, lphi_n
        n1 = 2*n + 1
        n2 = 2*n + 2
      do j = lzstart, lzend
      do i = lrstart, lrend
        energy_n(n) = energy_n(n) &
       + 0.50d0*grr(i,j,1) &
               *( vr_n(n1,i,j)**2 &
                 +vz_n(n1,i,j)**2 &
                 +vphi_n(n1,i,j)**2 )*rho_n(1,i,j) &
       + 0.50d0*grr(i,j,1) &
               *( br_n(n1,i,j)**2 &
                 +bz_n(n1,i,j)**2 &
                 +bphi_n(n1,i,j)**2 ) &
       + 0.50d0*grr(i,j,1) &
               *( vr_n(n2,i,j)**2 &
                 +vz_n(n2,i,j)**2 &
                 +vphi_n(n2,i,j)**2 )*rho_n(1,i,j) &
       + 0.50d0*grr(i,j,1) &
               *( br_n(n2,i,j)**2 &
                 +bz_n(n2,i,j)**2 &
                 +bphi_n(n2,i,j)**2 ) 
      end do
      end do
      end do

      factor_nrm = dr*dz*phileng/dble(mpi_proc_phi)
      energy_n(0) = energy_n(0)*factor_nrm
      energy_0(1) = energy_0(1)*factor_nrm
      energy_0(2) = energy_0(2)*factor_nrm
      energy_0(3) = energy_0(3)*factor_nrm

!$omp parallel do
      do n = 1, lphi_n
        energy_n(n) = energy_n(n)*factor_nrm*0.50d0
      end do

      n = lphi_n + 1
      call mpi_reduce(energy_n(0),energy_n_total(0),n,mpi_real8 &
                     ,mpi_sum,0,mhd_world,mpi_err)
      call mpi_reduce(energy_0(1),energy_0_total(1),3,mpi_real8 &
                     ,mpi_sum,0,mhd_world,mpi_err)


!-------------- start dissipation analysis --------------------------


! put perturbation of toroidal mode number 4*n
      do n = 0, lphi_n

        n1 = 2*n + 1
        n2 = 2*n + 2

!$omp parallel do
        do k = 1, lphi
        do j = 1, lz
        do i = 1, lr
          delta_vr(i,j,k)  = vr_n(n1,i,j)*cos_phi(n,k) + vr_n(n2,i,j)*sin_phi(n,k)
          delta_vz(i,j,k)  = vz_n(n1,i,j)*cos_phi(n,k) + vz_n(n2,i,j)*sin_phi(n,k)
          delta_vphi(i,j,k)= vphi_n(n1,i,j)*cos_phi(n,k)+vphi_n(n2,i,j)*sin_phi(n,k)
          delta_br(i,j,k)  = br0(i,j,k) + br_n(n1,i,j)*cos_phi(n,k) + br_n(n2,i,j)*sin_phi(n,k)
          delta_bz(i,j,k)  = bz0(i,j,k) + bz_n(n1,i,j)*cos_phi(n,k) + bz_n(n2,i,j)*sin_phi(n,k)
          delta_bphi(i,j,k)= bphi0(i,j,k)+bphi_n(n1,i,j)*cos_phi(n,k)+bphi_n(n2,i,j)*sin_phi(n,k)
          delta_rho(i,j,k) = rho_n(1,i,j)
!          delta_prs(i,j,k) = prs0(i,j,k)+ prs_n(n1,i,j)*cos_phi(n,k)+ prs_n(n2,i,j)*sin_phi(n,k)
          delta_er(i,j,k)  = er_n(n1,i,j)*cos_phi(n,k) + er_n(n2,i,j)*sin_phi(n,k)
          delta_ez(i,j,k)  = ez_n(n1,i,j)*cos_phi(n,k) + ez_n(n2,i,j)*sin_phi(n,k)
          delta_ephi(i,j,k)= ephi_n(n1,i,j)*cos_phi(n,k)+ephi_n(n2,i,j)*sin_phi(n,k)
        end do
        end do
        end do


! vorticity and velocity divergence
        call rotation(1,delta_br,delta_bz,delta_bphi,delta_cr,delta_cz,delta_cphi)
        call rotation(1,delta_vr,delta_vz,delta_vphi,domegar,domegaz,domegaphi)
        call divergence(1,delta_vr,delta_vz,delta_vphi,ddiv_v)


      dssp(:,n) = 0.0d0

!$omp parallel do
      do k = lphistart, lphiend
      do j = lzstart, lzend
      do i = lrstart, lrend
! omega^2
        dssp(2,n) = dssp(2,n) + nu_spt(i,j,k)*delta_rho(i,j,k)*grr(i,j,k) &
                   *(omegar(i,j,k)**2 + omegaz(i,j,k)**2 + omegaphi(i,j,k)**2 )
! divv^2
        dssp(3,n) = dssp(3,n) + nu_spt(i,j,k)*delta_rho(i,j,k)*grr(i,j,k) &
                   *(ddiv_v(i,j,k)**2*4.0d0/3.0d0 )
! joule, 2015-10-09s
        dssp(4,n) = dssp(4,n) + grr(i,j,k) &
                   *(delta_cr(i,j,k)*delta_er(i,j,k) &
                    +delta_cz(i,j,k)*delta_ez(i,j,k) &
                    +delta_cphi(i,j,k)*delta_ephi(i,j,k) &
                    )
! 2015-10-09e
!        dssp(4,n) = dssp(4,n) + eta_spt(i,j,k)*grr(i,j,k) &
!                   *(delta_cr(i,j,k)*(delta_cr(i,j,k) - cr0(i,j,k) ) &
!                    +delta_cz(i,j,k)*(delta_cz(i,j,k) - cz0(i,j,k) ) &
!                    +delta_cphi(i,j,k)*(delta_cphi(i,j,k) - cphi0(i,j,k) ) &
!                    )
      end do
      end do
      end do

      dssp(1,n) = dssp(2,n) + dssp(3,n) + dssp(4,n)
! mutiply by raxis for time normalizaion by Alfven frequency
      dssp(:,n) = dssp(:,n)*dr*dz*dphi*raxis

      end do

      n = 4*(lphi_n + 1)
      call mpi_reduce(dssp(1,0),dssp_total(1,0),n,mpi_real8 &
                     ,mpi_sum,0,mhd_world,mpi_err)


! n=0 mode free energy
      e_0_free = energy_n_total(0) - e_dissp
      
      if(my_rank.eq.0)then
#ifdef EXPEQ
        wat = t/omega_a*1.0d3
#else
        wat = t/raxis
#endif
        write(41,'(i8,31(1pe14.5))')kstep,wat,energy_n_total(0),energy_0_total(2),energy_0_total(3) &
                             ,e_dissp,e_0_free &
                             ,energy_0_total(1),(dssp_total(i,0),i=1,4) &
                             ,energy_n_total(1),(dssp_total(i,1),i=1,4) &
                             ,energy_n_total(2),(dssp_total(i,2),i=1,4) !&
!                             ,energy_n_total(3),(dssp_total(i,3),i=1,4) &
!                             ,energy_n_total(4),(dssp_total(i,4),i=1,4)

        flush(41)
      end if

end
!--------------------------------------------------------------------
#endif
!subroutines mhd_kti, e_field_kti, and write_check_kti
!--------------------------------------------------------------------


!--------------------------------------------------------------------
subroutine initial_mhd_balance(flag_HM) !2012-08-16
!--------------------------------------------------------------------
      use field
      implicit none
      logical::flag_HM !2012-08-16

! for electromagnetic field
! 2012-03-26
!$omp parallel
!$omp workshare
      dbr0 = 0.0d0
      dbz0 = 0.0d0
      dbphi0 = 0.0d0
      dvr0 = 0.0d0
      dvz0 = 0.0d0
      dvphi0 = 0.0d0
      drho0 = 0.0d0
      dprs0 = 0.0d0
      denrg0 = 0.0d0
!$omp end workshare
!$omp end parallel
! 2012-03-26 end

#ifdef KTI
      call mhd_kti(0,flag_HM) !subroutines hm and mhd are unified
#else      
      call mhd(0,flag_HM) !subroutines hm and mhd are unified
#endif
      
end
!--------------------------------------------------------------------
subroutine mhd_boundary(flag_HM)
!--------------------------------------------------------------------
      use parameters
      use mpiset
      use field, only:br,bz,bphi,br0,bz0,bphi0,vr,vz,vphi,rho,prs &
                     ,babs,babs0,vac,rho0,prs0,enrg,enrg0,gamma !2021-05-05
      implicit none
      logical::flag_HM !2013-01-12
      integer::i,j,k

#ifdef AMDGPU
!$omp target teams distribute parallel do private(j)
#else
!$omp parallel do private(j)
#endif
      do j = 1, lzphi
        br(1,j,1) = br0(1,j,1)
        vr(1,j,1) = 0.0d0
        vz(1,j,1) = 0.0d0
        vphi(1,j,1) = 0.0d0

        br(lr,j,1) = br0(lr,j,1)
        vr(lr,j,1) = 0.0d0
        vz(lr,j,1) = 0.0d0
        vphi(lr,j,1) = 0.0d0
      end do

#ifdef AMDGPU
!$omp target teams distribute parallel do collapse(2) private(k,i)
#else
!$omp parallel do private(k,i)
#endif
      do k = 1, lphi
      do i = 1, lr
        bz(i,1,k) = bz0(i,1,k)
        vr(i,1,k) = 0.0d0
        vz(i,1,k) = 0.0d0
        vphi(i,1,k) = 0.0d0

        bz(i,lz,k) = bz0(i,lz,k)
        vr(i,lz,k) = 0.0d0
        vz(i,lz,k) = 0.0d0
        vphi(i,lz,k) = 0.0d0
      end do
      end do

!      if(flag_HM)then
         call periodic_field_mlt8(rho,vr,vz,vphi,prs,br,bz,bphi)
!      else
!        call periodic_field_mlt8(rho,vr,vz,vphi,enrg,br,bz,bphi)
!      end if


! mask unused region

#ifdef AMDGPU
!$omp target teams distribute parallel do collapse(2) private(i,j)
#else
!$omp parallel do private(i,j)
#endif         
      do i = 1, lrstart-lr_shd-1
      do j = 1, lzphi
        rho(i,j,1) = rho0(i,j,1)
        prs(i,j,1) = prs0(i,j,1)
!        enrg(i,j,1) = enrg0(i,j,1)
        vr(i,j,1) = 0.0d0
        vz(i,j,1) = 0.0d0
        vphi(i,j,1) = 0.0d0
        br(i,j,1) = br0(i,j,1)
        bz(i,j,1) = bz0(i,j,1)
        bphi(i,j,1) = bphi0(i,j,1)
      end do
      end do

#ifdef AMDGPU
!$omp target teams distribute parallel do collapse(2) private(i,j)
#else
!$omp parallel do private(i,j)
#endif         
      do i = lrend+lr_shd+1, lr
      do j = 1, lzphi
        rho(i,j,1) = rho0(i,j,1)
        prs(i,j,1) = prs0(i,j,1)
!        enrg(i,j,1) = enrg0(i,j,1)
        vr(i,j,1) = 0.0d0
        vz(i,j,1) = 0.0d0
        vphi(i,j,1) = 0.0d0
        br(i,j,1) = br0(i,j,1)
        bz(i,j,1) = bz0(i,j,1)
        bphi(i,j,1) = bphi0(i,j,1)
      end do
      end do

#ifdef AMDGPU
!$omp target teams distribute parallel do collapse(3) private(k,j,i)
#else
!$omp parallel do private(k,j,i)
#endif         
      do k = 1, lphi
      do j = 1, lzstart-lz_shd-1
      do i = 1, lr
        rho(i,j,k) = rho0(i,j,k)
        prs(i,j,k) = prs0(i,j,k)
        enrg(i,j,k) = enrg0(i,j,k) !2021-05-05
        vr(i,j,k) = 0.0d0
        vz(i,j,k) = 0.0d0
        vphi(i,j,k) = 0.0d0
        br(i,j,k) = br0(i,j,k)
        bz(i,j,k) = bz0(i,j,k)
        bphi(i,j,k) = bphi0(i,j,k)
      end do
      end do
      end do

#ifdef AMDGPU
!$omp target teams distribute parallel do collapse(3) private(k,j,i)
#else
!$omp parallel do private(k,j,i)
#endif         
      do k = 1, lphi
      do j = lzend+lz_shd+1, lz
      do i = 1, lr
        rho(i,j,k) = rho0(i,j,k)
        prs(i,j,k) = prs0(i,j,k)
!        enrg(i,j,k) = enrg0(i,j,k) !2021-05-05
        vr(i,j,k) = 0.0d0
        vz(i,j,k) = 0.0d0
        vphi(i,j,k) = 0.0d0
        br(i,j,k) = br0(i,j,k)
        bz(i,j,k) = bz0(i,j,k)
        bphi(i,j,k) = bphi0(i,j,k)
      end do
      end do
      end do

#ifdef AMDGPU
!$omp target teams distribute parallel do private(i)
#else
!$omp parallel do private(i)
#endif      
      do i = 1, lrzphi
#ifndef KTI
        enrg(i,1,1) = 0.50d0*(vr(i,1,1)**2 + vz(i,1,1)**2 + vphi(i,1,1)**2)/rho(i,1,1) &
             + 0.50d0*(br(i,1,1)**2 + bz(i,1,1)**2 + bphi(i,1,1)**2) &
             + prs(i,1,1)/(gamma - 1.0d0)
!        prs(i,1,1) =(gamma - 1.0d0)*(enrg(i,1,1) &
!            - 0.50d0*(vr(i,1,1)**2 + vz(i,1,1)**2 + vphi(i,1,1)**2)/rho(i,1,1) &
!            - 0.50d0*(br(i,1,1)**2 + bz(i,1,1)**2 + bphi(i,1,1)**2) &
!                            )
#endif
      end do

      
end

!--------------------------------------------------------------------
subroutine mhd_smoothing(flag_HM)
!--------------------------------------------------------------------
      use parameters
      use mpiset
      use field, only:br,bz,bphi,dbr,dbz,dbphi,br_ref,bz_ref,bphi_ref &
                     ,vr,vz,vphi,dvr,dvz,dvphi,vr_ref,vz_ref,vphi_ref &
                     ,rho,prs,drho,dprs,rho_ref,prs_ref,enrg,denrg,enrg_ref & !2021-05-05
                     ,gamma
      use grid, only:grr,dr,dz,dphi
      use check, only:e_dissp !2019-06-30

      implicit none
      logical::flag_HM !2013-01-12
      real(8)::cwwa,cwwc
      real(8)::de_dissp_local !2019-06-30
      integer::i,j,k !2019-06-30


!------------------- 2019-06-30s -------------------------
! dissipated energy
      de_dissp_local = 0.0d0

#ifdef AMDGPU
!$omp target teams distribute parallel do collapse(3) private(k,j,i) reduction(+:de_dissp_local) 
#else
!$omp parallel do reduction(+:de_dissp_local)
#endif
      do k = lphistart, lphiend
      do j = lzstart, lzend
      do i = lrstart, lrend
        de_dissp_local = de_dissp_local &
       + grr(i,j,k)*( (vr(i,j,k)**2 + vz(i,j,k)**2 + vphi(i,j,k)**2 &
                      ) / rho(i,j,k) &
                    + (br(i,j,k)**2 + bz(i,j,k)**2 + bphi(i,j,k)**2 &
                      ) &
                    )
      end do
      end do
      end do

#ifdef AMDGPU
!$omp target teams distribute parallel do collapse(3) private(k,j,i)
#else
!$omp parallel do
#endif
      do k = 1, lphi
      do j = 1, lz
      do i = 1, lr
        dbr(i,j,k) = br(i,j,k) - br_ref(i,j,k)
        dbz(i,j,k) = bz(i,j,k) - bz_ref(i,j,k)
        dbphi(i,j,k) = bphi(i,j,k) - bphi_ref(i,j,k)

        dvr(i,j,k) = vr(i,j,k) - vr_ref(i,j,k)
        dvz(i,j,k) = vz(i,j,k) - vz_ref(i,j,k)
        dvphi(i,j,k) = vphi(i,j,k) - vphi_ref(i,j,k)
#ifndef KTI
        drho(i,j,k) = rho(i,j,k) - rho_ref(i,j,k)
#endif
        dprs(i,j,k) = prs(i,j,k) - prs_ref(i,j,k)
      end do
      end do
      end do
!------------------- 2019-06-30e -------------------------

      cwwa = 0.5d0
      cwwc =-1.d0/6.d0

      call partsm1(dbr,cwwa)
      call partsm1(dbr,cwwc)
      call partsm1(dbz,cwwa)
      call partsm1(dbz,cwwc)
      call partsm1(dbphi,cwwa)
      call partsm1(dbphi,cwwc)

      call partsm1(dvr,cwwa)
      call partsm1(dvr,cwwc)
      call partsm1(dvz,cwwa)
      call partsm1(dvz,cwwc)
      call partsm1(dvphi,cwwa)
      call partsm1(dvphi,cwwc)
#ifndef KTI
      call partsm1(drho,cwwa)
      call partsm1(drho,cwwc)
#endif
      call partsm1(dprs,cwwa)
      call partsm1(dprs,cwwc)

#ifdef AMDGPU
!$omp target teams distribute parallel do collapse(3) private(k,j,i)
#else
!$omp parallel do
#endif
      do k = 1, lphi
      do j = 1, lz
      do i = 1, lr
        br(i,j,k) = dbr(i,j,k) + br_ref(i,j,k)
        bz(i,j,k) = dbz(i,j,k) + bz_ref(i,j,k)
        bphi(i,j,k) = dbphi(i,j,k) + bphi_ref(i,j,k)

        vr(i,j,k) = dvr(i,j,k) + vr_ref(i,j,k)
        vz(i,j,k) = dvz(i,j,k) + vz_ref(i,j,k)
        vphi(i,j,k) = dvphi(i,j,k) + vphi_ref(i,j,k)
#ifndef KTI
        rho(i,j,k) = drho(i,j,k) + rho_ref(i,j,k)
#endif
        prs(i,j,k) = dprs(i,j,k) + prs_ref(i,j,k)
      end do
      end do
      end do

! boundary condition

      call mhd_boundary(flag_HM)

!------------------- 2019-06-30s -------------------------

#ifdef AMDGPU
!$omp target teams distribute parallel do collapse(3) private(k,j,i) reduction(+:de_dissp_local) 
#else
!$omp parallel do reduction(+:de_dissp_local)
#endif
      do k = lphistart, lphiend
      do j = lzstart, lzend
      do i = lrstart, lrend
        de_dissp_local = de_dissp_local &
       - grr(i,j,k)*( (vr(i,j,k)**2 + vz(i,j,k)**2 + vphi(i,j,k)**2 & ! negative sign
                      ) / rho(i,j,k) &
                    + (br(i,j,k)**2 + bz(i,j,k)**2 + bphi(i,j,k)**2 &
                      ) &
                    )     
      end do
      end do
      end do

      e_dissp = e_dissp + de_dissp_local*dr*dz*dphi*0.50d0

#ifdef AMDGPU
!$omp target teams distribute parallel do collapse(3) private(k,j,i)
#else
!$omp parallel do
#endif
      do k = 1, lphi
      do j = 1, lz
      do i = 1, lr
        br_ref(i,j,k) = br(i,j,k)
        bz_ref(i,j,k) = bz(i,j,k)
        bphi_ref(i,j,k) = bphi(i,j,k)
        
        vr_ref(i,j,k) = vr(i,j,k)
        vz_ref(i,j,k) = vz(i,j,k)
        vphi_ref(i,j,k) = vphi(i,j,k)
#ifndef KTI
        rho_ref(i,j,k) = rho(i,j,k)
#endif
        prs_ref(i,j,k) = prs(i,j,k)
      end do
      end do
      end do

      call gradmg 

end
!--------------------------------------------------------------------
subroutine mhd_lowpass(flag_HM)
!--------------------------------------------------------------------
      use parameters
      use mpiset
      use field, only:br,bz,bphi,dbr,dbz,dbphi,br0,bz0,bphi0 &
                     ,vr,vz,vphi &
                     ,rho,prs,drho,dprs,rho0,prs0,enrg,enrg0,denrg !2021-05-05
      use grid, only:grr,dr,dz,dphi !2025-04-30
      use check, only:e_dissp !2025-04-30
      implicit none
      logical::flag_HM
      real(8)::de_dissp_local !2025-04-30
      integer::i,j,k !2025-04-30

!------------------- 2025-04-30s -------------------------
! dissipated energy
      de_dissp_local = 0.0d0

#ifdef AMDGPU
!$omp target teams distribute parallel do collapse(3) private(k,j,i) reduction(+:de_dissp_local) 
#else
!$omp parallel do reduction(+:de_dissp_local)
#endif
      do k = lphistart, lphiend
      do j = lzstart, lzend
      do i = lrstart, lrend
        de_dissp_local = de_dissp_local &
       + grr(i,j,k)*( (vr(i,j,k)**2 + vz(i,j,k)**2 + vphi(i,j,k)**2 &
                      ) / rho(i,j,k) &
                    + (br(i,j,k)**2 + bz(i,j,k)**2 + bphi(i,j,k)**2 &
                      ) &
                    )
      end do
      end do
      end do

#ifdef AMDGPU
!$omp target teams distribute parallel do collapse(3) private(k,j,i)
#else
!$omp parallel do collapse(3) private(k,j,i)
#endif
      do k = 1, lphi
      do j = 1, lz
      do i = 1, lr
        dbr(i,j,k)  = br(i,j,k) - br0(i,j,k)
        dbz(i,j,k)  = bz(i,j,k) - bz0(i,j,k)
        dbphi(i,j,k)= bphi(i,j,k) - bphi0(i,j,k)
#ifndef KTI
        drho(i,j,k) = rho(i,j,k) - rho0(i,j,k)
#endif
        dprs(i,j,k) = prs(i,j,k) - prs0(i,j,k)
      end do
      end do
      end do
     
      call lowpass(dbr)
      call lowpass(dbz)
      call lowpass(dbphi)
      call lowpass(vr)
      call lowpass(vz)
      call lowpass(vphi)
#ifndef KTI
      call lowpass(drho)
#endif
      call lowpass(dprs)

#ifdef AMDGPU
!$omp target teams distribute parallel do collapse(3) private(k,j,i)
#else
!$omp parallel do collapse(3) private(k,j,i)
#endif
      do k = 1, lphi
      do j = 1, lz
      do i = 1, lr
        br(i,j,k)  = dbr(i,j,k)   + br0(i,j,k)
        bz(i,j,k)  = dbz(i,j,k)   + bz0(i,j,k)
        bphi(i,j,k)= dbphi(i,j,k) + bphi0(i,j,k)
#ifndef KTI
        rho(i,j,k) = drho(i,j,k)  + rho0(i,j,k)
#endif
        prs(i,j,k) = dprs(i,j,k)  + prs0(i,j,k)
      end do
      end do
      end do

! boundary condition

      call mhd_boundary(flag_HM)

!------------------- 2025-04-30s -------------------------

#ifdef AMDGPU
!$omp target teams distribute parallel do collapse(3) private(k,j,i) reduction(+:de_dissp_local) 
#else
!$omp parallel do collapse(3) reduction(+:de_dissp_local)
#endif
      do k = lphistart, lphiend
      do j = lzstart, lzend
      do i = lrstart, lrend
        de_dissp_local = de_dissp_local &
       - grr(i,j,k)*( (vr(i,j,k)**2 + vz(i,j,k)**2 + vphi(i,j,k)**2 & ! negative sign
                      ) / rho(i,j,k) &
                    + (br(i,j,k)**2 + bz(i,j,k)**2 + bphi(i,j,k)**2 &
                      ) &
                    )     
      end do
      end do
      end do

      e_dissp = e_dissp + de_dissp_local*dr*dz*dphi*0.50d0

!------------------- 2025-04-30e -------------------------

end
!--------------------------------------------------------------------
subroutine perturb_fluid_random(flag_HM)
!--------------------------------------------------------------------
      use parameters
      use mpiset
      use field, only:dvr,dvz,dvphi,dbr,dbz,dbphi,br,bz,bphi &
                     ,br0,bz0,bphi0,br_ref,bz_ref,bphi_ref &
                     ,rho,vr,vz,vphi,prs,enrg,enrg_ref,gamma !2021-05-05
      use grid
      use random_num

      implicit none
      logical::flag_HM !2013-01-12
      real(8)::dlt0,a_limit
      real(8)::rand_fld(lrzphi)
      integer::i,j,k,n,nrandom

      nrandom = lrzphi
      call ransuu(nrandom, rand_fld)

      dlt0 = 1.0d-6

      a_limit = 0.50d0*rleng*0.70d0

1000  continue

      dvr = 0.0d0
      dvz = 0.0d0
      dvphi = 0.0d0

!$omp parallel do private(n)
      do k = 1, lphi
      do j = 1, lz
      do i = 1, lr
        n = i + (j-1)*lr + (k-1)*lrz
        if(gaa(i,j,k).lt.a_limit)then
          dvphi(i,j,k) = dlt0*rand_fld(n)
        end if
      end do
      end do
      end do

      call periodic_field(dvphi)
      call rotation(1,dvr,dvz,dvphi,dbr,dbz,dbphi)

      br = br0 + dbr
      bz = bz0 + dbz
      bphi = bphi0 + dbphi
#ifndef KTI
      enrg = 0.50d0*(vr**2 + vz**2 + vphi**2)/rho &
           + 0.50d0*(br**2 + bz**2 + bphi**2) &
           + prs/(gamma - 1.0d0)
#endif
      
      call mhd_boundary(flag_HM)

      br_ref = br
      bz_ref = bz
      bphi_ref = bphi
#ifndef KTI
      enrg_ref = enrg
#endif

      call gradmg

#ifdef KTI
          call e_field_kti(1)
#else
          call e_field(1)
#endif


end
!--------------------------------------------------------------------
subroutine perturb_fluid(flag_HM)
!--------------------------------------------------------------------
      use parameters
      use field
      use grid
      implicit none
      logical::flag_HM !2013-01-12
      real(8)::mask,x,dn,dm,dlt0,dk,delta_b,delta_b0
      integer::i,j,k,iteration

      dn = phimode
      dm = dble(int(phimode*1.50d0))

      delta_b0 = 1.0d-5
      dlt0 = 3.0d-8
      dk = 1.0d0*pi/(0.60d0*0.50d0*rleng)
      iteration = 0

1000  continue

      dvr = 0.0d0
      dvz = 0.0d0
      dvphi = 0.0d0

!$omp parallel do private(x)
      do k = 1, lphi
      do j = 1, lz
      do i = 1, lr
        x = dk*gaa(i,j,k)
        if(x.lt.pi)then
          dvphi(i,j,k) = dlt0*sin(x)**3 &
                             *cos(dn*gphi(i,j,k) + dm*gtht(i,j,k) )
        end if
      end do
      end do
      end do

      call rotation(1,dvr,dvz,dvphi,dbr,dbz,dbphi)

      br = br0 + dbr
      bz = bz0 + dbz
      bphi = bphi0 + dbphi
#ifndef KTI
      enrg  = 0.50d0*(vr**2 + vz**2 + vphi**2)/rho &
           + 0.50d0*(br**2 + bz**2 + bphi**2) &
           + prs/(gamma - 1.0d0)
#endif
      call mhd_boundary(flag_HM)

      br_ref = br
      bz_ref = bz
      bphi_ref = bphi
#ifndef KTI
      enrg_ref = enrg
#endif

      call gradmg

#ifdef KTI
          call e_field_kti(1)
#else
          call e_field(1)
#endif

end
!--------------------------------------------------------------------
subroutine cc_density
!--------------------------------------------------------------------
      use parameters
      use field, only:br,bz,bphi,cr,cz,cphi
      implicit none
      integer::i

! current density

      call rotation(0,br,bz,bphi,cr,cz,cphi) !2013-05-26
!      call rotation(1,br,bz,bphi,cr,cz,cphi) !2012-08-27

end
!-------------- hereafter infrastructure subroutines ----------------
!                     date:   2012-06-17
! HITACHI version 2013-04-12
! shodow widths are defined by lr_shd, lz_shd, lphi_shd, 2015-06-28
!-------------- hereafter infrastructure subroutines ----------------

!--------------------------------------------------------------------
subroutine start_mpi
! initialization of mpi
! modified on 2011-01-29, and for particle parallelization, 2012-06-17
! modified on 2012-06-29 for new data type: mpi_satellite
! modified on 2015-06-27 for definition of shadow width
!--------------------------------------------------------------------
      use mpiset
      implicit none
      integer::n,k,j,i,icolor,ikey,nmhd !2012-06-17
!      integer::count,blocklength(2),displacements(2),types(2) !2012-06-29
!!!HITACHI
      integer array_of_blocklength(2)
      integer array_of_displacement(2)
      integer array_of_types(2)
!!!HITACHI

      call mpi_init(mpi_err)
      call mpi_comm_rank(mpi_comm_world,my_rank,mpi_err)
      call mpi_comm_size(mpi_comm_world,nprocess,mpi_err)

!!!HITACHI
      array_of_blocklength(1)=7
      array_of_blocklength(2)=2
      array_of_displacement(1)=0
      array_of_displacement(2)=7*8
      array_of_types(1)=mpi_real8
      array_of_types(2)=mpi_integer
      call mpi_type_struct(2,array_of_blocklength,array_of_displacement,array_of_types,new_d7i2_type,mpi_err)
      call mpi_type_commit(new_d7i2_type,mpi_err)
!!!HITACHI

      my_rank_ptcl = my_rank/mpi_proc_mhd
      my_rank_mhd = mod(my_rank,mpi_proc_mhd)

!      if(my_rank.eq.0)then
!        write(7,*)'my rank,size=',my_rank,nprocess
!      end if

      do n = 0, nprocess-1
        nmhd = mod(n, mpi_proc_mhd)
        k = nmhd/mpi_proc_pol
        kphi_offset(n)= lphinet*k/mpi_proc_phi

        j =(nmhd - mpi_proc_pol*k)/mpi_proc_r
        if(mpi_proc_z.eq.1)then
          kz_offset(n)= lznet*j/mpi_proc_z
          kz_s(n) = 1
          kz_e(n) = lz
        else
          if(j.eq.0)then
            kz_offset(n)= lznet*j/mpi_proc_z
            kz_s(n) = 1
            kz_e(n) = lz - 2*lz_shd
          else if(j.eq.mpi_proc_z-1)then
            kz_offset(n)= lznet*j/mpi_proc_z - 2*lz_shd
            kz_s(n) = 2*lz_shd + 1
            kz_e(n) = lz
          else
            kz_offset(n)= lznet*j/mpi_proc_z - lz_shd
            kz_s(n) = lz_shd + 1
            kz_e(n) = lz - lz_shd
          end if
        end if

        i = nmhd - mpi_proc_pol*k - mpi_proc_r*j !2012-06-17
        if(mpi_proc_r.eq.1)then
          kr_offset(n)= lrnet*i/mpi_proc_r
          kr_s(n) = 1
          kr_e(n) = lr
        else
          if(i.eq.0)then
            kr_offset(n)= lrnet*i/mpi_proc_r
            kr_s(n) = 1
            kr_e(n) = lr - 2*lr_shd
          else if(i.eq.mpi_proc_r-1)then
            kr_offset(n)= lrnet*i/mpi_proc_r - 2*lr_shd !2015-06-27
            kr_s(n) = 2*lr_shd + 1
            kr_e(n) = lr
          else
            kr_offset(n)= lrnet*i/mpi_proc_r - lr_shd !2015-06-27
            kr_s(n) = lr_shd + 1
            kr_e(n) = lr - lr_shd
          end if
        end if
      end do

      my_rank_phi = my_rank_mhd/mpi_proc_pol
      my_rank_z   =(my_rank_mhd - mpi_proc_pol*my_rank_phi) &
                  &/mpi_proc_r
      my_rank_r   = my_rank_mhd - mpi_proc_pol*my_rank_phi  &
                  - mpi_proc_r*my_rank_z
      my_rank_pol = my_rank_r + my_rank_z*mpi_proc_r !2015-07-17

      if(mpi_proc_r.eq.1)then
        lrstart = 1
        lrend   = lr
      else
        if(my_rank_r.eq.0)then
          lrstart = 1
          lrend   = lr - 2*lr_shd
        else if(my_rank_r.eq.(mpi_proc_r-1))then
          lrstart = 2*lr_shd + 1
          lrend   = lr
        else
          lrstart = lr_shd + 1
          lrend   = lr - lr_shd
        end if
      end if

      if(mpi_proc_z.eq.1)then
        lzstart = 1
        lzend   = lz
      else
        if(my_rank_z.eq.0)then
          lzstart = 1
          lzend   = lz - 2*lz_shd
        else if(my_rank_z.eq.(mpi_proc_z-1))then
          lzstart = 2*lz_shd + 1
          lzend   = lz
        else
          lzstart = lz_shd + 1
          lzend   = lz - lz_shd
        end if
      end if

      lphistart = lphi_shd + 1
      lphiend   = lphi - lphi_shd

!      lphiend   = 2 + lphinet*(my_rank_phi + 1)/mpi_proc_phi &
!                    - lphinet* my_rank_phi     /mpi_proc_phi

!-------------------------------------------
! definition of data type in z-communication

      isize(1) = lr
      isize(2) = lz
      isize(3) = lphi 

      isubsize(1) = lr
      isubsize(2) = lz_shd
      isubsize(3) = lphi 

! isend and irecv denote starting points of the subdata
! the first point is zero

      isend_up(1) = 0
      isend_up(2) = lzend - lz_shd
      isend_up(3) = 0
      irecv_up(1) = 0
      irecv_up(2) = lzend
      irecv_up(3) = 0

      isend_down(1) = 0
      isend_down(2) = lzstart - 1
      isend_down(3) = 0
      irecv_down(1) = 0
      irecv_down(2) = lzstart - lz_shd - 1
      irecv_down(3) = 0

      if(my_rank_z.ne.0)then
        call mpi_type_create_subarray(3,isize,isubsize,isend_down &
             &,mpi_order_fortran,mpi_real8,jz_senddown,mpi_err)
        call mpi_type_create_subarray(3,isize,isubsize,irecv_down &
             &,mpi_order_fortran,mpi_real8,jz_recvdown,mpi_err)
        call mpi_type_commit(jz_senddown,mpi_err)
        call mpi_type_commit(jz_recvdown,mpi_err)
      end if

      if(my_rank_z.ne.(mpi_proc_z-1))then
        call mpi_type_create_subarray(3,isize,isubsize,isend_up &
             &,mpi_order_fortran,mpi_real8,jz_sendup,mpi_err)
        call mpi_type_create_subarray(3,isize,isubsize,irecv_up &
             &,mpi_order_fortran,mpi_real8,jz_recvup,mpi_err)
        call mpi_type_commit(jz_sendup,mpi_err)
        call mpi_type_commit(jz_recvup,mpi_err)
      end if

!-------------------------------------------
! definition of data type in r-communication

      isize(1) = lr
      isize(2) = lz
      isize(3) = lphi 

      isubsize(1) = lr_shd
      isubsize(2) = lz
      isubsize(3) = lphi 

      isend_up(1) = lrend - lr_shd
      isend_up(2) = 0
      isend_up(3) = 0
      irecv_up(1) = lrend
      irecv_up(2) = 0
      irecv_up(3) = 0

      isend_down(1) = lrstart - 1
      isend_down(2) = 0
      isend_down(3) = 0
      irecv_down(1) = lrstart - lr_shd - 1
      irecv_down(2) = 0
      irecv_down(3) = 0

      if(my_rank_r.ne.0)then
        call mpi_type_create_subarray(3,isize,isubsize,isend_down &
             &,mpi_order_fortran,mpi_real8,ir_senddown,mpi_err)
        call mpi_type_create_subarray(3,isize,isubsize,irecv_down &
             &,mpi_order_fortran,mpi_real8,ir_recvdown,mpi_err)
        call mpi_type_commit(ir_senddown,mpi_err)
        call mpi_type_commit(ir_recvdown,mpi_err)
      end if

      if(my_rank_r.ne.(mpi_proc_r-1))then
        call mpi_type_create_subarray(3,isize,isubsize,isend_up &
             &,mpi_order_fortran,mpi_real8,ir_sendup,mpi_err)
        call mpi_type_create_subarray(3,isize,isubsize,irecv_up &
             &,mpi_order_fortran,mpi_real8,ir_recvup,mpi_err)
        call mpi_type_commit(ir_sendup,mpi_err)
        call mpi_type_commit(ir_recvup,mpi_err)
      end if

!-------------------------------------------
! definition of data type in particle parallelization

      isize(1) = lr
      isize(2) = lz
      isize(3) = lphi 

#ifdef PRM1MPI
      isubsize(1) = lr
      isubsize(2) = lz
#else
      isubsize(1) = lr - 2*lr_shd
      isubsize(2) = lz - 2*lz_shd
#endif

      isubsize(3) = lphi - 2*lphi_shd

      isend_up(1) = lrstart - 1
      isend_up(2) = lzstart - 1
      isend_up(3) = lphistart - 1

      call mpi_type_create_subarray(3,isize,isubsize,isend_up &
           &,mpi_order_fortran,mpi_real8,ptcl_allreduce,mpi_err)

      call mpi_type_commit(ptcl_allreduce,mpi_err)

!-------------------------------------------
! communicator for communication in r direction
! 2011-01-29
! 2012-06-17
      icolor = my_rank_ptcl*(mpi_proc_z*mpi_proc_phi) &
             + my_rank_phi*mpi_proc_z + my_rank_z
      ikey = my_rank_r
      call mpi_comm_split(mpi_comm_world,icolor,ikey,r_world,mpi_err)

!-------------------------------------------
! communicator for communication in z direction
! 2011-01-29
! 2012-06-17
      icolor = my_rank_ptcl*(mpi_proc_r*mpi_proc_phi) &
             + my_rank_phi*mpi_proc_r + my_rank_r
      ikey = my_rank_z
      call mpi_comm_split(mpi_comm_world,icolor,ikey,z_world,mpi_err)

!-------------------------------------------
! communicator for fourier decomposition in phi direction
! 2012-06-17
      icolor = my_rank_ptcl*(mpi_proc_r*mpi_proc_z) &
             + my_rank_z*mpi_proc_r + my_rank_r
      ikey = my_rank_phi
      call mpi_comm_split(mpi_comm_world,icolor,ikey,phi_world&
           &,mpi_err)

!-------------------------------------------
! communicator for snapshot on poloidal plane
! 2012-06-17
      icolor = my_rank_ptcl*mpi_proc_phi + my_rank_phi
      ikey = my_rank_z*mpi_proc_r + my_rank_r
      call mpi_comm_split(mpi_comm_world,icolor,ikey,poloidal_world&
           &,mpi_err)

!-------------------------------------------
! communicator for mhd
! 2012-06-17
      icolor = my_rank_ptcl
      ikey = my_rank_mhd
      call mpi_comm_split(mpi_comm_world,icolor,ikey,mhd_world&
           &,mpi_err)

!-------------------------------------------------
! communicator for particle parallelized processes
! 2012-06-17
      icolor = my_rank_mhd
      ikey = my_rank_ptcl
      call mpi_comm_split(mpi_comm_world,icolor,ikey,ptcl_world&
           &,mpi_err)


!-------------------------------------------------
! new data type for satellite communication
!            launch_satellite3, return_satellite3
! 2012-06-29

!      count = 2
!      blocklength(1) = 3
!      blocklength(2) = 7
!      displacements(1) = 0
!      displacements(2) = 3*4
!      types(1) = mpi_integer
!      types(2) = mpi_real8

!      call mpi_type_struct(count,blocklength,displacements,types &
!                          ,mpi_satellite,mpi_err)

!      call mpi_type_commit(mpi_satellite,mpi_err)

!-------------------------------------------------
! new data type for particle communication
!            com_particle3
! 2012-06-29

!      count = 2
!      blocklength(1) = 1
!      blocklength(2) = 10
!      displacements(1) = 0
!      displacements(2) = 1*4
!      types(1) = mpi_integer
!      types(2) = mpi_real8

!      call mpi_type_struct(count,blocklength,displacements,types &
!                          ,mpi_ptcl,mpi_err)

!      call mpi_type_commit(mpi_ptcl,mpi_err)


end
!--------------------------------------------------------------------
subroutine end_mpi
! finalize mpi
!--------------------------------------------------------------------
      use mpiset
      implicit none

      call mpi_finalize(mpi_err)

end
!--------------------------------------------------------------------
subroutine toroidal_mode_d(aaa,aaa_n)
!--------------------------------------------------------------------
      use check
      use mpiset
      implicit none
      integer::i,k,n,n1,n2

      real(8)::aaa(lr,lz,lphi)
      real(8)::work(lphi_n_size,lr,lz),aaa_n(lphi_n_size,lr,lz)

      work = 0.0d0

#ifdef AMDGPU
!$omp target teams distribute parallel do collapse(2) private(i,n,n1,n2,k)
#else      
!$omp parallel do private(n1,n2)
#endif
      do i = 1, lrz
        do n = 0, lphi_n
          n1 = 2*n + 1
          n2 = 2*n + 2
          do k = lphistart, lphiend
            work(n1,i,1) = work(n1,i,1) + aaa(i,1,k)*cos_phi(n,k)
            work(n2,i,1) = work(n2,i,1) + aaa(i,1,k)*sin_phi(n,k)
          end do
        end do
      end do

#ifdef AMDGPU
!$omp target teams distribute parallel do private(i,n)
#else
!$omp parallel do
#endif
      do i = 1, lrz
        work(1,i,1) = work(1,i,1)/dble(lphinet)
        work(2,i,1) = 0.0d0
        do n = 3, lphi_n_size
          work(n,i,1) = work(n,i,1)/dble(lphinet)*2.0d0
        end do
      end do


      n = lrz*lphi_n_size
! reduction in phi direction using communicator 'phi_world'
      call mpi_allreduce(work(1,1,1),aaa_n(1,1,1),n,mpi_real8 &
                        ,mpi_sum,phi_world,mpi_err)

end
!--------------------------------------------------------------------
subroutine toroidal_mode_d_mlt(imulti,aaa,aaa_n)
!--------------------------------------------------------------------
      use check
      use mpiset
      implicit none
      integer::i,k,n,n1,n2,imulti,j

      real(8)::aaa(imulti,lr,lz,lphi)
      real(8)::work(imulti,lr,lz,lphi_n_size),aaa_n(imulti,lr,lz,lphi_n_size)

      work = 0.0d0

#ifndef AMDGPU
!$omp parallel do private(n1,n2) reduction(+:work) !2024-12-16, correction suggested by Jialei Wang
      do k = lphistart, lphiend
        do n = 0, lphi_n
          n1 = 2*n + 1
          n2 = 2*n + 2
          do i = 1, lrz
            do j = 1, imulti
              work(j,i,1,n1) = work(j,i,1,n1) + aaa(j,i,1,k)*cos_phi(n,k)
              work(j,i,1,n2) = work(j,i,1,n2) + aaa(j,i,1,k)*sin_phi(n,k)
            end do
          end do
        end do
      end do
#else
!$omp target teams distribute parallel do collapse(3) private(n,i,j,n1,n2,k)
      do n = 0, lphi_n
        do i = 1, lrz
          do j = 1, imulti
            n1 = 2*n + 1
            n2 = 2*n + 2
            do k = lphistart, lphiend
              work(j,i,1,n1) = work(j,i,1,n1) + aaa(j,i,1,k)*cos_phi(n,k)
              work(j,i,1,n2) = work(j,i,1,n2) + aaa(j,i,1,k)*sin_phi(n,k)
            end do
          end do
        end do
      end do
#endif

#ifdef AMDGPU
!$omp target teams distribute parallel do collapse(2) private(i,j)
#else      
!$omp parallel do
#endif
      do i = 1, lrz
        do j = 1, imulti
          work(j,i,1,1) = work(j,i,1,1)/dble(lphinet)
          work(j,i,1,2) = 0.0d0
        end do
     end do
     
#ifdef AMDGPU
!$omp target teams distribute parallel do collapse(3) private(n,i,j)
#else
!$omp parallel do
#endif
      do n = 3, lphi_n_size
        do i = 1, lrz
          do j = 1, imulti
            work(j,i,1,n) = work(j,i,1,n)/dble(lphinet)*2.0d0
          end do
        end do
      end do

      n = lrz*lphi_n_size*imulti
! reduction in phi direction using communicator 'phi_world'
      call mpi_allreduce(work(1,1,1,1),aaa_n(1,1,1,1),n,mpi_real8 &
                        ,mpi_sum,phi_world,mpi_err)

end
!--------------------------------------------------------------------
subroutine lowpass(aaa)
!--------------------------------------------------------------------
      use check
      use mpiset
      implicit none
      integer::i,k,n,n1,n2

      real(8)::aaa(lr,lz,lphi)
      real(8)::aaa_n(lphi_n_size,lr,lz)

! fourier decomposition in toroidal direction

      call toroidal_mode_d(aaa,aaa_n)

      aaa = 0.0d0

#ifndef AMDGPU
      do k = 1, lphi
        do n = 0, lphi_n
          n1 = 2*n + 1
          n2 = 2*n + 2
!$omp parallel do
        do i = 1, lrz
          aaa(i,1,k) = aaa(i,1,k) + aaa_n(n1,i,1)*cos_phi(n,k) &
                                  + aaa_n(n2,i,1)*sin_phi(n,k)
        end do
        end do
      end do
#else
!$omp target teams distribute parallel do collapse(2) private(k,i,n,n1,n2)
      do k = 1, lphi
        do i = 1, lrz
          do n = 0, lphi_n
            n1 = 2*n + 1
            n2 = 2*n + 2
            aaa(i,1,k) = aaa(i,1,k) + aaa_n(n1,i,1)*cos_phi(n,k) &
                                    + aaa_n(n2,i,1)*sin_phi(n,k)
          end do
        end do
      end do
#endif      

end subroutine lowpass
!--------------------------------------------------------------------
subroutine lowpass2(n_min,n_max,aaa)
!--------------------------------------------------------------------
      use check
      use mpiset
      implicit none
      integer::i,k,n,n1,n2,n_min,n_max

      real(8)::aaa(lr,lz,lphi)
      real(8)::aaa_n(lphi_n_size,lr,lz)

! fourier decomposition in toroidal direction
      n_max = min(n_max, lphi_n) !2017-02-12

      call toroidal_mode_d(aaa,aaa_n)

      aaa = 0.0d0

#ifndef AMDGPU
      do k = 1, lphi
!        do n = 0, lphi_n
        do n = n_min, n_max
          n1 = 2*n + 1
          n2 = 2*n + 2
!$omp parallel do
        do i = 1, lrz
          aaa(i,1,k) = aaa(i,1,k) + aaa_n(n1,i,1)*cos_phi(n,k) &
                                  + aaa_n(n2,i,1)*sin_phi(n,k)
        end do
        end do
      end do
#else
!$omp target teams distribute parallel do collapse(2) private(k,i,n,n1,n2)
      do k = 1, lphi
        do i = 1, lrz
          do n = n_min, n_max
            n1 = 2*n + 1
            n2 = 2*n + 2
            aaa(i,1,k) = aaa(i,1,k) + aaa_n(n1,i,1)*cos_phi(n,k) &
                                    + aaa_n(n2,i,1)*sin_phi(n,k)
          end do
        end do
      end do
#endif

end
!--------------------------------------------------------------------
subroutine lowpass2_mlt2(n_min,n_max,aaa1,aaa2)
!--------------------------------------------------------------------
      use check
      use mpiset
      implicit none
      integer, parameter::lmulti=2
      integer::i,j,k,n,n1,n2,n_min,n_max,imulti=lmulti

      real(8)::aaa1(lr,lz,lphi),aaa2(lr,lz,lphi)
      real(8)::aaa(lmulti,lr,lz,lphi)
      real(8)::aaa_n(lmulti,lr,lz,lphi_n_size)

! fourier decomposition in toroidal direction
      n_max = min(n_max, lphi_n) !2017-02-12

#ifdef AMDGPU
!$omp target teams distribute parallel do collapse(3) private(k,j,i)
      do k=1,lphi
      do j=1,lz
      do i=1,lr
        aaa(1,i,j,k) = aaa1(i,j,k)
        aaa(2,i,j,k) = aaa2(i,j,k)
      end do
      end do
      end do
#else
!$omp parallel
!$omp workshare
      aaa(1,:,:,:) = aaa1(:,:,:)
      aaa(2,:,:,:) = aaa2(:,:,:)
!$omp end workshare
!$omp end parallel
#endif      
      
      call toroidal_mode_d_mlt(imulti,aaa,aaa_n)

#ifdef AMDGPU
!$omp target teams distribute parallel do collapse(3) private(k,j,i)
      do k=1,lphi
      do j=1,lz
      do i=1,lr
        aaa1(i,j,k) = 0.0d0
        aaa2(i,j,k) = 0.0d0
      end do
      end do
      end do
#else
!$omp parallel
!$omp workshare
      aaa1 = 0.0d0
      aaa2 = 0.0d0
!$omp end workshare
!$omp end parallel
#endif      


#ifndef AMDGPU
      do k = 1, lphi
!        do n = 0, lphi_n
        do n = n_min, n_max
          n1 = 2*n + 1
          n2 = 2*n + 2

!$omp parallel do
        do i = 1, lrz
          aaa1(i,1,k) = aaa1(i,1,k) + aaa_n(1,i,1,n1)*cos_phi(n,k) &
                                    + aaa_n(1,i,1,n2)*sin_phi(n,k)
          aaa2(i,1,k) = aaa2(i,1,k) + aaa_n(2,i,1,n1)*cos_phi(n,k) &
                                    + aaa_n(2,i,1,n2)*sin_phi(n,k)
        end do
        end do
      end do
#else
!$omp target teams distribute parallel do collapse(2) private(k,i,n,n1,n2)
      do k = 1, lphi
        do i = 1, lrz
          do n = n_min, n_max
            n1 = 2*n + 1
            n2 = 2*n + 2
            aaa1(i,1,k) = aaa1(i,1,k) + aaa_n(1,i,1,n1)*cos_phi(n,k) &
                                      + aaa_n(1,i,1,n2)*sin_phi(n,k)
            aaa2(i,1,k) = aaa2(i,1,k) + aaa_n(2,i,1,n1)*cos_phi(n,k) &
                                      + aaa_n(2,i,1,n2)*sin_phi(n,k)
          end do
        end do
      end do
#endif
      
end
!--------------------------------------------------------------------
subroutine n1(n,aaa)
! filter to choose toroidal mode number n
! modified on 2012-04-17
!--------------------------------------------------------------------
      use check
      use mpiset
      implicit none
      integer::i,k,n,n_size

      real(8)::aaa(lr,lz,lphi)
      real(8)::work(2,lr,lz),work_total(2,lr,lz)

#ifdef AMDGPU
!$omp target teams distribute parallel do private(i,k)
#else      
!$omp parallel do private(i,k)
#endif
      do i = 1, lrz
          work(1,i,1) = 0.0d0
          work(2,i,1) = 0.0d0
        do k = lphistart, lphiend
          work(1,i,1) = work(1,i,1) + aaa(i,1,k)*cos_phi(n,k)
          work(2,i,1) = work(2,i,1) + aaa(i,1,k)*sin_phi(n,k)
        end do
      end do

      if(n.eq.0)then
#ifdef AMDGPU
!$omp target teams distribute parallel do private(i)
#else      
!$omp parallel do private(i)
#endif
        do i = 1, lrz
          work(1,i,1) = work(1,i,1)/dble(lphinet)
          work(2,i,1) = 0.0d0
        end do
      else
#ifdef AMDGPU
!$omp target teams distribute parallel do private(i)
#else      
!$omp parallel do private(i)
#endif
        do i = 1, lrz
          work(1,i,1) = work(1,i,1)/dble(lphinet)*2.0d0
          work(2,i,1) = work(2,i,1)/dble(lphinet)*2.0d0
        end do
      end if

      n_size = lrz*2
! reduction in phi direction using communicator 'phi_world'
      call mpi_allreduce(work(1,1,1),work_total(1,1,1),n_size&
           &,mpi_real8 ,mpi_sum,phi_world,mpi_err)

#ifdef AMDGPU
!$omp target teams distribute parallel do collapse(2) private(k,i)
#else
!$omp parallel do private(k,i)
#endif
      do k = 1, lphi
        do i = 1, lrz
          aaa(i,1,k) = work_total(1,i,1)*cos_phi(n,k) &
                     + work_total(2,i,1)*sin_phi(n,k)
        end do
      end do

end
!--------------------------------------------------------------------
subroutine data_snapshot(aaa,aaa_pol,kp)
! gather data on a ploidal plane at kp (local grid number) 
! to rank=0 of poloidal world
!--------------------------------------------------------------------
      use mpiset
      implicit none
      integer::kp,n,node_send
      integer::i,j,ig,jg

      real(8)::aaa(lr,lz,lphi),aaa_pol(lrnet,lznet)
      real(8)::wpol(lr,lz),work(lr,lz,0:mpi_proc_pol1)

      do j = 1, lz
      do i = 1, lr
        wpol(i,j) = aaa(i,j,kp)
      end do
      end do

      call mpi_gather(wpol(1,1),lrz,mpi_real8 ,work(1,1,0),lrz&
           &,mpi_real8 ,0,poloidal_world,mpi_err)

      if(my_rank_r.eq.0.and.my_rank_z.eq.0)then
!$omp parallel do private(node_send,ig,jg)
        do n = 0, mpi_proc_pol1
          node_send = n + my_rank_phi*mpi_proc_pol
          do j = kz_s(node_send), kz_e(node_send)
          do i = kr_s(node_send), kr_e(node_send)
            ig = i + kr_offset(node_send)
            jg = j + kz_offset(node_send)
            aaa_pol(ig,jg) = work(i,j,n)
          end do
          end do
        end do

      end if
end
!--------------------------------------------------------------------
subroutine analyz_s(scal,fc)
! scalar analysis on the flux coordinates
!--------------------------------------------------------------------
      use mpiset
      use grid
      use flux_c
      implicit none

      real(8)::scal(lr,lz,lphi)
      real(8)::fc(0:mpol,-ntor:ntor,lpsi,2)
      integer::l,m,n,nb,ia,ia1,ja,ja1,ka,mode_num,i1
      real(8)::c_norm,ar,ar1,az,az1
      real(8)::aaa1,aaa2,aaa3,aaa4,value

#ifdef AMDGPU      
!$omp target teams distribute parallel do collapse(3) private(l,n,m)
#else
!$omp parallel do
#endif
      do l = 1, lpsi
      do n =-ntor, ntor
      do m = 0, mpol
        fc(m,n,l,1) = 0.0d0
        fc(m,n,l,2) = 0.0d0
      end do
      end do
      end do

      c_norm = 2.0d0/dble(ltheta*lphinet) !corrected, 2015-04-28

! #ifdef AMDGPU      
! !$omp target teams distribute parallel do collapse(2) &
! !$omp private(n,m,ka,nb,l,ia,ia1,ja,ja1,ar1,ar,az1,az,aaa1,aaa2,aaa3,aaa4,value,i1)
! #endif

      do ka = lphistart, lphiend
      do nb = 1, n_flux_grid
        l = l_flux_grid(nb)

        ia = max(1, min(lr-1, &
                 int( ( r_flux_grid(nb)-(major_r - minor_r) )/dr) + 1 &
                      - kr_offset(my_rank) ) )
        ia1= ia + 1
        ja = max(1, min(lz-1, &
                 int( z_flux_grid(nb)/dz ) + 1 - kz_offset(my_rank) ) )
        ja1= ja + 1

        ar1 = (r_flux_grid(nb)-grr(ia,ja,ka) ) /dr
        ar  = 1.0d0 - ar1
        az1 = (z_flux_grid(nb)-gzz(ia,ja,ka) ) /dz
        az  = 1.0d0 - az1

        aaa1 = ar *az
        aaa2 = ar1*az
        aaa3 = ar *az1
        aaa4 = ar1*az1

! interpolated value

        value = scal(ia, ja, ka )*aaa1  + scal(ia1,ja, ka )*aaa2 &
              + scal(ia, ja1,ka )*aaa3  + scal(ia1,ja1,ka )*aaa4

        i1 = i_flux_grid(nb)

!$omp parallel do private(n,m) collapse(2)
        do n= -ntor,ntor
        do m=     0,mpol
          fc(m,n,l,1) = fc(m,n,l,1) + c_norm*value*(cos_theta_flux(m,i1)*cos_phi_flux(n,ka) &
                                                   -sin_theta_flux(m,i1)*sin_phi_flux(n,ka) )
          fc(m,n,l,2) = fc(m,n,l,2) + c_norm*value*(sin_theta_flux(m,i1)*cos_phi_flux(n,ka) &
                                                   +cos_theta_flux(m,i1)*sin_phi_flux(n,ka) )
        end do
        end do

      end do
      end do

! adjust for m=0 and n=0
! corrected, 2015-04-28

#ifdef AMDGPU      
!$omp target teams distribute parallel do private(l)
#else
!$omp parallel do
#endif
      do l = 1, lpsi
          fc(0,0,l,1) = fc(0,0,l,1)*0.50d0
          fc(0,0,l,2) = fc(0,0,l,2)*0.50d0
      end do


end
!--------------------------------------------------------------------
subroutine analyz_v(wr,wz,wphi,fr,ft,fp)
! vector analysis on the flux coordinates
!--------------------------------------------------------------------
      use mpiset
      use grid
      use flux_c
      implicit none

      real(8)::wr(lr,lz,lphi),wz(lr,lz,lphi),wphi(lr,lz,lphi)
      real(8)::fr(0:mpol,-ntor:ntor,lpsi,2),ft(0:mpol,-ntor:ntor,lpsi,2)
      real(8)::fp(0:mpol,-ntor:ntor,lpsi,2)
      integer::l,m,n,nb,ia,ia1,ja,ja1,ka,mode_num,i1
      real(8)::c_norm,ar,ar1,az,az1
      real(8)::aaa1,aaa2,aaa3,aaa4
      real(8)::value_r,value_z,value_phi,val_r,val_t,val_p
 
#ifdef AMDGPU
!$omp target teams distribute parallel do collapse(3) private(l,n,m)
#else
!$omp parallel do
#endif
      do l = 1, lpsi
      do n =-ntor, ntor
      do m = 0, mpol
        fr(m,n,l,1) = 0.0d0
        fr(m,n,l,2) = 0.0d0
        ft(m,n,l,1) = 0.0d0
        ft(m,n,l,2) = 0.0d0
        fp(m,n,l,1) = 0.0d0
        fp(m,n,l,2) = 0.0d0
      end do
      end do
      end do

      c_norm = 2.0d0/dble(ltheta*lphinet) !corrected, 2015-04-28

! #ifdef AMDGPU
! !$omp target teams distribute parallel do collapse(2) &
! !$omp private(ka,nb,l,ia,ia1,ja,ja1,ar1,ar,az1,az,aaa1,aaa2,aaa3,aaa4,value_r,value_z,value_phi,&
! !$omp val_r,val_t,val_p,i1,n,m)
! #endif
      do ka = lphistart, lphiend
      do nb = 1, n_flux_grid
        l = l_flux_grid(nb)

        ia = max(1, min(lr-1, &
                 int( ( r_flux_grid(nb)-(major_r - minor_r) )/dr) + 1 &
                      - kr_offset(my_rank) ) )
        ia1= ia + 1
        ja = max(1, min(lz-1, &
                 int( z_flux_grid(nb)/dz ) + 1 - kz_offset(my_rank) ) )
        ja1= ja + 1

        ar1 = (r_flux_grid(nb)-grr(ia,ja,ka) ) /dr
        ar  = 1.0d0 - ar1
        az1 = (z_flux_grid(nb)-gzz(ia,ja,ka) ) /dz
        az  = 1.0d0 - az1

        aaa1 = ar *az
        aaa2 = ar1*az
        aaa3 = ar *az1
        aaa4 = ar1*az1

! interpolated value

        value_r = wr(ia, ja, ka )*aaa1  + wr(ia1,ja, ka )*aaa2 &
                + wr(ia, ja1,ka )*aaa3  + wr(ia1,ja1,ka )*aaa4
        value_z = wz(ia, ja, ka )*aaa1  + wz(ia1,ja, ka )*aaa2 &
                + wz(ia, ja1,ka )*aaa3  + wz(ia1,ja1,ka )*aaa4
        value_phi = wphi(ia, ja, ka )*aaa1  + wphi(ia1,ja, ka )*aaa2 &
                  + wphi(ia, ja1,ka )*aaa3  + wphi(ia1,ja1,ka )*aaa4

        val_r = value_r*dr_r(nb) + value_z*dr_z(nb) + value_phi*dr_phi(nb)
        val_t = value_r*dt_r(nb) + value_z*dt_z(nb) + value_phi*dt_phi(nb)
        val_p = value_r*dp_r(nb) + value_z*dp_z(nb) + value_phi*dp_phi(nb)

        i1 = i_flux_grid(nb)

!$omp parallel do private(n,m) collapse(2)
        do n= -ntor,ntor
        do m=     0,mpol
          fr(m,n,l,1) = fr(m,n,l,1) + c_norm*val_r*(cos_theta_flux(m,i1)*cos_phi_flux(n,ka) &
                                                   -sin_theta_flux(m,i1)*sin_phi_flux(n,ka) )
          fr(m,n,l,2) = fr(m,n,l,2) + c_norm*val_r*(sin_theta_flux(m,i1)*cos_phi_flux(n,ka) &
                                                   +cos_theta_flux(m,i1)*sin_phi_flux(n,ka) )
          ft(m,n,l,1) = ft(m,n,l,1) + c_norm*val_t*(cos_theta_flux(m,i1)*cos_phi_flux(n,ka) &
                                                   -sin_theta_flux(m,i1)*sin_phi_flux(n,ka) )
          ft(m,n,l,2) = ft(m,n,l,2) + c_norm*val_t*(sin_theta_flux(m,i1)*cos_phi_flux(n,ka) &
                                                   +cos_theta_flux(m,i1)*sin_phi_flux(n,ka) )
          fp(m,n,l,1) = fp(m,n,l,1) + c_norm*val_p*(cos_theta_flux(m,i1)*cos_phi_flux(n,ka) &
                                                   -sin_theta_flux(m,i1)*sin_phi_flux(n,ka) )
          fp(m,n,l,2) = fp(m,n,l,2) + c_norm*val_p*(sin_theta_flux(m,i1)*cos_phi_flux(n,ka) &
                                                   +cos_theta_flux(m,i1)*sin_phi_flux(n,ka) )
        end do
        end do

      end do
      end do


! adjust for m=0 and n=0
! corrected, 2015-04-28

#ifdef AMDGPU
!$omp target teams distribute parallel do private(l)
#else
!$omp parallel do
#endif
      do l = 1, lpsi
          fr(0,0,l,1) = fr(0,0,l,1)*0.50d0
          fr(0,0,l,2) = fr(0,0,l,2)*0.50d0
          ft(0,0,l,1) = ft(0,0,l,1)*0.50d0
          ft(0,0,l,2) = ft(0,0,l,2)*0.50d0
          fp(0,0,l,1) = fp(0,0,l,1)*0.50d0
          fp(0,0,l,2) = fp(0,0,l,2)*0.50d0
      end do

end
!--------------------------------------------------------------------
subroutine periodic_field(xxx)
! boundary condition with mpi communications
! 2015-06-28
!--------------------------------------------------------------------
      use parameters
      use mpiset
      implicit none
      real(8)::xxx(lr,lz,lphi)
      integer::i,j,k
      integer::node_up,node_down,lrz2,lrphi2,lzphi2

      real(8)::up_send_z(lr,lz_shd,lphi),down_send_z(lr,lz_shd,lphi)
      real(8)::up_recv_z(lr,lz_shd,lphi),down_recv_z(lr,lz_shd,lphi)
      real(8)::up_send_r(lr_shd,lz,lphi),down_send_r(lr_shd,lz,lphi)
      real(8)::up_recv_r(lr_shd,lz,lphi),down_recv_r(lr_shd,lz,lphi)


! periodic in phi-direction

      node_up = mod(my_rank_phi + 1, mpi_proc_phi)
      node_down = mod(my_rank_phi - 1 + mpi_proc_phi, mpi_proc_phi)

! communication of boundary data

      lrz2 = lr*lz*lphi_shd

      call mpi_isend(xxx(1,1,lphiend-lphi_shd+1),lrz2,mpi_real8,node_up,1 &
                    ,phi_world,my_request(1),mpi_err)
      call mpi_isend(xxx(1,1,lphistart),lrz2,mpi_real8,node_down,1 &
                    ,phi_world,my_request(2),mpi_err)
      call mpi_irecv(xxx(1,1,lphistart-lphi_shd),lrz2,mpi_real8,node_down,1 &
                    ,phi_world,my_request(3),mpi_err)
      call mpi_irecv(xxx(1,1,lphiend+1),lrz2,mpi_real8,node_up,1 &
                    ,phi_world,my_request(4),mpi_err)

      call mpi_waitall(4,my_request,my_status_all,mpi_err)

!------------------------------------------------------
! communication in z-direction
      if(mpi_proc_z.ge.2)then

      node_up = mod(my_rank_z + 1, mpi_proc_z)
      node_down = mod(my_rank_z - 1 + mpi_proc_z, mpi_proc_z)

      if(my_rank_z.ne.(mpi_proc_z-1))then
#ifdef AMDGPU
!$omp target teams distribute parallel do collapse(3) private(k,j,i)
#else
!$omp parallel do
#endif
        do k = 1, lphi
        do j = 1, lz_shd
        do i = 1, lr
          up_send_z(i,j,k) = xxx(i,lzend-lz_shd+j,k)
        end do
        end do
        end do
      end if

      if(my_rank_z.ne.0 )then
#ifdef AMDGPU
!$omp target teams distribute parallel do collapse(3) private(k,j,i)
#else
!$omp parallel do
#endif
        do k = 1, lphi
        do j = 1, lz_shd
        do i = 1, lr
          down_send_z(i,j,k) = xxx(i,lzstart+j-1,  k)
        end do
        end do
        end do
      end if

! communication of particle data

      lrphi2 = lr*lphi*lz_shd

      if(my_rank_z.eq.0)then
        call mpi_isend(up_send_z(1,1,1),lrphi2,mpi_real8,node_up,1 &
                      ,z_world,my_request(1),mpi_err)
        call mpi_irecv(up_recv_z(1,1,1),lrphi2,mpi_real8,node_up,1 &
                      ,z_world,my_request(2),mpi_err)

        call mpi_waitall(2,my_request,my_status_all,mpi_err)

#ifdef AMDGPU
!$omp target teams distribute parallel do collapse(3) private(k,j,i)
#else
!$omp parallel do
#endif
        do k = 1, lphi
        do j = 1, lz_shd
        do i = 1, lr
          xxx(i,lzend+j,k) = up_recv_z(i,j,k)
        end do
        end do
        end do

      else if(my_rank_z.eq.(mpi_proc_z-1) )then
        call mpi_isend(down_send_z(1,1,1),lrphi2,mpi_real8,node_down,1 &
                      ,z_world,my_request(1),mpi_err)
        call mpi_irecv(down_recv_z(1,1,1),lrphi2,mpi_real8,node_down,1 &
                      ,z_world,my_request(2),mpi_err)

        call mpi_waitall(2,my_request,my_status_all,mpi_err)

#ifdef AMDGPU
!$omp target teams distribute parallel do collapse(3) private(k,j,i)
#else
!$omp parallel do
#endif
        do k = 1, lphi
        do j = 1, lz_shd
        do i = 1, lr
          xxx(i,lzstart-lz_shd-1+j,k) = down_recv_z(i,j,k)
        end do
        end do
        end do

      else
        call mpi_isend(up_send_z(1,1,1),lrphi2,mpi_real8,node_up,1 &
                      ,z_world,my_request(1),mpi_err)
        call mpi_isend(down_send_z(1,1,1),lrphi2,mpi_real8,node_down,1 &
                      ,z_world,my_request(2),mpi_err)
        call mpi_irecv(up_recv_z(1,1,1),lrphi2,mpi_real8,node_up,1 &
                      ,z_world,my_request(3),mpi_err)
        call mpi_irecv(down_recv_z(1,1,1),lrphi2,mpi_real8,node_down,1 &
                      ,z_world,my_request(4),mpi_err)

        call mpi_waitall(4,my_request,my_status_all,mpi_err)

#ifdef AMDGPU
!$omp target teams distribute parallel do collapse(3) private(k,j,i)
#else
!$omp parallel do
#endif
        do k = 1, lphi
        do j = 1, lz_shd
        do i = 1, lr
          xxx(i,lzstart-lz_shd-1+j,k) = down_recv_z(i,j,k)
          xxx(i,lzend+j,k) = up_recv_z(i,j,k)
        end do
        end do
        end do

      end if


      end if    !end of z-communication
!------------------------------------------------------
! communication in r-direction
      if(mpi_proc_r.ge.2)then

      node_up = mod(my_rank + 1, mpi_proc_r)
      node_down = mod(my_rank - 1 + mpi_proc_r, mpi_proc_r)

      if(my_rank_r.ne.(mpi_proc_r-1))then
#ifdef AMDGPU
!$omp target teams distribute parallel do collapse(2) private(j,i)
#else
!$omp parallel do
#endif
        do j = 1, lzphi
        do i = 1, lr_shd
          up_send_r(i,j,1) = xxx(lrend-lr_shd+i,j,1)
        end do
        end do
      end if

      if(my_rank_r.ne.0 )then
#ifdef AMDGPU
!$omp target teams distribute parallel do collapse(2) private(j,i)
#else
!$omp parallel do
#endif
        do j = 1, lzphi
        do i = 1, lr_shd
          down_send_r(i,j,1) = xxx(lrstart+i-1,  j,1)
        end do
        end do
      end if


! communication of particle data

      lzphi2 = lz*lphi*lr_shd

      if(my_rank_r.eq.0)then
        call mpi_isend(up_send_r(1,1,1),lzphi2,mpi_real8,node_up,1 &
                      ,r_world,my_request(1),mpi_err)
        call mpi_irecv(up_recv_r(1,1,1),lzphi2,mpi_real8,node_up,1 &
                      ,r_world,my_request(2),mpi_err)

        call mpi_waitall(2,my_request,my_status_all,mpi_err)

#ifdef AMDGPU
!$omp target teams distribute parallel do collapse(2) private(j,i)
#else
!$omp parallel do
#endif
        do j = 1, lzphi
        do i = 1, lr_shd
          xxx(lrend+i,j,1) = up_recv_r(i,j,1)
        end do
        end do

      else if(my_rank_r.eq.(mpi_proc_r-1) )then
        call mpi_isend(down_send_r(1,1,1),lzphi2,mpi_real8,node_down,1 &
                      ,r_world,my_request(1),mpi_err)
        call mpi_irecv(down_recv_r(1,1,1),lzphi2,mpi_real8,node_down,1 &
                      ,r_world,my_request(2),mpi_err)

        call mpi_waitall(2,my_request,my_status_all,mpi_err)

#ifdef AMDGPU
!$omp target teams distribute parallel do collapse(2) private(j,i)
#else
!$omp parallel do
#endif
        do j = 1, lzphi
        do i = 1, lr_shd
          xxx(lrstart-lr_shd-1+i,j,1) = down_recv_r(i,j,1)
        end do
        end do

      else
        call mpi_isend(up_send_r(1,1,1),lzphi2,mpi_real8,node_up,1 &
                      ,r_world,my_request(1),mpi_err)
        call mpi_isend(down_send_r(1,1,1),lzphi2,mpi_real8,node_down,1 &
                      ,r_world,my_request(2),mpi_err)
        call mpi_irecv(up_recv_r(1,1,1),lzphi2,mpi_real8,node_up,1 &
                      ,r_world,my_request(3),mpi_err)
        call mpi_irecv(down_recv_r(1,1,1),lzphi2,mpi_real8,node_down,1 &
                      ,r_world,my_request(4),mpi_err)

        call mpi_waitall(4,my_request,my_status_all,mpi_err)

#ifdef AMDGPU
!$omp target teams distribute parallel do collapse(2) private(j,i)
#else
!$omp parallel do
#endif
        do j = 1, lzphi
        do i = 1, lr_shd
          xxx(lrstart-lr_shd-1+i,j,1) = down_recv_r(i,j,1)
          xxx(lrend+i,j,1) = up_recv_r(i,j,1)
        end do
        end do

      end if

      end if  !end of r-communication

end
!--------------------------------------------------------------------
subroutine periodic_field_mlt2(xxx1,xxx2)
! boundary condition with mpi communications for 2 variables
! 2015-06-28
!--------------------------------------------------------------------
      use parameters
      use mpiset
      implicit none
      integer, parameter::mlt=2
      real(8)::xxx1(lr,lz,lphi),xxx2(lr,lz,lphi)
      integer::i,j,k
      integer::node_up,node_down,lrz2,lrphi2,lzphi2

      Real(8)::up_send_phi(lr,lz,lphi_shd,mlt),down_send_phi(lr,lz,lphi_shd,mlt)
      real(8)::up_recv_phi(lr,lz,lphi_shd,mlt),down_recv_phi(lr,lz,lphi_shd,mlt)
      real(8)::up_send_z(lr,lz_shd,lphi,mlt),down_send_z(lr,lz_shd,lphi,mlt)
      real(8)::up_recv_z(lr,lz_shd,lphi,mlt),down_recv_z(lr,lz_shd,lphi,mlt)
      real(8)::up_send_r(lr_shd,lz,lphi,mlt),down_send_r(lr_shd,lz,lphi,mlt)
      real(8)::up_recv_r(lr_shd,lz,lphi,mlt),down_recv_r(lr_shd,lz,lphi,mlt)


! periodic in phi-direction

      node_up = mod(my_rank_phi + 1, mpi_proc_phi)
      node_down = mod(my_rank_phi - 1 + mpi_proc_phi, mpi_proc_phi)

!$omp parallel do
      do k = 1, lphi_shd
      do i = 1, lrz
        down_send_phi(i,1,k,1) = xxx1(i,1,lphistart+k-1)
        up_send_phi(i,1,k,1) = xxx1(i,1,lphiend-lphi_shd+k)

        down_send_phi(i,1,k,2) = xxx2(i,1,lphistart+k-1)
        up_send_phi(i,1,k,2) = xxx2(i,1,lphiend-lphi_shd+k)
      end do
      end do

! communication of boundary data

      lrz2 = lr*lz*lphi_shd*mlt

      call mpi_isend(up_send_phi(1,1,1,1),lrz2,mpi_real8,node_up,1 &
                    ,phi_world,my_request(1),mpi_err)
      call mpi_isend(down_send_phi(1,1,1,1),lrz2,mpi_real8,node_down,1 &
                    ,phi_world,my_request(2),mpi_err)
      call mpi_irecv(down_recv_phi(1,1,1,1),lrz2,mpi_real8,node_down,1 &
                    ,phi_world,my_request(3),mpi_err)
      call mpi_irecv(up_recv_phi(1,1,1,1),lrz2,mpi_real8,node_up,1 &
                    ,phi_world,my_request(4),mpi_err)

      call mpi_waitall(4,my_request,my_status_all,mpi_err)

!$omp parallel do
      do k = 1, lphi_shd
      do i = 1, lrz
        xxx1(i,1,lphistart-lphi_shd-1+k) = down_recv_phi(i,1,k,1)
        xxx1(i,1,lphiend+k) = up_recv_phi(i,1,k,1)

        xxx2(i,1,lphistart-lphi_shd-1+k) = down_recv_phi(i,1,k,2)
        xxx2(i,1,lphiend+k) = up_recv_phi(i,1,k,2)
      end do
      end do

!------------------------------------------------------
! communication in z-direction
      if(mpi_proc_z.ge.2)then

      node_up = mod(my_rank_z + 1, mpi_proc_z)
      node_down = mod(my_rank_z - 1 + mpi_proc_z, mpi_proc_z)

      if(my_rank_z.ne.(mpi_proc_z-1))then
!$omp parallel do
        do k = 1, lphi
        do j = 1, lz_shd
        do i = 1, lr
          up_send_z(i,j,k,1) = xxx1(i,lzend-lz_shd+j,k)
          up_send_z(i,j,k,2) = xxx2(i,lzend-lz_shd+j,k)
        end do
        end do
        end do
      end if

      if(my_rank_z.ne.0 )then
!$omp parallel do
        do k = 1, lphi
        do j = 1, lz_shd
        do i = 1, lr
          down_send_z(i,j,k,1) = xxx1(i,lzstart+j-1,k)
          down_send_z(i,j,k,2) = xxx2(i,lzstart+j-1,k)
        end do
        end do
        end do
      end if

! communication of boundary data

      lrphi2 = lr*lphi*lz_shd*mlt

      if(my_rank_z.eq.0)then
        call mpi_isend(up_send_z(1,1,1,1),lrphi2,mpi_real8,node_up,1 &
                      ,z_world,my_request(1),mpi_err)
        call mpi_irecv(up_recv_z(1,1,1,1),lrphi2,mpi_real8,node_up,1 &
                      ,z_world,my_request(2),mpi_err)

        call mpi_waitall(2,my_request,my_status_all,mpi_err)

!$omp parallel do
        do k = 1, lphi
        do j = 1, lz_shd
        do i = 1, lr
          xxx1(i,lzend+j,k) = up_recv_z(i,j,k,1)
          xxx2(i,lzend+j,k) = up_recv_z(i,j,k,2)
        end do
        end do
        end do

      else if(my_rank_z.eq.(mpi_proc_z-1) )then
        call mpi_isend(down_send_z(1,1,1,1),lrphi2,mpi_real8,node_down,1 &
                      ,z_world,my_request(1),mpi_err)
        call mpi_irecv(down_recv_z(1,1,1,1),lrphi2,mpi_real8,node_down,1 &
                      ,z_world,my_request(2),mpi_err)

        call mpi_waitall(2,my_request,my_status_all,mpi_err)

!$omp parallel do
        do k = 1, lphi
        do j = 1, lz_shd
        do i = 1, lr
          xxx1(i,lzstart-lz_shd-1+j,k) = down_recv_z(i,j,k,1)
          xxx2(i,lzstart-lz_shd-1+j,k) = down_recv_z(i,j,k,2)
        end do
        end do
        end do

      else
        call mpi_isend(up_send_z(1,1,1,1),lrphi2,mpi_real8,node_up,1 &
                      ,z_world,my_request(1),mpi_err)
        call mpi_isend(down_send_z(1,1,1,1),lrphi2,mpi_real8,node_down,1 &
                      ,z_world,my_request(2),mpi_err)
        call mpi_irecv(up_recv_z(1,1,1,1),lrphi2,mpi_real8,node_up,1 &
                      ,z_world,my_request(3),mpi_err)
        call mpi_irecv(down_recv_z(1,1,1,1),lrphi2,mpi_real8,node_down,1 &
                      ,z_world,my_request(4),mpi_err)

        call mpi_waitall(4,my_request,my_status_all,mpi_err)

!$omp parallel do
        do k = 1, lphi
        do j = 1, lz_shd
        do i = 1, lr
          xxx1(i,lzstart-lz_shd-1+j,k) = down_recv_z(i,j,k,1)
          xxx1(i,lzend+j,k) = up_recv_z(i,j,k,1)

          xxx2(i,lzstart-lz_shd-1+j,k) = down_recv_z(i,j,k,2)
          xxx2(i,lzend+j,k) = up_recv_z(i,j,k,2)
        end do
        end do
        end do

      end if


      end if    !end of z-communication
!------------------------------------------------------
! communication in r-direction
      if(mpi_proc_r.ge.2)then

      node_up = mod(my_rank + 1, mpi_proc_r)
      node_down = mod(my_rank - 1 + mpi_proc_r, mpi_proc_r)

      if(my_rank_r.ne.(mpi_proc_r-1))then
!$omp parallel do
        do j = 1, lzphi
        do i = 1, lr_shd
          up_send_r(i,j,1,1) = xxx1(lrend-lr_shd+i,j,1)
          up_send_r(i,j,1,2) = xxx2(lrend-lr_shd+i,j,1)
        end do
        end do
      end if

      if(my_rank_r.ne.0 )then
!$omp parallel do
        do j = 1, lzphi
        do i = 1, lr_shd
          down_send_r(i,j,1,1) = xxx1(lrstart+i-1,  j,1)
          down_send_r(i,j,1,2) = xxx2(lrstart+i-1,  j,1)
        end do
        end do
      end if


! communication of boundary data

      lzphi2 = lz*lphi*lr_shd*mlt

      if(my_rank_r.eq.0)then
        call mpi_isend(up_send_r(1,1,1,1),lzphi2,mpi_real8,node_up,1 &
                      ,r_world,my_request(1),mpi_err)
        call mpi_irecv(up_recv_r(1,1,1,1),lzphi2,mpi_real8,node_up,1 &
                      ,r_world,my_request(2),mpi_err)

        call mpi_waitall(2,my_request,my_status_all,mpi_err)

!$omp parallel do
        do j = 1, lzphi
        do i = 1, lr_shd
          xxx1(lrend+i,j,1) = up_recv_r(i,j,1,1)
          xxx2(lrend+i,j,1) = up_recv_r(i,j,1,2)
        end do
        end do

      else if(my_rank_r.eq.(mpi_proc_r-1) )then
        call mpi_isend(down_send_r(1,1,1,1),lzphi2,mpi_real8,node_down,1 &
                      ,r_world,my_request(1),mpi_err)
        call mpi_irecv(down_recv_r(1,1,1,1),lzphi2,mpi_real8,node_down,1 &
                      ,r_world,my_request(2),mpi_err)

        call mpi_waitall(2,my_request,my_status_all,mpi_err)

!$omp parallel do
        do j = 1, lzphi
        do i = 1, lr_shd
          xxx1(lrstart-lr_shd-1+i,j,1) = down_recv_r(i,j,1,1)
          xxx2(lrstart-lr_shd-1+i,j,1) = down_recv_r(i,j,1,2)
        end do
        end do

      else
        call mpi_isend(up_send_r(1,1,1,1),lzphi2,mpi_real8,node_up,1 &
                      ,r_world,my_request(1),mpi_err)
        call mpi_isend(down_send_r(1,1,1,1),lzphi2,mpi_real8,node_down,1 &
                      ,r_world,my_request(2),mpi_err)
        call mpi_irecv(up_recv_r(1,1,1,1),lzphi2,mpi_real8,node_up,1 &
                      ,r_world,my_request(3),mpi_err)
        call mpi_irecv(down_recv_r(1,1,1,1),lzphi2,mpi_real8,node_down,1 &
                      ,r_world,my_request(4),mpi_err)

        call mpi_waitall(4,my_request,my_status_all,mpi_err)

!$omp parallel do
        do j = 1, lzphi
        do i = 1, lr_shd
          xxx1(lrstart-lr_shd-1+i,j,1) = down_recv_r(i,j,1,1)
          xxx1(lrend+i,j,1) = up_recv_r(i,j,1,1)

          xxx2(lrstart-lr_shd-1+i,j,1) = down_recv_r(i,j,1,2)
          xxx2(lrend+i,j,1) = up_recv_r(i,j,1,2)
        end do
        end do

      end if

      end if  !end of r-communication

end
!--------------------------------------------------------------------
subroutine periodic_field_mlt3(xxx1,xxx2,xxx3)
! boundary condition with mpi communications for 3 variables
! 2015-06-28
!--------------------------------------------------------------------
      use parameters
      use mpiset
      implicit none
      integer, parameter::mlt=3
      real(8)::xxx1(lr,lz,lphi),xxx2(lr,lz,lphi),xxx3(lr,lz,lphi)
      integer::i,j,k
      integer::node_up,node_down,lrz2,lrphi2,lzphi2

      Real(8)::up_send_phi(lr,lz,lphi_shd,mlt),down_send_phi(lr,lz,lphi_shd,mlt)
      real(8)::up_recv_phi(lr,lz,lphi_shd,mlt),down_recv_phi(lr,lz,lphi_shd,mlt)
      real(8)::up_send_z(lr,lz_shd,lphi,mlt),down_send_z(lr,lz_shd,lphi,mlt)
      real(8)::up_recv_z(lr,lz_shd,lphi,mlt),down_recv_z(lr,lz_shd,lphi,mlt)
      real(8)::up_send_r(lr_shd,lz,lphi,mlt),down_send_r(lr_shd,lz,lphi,mlt)
      real(8)::up_recv_r(lr_shd,lz,lphi,mlt),down_recv_r(lr_shd,lz,lphi,mlt)


! periodic in phi-direction

      node_up = mod(my_rank_phi + 1, mpi_proc_phi)
      node_down = mod(my_rank_phi - 1 + mpi_proc_phi, mpi_proc_phi)

#ifdef AMDGPU
!$omp target teams distribute parallel do collapse(2) private(k,i)
#else
!$omp parallel do
#endif
      do k = 1, lphi_shd
      do i = 1, lrz
        down_send_phi(i,1,k,1) = xxx1(i,1,lphistart+k-1)
        up_send_phi(i,1,k,1) = xxx1(i,1,lphiend-lphi_shd+k)

        down_send_phi(i,1,k,2) = xxx2(i,1,lphistart+k-1)
        up_send_phi(i,1,k,2) = xxx2(i,1,lphiend-lphi_shd+k)

        down_send_phi(i,1,k,3) = xxx3(i,1,lphistart+k-1)
        up_send_phi(i,1,k,3) = xxx3(i,1,lphiend-lphi_shd+k)
      end do
      end do

! communication of boundary data

      lrz2 = lr*lz*lphi_shd*mlt

      call mpi_isend(up_send_phi(1,1,1,1),lrz2,mpi_real8,node_up,1 &
                    ,phi_world,my_request(1),mpi_err)
      call mpi_isend(down_send_phi(1,1,1,1),lrz2,mpi_real8,node_down,1 &
                    ,phi_world,my_request(2),mpi_err)
      call mpi_irecv(down_recv_phi(1,1,1,1),lrz2,mpi_real8,node_down,1 &
                    ,phi_world,my_request(3),mpi_err)
      call mpi_irecv(up_recv_phi(1,1,1,1),lrz2,mpi_real8,node_up,1 &
                    ,phi_world,my_request(4),mpi_err)

      call mpi_waitall(4,my_request,my_status_all,mpi_err)

#ifdef AMDGPU
!$omp target teams distribute parallel do collapse(2) private(k,i)
#else
!$omp parallel do
#endif
      do k = 1, lphi_shd
      do i = 1, lrz
        xxx1(i,1,lphistart-lphi_shd-1+k) = down_recv_phi(i,1,k,1)
        xxx1(i,1,lphiend+k) = up_recv_phi(i,1,k,1)

        xxx2(i,1,lphistart-lphi_shd-1+k) = down_recv_phi(i,1,k,2)
        xxx2(i,1,lphiend+k) = up_recv_phi(i,1,k,2)

        xxx3(i,1,lphistart-lphi_shd-1+k) = down_recv_phi(i,1,k,3)
        xxx3(i,1,lphiend+k) = up_recv_phi(i,1,k,3)
      end do
      end do

!------------------------------------------------------
! communication in z-direction
      if(mpi_proc_z.ge.2)then

      node_up = mod(my_rank_z + 1, mpi_proc_z)
      node_down = mod(my_rank_z - 1 + mpi_proc_z, mpi_proc_z)

      if(my_rank_z.ne.(mpi_proc_z-1))then
#ifdef AMDGPU
!$omp target teams distribute parallel do collapse(3) private(k,j,i)
#else
!$omp parallel do
#endif
        do k = 1, lphi
        do j = 1, lz_shd
        do i = 1, lr
          up_send_z(i,j,k,1) = xxx1(i,lzend-lz_shd+j,k)
          up_send_z(i,j,k,2) = xxx2(i,lzend-lz_shd+j,k)
          up_send_z(i,j,k,3) = xxx3(i,lzend-lz_shd+j,k)
        end do
        end do
        end do
      end if

      if(my_rank_z.ne.0 )then
#ifdef AMDGPU
!$omp target teams distribute parallel do collapse(3) private(k,j,i)
#else
!$omp parallel do
#endif
        do k = 1, lphi
        do j = 1, lz_shd
        do i = 1, lr
          down_send_z(i,j,k,1) = xxx1(i,lzstart+j-1,k)
          down_send_z(i,j,k,2) = xxx2(i,lzstart+j-1,k)
          down_send_z(i,j,k,3) = xxx3(i,lzstart+j-1,k)
        end do
        end do
        end do
      end if

! communication of boundary data

      lrphi2 = lr*lphi*lz_shd*mlt

      if(my_rank_z.eq.0)then
        call mpi_isend(up_send_z(1,1,1,1),lrphi2,mpi_real8,node_up,1 &
                      ,z_world,my_request(1),mpi_err)
        call mpi_irecv(up_recv_z(1,1,1,1),lrphi2,mpi_real8,node_up,1 &
                      ,z_world,my_request(2),mpi_err)

        call mpi_waitall(2,my_request,my_status_all,mpi_err)

#ifdef AMDGPU
!$omp target teams distribute parallel do collapse(3) private(k,j,i)
#else
!$omp parallel do
#endif
        do k = 1, lphi
        do j = 1, lz_shd
        do i = 1, lr
          xxx1(i,lzend+j,k) = up_recv_z(i,j,k,1)
          xxx2(i,lzend+j,k) = up_recv_z(i,j,k,2)
          xxx3(i,lzend+j,k) = up_recv_z(i,j,k,3)
        end do
        end do
        end do

      else if(my_rank_z.eq.(mpi_proc_z-1) )then
        call mpi_isend(down_send_z(1,1,1,1),lrphi2,mpi_real8,node_down,1 &
                      ,z_world,my_request(1),mpi_err)
        call mpi_irecv(down_recv_z(1,1,1,1),lrphi2,mpi_real8,node_down,1 &
                      ,z_world,my_request(2),mpi_err)

        call mpi_waitall(2,my_request,my_status_all,mpi_err)

#ifdef AMDGPU
!$omp target teams distribute parallel do collapse(3) private(k,j,i)
#else
!$omp parallel do
#endif
        do k = 1, lphi
        do j = 1, lz_shd
        do i = 1, lr
          xxx1(i,lzstart-lz_shd-1+j,k) = down_recv_z(i,j,k,1)
          xxx2(i,lzstart-lz_shd-1+j,k) = down_recv_z(i,j,k,2)
          xxx3(i,lzstart-lz_shd-1+j,k) = down_recv_z(i,j,k,3)
        end do
        end do
        end do

      else
        call mpi_isend(up_send_z(1,1,1,1),lrphi2,mpi_real8,node_up,1 &
                      ,z_world,my_request(1),mpi_err)
        call mpi_isend(down_send_z(1,1,1,1),lrphi2,mpi_real8,node_down,1 &
                      ,z_world,my_request(2),mpi_err)
        call mpi_irecv(up_recv_z(1,1,1,1),lrphi2,mpi_real8,node_up,1 &
                      ,z_world,my_request(3),mpi_err)
        call mpi_irecv(down_recv_z(1,1,1,1),lrphi2,mpi_real8,node_down,1 &
                      ,z_world,my_request(4),mpi_err)

        call mpi_waitall(4,my_request,my_status_all,mpi_err)

#ifdef AMDGPU
!$omp target teams distribute parallel do collapse(3) private(k,j,i)
#else
!$omp parallel do
#endif
        do k = 1, lphi
        do j = 1, lz_shd
        do i = 1, lr
          xxx1(i,lzstart-lz_shd-1+j,k) = down_recv_z(i,j,k,1)
          xxx1(i,lzend+j,k) = up_recv_z(i,j,k,1)

          xxx2(i,lzstart-lz_shd-1+j,k) = down_recv_z(i,j,k,2)
          xxx2(i,lzend+j,k) = up_recv_z(i,j,k,2)

          xxx3(i,lzstart-lz_shd-1+j,k) = down_recv_z(i,j,k,3)
          xxx3(i,lzend+j,k) = up_recv_z(i,j,k,3)
        end do
        end do
        end do

      end if


      end if    !end of z-communication
!------------------------------------------------------
! communication in r-direction
      if(mpi_proc_r.ge.2)then

      node_up = mod(my_rank + 1, mpi_proc_r)
      node_down = mod(my_rank - 1 + mpi_proc_r, mpi_proc_r)

      if(my_rank_r.ne.(mpi_proc_r-1))then
#ifdef AMDGPU
!$omp target teams distribute parallel do collapse(2) private(j,i)
#else
!$omp parallel do
#endif
        do j = 1, lzphi
        do i = 1, lr_shd
          up_send_r(i,j,1,1) = xxx1(lrend-lr_shd+i,j,1)
          up_send_r(i,j,1,2) = xxx2(lrend-lr_shd+i,j,1)
          up_send_r(i,j,1,3) = xxx3(lrend-lr_shd+i,j,1)
        end do
        end do
      end if

      if(my_rank_r.ne.0 )then
#ifdef AMDGPU
!$omp target teams distribute parallel do collapse(2) private(j,i)
#else
!$omp parallel do
#endif
        do j = 1, lzphi
        do i = 1, lr_shd
          down_send_r(i,j,1,1) = xxx1(lrstart+i-1,  j,1)
          down_send_r(i,j,1,2) = xxx2(lrstart+i-1,  j,1)
          down_send_r(i,j,1,3) = xxx3(lrstart+i-1,  j,1)
        end do
        end do
      end if


! communication of boundary data

      lzphi2 = lz*lphi*lr_shd*mlt

      if(my_rank_r.eq.0)then
        call mpi_isend(up_send_r(1,1,1,1),lzphi2,mpi_real8,node_up,1 &
                      ,r_world,my_request(1),mpi_err)
        call mpi_irecv(up_recv_r(1,1,1,1),lzphi2,mpi_real8,node_up,1 &
                      ,r_world,my_request(2),mpi_err)

        call mpi_waitall(2,my_request,my_status_all,mpi_err)

#ifdef AMDGPU
!$omp target teams distribute parallel do collapse(2) private(j,i)
#else
!$omp parallel do
#endif
        do j = 1, lzphi
        do i = 1, lr_shd
          xxx1(lrend+i,j,1) = up_recv_r(i,j,1,1)
          xxx2(lrend+i,j,1) = up_recv_r(i,j,1,2)
          xxx3(lrend+i,j,1) = up_recv_r(i,j,1,3)
        end do
        end do

      else if(my_rank_r.eq.(mpi_proc_r-1) )then
        call mpi_isend(down_send_r(1,1,1,1),lzphi2,mpi_real8,node_down,1 &
                      ,r_world,my_request(1),mpi_err)
        call mpi_irecv(down_recv_r(1,1,1,1),lzphi2,mpi_real8,node_down,1 &
                      ,r_world,my_request(2),mpi_err)

        call mpi_waitall(2,my_request,my_status_all,mpi_err)

#ifdef AMDGPU
!$omp target teams distribute parallel do collapse(2) private(j,i)
#else
!$omp parallel do
#endif
        do j = 1, lzphi
        do i = 1, lr_shd
          xxx1(lrstart-lr_shd-1+i,j,1) = down_recv_r(i,j,1,1)
          xxx2(lrstart-lr_shd-1+i,j,1) = down_recv_r(i,j,1,2)
          xxx3(lrstart-lr_shd-1+i,j,1) = down_recv_r(i,j,1,3)
        end do
        end do

      else
        call mpi_isend(up_send_r(1,1,1,1),lzphi2,mpi_real8,node_up,1 &
                      ,r_world,my_request(1),mpi_err)
        call mpi_isend(down_send_r(1,1,1,1),lzphi2,mpi_real8,node_down,1 &
                      ,r_world,my_request(2),mpi_err)
        call mpi_irecv(up_recv_r(1,1,1,1),lzphi2,mpi_real8,node_up,1 &
                      ,r_world,my_request(3),mpi_err)
        call mpi_irecv(down_recv_r(1,1,1,1),lzphi2,mpi_real8,node_down,1 &
                      ,r_world,my_request(4),mpi_err)

        call mpi_waitall(4,my_request,my_status_all,mpi_err)

#ifdef AMDGPU
!$omp target teams distribute parallel do collapse(2) private(j,i)
#else
!$omp parallel do
#endif
        do j = 1, lzphi
        do i = 1, lr_shd
          xxx1(lrstart-lr_shd-1+i,j,1) = down_recv_r(i,j,1,1)
          xxx1(lrend+i,j,1) = up_recv_r(i,j,1,1)

          xxx2(lrstart-lr_shd-1+i,j,1) = down_recv_r(i,j,1,2)
          xxx2(lrend+i,j,1) = up_recv_r(i,j,1,2)

          xxx3(lrstart-lr_shd-1+i,j,1) = down_recv_r(i,j,1,3)
          xxx3(lrend+i,j,1) = up_recv_r(i,j,1,3)
        end do
        end do

      end if

      end if  !end of r-communication

end
!--------------------------------------------------------------------
subroutine periodic_field_mlt7(xxx1,xxx2,xxx3,xxx4,xxx5,xxx6,xxx7)
! boundary condition with mpi communications for 7 variables
! 2015-06-28
!--------------------------------------------------------------------
      use parameters
      use mpiset
      implicit none
      integer, parameter::mlt=7
      real(8)::xxx1(lr,lz,lphi),xxx2(lr,lz,lphi),xxx3(lr,lz,lphi),xxx4(lr,lz,lphi)
      real(8)::xxx5(lr,lz,lphi),xxx6(lr,lz,lphi),xxx7(lr,lz,lphi)
      integer::i,j,k
      integer::node_up,node_down,lrz2,lrphi2,lzphi2

      Real(8)::up_send_phi(lr,lz,lphi_shd,mlt),down_send_phi(lr,lz,lphi_shd,mlt)
      real(8)::up_recv_phi(lr,lz,lphi_shd,mlt),down_recv_phi(lr,lz,lphi_shd,mlt)
      real(8)::up_send_z(lr,lz_shd,lphi,mlt),down_send_z(lr,lz_shd,lphi,mlt)
      real(8)::up_recv_z(lr,lz_shd,lphi,mlt),down_recv_z(lr,lz_shd,lphi,mlt)
      real(8)::up_send_r(lr_shd,lz,lphi,mlt),down_send_r(lr_shd,lz,lphi,mlt)
      real(8)::up_recv_r(lr_shd,lz,lphi,mlt),down_recv_r(lr_shd,lz,lphi,mlt)


! periodic in phi-direction

      node_up = mod(my_rank_phi + 1, mpi_proc_phi)
      node_down = mod(my_rank_phi - 1 + mpi_proc_phi, mpi_proc_phi)

#ifdef AMDGPU      
!$omp target teams distribute parallel do collapse(2) private(k,i)
#else
!$omp parallel do
#endif
      do k = 1, lphi_shd
      do i = 1, lrz
        down_send_phi(i,1,k,1) = xxx1(i,1,lphistart+k-1)
        up_send_phi(i,1,k,1) = xxx1(i,1,lphiend-lphi_shd+k)

        down_send_phi(i,1,k,2) = xxx2(i,1,lphistart+k-1)
        up_send_phi(i,1,k,2) = xxx2(i,1,lphiend-lphi_shd+k)

        down_send_phi(i,1,k,3) = xxx3(i,1,lphistart+k-1)
        up_send_phi(i,1,k,3) = xxx3(i,1,lphiend-lphi_shd+k)

        down_send_phi(i,1,k,4) = xxx4(i,1,lphistart+k-1)
        up_send_phi(i,1,k,4) = xxx4(i,1,lphiend-lphi_shd+k)

        down_send_phi(i,1,k,5) = xxx5(i,1,lphistart+k-1)
        up_send_phi(i,1,k,5) = xxx5(i,1,lphiend-lphi_shd+k)

        down_send_phi(i,1,k,6) = xxx6(i,1,lphistart+k-1)
        up_send_phi(i,1,k,6) = xxx6(i,1,lphiend-lphi_shd+k)

        down_send_phi(i,1,k,7) = xxx7(i,1,lphistart+k-1)
        up_send_phi(i,1,k,7) = xxx7(i,1,lphiend-lphi_shd+k)
      end do
      end do

! communication of boundary data

      lrz2 = lr*lz*lphi_shd*mlt

      call mpi_isend(up_send_phi(1,1,1,1),lrz2,mpi_real8,node_up,1 &
                    ,phi_world,my_request(1),mpi_err)
      call mpi_isend(down_send_phi(1,1,1,1),lrz2,mpi_real8,node_down,1 &
                    ,phi_world,my_request(2),mpi_err)
      call mpi_irecv(down_recv_phi(1,1,1,1),lrz2,mpi_real8,node_down,1 &
                    ,phi_world,my_request(3),mpi_err)
      call mpi_irecv(up_recv_phi(1,1,1,1),lrz2,mpi_real8,node_up,1 &
                    ,phi_world,my_request(4),mpi_err)

      call mpi_waitall(4,my_request,my_status_all,mpi_err)

#ifdef AMDGPU      
!$omp target teams distribute parallel do collapse(2) private(k,i)
#else
!$omp parallel do
#endif
      do k = 1, lphi_shd
      do i = 1, lrz
        xxx1(i,1,lphistart-lphi_shd-1+k) = down_recv_phi(i,1,k,1)
        xxx1(i,1,lphiend+k) = up_recv_phi(i,1,k,1)

        xxx2(i,1,lphistart-lphi_shd-1+k) = down_recv_phi(i,1,k,2)
        xxx2(i,1,lphiend+k) = up_recv_phi(i,1,k,2)

        xxx3(i,1,lphistart-lphi_shd-1+k) = down_recv_phi(i,1,k,3)
        xxx3(i,1,lphiend+k) = up_recv_phi(i,1,k,3)

        xxx4(i,1,lphistart-lphi_shd-1+k) = down_recv_phi(i,1,k,4)
        xxx4(i,1,lphiend+k) = up_recv_phi(i,1,k,4)

        xxx5(i,1,lphistart-lphi_shd-1+k) = down_recv_phi(i,1,k,5)
        xxx5(i,1,lphiend+k) = up_recv_phi(i,1,k,5)

        xxx6(i,1,lphistart-lphi_shd-1+k) = down_recv_phi(i,1,k,6)
        xxx6(i,1,lphiend+k) = up_recv_phi(i,1,k,6)

        xxx7(i,1,lphistart-lphi_shd-1+k) = down_recv_phi(i,1,k,7)
        xxx7(i,1,lphiend+k) = up_recv_phi(i,1,k,7)
      end do
      end do

!------------------------------------------------------
! communication in z-direction
      if(mpi_proc_z.ge.2)then

      node_up = mod(my_rank_z + 1, mpi_proc_z)
      node_down = mod(my_rank_z - 1 + mpi_proc_z, mpi_proc_z)

      if(my_rank_z.ne.(mpi_proc_z-1))then
#ifdef AMDGPU      
!$omp target teams distribute parallel do collapse(3) private(k,j,i)
#else
!$omp parallel do
#endif
         do k = 1, lphi
        do j = 1, lz_shd
        do i = 1, lr
          up_send_z(i,j,k,1) = xxx1(i,lzend-lz_shd+j,k)
          up_send_z(i,j,k,2) = xxx2(i,lzend-lz_shd+j,k)
          up_send_z(i,j,k,3) = xxx3(i,lzend-lz_shd+j,k)
          up_send_z(i,j,k,4) = xxx4(i,lzend-lz_shd+j,k)
          up_send_z(i,j,k,5) = xxx5(i,lzend-lz_shd+j,k)
          up_send_z(i,j,k,6) = xxx6(i,lzend-lz_shd+j,k)
          up_send_z(i,j,k,7) = xxx7(i,lzend-lz_shd+j,k)
        end do
        end do
        end do
      end if

      if(my_rank_z.ne.0 )then
#ifdef AMDGPU      
!$omp target teams distribute parallel do collapse(3) private(k,j,i)
#else
!$omp parallel do
#endif
        do k = 1, lphi
        do j = 1, lz_shd
        do i = 1, lr
          down_send_z(i,j,k,1) = xxx1(i,lzstart+j-1,k)
          down_send_z(i,j,k,2) = xxx2(i,lzstart+j-1,k)
          down_send_z(i,j,k,3) = xxx3(i,lzstart+j-1,k)
          down_send_z(i,j,k,4) = xxx4(i,lzstart+j-1,k)
          down_send_z(i,j,k,5) = xxx5(i,lzstart+j-1,k)
          down_send_z(i,j,k,6) = xxx6(i,lzstart+j-1,k)
          down_send_z(i,j,k,7) = xxx7(i,lzstart+j-1,k)
        end do
        end do
        end do
      end if

! communication of boundary data

      lrphi2 = lr*lphi*lz_shd*mlt

      if(my_rank_z.eq.0)then
        call mpi_isend(up_send_z(1,1,1,1),lrphi2,mpi_real8,node_up,1 &
                      ,z_world,my_request(1),mpi_err)
        call mpi_irecv(up_recv_z(1,1,1,1),lrphi2,mpi_real8,node_up,1 &
                      ,z_world,my_request(2),mpi_err)

        call mpi_waitall(2,my_request,my_status_all,mpi_err)

#ifdef AMDGPU      
!$omp target teams distribute parallel do collapse(3) private(k,j,i)
#else
!$omp parallel do
#endif
        do k = 1, lphi
        do j = 1, lz_shd
        do i = 1, lr
          xxx1(i,lzend+j,k) = up_recv_z(i,j,k,1)
          xxx2(i,lzend+j,k) = up_recv_z(i,j,k,2)
          xxx3(i,lzend+j,k) = up_recv_z(i,j,k,3)
          xxx4(i,lzend+j,k) = up_recv_z(i,j,k,4)
          xxx5(i,lzend+j,k) = up_recv_z(i,j,k,5)
          xxx6(i,lzend+j,k) = up_recv_z(i,j,k,6)
          xxx7(i,lzend+j,k) = up_recv_z(i,j,k,7)
        end do
        end do
        end do

      else if(my_rank_z.eq.(mpi_proc_z-1) )then
        call mpi_isend(down_send_z(1,1,1,1),lrphi2,mpi_real8,node_down,1 &
                      ,z_world,my_request(1),mpi_err)
        call mpi_irecv(down_recv_z(1,1,1,1),lrphi2,mpi_real8,node_down,1 &
                      ,z_world,my_request(2),mpi_err)

        call mpi_waitall(2,my_request,my_status_all,mpi_err)

#ifdef AMDGPU      
!$omp target teams distribute parallel do collapse(3) private(k,j,i)
#else
!$omp parallel do
#endif
        do k = 1, lphi
        do j = 1, lz_shd
        do i = 1, lr
          xxx1(i,lzstart-lz_shd-1+j,k) = down_recv_z(i,j,k,1)
          xxx2(i,lzstart-lz_shd-1+j,k) = down_recv_z(i,j,k,2)
          xxx3(i,lzstart-lz_shd-1+j,k) = down_recv_z(i,j,k,3)
          xxx4(i,lzstart-lz_shd-1+j,k) = down_recv_z(i,j,k,4)
          xxx5(i,lzstart-lz_shd-1+j,k) = down_recv_z(i,j,k,5)
          xxx6(i,lzstart-lz_shd-1+j,k) = down_recv_z(i,j,k,6)
          xxx7(i,lzstart-lz_shd-1+j,k) = down_recv_z(i,j,k,7)
        end do
        end do
        end do

      else
        call mpi_isend(up_send_z(1,1,1,1),lrphi2,mpi_real8,node_up,1 &
                      ,z_world,my_request(1),mpi_err)
        call mpi_isend(down_send_z(1,1,1,1),lrphi2,mpi_real8,node_down,1 &
                      ,z_world,my_request(2),mpi_err)
        call mpi_irecv(up_recv_z(1,1,1,1),lrphi2,mpi_real8,node_up,1 &
                      ,z_world,my_request(3),mpi_err)
        call mpi_irecv(down_recv_z(1,1,1,1),lrphi2,mpi_real8,node_down,1 &
                      ,z_world,my_request(4),mpi_err)

        call mpi_waitall(4,my_request,my_status_all,mpi_err)

#ifdef AMDGPU      
!$omp target teams distribute parallel do collapse(3) private(k,j,i)
#else
!$omp parallel do
#endif
        do k = 1, lphi
        do j = 1, lz_shd
        do i = 1, lr
          xxx1(i,lzstart-lz_shd-1+j,k) = down_recv_z(i,j,k,1)
          xxx1(i,lzend+j,k) = up_recv_z(i,j,k,1)

          xxx2(i,lzstart-lz_shd-1+j,k) = down_recv_z(i,j,k,2)
          xxx2(i,lzend+j,k) = up_recv_z(i,j,k,2)

          xxx3(i,lzstart-lz_shd-1+j,k) = down_recv_z(i,j,k,3)
          xxx3(i,lzend+j,k) = up_recv_z(i,j,k,3)

          xxx4(i,lzstart-lz_shd-1+j,k) = down_recv_z(i,j,k,4)
          xxx4(i,lzend+j,k) = up_recv_z(i,j,k,4)

          xxx5(i,lzstart-lz_shd-1+j,k) = down_recv_z(i,j,k,5)
          xxx5(i,lzend+j,k) = up_recv_z(i,j,k,5)

          xxx6(i,lzstart-lz_shd-1+j,k) = down_recv_z(i,j,k,6)
          xxx6(i,lzend+j,k) = up_recv_z(i,j,k,6)

          xxx7(i,lzstart-lz_shd-1+j,k) = down_recv_z(i,j,k,7)
          xxx7(i,lzend+j,k) = up_recv_z(i,j,k,7)
        end do
        end do
        end do

      end if


      end if    !end of z-communication
!------------------------------------------------------
! communication in r-direction
      if(mpi_proc_r.ge.2)then

      node_up = mod(my_rank + 1, mpi_proc_r)
      node_down = mod(my_rank - 1 + mpi_proc_r, mpi_proc_r)

      if(my_rank_r.ne.(mpi_proc_r-1))then
#ifdef AMDGPU      
!$omp target teams distribute parallel do collapse(2) private(j,i)
#else
!$omp parallel do
#endif
        do j = 1, lzphi
        do i = 1, lr_shd
          up_send_r(i,j,1,1) = xxx1(lrend-lr_shd+i,j,1)
          up_send_r(i,j,1,2) = xxx2(lrend-lr_shd+i,j,1)
          up_send_r(i,j,1,3) = xxx3(lrend-lr_shd+i,j,1)
          up_send_r(i,j,1,4) = xxx4(lrend-lr_shd+i,j,1)
          up_send_r(i,j,1,5) = xxx5(lrend-lr_shd+i,j,1)
          up_send_r(i,j,1,6) = xxx6(lrend-lr_shd+i,j,1)
          up_send_r(i,j,1,7) = xxx7(lrend-lr_shd+i,j,1)
        end do
        end do
      end if

      if(my_rank_r.ne.0 )then
#ifdef AMDGPU      
!$omp target teams distribute parallel do collapse(2) private(j,i)
#else
!$omp parallel do
#endif
        do j = 1, lzphi
        do i = 1, lr_shd
          down_send_r(i,j,1,1) = xxx1(lrstart+i-1,  j,1)
          down_send_r(i,j,1,2) = xxx2(lrstart+i-1,  j,1)
          down_send_r(i,j,1,3) = xxx3(lrstart+i-1,  j,1)
          down_send_r(i,j,1,4) = xxx4(lrstart+i-1,  j,1)
          down_send_r(i,j,1,5) = xxx5(lrstart+i-1,  j,1)
          down_send_r(i,j,1,6) = xxx6(lrstart+i-1,  j,1)
          down_send_r(i,j,1,7) = xxx7(lrstart+i-1,  j,1)
        end do
        end do
      end if


! communication of boundary data

      lzphi2 = lz*lphi*lr_shd*mlt

      if(my_rank_r.eq.0)then
        call mpi_isend(up_send_r(1,1,1,1),lzphi2,mpi_real8,node_up,1 &
                      ,r_world,my_request(1),mpi_err)
        call mpi_irecv(up_recv_r(1,1,1,1),lzphi2,mpi_real8,node_up,1 &
                      ,r_world,my_request(2),mpi_err)

        call mpi_waitall(2,my_request,my_status_all,mpi_err)

#ifdef AMDGPU      
!$omp target teams distribute parallel do collapse(2) private(j,i)
#else
!$omp parallel do
#endif
        do j = 1, lzphi
        do i = 1, lr_shd
          xxx1(lrend+i,j,1) = up_recv_r(i,j,1,1)
          xxx2(lrend+i,j,1) = up_recv_r(i,j,1,2)
          xxx3(lrend+i,j,1) = up_recv_r(i,j,1,3)
          xxx4(lrend+i,j,1) = up_recv_r(i,j,1,4)
          xxx5(lrend+i,j,1) = up_recv_r(i,j,1,5)
          xxx6(lrend+i,j,1) = up_recv_r(i,j,1,6)
          xxx7(lrend+i,j,1) = up_recv_r(i,j,1,7)
        end do
        end do

      else if(my_rank_r.eq.(mpi_proc_r-1) )then
        call mpi_isend(down_send_r(1,1,1,1),lzphi2,mpi_real8,node_down,1 &
                      ,r_world,my_request(1),mpi_err)
        call mpi_irecv(down_recv_r(1,1,1,1),lzphi2,mpi_real8,node_down,1 &
                      ,r_world,my_request(2),mpi_err)

        call mpi_waitall(2,my_request,my_status_all,mpi_err)

#ifdef AMDGPU      
!$omp target teams distribute parallel do collapse(2) private(j,i)
#else
!$omp parallel do
#endif
        do j = 1, lzphi
        do i = 1, lr_shd
          xxx1(lrstart-lr_shd-1+i,j,1) = down_recv_r(i,j,1,1)
          xxx2(lrstart-lr_shd-1+i,j,1) = down_recv_r(i,j,1,2)
          xxx3(lrstart-lr_shd-1+i,j,1) = down_recv_r(i,j,1,3)
          xxx4(lrstart-lr_shd-1+i,j,1) = down_recv_r(i,j,1,4)
          xxx5(lrstart-lr_shd-1+i,j,1) = down_recv_r(i,j,1,5)
          xxx6(lrstart-lr_shd-1+i,j,1) = down_recv_r(i,j,1,6)
          xxx7(lrstart-lr_shd-1+i,j,1) = down_recv_r(i,j,1,7)
        end do
        end do

      else
        call mpi_isend(up_send_r(1,1,1,1),lzphi2,mpi_real8,node_up,1 &
                      ,r_world,my_request(1),mpi_err)
        call mpi_isend(down_send_r(1,1,1,1),lzphi2,mpi_real8,node_down,1 &
                      ,r_world,my_request(2),mpi_err)
        call mpi_irecv(up_recv_r(1,1,1,1),lzphi2,mpi_real8,node_up,1 &
                      ,r_world,my_request(3),mpi_err)
        call mpi_irecv(down_recv_r(1,1,1,1),lzphi2,mpi_real8,node_down,1 &
                      ,r_world,my_request(4),mpi_err)

        call mpi_waitall(4,my_request,my_status_all,mpi_err)

#ifdef AMDGPU      
!$omp target teams distribute parallel do collapse(2) private(j,i)
#else
!$omp parallel do
#endif
        do j = 1, lzphi
        do i = 1, lr_shd
          xxx1(lrstart-lr_shd-1+i,j,1) = down_recv_r(i,j,1,1)
          xxx1(lrend+i,j,1) = up_recv_r(i,j,1,1)

          xxx2(lrstart-lr_shd-1+i,j,1) = down_recv_r(i,j,1,2)
          xxx2(lrend+i,j,1) = up_recv_r(i,j,1,2)

          xxx3(lrstart-lr_shd-1+i,j,1) = down_recv_r(i,j,1,3)
          xxx3(lrend+i,j,1) = up_recv_r(i,j,1,3)

          xxx4(lrstart-lr_shd-1+i,j,1) = down_recv_r(i,j,1,4)
          xxx4(lrend+i,j,1) = up_recv_r(i,j,1,4)

          xxx5(lrstart-lr_shd-1+i,j,1) = down_recv_r(i,j,1,5)
          xxx5(lrend+i,j,1) = up_recv_r(i,j,1,5)

          xxx6(lrstart-lr_shd-1+i,j,1) = down_recv_r(i,j,1,6)
          xxx6(lrend+i,j,1) = up_recv_r(i,j,1,6)

          xxx7(lrstart-lr_shd-1+i,j,1) = down_recv_r(i,j,1,7)
          xxx7(lrend+i,j,1) = up_recv_r(i,j,1,7)
        end do
        end do

      end if

      end if  !end of r-communication

end
!--------------------------------------------------------------------
subroutine periodic_field_mlt8(xxx1,xxx2,xxx3,xxx4,xxx5,xxx6,xxx7,xxx8)
! boundary condition with mpi communications for 8 variables
! 2015-06-28
!--------------------------------------------------------------------
      use parameters
      use mpiset
      implicit none
      integer, parameter::mlt=8
      real(8)::xxx1(lr,lz,lphi),xxx2(lr,lz,lphi),xxx3(lr,lz,lphi),xxx4(lr,lz,lphi)
      real(8)::xxx5(lr,lz,lphi),xxx6(lr,lz,lphi),xxx7(lr,lz,lphi),xxx8(lr,lz,lphi)
      integer::i,j,k
      integer::node_up,node_down,lrz2,lrphi2,lzphi2

      Real(8)::up_send_phi(lr,lz,lphi_shd,mlt),down_send_phi(lr,lz,lphi_shd,mlt)
      real(8)::up_recv_phi(lr,lz,lphi_shd,mlt),down_recv_phi(lr,lz,lphi_shd,mlt)
      real(8)::up_send_z(lr,lz_shd,lphi,mlt),down_send_z(lr,lz_shd,lphi,mlt)
      real(8)::up_recv_z(lr,lz_shd,lphi,mlt),down_recv_z(lr,lz_shd,lphi,mlt)
      real(8)::up_send_r(lr_shd,lz,lphi,mlt),down_send_r(lr_shd,lz,lphi,mlt)
      real(8)::up_recv_r(lr_shd,lz,lphi,mlt),down_recv_r(lr_shd,lz,lphi,mlt)


! periodic in phi-direction

      node_up = mod(my_rank_phi + 1, mpi_proc_phi)
      node_down = mod(my_rank_phi - 1 + mpi_proc_phi, mpi_proc_phi)

#ifdef AMDGPU      
!$omp target teams distribute parallel do collapse(2) private(k,i)
#else
!$omp parallel do
#endif
      do k = 1, lphi_shd
      do i = 1, lrz
        down_send_phi(i,1,k,1) = xxx1(i,1,lphistart+k-1)
        up_send_phi(i,1,k,1) = xxx1(i,1,lphiend-lphi_shd+k)

        down_send_phi(i,1,k,2) = xxx2(i,1,lphistart+k-1)
        up_send_phi(i,1,k,2) = xxx2(i,1,lphiend-lphi_shd+k)

        down_send_phi(i,1,k,3) = xxx3(i,1,lphistart+k-1)
        up_send_phi(i,1,k,3) = xxx3(i,1,lphiend-lphi_shd+k)

        down_send_phi(i,1,k,4) = xxx4(i,1,lphistart+k-1)
        up_send_phi(i,1,k,4) = xxx4(i,1,lphiend-lphi_shd+k)

        down_send_phi(i,1,k,5) = xxx5(i,1,lphistart+k-1)
        up_send_phi(i,1,k,5) = xxx5(i,1,lphiend-lphi_shd+k)

        down_send_phi(i,1,k,6) = xxx6(i,1,lphistart+k-1)
        up_send_phi(i,1,k,6) = xxx6(i,1,lphiend-lphi_shd+k)

        down_send_phi(i,1,k,7) = xxx7(i,1,lphistart+k-1)
        up_send_phi(i,1,k,7) = xxx7(i,1,lphiend-lphi_shd+k)

        down_send_phi(i,1,k,8) = xxx8(i,1,lphistart+k-1)
        up_send_phi(i,1,k,8) = xxx8(i,1,lphiend-lphi_shd+k)
      end do
      end do

! communication of boundary data

      lrz2 = lr*lz*lphi_shd*mlt

      call mpi_isend(up_send_phi(1,1,1,1),lrz2,mpi_real8,node_up,1 &
                    ,phi_world,my_request(1),mpi_err)
      call mpi_isend(down_send_phi(1,1,1,1),lrz2,mpi_real8,node_down,1 &
                    ,phi_world,my_request(2),mpi_err)
      call mpi_irecv(down_recv_phi(1,1,1,1),lrz2,mpi_real8,node_down,1 &
                    ,phi_world,my_request(3),mpi_err)
      call mpi_irecv(up_recv_phi(1,1,1,1),lrz2,mpi_real8,node_up,1 &
                    ,phi_world,my_request(4),mpi_err)

      call mpi_waitall(4,my_request,my_status_all,mpi_err)

#ifdef AMDGPU      
!$omp target teams distribute parallel do collapse(2) private(k,i)
#else
!$omp parallel do
#endif
      do k = 1, lphi_shd
      do i = 1, lrz
        xxx1(i,1,lphistart-lphi_shd-1+k) = down_recv_phi(i,1,k,1)
        xxx1(i,1,lphiend+k) = up_recv_phi(i,1,k,1)

        xxx2(i,1,lphistart-lphi_shd-1+k) = down_recv_phi(i,1,k,2)
        xxx2(i,1,lphiend+k) = up_recv_phi(i,1,k,2)

        xxx3(i,1,lphistart-lphi_shd-1+k) = down_recv_phi(i,1,k,3)
        xxx3(i,1,lphiend+k) = up_recv_phi(i,1,k,3)

        xxx4(i,1,lphistart-lphi_shd-1+k) = down_recv_phi(i,1,k,4)
        xxx4(i,1,lphiend+k) = up_recv_phi(i,1,k,4)

        xxx5(i,1,lphistart-lphi_shd-1+k) = down_recv_phi(i,1,k,5)
        xxx5(i,1,lphiend+k) = up_recv_phi(i,1,k,5)

        xxx6(i,1,lphistart-lphi_shd-1+k) = down_recv_phi(i,1,k,6)
        xxx6(i,1,lphiend+k) = up_recv_phi(i,1,k,6)

        xxx7(i,1,lphistart-lphi_shd-1+k) = down_recv_phi(i,1,k,7)
        xxx7(i,1,lphiend+k) = up_recv_phi(i,1,k,7)

        xxx8(i,1,lphistart-lphi_shd-1+k) = down_recv_phi(i,1,k,8)
        xxx8(i,1,lphiend+k) = up_recv_phi(i,1,k,8)
      end do
      end do

!------------------------------------------------------
! communication in z-direction
      if(mpi_proc_z.ge.2)then

      node_up = mod(my_rank_z + 1, mpi_proc_z)
      node_down = mod(my_rank_z - 1 + mpi_proc_z, mpi_proc_z)

      if(my_rank_z.ne.(mpi_proc_z-1))then
#ifdef AMDGPU      
!$omp target teams distribute parallel do collapse(3) private(k,j,i)
#else
!$omp parallel do
#endif
        do k = 1, lphi
        do j = 1, lz_shd
        do i = 1, lr
          up_send_z(i,j,k,1) = xxx1(i,lzend-lz_shd+j,k)
          up_send_z(i,j,k,2) = xxx2(i,lzend-lz_shd+j,k)
          up_send_z(i,j,k,3) = xxx3(i,lzend-lz_shd+j,k)
          up_send_z(i,j,k,4) = xxx4(i,lzend-lz_shd+j,k)
          up_send_z(i,j,k,5) = xxx5(i,lzend-lz_shd+j,k)
          up_send_z(i,j,k,6) = xxx6(i,lzend-lz_shd+j,k)
          up_send_z(i,j,k,7) = xxx7(i,lzend-lz_shd+j,k)
          up_send_z(i,j,k,8) = xxx8(i,lzend-lz_shd+j,k)
        end do
        end do
        end do
      end if

      if(my_rank_z.ne.0 )then
#ifdef AMDGPU      
!$omp target teams distribute parallel do collapse(3) private(k,j,i)
#else
!$omp parallel do
#endif
        do k = 1, lphi
        do j = 1, lz_shd
        do i = 1, lr
          down_send_z(i,j,k,1) = xxx1(i,lzstart+j-1,k)
          down_send_z(i,j,k,2) = xxx2(i,lzstart+j-1,k)
          down_send_z(i,j,k,3) = xxx3(i,lzstart+j-1,k)
          down_send_z(i,j,k,4) = xxx4(i,lzstart+j-1,k)
          down_send_z(i,j,k,5) = xxx5(i,lzstart+j-1,k)
          down_send_z(i,j,k,6) = xxx6(i,lzstart+j-1,k)
          down_send_z(i,j,k,7) = xxx7(i,lzstart+j-1,k)
          down_send_z(i,j,k,8) = xxx8(i,lzstart+j-1,k)
        end do
        end do
        end do
      end if

! communication of boundary data

      lrphi2 = lr*lphi*lz_shd*mlt

      if(my_rank_z.eq.0)then
        call mpi_isend(up_send_z(1,1,1,1),lrphi2,mpi_real8,node_up,1 &
                      ,z_world,my_request(1),mpi_err)
        call mpi_irecv(up_recv_z(1,1,1,1),lrphi2,mpi_real8,node_up,1 &
                      ,z_world,my_request(2),mpi_err)

        call mpi_waitall(2,my_request,my_status_all,mpi_err)

#ifdef AMDGPU      
!$omp target teams distribute parallel do collapse(3) private(k,j,i)
#else
!$omp parallel do
#endif
        do k = 1, lphi
        do j = 1, lz_shd
        do i = 1, lr
          xxx1(i,lzend+j,k) = up_recv_z(i,j,k,1)
          xxx2(i,lzend+j,k) = up_recv_z(i,j,k,2)
          xxx3(i,lzend+j,k) = up_recv_z(i,j,k,3)
          xxx4(i,lzend+j,k) = up_recv_z(i,j,k,4)
          xxx5(i,lzend+j,k) = up_recv_z(i,j,k,5)
          xxx6(i,lzend+j,k) = up_recv_z(i,j,k,6)
          xxx7(i,lzend+j,k) = up_recv_z(i,j,k,7)
          xxx8(i,lzend+j,k) = up_recv_z(i,j,k,8)
        end do
        end do
        end do

      else if(my_rank_z.eq.(mpi_proc_z-1) )then
        call mpi_isend(down_send_z(1,1,1,1),lrphi2,mpi_real8,node_down,1 &
                      ,z_world,my_request(1),mpi_err)
        call mpi_irecv(down_recv_z(1,1,1,1),lrphi2,mpi_real8,node_down,1 &
                      ,z_world,my_request(2),mpi_err)

        call mpi_waitall(2,my_request,my_status_all,mpi_err)

#ifdef AMDGPU      
!$omp target teams distribute parallel do collapse(3) private(k,j,i)
#else
!$omp parallel do
#endif
        do k = 1, lphi
        do j = 1, lz_shd
        do i = 1, lr
          xxx1(i,lzstart-lz_shd-1+j,k) = down_recv_z(i,j,k,1)
          xxx2(i,lzstart-lz_shd-1+j,k) = down_recv_z(i,j,k,2)
          xxx3(i,lzstart-lz_shd-1+j,k) = down_recv_z(i,j,k,3)
          xxx4(i,lzstart-lz_shd-1+j,k) = down_recv_z(i,j,k,4)
          xxx5(i,lzstart-lz_shd-1+j,k) = down_recv_z(i,j,k,5)
          xxx6(i,lzstart-lz_shd-1+j,k) = down_recv_z(i,j,k,6)
          xxx7(i,lzstart-lz_shd-1+j,k) = down_recv_z(i,j,k,7)
          xxx8(i,lzstart-lz_shd-1+j,k) = down_recv_z(i,j,k,8)
        end do
        end do
        end do

      else
        call mpi_isend(up_send_z(1,1,1,1),lrphi2,mpi_real8,node_up,1 &
                      ,z_world,my_request(1),mpi_err)
        call mpi_isend(down_send_z(1,1,1,1),lrphi2,mpi_real8,node_down,1 &
                      ,z_world,my_request(2),mpi_err)
        call mpi_irecv(up_recv_z(1,1,1,1),lrphi2,mpi_real8,node_up,1 &
                      ,z_world,my_request(3),mpi_err)
        call mpi_irecv(down_recv_z(1,1,1,1),lrphi2,mpi_real8,node_down,1 &
                      ,z_world,my_request(4),mpi_err)

        call mpi_waitall(4,my_request,my_status_all,mpi_err)

#ifdef AMDGPU      
!$omp target teams distribute parallel do collapse(3) private(k,j,i)
#else
!$omp parallel do
#endif
        do k = 1, lphi
        do j = 1, lz_shd
        do i = 1, lr
          xxx1(i,lzstart-lz_shd-1+j,k) = down_recv_z(i,j,k,1)
          xxx1(i,lzend+j,k) = up_recv_z(i,j,k,1)

          xxx2(i,lzstart-lz_shd-1+j,k) = down_recv_z(i,j,k,2)
          xxx2(i,lzend+j,k) = up_recv_z(i,j,k,2)

          xxx3(i,lzstart-lz_shd-1+j,k) = down_recv_z(i,j,k,3)
          xxx3(i,lzend+j,k) = up_recv_z(i,j,k,3)

          xxx4(i,lzstart-lz_shd-1+j,k) = down_recv_z(i,j,k,4)
          xxx4(i,lzend+j,k) = up_recv_z(i,j,k,4)

          xxx5(i,lzstart-lz_shd-1+j,k) = down_recv_z(i,j,k,5)
          xxx5(i,lzend+j,k) = up_recv_z(i,j,k,5)

          xxx6(i,lzstart-lz_shd-1+j,k) = down_recv_z(i,j,k,6)
          xxx6(i,lzend+j,k) = up_recv_z(i,j,k,6)

          xxx7(i,lzstart-lz_shd-1+j,k) = down_recv_z(i,j,k,7)
          xxx7(i,lzend+j,k) = up_recv_z(i,j,k,7)

          xxx8(i,lzstart-lz_shd-1+j,k) = down_recv_z(i,j,k,8)
          xxx8(i,lzend+j,k) = up_recv_z(i,j,k,8)
        end do
        end do
        end do

      end if


      end if    !end of z-communication
!------------------------------------------------------
! communication in r-direction
      if(mpi_proc_r.ge.2)then

      node_up = mod(my_rank + 1, mpi_proc_r)
      node_down = mod(my_rank - 1 + mpi_proc_r, mpi_proc_r)

      if(my_rank_r.ne.(mpi_proc_r-1))then
#ifdef AMDGPU      
!$omp target teams distribute parallel do collapse(2) private(j,i)
#else
!$omp parallel do
#endif
        do j = 1, lzphi
        do i = 1, lr_shd
          up_send_r(i,j,1,1) = xxx1(lrend-lr_shd+i,j,1)
          up_send_r(i,j,1,2) = xxx2(lrend-lr_shd+i,j,1)
          up_send_r(i,j,1,3) = xxx3(lrend-lr_shd+i,j,1)
          up_send_r(i,j,1,4) = xxx4(lrend-lr_shd+i,j,1)
          up_send_r(i,j,1,5) = xxx5(lrend-lr_shd+i,j,1)
          up_send_r(i,j,1,6) = xxx6(lrend-lr_shd+i,j,1)
          up_send_r(i,j,1,7) = xxx7(lrend-lr_shd+i,j,1)
          up_send_r(i,j,1,8) = xxx8(lrend-lr_shd+i,j,1)
        end do
        end do
      end if

      if(my_rank_r.ne.0 )then
#ifdef AMDGPU      
!$omp target teams distribute parallel do collapse(2) private(j,i)
#else
!$omp parallel do
#endif
        do j = 1, lzphi
        do i = 1, lr_shd
          down_send_r(i,j,1,1) = xxx1(lrstart+i-1,  j,1)
          down_send_r(i,j,1,2) = xxx2(lrstart+i-1,  j,1)
          down_send_r(i,j,1,3) = xxx3(lrstart+i-1,  j,1)
          down_send_r(i,j,1,4) = xxx4(lrstart+i-1,  j,1)
          down_send_r(i,j,1,5) = xxx5(lrstart+i-1,  j,1)
          down_send_r(i,j,1,6) = xxx6(lrstart+i-1,  j,1)
          down_send_r(i,j,1,7) = xxx7(lrstart+i-1,  j,1)
          down_send_r(i,j,1,8) = xxx8(lrstart+i-1,  j,1)
        end do
        end do
      end if


! communication of boundary data

      lzphi2 = lz*lphi*lr_shd*mlt

      if(my_rank_r.eq.0)then
        call mpi_isend(up_send_r(1,1,1,1),lzphi2,mpi_real8,node_up,1 &
                      ,r_world,my_request(1),mpi_err)
        call mpi_irecv(up_recv_r(1,1,1,1),lzphi2,mpi_real8,node_up,1 &
                      ,r_world,my_request(2),mpi_err)

        call mpi_waitall(2,my_request,my_status_all,mpi_err)

#ifdef AMDGPU      
!$omp target teams distribute parallel do collapse(2) private(j,i)
#else
!$omp parallel do
#endif
        do j = 1, lzphi
        do i = 1, lr_shd
          xxx1(lrend+i,j,1) = up_recv_r(i,j,1,1)
          xxx2(lrend+i,j,1) = up_recv_r(i,j,1,2)
          xxx3(lrend+i,j,1) = up_recv_r(i,j,1,3)
          xxx4(lrend+i,j,1) = up_recv_r(i,j,1,4)
          xxx5(lrend+i,j,1) = up_recv_r(i,j,1,5)
          xxx6(lrend+i,j,1) = up_recv_r(i,j,1,6)
          xxx7(lrend+i,j,1) = up_recv_r(i,j,1,7)
          xxx8(lrend+i,j,1) = up_recv_r(i,j,1,8)
        end do
        end do

      else if(my_rank_r.eq.(mpi_proc_r-1) )then
        call mpi_isend(down_send_r(1,1,1,1),lzphi2,mpi_real8,node_down,1 &
                      ,r_world,my_request(1),mpi_err)
        call mpi_irecv(down_recv_r(1,1,1,1),lzphi2,mpi_real8,node_down,1 &
                      ,r_world,my_request(2),mpi_err)

        call mpi_waitall(2,my_request,my_status_all,mpi_err)

#ifdef AMDGPU      
!$omp target teams distribute parallel do collapse(2) private(j,i)
#else
!$omp parallel do
#endif
        do j = 1, lzphi
        do i = 1, lr_shd
          xxx1(lrstart-lr_shd-1+i,j,1) = down_recv_r(i,j,1,1)
          xxx2(lrstart-lr_shd-1+i,j,1) = down_recv_r(i,j,1,2)
          xxx3(lrstart-lr_shd-1+i,j,1) = down_recv_r(i,j,1,3)
          xxx4(lrstart-lr_shd-1+i,j,1) = down_recv_r(i,j,1,4)
          xxx5(lrstart-lr_shd-1+i,j,1) = down_recv_r(i,j,1,5)
          xxx6(lrstart-lr_shd-1+i,j,1) = down_recv_r(i,j,1,6)
          xxx7(lrstart-lr_shd-1+i,j,1) = down_recv_r(i,j,1,7)
          xxx8(lrstart-lr_shd-1+i,j,1) = down_recv_r(i,j,1,8)
        end do
        end do

      else
        call mpi_isend(up_send_r(1,1,1,1),lzphi2,mpi_real8,node_up,1 &
                      ,r_world,my_request(1),mpi_err)
        call mpi_isend(down_send_r(1,1,1,1),lzphi2,mpi_real8,node_down,1 &
                      ,r_world,my_request(2),mpi_err)
        call mpi_irecv(up_recv_r(1,1,1,1),lzphi2,mpi_real8,node_up,1 &
                      ,r_world,my_request(3),mpi_err)
        call mpi_irecv(down_recv_r(1,1,1,1),lzphi2,mpi_real8,node_down,1 &
                      ,r_world,my_request(4),mpi_err)

        call mpi_waitall(4,my_request,my_status_all,mpi_err)

#ifdef AMDGPU      
!$omp target teams distribute parallel do collapse(2) private(j,i)
#else
!$omp parallel do
#endif
        do j = 1, lzphi
        do i = 1, lr_shd
          xxx1(lrstart-lr_shd-1+i,j,1) = down_recv_r(i,j,1,1)
          xxx1(lrend+i,j,1) = up_recv_r(i,j,1,1)

          xxx2(lrstart-lr_shd-1+i,j,1) = down_recv_r(i,j,1,2)
          xxx2(lrend+i,j,1) = up_recv_r(i,j,1,2)

          xxx3(lrstart-lr_shd-1+i,j,1) = down_recv_r(i,j,1,3)
          xxx3(lrend+i,j,1) = up_recv_r(i,j,1,3)

          xxx4(lrstart-lr_shd-1+i,j,1) = down_recv_r(i,j,1,4)
          xxx4(lrend+i,j,1) = up_recv_r(i,j,1,4)

          xxx5(lrstart-lr_shd-1+i,j,1) = down_recv_r(i,j,1,5)
          xxx5(lrend+i,j,1) = up_recv_r(i,j,1,5)

          xxx6(lrstart-lr_shd-1+i,j,1) = down_recv_r(i,j,1,6)
          xxx6(lrend+i,j,1) = up_recv_r(i,j,1,6)

          xxx7(lrstart-lr_shd-1+i,j,1) = down_recv_r(i,j,1,7)
          xxx7(lrend+i,j,1) = up_recv_r(i,j,1,7)

          xxx8(lrstart-lr_shd-1+i,j,1) = down_recv_r(i,j,1,8)
          xxx8(lrend+i,j,1) = up_recv_r(i,j,1,8)
        end do
        end do

      end if

      end if  !end of r-communication

end
!--------------------------------------------------------------------
subroutine periodic_particle_mlt2b(xxx1,xxx2)
! boundary condition for particle density with mpi communications
! summation over ptcl_world: 2012-06-21
! for 2 variables
! communicate for lx_shd * 2: 2016-12-23
! correction: 2017-02-12
!--------------------------------------------------------------------
      use parameters
      use mpiset
      implicit none
      integer, parameter::mlt=2
      integer, parameter::lrcom=lrnet/mpi_proc_r,lzcom=lznet/mpi_proc_z &
                         ,lphicom=lphinet/mpi_proc_phi
      integer, parameter::lr_shd2=lr_shd*2, lz_shd2=lz_shd*2, lphi_shd2=lphi_shd*2

      real(8)::xxx1(lr,lz,lphi),xxx2(lr,lz,lphi)
      real(8)::xxx_com(lr,lz,lphi,mlt)
      integer::i,j,k,lsize                      !2012-06-17
      integer::node_up,node_down,lrz2,lrphi2,lzphi2
      integer::lr_down,lr_up,lz_down,lz_up,lphi_down,lphi_up

      real(8)::up_send_phi(lr,lz,lphi_shd2,mlt),down_send_phi(lr,lz,lphi_shd2,mlt)
      real(8)::up_recv_phi(lr,lz,lphi_shd2,mlt),down_recv_phi(lr,lz,lphi_shd2,mlt)
      real(8)::up_send_z(lr,lz_shd2,lphi,mlt),down_send_z(lr,lz_shd2,lphi,mlt)
      real(8)::up_recv_z(lr,lz_shd2,lphi,mlt),down_recv_z(lr,lz_shd2,lphi,mlt)
      real(8)::up_send_r(lr_shd2,lz,lphi,mlt),down_send_r(lr_shd2,lz,lphi,mlt)
      real(8)::up_recv_r(lr_shd2,lz,lphi,mlt),down_recv_r(lr_shd2,lz,lphi,mlt)


! summation on phi-boundary

      node_up = mod(my_rank_phi + 1, mpi_proc_phi)
      node_down = mod(my_rank_phi - 1 + mpi_proc_phi, mpi_proc_phi)

      lr_down = lrstart - lr_shd - 1
      lr_up   = lrend   - lr_shd

      lz_down = lzstart - lz_shd - 1
      lz_up   = lzend   - lz_shd

      lphi_down = lphistart - lphi_shd - 1
      lphi_up   = lphiend   - lphi_shd

#ifdef AMDGPU
!$omp target teams distribute parallel do collapse(3) private(k,j,i)
#else
!$omp parallel do
#endif
      do k = 1, lphi_shd2
      do j = 1, lz
      do i = 1, lr
        down_send_phi(i,j,k,1)= xxx1(i,j,lphi_down + k)
        up_send_phi(i,j,k,1)  = xxx1(i,j,lphi_up   + k)

        down_send_phi(i,j,k,2)= xxx2(i,j,lphi_down + k)
        up_send_phi(i,j,k,2)  = xxx2(i,j,lphi_up   + k)
      end do
      end do
      end do

! communication of particle data

      lrz2 = lr*lz*lphi_shd2*mlt

      call mpi_isend(up_send_phi(1,1,1,1),lrz2,mpi_real8,node_up,1 &
                    ,phi_world,my_request(1),mpi_err)
      call mpi_isend(down_send_phi(1,1,1,1),lrz2,mpi_real8,node_down,1 &
                    ,phi_world,my_request(2),mpi_err)
      call mpi_irecv(down_recv_phi(1,1,1,1),lrz2,mpi_real8,node_down,1 &
                    ,phi_world,my_request(3),mpi_err)
      call mpi_irecv(up_recv_phi(1,1,1,1),lrz2,mpi_real8,node_up,1 &
                    ,phi_world,my_request(4),mpi_err)

      call mpi_waitall(4,my_request,my_status_all,mpi_err)

#ifdef AMDGPU
!$omp target teams distribute parallel do collapse(2) private(k,i)
#else
!$omp parallel do
#endif
      do k = 1, lphi_shd2
      do i = 1, lrz
        xxx1(i,1,lphi_down+ k)  = xxx1(i,1,lphi_down + k)  + down_recv_phi(i,1,k,1)

        xxx2(i,1,lphi_down+ k)  = xxx2(i,1,lphi_down + k)  + down_recv_phi(i,1,k,2)
      end do
      end do

#ifdef AMDGPU
!$omp target teams distribute parallel do collapse(2) private(k,i)
#else
!$omp parallel do
#endif
      do k = 1, lphi_shd2
      do i = 1, lrz
        xxx1(i,1,lphi_up  + k)  = xxx1(i,1,lphi_up   + k)  + up_recv_phi(i,1,k,1)

        xxx2(i,1,lphi_up  + k)  = xxx2(i,1,lphi_up   + k)  + up_recv_phi(i,1,k,2)
      end do
      end do


!------------------------------------------------------
! communication in z-direction
      if(mpi_proc_z.ge.2)then

      node_up = mod(my_rank_z + 1, mpi_proc_z)
      node_down = mod(my_rank_z - 1 + mpi_proc_z, mpi_proc_z)

      if(my_rank_z.ne.(mpi_proc_z-1))then
#ifdef AMDGPU
!$omp target teams distribute parallel do collapse(3) private(k,j,i)
#else
!$omp parallel do
#endif
        do k = 1, lphi
        do j = 1, lz_shd2
        do i = 1, lr
          up_send_z(i,j,k,1) = xxx1(i,lz_up + j,k)
          up_send_z(i,j,k,2) = xxx2(i,lz_up + j,k)
        end do
        end do
        end do
      end if

      if(my_rank_z.ne.0 )then
#ifdef AMDGPU
!$omp target teams distribute parallel do collapse(3) private(k,j,i)
#else
!$omp parallel do
#endif
        do k = 1, lphi
        do j = 1, lz_shd2
        do i = 1, lr
          down_send_z(i,j,k,1) = xxx1(i,lz_down + j,k)
          down_send_z(i,j,k,2) = xxx2(i,lz_down + j,k)
        end do
        end do
        end do
      end if

! communication of particle data

      lrphi2 = lr*lphi*lz_shd2*mlt

      if(my_rank_z.eq.0)then
        call mpi_isend(up_send_z(1,1,1,1),lrphi2,mpi_real8,node_up,1 &
                      ,z_world,my_request(1),mpi_err)
        call mpi_irecv(up_recv_z(1,1,1,1),lrphi2,mpi_real8,node_up,1 &
                      ,z_world,my_request(2),mpi_err)

        call mpi_waitall(2,my_request,my_status_all,mpi_err)

#ifdef AMDGPU
!$omp target teams distribute parallel do collapse(3) private(k,j,i)
#else
!$omp parallel do
#endif
        do k = 1, lphi
        do j = 1, lz_shd2
        do i = 1, lr
          xxx1(i,lz_up+j,k)= xxx1(i,lz_up+j,k) + up_recv_z(i,j,k,1)
          xxx2(i,lz_up+j,k)= xxx2(i,lz_up+j,k) + up_recv_z(i,j,k,2)
        end do
        end do
        end do

      else if(my_rank_z.eq.(mpi_proc_z-1) )then
        call mpi_isend(down_send_z(1,1,1,1),lrphi2,mpi_real8,node_down,1 &
                      ,z_world,my_request(1),mpi_err)
        call mpi_irecv(down_recv_z(1,1,1,1),lrphi2,mpi_real8,node_down,1 &
                      ,z_world,my_request(2),mpi_err)

        call mpi_waitall(2,my_request,my_status_all,mpi_err)

#ifdef AMDGPU
!$omp target teams distribute parallel do collapse(3) private(k,j,i)
#else
!$omp parallel do
#endif
        do k = 1, lphi
        do j = 1, lz_shd2
        do i = 1, lr
          xxx1(i,lz_down+j,k)  = xxx1(i,lz_down+j,k)  + down_recv_z(i,j,k,1)
          xxx2(i,lz_down+j,k)  = xxx2(i,lz_down+j,k)  + down_recv_z(i,j,k,2)
        end do
        end do
        end do

      else
        call mpi_isend(up_send_z(1,1,1,1),lrphi2,mpi_real8,node_up,1 &
                      ,z_world,my_request(1),mpi_err)
        call mpi_isend(down_send_z(1,1,1,1),lrphi2,mpi_real8,node_down,1 &
                      ,z_world,my_request(2),mpi_err)
        call mpi_irecv(up_recv_z(1,1,1,1),lrphi2,mpi_real8,node_up,1 &
                      ,z_world,my_request(3),mpi_err)
        call mpi_irecv(down_recv_z(1,1,1,1),lrphi2,mpi_real8,node_down,1 &
                      ,z_world,my_request(4),mpi_err)

        call mpi_waitall(4,my_request,my_status_all,mpi_err)


#ifdef AMDGPU
!$omp target teams distribute parallel do collapse(3) private(k,j,i)
#else
!$omp parallel do
#endif
        do k = 1, lphi
        do j = 1, lz_shd2
        do i = 1, lr
          xxx1(i,lz_down+j,k)  = xxx1(i,lz_down+j,k)  + down_recv_z(i,j,k,1)
          xxx2(i,lz_down+j,k)  = xxx2(i,lz_down+j,k)  + down_recv_z(i,j,k,2)
        end do
        end do
        end do

#ifdef AMDGPU
!$omp target teams distribute parallel do collapse(3) private(k,j,i)
#else
!$omp parallel do
#endif
        do k = 1, lphi
        do j = 1, lz_shd2
        do i = 1, lr
          xxx1(i,lz_up+j,k)= xxx1(i,lz_up+j,k) + up_recv_z(i,j,k,1)
          xxx2(i,lz_up+j,k)= xxx2(i,lz_up+j,k) + up_recv_z(i,j,k,2)
        end do
        end do
        end do

      end if

      end if  !end of z-communication

!------------------------------------------------------
! communication in r-direction
      if(mpi_proc_r.ge.2)then

      node_up = mod(my_rank_r + 1, mpi_proc_r)
      node_down = mod(my_rank_r - 1 + mpi_proc_r, mpi_proc_r)

      if(my_rank_r.ne.(mpi_proc_r-1))then
#ifdef AMDGPU
!$omp target teams distribute parallel do collapse(2) private(j,i)
#else
!$omp parallel do
#endif         
        do j = 1, lzphi
        do i = 1, lr_shd2
          up_send_r(i,j,1,1) = xxx1(lr_up+i,j,1)
          up_send_r(i,j,1,2) = xxx2(lr_up+i,j,1)
        end do
        end do
      end if

      if(my_rank_r.ne.0 )then
#ifdef AMDGPU
!$omp target teams distribute parallel do collapse(2) private(j,i)
#else
!$omp parallel do
#endif
        do j = 1, lzphi
        do i = 1, lr_shd2
          down_send_r(i,j,1,1) = xxx1(lr_down+i,j,1)
          down_send_r(i,j,1,2) = xxx2(lr_down+i,j,1)
        end do
        end do
      end if


! communication of particle data

      lzphi2 = lz*lphi*lr_shd2*mlt

      if(my_rank_r.eq.0)then
        call mpi_isend(up_send_r(1,1,1,1),lzphi2,mpi_real8,node_up,1 &
                      ,r_world,my_request(1),mpi_err)
        call mpi_irecv(up_recv_r(1,1,1,1),lzphi2,mpi_real8,node_up,1 &
                      ,r_world,my_request(2),mpi_err)

        call mpi_waitall(2,my_request,my_status_all,mpi_err)

#ifdef AMDGPU
!$omp target teams distribute parallel do collapse(2) private(j,i)
#else
!$omp parallel do
#endif
        do j = 1, lzphi
        do i = 1, lr_shd2
          xxx1(lr_up+i,j,1)= xxx1(lr_up+i,j,1) + up_recv_r(i,j,1,1)
          xxx2(lr_up+i,j,1)= xxx2(lr_up+i,j,1) + up_recv_r(i,j,1,2)
        end do
        end do

      else if(my_rank_r.eq.(mpi_proc_r-1) )then
        call mpi_isend(down_send_r(1,1,1,1),lzphi2,mpi_real8,node_down,1 &
                      ,r_world,my_request(1),mpi_err)
        call mpi_irecv(down_recv_r(1,1,1,1),lzphi2,mpi_real8,node_down,1 &
                      ,r_world,my_request(2),mpi_err)

        call mpi_waitall(2,my_request,my_status_all,mpi_err)

#ifdef AMDGPU
!$omp target teams distribute parallel do collapse(2) private(j,i)
#else
!$omp parallel do
#endif
        do j = 1, lzphi
        do i = 1, lr_shd2
          xxx1(lr_down+i,j,1)  = xxx1(lr_down+i,j,1)  + down_recv_r(i,j,1,1)
          xxx2(lr_down+i,j,1)  = xxx2(lr_down+i,j,1)  + down_recv_r(i,j,1,2)
        end do
        end do

      else
        call mpi_isend(up_send_r(1,1,1,1),lzphi2,mpi_real8,node_up,1 &
                      ,r_world,my_request(1),mpi_err)
        call mpi_isend(down_send_r(1,1,1,1),lzphi2,mpi_real8,node_down,1 &
                      ,r_world,my_request(2),mpi_err)
        call mpi_irecv(up_recv_r(1,1,1,1),lzphi2,mpi_real8,node_up,1 &
                      ,r_world,my_request(3),mpi_err)
        call mpi_irecv(down_recv_r(1,1,1,1),lzphi2,mpi_real8,node_down,1 &
                      ,r_world,my_request(4),mpi_err)

        call mpi_waitall(4,my_request,my_status_all,mpi_err)

#ifdef AMDGPU
!$omp target teams distribute parallel do collapse(2) private(j,i)
#else
!$omp parallel do
#endif
        do j = 1, lzphi
        do i = 1, lr_shd2
          xxx1(lr_down+i,j,1)  = xxx1(lr_down+i,j,1)  + down_recv_r(i,j,1,1)
          xxx2(lr_down+i,j,1)  = xxx2(lr_down+i,j,1)  + down_recv_r(i,j,1,2)
        end do
        end do

#ifdef AMDGPU
!$omp target teams distribute parallel do collapse(2) private(j,i)
#else
!$omp parallel do
#endif
        do j = 1, lzphi
        do i = 1, lr_shd2
          xxx1(lr_up+i,j,1)= xxx1(lr_up+i,j,1) + up_recv_r(i,j,1,1)
          xxx2(lr_up+i,j,1)= xxx2(lr_up+i,j,1) + up_recv_r(i,j,1,2)
        end do
        end do

      end if

      end if  !end of r-communication

!--------------------------------

! summation over ptcl_world
!2012-06-17
      if(mpi_proc_ptcl.ge.2)then
#ifdef AMDGPU
!$omp target teams distribute parallel do collapse(3) private(k,j,i)
#else
!$omp parallel do
#endif
        do k = 1, lphi
        do j = 1, lz
        do i = 1, lr
          xxx_com(i,j,k,1) = xxx1(i,j,k)
          xxx_com(i,j,k,2) = xxx2(i,j,k)
        end do
        end do
        end do

        lsize = lrzphi*mlt
        call mpi_allreduce(mpi_in_place,xxx_com(1,1,1,1),lsize,mpi_real8 &
                          ,mpi_sum,ptcl_world,mpi_err)

#ifdef AMDGPU
!$omp target teams distribute parallel do collapse(3) private(k,j,i)
#else
!$omp parallel do
#endif
        do k = 1, lphi
        do j = 1, lz
        do i = 1, lr
          xxx1(i,j,k) = xxx_com(i,j,k,1)
          xxx2(i,j,k) = xxx_com(i,j,k,2)
        end do
        end do
        end do
      end if
!2012-06-17 end


end
!--------------------------------------------------------------------
subroutine periodic_particle_mlt4b(xxx1,xxx2,xxx3,xxx4)
! boundary condition for particle density with mpi communications
! summation over ptcl_world: 2012-06-21
! for 4 variables: 2017-02-02
! communicate for lx_shd * 2: 2016-12-23
! correction: 2017-02-12
!--------------------------------------------------------------------
      use parameters
      use mpiset
      implicit none
      integer, parameter::mlt=4
      integer, parameter::lrcom=lrnet/mpi_proc_r,lzcom=lznet/mpi_proc_z &
                         ,lphicom=lphinet/mpi_proc_phi
      integer, parameter::lr_shd2=lr_shd*2, lz_shd2=lz_shd*2, lphi_shd2=lphi_shd*2

      real(8)::xxx1(lr,lz,lphi),xxx2(lr,lz,lphi),xxx3(lr,lz,lphi)
      real(8)::xxx4(lr,lz,lphi)
      real(8)::xxx_com(lr,lz,lphi,mlt)
      integer::i,j,k,lsize                      !2012-06-17
      integer::node_up,node_down,lrz2,lrphi2,lzphi2
      integer::lr_down,lr_up,lz_down,lz_up,lphi_down,lphi_up

      real(8)::up_send_phi(lr,lz,lphi_shd2,mlt),down_send_phi(lr,lz,lphi_shd2,mlt)
      real(8)::up_recv_phi(lr,lz,lphi_shd2,mlt),down_recv_phi(lr,lz,lphi_shd2,mlt)
      real(8)::up_send_z(lr,lz_shd2,lphi,mlt),down_send_z(lr,lz_shd2,lphi,mlt)
      real(8)::up_recv_z(lr,lz_shd2,lphi,mlt),down_recv_z(lr,lz_shd2,lphi,mlt)
      real(8)::up_send_r(lr_shd2,lz,lphi,mlt),down_send_r(lr_shd2,lz,lphi,mlt)
      real(8)::up_recv_r(lr_shd2,lz,lphi,mlt),down_recv_r(lr_shd2,lz,lphi,mlt)


! summation on phi-boundary

      node_up = mod(my_rank_phi + 1, mpi_proc_phi)
      node_down = mod(my_rank_phi - 1 + mpi_proc_phi, mpi_proc_phi)

      lr_down = lrstart - lr_shd - 1
      lr_up   = lrend   - lr_shd

      lz_down = lzstart - lz_shd - 1
      lz_up   = lzend   - lz_shd

      lphi_down = lphistart - lphi_shd - 1
      lphi_up   = lphiend   - lphi_shd

#ifdef AMDGPU
!$omp target teams distribute parallel do collapse(3) private(k,j,i)
#else
!$omp parallel do
#endif
      do k = 1, lphi_shd2
      do j = 1, lz
      do i = 1, lr
        down_send_phi(i,j,k,1)= xxx1(i,j,lphi_down + k)
        up_send_phi(i,j,k,1)  = xxx1(i,j,lphi_up   + k)

        down_send_phi(i,j,k,2)= xxx2(i,j,lphi_down + k)
        up_send_phi(i,j,k,2)  = xxx2(i,j,lphi_up   + k)

        down_send_phi(i,j,k,3)= xxx3(i,j,lphi_down + k)
        up_send_phi(i,j,k,3)  = xxx3(i,j,lphi_up   + k)

        down_send_phi(i,j,k,4)= xxx4(i,j,lphi_down + k)
        up_send_phi(i,j,k,4)  = xxx4(i,j,lphi_up   + k)
      end do
      end do
      end do

! communication of particle data

      lrz2 = lr*lz*lphi_shd2*mlt

      call mpi_isend(up_send_phi(1,1,1,1),lrz2,mpi_real8,node_up,1 &
                    ,phi_world,my_request(1),mpi_err)
      call mpi_isend(down_send_phi(1,1,1,1),lrz2,mpi_real8,node_down,1 &
                    ,phi_world,my_request(2),mpi_err)
      call mpi_irecv(down_recv_phi(1,1,1,1),lrz2,mpi_real8,node_down,1 &
                    ,phi_world,my_request(3),mpi_err)
      call mpi_irecv(up_recv_phi(1,1,1,1),lrz2,mpi_real8,node_up,1 &
                    ,phi_world,my_request(4),mpi_err)

      call mpi_waitall(4,my_request,my_status_all,mpi_err)

#ifdef AMDGPU
!$omp target teams distribute parallel do collapse(2) private(k,i)
#else
!$omp parallel do
#endif
      do k = 1, lphi_shd2
      do i = 1, lrz
        xxx1(i,1,lphi_down+ k) = xxx1(i,1,lphi_down + k) + down_recv_phi(i,1,k,1)

        xxx2(i,1,lphi_down+ k) = xxx2(i,1,lphi_down + k) + down_recv_phi(i,1,k,2)

        xxx3(i,1,lphi_down+ k) = xxx3(i,1,lphi_down + k) + down_recv_phi(i,1,k,3)

        xxx4(i,1,lphi_down+ k) = xxx4(i,1,lphi_down + k) + down_recv_phi(i,1,k,4)
      end do
      end do

#ifdef AMDGPU
!$omp target teams distribute parallel do collapse(2) private(k,i)
#else
!$omp parallel do
#endif
      do k = 1, lphi_shd2
      do i = 1, lrz
        xxx1(i,1,lphi_up  + k) = xxx1(i,1,lphi_up   + k) + up_recv_phi(i,1,k,1)

        xxx2(i,1,lphi_up  + k) = xxx2(i,1,lphi_up   + k) + up_recv_phi(i,1,k,2)

        xxx3(i,1,lphi_up  + k) = xxx3(i,1,lphi_up   + k) + up_recv_phi(i,1,k,3)

        xxx4(i,1,lphi_up  + k) = xxx4(i,1,lphi_up   + k) + up_recv_phi(i,1,k,4)
      end do
      end do


!------------------------------------------------------
! communication in z-direction
      if(mpi_proc_z.ge.2)then

      node_up = mod(my_rank_z + 1, mpi_proc_z)
      node_down = mod(my_rank_z - 1 + mpi_proc_z, mpi_proc_z)

      if(my_rank_z.ne.(mpi_proc_z-1))then
#ifdef AMDGPU
!$omp target teams distribute parallel do collapse(3) private(k,j,i)
#else
!$omp parallel do
#endif
        do k = 1, lphi
        do j = 1, lz_shd2
        do i = 1, lr
          up_send_z(i,j,k,1) = xxx1(i,lz_up + j,k)
          up_send_z(i,j,k,2) = xxx2(i,lz_up + j,k)
          up_send_z(i,j,k,3) = xxx3(i,lz_up + j,k)
          up_send_z(i,j,k,4) = xxx4(i,lz_up + j,k)
        end do
        end do
        end do
      end if

      if(my_rank_z.ne.0 )then
#ifdef AMDGPU
!$omp target teams distribute parallel do collapse(3) private(k,j,i)
#else
!$omp parallel do
#endif
        do k = 1, lphi
        do j = 1, lz_shd2
        do i = 1, lr
          down_send_z(i,j,k,1) = xxx1(i,lz_down + j,k)
          down_send_z(i,j,k,2) = xxx2(i,lz_down + j,k)
          down_send_z(i,j,k,3) = xxx3(i,lz_down + j,k)
          down_send_z(i,j,k,4) = xxx4(i,lz_down + j,k)
        end do
        end do
        end do
      end if

! communication of particle data

      lrphi2 = lr*lphi*lz_shd2*mlt

      if(my_rank_z.eq.0)then
        call mpi_isend(up_send_z(1,1,1,1),lrphi2,mpi_real8,node_up,1 &
                      ,z_world,my_request(1),mpi_err)
        call mpi_irecv(up_recv_z(1,1,1,1),lrphi2,mpi_real8,node_up,1 &
                      ,z_world,my_request(2),mpi_err)

        call mpi_waitall(2,my_request,my_status_all,mpi_err)

#ifdef AMDGPU
!$omp target teams distribute parallel do collapse(3) private(k,j,i)
#else
!$omp parallel do
#endif
        do k = 1, lphi
        do j = 1, lz_shd2
        do i = 1, lr
          xxx1(i,lz_up+j,k) = xxx1(i,lz_up+j,k) + up_recv_z(i,j,k,1)
          xxx2(i,lz_up+j,k) = xxx2(i,lz_up+j,k) + up_recv_z(i,j,k,2)
          xxx3(i,lz_up+j,k) = xxx3(i,lz_up+j,k) + up_recv_z(i,j,k,3)
          xxx4(i,lz_up+j,k) = xxx4(i,lz_up+j,k) + up_recv_z(i,j,k,4)
        end do
        end do
        end do

      else if(my_rank_z.eq.(mpi_proc_z-1) )then
        call mpi_isend(down_send_z(1,1,1,1),lrphi2,mpi_real8,node_down,1 &
                      ,z_world,my_request(1),mpi_err)
        call mpi_irecv(down_recv_z(1,1,1,1),lrphi2,mpi_real8,node_down,1 &
                      ,z_world,my_request(2),mpi_err)

        call mpi_waitall(2,my_request,my_status_all,mpi_err)

#ifdef AMDGPU
!$omp target teams distribute parallel do collapse(3) private(k,j,i)
#else
!$omp parallel do
#endif
        do k = 1, lphi
        do j = 1, lz_shd2
        do i = 1, lr
          xxx1(i,lz_down+j,k) = xxx1(i,lz_down+j,k) + down_recv_z(i,j,k,1)
          xxx2(i,lz_down+j,k) = xxx2(i,lz_down+j,k) + down_recv_z(i,j,k,2)
          xxx3(i,lz_down+j,k) = xxx3(i,lz_down+j,k) + down_recv_z(i,j,k,3)
          xxx4(i,lz_down+j,k) = xxx4(i,lz_down+j,k) + down_recv_z(i,j,k,4)
        end do
        end do
        end do

      else
        call mpi_isend(up_send_z(1,1,1,1),lrphi2,mpi_real8,node_up,1 &
                      ,z_world,my_request(1),mpi_err)
        call mpi_isend(down_send_z(1,1,1,1),lrphi2,mpi_real8,node_down,1 &
                      ,z_world,my_request(2),mpi_err)
        call mpi_irecv(up_recv_z(1,1,1,1),lrphi2,mpi_real8,node_up,1 &
                      ,z_world,my_request(3),mpi_err)
        call mpi_irecv(down_recv_z(1,1,1,1),lrphi2,mpi_real8,node_down,1 &
                      ,z_world,my_request(4),mpi_err)

        call mpi_waitall(4,my_request,my_status_all,mpi_err)

#ifdef AMDGPU
!$omp target teams distribute parallel do collapse(3) private(k,j,i)
#else
!$omp parallel do
#endif
        do k = 1, lphi
        do j = 1, lz_shd2
        do i = 1, lr
          xxx1(i,lz_down+j,k)= xxx1(i,lz_down+j,k) + down_recv_z(i,j,k,1)

          xxx2(i,lz_down+j,k)= xxx2(i,lz_down+j,k) + down_recv_z(i,j,k,2)

          xxx3(i,lz_down+j,k)= xxx3(i,lz_down+j,k) + down_recv_z(i,j,k,3)

          xxx4(i,lz_down+j,k)= xxx4(i,lz_down+j,k) + down_recv_z(i,j,k,4)
        end do
        end do
        end do

#ifdef AMDGPU
!$omp target teams distribute parallel do collapse(3) private(k,j,i)
#else
!$omp parallel do
#endif
        do k = 1, lphi
        do j = 1, lz_shd2
        do i = 1, lr
          xxx1(i,lz_up+j,k)  = xxx1(i,lz_up+j,k)   + up_recv_z(i,j,k,1)

          xxx2(i,lz_up+j,k)  = xxx2(i,lz_up+j,k)   + up_recv_z(i,j,k,2)

          xxx3(i,lz_up+j,k)  = xxx3(i,lz_up+j,k)   + up_recv_z(i,j,k,3)

          xxx4(i,lz_up+j,k)  = xxx4(i,lz_up+j,k)   + up_recv_z(i,j,k,4)
        end do
        end do
        end do

      end if

      end if  !end of z-communication

!------------------------------------------------------
! communication in r-direction
      if(mpi_proc_r.ge.2)then

      node_up = mod(my_rank_r + 1, mpi_proc_r)
      node_down = mod(my_rank_r - 1 + mpi_proc_r, mpi_proc_r)

      if(my_rank_r.ne.(mpi_proc_r-1))then
#ifdef AMDGPU
!$omp target teams distribute parallel do collapse(2) private(j,i)
#else
!$omp parallel do
#endif         
        do j = 1, lzphi
        do i = 1, lr_shd2
          up_send_r(i,j,1,1) = xxx1(lr_up+i,j,1)
          up_send_r(i,j,1,2) = xxx2(lr_up+i,j,1)
          up_send_r(i,j,1,3) = xxx3(lr_up+i,j,1)
          up_send_r(i,j,1,4) = xxx4(lr_up+i,j,1)
        end do
        end do
      end if

      if(my_rank_r.ne.0 )then
#ifdef AMDGPU
!$omp target teams distribute parallel do collapse(2) private(j,i)
#else
!$omp parallel do
#endif
        do j = 1, lzphi
        do i = 1, lr_shd2
          down_send_r(i,j,1,1) = xxx1(lr_down+i,j,1)
          down_send_r(i,j,1,2) = xxx2(lr_down+i,j,1)
          down_send_r(i,j,1,3) = xxx3(lr_down+i,j,1)
          down_send_r(i,j,1,4) = xxx4(lr_down+i,j,1)
        end do
        end do
      end if


! communication of particle data

      lzphi2 = lz*lphi*lr_shd2*mlt

      if(my_rank_r.eq.0)then
        call mpi_isend(up_send_r(1,1,1,1),lzphi2,mpi_real8,node_up,1 &
                      ,r_world,my_request(1),mpi_err)
        call mpi_irecv(up_recv_r(1,1,1,1),lzphi2,mpi_real8,node_up,1 &
                      ,r_world,my_request(2),mpi_err)

        call mpi_waitall(2,my_request,my_status_all,mpi_err)

#ifdef AMDGPU
!$omp target teams distribute parallel do collapse(2) private(j,i)
#else
!$omp parallel do
#endif
        do j = 1, lzphi
        do i = 1, lr_shd2
          xxx1(lr_up+i,j,1) = xxx1(lr_up+i,j,1) + up_recv_r(i,j,1,1)
          xxx2(lr_up+i,j,1) = xxx2(lr_up+i,j,1) + up_recv_r(i,j,1,2)
          xxx3(lr_up+i,j,1) = xxx3(lr_up+i,j,1) + up_recv_r(i,j,1,3)
          xxx4(lr_up+i,j,1) = xxx4(lr_up+i,j,1) + up_recv_r(i,j,1,4)
        end do
        end do

      else if(my_rank_r.eq.(mpi_proc_r-1) )then
        call mpi_isend(down_send_r(1,1,1,1),lzphi2,mpi_real8,node_down,1 &
                      ,r_world,my_request(1),mpi_err)
        call mpi_irecv(down_recv_r(1,1,1,1),lzphi2,mpi_real8,node_down,1 &
                      ,r_world,my_request(2),mpi_err)

        call mpi_waitall(2,my_request,my_status_all,mpi_err)

#ifdef AMDGPU
!$omp target teams distribute parallel do collapse(2) private(j,i)
#else
!$omp parallel do
#endif
        do j = 1, lzphi
        do i = 1, lr_shd2
          xxx1(lr_down+i,j,1) = xxx1(lr_down+i,j,1) + down_recv_r(i,j,1,1)
          xxx2(lr_down+i,j,1) = xxx2(lr_down+i,j,1) + down_recv_r(i,j,1,2)
          xxx3(lr_down+i,j,1) = xxx3(lr_down+i,j,1) + down_recv_r(i,j,1,3)
          xxx4(lr_down+i,j,1) = xxx4(lr_down+i,j,1) + down_recv_r(i,j,1,4)
        end do
        end do

      else
        call mpi_isend(up_send_r(1,1,1,1),lzphi2,mpi_real8,node_up,1 &
                      ,r_world,my_request(1),mpi_err)
        call mpi_isend(down_send_r(1,1,1,1),lzphi2,mpi_real8,node_down,1 &
                      ,r_world,my_request(2),mpi_err)
        call mpi_irecv(up_recv_r(1,1,1,1),lzphi2,mpi_real8,node_up,1 &
                      ,r_world,my_request(3),mpi_err)
        call mpi_irecv(down_recv_r(1,1,1,1),lzphi2,mpi_real8,node_down,1 &
                      ,r_world,my_request(4),mpi_err)

        call mpi_waitall(4,my_request,my_status_all,mpi_err)

#ifdef AMDGPU
!$omp target teams distribute parallel do collapse(2) private(j,i)
#else
!$omp parallel do
#endif
        do j = 1, lzphi
        do i = 1, lr_shd2
          xxx1(lr_down+i,j,1)= xxx1(lr_down+i,j,1) + down_recv_r(i,j,1,1)

          xxx2(lr_down+i,j,1)= xxx2(lr_down+i,j,1) + down_recv_r(i,j,1,2)

          xxx3(lr_down+i,j,1)= xxx3(lr_down+i,j,1) + down_recv_r(i,j,1,3)

          xxx4(lr_down+i,j,1)= xxx4(lr_down+i,j,1) + down_recv_r(i,j,1,4)
        end do
        end do

#ifdef AMDGPU
!$omp target teams distribute parallel do collapse(2) private(j,i)
#else
!$omp parallel do
#endif
        do j = 1, lzphi
        do i = 1, lr_shd2
          xxx1(lr_up+i,j,1)  = xxx1(lr_up+i,j,1)   + up_recv_r(i,j,1,1)

          xxx2(lr_up+i,j,1)  = xxx2(lr_up+i,j,1)   + up_recv_r(i,j,1,2)

          xxx3(lr_up+i,j,1)  = xxx3(lr_up+i,j,1)   + up_recv_r(i,j,1,3)

          xxx4(lr_up+i,j,1)  = xxx4(lr_up+i,j,1)   + up_recv_r(i,j,1,4)
        end do
        end do

      end if

      end if  !end of r-communication

!--------------------------------

! summation over ptcl_world
!2012-06-17
      if(mpi_proc_ptcl.ge.2)then
#ifdef AMDGPU
!$omp target teams distribute parallel do collapse(3) private(k,j,i)
#else
!$omp parallel do
#endif
        do k = 1, lphi
        do j = 1, lz
        do i = 1, lr
          xxx_com(i,j,k,1) = xxx1(i,j,k)
          xxx_com(i,j,k,2) = xxx2(i,j,k)
          xxx_com(i,j,k,3) = xxx3(i,j,k)
          xxx_com(i,j,k,4) = xxx4(i,j,k)
        end do
        end do
        end do

        lsize = lrzphi*mlt
        call mpi_allreduce(mpi_in_place,xxx_com(1,1,1,1),lsize,mpi_real8 &
                          ,mpi_sum,ptcl_world,mpi_err)

#ifdef AMDGPU
!$omp target teams distribute parallel do collapse(3) private(k,j,i)
#else
!$omp parallel do
#endif
        do k = 1, lphi
        do j = 1, lz
        do i = 1, lr
          xxx1(i,j,k) = xxx_com(i,j,k,1)
          xxx2(i,j,k) = xxx_com(i,j,k,2)
          xxx3(i,j,k) = xxx_com(i,j,k,3)
          xxx4(i,j,k) = xxx_com(i,j,k,4)
        end do
        end do
        end do
      end if
!2012-06-17 end


end
!--------------------------------------------------------------------
subroutine periodic_particle_mlt6b(xxx1,xxx2,xxx3,xxx4,xxx5,xxx6)
! boundary condition for particle density with mpi communications
! summation over ptcl_world: 2012-06-21
! for 6 variables: 2017-02-02
! communicate for lx_shd * 2: 2016-12-23
! correction: 2017-02-12
!--------------------------------------------------------------------
      use parameters
      use mpiset
      implicit none
      integer, parameter::mlt=6
      integer, parameter::lrcom=lrnet/mpi_proc_r,lzcom=lznet/mpi_proc_z &
                         ,lphicom=lphinet/mpi_proc_phi
      integer, parameter::lr_shd2=lr_shd*2, lz_shd2=lz_shd*2, lphi_shd2=lphi_shd*2

      real(8)::xxx1(lr,lz,lphi),xxx2(lr,lz,lphi),xxx3(lr,lz,lphi)
      real(8)::xxx4(lr,lz,lphi),xxx5(lr,lz,lphi),xxx6(lr,lz,lphi)
      real(8)::xxx_com(lr,lz,lphi,mlt)
      integer::i,j,k,lsize                      !2012-06-17
      integer::node_up,node_down,lrz2,lrphi2,lzphi2
      integer::lr_down,lr_up,lz_down,lz_up,lphi_down,lphi_up

      real(8)::up_send_phi(lr,lz,lphi_shd2,mlt),down_send_phi(lr,lz,lphi_shd2,mlt)
      real(8)::up_recv_phi(lr,lz,lphi_shd2,mlt),down_recv_phi(lr,lz,lphi_shd2,mlt)
      real(8)::up_send_z(lr,lz_shd2,lphi,mlt),down_send_z(lr,lz_shd2,lphi,mlt)
      real(8)::up_recv_z(lr,lz_shd2,lphi,mlt),down_recv_z(lr,lz_shd2,lphi,mlt)
      real(8)::up_send_r(lr_shd2,lz,lphi,mlt),down_send_r(lr_shd2,lz,lphi,mlt)
      real(8)::up_recv_r(lr_shd2,lz,lphi,mlt),down_recv_r(lr_shd2,lz,lphi,mlt)


! summation on phi-boundary

      node_up = mod(my_rank_phi + 1, mpi_proc_phi)
      node_down = mod(my_rank_phi - 1 + mpi_proc_phi, mpi_proc_phi)

      lr_down = lrstart - lr_shd - 1
      lr_up   = lrend   - lr_shd

      lz_down = lzstart - lz_shd - 1
      lz_up   = lzend   - lz_shd

      lphi_down = lphistart - lphi_shd - 1
      lphi_up   = lphiend   - lphi_shd

#ifdef AMDGPU
!$omp target teams distribute parallel do collapse(3) private(k,j,i)
#else
!$omp parallel do
#endif
      do k = 1, lphi_shd2
      do j = 1, lz
      do i = 1, lr
        down_send_phi(i,j,k,1)= xxx1(i,j,lphi_down + k)
        up_send_phi(i,j,k,1)  = xxx1(i,j,lphi_up   + k)

        down_send_phi(i,j,k,2)= xxx2(i,j,lphi_down + k)
        up_send_phi(i,j,k,2)  = xxx2(i,j,lphi_up   + k)

        down_send_phi(i,j,k,3)= xxx3(i,j,lphi_down + k)
        up_send_phi(i,j,k,3)  = xxx3(i,j,lphi_up   + k)

        down_send_phi(i,j,k,4)= xxx4(i,j,lphi_down + k)
        up_send_phi(i,j,k,4)  = xxx4(i,j,lphi_up   + k)

        down_send_phi(i,j,k,5)= xxx5(i,j,lphi_down + k)
        up_send_phi(i,j,k,5)  = xxx5(i,j,lphi_up   + k)

        down_send_phi(i,j,k,6)= xxx6(i,j,lphi_down + k)
        up_send_phi(i,j,k,6)  = xxx6(i,j,lphi_up   + k)
      end do
      end do
      end do

! communication of particle data

      lrz2 = lr*lz*lphi_shd2*mlt

      call mpi_isend(up_send_phi(1,1,1,1),lrz2,mpi_real8,node_up,1 &
                    ,phi_world,my_request(1),mpi_err)
      call mpi_isend(down_send_phi(1,1,1,1),lrz2,mpi_real8,node_down,1 &
                    ,phi_world,my_request(2),mpi_err)
      call mpi_irecv(down_recv_phi(1,1,1,1),lrz2,mpi_real8,node_down,1 &
                    ,phi_world,my_request(3),mpi_err)
      call mpi_irecv(up_recv_phi(1,1,1,1),lrz2,mpi_real8,node_up,1 &
                    ,phi_world,my_request(4),mpi_err)

      call mpi_waitall(4,my_request,my_status_all,mpi_err)

#ifdef AMDGPU
!$omp target teams distribute parallel do collapse(2) private(k,i)
#else
!$omp parallel do
#endif
      do k = 1, lphi_shd2
      do i = 1, lrz
        xxx1(i,1,lphi_down+ k) = xxx1(i,1,lphi_down + k) + down_recv_phi(i,1,k,1)

        xxx2(i,1,lphi_down+ k) = xxx2(i,1,lphi_down + k) + down_recv_phi(i,1,k,2)

        xxx3(i,1,lphi_down+ k) = xxx3(i,1,lphi_down + k) + down_recv_phi(i,1,k,3)

        xxx4(i,1,lphi_down+ k) = xxx4(i,1,lphi_down + k) + down_recv_phi(i,1,k,4)

        xxx5(i,1,lphi_down+ k) = xxx5(i,1,lphi_down + k) + down_recv_phi(i,1,k,5)

        xxx6(i,1,lphi_down+ k) = xxx6(i,1,lphi_down + k) + down_recv_phi(i,1,k,6)
      end do
      end do

#ifdef AMDGPU
!$omp target teams distribute parallel do collapse(2) private(k,i)
#else
!$omp parallel do
#endif
      do k = 1, lphi_shd2
      do i = 1, lrz
        xxx1(i,1,lphi_up  + k) = xxx1(i,1,lphi_up   + k) + up_recv_phi(i,1,k,1)

        xxx2(i,1,lphi_up  + k) = xxx2(i,1,lphi_up   + k) + up_recv_phi(i,1,k,2)

        xxx3(i,1,lphi_up  + k) = xxx3(i,1,lphi_up   + k) + up_recv_phi(i,1,k,3)

        xxx4(i,1,lphi_up  + k) = xxx4(i,1,lphi_up   + k) + up_recv_phi(i,1,k,4)

        xxx5(i,1,lphi_up  + k) = xxx5(i,1,lphi_up   + k) + up_recv_phi(i,1,k,5)

        xxx6(i,1,lphi_up  + k) = xxx6(i,1,lphi_up   + k) + up_recv_phi(i,1,k,6)
      end do
      end do


!------------------------------------------------------
! communication in z-direction
      if(mpi_proc_z.ge.2)then

      node_up = mod(my_rank_z + 1, mpi_proc_z)
      node_down = mod(my_rank_z - 1 + mpi_proc_z, mpi_proc_z)

      if(my_rank_z.ne.(mpi_proc_z-1))then
#ifdef AMDGPU
!$omp target teams distribute parallel do collapse(3) private(k,j,i)
#else
!$omp parallel do
#endif
        do k = 1, lphi
        do j = 1, lz_shd2
        do i = 1, lr
          up_send_z(i,j,k,1) = xxx1(i,lz_up + j,k)
          up_send_z(i,j,k,2) = xxx2(i,lz_up + j,k)
          up_send_z(i,j,k,3) = xxx3(i,lz_up + j,k)
          up_send_z(i,j,k,4) = xxx4(i,lz_up + j,k)
          up_send_z(i,j,k,5) = xxx5(i,lz_up + j,k)
          up_send_z(i,j,k,6) = xxx6(i,lz_up + j,k)
        end do
        end do
        end do
      end if

      if(my_rank_z.ne.0 )then
#ifdef AMDGPU
!$omp target teams distribute parallel do collapse(3) private(k,j,i)
#else
!$omp parallel do
#endif
        do k = 1, lphi
        do j = 1, lz_shd2
        do i = 1, lr
          down_send_z(i,j,k,1) = xxx1(i,lz_down + j,k)
          down_send_z(i,j,k,2) = xxx2(i,lz_down + j,k)
          down_send_z(i,j,k,3) = xxx3(i,lz_down + j,k)
          down_send_z(i,j,k,4) = xxx4(i,lz_down + j,k)
          down_send_z(i,j,k,5) = xxx5(i,lz_down + j,k)
          down_send_z(i,j,k,6) = xxx6(i,lz_down + j,k)
        end do
        end do
        end do
      end if

! communication of particle data

      lrphi2 = lr*lphi*lz_shd2*mlt

      if(my_rank_z.eq.0)then
        call mpi_isend(up_send_z(1,1,1,1),lrphi2,mpi_real8,node_up,1 &
                      ,z_world,my_request(1),mpi_err)
        call mpi_irecv(up_recv_z(1,1,1,1),lrphi2,mpi_real8,node_up,1 &
                      ,z_world,my_request(2),mpi_err)

        call mpi_waitall(2,my_request,my_status_all,mpi_err)

#ifdef AMDGPU
!$omp target teams distribute parallel do collapse(3) private(k,j,i)
#else
!$omp parallel do
#endif
        do k = 1, lphi
        do j = 1, lz_shd2
        do i = 1, lr
          xxx1(i,lz_up+j,k) = xxx1(i,lz_up+j,k) + up_recv_z(i,j,k,1)
          xxx2(i,lz_up+j,k) = xxx2(i,lz_up+j,k) + up_recv_z(i,j,k,2)
          xxx3(i,lz_up+j,k) = xxx3(i,lz_up+j,k) + up_recv_z(i,j,k,3)
          xxx4(i,lz_up+j,k) = xxx4(i,lz_up+j,k) + up_recv_z(i,j,k,4)
          xxx5(i,lz_up+j,k) = xxx5(i,lz_up+j,k) + up_recv_z(i,j,k,5)
          xxx6(i,lz_up+j,k) = xxx6(i,lz_up+j,k) + up_recv_z(i,j,k,6)
        end do
        end do
        end do

      else if(my_rank_z.eq.(mpi_proc_z-1) )then
        call mpi_isend(down_send_z(1,1,1,1),lrphi2,mpi_real8,node_down,1 &
                      ,z_world,my_request(1),mpi_err)
        call mpi_irecv(down_recv_z(1,1,1,1),lrphi2,mpi_real8,node_down,1 &
                      ,z_world,my_request(2),mpi_err)

        call mpi_waitall(2,my_request,my_status_all,mpi_err)

#ifdef AMDGPU
!$omp target teams distribute parallel do collapse(3) private(k,j,i)
#else
!$omp parallel do
#endif
        do k = 1, lphi
        do j = 1, lz_shd2
        do i = 1, lr
          xxx1(i,lz_down+j,k) = xxx1(i,lz_down+j,k) + down_recv_z(i,j,k,1)
          xxx2(i,lz_down+j,k) = xxx2(i,lz_down+j,k) + down_recv_z(i,j,k,2)
          xxx3(i,lz_down+j,k) = xxx3(i,lz_down+j,k) + down_recv_z(i,j,k,3)
          xxx4(i,lz_down+j,k) = xxx4(i,lz_down+j,k) + down_recv_z(i,j,k,4)
          xxx5(i,lz_down+j,k) = xxx5(i,lz_down+j,k) + down_recv_z(i,j,k,5)
          xxx6(i,lz_down+j,k) = xxx6(i,lz_down+j,k) + down_recv_z(i,j,k,6)
        end do
        end do
        end do

      else
        call mpi_isend(up_send_z(1,1,1,1),lrphi2,mpi_real8,node_up,1 &
                      ,z_world,my_request(1),mpi_err)
        call mpi_isend(down_send_z(1,1,1,1),lrphi2,mpi_real8,node_down,1 &
                      ,z_world,my_request(2),mpi_err)
        call mpi_irecv(up_recv_z(1,1,1,1),lrphi2,mpi_real8,node_up,1 &
                      ,z_world,my_request(3),mpi_err)
        call mpi_irecv(down_recv_z(1,1,1,1),lrphi2,mpi_real8,node_down,1 &
                      ,z_world,my_request(4),mpi_err)

        call mpi_waitall(4,my_request,my_status_all,mpi_err)

#ifdef AMDGPU
!$omp target teams distribute parallel do collapse(3) private(k,j,i)
#else
!$omp parallel do
#endif
        do k = 1, lphi
        do j = 1, lz_shd2
        do i = 1, lr
          xxx1(i,lz_down+j,k)= xxx1(i,lz_down+j,k) + down_recv_z(i,j,k,1)

          xxx2(i,lz_down+j,k)= xxx2(i,lz_down+j,k) + down_recv_z(i,j,k,2)

          xxx3(i,lz_down+j,k)= xxx3(i,lz_down+j,k) + down_recv_z(i,j,k,3)

          xxx4(i,lz_down+j,k)= xxx4(i,lz_down+j,k) + down_recv_z(i,j,k,4)

          xxx5(i,lz_down+j,k)= xxx5(i,lz_down+j,k) + down_recv_z(i,j,k,5)

          xxx6(i,lz_down+j,k)= xxx6(i,lz_down+j,k) + down_recv_z(i,j,k,6)
        end do
        end do
        end do

#ifdef AMDGPU
!$omp target teams distribute parallel do collapse(3) private(k,j,i)
#else
!$omp parallel do
#endif
        do k = 1, lphi
        do j = 1, lz_shd2
        do i = 1, lr
          xxx1(i,lz_up+j,k)  = xxx1(i,lz_up+j,k)   + up_recv_z(i,j,k,1)

          xxx2(i,lz_up+j,k)  = xxx2(i,lz_up+j,k)   + up_recv_z(i,j,k,2)

          xxx3(i,lz_up+j,k)  = xxx3(i,lz_up+j,k)   + up_recv_z(i,j,k,3)

          xxx4(i,lz_up+j,k)  = xxx4(i,lz_up+j,k)   + up_recv_z(i,j,k,4)

          xxx5(i,lz_up+j,k)  = xxx5(i,lz_up+j,k)   + up_recv_z(i,j,k,5)

          xxx6(i,lz_up+j,k)  = xxx6(i,lz_up+j,k)   + up_recv_z(i,j,k,6)
        end do
        end do
        end do

      end if

      end if  !end of z-communication

!------------------------------------------------------
! communication in r-direction
      if(mpi_proc_r.ge.2)then

      node_up = mod(my_rank_r + 1, mpi_proc_r)
      node_down = mod(my_rank_r - 1 + mpi_proc_r, mpi_proc_r)

      if(my_rank_r.ne.(mpi_proc_r-1))then
#ifdef AMDGPU
!$omp target teams distribute parallel do collapse(2) private(j,i)
#else
!$omp parallel do
#endif         
        do j = 1, lzphi
        do i = 1, lr_shd2
          up_send_r(i,j,1,1) = xxx1(lr_up+i,j,1)
          up_send_r(i,j,1,2) = xxx2(lr_up+i,j,1)
          up_send_r(i,j,1,3) = xxx3(lr_up+i,j,1)
          up_send_r(i,j,1,4) = xxx4(lr_up+i,j,1)
          up_send_r(i,j,1,5) = xxx5(lr_up+i,j,1)
          up_send_r(i,j,1,6) = xxx6(lr_up+i,j,1)
        end do
        end do
      end if

      if(my_rank_r.ne.0 )then
#ifdef AMDGPU
!$omp target teams distribute parallel do collapse(2) private(j,i)
#else
!$omp parallel do
#endif
        do j = 1, lzphi
        do i = 1, lr_shd2
          down_send_r(i,j,1,1) = xxx1(lr_down+i,j,1)
          down_send_r(i,j,1,2) = xxx2(lr_down+i,j,1)
          down_send_r(i,j,1,3) = xxx3(lr_down+i,j,1)
          down_send_r(i,j,1,4) = xxx4(lr_down+i,j,1)
          down_send_r(i,j,1,5) = xxx5(lr_down+i,j,1)
          down_send_r(i,j,1,6) = xxx6(lr_down+i,j,1)
        end do
        end do
      end if


! communication of particle data

      lzphi2 = lz*lphi*lr_shd2*mlt

      if(my_rank_r.eq.0)then
        call mpi_isend(up_send_r(1,1,1,1),lzphi2,mpi_real8,node_up,1 &
                      ,r_world,my_request(1),mpi_err)
        call mpi_irecv(up_recv_r(1,1,1,1),lzphi2,mpi_real8,node_up,1 &
                      ,r_world,my_request(2),mpi_err)

        call mpi_waitall(2,my_request,my_status_all,mpi_err)

#ifdef AMDGPU
!$omp target teams distribute parallel do collapse(2) private(j,i)
#else
!$omp parallel do
#endif
        do j = 1, lzphi
        do i = 1, lr_shd2
          xxx1(lr_up+i,j,1) = xxx1(lr_up+i,j,1) + up_recv_r(i,j,1,1)
          xxx2(lr_up+i,j,1) = xxx2(lr_up+i,j,1) + up_recv_r(i,j,1,2)
          xxx3(lr_up+i,j,1) = xxx3(lr_up+i,j,1) + up_recv_r(i,j,1,3)
          xxx4(lr_up+i,j,1) = xxx4(lr_up+i,j,1) + up_recv_r(i,j,1,4)
          xxx5(lr_up+i,j,1) = xxx5(lr_up+i,j,1) + up_recv_r(i,j,1,5)
          xxx6(lr_up+i,j,1) = xxx6(lr_up+i,j,1) + up_recv_r(i,j,1,6)
        end do
        end do

      else if(my_rank_r.eq.(mpi_proc_r-1) )then
        call mpi_isend(down_send_r(1,1,1,1),lzphi2,mpi_real8,node_down,1 &
                      ,r_world,my_request(1),mpi_err)
        call mpi_irecv(down_recv_r(1,1,1,1),lzphi2,mpi_real8,node_down,1 &
                      ,r_world,my_request(2),mpi_err)

        call mpi_waitall(2,my_request,my_status_all,mpi_err)

#ifdef AMDGPU
!$omp target teams distribute parallel do collapse(2) private(j,i)
#else
!$omp parallel do
#endif
        do j = 1, lzphi
        do i = 1, lr_shd2
          xxx1(lr_down+i,j,1) = xxx1(lr_down+i,j,1) + down_recv_r(i,j,1,1)
          xxx2(lr_down+i,j,1) = xxx2(lr_down+i,j,1) + down_recv_r(i,j,1,2)
          xxx3(lr_down+i,j,1) = xxx3(lr_down+i,j,1) + down_recv_r(i,j,1,3)
          xxx4(lr_down+i,j,1) = xxx4(lr_down+i,j,1) + down_recv_r(i,j,1,4)
          xxx5(lr_down+i,j,1) = xxx5(lr_down+i,j,1) + down_recv_r(i,j,1,5)
          xxx6(lr_down+i,j,1) = xxx6(lr_down+i,j,1) + down_recv_r(i,j,1,6)
        end do
        end do

      else
        call mpi_isend(up_send_r(1,1,1,1),lzphi2,mpi_real8,node_up,1 &
                      ,r_world,my_request(1),mpi_err)
        call mpi_isend(down_send_r(1,1,1,1),lzphi2,mpi_real8,node_down,1 &
                      ,r_world,my_request(2),mpi_err)
        call mpi_irecv(up_recv_r(1,1,1,1),lzphi2,mpi_real8,node_up,1 &
                      ,r_world,my_request(3),mpi_err)
        call mpi_irecv(down_recv_r(1,1,1,1),lzphi2,mpi_real8,node_down,1 &
                      ,r_world,my_request(4),mpi_err)

        call mpi_waitall(4,my_request,my_status_all,mpi_err)

#ifdef AMDGPU
!$omp target teams distribute parallel do collapse(2) private(j,i)
#else
!$omp parallel do
#endif
        do j = 1, lzphi
        do i = 1, lr_shd2
          xxx1(lr_down+i,j,1)= xxx1(lr_down+i,j,1) + down_recv_r(i,j,1,1)

          xxx2(lr_down+i,j,1)= xxx2(lr_down+i,j,1) + down_recv_r(i,j,1,2)

          xxx3(lr_down+i,j,1)= xxx3(lr_down+i,j,1) + down_recv_r(i,j,1,3)

          xxx4(lr_down+i,j,1)= xxx4(lr_down+i,j,1) + down_recv_r(i,j,1,4)

          xxx5(lr_down+i,j,1)= xxx5(lr_down+i,j,1) + down_recv_r(i,j,1,5)

          xxx6(lr_down+i,j,1)= xxx6(lr_down+i,j,1) + down_recv_r(i,j,1,6)
        end do
        end do

#ifdef AMDGPU
!$omp target teams distribute parallel do collapse(2) private(j,i)
#else
!$omp parallel do
#endif
        do j = 1, lzphi
        do i = 1, lr_shd2
          xxx1(lr_up+i,j,1)  = xxx1(lr_up+i,j,1)   + up_recv_r(i,j,1,1)

          xxx2(lr_up+i,j,1)  = xxx2(lr_up+i,j,1)   + up_recv_r(i,j,1,2)

          xxx3(lr_up+i,j,1)  = xxx3(lr_up+i,j,1)   + up_recv_r(i,j,1,3)

          xxx4(lr_up+i,j,1)  = xxx4(lr_up+i,j,1)   + up_recv_r(i,j,1,4)

          xxx5(lr_up+i,j,1)  = xxx5(lr_up+i,j,1)   + up_recv_r(i,j,1,5)

          xxx6(lr_up+i,j,1)  = xxx6(lr_up+i,j,1)   + up_recv_r(i,j,1,6)
        end do
        end do

      end if

      end if  !end of r-communication

!--------------------------------

! summation over ptcl_world
!2012-06-17
      if(mpi_proc_ptcl.ge.2)then
#ifdef AMDGPU
!$omp target teams distribute parallel do collapse(3) private(k,j,i)
#else
!$omp parallel do
#endif
        do k = 1, lphi
        do j = 1, lz
        do i = 1, lr
          xxx_com(i,j,k,1) = xxx1(i,j,k)
          xxx_com(i,j,k,2) = xxx2(i,j,k)
          xxx_com(i,j,k,3) = xxx3(i,j,k)
          xxx_com(i,j,k,4) = xxx4(i,j,k)
          xxx_com(i,j,k,5) = xxx5(i,j,k)
          xxx_com(i,j,k,6) = xxx6(i,j,k)
        end do
        end do
        end do

        lsize = lrzphi*mlt
        call mpi_allreduce(mpi_in_place,xxx_com(1,1,1,1),lsize,mpi_real8 &
                          ,mpi_sum,ptcl_world,mpi_err)

#ifdef AMDGPU
!$omp target teams distribute parallel do collapse(3) private(k,j,i)
#else
!$omp parallel do
#endif
        do k = 1, lphi
        do j = 1, lz
        do i = 1, lr
          xxx1(i,j,k) = xxx_com(i,j,k,1)
          xxx2(i,j,k) = xxx_com(i,j,k,2)
          xxx3(i,j,k) = xxx_com(i,j,k,3)
          xxx4(i,j,k) = xxx_com(i,j,k,4)
          xxx5(i,j,k) = xxx_com(i,j,k,5)
          xxx6(i,j,k) = xxx_com(i,j,k,6)
        end do
        end do
        end do
      end if
!2012-06-17 end


end
!--------------------------------------------------------------------
subroutine end_ransuu !2020-07-28
! end random number generator
!--------------------------------------------------------------------
      use mpiset
      use random_num

#ifdef HITACHI
      call hmatexit(ierr)
      deallocate(lis)
      deallocate(lstpu)

!2020-07-28s
#elif AURORA 
      ! Generator Finalization
      call asl_random_destroy(rng)

      ! Library Finalization
      call asl_library_finalize()
!2020-07-28e

#endif

end
!--------------------------------------------------------------------
subroutine ransuu(nrandom, rand_no)
!2012-01-28 with random_number routine of fortran95
!2013-01-30 with mt_stream (Multiple Stream Mersenne Twister):MSMT
! random number generation
!--------------------------------------------------------------------
      use mpiset
      use random_num

      implicit none
      integer::nrandom
      real(8)::rand_no(*)
      integer::node_up,node_down,i
      integer::m,n !2021-02-06
      integer::k_rank !2017-01-23
      integer::ierr

#ifdef HITACHI
      if (iflag==0) then
        allocate(lis(nprocess))
        allocate(lstpu(nprocess))
        iwksize = max(16,(nprocess+7))*8

        ix = 0
        is = -1
        lis(1) = -1
        lstpu(1) = -1
        iopt1 = 0
        iopt2(1) = 1
        iopt2(2) = 1
        call hmatinit(iwksize,lstpu,nprocess,ierr)
        call hdru4mdp(nrandom,ix,is,lis,lstpu,nprocess,iopt1,iopt2 &
                     ,mtbl,mtblp,rand_no,ierr)
        iopt1 = 2

        iflag=1
      end if

      call hdru4mdp(nrandom,ix,is,lis,lstpu,nprocess,iopt1,iopt2 &
                   ,mtbl,mtblp,rand_no,ierr)

#elif MKL
        if (iflag==0) then
          my_rank_seed = mod(my_rank, 4096)
          brng = VSL_BRNG_MT2203 + my_rank_seed
          ierr = vslNewStream(stream, brng, iseed)
          iflag=1
        end if

        lb=0.d0
        rb=1.d0
        ierr = vdrnguniform(method, stream, nrandom, rand_no, lb, rb) !2020-02-15

#elif MSMT
        if (iflag==0) then
          call set_mt19937
          call new(mts)
          call init(mts,iseed)
          id = my_rank

          if(my_rank/=0)then
            call create_stream(mts,mts_new,id)
            mts = mts_new
          end if

          iflag=1
        end if

!$omp parallel do
        do i=1,nrandom
          rand_no(i) = genrand_double1(mts)
        end do

#elif DCMT
      w = 32
      p = 521
!     p = 607
!     p = 1279
!     p = 2203
      seed0 = abs(4172*my_rank)
      seed1 = abs(1234*my_rank)
      id = my_rank
      if (iflag==0) call get_mt_parameter_id_st_f(w,p,id,seed0)
      if (iflag==0) call sgenrand_mt_f(seed1)
      iflag=1
!$omp parallel do
      do i=1,nrandom
         call genrand_mt_f(rand_no(i))
      end do

#elif FUJITSU
! Fujitsu START
!2021-02-06s
      nbuf = nrandom / NUMT + min(1, mod(nrandom, NUMT) ) 

      NWORK = NMAX
      if (iflag==0) then
        IX = iseed+my_rank
        CALL DM_VRANU4(IX,DA,NRNDMAX,NBUF,DWORK,NWORK,ICON)
        iflag=1
      else
        CALL DM_VRANU4(IX,DA,NRNDMAX,NBUF,DWORK,NWORK,ICON)
      end if

!$omp parallel do private(m,i)
      do n = 1, nrandom
        m = n / nbuf + 1
        i = mod(n, nbuf) + 1
        rand_no(n) = DA(i, m)
      end do

!2021-02-06e
! Fujitsu END

#elif AURORA
!2020-07-28s
      if (iflag==0) then

        ! Library Initialization
        call asl_library_initialize()

        ! Generator Preparation
        call asl_random_parallel_create(rng, ASL_RANDOMMETHOD_MT19937_64, mpi_comm_world)

        iflag = 1
      end if
!2020-07-28e

      ! Generation
      call asl_random_generate_d(rng, nrandom, rand_no)

!2020-07-28s, moved to end_ransuu
      ! Generator Finalization
!      call asl_random_destroy(rng)

      ! Library Finalization
!      call asl_library_finalize()

!2020-07-28e
      
#elif AMDGPU
      if (iflag==0) then
        call hiprand_init()
        iflag = 1

        allocate(rand_no_f95(nrandom)) !corrected 2025-09-06
      end if


      do k_rank = 0, my_rank
        call hiprand_gen(nrandom, rand_no_f95)
      end do

!$omp target teams distribute parallel do private(i)
      do i = 1, nrandom
        rand_no(i) = rand_no_f95(i)
      end do

      do k_rank = my_rank + 1, nprocess-1
        call hiprand_gen(nrandom, rand_no_f95)
      end do

!Fortran95 version, 2017-01-23
#else
!2017-09-13s
      if (iflag==0) then
        call random_seed(size=nsize)
        allocate(seed(nsize))

        iflag = 1
      end if

      allocate(rand_no_f95(nrandom))

      do k_rank = 0, my_rank
        call random_number(rand_no_f95)
!        call random_seed(get=seed)
!        call random_seed(put=seed)
      end do

!$omp parallel do
      do i = 1, nrandom
        rand_no(i) = rand_no_f95(i)
      end do

      do k_rank = my_rank + 1, nprocess-1
        call random_number(rand_no_f95)
!        call random_seed(get=seed)
!        call random_seed(put=seed)
      end do

      deallocate(rand_no_f95)
!2017-09-13e

#endif

end
!--------------------------------------------------------------------
subroutine order_particle(marker_num,gc,node)
! 2015-06-23
!--------------------------------------------------------------------
      use mpiset
      use field
      use grid
      implicit none

#ifdef OPEN_MP
      integer, parameter::marker_thread_max=marker_local_max/lpara
      integer, parameter::marker_each_thread=marker_each/lpara
#endif

      integer::marker_num
      real(8)::gc(ngc2,marker_num)
      integer::node(marker_num)

      real(8)::gc_buf(ngc2,marker_each)
      integer::node_buf(marker_each)

      real(8)::dr1,dz1,dphi1
      integer::n,ia,ja,ka,m,nb
      integer::kr,kz,kphi
      real(8)::ma_mi_r
      integer::igc

#ifdef OPEN_MP
      integer(4) :: me,nthreads
      integer(4),external :: omp_get_thread_num,omp_get_num_threads
      integer::marker_index(marker_thread_max,0:lpara-1,lr,lz,lphi),imrk_grid(0:lpara-1,lr,lz,lphi)
      integer::marker_index_dead(marker_each_thread,0:lpara-1),imrk_dead(0:lpara-1)
      integer::nb_thread(0:lpara-1)
      integer(4) :: nmod,nvec0,nvec,nn,l
#else
      integer::marker_index(marker_local_max,lr,lz,lphi),imrk_grid(lr,lz,lphi)
      integer::marker_index_dead(marker_each),imrk_dead
#endif

      dr1 = 1.0d0/dr
      dz1 = 1.0d0/dz
      dphi1 = 1.0d0/dphi

      kr   = 1 - kr_offset(my_rank)
      kz   = 1 - kz_offset(my_rank)
      kphi = 1 + lphi_shd - kphi_offset(my_rank)
!      kphi = 3 - kphi_offset(my_rank)
      ma_mi_r = major_r-minor_r

!$omp parallel
!$omp workshare
      imrk_grid = 0
      imrk_dead = 0
!$omp end workshare
!$omp end parallel

#ifdef OPEN_MP
      nb_thread = 0

! !$omp parallel private(nthreads,me,nvec,nn,n,ia,ja,ka,nb,l)
!$omp parallel private(nthreads,me,nmod,nvec0,nvec,nn,n,ia,ja,ka,nb,l,m)
      nthreads=omp_get_num_threads()
      me=omp_get_thread_num()

      nmod=mod(marker_num,nthreads)
      nvec0=marker_num/nthreads
      nvec=nvec0+min(1,nmod/(me+1))
      do nn = 1, nvec
        n = nvec0*me+min(nmod,me)+nn

        ia=max(1,min(lr  -1,int((gc(1,n)-ma_mi_r)*dr1  ) + kr  ))
        ja=max(1,min(lz  -1,int(gc(2,n)          *dz1  ) + kz  ))
        ka=max(1,min(lphi-1,int(gc(3,n)          *dphi1) + kphi))

        if(imrk_grid(me,ia,ja,ka).lt.marker_thread_max)then
          imrk_grid(me,ia,ja,ka) = imrk_grid(me,ia,ja,ka) + 1
          marker_index(imrk_grid(me,ia,ja,ka),me,ia,ja,ka) = n
        else
          imrk_dead(me) = imrk_dead(me) + 1
          marker_index_dead(imrk_dead(me),me) = n
        end if
      end do

      nb_thread(me) = nvec
!$omp barrier

      nb = 0
      do l = 1, me
        nb = nb + nb_thread(l-1)
      end do

      do ka = 1, lphi-1
      do ja = 1, lz-1
      do ia = 1, lr-1
      do m = 1, imrk_grid(me,ia,ja,ka)
        nb = nb + 1
        n = marker_index(m,me,ia,ja,ka)

        do igc = 1, ngc2
          gc_buf(igc,nb) = gc(igc,n)
        end do
        node_buf(nb) = node(n)

      end do
      end do
      end do
      end do

! for inactive particles
      do m = 1, imrk_dead(me)
        nb = nb + 1
        n = marker_index_dead(m,me)

        do igc = 1, ngc2
          gc_buf(igc,nb) = gc(igc,n)
        end do
        node_buf(nb) = node(n)

      end do

      do nn = 1, nvec
        n = nvec0*me+min(nmod,me)+nn

        do igc = 1, ngc2
          gc(igc,n) = gc_buf(igc,n)
        end do
        node(n) = node_buf(n)

      end do
!$omp end parallel

#else
      do n = 1, marker_num
        ia=max(1,min(lr  -1,int((gc(1,n)-ma_mi_r)*dr1  ) + kr  ))
        ja=max(1,min(lz  -1,int(gc(2,n)          *dz1  ) + kz  ))
        ka=max(1,min(lphi-1,int(gc(3,n)          *dphi1) + kphi))

!        if(node(n).eq.my_rank)then
          if(imrk_grid(ia,ja,ka).lt.marker_local_max)then
            imrk_grid(ia,ja,ka) = imrk_grid(ia,ja,ka) + 1
            marker_index(imrk_grid(ia,ja,ka),ia,ja,ka) = n
          else
            imrk_dead = imrk_dead + 1          !2011-11-24
            marker_index_dead(imrk_dead) = n   !2011-11-24
          end if
!        else
!          imrk_dead = imrk_dead + 1 
!          marker_index_dead(imrk_dead) = n
!        end if
      end do

      nb = 0
      do ka = 1, lphi-1
      do ja = 1, lz-1
      do ia = 1, lr-1
      do m = 1, imrk_grid(ia,ja,ka)
        nb = nb + 1
        n = marker_index(m,ia,ja,ka)

        do igc = 1, ngc2
          gc_buf(igc,nb) = gc(igc,n)
        end do
        node_buf(nb) = node(n)

      end do
      end do
      end do
      end do

! for inactive particles
      do m = 1, imrk_dead
        nb = nb + 1
        n = marker_index_dead(m)

        do igc = 1, ngc2
          gc_buf(igc,nb) = gc(igc,n)
        end do
        node_buf(nb) = node(n)

      end do

!$omp parallel do
      do n = 1, marker_num

        do igc = 1, ngc2
          gc(igc,n) = gc_buf(igc,n)
        end do
        node(n) = node_buf(n)

      end do
#endif

end
!--------------------------------------------------------------------
subroutine com_particle3(idirection, &
                 marker_num,gc,node,node_now)
!     communication of particle data in 3-dimensional domain-decomp.
!     idirection=0: r-direction
!     idirection=1: z-direction
!     idirection=2: phi-direction
!     OPEN_MP: 2013-07-27
! 2015-06-23
!--------------------------------------------------------------------
      use mpiset
      use grid
      implicit none

      integer::marker_num
      real(8)::gc(ngc2,marker_num)
      integer::node(marker_num),node_now(marker_num)

      integer::n,idirection,node_up,node_down,i,nloop
      integer::nup,ndown,nup_r,ndown_r
      integer::nup8,ndown8,nup_r8,ndown_r8
      integer::nsend_up(ncomm),nsend_down(ncomm)
      real(8)::send_up(ngc2,ncomm),send_down(ngc2,ncomm) &
              ,recv_up(ngc2,ncomm),recv_down(ngc2,ncomm) !2016-01-05
      integer::each_world,node_up_each,node_down_each !2012-06-22
      integer::igc

#ifdef OPEN_MP
      integer(4) :: me,nthreads
      integer(4),external :: omp_get_thread_num,omp_get_num_threads
      integer(4) :: l,nn,nup0,ndown0,nvec,nvec0,nmod
!      integer(4),parameter :: maxthreads=16
      integer(4),parameter :: maxthreads=lpara !2013-06-26
      integer::nsend_up0(ncomm,0:maxthreads-1) &
              ,nsend_down0(ncomm,0:maxthreads-1)
      integer::ia(marker_num),na(marker_num)
      integer::nnmax

#elif AURORA
      integer(4) :: nn

#endif

! check particles communicated

      if(idirection.eq.0)then
!2012-06-17
!2012-06-22
        each_world = r_world
        node_up_each = mod(my_rank_r + 1, mpi_proc_r)
        node_down_each = mod(my_rank_r - 1 + mpi_proc_r, mpi_proc_r)
!2012-06-22 end

        node_up = mod(my_rank_mhd + 1, mpi_proc_mhd) &
                + my_rank_ptcl*mpi_proc_mhd
        node_down = mod(my_rank_mhd - 1 + mpi_proc_mhd, mpi_proc_mhd) &
                  + my_rank_ptcl*mpi_proc_mhd
!2012-06-17 end

#ifdef AMDGPU        
!$omp target teams distribute parallel do private(n)
#else
!$omp parallel do
#endif
        do n = 1, marker_num
          node_now(n) = max(0, min(mpi_proc_r-1, &
               int(dble(mpi_proc_r)*(gc(1,n)-(major_r-minor_r) )/rleng ) &
                           ) ) &
                      + my_rank_z*mpi_proc_r &
                      + my_rank_phi*mpi_proc_pol &
                      + my_rank_ptcl*mpi_proc_mhd  !2012-06-17
        end do
      else if(idirection.eq.1)then
!2012-06-17
!2012-06-22
        each_world = z_world
        node_up_each = mod(my_rank_z + 1, mpi_proc_z)
        node_down_each = mod(my_rank_z - 1 + mpi_proc_z, mpi_proc_z)
!2012-06-22 end

        node_up = mod(my_rank_mhd + mpi_proc_r, mpi_proc_mhd) &
                + my_rank_ptcl*mpi_proc_mhd
        node_down = mod(my_rank_mhd - mpi_proc_r + mpi_proc_mhd, mpi_proc_mhd) &
                  + my_rank_ptcl*mpi_proc_mhd
!2012-06-17 end

#ifdef AMDGPU        
!$omp target teams distribute parallel do private(n)
#else
!$omp parallel do
#endif
        do n = 1, marker_num
          node_now(n) = max(0, min(mpi_proc_z-1, &
               int(dble(mpi_proc_z)*gc(2,n)/zleng ) &
                           ) )*mpi_proc_r &
                      + my_rank_r &
                      + my_rank_phi*mpi_proc_pol &
                      + my_rank_ptcl*mpi_proc_mhd  !2012-06-17
        end do
      else if(idirection.eq.2)then
!2012-06-17
!2012-06-22
        each_world = phi_world
        node_up_each = mod(my_rank_phi + 1, mpi_proc_phi)
        node_down_each = mod(my_rank_phi - 1 + mpi_proc_phi &
                            ,mpi_proc_phi)
!2012-06-22 end

        node_up = mod(my_rank_mhd + mpi_proc_pol, mpi_proc_mhd) &
                + my_rank_ptcl*mpi_proc_mhd
        node_down = mod(my_rank_mhd - mpi_proc_pol + mpi_proc_mhd &
                       ,mpi_proc_mhd) &
                  + my_rank_ptcl*mpi_proc_mhd
!2012-06-17 end

#ifdef AMDGPU        
!$omp target teams distribute parallel do private(n)
#else
!$omp parallel do
#endif
        do n = 1, marker_num
          node_now(n) = max(0, min(mpi_proc_phi-1, & 
               int(dble(mpi_proc_phi)*gc(3,n)/phileng ) &
                           ) )*mpi_proc_pol &
                      + my_rank_z*mpi_proc_r &
                      + my_rank_r            &
                      + my_rank_ptcl*mpi_proc_mhd  !2012-06-17
        end do
      end if


      nup = 0
      ndown = 0
      nup_r = 0
      ndown_r = 0



#ifdef AMDGPU
      do n = 1, marker_num
        if(node(n).eq.my_rank.and.node_now(n).eq.node_up)then
          node(n) = node_up
          nup = nup + 1
          nsend_up(nup) = n

!!$omp target teams distribute parallel do private(igc)
          do igc = 1, ngc2
            send_up(igc,nup) = gc(igc,n) 
          end do
          gc(7,n) = 0.0d0

       else if(node(n).eq.my_rank.and.node_now(n).eq.node_down)then
          node(n) = node_down
          ndown = ndown + 1
          nsend_down(ndown) = n

!!$omp target teams distribute parallel do private(igc)
          do igc = 1, ngc2
            send_down(igc,ndown) = gc(igc,n) 
          end do
          gc(7,n) = 0.0d0

        end if
      end do

#elif OPEN_MP
!$omp parallel private(me,nup0,ndown0,nvec,nvec0,nmod,n,l,nthreads)
      nthreads=omp_get_num_threads()
      me=omp_get_thread_num()
      nup0=0
      ndown0=0
      nmod=mod(marker_num,nthreads)
      nvec0=marker_num/nthreads
      nvec=nvec0+min(1,nmod/(me+1))
      do nn = 1, nvec
        n = nvec0*me+min(nmod,me)+nn
        if(node(n).eq.my_rank.and. &
          (node_now(n).eq.node_up) )then
          nup0 = nup0 + 1
          nsend_up0(nup0,me) = n
        else if(node(n).eq.my_rank.and. &
               (node_now(n).eq.node_down) )then
          ndown0 = ndown0 + 1
          nsend_down0(ndown0,me) = n
        end if
      end do
      do l=0,nthreads-1
         if (l==me) then
            do n=1,nup0
               nup=nup+1
               nsend_up(nup)=nsend_up0(n,l)
            end do
            do n=1,ndown0
               ndown=ndown+1
               nsend_down(ndown)=nsend_down0(n,l)
            end do
          end if
!$omp barrier
      end do
!$omp end parallel

!2016-01-10s
      if(nup.gt.ncomm)then
        write(7,*)'nup > ncomm, rank=',my_rank
        go to 9999
      end if

      if(ndown.gt.ncomm)then
        write(7,*)'ndown > ncomm, rank=',my_rank
        go to 9999
      end if
!2016-01-10e


!$omp parallel do private(n)
      do nn = 1, nup
        n = nsend_up(nn)
        node(n) = node_up

        do igc = 1, ngc2
          send_up(igc,nn) = gc(igc,n) 
        end do
        gc(7,n) = 0.0d0

      end do
!$omp parallel do private(n)
      do nn = 1, ndown
        n = nsend_down(nn)
        node(n) = node_down

        do igc = 1, ngc2
          send_down(igc,nn) = gc(igc,n) 
        end do
        gc(7,n) = 0.0d0

      end do

!NEC start
#elif AURORA
      do n = 1, marker_num
        if(node(n).eq.my_rank.and. &
          (node_now(n).eq.node_up) )then
          nup = nup + 1
          nsend_up(nup) = n
        else if(node(n).eq.my_rank.and. &
               (node_now(n).eq.node_down) )then
          ndown = ndown + 1
          nsend_down(ndown) = n
        end if
      end do

      if(nup.gt.ncomm)then
        write(7,*)'nup > ncomm, rank=',my_rank
        go to 9999
      end if

      if(ndown.gt.ncomm)then
        write(7,*)'ndown > ncomm, rank=',my_rank
        go to 9999
      end if

      do nn = 1, nup
        n = nsend_up(nn)
        node(n) = node_up

        do igc = 1, ngc2
          send_up(igc,nn) = gc(igc,n) 
        end do
        gc(7,n) = 0.0d0

      end do

      do nn = 1, ndown
        n = nsend_down(nn)
        node(n) = node_down

        do igc = 1, ngc2
          send_down(igc,nn) = gc(igc,n) 
        end do
        gc(7,n) = 0.0d0

      end do
!NEC end

#else
      do n = 1, marker_num
        if(node(n).eq.my_rank.and.node_now(n).eq.node_up)then
          node(n) = node_up
          nup = nup + 1
          nsend_up(nup) = n

          do igc = 1, ngc2
            send_up(igc,nup) = gc(igc,n) 
          end do
          gc(7,n) = 0.0d0

       else if(node(n).eq.my_rank.and.node_now(n).eq.node_down)then
          node(n) = node_down
          ndown = ndown + 1
          nsend_down(ndown) = n

          do igc = 1, ngc2
            send_down(igc,ndown) = gc(igc,n) 
          end do
          gc(7,n) = 0.0d0

        end if
      end do
#endif


      if((idirection.eq.0.and.(my_rank_r.eq.0.or.mpi_proc_mhd.eq.2)).or. &
         (idirection.eq.1.and.(my_rank_z.eq.0.or.mpi_proc_z*mpi_proc_phi.eq.2)) &
        )then !2013-07-17
!2012-06-22
          call mpi_isend(nup,1,mpi_integer,node_up_each,1,each_world &
                        ,my_request(1),mpi_err)
          call mpi_irecv(nup_r,1,mpi_integer,node_up_each,1,each_world &
                        ,my_request(2),mpi_err)
!2012-06-22 end
          call mpi_waitall(2,my_request,my_status_all,mpi_err)

      else if((idirection.eq.0.and.my_rank_r.eq.(mpi_proc_r-1)).or. &
              (idirection.eq.1.and.my_rank_z.eq.(mpi_proc_z-1)) )then
!2012-06-22
          call mpi_isend(ndown,1,mpi_integer,node_down_each,1 &
                        ,each_world &
                        ,my_request(1),mpi_err)
          call mpi_irecv(ndown_r,1,mpi_integer,node_down_each,1 &
                        ,each_world &
                        ,my_request(2),mpi_err)
!2012-06-22 end

          call mpi_waitall(2,my_request,my_status_all,mpi_err)
      else
!2012-06-22
          call mpi_isend(nup,1,mpi_integer,node_up_each,1,each_world &
                        ,my_request(1),mpi_err)
          call mpi_isend(ndown,1,mpi_integer,node_down_each,1 &
                        ,each_world &
                        ,my_request(2),mpi_err)
          call mpi_irecv(nup_r,1,mpi_integer,node_up_each,1,each_world &
                        ,my_request(3),mpi_err)
          call mpi_irecv(ndown_r,1,mpi_integer,node_down_each,1 &
                        ,each_world &
                        ,my_request(4),mpi_err)
!2012-06-22 end

          call mpi_waitall(4,my_request,my_status_all,mpi_err)
      end if


! communication of particle data

!2016-01-05s
      nup8 = nup*ngc2
      ndown8 = ndown*ngc2
      nup_r8 = nup_r*ngc2
      ndown_r8 = ndown_r*ngc2
!2016-01-05e

      if((idirection.eq.0.and.(my_rank_r.eq.0.or.mpi_proc_mhd.eq.2)).or. &
         (idirection.eq.1.and.(my_rank_z.eq.0.or.mpi_proc_z*mpi_proc_phi.eq.2)) &
        )then !2013-07-17
!2012-06-22
        call mpi_isend(send_up(1,1),nup8,mpi_real8,node_up_each,1 &
                      ,each_world,my_request(1),mpi_err)
        call mpi_irecv(recv_up(1,1),nup_r8,mpi_real8,node_up_each,1 &
                      ,each_world,my_request(2),mpi_err)
!2012-06-22 end
        call mpi_waitall(2,my_request,my_status_all,mpi_err)

      else if((idirection.eq.0.and.my_rank_r.eq.(mpi_proc_r-1)).or. &
              (idirection.eq.1.and.my_rank_z.eq.(mpi_proc_z-1)) )then
!2012-06-22
        call mpi_isend(send_down(1,1),ndown8,mpi_real8,node_down_each,1 &
                      ,each_world,my_request(1),mpi_err)
        call mpi_irecv(recv_down(1,1),ndown_r8,mpi_real8,node_down_each,1 &
                      ,each_world,my_request(2),mpi_err)
!2012-06-22 end
        call mpi_waitall(2,my_request,my_status_all,mpi_err)

      else
!2012-06-22
        call mpi_isend(send_up(1,1),nup8,mpi_real8,node_up_each,1 &
                      ,each_world,my_request(1),mpi_err)
        call mpi_isend(send_down(1,1),ndown8,mpi_real8,node_down_each,1 &
                      ,each_world,my_request(2),mpi_err)
        call mpi_irecv(recv_up(1,1),nup_r8,mpi_real8,node_up_each,1 &
                      ,each_world,my_request(3),mpi_err)
        call mpi_irecv(recv_down(1,1),ndown_r8,mpi_real8,node_down_each,1 &
                      ,each_world,my_request(4),mpi_err)
!2012-06-22 end
        call mpi_waitall(4,my_request,my_status_all,mpi_err)

      end if


!----------------------------------------------
! arrange particles received from upward domain

!check      if(my_rank.eq.5)then
!        write(7,*)'kstep, marker_num=',kstep,marker_num
!        write(7,*)'nup, nup_r=',nup,nup_r
!        write(7,*)'ndown, ndown_r=',ndown,ndown_r
!check      end if

      nloop = min(nup, nup_r)

! #ifdef AMDGPU
!!$omp target teams distribute parallel do private(i,n,igc)
! #else      
!$omp parallel do private(n)
! #endif
      do i = 1, nloop
        n = nsend_up(i)
        node(n) = my_rank

        do igc = 1, ngc2
          gc(igc,n) = recv_up(igc,i)
        end do

      end do


!----------------------------------------------
! arrange particles received from downward domain

      nloop = min(ndown, ndown_r)

! #ifdef AMDGPU
!!$omp target teams distribute parallel do private(i,n,igc)
! #else      
!$omp parallel do private(n)
! #endif
      do i = 1, nloop
        n = nsend_down(i)
        node(n) = my_rank

        do igc = 1, ngc2
          gc(igc,n) = recv_down(igc,i)
        end do

      end do

!----------------------------------------------

      if(nup_r.gt.nup)then

#ifdef OPEN_MP
        nn=0
        n = 0
      do i = nup+1, nup_r
 1000 continue
        n = n + 1
        if(n.gt.marker_num)then
          write(7,*)'overflow in particle buffer-up, rank=',my_rank
          go to 9999
        end if
        if(node(n).eq.my_rank)then
          go to 1000
        end if
        nn=nn+1
        ia(nn)=i
        na(nn)=n
      end do
      nnmax=nn

!$omp parallel do private(i,n)
      do nn = 1,nnmax
        i=ia(nn)
        n=na(nn)
        node(n) = my_rank

        do igc = 1, ngc2
          gc(igc,n) = recv_up(igc,i)
        end do

      end do

#else

        n = 0
      do i = nup+1, nup_r
 1000 continue
        n = n + 1
        if(n.gt.marker_num)then
          write(7,*)'overflow in particle buffer-up, rank=',my_rank
          go to 9999
        end if
        if(node(n).eq.my_rank)then
          go to 1000
        end if

        node(n) = my_rank

        do igc = 1, ngc2
          gc(igc,n) = recv_up(igc,i)
        end do

      end do

#endif
      end if


      if(ndown_r.gt.ndown)then
#ifdef OPEN_MP
        nn=0
        n = 0
      do i = ndown+1, ndown_r
 2000 continue
        n = n + 1
        if(n.gt.marker_num)then
          write(7,*)'overflow in particle buffer-down, rank=',my_rank
          go to 9999
        end if
        if(node(n).eq.my_rank)then
          go to 2000
        end if
        nn=nn+1
        ia(nn)=i
        na(nn)=n
      end do
      nnmax=nn

!$omp parallel do private(i,n)
      do nn = 1,nnmax
        i=ia(nn)
        n=na(nn)
        node(n) = my_rank

        do igc = 1, ngc2
          gc(igc,n) = recv_down(igc,i)
        end do

      end do

#else

        n = 0
      do i = ndown+1, ndown_r
 2000 continue
        n = n + 1
        if(n.gt.marker_num)then
          write(7,*)'overflow in particle buffer-down, rank=',my_rank
          go to 9999
        end if
        if(node(n).eq.my_rank)then
          go to 2000
        end if

        node(n) = my_rank

        do igc = 1, ngc2
          gc(igc,n) = recv_down(igc,i)
        end do

      end do
#endif

      end if

 9999 continue


end
!--------------------------------------------------------------------
subroutine partsm1(xxx,cww)
!--------------------------------------------------------------------
      use grid
      implicit none

      real(8)::xxx(lr,lz,lphi),aaa(lr,lz,lphi)
      real(8)::cww,c1,c0
      integer::i,j,k,i1,i2,j1,j2,k1,k2

      c1 = cww/(1.0d0 + 2.0d0*cww)
      c0 = 1.0d0/(1.0d0 + 2.0d0*cww)

#ifdef AMDGPU
!$omp target teams distribute parallel do collapse(3) private(k,j,i)
#else
!$omp parallel do
#endif
      do k = 1, lphi
      do j = 1, lz
      do i = 1, lr
        aaa(i,j,k) = xxx(i,j,k)
      end do
      end do
      end do

#ifdef AMDGPU
!$omp target teams distribute parallel do collapse(3) private(k,j,i,i2,i1)
#else
!$omp parallel do private(i2,i1)
#endif
      do k = 1, lphi
      do j = 1, lz
      do i = 2, lr-1
        i2 = i + 1 
        i1 = i - 1 
        xxx(i,j,k) = c1*aaa(i1,j,k) + c0*aaa(i,j,k) + c1*aaa(i2,j,k) 
      end do
      end do
      end do

#ifdef AMDGPU
!$omp target teams distribute parallel do collapse(3) private(k,j,i,j2,j1)
#else
!$omp parallel do private(j2,j1)
#endif
      do k = 1, lphi
      do j = 2, lz-1
      do i = 1, lr
        j2 = j + 1 
        j1 = j - 1 
        aaa(i,j,k) = c1*xxx(i,j1,k) + c0*xxx(i,j,k) + c1*xxx(i,j2,k) 
      end do
      end do
      end do

#ifdef AMDGPU
!$omp target teams distribute parallel do collapse(3) private(k,j,i,k2,k1)
#else
!$omp parallel do private(k2,k1)
#endif
      do k = 3, lphi-2
      do j = 1, lz
      do i = 1, lr
        k2 = k + 1 
        k1 = k - 1 
        xxx(i,j,k) = c1*aaa(i,j,k1) + c0*aaa(i,j,k) + c1*aaa(i,j,k2) 
      end do
      end do
      end do

! periodic condition

      call periodic_field(xxx)


end
!--------------------------------------------------------------------
subroutine gradient(iflag,a,a_r,a_z,a_phi)
!--------------------------------------------------------------------
      use parameters
      use grid
      implicit none

      real(8)::a(lr,lz,lphi),a_r(lr,lz,lphi),a_z(lr,lz,lphi),a_phi(lr,lz,lphi)
      real(8)::dr1,dz1,dphi1
      integer::i,j,k,i0,i1,i2,i3,k0,k1,k2,k3
      integer::iflag
      real(8)::cc0=4.0d0/3.0d0,cc1=-1.0d0/6.0d0

      dr1= 0.50d0/dr
      dz1= 0.50d0/dz
      dphi1= 0.50d0/dphi
!      cc0= 4.0d0/3.0d0
!      cc1=-1.0d0/6.0d0 


! r direction
#ifdef AMDGPU
!$omp target teams distribute parallel do collapse(3) private(k,j,i,i3,i2,i1,i0)
#else
!$omp parallel do private(i0,i1,i2,i3)
#endif
      do k = 3, lphi-2
      do j = 1, lz
      do i = 3, lr-2
        i3 = i + 2
        i2 = i + 1
        i1 = i - 1
        i0 = i - 2

        a_r(i,j,k) = (cc0*a(i2,j,k) + cc1*a(i3,j,k)        &
                     -cc0*a(i1,j,k) - cc1*a(i0,j,k) )*dr1
      end do
      end do
      end do

#ifdef AMDGPU
!$omp target teams distribute parallel do collapse(2) private(k,j)
#else
!$omp parallel do
#endif
      do k = 3, lphi-2
      do j = 1, lz
        a_r(1,   j,k) = (a(2, j,k)-a(1,   j,k) )*dr1*2.0d0
        a_r(2,   j,k) = (a(3, j,k)-a(1,   j,k) )*dr1
        a_r(lr-1,j,k) = (a(lr,j,k)-a(lr-2,j,k) )*dr1
        a_r(lr,  j,k) = (a(lr,j,k)-a(lr-1,j,k) )*dr1*2.0d0
      end do
      end do

! z direction

#ifdef AMDGPU
!$omp target teams distribute parallel do collapse(2) private(k,i,i0,i1,i2,i3)
#else
!$omp parallel do private(i0,i1,i2,i3)
#endif
      do k = 3, lphi-2
      do i = 2*lr+1, lrz-2*lr
        i3 = i + 2*lr
        i2 = i + lr
        i1 = i - lr
        i0 = i - 2*lr

        a_z(i,1,k) = (a(i2,1,k)*cc0 + a(i3,1,k)*cc1       &
                     -a(i1,1,k)*cc0 - a(i0,1,k)*cc1 )*dz1

!      end do
      end do
      end do

#ifdef AMDGPU  
!$omp target teams distribute parallel do collapse(2) private(k,i)
#else
!$omp parallel do
#endif
      do k = 3, lphi-2
      do i = 1, lr
        a_z(i,1,   k) = (a(i,2, k) - a(i,1,   k) )*dz1*2.0d0
        a_z(i,2,   k) = (a(i,3, k) - a(i,1,   k) )*dz1
        a_z(i,lz-1,k) = (a(i,lz,k) - a(i,lz-2,k) )*dz1
        a_z(i,lz,  k) = (a(i,lz,k) - a(i,lz-1,k) )*dz1*2.0d0
      end do
      end do

! phi direction

#ifdef AMDGPU  
!$omp target teams distribute parallel do collapse(2) private(k,i,k3,k2,k1,k0)
#else
!$omp parallel do private(k0,k1,k2,k3)
#endif
      do k = 3, lphi-2
      do i = 1, lrz
        k3 = k + 2
        k2 = k + 1
        k1 = k - 1
        k0 = k - 2

        a_phi(i,1,k) = (a(i,1,k2)*cc0 + a(i,1,k3)*cc1  &
                       -a(i,1,k1)*cc0 - a(i,1,k0)*cc1  &
                       )*dphi1/grr(i,1,k)

      end do
      end do

      if(iflag.eq.1)then
        call periodic_field_mlt3(a_r,a_z,a_phi) !2012-04-13
!        call periodic_field(a_r)
!        call periodic_field(a_z)
!        call periodic_field(a_phi)
      end if

      end
!--------------------------------------------------------------------
subroutine gradient5(iflag,cf,sc,a,a_r,a_z,a_phi)
!--------------------------------------------------------------------
      use parameters
      use grid
      implicit none

      real(8)::cf(lr,lz,lphi),sc(lr,lz,lphi)
      real(8)::a(lr,lz,lphi),a_r(lr,lz,lphi),a_z(lr,lz,lphi),a_phi(lr,lz,lphi)
      real(8)::dr1,dz1,dphi1
      integer::i,j,k,i0,i1,i2,i3,k0,k1,k2,k3
      integer::iflag
      real(8)::cc0=4.0d0/3.0d0,cc1=-1.0d0/6.0d0

      dr1= 0.50d0/dr
      dz1= 0.50d0/dz
      dphi1= 0.50d0/dphi
!      cc0= 4.0d0/3.0d0
!      cc1=-1.0d0/6.0d0 


! r direction

#ifdef AMDGPU  
!$omp target teams distribute parallel do collapse(3) private(k,j,i,i3,i2,i1,i0)
#else
!$omp parallel do private(i0,i1,i2,i3)
#endif
      do k = 3, lphi-2
      do j = 1, lz
      do i = 3, lr-2
        i3 = i + 2
        i2 = i + 1
        i1 = i - 1
        i0 = i - 2

        a_r(i,j,k) = (cc0*cf(i2,j,k)*sc(i2,j,k)*a(i2,j,k) &
                    + cc1*cf(i3,j,k)*sc(i3,j,k)*a(i3,j,k) &
                     -cc0*cf(i1,j,k)*sc(i1,j,k)*a(i1,j,k) &
                     -cc1*cf(i0,j,k)*sc(i0,j,k)*a(i0,j,k) )*dr1
      end do
      end do
      end do

#ifdef AMDGPU  
!$omp target teams distribute parallel do collapse(2) private(k,j)
#else
!$omp parallel do
#endif
      do k = 3, lphi-2
      do j = 1, lz
        a_r(1,   j,k) = (cf(2,j,k)*sc(2,j,k)*a(2, j,k)-cf(1,j,k)*sc(1,j,k)*a(1,   j,k) )*dr1*2.0d0
        a_r(2,   j,k) = (cf(3,j,k)*sc(3,j,k)*a(3, j,k)-cf(1,j,k)*sc(1,j,k)*a(1,   j,k) )*dr1
        a_r(lr-1,j,k) = (cf(lr,j,k)*sc(lr,j,k)*a(lr,j,k)-cf(lr-2,j,k)*sc(lr-2,j,k)*a(lr-2,j,k) )*dr1
        a_r(lr,  j,k) = (cf(lr,j,k)*sc(lr,j,k)*a(lr,j,k)-cf(lr-1,j,k)*sc(lr-1,j,k)*a(lr-1,j,k) )*dr1*2.0d0
      end do
      end do

! z direction

#ifdef AMDGPU  
!$omp target teams distribute parallel do collapse(2) private(k,i,i3,i2,i1,i0)
#else
!$omp parallel do private(i0,i1,i2,i3)
#endif
      do k = 3, lphi-2
      do i = 2*lr+1, lrz-2*lr
        i3 = i + 2*lr
        i2 = i + lr
        i1 = i - lr
        i0 = i - 2*lr

        a_z(i,1,k) = (cf(i2,1,k)*sc(i2,1,k)*a(i2,1,k)*cc0 &
                    + cf(i3,1,k)*sc(i3,1,k)*a(i3,1,k)*cc1 &
                     -cf(i1,1,k)*sc(i1,1,k)*a(i1,1,k)*cc0 &
                     -cf(i0,1,k)*sc(i0,1,k)*a(i0,1,k)*cc1 )*dz1
      end do
      end do

#ifdef AMDGPU  
!$omp target teams distribute parallel do collapse(2) private(k,i)
#else
!$omp parallel do
#endif
      do k = 3, lphi-2
      do i = 1, lr
        a_z(i,1,   k) = (cf(i,2,k)*sc(i,2,k)*a(i,2, k) - cf(i,1,k)*sc(i,1,k)*a(i,1,   k) )*dz1*2.0d0
        a_z(i,2,   k) = (cf(i,3,k)*sc(i,3,k)*a(i,3, k) - cf(i,1,k)*sc(i,1,k)*a(i,1,   k) )*dz1
        a_z(i,lz-1,k) = (cf(i,lz,k)*sc(i,lz,k)*a(i,lz,k) - cf(i,lz-2,k)*sc(i,lz-2,k)*a(i,lz-2,k) )*dz1
        a_z(i,lz,  k) = (cf(i,lz,k)*sc(i,lz,k)*a(i,lz,k) - cf(i,lz-1,k)*sc(i,lz-1,k)*a(i,lz-1,k) )*dz1*2.0d0
      end do
      end do

! phi direction

#ifdef AMDGPU  
!$omp target teams distribute parallel do collapse(2) private(k,i,k0,k1,k2,k3)
#else
!$omp parallel do private(k0,k1,k2,k3)
#endif
      do k = 3, lphi-2
      do i = 1, lrz
        k3 = k + 2
        k2 = k + 1
        k1 = k - 1
        k0 = k - 2

        a_phi(i,1,k) = (cf(i,1,k2)*sc(i,1,k2)*a(i,1,k2)*cc0 &
                      + cf(i,1,k3)*sc(i,1,k3)*a(i,1,k3)*cc1 &
                       -cf(i,1,k1)*sc(i,1,k1)*a(i,1,k1)*cc0 &
                       -cf(i,1,k0)*sc(i,1,k0)*a(i,1,k0)*cc1  &
                       )*dphi1/grr(i,1,k)

      end do
      end do

      if(iflag.eq.1)then
        call periodic_field_mlt3(a_r,a_z,a_phi) !2012-04-13
!        call periodic_field(a_r)
!        call periodic_field(a_z)
!        call periodic_field(a_phi)
      end if

      end
!--------------------------------------------------------------------
subroutine divergence(iflag,ar,az,aphi,div)
!--------------------------------------------------------------------
      use parameters
      use grid
      implicit none

      real(8)::ar(lr,lz,lphi),az(lr,lz,lphi),aphi(lr,lz,lphi),div(lr,lz,lphi)
      real(8)::dr1,dz1,dphi1
      integer::i,j,k,i0,i1,i2,i3,k0,k1,k2,k3
      integer::iflag
      real(8)::cc0=4.0d0/3.0d0,cc1=-1.0d0/6.0d0,cc2=7.0d0/6.0d0

      dr1= 0.50d0/dr
      dz1= 0.50d0/dz
      dphi1= 0.50d0/dphi
!      cc0= 4.0d0/3.0d0
!      cc1=-1.0d0/6.0d0 
!      cc2= 7.0d0/6.0d0 

! r direction

#ifdef AMDGPU
!$omp target teams distribute parallel do collapse(3) private(k,j,i,i3,i2,i1,i0)
#else
!$omp parallel do private(i0,i1,i2,i3)
#endif
      do k = 3, lphi-2
      do j = 1, lz
      do i = 3, lr-2
        i3 = i + 2
        i2 = i + 1
        i1 = i - 1
        i0 = i - 2

        div(i,j,k) = (cc0*ar(i2,j,k)*grr(i2,j,k) &
                     +cc1*ar(i3,j,k)*grr(i3,j,k) &
                     -cc0*ar(i1,j,k)*grr(i1,j,k) &
                     -cc1*ar(i0,j,k)*grr(i0,j,k) &
                     )*dr1/grr(i,j,k)
      end do
      end do
      end do

#ifdef AMDGPU
!$omp target teams distribute parallel do collapse(2) private(k,j)
#else
!$omp parallel do
#endif
      do k = 3, lphi-2
      do j = 1, lz
        div(1,j,k) = ( ar(3,j,k)*grr(3,j,k)*cc1 &
                      +ar(2,j,k)*grr(2,j,k)*cc0 &
                      -ar(1,j,k)*grr(1,j,k)*cc2 &
                     )*dr1/grr(1,j,k)
        div(2,j,k) = (ar(4,j,k)*grr(4,j,k)*cc1 &
                     +ar(3,j,k)*grr(3,j,k)*cc0 &
                     +ar(2,j,k)*grr(2,j,k)*cc1 &
                     -ar(1,j,k)*grr(1,j,k)     &
                     )*dr1/grr(2,j,k)
        div(lr-1,j,k) = (ar(lr,  j,k)*grr(lr,  j,k)     &
                        -ar(lr-1,j,k)*grr(lr-1,j,k)*cc1 &
                        -ar(lr-2,j,k)*grr(lr-2,j,k)*cc0 &
                        -ar(lr-3,j,k)*grr(lr-3,j,k)*cc1 &
                        )*dr1/grr(lr-1,j,k)
        div(lr,  j,k) = (ar(lr,  j,k)*grr(lr,  j,k)*cc2 &
                        -ar(lr-1,j,k)*grr(lr-1,j,k)*cc0 &
                        -ar(lr-2,j,k)*grr(lr-2,j,k)*cc1 &
                        )*dr1/grr(lr,j,k)
      end do
      end do

! z direction

#ifdef AMDGPU
!$omp target teams distribute parallel do collapse(2) private(k,i,i3,i2,i1,i0)
#else
!$omp parallel do private(i0,i1,i2,i3)
#endif
      do k = 3, lphi-2
      do i = 2*lr+1, lrz-2*lr

        i3 = i + 2*lr
        i2 = i + lr
        i1 = i - lr
        i0 = i - 2*lr
        div(i,1,k) = div(i,1,k) &
                   + (az(i2,1,k)*cc0 + az(i3,1,k)*cc1 &
                     -az(i1,1,k)*cc0 - az(i0,1,k)*cc1 )*dz1

      end do
      end do

#ifdef AMDGPU
!$omp target teams distribute parallel do collapse(2) private(k,i)
#else
!$omp parallel do
#endif
      do k = 3, lphi-2
      do i = 1, lr
        div(i,1,   k) = div(i,1,   k) &
              + (az(i,3,k)*cc1 + az(i,2,k)*cc0 - az(i,1,k)*cc2 )*dz1
        div(i,2,   k) = div(i,2,   k) &
              + (az(i,4,k)*cc1 + az(i,3,k)*cc0 + az(i,2,k)*cc1 -az(i,1,k))*dz1
        div(i,lz-1,k) = div(i,lz-1,k) &
              + (az(i,lz,k)  &
               - az(i,lz-1,k)*cc1 - az(i,lz-2,k)*cc0 - az(i,lz-3,k)*cc1 )*dz1
        div(i,lz,  k) = div(i,lz,  k) &
              + (az(i,lz,k)*cc2 - az(i,lz-1,k)*cc0 - az(i,lz-2,k)*cc1 )*dz1
      end do
      end do

! phi direction

#ifdef AMDGPU
!$omp target teams distribute parallel do collapse(2) private(k,i,k0,k1,k2,k3)
#else
!$omp parallel do private(k0,k1,k2,k3)
#endif
      do k = 3, lphi-2
      do i = 1, lrz
        k3 = k + 2
        k2 = k + 1
        k1 = k - 1
        k0 = k - 2

        div(i,1,k) = div(i,1,k) &
                   + (aphi(i,1,k2)*cc0 + aphi(i,1,k3)*cc1  &
                     -aphi(i,1,k1)*cc0 - aphi(i,1,k0)*cc1  &
                     )*dphi1/grr(i,1,k)
      end do
      end do

      if(iflag.eq.1)then
        call periodic_field(div)
      end if

      end
!--------------------------------------------------------------------
subroutine divergence4
! 2012-04-11 for dvr,dvz,dvphi,dprs in subroutine mhd
!--------------------------------------------------------------------
      use parameters
      use grid
      use field, only:vr,vz,vphi,prs,ur,uz,uphi,div_u,dvr,dvz,dvphi,dprs,gamma
      implicit none

      real(8)::dr1,dz1,dphi1
      real(8)::cc0=4.0d0/3.0d0,cc1=-1.0d0/6.0d0,cc2=7.0d0/6.0d0
      integer::i,j,k,i0,i1,i2,i3,k0,k1,k2,k3

      dr1= 0.50d0/dr
      dz1= 0.50d0/dz
      dphi1= 0.50d0/dphi
!      cc0= 4.0d0/3.0d0
!      cc1=-1.0d0/6.0d0 
!      cc2= 7.0d0/6.0d0 

! r direction

#ifdef AMDGPU
!$omp target teams distribute parallel do collapse(3) private(k,j,i,i0,i1,i2,i3)
#else
!$omp parallel do private(i0,i1,i2,i3)
#endif
      do k = 3, lphi-2
      do j = 1, lz
      do i = 3, lr-2
        i3 = i + 2
        i2 = i + 1
        i1 = i - 1
        i0 = i - 2

        dvr(i,j,k) = (cc0*vr(i2,j,k)*ur(i2,j,k)*grr(i2,j,k) &
                     +cc1*vr(i3,j,k)*ur(i3,j,k)*grr(i3,j,k) &
                     -cc0*vr(i1,j,k)*ur(i1,j,k)*grr(i1,j,k) &
                     -cc1*vr(i0,j,k)*ur(i0,j,k)*grr(i0,j,k) &
                     )*dr1/grr(i,j,k) &
                   + (cc0*prs(i2,j,k) &
                     +cc1*prs(i3,j,k) &
                     -cc0*prs(i1,j,k) &
                     -cc1*prs(i0,j,k) &
                     )*dr1
        dvz(i,j,k) = (cc0*vz(i2,j,k)*ur(i2,j,k)*grr(i2,j,k) &
                     +cc1*vz(i3,j,k)*ur(i3,j,k)*grr(i3,j,k) &
                     -cc0*vz(i1,j,k)*ur(i1,j,k)*grr(i1,j,k) &
                     -cc1*vz(i0,j,k)*ur(i0,j,k)*grr(i0,j,k) &
                     )*dr1/grr(i,j,k)
        dvphi(i,j,k) = (cc0*vphi(i2,j,k)*ur(i2,j,k)*grr(i2,j,k) &
                       +cc1*vphi(i3,j,k)*ur(i3,j,k)*grr(i3,j,k) &
                       -cc0*vphi(i1,j,k)*ur(i1,j,k)*grr(i1,j,k) &
                       -cc1*vphi(i0,j,k)*ur(i0,j,k)*grr(i0,j,k) &
                       )*dr1/grr(i,j,k)
        dprs(i,j,k) = (cc0*prs(i2,j,k)*ur(i2,j,k)*grr(i2,j,k) &
                      +cc1*prs(i3,j,k)*ur(i3,j,k)*grr(i3,j,k) &
                      -cc0*prs(i1,j,k)*ur(i1,j,k)*grr(i1,j,k) &
                      -cc1*prs(i0,j,k)*ur(i0,j,k)*grr(i0,j,k) &
                      )*dr1/grr(i,j,k)
      end do
      end do
      end do

#ifdef AMDGPU
!$omp target teams distribute parallel do collapse(2) private(k,j)
#else
!$omp parallel do
#endif
      do k = 3, lphi-2
      do j = 1, lz
        dvr(1,j,k) = ( vr(3,j,k)*ur(3,j,k)*grr(3,j,k)*cc1 &
                      +vr(2,j,k)*ur(2,j,k)*grr(2,j,k)*cc0 &
                      -vr(1,j,k)*ur(1,j,k)*grr(1,j,k)*cc2 &
                     )*dr1/grr(1,j,k) &
                   + (cc1*prs(3,j,k) &
                     +cc0*prs(2,j,k) &
                     -cc2*prs(1,j,k) &
                     )*dr1
        dvz(1,j,k) = ( vz(3,j,k)*ur(3,j,k)*grr(3,j,k)*cc1 &
                      +vz(2,j,k)*ur(2,j,k)*grr(2,j,k)*cc0 &
                      -vz(1,j,k)*ur(1,j,k)*grr(1,j,k)*cc2 &
                     )*dr1/grr(1,j,k)
        dvphi(1,j,k) = ( vphi(3,j,k)*ur(3,j,k)*grr(3,j,k)*cc1 &
                        +vphi(2,j,k)*ur(2,j,k)*grr(2,j,k)*cc0 &
                        -vphi(1,j,k)*ur(1,j,k)*grr(1,j,k)*cc2 &
                       )*dr1/grr(1,j,k)
        dprs(1,j,k) = ( prs(3,j,k)*ur(3,j,k)*grr(3,j,k)*cc1 &
                       +prs(2,j,k)*ur(2,j,k)*grr(2,j,k)*cc0 &
                       -prs(1,j,k)*ur(1,j,k)*grr(1,j,k)*cc2 &
                      )*dr1/grr(1,j,k)

        dvr(2,j,k) = (vr(4,j,k)*ur(4,j,k)*grr(4,j,k)*cc1 &
                     +vr(3,j,k)*ur(3,j,k)*grr(3,j,k)*cc0 &
                     +vr(2,j,k)*ur(2,j,k)*grr(2,j,k)*cc1 &
                     -vr(1,j,k)*ur(1,j,k)*grr(1,j,k)     &
                     )*dr1/grr(2,j,k) &
                   + (cc1*prs(4,j,k) &
                     +cc0*prs(3,j,k) &
                     +cc1*prs(2,j,k) &
                     -    prs(1,j,k) &
                     )*dr1
        dvz(2,j,k) = (vz(4,j,k)*ur(4,j,k)*grr(4,j,k)*cc1 &
                     +vz(3,j,k)*ur(3,j,k)*grr(3,j,k)*cc0 &
                     +vz(2,j,k)*ur(2,j,k)*grr(2,j,k)*cc1 &
                     -vz(1,j,k)*ur(1,j,k)*grr(1,j,k)     &
                     )*dr1/grr(2,j,k)
        dvphi(2,j,k) = (vphi(4,j,k)*ur(4,j,k)*grr(4,j,k)*cc1 &
                       +vphi(3,j,k)*ur(3,j,k)*grr(3,j,k)*cc0 &
                       +vphi(2,j,k)*ur(2,j,k)*grr(2,j,k)*cc1 &
                       -vphi(1,j,k)*ur(1,j,k)*grr(1,j,k)     &
                       )*dr1/grr(2,j,k)
        dprs(2,j,k) = (prs(4,j,k)*ur(4,j,k)*grr(4,j,k)*cc1 &
                      +prs(3,j,k)*ur(3,j,k)*grr(3,j,k)*cc0 &
                      +prs(2,j,k)*ur(2,j,k)*grr(2,j,k)*cc1 &
                      -prs(1,j,k)*ur(1,j,k)*grr(1,j,k)     &
                      )*dr1/grr(2,j,k)

        dvr(lr-1,j,k) = (vr(lr,  j,k)*ur(lr,  j,k)*grr(lr,  j,k)     &
                        -vr(lr-1,j,k)*ur(lr-1,j,k)*grr(lr-1,j,k)*cc1 &
                        -vr(lr-2,j,k)*ur(lr-2,j,k)*grr(lr-2,j,k)*cc0 &
                        -vr(lr-3,j,k)*ur(lr-3,j,k)*grr(lr-3,j,k)*cc1 &
                        )*dr1/grr(lr-1,j,k) &
                   + (    prs(lr,j,k) &
                     -cc1*prs(lr-1,j,k) &
                     -cc0*prs(lr-2,j,k) &
                     -cc1*prs(lr-3,j,k) &
                     )*dr1
        dvz(lr-1,j,k) = (vz(lr,  j,k)*ur(lr,  j,k)*grr(lr,  j,k)     &
                        -vz(lr-1,j,k)*ur(lr-1,j,k)*grr(lr-1,j,k)*cc1 &
                        -vz(lr-2,j,k)*ur(lr-2,j,k)*grr(lr-2,j,k)*cc0 &
                        -vz(lr-3,j,k)*ur(lr-3,j,k)*grr(lr-3,j,k)*cc1 &
                        )*dr1/grr(lr-1,j,k)
        dvphi(lr-1,j,k) = (vphi(lr,  j,k)*ur(lr,  j,k)*grr(lr,  j,k)     &
                          -vphi(lr-1,j,k)*ur(lr-1,j,k)*grr(lr-1,j,k)*cc1 &
                          -vphi(lr-2,j,k)*ur(lr-2,j,k)*grr(lr-2,j,k)*cc0 &
                          -vphi(lr-3,j,k)*ur(lr-3,j,k)*grr(lr-3,j,k)*cc1 &
                          )*dr1/grr(lr-1,j,k)
        dprs(lr-1,j,k) = (prs(lr,  j,k)*ur(lr,  j,k)*grr(lr,  j,k)     &
                         -prs(lr-1,j,k)*ur(lr-1,j,k)*grr(lr-1,j,k)*cc1 &
                         -prs(lr-2,j,k)*ur(lr-2,j,k)*grr(lr-2,j,k)*cc0 &
                         -prs(lr-3,j,k)*ur(lr-3,j,k)*grr(lr-3,j,k)*cc1 &
                         )*dr1/grr(lr-1,j,k)

        dvr(lr,  j,k) = (vr(lr,  j,k)*ur(lr,  j,k)*grr(lr,  j,k)*cc2 &
                        -vr(lr-1,j,k)*ur(lr-1,j,k)*grr(lr-1,j,k)*cc0 &
                        -vr(lr-2,j,k)*ur(lr-2,j,k)*grr(lr-2,j,k)*cc1 &
                        )*dr1/grr(lr,j,k) &
                   + (cc2*prs(lr,j,k) &
                     -cc0*prs(lr-1,j,k) &
                     -cc1*prs(lr-2,j,k) &
                     )*dr1
        dvz(lr,  j,k) = (vz(lr,  j,k)*ur(lr,  j,k)*grr(lr,  j,k)*cc2 &
                        -vz(lr-1,j,k)*ur(lr-1,j,k)*grr(lr-1,j,k)*cc0 &
                        -vz(lr-2,j,k)*ur(lr-2,j,k)*grr(lr-2,j,k)*cc1 &
                        )*dr1/grr(lr,j,k)
        dvphi(lr,  j,k) = (vphi(lr,  j,k)*ur(lr,  j,k)*grr(lr,  j,k)*cc2 &
                          -vphi(lr-1,j,k)*ur(lr-1,j,k)*grr(lr-1,j,k)*cc0 &
                          -vphi(lr-2,j,k)*ur(lr-2,j,k)*grr(lr-2,j,k)*cc1 &
                          )*dr1/grr(lr,j,k)
        dprs(lr,  j,k) = (prs(lr,  j,k)*ur(lr,  j,k)*grr(lr,  j,k)*cc2 &
                         -prs(lr-1,j,k)*ur(lr-1,j,k)*grr(lr-1,j,k)*cc0 &
                         -prs(lr-2,j,k)*ur(lr-2,j,k)*grr(lr-2,j,k)*cc1 &
                         )*dr1/grr(lr,j,k)
      end do
      end do

! z direction

#ifdef AMDGPU
!$omp target teams distribute parallel do collapse(2) private(k,i,i0,i1,i2,i3)
#else
!$omp parallel do private(i0,i1,i2,i3)
#endif
      do k = 3, lphi-2
      do i = 2*lr+1, lrz-2*lr

        i3 = i + 2*lr
        i2 = i + lr
        i1 = i - lr
        i0 = i - 2*lr
        dvr(i,1,k) = dvr(i,1,k) &
                   + (vr(i2,1,k)*uz(i2,1,k)*cc0 + vr(i3,1,k)*uz(i3,1,k)*cc1 &
                     -vr(i1,1,k)*uz(i1,1,k)*cc0 - vr(i0,1,k)*uz(i0,1,k)*cc1 )*dz1
        dvz(i,1,k) = dvz(i,1,k) &
                   + (vz(i2,1,k)*uz(i2,1,k)*cc0 + vz(i3,1,k)*uz(i3,1,k)*cc1 &
                     -vz(i1,1,k)*uz(i1,1,k)*cc0 - vz(i0,1,k)*uz(i0,1,k)*cc1 )*dz1 &
                   + (prs(i2,1,k)*cc0 + prs(i3,1,k)*cc1 &
                     -prs(i1,1,k)*cc0 - prs(i0,1,k)*cc1 )*dz1
        dvphi(i,1,k) = dvphi(i,1,k) &
                   + (vphi(i2,1,k)*uz(i2,1,k)*cc0 + vphi(i3,1,k)*uz(i3,1,k)*cc1 &
                     -vphi(i1,1,k)*uz(i1,1,k)*cc0 - vphi(i0,1,k)*uz(i0,1,k)*cc1 )*dz1
        dprs(i,1,k) = dprs(i,1,k) &
                   + (prs(i2,1,k)*uz(i2,1,k)*cc0 + prs(i3,1,k)*uz(i3,1,k)*cc1 &
                     -prs(i1,1,k)*uz(i1,1,k)*cc0 - prs(i0,1,k)*uz(i0,1,k)*cc1 )*dz1

      end do
      end do

#ifdef AMDGPU
!$omp target teams distribute parallel do collapse(2) private(k,i)
#else
!$omp parallel do
#endif
      do k = 3, lphi-2
      do i = 1, lr
        dvr(i,1,   k) = dvr(i,1,   k) &
              + (vr(i,3,k)*uz(i,3,k)*cc1 + vr(i,2,k)*uz(i,2,k)*cc0 - vr(i,1,1)*uz(i,1,k)*cc2 )*dz1
        dvz(i,1,   k) = dvz(i,1,   k) &
              + (vz(i,3,k)*uz(i,3,k)*cc1 + vz(i,2,k)*uz(i,2,k)*cc0 - vz(i,1,1)*uz(i,1,k)*cc2 )*dz1 &
              + (prs(i,3,k)*cc1 + prs(i,2,k)*cc0 - prs(i,1,1)*cc2 )*dz1
        dvphi(i,1,   k) = dvphi(i,1,   k) &
              + (vphi(i,3,k)*uz(i,3,k)*cc1 + vphi(i,2,k)*uz(i,2,k)*cc0 - vphi(i,1,1)*uz(i,1,k)*cc2 )*dz1
        dprs(i,1,   k) = dprs(i,1,   k) &
              + (prs(i,3,k)*uz(i,3,k)*cc1 + prs(i,2,k)*uz(i,2,k)*cc0 - prs(i,1,1)*uz(i,1,k)*cc2 )*dz1


        dvr(i,2,   k) = dvr(i,2,   k) &
              + (vr(i,4,k)*uz(i,4,k)*cc1 + vr(i,3,k)*uz(i,3,k)*cc0 &
                +vr(i,2,k)*uz(i,2,k)*cc1 - vr(i,1,k)*uz(i,1,k))*dz1
        dvz(i,2,   k) = dvz(i,2,   k) &
              + (vz(i,4,k)*uz(i,4,k)*cc1 + vz(i,3,k)*uz(i,3,k)*cc0 &
                +vz(i,2,k)*uz(i,2,k)*cc1 - vz(i,1,k)*uz(i,1,k))*dz1 &
              + (prs(i,4,k)*cc1 + prs(i,3,k)*cc0 &
                +prs(i,2,k)*cc1 - prs(i,1,k))*dz1
        dvphi(i,2,   k) = dvphi(i,2,   k) &
              + (vphi(i,4,k)*uz(i,4,k)*cc1 + vphi(i,3,k)*uz(i,3,k)*cc0 &
                +vphi(i,2,k)*uz(i,2,k)*cc1 - vphi(i,1,k)*uz(i,1,k))*dz1
        dprs(i,2,   k) = dprs(i,2,   k) &
              + (prs(i,4,k)*uz(i,4,k)*cc1 + prs(i,3,k)*uz(i,3,k)*cc0 &
                +prs(i,2,k)*uz(i,2,k)*cc1 - prs(i,1,k)*uz(i,1,k))*dz1

        dvr(i,lz-1,k) = dvr(i,lz-1,k) &
              + (vr(i,lz,k)*uz(i,lz,k) - vr(i,lz-1,k)*uz(i,lz-1,k)*cc1 &
                -vr(i,lz-2,k)*uz(i,lz-2,k)*cc0 - vr(i,lz-3,k)*uz(i,lz-3,k)*cc1 )*dz1
        dvz(i,lz-1,k) = dvz(i,lz-1,k) &
              + (vz(i,lz,k)*uz(i,lz,k) - vz(i,lz-1,k)*uz(i,lz-1,k)*cc1 &
                -vz(i,lz-2,k)*uz(i,lz-2,k)*cc0 - vz(i,lz-3,k)*uz(i,lz-3,k)*cc1 )*dz1 &
              + (prs(i,lz,k) - prs(i,lz-1,k)*cc1 &
                -prs(i,lz-2,k)*cc0 - prs(i,lz-3,k)*cc1 )*dz1
        dvphi(i,lz-1,k) = dvphi(i,lz-1,k) &
              + (vphi(i,lz,k)*uz(i,lz,k) - vphi(i,lz-1,k)*uz(i,lz-1,k)*cc1 &
                -vphi(i,lz-2,k)*uz(i,lz-2,k)*cc0 - vphi(i,lz-3,k)*uz(i,lz-3,k)*cc1 )*dz1
        dprs(i,lz-1,k) = dprs(i,lz-1,k) &
              + (prs(i,lz,k)*uz(i,lz,k) - prs(i,lz-1,k)*uz(i,lz-1,k)*cc1 &
                -prs(i,lz-2,k)*uz(i,lz-2,k)*cc0 - prs(i,lz-3,k)*uz(i,lz-3,k)*cc1 )*dz1

        dvr(i,lz,  k) = dvr(i,lz,  k) &
              + (vr(i,lz,k)*uz(i,lz,k)*cc2 - vr(i,lz-1,k)*uz(i,lz-1,k)*cc0 - vr(i,lz-2,k)*uz(i,lz-2,k)*cc1 )*dz1
        dvz(i,lz,  k) = dvz(i,lz,  k) &
              + (vz(i,lz,k)*uz(i,lz,k)*cc2 - vz(i,lz-1,k)*uz(i,lz-1,k)*cc0 - vz(i,lz-2,k)*uz(i,lz-2,k)*cc1 )*dz1 &
              + (prs(i,lz,k)*cc2 - prs(i,lz-1,k)*cc0 - prs(i,lz-2,k)*cc1 )*dz1
        dvphi(i,lz,  k) = dvphi(i,lz,  k) &
              + (vphi(i,lz,k)*uz(i,lz,k)*cc2 - vphi(i,lz-1,k)*uz(i,lz-1,k)*cc0 - vphi(i,lz-2,k)*uz(i,lz-2,k)*cc1 )*dz1
        dprs(i,lz,  k) = dprs(i,lz,  k) &
              + (prs(i,lz,k)*uz(i,lz,k)*cc2 - prs(i,lz-1,k)*uz(i,lz-1,k)*cc0 - prs(i,lz-2,k)*uz(i,lz-2,k)*cc1 )*dz1

      end do
      end do

! phi direction

#ifdef AMDGPU
!$omp target teams distribute parallel do collapse(2) private(k,i,k0,k1,k2,k3)
#else
!$omp parallel do private(k0,k1,k2,k3)
#endif
      do k = 3, lphi-2
      do i = 1, lrz
        k3 = k + 2
        k2 = k + 1
        k1 = k - 1
        k0 = k - 2

        dvr(i,1,k) = dvr(i,1,k) &
                   + (vr(i,1,k2)*uphi(i,1,k2)*cc0 + vr(i,1,k3)*uphi(i,1,k3)*cc1  &
                     -vr(i,1,k1)*uphi(i,1,k1)*cc0 - vr(i,1,k0)*uphi(i,1,k0)*cc1  &
                     )*dphi1/grr(i,1,k)
        dvz(i,1,k) = dvz(i,1,k) &
                   + (vz(i,1,k2)*uphi(i,1,k2)*cc0 + vz(i,1,k3)*uphi(i,1,k3)*cc1  &
                     -vz(i,1,k1)*uphi(i,1,k1)*cc0 - vz(i,1,k0)*uphi(i,1,k0)*cc1  &
                     )*dphi1/grr(i,1,k)
        dvphi(i,1,k) = dvphi(i,1,k) &
                   + (vphi(i,1,k2)*uphi(i,1,k2)*cc0 + vphi(i,1,k3)*uphi(i,1,k3)*cc1  &
                     -vphi(i,1,k1)*uphi(i,1,k1)*cc0 - vphi(i,1,k0)*uphi(i,1,k0)*cc1  &
                   +  prs(i,1,k2)*cc0 + prs(i,1,k3)*cc1  &
                     -prs(i,1,k1)*cc0 - prs(i,1,k0)*cc1  &
                     )*dphi1/grr(i,1,k)
        dprs(i,1,k) = dprs(i,1,k) &
                   + (prs(i,1,k2)*uphi(i,1,k2)*cc0 + prs(i,1,k3)*uphi(i,1,k3)*cc1  &
                     -prs(i,1,k1)*uphi(i,1,k1)*cc0 - prs(i,1,k0)*uphi(i,1,k0)*cc1  &
                     )*dphi1/grr(i,1,k) &
                    + (gamma-1.0d0)*prs(i,1,k)*div_u(i,1,k)
      end do
      end do


      end
!--------------------------------------------------------------------
subroutine divergence3hm
! 2012-08-16 for dvr,dvz,dvphi,dprs in subroutine hm
! 2014-01-10 toroidal rotation is introduced, in utotphi
! 2017-05-04 for dvr,dvz,dvphi in subroutine hm
!--------------------------------------------------------------------
      use parameters
      use grid
      use field, only:vr,vz,vphi,prs,vdrftr,vdrftz,vdrftphi &
                     ,dvr,dvz,dvphi
      implicit none

      real(8)::dr1,dz1,dphi1
      real(8)::cc0=4.0d0/3.0d0,cc1=-1.0d0/6.0d0,cc2=7.0d0/6.0d0
      integer::i,j,k,i0,i1,i2,i3,k0,k1,k2,k3

      dr1= 0.50d0/dr
      dz1= 0.50d0/dz
      dphi1= 0.50d0/dphi
!      cc0= 4.0d0/3.0d0
!      cc1=-1.0d0/6.0d0 
!      cc2= 7.0d0/6.0d0 

! r direction

#ifdef AMDGPU
!$omp target teams distribute parallel do collapse(3) private(k,j,i,i0,i1,i2,i3)
#else
!$omp parallel do private(i0,i1,i2,i3)
#endif
      do k = 3, lphi-2
      do j = 1, lz
      do i = 3, lr-2
        i3 = i + 2
        i2 = i + 1
        i1 = i - 1
        i0 = i - 2

        dvr(i,j,k) = (cc0*vr(i2,j,k)*vdrftr(i2,j,k)*grr(i2,j,k) &
                     +cc1*vr(i3,j,k)*vdrftr(i3,j,k)*grr(i3,j,k) &
                     -cc0*vr(i1,j,k)*vdrftr(i1,j,k)*grr(i1,j,k) &
                     -cc1*vr(i0,j,k)*vdrftr(i0,j,k)*grr(i0,j,k) &
                     )*dr1/grr(i,j,k) &
                   + (cc0*prs(i2,j,k) &
                     +cc1*prs(i3,j,k) &
                     -cc0*prs(i1,j,k) &
                     -cc1*prs(i0,j,k) &
                     )*dr1
        dvz(i,j,k) = (cc0*vz(i2,j,k)*vdrftr(i2,j,k)*grr(i2,j,k) &
                     +cc1*vz(i3,j,k)*vdrftr(i3,j,k)*grr(i3,j,k) &
                     -cc0*vz(i1,j,k)*vdrftr(i1,j,k)*grr(i1,j,k) &
                     -cc1*vz(i0,j,k)*vdrftr(i0,j,k)*grr(i0,j,k) &
                     )*dr1/grr(i,j,k)
        dvphi(i,j,k) = (cc0*vphi(i2,j,k)*vdrftr(i2,j,k)*grr(i2,j,k) &
                       +cc1*vphi(i3,j,k)*vdrftr(i3,j,k)*grr(i3,j,k) &
                       -cc0*vphi(i1,j,k)*vdrftr(i1,j,k)*grr(i1,j,k) &
                       -cc1*vphi(i0,j,k)*vdrftr(i0,j,k)*grr(i0,j,k) &
                       )*dr1/grr(i,j,k)
      end do
      end do
      end do

#ifdef AMDGPU
!$omp target teams distribute parallel do collapse(2) private(k,j)
#else
!$omp parallel do
#endif
      do k = 3, lphi-2
      do j = 1, lz
        dvr(1,j,k) = ( vr(3,j,k)*vdrftr(3,j,k)*grr(3,j,k)*cc1 &
                      +vr(2,j,k)*vdrftr(2,j,k)*grr(2,j,k)*cc0 &
                      -vr(1,j,k)*vdrftr(1,j,k)*grr(1,j,k)*cc2 &
                     )*dr1/grr(1,j,k) &
                   + (cc1*prs(3,j,k) &
                     +cc0*prs(2,j,k) &
                     -cc2*prs(1,j,k) &
                     )*dr1
        dvz(1,j,k) = ( vz(3,j,k)*vdrftr(3,j,k)*grr(3,j,k)*cc1 &
                      +vz(2,j,k)*vdrftr(2,j,k)*grr(2,j,k)*cc0 &
                      -vz(1,j,k)*vdrftr(1,j,k)*grr(1,j,k)*cc2 &
                     )*dr1/grr(1,j,k)
        dvphi(1,j,k) = ( vphi(3,j,k)*vdrftr(3,j,k)*grr(3,j,k)*cc1 &
                        +vphi(2,j,k)*vdrftr(2,j,k)*grr(2,j,k)*cc0 &
                        -vphi(1,j,k)*vdrftr(1,j,k)*grr(1,j,k)*cc2 &
                       )*dr1/grr(1,j,k)

        dvr(2,j,k) = (vr(4,j,k)*vdrftr(4,j,k)*grr(4,j,k)*cc1 &
                     +vr(3,j,k)*vdrftr(3,j,k)*grr(3,j,k)*cc0 &
                     +vr(2,j,k)*vdrftr(2,j,k)*grr(2,j,k)*cc1 &
                     -vr(1,j,k)*vdrftr(1,j,k)*grr(1,j,k)     &
                     )*dr1/grr(2,j,k) &
                   + (cc1*prs(4,j,k) &
                     +cc0*prs(3,j,k) &
                     +cc1*prs(2,j,k) &
                     -    prs(1,j,k) &
                     )*dr1
        dvz(2,j,k) = (vz(4,j,k)*vdrftr(4,j,k)*grr(4,j,k)*cc1 &
                     +vz(3,j,k)*vdrftr(3,j,k)*grr(3,j,k)*cc0 &
                     +vz(2,j,k)*vdrftr(2,j,k)*grr(2,j,k)*cc1 &
                     -vz(1,j,k)*vdrftr(1,j,k)*grr(1,j,k)     &
                     )*dr1/grr(2,j,k)
        dvphi(2,j,k) = (vphi(4,j,k)*vdrftr(4,j,k)*grr(4,j,k)*cc1 &
                       +vphi(3,j,k)*vdrftr(3,j,k)*grr(3,j,k)*cc0 &
                       +vphi(2,j,k)*vdrftr(2,j,k)*grr(2,j,k)*cc1 &
                       -vphi(1,j,k)*vdrftr(1,j,k)*grr(1,j,k)     &
                       )*dr1/grr(2,j,k)

        dvr(lr-1,j,k) = (vr(lr,  j,k)*vdrftr(lr,  j,k)*grr(lr,  j,k)     &
                        -vr(lr-1,j,k)*vdrftr(lr-1,j,k)*grr(lr-1,j,k)*cc1 &
                        -vr(lr-2,j,k)*vdrftr(lr-2,j,k)*grr(lr-2,j,k)*cc0 &
                        -vr(lr-3,j,k)*vdrftr(lr-3,j,k)*grr(lr-3,j,k)*cc1 &
                        )*dr1/grr(lr-1,j,k) &
                   + (    prs(lr,j,k) &
                     -cc1*prs(lr-1,j,k) &
                     -cc0*prs(lr-2,j,k) &
                     -cc1*prs(lr-3,j,k) &
                     )*dr1
        dvz(lr-1,j,k) = (vz(lr,  j,k)*vdrftr(lr,  j,k)*grr(lr,  j,k)     &
                        -vz(lr-1,j,k)*vdrftr(lr-1,j,k)*grr(lr-1,j,k)*cc1 &
                        -vz(lr-2,j,k)*vdrftr(lr-2,j,k)*grr(lr-2,j,k)*cc0 &
                        -vz(lr-3,j,k)*vdrftr(lr-3,j,k)*grr(lr-3,j,k)*cc1 &
                        )*dr1/grr(lr-1,j,k)
        dvphi(lr-1,j,k) = (vphi(lr,  j,k)*vdrftr(lr,  j,k)*grr(lr,  j,k)     &
                          -vphi(lr-1,j,k)*vdrftr(lr-1,j,k)*grr(lr-1,j,k)*cc1 &
                          -vphi(lr-2,j,k)*vdrftr(lr-2,j,k)*grr(lr-2,j,k)*cc0 &
                          -vphi(lr-3,j,k)*vdrftr(lr-3,j,k)*grr(lr-3,j,k)*cc1 &
                          )*dr1/grr(lr-1,j,k)

        dvr(lr,  j,k) = (vr(lr,  j,k)*vdrftr(lr,  j,k)*grr(lr,  j,k)*cc2 &
                        -vr(lr-1,j,k)*vdrftr(lr-1,j,k)*grr(lr-1,j,k)*cc0 &
                        -vr(lr-2,j,k)*vdrftr(lr-2,j,k)*grr(lr-2,j,k)*cc1 &
                        )*dr1/grr(lr,j,k) &
                   + (cc2*prs(lr,j,k) &
                     -cc0*prs(lr-1,j,k) &
                     -cc1*prs(lr-2,j,k) &
                     )*dr1
        dvz(lr,  j,k) = (vz(lr,  j,k)*vdrftr(lr,  j,k)*grr(lr,  j,k)*cc2 &
                        -vz(lr-1,j,k)*vdrftr(lr-1,j,k)*grr(lr-1,j,k)*cc0 &
                        -vz(lr-2,j,k)*vdrftr(lr-2,j,k)*grr(lr-2,j,k)*cc1 &
                        )*dr1/grr(lr,j,k)
        dvphi(lr,  j,k) = (vphi(lr,  j,k)*vdrftr(lr,  j,k)*grr(lr,  j,k)*cc2 &
                          -vphi(lr-1,j,k)*vdrftr(lr-1,j,k)*grr(lr-1,j,k)*cc0 &
                          -vphi(lr-2,j,k)*vdrftr(lr-2,j,k)*grr(lr-2,j,k)*cc1 &
                          )*dr1/grr(lr,j,k)
      end do
      end do

! z direction

#ifdef AMDGPU
!$omp target teams distribute parallel do collapse(2) private(k,i,i0,i1,i2,i3)
#else
!$omp parallel do private(i0,i1,i2,i3)
#endif
      do k = 3, lphi-2
      do i = 2*lr+1, lrz-2*lr

        i3 = i + 2*lr
        i2 = i + lr
        i1 = i - lr
        i0 = i - 2*lr
        dvr(i,1,k) = dvr(i,1,k) &
                   + (vr(i2,1,k)*vdrftz(i2,1,k)*cc0 + vr(i3,1,k)*vdrftz(i3,1,k)*cc1 &
                     -vr(i1,1,k)*vdrftz(i1,1,k)*cc0 - vr(i0,1,k)*vdrftz(i0,1,k)*cc1 )*dz1
        dvz(i,1,k) = dvz(i,1,k) &
                   + (vz(i2,1,k)*vdrftz(i2,1,k)*cc0 + vz(i3,1,k)*vdrftz(i3,1,k)*cc1 &
                     -vz(i1,1,k)*vdrftz(i1,1,k)*cc0 - vz(i0,1,k)*vdrftz(i0,1,k)*cc1 )*dz1 &
                   + (prs(i2,1,k)*cc0 + prs(i3,1,k)*cc1 &
                     -prs(i1,1,k)*cc0 - prs(i0,1,k)*cc1 )*dz1
        dvphi(i,1,k) = dvphi(i,1,k) &
                   + (vphi(i2,1,k)*vdrftz(i2,1,k)*cc0 + vphi(i3,1,k)*vdrftz(i3,1,k)*cc1 &
                     -vphi(i1,1,k)*vdrftz(i1,1,k)*cc0 - vphi(i0,1,k)*vdrftz(i0,1,k)*cc1 )*dz1
      end do
      end do

#ifdef AMDGPU
!$omp target teams distribute parallel do collapse(2) private(k,i)
#else
!$omp parallel do
#endif
      do k = 3, lphi-2
      do i = 1, lr
        dvr(i,1,   k) = dvr(i,1,   k) &
              + (vr(i,3,k)*vdrftz(i,3,k)*cc1 + vr(i,2,k)*vdrftz(i,2,k)*cc0 &
               - vr(i,1,1)*vdrftz(i,1,k)*cc2 )*dz1
        dvz(i,1,   k) = dvz(i,1,   k) &
              + (vz(i,3,k)*vdrftz(i,3,k)*cc1 + vz(i,2,k)*vdrftz(i,2,k)*cc0 &
               - vz(i,1,1)*vdrftz(i,1,k)*cc2 )*dz1 &
              + (prs(i,3,k)*cc1 + prs(i,2,k)*cc0 - prs(i,1,1)*cc2 )*dz1
        dvphi(i,1,   k) = dvphi(i,1,   k) &
              + (vphi(i,3,k)*vdrftz(i,3,k)*cc1 + vphi(i,2,k)*vdrftz(i,2,k)*cc0 &
               - vphi(i,1,1)*vdrftz(i,1,k)*cc2 )*dz1

        dvr(i,2,   k) = dvr(i,2,   k) &
              + (vr(i,4,k)*vdrftz(i,4,k)*cc1 + vr(i,3,k)*vdrftz(i,3,k)*cc0 &
                +vr(i,2,k)*vdrftz(i,2,k)*cc1 - vr(i,1,k)*vdrftz(i,1,k))*dz1
        dvz(i,2,   k) = dvz(i,2,   k) &
              + (vz(i,4,k)*vdrftz(i,4,k)*cc1 + vz(i,3,k)*vdrftz(i,3,k)*cc0 &
                +vz(i,2,k)*vdrftz(i,2,k)*cc1 - vz(i,1,k)*vdrftz(i,1,k))*dz1 &
              + (prs(i,4,k)*cc1 + prs(i,3,k)*cc0 &
                +prs(i,2,k)*cc1 - prs(i,1,k))*dz1
        dvphi(i,2,   k) = dvphi(i,2,   k) &
              + (vphi(i,4,k)*vdrftz(i,4,k)*cc1 + vphi(i,3,k)*vdrftz(i,3,k)*cc0 &
                +vphi(i,2,k)*vdrftz(i,2,k)*cc1 - vphi(i,1,k)*vdrftz(i,1,k))*dz1

        dvr(i,lz-1,k) = dvr(i,lz-1,k) &
              + (vr(i,lz,k)*vdrftz(i,lz,k) - vr(i,lz-1,k)*vdrftz(i,lz-1,k)*cc1 &
                -vr(i,lz-2,k)*vdrftz(i,lz-2,k)*cc0 - vr(i,lz-3,k)*vdrftz(i,lz-3,k)*cc1 )*dz1
        dvz(i,lz-1,k) = dvz(i,lz-1,k) &
              + (vz(i,lz,k)*vdrftz(i,lz,k) - vz(i,lz-1,k)*vdrftz(i,lz-1,k)*cc1 &
                -vz(i,lz-2,k)*vdrftz(i,lz-2,k)*cc0 - vz(i,lz-3,k)*vdrftz(i,lz-3,k)*cc1 )*dz1 &
              + (prs(i,lz,k) - prs(i,lz-1,k)*cc1 &
                -prs(i,lz-2,k)*cc0 - prs(i,lz-3,k)*cc1 )*dz1
        dvphi(i,lz-1,k) = dvphi(i,lz-1,k) &
              + (vphi(i,lz,k)*vdrftz(i,lz,k) - vphi(i,lz-1,k)*vdrftz(i,lz-1,k)*cc1 &
                -vphi(i,lz-2,k)*vdrftz(i,lz-2,k)*cc0 - vphi(i,lz-3,k)*vdrftz(i,lz-3,k)*cc1 )*dz1

        dvr(i,lz,  k) = dvr(i,lz,  k) &
              + (vr(i,lz,k)*vdrftz(i,lz,k)*cc2 - vr(i,lz-1,k)*vdrftz(i,lz-1,k)*cc0 &
               - vr(i,lz-2,k)*vdrftz(i,lz-2,k)*cc1 )*dz1
        dvz(i,lz,  k) = dvz(i,lz,  k) &
              + (vz(i,lz,k)*vdrftz(i,lz,k)*cc2 - vz(i,lz-1,k)*vdrftz(i,lz-1,k)*cc0 &
               - vz(i,lz-2,k)*vdrftz(i,lz-2,k)*cc1 )*dz1 &
              + (prs(i,lz,k)*cc2 - prs(i,lz-1,k)*cc0 - prs(i,lz-2,k)*cc1 )*dz1
        dvphi(i,lz,  k) = dvphi(i,lz,  k) &
              + (vphi(i,lz,k)*vdrftz(i,lz,k)*cc2 - vphi(i,lz-1,k)*vdrftz(i,lz-1,k)*cc0 &
               - vphi(i,lz-2,k)*vdrftz(i,lz-2,k)*cc1 )*dz1
      end do
      end do

! phi direction

#ifdef AMDGPU
!$omp target teams distribute parallel do collapse(2) private(k,i,k3,k2,k1,k0)
#else
!$omp parallel do private(k0,k1,k2,k3)
#endif
      do k = 3, lphi-2
      do i = 1, lrz
        k3 = k + 2
        k2 = k + 1
        k1 = k - 1
        k0 = k - 2

        dvr(i,1,k) = dvr(i,1,k) &
                   + (vr(i,1,k2)*vdrftphi(i,1,k2)*cc0 + vr(i,1,k3)*vdrftphi(i,1,k3)*cc1  &
                     -vr(i,1,k1)*vdrftphi(i,1,k1)*cc0 - vr(i,1,k0)*vdrftphi(i,1,k0)*cc1  &
                     )*dphi1/grr(i,1,k)
        dvz(i,1,k) = dvz(i,1,k) &
                   + (vz(i,1,k2)*vdrftphi(i,1,k2)*cc0 + vz(i,1,k3)*vdrftphi(i,1,k3)*cc1  &
                     -vz(i,1,k1)*vdrftphi(i,1,k1)*cc0 - vz(i,1,k0)*vdrftphi(i,1,k0)*cc1  &
                     )*dphi1/grr(i,1,k)
        dvphi(i,1,k) = dvphi(i,1,k) &
                   + (vphi(i,1,k2)*vdrftphi(i,1,k2)*cc0 + vphi(i,1,k3)*vdrftphi(i,1,k3)*cc1  &
                     -vphi(i,1,k1)*vdrftphi(i,1,k1)*cc0 - vphi(i,1,k0)*vdrftphi(i,1,k0)*cc1  &
                   +  prs(i,1,k2)*cc0 + prs(i,1,k3)*cc1  &
                     -prs(i,1,k1)*cc0 - prs(i,1,k0)*cc1  &
                     )*dphi1/grr(i,1,k)
      end do
      end do


      end
!--------------------------------------------------------------------
subroutine divergence4hm
! 2012-08-16 for dvr,dvz,dvphi,dprs in subroutine hm
! 2014-01-10 toroidal rotation is introduced, in utotphi
!--------------------------------------------------------------------
      use parameters
      use grid
      use field, only:vr,vz,vphi,prs,vdrftr,vdrftz,vdrftphi &
                     ,ur,uz,utotphi,div_u,dvr,dvz,dvphi,dprs,gamma
      implicit none

      real(8)::dr1,dz1,dphi1
      real(8)::cc0=4.0d0/3.0d0,cc1=-1.0d0/6.0d0,cc2=7.0d0/6.0d0
      integer::i,j,k,i0,i1,i2,i3,k0,k1,k2,k3

      dr1= 0.50d0/dr
      dz1= 0.50d0/dz
      dphi1= 0.50d0/dphi
!      cc0= 4.0d0/3.0d0
!      cc1=-1.0d0/6.0d0 
!      cc2= 7.0d0/6.0d0 

! r direction

#ifdef AMDGPU
!$omp target teams distribute parallel do collapse(3) private(k,j,i,i0,i1,i2,i3)
#else
!$omp parallel do private(i0,i1,i2,i3)
#endif
      do k = 3, lphi-2
      do j = 1, lz
      do i = 3, lr-2
        i3 = i + 2
        i2 = i + 1
        i1 = i - 1
        i0 = i - 2

        dvr(i,j,k) = (cc0*vr(i2,j,k)*vdrftr(i2,j,k)*grr(i2,j,k) &
                     +cc1*vr(i3,j,k)*vdrftr(i3,j,k)*grr(i3,j,k) &
                     -cc0*vr(i1,j,k)*vdrftr(i1,j,k)*grr(i1,j,k) &
                     -cc1*vr(i0,j,k)*vdrftr(i0,j,k)*grr(i0,j,k) &
                     )*dr1/grr(i,j,k) &
                   + (cc0*prs(i2,j,k) &
                     +cc1*prs(i3,j,k) &
                     -cc0*prs(i1,j,k) &
                     -cc1*prs(i0,j,k) &
                     )*dr1
        dvz(i,j,k) = (cc0*vz(i2,j,k)*vdrftr(i2,j,k)*grr(i2,j,k) &
                     +cc1*vz(i3,j,k)*vdrftr(i3,j,k)*grr(i3,j,k) &
                     -cc0*vz(i1,j,k)*vdrftr(i1,j,k)*grr(i1,j,k) &
                     -cc1*vz(i0,j,k)*vdrftr(i0,j,k)*grr(i0,j,k) &
                     )*dr1/grr(i,j,k)
        dvphi(i,j,k) = (cc0*vphi(i2,j,k)*vdrftr(i2,j,k)*grr(i2,j,k) &
                       +cc1*vphi(i3,j,k)*vdrftr(i3,j,k)*grr(i3,j,k) &
                       -cc0*vphi(i1,j,k)*vdrftr(i1,j,k)*grr(i1,j,k) &
                       -cc1*vphi(i0,j,k)*vdrftr(i0,j,k)*grr(i0,j,k) &
                       )*dr1/grr(i,j,k)
        dprs(i,j,k) = (cc0*prs(i2,j,k)*ur(i2,j,k)*grr(i2,j,k) &
                      +cc1*prs(i3,j,k)*ur(i3,j,k)*grr(i3,j,k) &
                      -cc0*prs(i1,j,k)*ur(i1,j,k)*grr(i1,j,k) &
                      -cc1*prs(i0,j,k)*ur(i0,j,k)*grr(i0,j,k) &
                      )*dr1/grr(i,j,k)
      end do
      end do
      end do

#ifdef AMDGPU
!$omp target teams distribute parallel do collapse(2) private(k,j)
#else
!$omp parallel do
#endif
      do k = 3, lphi-2
      do j = 1, lz
        dvr(1,j,k) = ( vr(3,j,k)*vdrftr(3,j,k)*grr(3,j,k)*cc1 &
                      +vr(2,j,k)*vdrftr(2,j,k)*grr(2,j,k)*cc0 &
                      -vr(1,j,k)*vdrftr(1,j,k)*grr(1,j,k)*cc2 &
                     )*dr1/grr(1,j,k) &
                   + (cc1*prs(3,j,k) &
                     +cc0*prs(2,j,k) &
                     -cc2*prs(1,j,k) &
                     )*dr1
        dvz(1,j,k) = ( vz(3,j,k)*vdrftr(3,j,k)*grr(3,j,k)*cc1 &
                      +vz(2,j,k)*vdrftr(2,j,k)*grr(2,j,k)*cc0 &
                      -vz(1,j,k)*vdrftr(1,j,k)*grr(1,j,k)*cc2 &
                     )*dr1/grr(1,j,k)
        dvphi(1,j,k) = ( vphi(3,j,k)*vdrftr(3,j,k)*grr(3,j,k)*cc1 &
                        +vphi(2,j,k)*vdrftr(2,j,k)*grr(2,j,k)*cc0 &
                        -vphi(1,j,k)*vdrftr(1,j,k)*grr(1,j,k)*cc2 &
                       )*dr1/grr(1,j,k)
        dprs(1,j,k) = ( prs(3,j,k)*ur(3,j,k)*grr(3,j,k)*cc1 &
                       +prs(2,j,k)*ur(2,j,k)*grr(2,j,k)*cc0 &
                       -prs(1,j,k)*ur(1,j,k)*grr(1,j,k)*cc2 &
                      )*dr1/grr(1,j,k)

        dvr(2,j,k) = (vr(4,j,k)*vdrftr(4,j,k)*grr(4,j,k)*cc1 &
                     +vr(3,j,k)*vdrftr(3,j,k)*grr(3,j,k)*cc0 &
                     +vr(2,j,k)*vdrftr(2,j,k)*grr(2,j,k)*cc1 &
                     -vr(1,j,k)*vdrftr(1,j,k)*grr(1,j,k)     &
                     )*dr1/grr(2,j,k) &
                   + (cc1*prs(4,j,k) &
                     +cc0*prs(3,j,k) &
                     +cc1*prs(2,j,k) &
                     -    prs(1,j,k) &
                     )*dr1
        dvz(2,j,k) = (vz(4,j,k)*vdrftr(4,j,k)*grr(4,j,k)*cc1 &
                     +vz(3,j,k)*vdrftr(3,j,k)*grr(3,j,k)*cc0 &
                     +vz(2,j,k)*vdrftr(2,j,k)*grr(2,j,k)*cc1 &
                     -vz(1,j,k)*vdrftr(1,j,k)*grr(1,j,k)     &
                     )*dr1/grr(2,j,k)
        dvphi(2,j,k) = (vphi(4,j,k)*vdrftr(4,j,k)*grr(4,j,k)*cc1 &
                       +vphi(3,j,k)*vdrftr(3,j,k)*grr(3,j,k)*cc0 &
                       +vphi(2,j,k)*vdrftr(2,j,k)*grr(2,j,k)*cc1 &
                       -vphi(1,j,k)*vdrftr(1,j,k)*grr(1,j,k)     &
                       )*dr1/grr(2,j,k)
        dprs(2,j,k) = (prs(4,j,k)*ur(4,j,k)*grr(4,j,k)*cc1 &
                      +prs(3,j,k)*ur(3,j,k)*grr(3,j,k)*cc0 &
                      +prs(2,j,k)*ur(2,j,k)*grr(2,j,k)*cc1 &
                      -prs(1,j,k)*ur(1,j,k)*grr(1,j,k)     &
                      )*dr1/grr(2,j,k)

        dvr(lr-1,j,k) = (vr(lr,  j,k)*vdrftr(lr,  j,k)*grr(lr,  j,k)     &
                        -vr(lr-1,j,k)*vdrftr(lr-1,j,k)*grr(lr-1,j,k)*cc1 &
                        -vr(lr-2,j,k)*vdrftr(lr-2,j,k)*grr(lr-2,j,k)*cc0 &
                        -vr(lr-3,j,k)*vdrftr(lr-3,j,k)*grr(lr-3,j,k)*cc1 &
                        )*dr1/grr(lr-1,j,k) &
                   + (    prs(lr,j,k) &
                     -cc1*prs(lr-1,j,k) &
                     -cc0*prs(lr-2,j,k) &
                     -cc1*prs(lr-3,j,k) &
                     )*dr1
        dvz(lr-1,j,k) = (vz(lr,  j,k)*vdrftr(lr,  j,k)*grr(lr,  j,k)     &
                        -vz(lr-1,j,k)*vdrftr(lr-1,j,k)*grr(lr-1,j,k)*cc1 &
                        -vz(lr-2,j,k)*vdrftr(lr-2,j,k)*grr(lr-2,j,k)*cc0 &
                        -vz(lr-3,j,k)*vdrftr(lr-3,j,k)*grr(lr-3,j,k)*cc1 &
                        )*dr1/grr(lr-1,j,k)
        dvphi(lr-1,j,k) = (vphi(lr,  j,k)*vdrftr(lr,  j,k)*grr(lr,  j,k)     &
                          -vphi(lr-1,j,k)*vdrftr(lr-1,j,k)*grr(lr-1,j,k)*cc1 &
                          -vphi(lr-2,j,k)*vdrftr(lr-2,j,k)*grr(lr-2,j,k)*cc0 &
                          -vphi(lr-3,j,k)*vdrftr(lr-3,j,k)*grr(lr-3,j,k)*cc1 &
                          )*dr1/grr(lr-1,j,k)
        dprs(lr-1,j,k) = (prs(lr,  j,k)*ur(lr,  j,k)*grr(lr,  j,k)     &
                         -prs(lr-1,j,k)*ur(lr-1,j,k)*grr(lr-1,j,k)*cc1 &
                         -prs(lr-2,j,k)*ur(lr-2,j,k)*grr(lr-2,j,k)*cc0 &
                         -prs(lr-3,j,k)*ur(lr-3,j,k)*grr(lr-3,j,k)*cc1 &
                         )*dr1/grr(lr-1,j,k)

        dvr(lr,  j,k) = (vr(lr,  j,k)*vdrftr(lr,  j,k)*grr(lr,  j,k)*cc2 &
                        -vr(lr-1,j,k)*vdrftr(lr-1,j,k)*grr(lr-1,j,k)*cc0 &
                        -vr(lr-2,j,k)*vdrftr(lr-2,j,k)*grr(lr-2,j,k)*cc1 &
                        )*dr1/grr(lr,j,k) &
                   + (cc2*prs(lr,j,k) &
                     -cc0*prs(lr-1,j,k) &
                     -cc1*prs(lr-2,j,k) &
                     )*dr1
        dvz(lr,  j,k) = (vz(lr,  j,k)*vdrftr(lr,  j,k)*grr(lr,  j,k)*cc2 &
                        -vz(lr-1,j,k)*vdrftr(lr-1,j,k)*grr(lr-1,j,k)*cc0 &
                        -vz(lr-2,j,k)*vdrftr(lr-2,j,k)*grr(lr-2,j,k)*cc1 &
                        )*dr1/grr(lr,j,k)
        dvphi(lr,  j,k) = (vphi(lr,  j,k)*vdrftr(lr,  j,k)*grr(lr,  j,k)*cc2 &
                          -vphi(lr-1,j,k)*vdrftr(lr-1,j,k)*grr(lr-1,j,k)*cc0 &
                          -vphi(lr-2,j,k)*vdrftr(lr-2,j,k)*grr(lr-2,j,k)*cc1 &
                          )*dr1/grr(lr,j,k)
        dprs(lr,  j,k) = (prs(lr,  j,k)*ur(lr,  j,k)*grr(lr,  j,k)*cc2 &
                         -prs(lr-1,j,k)*ur(lr-1,j,k)*grr(lr-1,j,k)*cc0 &
                         -prs(lr-2,j,k)*ur(lr-2,j,k)*grr(lr-2,j,k)*cc1 &
                         )*dr1/grr(lr,j,k)
      end do
      end do

! z direction

#ifdef AMDGPU
!$omp target teams distribute parallel do collapse(2) private(k,i,i0,i1,i2,i3)
#else
!$omp parallel do private(i0,i1,i2,i3)
#endif
      do k = 3, lphi-2
      do i = 2*lr+1, lrz-2*lr

        i3 = i + 2*lr
        i2 = i + lr
        i1 = i - lr
        i0 = i - 2*lr
        dvr(i,1,k) = dvr(i,1,k) &
                   + (vr(i2,1,k)*vdrftz(i2,1,k)*cc0 + vr(i3,1,k)*vdrftz(i3,1,k)*cc1 &
                     -vr(i1,1,k)*vdrftz(i1,1,k)*cc0 - vr(i0,1,k)*vdrftz(i0,1,k)*cc1 )*dz1
        dvz(i,1,k) = dvz(i,1,k) &
                   + (vz(i2,1,k)*vdrftz(i2,1,k)*cc0 + vz(i3,1,k)*vdrftz(i3,1,k)*cc1 &
                     -vz(i1,1,k)*vdrftz(i1,1,k)*cc0 - vz(i0,1,k)*vdrftz(i0,1,k)*cc1 )*dz1 &
                   + (prs(i2,1,k)*cc0 + prs(i3,1,k)*cc1 &
                     -prs(i1,1,k)*cc0 - prs(i0,1,k)*cc1 )*dz1
        dvphi(i,1,k) = dvphi(i,1,k) &
                   + (vphi(i2,1,k)*vdrftz(i2,1,k)*cc0 + vphi(i3,1,k)*vdrftz(i3,1,k)*cc1 &
                     -vphi(i1,1,k)*vdrftz(i1,1,k)*cc0 - vphi(i0,1,k)*vdrftz(i0,1,k)*cc1 )*dz1
        dprs(i,1,k) = dprs(i,1,k) &
                   + (prs(i2,1,k)*uz(i2,1,k)*cc0 + prs(i3,1,k)*uz(i3,1,k)*cc1 &
                     -prs(i1,1,k)*uz(i1,1,k)*cc0 - prs(i0,1,k)*uz(i0,1,k)*cc1 )*dz1

      end do
      end do

#ifdef AMDGPU
!$omp target teams distribute parallel do collapse(2) private(k,i)
#else
!$omp parallel do
#endif
      do k = 3, lphi-2
      do i = 1, lr
        dvr(i,1,   k) = dvr(i,1,   k) &
              + (vr(i,3,k)*vdrftz(i,3,k)*cc1 + vr(i,2,k)*vdrftz(i,2,k)*cc0 &
               - vr(i,1,1)*vdrftz(i,1,k)*cc2 )*dz1
        dvz(i,1,   k) = dvz(i,1,   k) &
              + (vz(i,3,k)*vdrftz(i,3,k)*cc1 + vz(i,2,k)*vdrftz(i,2,k)*cc0 &
               - vz(i,1,1)*vdrftz(i,1,k)*cc2 )*dz1 &
              + (prs(i,3,k)*cc1 + prs(i,2,k)*cc0 - prs(i,1,1)*cc2 )*dz1
        dvphi(i,1,   k) = dvphi(i,1,   k) &
              + (vphi(i,3,k)*vdrftz(i,3,k)*cc1 + vphi(i,2,k)*vdrftz(i,2,k)*cc0 &
               - vphi(i,1,1)*vdrftz(i,1,k)*cc2 )*dz1
        dprs(i,1,   k) = dprs(i,1,   k) &
              + (prs(i,3,k)*uz(i,3,k)*cc1 + prs(i,2,k)*uz(i,2,k)*cc0 - prs(i,1,1)*uz(i,1,k)*cc2 )*dz1


        dvr(i,2,   k) = dvr(i,2,   k) &
              + (vr(i,4,k)*vdrftz(i,4,k)*cc1 + vr(i,3,k)*vdrftz(i,3,k)*cc0 &
                +vr(i,2,k)*vdrftz(i,2,k)*cc1 - vr(i,1,k)*vdrftz(i,1,k))*dz1
        dvz(i,2,   k) = dvz(i,2,   k) &
              + (vz(i,4,k)*vdrftz(i,4,k)*cc1 + vz(i,3,k)*vdrftz(i,3,k)*cc0 &
                +vz(i,2,k)*vdrftz(i,2,k)*cc1 - vz(i,1,k)*vdrftz(i,1,k))*dz1 &
              + (prs(i,4,k)*cc1 + prs(i,3,k)*cc0 &
                +prs(i,2,k)*cc1 - prs(i,1,k))*dz1
        dvphi(i,2,   k) = dvphi(i,2,   k) &
              + (vphi(i,4,k)*vdrftz(i,4,k)*cc1 + vphi(i,3,k)*vdrftz(i,3,k)*cc0 &
                +vphi(i,2,k)*vdrftz(i,2,k)*cc1 - vphi(i,1,k)*vdrftz(i,1,k))*dz1
        dprs(i,2,   k) = dprs(i,2,   k) &
              + (prs(i,4,k)*uz(i,4,k)*cc1 + prs(i,3,k)*uz(i,3,k)*cc0 &
                +prs(i,2,k)*uz(i,2,k)*cc1 - prs(i,1,k)*uz(i,1,k))*dz1

        dvr(i,lz-1,k) = dvr(i,lz-1,k) &
              + (vr(i,lz,k)*vdrftz(i,lz,k) - vr(i,lz-1,k)*vdrftz(i,lz-1,k)*cc1 &
                -vr(i,lz-2,k)*vdrftz(i,lz-2,k)*cc0 - vr(i,lz-3,k)*vdrftz(i,lz-3,k)*cc1 )*dz1
        dvz(i,lz-1,k) = dvz(i,lz-1,k) &
              + (vz(i,lz,k)*vdrftz(i,lz,k) - vz(i,lz-1,k)*vdrftz(i,lz-1,k)*cc1 &
                -vz(i,lz-2,k)*vdrftz(i,lz-2,k)*cc0 - vz(i,lz-3,k)*vdrftz(i,lz-3,k)*cc1 )*dz1 &
              + (prs(i,lz,k) - prs(i,lz-1,k)*cc1 &
                -prs(i,lz-2,k)*cc0 - prs(i,lz-3,k)*cc1 )*dz1
        dvphi(i,lz-1,k) = dvphi(i,lz-1,k) &
              + (vphi(i,lz,k)*vdrftz(i,lz,k) - vphi(i,lz-1,k)*vdrftz(i,lz-1,k)*cc1 &
                -vphi(i,lz-2,k)*vdrftz(i,lz-2,k)*cc0 - vphi(i,lz-3,k)*vdrftz(i,lz-3,k)*cc1 )*dz1
        dprs(i,lz-1,k) = dprs(i,lz-1,k) &
              + (prs(i,lz,k)*uz(i,lz,k) - prs(i,lz-1,k)*uz(i,lz-1,k)*cc1 &
                -prs(i,lz-2,k)*uz(i,lz-2,k)*cc0 - prs(i,lz-3,k)*uz(i,lz-3,k)*cc1 )*dz1

        dvr(i,lz,  k) = dvr(i,lz,  k) &
              + (vr(i,lz,k)*vdrftz(i,lz,k)*cc2 - vr(i,lz-1,k)*vdrftz(i,lz-1,k)*cc0 &
               - vr(i,lz-2,k)*vdrftz(i,lz-2,k)*cc1 )*dz1
        dvz(i,lz,  k) = dvz(i,lz,  k) &
              + (vz(i,lz,k)*vdrftz(i,lz,k)*cc2 - vz(i,lz-1,k)*vdrftz(i,lz-1,k)*cc0 &
               - vz(i,lz-2,k)*vdrftz(i,lz-2,k)*cc1 )*dz1 &
              + (prs(i,lz,k)*cc2 - prs(i,lz-1,k)*cc0 - prs(i,lz-2,k)*cc1 )*dz1
        dvphi(i,lz,  k) = dvphi(i,lz,  k) &
              + (vphi(i,lz,k)*vdrftz(i,lz,k)*cc2 - vphi(i,lz-1,k)*vdrftz(i,lz-1,k)*cc0 &
               - vphi(i,lz-2,k)*vdrftz(i,lz-2,k)*cc1 )*dz1
        dprs(i,lz,  k) = dprs(i,lz,  k) &
              + (prs(i,lz,k)*uz(i,lz,k)*cc2 - prs(i,lz-1,k)*uz(i,lz-1,k)*cc0 &
               - prs(i,lz-2,k)*uz(i,lz-2,k)*cc1 )*dz1

      end do
      end do

! phi direction

#ifdef AMDGPU
!$omp target teams distribute parallel do collapse(2) private(k,i,k3,k2,k1,k0)
#else
!$omp parallel do private(k0,k1,k2,k3)
#endif
      do k = 3, lphi-2
      do i = 1, lrz
        k3 = k + 2
        k2 = k + 1
        k1 = k - 1
        k0 = k - 2

        dvr(i,1,k) = dvr(i,1,k) &
                   + (vr(i,1,k2)*vdrftphi(i,1,k2)*cc0 + vr(i,1,k3)*vdrftphi(i,1,k3)*cc1  &
                     -vr(i,1,k1)*vdrftphi(i,1,k1)*cc0 - vr(i,1,k0)*vdrftphi(i,1,k0)*cc1  &
                     )*dphi1/grr(i,1,k)
        dvz(i,1,k) = dvz(i,1,k) &
                   + (vz(i,1,k2)*vdrftphi(i,1,k2)*cc0 + vz(i,1,k3)*vdrftphi(i,1,k3)*cc1  &
                     -vz(i,1,k1)*vdrftphi(i,1,k1)*cc0 - vz(i,1,k0)*vdrftphi(i,1,k0)*cc1  &
                     )*dphi1/grr(i,1,k)
        dvphi(i,1,k) = dvphi(i,1,k) &
                   + (vphi(i,1,k2)*vdrftphi(i,1,k2)*cc0 + vphi(i,1,k3)*vdrftphi(i,1,k3)*cc1  &
                     -vphi(i,1,k1)*vdrftphi(i,1,k1)*cc0 - vphi(i,1,k0)*vdrftphi(i,1,k0)*cc1  &
                   +  prs(i,1,k2)*cc0 + prs(i,1,k3)*cc1  &
                     -prs(i,1,k1)*cc0 - prs(i,1,k0)*cc1  &
                     )*dphi1/grr(i,1,k)
        dprs(i,1,k) = dprs(i,1,k) &
                   + (prs(i,1,k2)*utotphi(i,1,k2)*cc0 + prs(i,1,k3)*utotphi(i,1,k3)*cc1  &
                     -prs(i,1,k1)*utotphi(i,1,k1)*cc0 - prs(i,1,k0)*utotphi(i,1,k0)*cc1  &
                     )*dphi1/grr(i,1,k) &
                    + (gamma-1.0d0)*prs(i,1,k)*div_u(i,1,k)
      end do
      end do


      end
!--------------------------------------------------------------------
subroutine rotation(iflag,a_r,a_z,a_phi,b_r,b_z,b_phi)
!--------------------------------------------------------------------
      use parameters
      use grid
      implicit none

      real(8)::a_r(lr,lz,lphi),a_z(lr,lz,lphi),a_phi(lr,lz,lphi)
      real(8)::b_r(lr,lz,lphi),b_z(lr,lz,lphi),b_phi(lr,lz,lphi)
      real(8)::dr1,dz1,dphi1
      real(8)::cc0=4.0d0/3.0d0,cc1=-1.0d0/6.0d0,cc2=7.0d0/6.0d0
      integer::i,j,k,i0,i1,i2,i3,k0,k1,k2,k3
      integer::iflag

      dr1= 0.50d0/dr
      dz1= 0.50d0/dz
      dphi1= 0.50d0/dphi
!      cc0= 4.0d0/3.0d0
!      cc1=-1.0d0/6.0d0
!      cc2= 7.0d0/6.0d0


! derivative in r direction

#ifdef AMDGPU
!$omp target teams distribute parallel do collapse(3) private(k,j,i,i0,i1,i2,i3)
#else
!$omp parallel do private(i0,i1,i2,i3)
#endif
      do k = 3, lphi-2
      do j = 1, lz
      do i = 3, lr-2
        i3 = i + 2
        i2 = i + 1
        i1 = i - 1
        i0 = i - 2

        b_phi(i,j,k) = -(cc0*a_z(i2,j,k) + cc1*a_z(i3,j,k)      &
                        -cc0*a_z(i1,j,k) - cc1*a_z(i0,j,k) )*dr1
        b_z(i,j,k) = (cc0*grr(i2,j,k)*a_phi(i2,j,k) &
                     +cc1*grr(i3,j,k)*a_phi(i3,j,k) &
                     -cc0*grr(i1,j,k)*a_phi(i1,j,k) &
                     -cc1*grr(i0,j,k)*a_phi(i0,j,k) &
                     )*dr1/grr(i,j,k)
      end do
      end do
      end do

#ifdef AMDGPU
!$omp target teams distribute parallel do collapse(2) private(k,j)
#else
!$omp parallel do
#endif
      do k = 3, lphi-2
      do j = 1, lz
        b_phi(1,j,k) = -(a_z(3,j,k)*cc1 + a_z(2,j,k)*cc0 - a_z(1,j,k)*cc2)*dr1
        b_z(1,j,k) = (grr(3,j,k)*a_phi(3,j,k)*cc1 &
                     +grr(2,j,k)*a_phi(2,j,k)*cc0 &
                     -grr(1,j,k)*a_phi(1,j,k)*cc2 &
                     )*dr1/grr(1,j,k)

        b_phi(2,j,k) = -(a_z(4,j,k)*cc1 + a_z(3,j,k)*cc0 &
                        +a_z(2,j,k)*cc1 - a_z(1,j,k))*dr1
        b_z(2,j,k) = (grr(4,j,k)*a_phi(4,j,k)*cc1 + grr(3,j,k)*a_phi(3,j,k)*cc0 &
                     +grr(2,j,k)*a_phi(2,j,k)*cc1 - grr(1,j,k)*a_phi(1,j,k) &
                     )*dr1/grr(2,j,k)

        b_phi(lr-1,j,k) = -(a_z(lr,j,k) - a_z(lr-1,j,k)*cc1 &
                           -a_z(lr-2,j,k)*cc0 - a_z(lr-3,j,k)*cc1 )*dr1
        b_z(lr-1,j,k) = (grr(lr,j,k)  *a_phi(lr,j,k) &
                        -grr(lr-1,j,k)*a_phi(lr-1,j,k)*cc1 &
                        -grr(lr-2,j,k)*a_phi(lr-2,j,k)*cc0 &
                        -grr(lr-3,j,k)*a_phi(lr-3,j,k)*cc1 &
                        )*dr1/grr(lr-1,j,k)

        b_phi(lr,j,k) = -(a_z(lr,j,k)*cc2 &
                         -a_z(lr-1,j,k)*cc0 - a_z(lr-2,j,k)*cc1)*dr1
        b_z(lr,j,k) = (grr(lr,j,k)  *a_phi(lr,j,k)*cc2 &
                      -grr(lr-1,j,k)*a_phi(lr-1,j,k)*cc0 &
                      -grr(lr-2,j,k)*a_phi(lr-2,j,k)*cc1 &
                       )*dr1/grr(lr,j,k)
      end do
      end do

! derivative in z direction

#ifdef AMDGPU
!$omp target teams distribute parallel do collapse(2) private(k,i,i0,i1,i2,i3)
#else
!$omp parallel do private(i0,i1,i2,i3)
#endif
      do k = 3, lphi-2
      do i = 2*lr+1, lrz-2*lr
        i3 = i + 2*lr
        i2 = i + lr
        i1 = i - lr
        i0 = i - 2*lr

        b_r(i,1,k) =-(a_phi(i2,1,k)*cc0 + a_phi(i3,1,k)*cc1 &
                     -a_phi(i1,1,k)*cc0 - a_phi(i0,1,k)*cc1 )*dz1
        b_phi(i,1,k) = (a_r(i2,1,k)*cc0 + a_r(i3,1,k)*cc1   &
                       -a_r(i1,1,k)*cc0 - a_r(i0,1,k)*cc1 )*dz1 &
                     + b_phi(i,1,k)
      end do
      end do

#ifdef AMDGPU
!$omp target teams distribute parallel do collapse(2) private(k,i)
#else
!$omp parallel do
#endif
      do k = 3, lphi-2
      do i = 1, lr
        b_r(i,1,k)   =-(a_phi(i,3,k)*cc1 + a_phi(i,2,k)*cc0 &
                       -a_phi(i,1,k)*cc2 )*dz1
        b_phi(i,1,k) = (a_r(i,3,k)*cc1 + a_r(i,2,k)*cc0 &
                       -a_r(i,1,k)*cc2   )*dz1 &
                     + b_phi(i,1,k)

        b_r(i,2,k)   =-(a_phi(i,4,k)*cc1 + a_phi(i,3,k)*cc0 + a_phi(i,2,k)*cc1 &
                       -a_phi(i,1,k) )*dz1
        b_phi(i,2,k) = (a_r(i,4,k)*cc1 + a_r(i,3,k)*cc0 + a_r(i,2,k)*cc1 &
                       -a_r(i,1,k)   )*dz1 &
                     + b_phi(i,2,k)

        b_r(i,lz-1,k)   =-(a_phi(i,lz,k) - a_phi(i,lz-1,k)*cc1 &
                          -a_phi(i,lz-2,k)*cc0 - a_phi(i,lz-3,k)*cc1 )*dz1
        b_phi(i,lz-1,k) = (a_r(i,lz,k) - a_r(i,lz-1,k)*cc1 &
                          -a_r(i,lz-2,k)*cc0 - a_r(i,lz-3,k)*cc1 )*dz1 &
                        + b_phi(i,lz-1,k)

        b_r(i,lz,k)   =-(a_phi(i,lz,k)*cc2 &
                        -a_phi(i,lz-1,k)*cc0 - a_phi(i,lz-2,k)*cc1 )*dz1
        b_phi(i,lz,k) = (a_r(i,lz,k)*cc2 &
                        -a_r(i,lz-1,k)*cc0 - a_r(i,lz-2,k)*cc1 )*dz1 &
                      + b_phi(i,lz,k)
      end do
      end do

! derivative in phi direction

#ifdef AMDGPU
!$omp target teams distribute parallel do collapse(2) private(k,i,k0,k1,k2,k3)
#else
!$omp parallel do private(k0,k1,k2,k3)
#endif
      do k = 3, lphi-2
      do i = 1, lrz
        k3 = k + 2
        k2 = k + 1
        k1 = k - 1
        k0 = k - 2

        b_r(i,1,k) = (a_z(i,1,k2)*cc0 + a_z(i,1,k3)*cc1 &
                     -a_z(i,1,k1)*cc0 - a_z(i,1,k0)*cc1 )*dphi1/grr(i,1,k) &
                   + b_r(i,1,k)
        b_z(i,1,k) =-(a_r(i,1,k2)*cc0 + a_r(i,1,k3)*cc1 &
                     -a_r(i,1,k1)*cc0 - a_r(i,1,k0)*cc1 )*dphi1/grr(i,1,k) &
                   + b_z(i,1,k)
      end do
      end do

      if(iflag.eq.1)then
        call periodic_field_mlt3(b_r,b_z,b_phi) !2012-04-13
!        call periodic_field(b_r)
!        call periodic_field(b_z)
!        call periodic_field(b_phi)
      end if

end
!--------------------------------------------------------------------
subroutine rotation5(iflag,cf,sc,a_r,a_z,a_phi,b_r,b_z,b_phi)
! 2012-04-16 for subroutine mhd
!--------------------------------------------------------------------
      use parameters
      use grid
      implicit none

      real(8)::cf(lr,lz,lphi),sc(lr,lz,lphi)
      real(8)::a_r(lr,lz,lphi),a_z(lr,lz,lphi),a_phi(lr,lz,lphi)
      real(8)::b_r(lr,lz,lphi),b_z(lr,lz,lphi),b_phi(lr,lz,lphi)
      real(8)::dr1,dz1,dphi1
      real(8)::cc0=4.0d0/3.0d0,cc1=-1.0d0/6.0d0,cc2=7.0d0/6.0d0
      integer::i,j,k,i0,i1,i2,i3,k0,k1,k2,k3
      integer::iflag

      dr1= 0.50d0/dr
      dz1= 0.50d0/dz
      dphi1= 0.50d0/dphi
!      cc0= 4.0d0/3.0d0
!      cc1=-1.0d0/6.0d0
!      cc2= 7.0d0/6.0d0


! derivative in r direction

#ifdef AMDGPU
!$omp target teams distribute parallel do collapse(3) private(k,j,i,i0,i1,i2,i3)
#else
!$omp parallel do private(i0,i1,i2,i3)
#endif
      do k = 3, lphi-2
      do j = 1, lz
      do i = 3, lr-2
        i3 = i + 2
        i2 = i + 1
        i1 = i - 1
        i0 = i - 2

        b_phi(i,j,k) = -(cc0*cf(i2,j,k)*sc(i2,j,k)*a_z(i2,j,k) &
                       + cc1*cf(i3,j,k)*sc(i3,j,k)*a_z(i3,j,k) &
                       - cc0*cf(i1,j,k)*sc(i1,j,k)*a_z(i1,j,k) &
                       - cc1*cf(i0,j,k)*sc(i0,j,k)*a_z(i0,j,k) )*dr1
        b_z(i,j,k) = (cc0*grr(i2,j,k)*cf(i2,j,k)*sc(i2,j,k)*a_phi(i2,j,k) &
                     +cc1*grr(i3,j,k)*cf(i3,j,k)*sc(i3,j,k)*a_phi(i3,j,k) &
                     -cc0*grr(i1,j,k)*cf(i1,j,k)*sc(i1,j,k)*a_phi(i1,j,k) &
                     -cc1*grr(i0,j,k)*cf(i0,j,k)*sc(i0,j,k)*a_phi(i0,j,k) &
                     )*dr1/grr(i,j,k)
      end do
      end do
      end do

#ifdef AMDGPU
!$omp target teams distribute parallel do collapse(2) private(k,j)
#else
!$omp parallel do
#endif
      do k = 3, lphi-2
      do j = 1, lz
        b_phi(1,j,k) = -(cf(3,j,k)*sc(3,j,k)*a_z(3,j,k)*cc1 &
                       + cf(2,j,k)*sc(2,j,k)*a_z(2,j,k)*cc0 &
                       - cf(1,j,k)*sc(1,j,k)*a_z(1,j,k)*cc2)*dr1
        b_z(1,j,k) = (grr(3,j,k)*cf(3,j,k)*sc(3,j,k)*a_phi(3,j,k)*cc1 &
                     +grr(2,j,k)*cf(2,j,k)*sc(2,j,k)*a_phi(2,j,k)*cc0 &
                     -grr(1,j,k)*cf(1,j,k)*sc(1,j,k)*a_phi(1,j,k)*cc2 &
                     )*dr1/grr(1,j,k)

        b_phi(2,j,k) = -(cf(4,j,k)*sc(4,j,k)*a_z(4,j,k)*cc1 &
                       + cf(3,j,k)*sc(3,j,k)*a_z(3,j,k)*cc0 &
                       + cf(2,j,k)*sc(2,j,k)*a_z(2,j,k)*cc1 &
                       - cf(1,j,k)*sc(1,j,k)*a_z(1,j,k))*dr1
        b_z(2,j,k) = (grr(4,j,k)*cf(4,j,k)*sc(4,j,k)*a_phi(4,j,k)*cc1 &
                    + grr(3,j,k)*cf(3,j,k)*sc(3,j,k)*a_phi(3,j,k)*cc0 &
                    + grr(2,j,k)*cf(2,j,k)*sc(2,j,k)*a_phi(2,j,k)*cc1 &
                    - grr(1,j,k)*cf(1,j,k)*sc(1,j,k)*a_phi(1,j,k) &
                     )*dr1/grr(2,j,k)

        b_phi(lr-1,j,k) = -(cf(lr,j,k)*sc(lr,j,k)*a_z(lr,j,k) &
                          - cf(lr-1,j,k)*sc(lr-1,j,k)*a_z(lr-1,j,k)*cc1 &
                          - cf(lr-2,j,k)*sc(lr-2,j,k)*a_z(lr-2,j,k)*cc0 &
                          - cf(lr-3,j,k)*sc(lr-3,j,k)*a_z(lr-3,j,k)*cc1 )*dr1
        b_z(lr-1,j,k) = (grr(lr,j,k)  *cf(lr,j,k)*sc(lr,j,k)*a_phi(lr,j,k) &
                        -grr(lr-1,j,k)*cf(lr-1,j,k)*sc(lr-1,j,k)*a_phi(lr-1,j,k)*cc1 &
                        -grr(lr-2,j,k)*cf(lr-2,j,k)*sc(lr-2,j,k)*a_phi(lr-2,j,k)*cc0 &
                        -grr(lr-3,j,k)*cf(lr-3,j,k)*sc(lr-3,j,k)*a_phi(lr-3,j,k)*cc1 &
                        )*dr1/grr(lr-1,j,k)

        b_phi(lr,j,k) = -(cf(lr,j,k)*sc(lr,j,k)*a_z(lr,j,k)*cc2 &
                         -cf(lr-1,j,k)*sc(lr-1,j,k)*a_z(lr-1,j,k)*cc0 &
                         -cf(lr-2,j,k)*sc(lr-2,j,k)*a_z(lr-2,j,k)*cc1)*dr1
        b_z(lr,j,k) = (grr(lr,j,k)  *cf(lr,j,k)*sc(lr,j,k)*a_phi(lr,j,k)*cc2 &
                      -grr(lr-1,j,k)*cf(lr-1,j,k)*sc(lr-1,j,k)*a_phi(lr-1,j,k)*cc0 &
                      -grr(lr-2,j,k)*cf(lr-2,j,k)*sc(lr-2,j,k)*a_phi(lr-2,j,k)*cc1 &
                       )*dr1/grr(lr,j,k)
      end do
      end do

! derivative in z direction

#ifdef AMDGPU
!$omp target teams distribute parallel do collapse(2) private(k,i,i0,i1,i2,i3)
#else
!$omp parallel do private(i0,i1,i2,i3)
#endif
      do k = 3, lphi-2
      do i = 2*lr+1, lrz-2*lr
        i3 = i + 2*lr
        i2 = i + lr
        i1 = i - lr
        i0 = i - 2*lr

        b_r(i,1,k) =-(cf(i2,1,k)*sc(i2,1,k)*a_phi(i2,1,k)*cc0 &
                    + cf(i3,1,k)*sc(i3,1,k)*a_phi(i3,1,k)*cc1 &
                    - cf(i1,1,k)*sc(i1,1,k)*a_phi(i1,1,k)*cc0 &
                    - cf(i0,1,k)*sc(i0,1,k)*a_phi(i0,1,k)*cc1 )*dz1
        b_phi(i,1,k) = (cf(i2,1,k)*sc(i2,1,k)*a_r(i2,1,k)*cc0 &
                      + cf(i3,1,k)*sc(i3,1,k)*a_r(i3,1,k)*cc1 &
                      - cf(i1,1,k)*sc(i1,1,k)*a_r(i1,1,k)*cc0 &
                      - cf(i0,1,k)*sc(i0,1,k)*a_r(i0,1,k)*cc1 )*dz1 &
                     + b_phi(i,1,k)
      end do
      end do

#ifdef AMDGPU
!$omp target teams distribute parallel do collapse(2) private(k,i)
#else
!$omp parallel do
#endif
      do k = 3, lphi-2
      do i = 1, lr
        b_r(i,1,k)   =-(cf(i,3,k)*sc(i,3,k)*a_phi(i,3,k)*cc1 &
                      + cf(i,2,k)*sc(i,2,k)*a_phi(i,2,k)*cc0 &
                      - cf(i,1,k)*sc(i,1,k)*a_phi(i,1,k)*cc2 )*dz1
        b_phi(i,1,k) = (cf(i,3,k)*sc(i,3,k)*a_r(i,3,k)*cc1 &
                      + cf(i,2,k)*sc(i,2,k)*a_r(i,2,k)*cc0 &
                      - cf(i,1,k)*sc(i,1,k)*a_r(i,1,k)*cc2   )*dz1 &
                     + b_phi(i,1,k)

        b_r(i,2,k)   =-(cf(i,4,k)*sc(i,4,k)*a_phi(i,4,k)*cc1 &
                      + cf(i,3,k)*sc(i,3,k)*a_phi(i,3,k)*cc0 &
                      + cf(i,2,k)*sc(i,2,k)*a_phi(i,2,k)*cc1 &
                      - cf(i,1,k)*sc(i,1,k)*a_phi(i,1,k) )*dz1
        b_phi(i,2,k) = (cf(i,4,k)*sc(i,4,k)*a_r(i,4,k)*cc1 &
                      + cf(i,3,k)*sc(i,3,k)*a_r(i,3,k)*cc0 &
                      + cf(i,2,k)*sc(i,2,k)*a_r(i,2,k)*cc1 &
                      - cf(i,1,k)*sc(i,1,k)*a_r(i,1,k)   )*dz1 &
                     + b_phi(i,2,k)

        b_r(i,lz-1,k)   =-(cf(i,lz,k)*sc(i,lz,k)*a_phi(i,lz,k) &
                         - cf(i,lz-1,k)*sc(i,lz-1,k)*a_phi(i,lz-1,k)*cc1 &
                         - cf(i,lz-2,k)*sc(i,lz-2,k)*a_phi(i,lz-2,k)*cc0 &
                         - cf(i,lz-3,k)*sc(i,lz-3,k)*a_phi(i,lz-3,k)*cc1 )*dz1
        b_phi(i,lz-1,k) = (cf(i,lz,k)*sc(i,lz,k)*a_r(i,lz,k) &
                         - cf(i,lz-1,k)*sc(i,lz-1,k)*a_r(i,lz-1,k)*cc1 &
                         - cf(i,lz-2,k)*sc(i,lz-2,k)*a_r(i,lz-2,k)*cc0 &
                         - cf(i,lz-3,k)*sc(i,lz-3,k)*a_r(i,lz-3,k)*cc1 )*dz1 &
                        + b_phi(i,lz-1,k)

        b_r(i,lz,k)   =-(cf(i,lz,k)*sc(i,lz,k)*a_phi(i,lz,k)*cc2 &
                       - cf(i,lz-1,k)*sc(i,lz-1,k)*a_phi(i,lz-1,k)*cc0 &
                       - cf(i,lz-2,k)*sc(i,lz-2,k)*a_phi(i,lz-2,k)*cc1 )*dz1
        b_phi(i,lz,k) = (cf(i,lz,k)*sc(i,lz,k)*a_r(i,lz,k)*cc2 &
                       - cf(i,lz-1,k)*sc(i,lz-1,k)*a_r(i,lz-1,k)*cc0 &
                       - cf(i,lz-2,k)*sc(i,lz-2,k)*a_r(i,lz-2,k)*cc1 )*dz1 &
                      + b_phi(i,lz,k)
      end do
      end do

! derivative in phi direction

#ifdef AMDGPU
!$omp target teams distribute parallel do collapse(2) private(k,i,k0,k1,k2,k3)
#else
!$omp parallel do private(k0,k1,k2,k3)
#endif
      do k = 3, lphi-2
      do i = 1, lrz
        k3 = k + 2
        k2 = k + 1
        k1 = k - 1
        k0 = k - 2

        b_r(i,1,k) = (cf(i,1,k2)*sc(i,1,k2)*a_z(i,1,k2)*cc0 &
                    + cf(i,1,k3)*sc(i,1,k3)*a_z(i,1,k3)*cc1 &
                    - cf(i,1,k1)*sc(i,1,k1)*a_z(i,1,k1)*cc0 &
                    - cf(i,1,k0)*sc(i,1,k0)*a_z(i,1,k0)*cc1 )*dphi1/grr(i,1,k) &
                   + b_r(i,1,k)
        b_z(i,1,k) =-(cf(i,1,k2)*sc(i,1,k2)*a_r(i,1,k2)*cc0 &
                    + cf(i,1,k3)*sc(i,1,k3)*a_r(i,1,k3)*cc1 &
                    - cf(i,1,k1)*sc(i,1,k1)*a_r(i,1,k1)*cc0 &
                    - cf(i,1,k0)*sc(i,1,k0)*a_r(i,1,k0)*cc1 )*dphi1/grr(i,1,k) &
                   + b_z(i,1,k)
      end do
      end do

      if(iflag.eq.1)then
        call periodic_field_mlt3(b_r,b_z,b_phi) !2012-04-13
!        call periodic_field(b_r)
!        call periodic_field(b_z)
!        call periodic_field(b_phi)
      end if

end
!--------------------------------------------------------------------
! Laplacian with conservation property, modified on 2013-02-17
subroutine laplacian5(iflag,a,a0,lpl3_a,coef)
!--------------------------------------------------------------------
      use grid
      implicit none

      real(8)::a(lr,lz,lphi),a0(lr,lz,lphi)
      real(8)::lpl3_a(lr,lz,lphi),coef(lr,lz,lphi)
      real(8)::dri2,dzi2,dphii2,r0,r1,z0,z1,phi0,phi1
      integer::i,j,k,i0,i1,j0,j1,k0,k1
      integer::iflag

      dri2 = 1.0d0/dr**2
      dzi2 = 1.0d0/dz**2
      dphii2 = 1.0d0/dphi**2


! laplacian
! r-direction      

#ifdef AMDGPU  
!$omp target teams distribute parallel do collapse(3) private(k,j,i,i0,i1,r0,r1)
#else
!$omp parallel do private(i0,i1,r0,r1)
#endif
      do k = 3, lphi-2
        do j = 1, lz
        do i = 2, lr-1
          i0 = i - 1
          i1 = i + 1

          r1 =(grr(i,j,k) + 0.50d0*dr)*(coef(i1,j,k)+coef(i,j,k) )*0.50d0
          r0 =(grr(i,j,k) - 0.50d0*dr)*(coef(i0,j,k)+coef(i,j,k) )*0.50d0

          lpl3_a(i,j,k) = dri2*(r0*(a(i0,j,k) - a0(i0,j,k) - a(i,j,k) + a0(i,j,k) ) &
                              + r1*(a(i1,j,k) - a0(i1,j,k) - a(i,j,k) + a0(i,j,k) ) &
                                )/grr(i,j,k)
        end do
        end do
        end do
     
#ifdef AMDGPU  
!$omp target teams distribute parallel do collapse(2) private(k,j,i,i0,i1,r0,r1)
#else
!$omp parallel do private(i,i0,i1,r0,r1)
#endif
        do k = 3, lphi-2
        do j = 1, lz
          i = 1
          i1 = i + 1

          r1 =(grr(i,j,k) + 0.50d0*dr)*(coef(i1,j,k)+coef(i,j,k) )*0.50d0
          lpl3_a(i,j,k) = dri2*r1/grr(i,j,k) &
                              *(a(i1,j,k) - a0(i1,j,k) - a(i,j,k) + a0(i,j,k) )
          
          i = lr
          i0 = i - 1

          r0 =(grr(i,j,k) - 0.50d0*dr)*(coef(i0,j,k)+coef(i,j,k) )*0.50d0
          lpl3_a(i,j,k) = dri2*r0/grr(i,j,k) &
                              *(a(i0,j,k) - a0(i0,j,k) - a(i,j,k) + a0(i,j,k) )
        end do
        end do

! z-direction      
        
#ifdef AMDGPU  
!$omp target teams distribute parallel do collapse(3) private(k,j,i,j0,j1,z0,z1)
#else
!$omp parallel do private(j0,j1,z0,z1)
#endif
        do k = 3, lphi-2
        do j = 2, lz-1
        do i = 1, lr
          j0 = j - 1
          j1 = j + 1

          z1 = (coef(i,j1,k)+coef(i,j,k) )*0.50d0
          z0 = (coef(i,j0,k)+coef(i,j,k) )*0.50d0

          lpl3_a(i,j,k) = lpl3_a(i,j,k) &
                        + dzi2*( z0*(a(i,j0,k) - a0(i,j0,k) - a(i,j,k) + a0(i,j,k) ) &
                               + z1*(a(i,j1,k) - a0(i,j1,k) - a(i,j,k) + a0(i,j,k) ) &
                               )
        end do
        end do
        end do
     
#ifdef AMDGPU  
!$omp target teams distribute parallel do collapse(2) private(k,j,i,j0,j1,z0,z1)
#else
!$omp parallel do private(j,j0,j1,z0,z1)
#endif
        do k = 3, lphi-2
        do i = 1, lr
          j = 1 
          j1 = j + 1

          z1 = (coef(i,j1,k)+coef(i,j,k) )*0.50d0
          lpl3_a(i,j,k) = lpl3_a(i,j,k) &
                        + dzi2*(        &
                               + z1*(a(i,j1,k) - a0(i,j1,k) - a(i,j,k) + a0(i,j,k) ) &
                               )
          j = lz  
          j0 = j - 1

          z0 = (coef(i,j0,k)+coef(i,j,k) )*0.50d0
          lpl3_a(i,j,k) = lpl3_a(i,j,k) &
                        + dzi2*( z0*(a(i,j0,k) - a0(i,j0,k) - a(i,j,k) + a0(i,j,k) ) &
                               )
        end do
        end do

! phi-direction      

#ifdef AMDGPU  
!$omp target teams distribute parallel do collapse(2) private(k,i,k0,k1,phi0,phi1)
#else
!$omp parallel do private(k0,k1,phi0,phi1)
#endif
        do k = 3, lphi-2
        do i = 1, lrz
          k0 = k - 1
          k1 = k + 1

          phi1 = (coef(i,1,k1)+coef(i,1,k) )*0.50d0
          phi0 = (coef(i,1,k0)+coef(i,1,k) )*0.50d0

          lpl3_a(i,1,k) = lpl3_a(i,1,k) &
                        + dphii2*( phi0*(a(i,1,k0) - a0(i,1,k0) - a(i,1,k) + a0(i,1,k) ) &
                                 + phi1*(a(i,1,k1) - a0(i,1,k1) - a(i,1,k) + a0(i,1,k) ) &
                                 )/grr(i,1,k)**2
        end do
        end do


      if(iflag.eq.1)then
        call periodic_field(lpl3_a)
      end if

end
!--------------------------------------------------------------------
subroutine wall_clock(wtime)
! wall clock
!--------------------------------------------------------------------
      use mpiset
      real(8)::wtime

      call mpi_barrier(mpi_comm_world,mpi_err)
      wtime = mpi_wtime()

end
