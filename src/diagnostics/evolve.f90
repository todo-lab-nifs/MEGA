module harmonics_array
  integer, parameter::lpsi=201,mpol=64,ntor=2
  real(8)::r_psi(lpsi),q_psi(lpsi)
  real(8)::gpsi_nrm(lpsi),gpsi_sqrt(lpsi)
  real(8)::vrad(0:mpol,-ntor:ntor,lpsi,2)
  real(8)::vtheta(0:mpol,-ntor:ntor,lpsi,2)
  real(8)::vphi(0:mpol,-ntor:ntor,lpsi,2)
  real(8)::brad(0:mpol,-ntor:ntor,lpsi,2)
  real(8)::btheta(0:mpol,-ntor:ntor,lpsi,2)
  real(8)::bphi(0:mpol,-ntor:ntor,lpsi,2)
  real(8)::erad(0:mpol,-ntor:ntor,lpsi,2)
  real(8)::etheta(0:mpol,-ntor:ntor,lpsi,2)
  real(8)::ephi(0:mpol,-ntor:ntor,lpsi,2)
  real(8)::rho(0:mpol,-ntor:ntor,lpsi,2)
  real(8)::prs(0:mpol,-ntor:ntor,lpsi,2)

  real(8)::dns_e(0:mpol,-ntor:ntor,lpsi,2)
  real(8)::mom_e(0:mpol,-ntor:ntor,lpsi,2)
  real(8)::ppara_e(0:mpol,-ntor:ntor,lpsi,2)
  real(8)::pperp_e(0:mpol,-ntor:ntor,lpsi,2)
  real(8)::dns_i(0:mpol,-ntor:ntor,lpsi,2)
  real(8)::mom_i(0:mpol,-ntor:ntor,lpsi,2)
  real(8)::ppara_i(0:mpol,-ntor:ntor,lpsi,2)
  real(8)::pperp_i(0:mpol,-ntor:ntor,lpsi,2)
  real(8)::dns_a(0:mpol,-ntor:ntor,lpsi,2)
  real(8)::mom_a(0:mpol,-ntor:ntor,lpsi,2)
  real(8)::ppara_a(0:mpol,-ntor:ntor,lpsi,2)
  real(8)::pperp_a(0:mpol,-ntor:ntor,lpsi,2)
!  real(8)::qpara_a(0:mpol,-ntor:ntor,lpsi,2)
!  real(8)::qperp_a(0:mpol,-ntor:ntor,lpsi,2)
end module

!--------------------------------------------------
program evolve
! evolution analysis program
!--------------------------------------------------
      use harmonics_array
      implicit none
      integer, parameter::ltime=5000
      real(8)::time(ltime),amplitude(ltime),phase(ltime),amp_cos(ltime),amp_sin(ltime)
      character(6)::job_no
      character(50)::file_in
      character(50)::file_out
      integer::no_serial,no_starting,no_ending
      integer::kst,kst0,ievol,i,n
      integer::m_poloidal,n_toroidal,l_peak,n_fundamental
      real(8)::pi,twopi,t,wa,omega

      file_in ='/data/usr1/itp001/itp001_001.harmonics'
      file_out ='itp001_evol_m=32_n=+19_l=117-vrad.txt'

      write(6,*)'input jobname in a form like itp001'
      read(5,'(a)')job_no
      write(6,*)'input starting jobno in a form like 1'
      read(5,*)no_starting
      write(6,*)'input ending jobno in a form like 99'
      read(5,*)no_ending

      write(file_in(12:17),'(a)')job_no
      write(file_in(19:24),'(a)')job_no
      write(file_out(1:6),'(a)')job_no


      print *, 'input poloidal mode number (>=0) like 15'
      read *, m_poloidal

      print *, 'input toroidal mode number like -11 or 5'
      read *, n_toroidal

      print *, 'input 2*pi / toroidal length'
      read *, n_fundamental

      print *, 'input peak location (>=1) like 49'
      read *, l_peak

      write(file_out(15:16),'(i2.2)')m_poloidal
      write(file_out(26:28),'(i3.3)')l_peak
      if(n_toroidal.ge.0)then
        write(file_out(21:22),'(i2.2)')n_toroidal
      else
        write(file_out(20:22),'(i3.2)')n_toroidal
      end if     

      pi = 4.0d0*atan(1.0d0)
      twopi = pi*2.d0
!      wa = 1.0d0/(1.6006d0*1.22d1) !Raxis=1.6006 in eql066, minor_r=12.2 in mega2025.f90
      wa = 1.0d0/(3.30428d0*1.6d1) !Raxis=3.30428 in eql045, minor_r=16 in mega2025_v6.f90

      ievol = 0
      kst0 = -1
loop_file: do no_serial = no_starting, no_ending
        write(file_in(26:28),'(i3.3)')no_serial

!        open(10,file=file_in,form='unformatted',convert='BIG_ENDIAN')
        open(10,file=file_in,form='unformatted')

      do while(ievol.le.ltime)
        read(10,end=9999)kst,t,r_psi,gpsi_nrm,q_psi &
                ,vrad,vtheta,vphi,brad,btheta,bphi &
                ,erad,etheta,ephi,prs,rho &
!                ,dns_e,mom_e,ppara_e,pperp_e &
                ,dns_i,mom_i,ppara_i,pperp_i &
                ,dns_a,mom_a,ppara_a,pperp_a
!                ,dns_a,mom_a,ppara_a,pperp_a,qpara_a,qperp_a
        print *, 'kstep=',kst

        if(kst.gt.kst0)then
          kst0 = kst
          ievol = ievol + 1
          time(ievol) = t*wa
          n = n_toroidal/n_fundamental

          amplitude(ievol) = sqrt(vrad(m_poloidal,n,l_peak,1)**2 &
                                 +vrad(m_poloidal,n,l_peak,2)**2 &
                                 )
          phase(ievol) = atan2(vrad(m_poloidal,n,l_peak,2) &
                              ,max(1.0d-20, abs(vrad(m_poloidal,n,l_peak,1))) &
                               *sign(1.0d0, vrad(m_poloidal,n,l_peak,1) ) &
                              )
          amp_cos(ievol) = vrad(m_poloidal,n,l_peak,1)
          amp_sin(ievol) = vrad(m_poloidal,n,l_peak,2)
        end if
      end do
 9999 continue
      close(10)

      end do loop_file

      open(7,file=file_out)

      write(7,'(5(2x,a))')'time','amplitude','omega','amp_cos','amp_sin'
      do i = 2, ievol-1
        omega = phase(i) - phase(i-1)
        if(abs(omega).gt.pi)then
          omega = -twopi*sign(1.0d0, omega) + omega
        end if
        omega = omega/(time(i) - time(i-1) )
        write(7,'(5(1pe14.5))')time(i),amplitude(i),omega,amp_cos(i),amp_sin(i)
      end do

end
