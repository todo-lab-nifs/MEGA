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
program profile
! profile analysis program
!--------------------------------------------------
      use harmonics_array
      implicit none
      real(8)::fm(0:mpol,lpsi,2)
      character(6)::job_name
      character(3)::job_no
      character(50)::file_in
      character(50)::file_out
      integer::kstep,n_toroidal,kst,l,m,n,l_max,l_max2,m_max,m1,n_max,i
      integer::n_fundamental
      real(8)::pi,twopi,t,amp_max,amp,alpha1,alpha2,alpha,theta0,phi0,sgn

      file_in ='/data/usr1/itp023/itp023_001.harmonics'
      file_out ='itp023_001_kstep=1000000_n=+01-vrad.txt'

      print *, 'input job name like itp023'
      read *, job_name

      print *, 'input job number like 001'
      read *, job_no

      print *, 'input step number to anlyze like 20000'
      read *, kstep

      print *, 'input toroidal mode number like -11 or 5'
      read *, n_toroidal

      print *, 'input 2*pi / toroidal length'
      read *, n_fundamental

      write(file_in(12:17),'(a)')job_name
      write(file_in(19:24),'(a)')job_name
      write(file_in(26:28),'(a)')job_no
      write(file_out(1:6),'(a)')job_name
      write(file_out(8:10),'(a)')job_no
      write(file_out(18:24),'(i7.7)')kstep

      if(n_toroidal.ge.0)then
        write(file_out(29:30),'(i2.2)')n_toroidal
      else
        write(file_out(28:30),'(i3.2)')n_toroidal
      end if     
      pi = 4.0d0*atan(1.0d0)
      twopi = pi*2.d0

!      open(10,file=file_in,form='unformatted',convert='BIG_ENDIAN')
      open(10,file=file_in,form='unformatted')

      kst = -1
      do while (kst.lt.kstep)
        read(10,end=9999)kst,t,r_psi,gpsi_nrm,q_psi &
                ,vrad,vtheta,vphi,brad,btheta,bphi &
                ,erad,etheta,ephi,prs,rho &
!                ,dns_e,mom_e,ppara_e,pperp_e &
                ,dns_i,mom_i,ppara_i,pperp_i &
                ,dns_a,mom_a,ppara_a,pperp_a
!                ,dns_a,mom_a,ppara_a,pperp_a,qpara_a,qperp_a
        print *, 'kstep=',kst
      end do
 9999 continue
      close(10)
      kstep = kst
      print *, 'kstep=',kstep

      open(7,file=file_out)

      write(7,7000)'r/a','gpsi_sqrt','q','m=0cos','m=0sin' &
                  ,'m=1cos','m=1sin','m=2cos','m=2sin' &
                  ,'m=3cos','m=3sin','m=4cos','m=4sin' &
                  ,'m=5cos','m=5sin','m=6cos','m=6sin' &
                  ,'m=7cos','m=7sin','m=8cos','m=8sin' &
                  ,'m=9cos','m=9sin','m=10cos','m=10sin' &
                  ,'m=11cos','m=11sin','m=12cos','m=12sin' &
                  ,'m=13cos','m=13sin','m=14cos','m=14sin' &
                  ,'m=15cos','m=15sin','m=16cos','m=16sin' &
                  ,'m=17cos','m=17sin','m=18cos','m=18sin' &
                  ,'m=19cos','m=19sin','m=20cos','m=20sin' &
                  ,'m=21cos','m=21sin','m=22cos','m=22sin' &
                  ,'m=23cos','m=23sin','m=24cos','m=24sin' &
                  ,'m=25cos','m=25sin','m=26cos','m=26sin' &
                  ,'m=27cos','m=27sin','m=28cos','m=28sin' &
                  ,'m=29cos','m=29sin','m=30cos','m=30sin' &
                  ,'m=31cos','m=31sin','m=32cos','m=32sin' &
                  ,'m=33cos','m=33sin','m=34cos','m=34sin' &
                  ,'m=35cos','m=35sin','m=36cos','m=36sin' &
                  ,'m=37cos','m=37sin','m=38cos','m=38sin' &
                  ,'m=39cos','m=39sin','m=40cos','m=40sin' &
                  ,'m=41cos','m=41sin','m=42cos','m=42sin' &
                  ,'m=43cos','m=43sin','m=44cos','m=44sin' &
                  ,'m=45cos','m=45sin','m=46cos','m=46sin' &
                  ,'m=47cos','m=47sin','m=48cos','m=48sin' &
                  ,'m=49cos','m=49sin','m=50cos','m=50sin' &
                  ,'m=51cos','m=51sin','m=52cos','m=52sin' &
                  ,'m=53cos','m=53sin','m=54cos','m=54sin' &
                  ,'m=55cos','m=55sin','m=56cos','m=56sin' &
                  ,'m=57cos','m=57sin','m=58cos','m=58sin' &
                  ,'m=59cos','m=59sin','m=60cos','m=60sin' &
                  ,'m=61cos','m=61sin','m=62cos','m=62sin' &
                  ,'m=63cos','m=63sin','m=64cos','m=64sin'
 7000 format(150(2x,a))

      amp_max = 0.0d0
      l_max = 1
      m_max = 0
      n_max = -ntor

      do l = 1, lpsi
      do n = -ntor, ntor
      do m = 0, mpol 
        amp = sqrt( vrad(m,n,l,1)**2 &
                   +vrad(m,n,l,2)**2 &
                   )
        if(amp.gt.amp_max)then
          m_max = m
          n_max = n
          l_max = l
          amp_max = amp
        end if
      end do
      end do
      end do

      print *, 'm_max,n_max,l_max=',m_max,n_max*n_fundamental,l_max

      amp_max = 0.0d0
      l_max = 1
      m_max = 0
      n = n_toroidal/n_fundamental

!      do l = 1, lpsi
      do l = 3, lpsi
      do m = 0, mpol 
        amp = sqrt( vrad(m,n,l,1)**2 &
                   +vrad(m,n,l,2)**2 &
                   )
        if(amp.gt.amp_max)then
          m_max = m
          l_max = l
          amp_max = amp
        end if
      end do
      end do

      print *, 'm_max,n_toroidal,l_max=',m_max,n_toroidal,l_max

      alpha1 = atan2(vrad(m_max,n,l_max,2), &
                     vrad(m_max,n,l_max,1)  &
                    )

      amp_max = 0.0d0
      l_max2 = 1
      m1 = m_max - 1
      n = n_toroidal/n_fundamental

      do l = 1, lpsi
        amp = sqrt( vrad(m1,n,l,1)**2 &
                   +vrad(m1,n,l,2)**2 &
                   )
        if(amp.gt.amp_max)then
          l_max2 = l
          amp_max = amp
        end if
      end do

      alpha2 = atan2(vrad(m1,n,l_max2,2), &
                     vrad(m1,n,l_max2,1)  &
                    )

! alpha = m*theta0 + n*phi0
! change phase by alpha
! x = a'*cos(n*phi+m*theta) + b'*sin(n*phi+m*theta)
!   = a*cos(n*phi+m*theta - alpha) + b*sin(n*phi+m*theta - alpha)
!   = a*( cos()*cos(alpha) + sin()*sin(alpha) )
!    +b*( sin()*cos(alpha) - cos()*sin(alpha) )
! a' = a*cos(alpha) - b*sin(alpha)
! b' = a*sin(alpha) + b*cos(alpha)
! if tan(alpha)= b'/a', b=0

! cos(m*(theta-theta0) + n*(phi-phi0))
! m*theta0 + n*phi0 = alpha1
! m1*theta0 + n*phi0 = alpha2
! theta0 =(alpha2 - alpha1)/(m1-m)

!      theta0 =(alpha2 - alpha1)/dble(m1 - m_max)
      theta0 = 0.0d0
      if(n_toroidal.ne.0)then
        phi0 =(-dble(m_max)*theta0 + alpha1)/dble(n_toroidal)
      else
        phi0 = 0.0d0
      end if

      print *, 'alpha1,alpha2=',alpha1,alpha2
      print *, 'theta0,phi0=',theta0,phi0

      n = n_toroidal/n_fundamental
      alpha = dble(m_max)*theta0 + dble(n_toroidal)*phi0
      sgn = sign(1.0d0, vrad(m_max,n,l_max,1)*cos(alpha) &
                       +vrad(m_max,n,l_max,2)*sin(alpha) )

      do l = 1, lpsi
        gpsi_sqrt(l) = sqrt(max(0.0d0,gpsi_nrm(l)))
      do m = 0, mpol
        alpha = dble(m)*theta0 + dble(n_toroidal)*phi0
        n = n_toroidal/n_fundamental

        fm(m,l,1) =(vrad(m,n,l,1)*cos(alpha) &
                   +vrad(m,n,l,2)*sin(alpha) )*sgn
!                   +vrad(m,n,l,2)*sin(alpha) )*gpsi_sqrt(l)*sgn
        fm(m,l,2) =(-vrad(m,n,l,1)*sin(alpha) &
                    +vrad(m,n,l,2)*cos(alpha) )*sgn
!                    +vrad(m,n,l,2)*cos(alpha) )*gpsi_sqrt(l)*sgn
      end do 
      end do 

      do l = 1, lpsi
        write(7,7100)r_psi(l),gpsi_sqrt(l),q_psi(l),((fm(m,l,i),i=1,2),m=0,mpol)
      end do
 7100 format(150(1PE14.5) )

end
