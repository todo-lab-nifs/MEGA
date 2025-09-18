!---------------------------------------------------------------------------------
! modules: aurora, fx100, fugaku, xeon6900p, mi300-A
! contain optimized subroutines: push, density, moments, density_gyro, emf_gyro, moments_gyro
! 2020-07-11
! 2020-11-22: module fugaku is created
! 2025-07-11: module xeon6900p is created for subsystem A of new Plasma Simulator
! 2025-09-03: module amdgpu is created for subsystem B of new Plasma Simulator
!---------------------------------------------------------------------------------

!----------------------------------------------------------------------------
module aurora
contains
!--------------------------------------------------------------------
subroutine push(marker_num,amassp,ep &
               ,type,temp,valpha,deltav,clambda0,dclambda & !2012-06-17
               ,gc,dgc,v &
               ,cf_pphi,pphi_min,pphi_max &
               ,flp,ispecies,igyro)
! type=0: maxwellian, type=1: slowing down, type=2: beam
! type=3: beam with finite pitch angle width !2012-06-17
! igyro=0: w/o FLR, igyro=1: w/ FLR
! modified for NEC SX-Aurora TSUBASA 2020-07-10
!   gc_x(:,n) = 1:r, 2:z, 3:phi, 4:p, 5:mu, 6:weight, 7:active, 8:p_phi, 9:fff, 10:fnrml
!--------------------------------------------------------------------
      use mpiset
      use field
      use grid
!      use equi_sol, only:raxis
      use particle, only:nu_krook !2024-04-23
      implicit none

      integer::marker_num,type
      integer::ispecies,igyro
      real(8)::amassp,ep
      real(8)::gc(ngc2,marker_num)
      real(8)::dgc(ngc1,marker_num)
      real(8)::v(marker_num)
      real(8)::cf_pphi(0:lcfpphi),pphi_min,pphi_max !2016-02-04
      real(8)::temp,valpha,deltav
      real(8)::dr1,dz1,dphi1,ep1,amsp1
      integer::n,ia,ia1,ja,ja1,ka,ka1,i
      real(8)::ar,ar1,az,az1,aphi,aphi1
!      real(8)::aaa1,aaa2,aaa3,aaa4,aaa5,aaa6,aaa7,aaa8

!2015-06-04s
!      integer::ijk_a(3,marker_each)
!      real(8)::aaa(8,marker_each)
!2015-06-04e

      real(8)::flp(nflp,marker_num),detotal
      real(8)::b2,babse,babs0e,b21,bre,bze,bphie
      real(8)::b1,br1a,bz1a,bphi1a,b10,br10,bz10,bphi10,dbpar,dbpr,dbpz,dbpphi
      real(8)::rhopar,orbpr,orbpz,orbpphi,denom1
      real(8)::wvnmlr,wvnmlz,wvnmlphi
      real(8)::dgrdbr,dgrdbz,dgrdbphi,dppar,dppar1,dppar2 !2012-08-31
      real(8)::pvpar,w1nmlr,w1nmlz,w1nmlphi
      real(8)::psinrm,pactive,prof,dprofdpsi,energy0,pphi_n,rminor,bmax
      real(8)::dpphi_n,vpara,vpara0,dvpara,sigma
      real(8)::vlcpart,dvlcpart,vnrm,sd_factor
      real(8)::sqr_pi1
      real(8)::nudt,coef,vrat !2016-01-09
      real(8)::nu_krook_dt !2024-04-23

      integer::kr,kz,kphi
      real(8)::ma_mi_r
      integer::kfl_start,kfl

      real(8)::pt_factor,energy,clambda0,dclambda,clambda(marker_each) !2012-06-17
      real(8)::dwpsi(marker_each),dwenrc(marker_each) !2015-06-23
!      real(8)::psip(marker_each)
!      real(8)::dpsi_dr(marker_each),dpsi_dz(marker_each),dpsi_dphi(marker_each)
      real(8)::bphi1e(marker_each) !2016-01-09

      integer,parameter::nblkp=1024 !NEC
!2021-01-11s
      integer::ijk_a(nblkp,3)
      real(8)::aaa(nblkp,8)
      integer::nn,nsta,nend,ivect
!2021-01-11e


! time derivative of particle position and velocity

      dr1 = 1.0d0/dr
      dz1 = 1.0d0/dz
      dphi1 = 1.0d0/dphi
      ep1 = 1.0d0/ep
      amsp1 = 1.0d0/amassp
      sqr_pi1=1.0d0/sqrt(pi)

      nu_krook_dt = nu_krook * dt !2024-04-23

      kr   = 1 - kr_offset(my_rank)
      kz   = 1 - kz_offset(my_rank)
      kphi = 1 + lphi_shd - kphi_offset(my_rank)
      ma_mi_r = major_r-minor_r

! igyro=0: w/o FLR, igyro/=0: flp(1:7,:) is given in subr. extract_em
      if(igyro.eq.0)then
        kfl_start = 1
      else
        kfl_start = nflp_gyro + 1
      end if

!      call ftrace_region_begin('push0')

!$omp parallel
      if(igyro.eq.0)then !w/o FLR
!$omp workshare
        fld(1,:,:,:) = er(:,:,:)
        fld(2,:,:,:) = ez(:,:,:)
        fld(3,:,:,:) = ephi(:,:,:)
        fld(4,:,:,:) = epara(:,:,:)
        fld(5,:,:,:) = br(:,:,:) - br0(:,:,:)
        fld(6,:,:,:) = bz(:,:,:) - bz0(:,:,:)
        fld(7,:,:,:) = bphi(:,:,:) - bphi0(:,:,:)
!$omp end workshare
      end if

!2025-02-05      
!$omp workshare
      fld(11,:,:,:)= gradbr(:,:,:)
      fld(12,:,:,:)= gradbz(:,:,:)
      fld(13,:,:,:)= gradbphi(:,:,:)
      fld(14,:,:,:)= curvbr(:,:,:)
      fld(15,:,:,:)= curvbz(:,:,:)
      fld(16,:,:,:)= curvbphi(:,:,:)
!$omp end workshare

      if(ispecies.eq.0)then
!$omp workshare
        fld(23,:,:,:)= dns_e0(:,:,:)
        fld(24,:,:,:)= dns_e0_r(:,:,:)
        fld(25,:,:,:)= dns_e0_z(:,:,:)
        fld(26,:,:,:)= dns_e0_phi(:,:,:)
        fld(27,:,:,:)= temp_e0(:,:,:)
        fld(28,:,:,:)= temp_e0_r(:,:,:)
        fld(29,:,:,:)= temp_e0_z(:,:,:)
        fld(30,:,:,:)= temp_e0_phi(:,:,:)
!$omp end workshare
      else if(ispecies.eq.1)then
!$omp workshare
        fld(23,:,:,:)= dns_i0(:,:,:)
        fld(24,:,:,:)= dns_i0_r(:,:,:)
        fld(25,:,:,:)= dns_i0_z(:,:,:)
        fld(26,:,:,:)= dns_i0_phi(:,:,:)
        fld(27,:,:,:)= temp_i0(:,:,:)
        fld(28,:,:,:)= temp_i0_r(:,:,:)
        fld(29,:,:,:)= temp_i0_z(:,:,:)
        fld(30,:,:,:)= temp_i0_phi(:,:,:)
!$omp end workshare
      else if(ispecies.eq.2)then
!$omp workshare
        fld(23,:,:,:)= dns_a0(:,:,:)
        fld(24,:,:,:)= dns_a0_r(:,:,:)
        fld(25,:,:,:)= dns_a0_z(:,:,:)
        fld(26,:,:,:)= dns_a0_phi(:,:,:)
        fld(27,:,:,:)= temp_a0(:,:,:)
        fld(28,:,:,:)= temp_a0_r(:,:,:)
        fld(29,:,:,:)= temp_a0_z(:,:,:)
        fld(30,:,:,:)= temp_a0_phi(:,:,:)
!$omp end workshare
      end if
!$omp end parallel

!      call ftrace_region_end('push0')

!$omp parallel private(nsta,nend,n,ivect,kfl,ijk_a,ia,ja,ka &
!$omp& ,aaa,ar1,ar,az1,az,aphi1,aphi &
!$omp& ,bre,bze,bphie,b2,babse,babs0e,b1,b21,br1a,bz1a,bphi1a &
!$omp& ,b10,br10,bz10,bphi10,dbpar,dbpr,dbpz,dbpphi,rhopar,denom1,orbpr,orbpz,orbpphi &
!$omp& ,wvnmlr,wvnmlz,wvnmlphi,dgrdbr,dgrdbz,dgrdbphi,dppar,dppar2,dppar1,pvpar &
!$omp& ,w1nmlr,w1nmlz,w1nmlphi,detotal,energy &
!$omp& ,rminor,bmax,energy0,sigma,vpara0,vpara,dvpara &
!$omp& ,pphi_n,dpphi_n,prof,dprofdpsi &
!$omp& ,vlcpart,dvlcpart,vnrm,sd_factor,pt_factor,nudt,vrat,coef)
!$omp do schedule(dynamic,1000)
      do nn = 1, marker_num, nblkp
        nsta = nn
        nend = min(nn+nblkp-1,marker_num)

!      call ftrace_region_begin('push1')

      do n = nsta, nend
        ivect = n - nsta + 1

        ijk_a(ivect,1)=max(1,min(lr  -1,int((gc(1,n)-ma_mi_r)*dr1  ) + kr  ))
        ijk_a(ivect,2)=max(1,min(lz  -1,int(gc(2,n)          *dz1  ) + kz  ))
        ijk_a(ivect,3)=max(1,min(lphi-1,int(gc(3,n)          *dphi1) + kphi))

        ar1  = max(0.0d0, min(1.0d0, (gc(1,n) - ma_mi_r)*dr1 -  dble(ijk_a(ivect,1) - kr)  ) )
        ar   = 1.0d0 - ar1
        az1  = max(0.0d0, min(1.0d0, gc(2,n)*dz1 - dble(ijk_a(ivect,2) - kz) ) )
        az   = 1.0d0 - az1
        aphi1 = max(0.0d0, min(1.0d0, gc(3,n)*dphi1 - dble(ijk_a(ivect,3) - kphi) ) )
        aphi = 1.0d0 - aphi1

        aaa(ivect,1) = ar *az *aphi
        aaa(ivect,2) = ar1*az *aphi
        aaa(ivect,3) = ar *az1*aphi
        aaa(ivect,4) = ar1*az1*aphi
        aaa(ivect,5) = ar *az *aphi1
        aaa(ivect,6) = ar1*az *aphi1
        aaa(ivect,7) = ar *az1*aphi1
        aaa(ivect,8) = ar1*az1*aphi1
      end do


! fields at each particle position
! igyro=0: w/o FLR, igyro/=0: flp(1:7,:) is given in subr. extract_em

!      call ftrace_region_end('push1')
!      call ftrace_region_begin('push2')

!NEC$ outerloop_unroll(4)
      do kfl = kfl_start, nflp

      do n = nsta, nend
        ivect = n - nsta + 1

        ia=ijk_a(ivect,1)
        ja=ijk_a(ivect,2)
        ka=ijk_a(ivect,3)

        flp(kfl,n) = fld(kfl, ia, ja,  ka  )*aaa(ivect,1) + fld(kfl, ia+1,ja,  ka  )*aaa(ivect,2) &
                   + fld(kfl, ia, ja+1,ka  )*aaa(ivect,3) + fld(kfl, ia+1,ja+1,ka  )*aaa(ivect,4) &
                   + fld(kfl, ia, ja,  ka+1)*aaa(ivect,5) + fld(kfl, ia+1,ja,  ka+1)*aaa(ivect,6) &
                   + fld(kfl, ia, ja+1,ka+1)*aaa(ivect,7) + fld(kfl, ia+1,ja+1,ka+1)*aaa(ivect,8)
        end do

      end do

!      call ftrace_region_end('push2')
!      call ftrace_region_begin('push3')

      do n = nsta, nend

! flp(5:7,n): delta_br(z,phi)
        bre = flp(5,n) + flp(8,n)
        bze = flp(6,n) + flp(9,n)
        bphie = flp(7,n) + flp(10,n)

        b2  = bre**2 + bze**2 + bphie**2
        babse= max(eps_b, sqrt(b2) )
        babs0e= max(eps_b, sqrt(flp(8,n)**2 + flp(9,n)**2 + flp(10,n)**2) )
        b1 = 1.0d0/babse
        b21= 1.0d0/b2
        br1a = bre*b1
        bz1a = bze*b1
        bphi1a = bphie*b1
        bphi1e(n) = bphi1a !2016-02-4

        b10 = 1.0d0/babs0e
        br10 = flp(8,n)*b10
        bz10 = flp(9,n)*b10
        bphi10 = flp(10,n)*b10

        dbpar = br1a*br10 + bz1a*bz10 + bphi1a*bphi10
        dbpr  = br1a  - br10*dbpar
        dbpz  = bz1a  - bz10*dbpar
        dbpphi= bphi1a- bphi10*dbpar

! guiding center motion

        rhopar = gc(4,n)*ep1*b1

        denom1 = 1.0d0/(1.0d0 + rhopar*(br1a*flp(14,n) + bz1a*flp(15,n) + bphi1a*flp(16,n)))

        orbpr = (br1a + rhopar*flp(14,n))*denom1
        orbpz = (bz1a + rhopar*flp(15,n))*denom1
        orbpphi = (bphi1a + rhopar*flp(16,n))*denom1

! e x b drift

        wvnmlr = (flp(3,n)*bze-flp(2,n)*bphie)*b21*denom1
        wvnmlz = (flp(1,n)*bphie-flp(3,n)*bre)*b21*denom1
        wvnmlphi = (flp(2,n)*bre-flp(1,n)*bze)*b21*denom1

! grad-b drift

        dgrdbr = gc(5,n)*(bphie*flp(12,n) - bze*flp(13,n))*ep1*b21*denom1
        dgrdbz = gc(5,n)*(bre*flp(13,n) - bphie*flp(11,n))*ep1*b21*denom1
        dgrdbphi = gc(5,n)*(bze*flp(11,n) - bre*flp(12,n))*ep1*b21*denom1

! mirror force

        dppar =-gc(5,n)*( flp(11,n)*orbpr &
                      + flp(12,n)*orbpz &
                      + flp(13,n)*orbpphi &
                      )*dt

!2012-08-31
        dppar2=-gc(5,n)*( flp(11,n)*dbpr*denom1 &
                      + flp(12,n)*dbpz*denom1 &
                      + flp(13,n)*orbpphi &
                      )*dt
!2012-08-31 end

! aceeleration due to electric field and curvature drift

!2012-07-07
        dppar1 = ep*((flp(1,n)*flp(14,n) + flp(2,n)*flp(15,n) + flp(3,n)*flp(16,n) &
                      )*rhopar &
                    + flp(4,n) &
                     )*dt*denom1
!2012-07-07 end

! total drift velocity

        pvpar = gc(4,n)*amsp1
        dgc(1,n) = dt*(pvpar*orbpr + wvnmlr + dgrdbr)*gc(7,n)
        dgc(2,n) = dt*(pvpar*orbpz + wvnmlz + dgrdbz)*gc(7,n)
        dgc(3,n) = dt*(pvpar*orbpphi + wvnmlphi + dgrdbphi)/gc(1,n)*gc(7,n)
        dgc(4,n) =(dppar + dppar1)*gc(7,n)

! temporal evolution of weight of high-energy ion particle

        w1nmlr = wvnmlr + pvpar*dbpr*denom1
        w1nmlz = wvnmlz + pvpar*dbpz*denom1
        w1nmlphi = wvnmlphi + pvpar*dbpphi*denom1

        dgc(8,n) =(ep*(w1nmlr*flp(18,n) + w1nmlz*flp(19,n) + w1nmlphi*flp(20,n))*dt &
                  + bphi1a*(w1nmlr*gc(4,n)*dt + gc(1,n)*(dppar1 + dppar2) ) & !2012-08-31
                  )*gc(7,n)

        detotal =(  ep*(flp(1,n)*dgrdbr + flp(2,n)*dgrdbz + flp(3,n)*dgrdbphi)*dt &
                + dppar1*pvpar )*gc(7,n) &
! the following term considers 'mu* v * grad(B0 - B)', 2025-04-04
                + gc(5,n)*(dgc(1,n)*(flp(31,n)-flp(11,n) ) &
                          +dgc(2,n)*(flp(32,n)-flp(12,n) ) &
                          +dgc(3,n)*(flp(33,n)-flp(13,n) )*gc(1,n) &
                          )
!2025-04-27s
        energy0 = max(1.0d-30, 0.50d0*amassp*pvpar**2 + gc(5,n)*babs0e)
        clambda(n) = gc(5,n)*b0/energy0
        v(n) = sqrt(2.0d0*energy0*amsp1)

!        energy = max(1.0d-30, 0.50d0*amassp*pvpar**2 + gc(5,n)*babse)
!        clambda(n) = gc(5,n)*b0/energy
!        v(n) = sqrt(2.0d0*energy*amsp1)
!2025-04-27e

! weight evolution : weight = f - f0
! d weight /dt = - d f0 /dt
! f0 = f_nrml*prof(psinrm)/(v**3 + flp(21,n)**3)*0.50*erfc((v(n)-valpha)/deltav)

!        rminor = sqrt( (gc(1,n)-raxis)**2 &
!                     + ((gc(2,n)-0.50d0*zleng) )**2 &
!                     )
!        bmax = b0*raxis/(raxis-rminor)
!        energy0 = 0.50d0*amsp1*gc(4,n)**2 + gc(5,n)*babs0e
!        sigma = 0.50d0*(1.0d0 + sign(1.0d0, energy0-gc(5,n)*bmax) ) &
!                      *sign(1.0d0,gc(4,n) )
!        vpara0 = sqrt(2.0d0*(energy0-gc(5,n)*bmin)*amsp1)
!        vpara = vpara0*sigma
!        dvpara = sigma/vpara0*gc(4,n)*dppar1*amsp1**2

!        pphi_n = gc(8,n) - amassp*raxis*vpara
!        dpphi_n= dgc(8,n) - amassp*raxis*dvpara

!        prof = exp(pphi_n/(ep*psimax*0.37d0) )
!        dprofdpsi = prof/(ep*psimax*0.37d0)

!        dwpsi(n) = dpphi_n*dprofdpsi*gc(10,n)
!        dwenrc(n)= detotal/(amassp*v(n) )*prof*gc(10,n)

!2016-08-05s
        dwpsi(n) = ( w1nmlr*flp(24,n) + w1nmlz*flp(25,n) + w1nmlphi*flp(26,n) &
                    +(w1nmlr*flp(28,n) + w1nmlz*flp(29,n) + w1nmlphi*flp(30,n) ) &
                    *0.50d0*(amassp*v(n)**2/flp(27,n) - 3.0d0)*flp(23,n)/flp(27,n) &
                   )*dt*gc(7,n)*gc(10,n)

!        dwpsi(n) = (w1nmlr*flp(24,n) + w1nmlz*flp(25,n) + w1nmlphi*flp(26,n))*dt & !2016-01-09
!                   *gc(7,n)*gc(10,n)
!2016-08-05e

        dwenrc(n)= detotal/(amassp*v(n) )*flp(23,n)*gc(10,n) !2016-01-09

      end do

!      call ftrace_region_end('push3')
!      call ftrace_region_begin('push4')

      if(type.eq.0)then

        do n = nsta, nend
!2016-08-05s
          vlcpart = exp(-0.50d0*amassp*v(n)**2/flp(27,n) )*flp(27,n)**(-1.5d0) !2013-07-17
          dvlcpart = -amassp*v(n)/flp(27,n)*vlcpart 
!          vlcpart = exp(-0.50d0*amassp*v(n)**2/temp)
!          dvlcpart = -amassp*v(n)/temp*vlcpart 
!2016-08-05e
          dwpsi(n) = dwpsi(n)*vlcpart
          dwenrc(n) = dwenrc(n)*dvlcpart
          dgc(6,n) = - dwpsi(n) - dwenrc(n) &
                     - gc(6,n)*nu_krook_dt & !Krook operator at normalized psi > psi_edge, flp(17)=psi
                     *(0.50d0 + sign(0.50d0, 1.0d0 - flp(17,n)/psimax - psi_edge) ) !2024-04-23
        end do

      else if(type.eq.1.or.type.eq.2)then

        do n = nsta, nend
          vnrm =(v(n)-valpha)/deltav
          sd_factor = 1.0d0/(v(n)**3 + flp(21,n)**3)
          vlcpart = 0.50d0*erfc(vnrm)*sd_factor
          dvlcpart= -( 3.0d0*v(n)**2*vlcpart &
                     + exp(-vnrm**2)*sqr_pi1/deltav &
                     )*sd_factor

          dwpsi(n) = dwpsi(n)*vlcpart
          dwenrc(n) = dwenrc(n)*dvlcpart
          dgc(6,n) = - dwpsi(n) - dwenrc(n) &
                     - gc(6,n)*nu_krook_dt & !Krook operator at normalized psi > psi_edge, flp(17)=psi
                     *(0.50d0 + sign(0.50d0, 1.0d0 - flp(17,n)/psimax - psi_edge) ) !2024-04-23
        end do

!2012-06-17
      else if(type.eq.3)then

        do n = nsta, nend
          vnrm =(v(n)-valpha)/deltav
          sd_factor = 1.0d0/(v(n)**3 + flp(21,n)**3)
          pt_factor = exp(-(clambda(n)-clambda0)**2/dclambda**2)

          vlcpart = 0.50d0*erfc(vnrm)*sd_factor*pt_factor
          dvlcpart= -( 3.0d0*v(n)**2*vlcpart &
                     + exp(-vnrm**2)*sqr_pi1/deltav*pt_factor &
                     )*sd_factor &
                  + 4.0d0*vlcpart*clambda(n)*(clambda(n)-clambda0) &
                         /(v(n)*dclambda**2)

          dwpsi(n) = dwpsi(n)*vlcpart
          dwenrc(n) = dwenrc(n)*dvlcpart
          dgc(6,n) = - dwpsi(n) - dwenrc(n) &
                     - gc(6,n)*nu_krook_dt & !Krook operator at normalized psi > psi_edge, flp(17)=psi
                     *(0.50d0 + sign(0.50d0, 1.0d0 - flp(17,n)/psimax - psi_edge) ) !2024-04-23
        end do
!2012-06-17 end

      end if

! slowing down
      if(type.eq.-5)then
        do n = nsta, nend
            dgc(6,n) = 0.0d0 !full-f
            nudt = flp(22,n)*dt*gc(7,n)
            vrat = flp(21,n)/v(n)
            coef =(1.0d0 + vrat**3)*nudt
            dgc(10,n) =-3.0d0*nudt*gc(10,n)
            dgc(5,n)= -2.0d0*gc(5,n)*coef
            dgc(4,n) = dgc(4,n)- gc(4,n)*coef
            dgc(8,n) = dgc(8,n) - gc(4,n)*coef*gc(1,n)*bphi1e(n)
        end do
      end if

!      call ftrace_region_end('push4')

      end do !nn
!$omp end parallel

end subroutine push
!--------------------------------------------------------------------
subroutine density(marker_num,mass,gc &
                  ,dns,mom,ppara,pperp,dns0,mom0,ppara0,pperp0)
! calculate pressure
! modified for NEC SX-Aurora TSUBASA 2020-07-10
! gc_x(:,n) = 1:r, 2:z, 3:phi, 4:p, 5:mu, 6:weight, 7:active, 8:p_phi, 9:fff, 10:fnrml
!--------------------------------------------------------------------
      use mpiset
      use grid
      use field, only:babs,babs0
      implicit none

      integer::marker_num,nvec
      real(8)::gc(ngc2,marker_num)
      real(8)::mass
      real(8)::dns(lr,lz,lphi),mom(lr,lz,lphi)
      real(8)::ppara(lr,lz,lphi),pperp(lr,lz,lphi)
      real(8)::dns0(lr,lz,lphi),mom0(lr,lz,lphi)
      real(8)::ppara0(lr,lz,lphi),pperp0(lr,lz,lphi)
!NEC      real(8)::wdns(lr,lz,lphi,lpara),wmom(lr,lz,lphi,lpara)
!NEC      real(8)::wpar(lr,lz,lphi,lpara),wprp(lr,lz,lphi,lpara)
!NEC      integer::ijk_a(3,marker_each)
!NEC      real(8)::aaa(8,marker_each)
      integer::i,j,k,l,m,n,ia,ia1,ja,ja1,ka,ka1,lr1,lz1
      real(8)::ar,ar1,az,az1,aphi,aphi1
      real(8)::aaa1,aaa2,aaa3,aaa4,aaa5,aaa6,aaa7,aaa8
      real(8)::vol,cwwa,cwwc
!      real(8)::d0,d1,d2,d3,d4,d5,t1,t2,t3,t4,t5
      real(8)::dr1,dz1,dphi1,ma_mi_r
      integer::kr,kz,kphi
      real(8)::p0,p1,p2,mu1
      integer::n_min,n_max !2013-05-22

!NEC start
!      integer,parameter::nblkd=255 !20NOV
      integer,parameter::nblkd=15
      real(8)::wkdns(nblkd,lr,lz,lphi,lpara),wkmom(nblkd,lr,lz,lphi,lpara)
      real(8)::wkpar(nblkd,lr,lz,lphi,lpara),wkprp(nblkd,lr,lz,lphi,lpara)
      integer::nvec0,nmod

!$omp parallel do
        do l = 1, lpara
          do i = 1, lrzphi
            do k = 1, nblkd
              wkdns(k,i,1,1,l)=0.0d0
              wkmom(k,i,1,1,l)=0.0d0
              wkpar(k,i,1,1,l)=0.0d0
              wkprp(k,i,1,1,l)=0.0d0
            end do
          end do
        end do
!NEC end

      dr1 = 1.0d0/dr
      dz1 = 1.0d0/dz
      dphi1 = 1.0d0/dphi

      kr   = 1 - kr_offset(my_rank)
      kz   = 1 - kz_offset(my_rank)
      kphi = 1 + lphi_shd - kphi_offset(my_rank)
!      kphi = 3 - kphi_offset(my_rank)
      ma_mi_r = major_r-minor_r

!$omp parallel do
        do i = 1, lrzphi
          dns(i,1,1) = 0.0d0
          mom(i,1,1) = 0.0d0
          ppara(i,1,1) = 0.0d0
          pperp(i,1,1) = 0.0d0
        end do

      nvec = marker_num/lpara
      nvec0 = marker_num/lpara
      nmod  = mod(marker_num,lpara)

!$omp parallel private(ar1,ar,az1,az,aphi1,aphi &
!$omp& ,nvec,n,m,k,ia,ja,ka,p0,p1,p2,mu1 &
!$omp& ,aaa1,aaa2,aaa3,aaa4,aaa5,aaa6,aaa7,aaa8)
!$omp do
      do l = 1, lpara
        nvec = nvec0 + min(1,nmod/l)

      do m = 1, nvec, nblkd

       do k = 1, min(nblkd, nvec-m+1)
        n = m + k - 1 + nvec0*(l-1) + min(l-1,nmod)

! gc_x(:,n) = 1:r, 2:z, 3:phi, 4:p, 5:mu, 6:weight, 7:active, 8:p_phi, 9:fff, 10:fnrml
        ia=max(1,min(lr  -1,int((gc(1,n)-ma_mi_r)*dr1  ) + kr  ))
        ja=max(1,min(lz  -1,int(gc(2,n)          *dz1  ) + kz  ))
        ka=max(1,min(lphi-1,int(gc(3,n)        *dphi1) + kphi))

        ar1 = (gc(1,n)-grr(ia,ja,ka ) ) *dr1
        ar  = 1.0d0 - ar1
        az1 = (gc(2,n)-gzz(ia,ja,ka ) ) *dz1
        az  = 1.0d0 - az1
        aphi1 = (gc(3,n)-gphi(ia,ja,ka ) ) *dphi1
        aphi  = 1.0d0 - aphi1

        aaa1 = ar *az *aphi
        aaa2 = ar1*az *aphi
        aaa3 = ar *az1*aphi
        aaa4 = ar1*az1*aphi
        aaa5 = ar *az *aphi1
        aaa6 = ar1*az *aphi1
        aaa7 = ar *az1*aphi1
        aaa8 = ar1*az1*aphi1

        p0 =  gc(6,n)*gc(7,n)
        p1 =  gc(4,n)*p0
        p2 =  gc(4,n)*p1
        mu1 = gc(5,n)*p0

        wkdns(k,ia,  ja,  ka,  l) = wkdns(k,ia,  ja,  ka,  l) + aaa1*p0
        wkdns(k,ia+1,ja,  ka,  l) = wkdns(k,ia+1,ja,  ka,  l) + aaa2*p0
        wkdns(k,ia,  ja+1,ka,  l) = wkdns(k,ia,  ja+1,ka,  l) + aaa3*p0
        wkdns(k,ia+1,ja+1,ka,  l) = wkdns(k,ia+1,ja+1,ka,  l) + aaa4*p0
        wkdns(k,ia,  ja,  ka+1,l) = wkdns(k,ia,  ja,  ka+1,l) + aaa5*p0
        wkdns(k,ia+1,ja,  ka+1,l) = wkdns(k,ia+1,ja,  ka+1,l) + aaa6*p0
        wkdns(k,ia,  ja+1,ka+1,l) = wkdns(k,ia,  ja+1,ka+1,l) + aaa7*p0
        wkdns(k,ia+1,ja+1,ka+1,l) = wkdns(k,ia+1,ja+1,ka+1,l) + aaa8*p0

        wkmom(k,ia,  ja,  ka,  l) = wkmom(k,ia,  ja,  ka,  l) + aaa1*p1
        wkmom(k,ia+1,ja,  ka,  l) = wkmom(k,ia+1,ja,  ka,  l) + aaa2*p1
        wkmom(k,ia,  ja+1,ka,  l) = wkmom(k,ia,  ja+1,ka,  l) + aaa3*p1
        wkmom(k,ia+1,ja+1,ka,  l) = wkmom(k,ia+1,ja+1,ka,  l) + aaa4*p1
        wkmom(k,ia,  ja,  ka+1,l) = wkmom(k,ia,  ja,  ka+1,l) + aaa5*p1
        wkmom(k,ia+1,ja,  ka+1,l) = wkmom(k,ia+1,ja,  ka+1,l) + aaa6*p1
        wkmom(k,ia,  ja+1,ka+1,l) = wkmom(k,ia,  ja+1,ka+1,l) + aaa7*p1
        wkmom(k,ia+1,ja+1,ka+1,l) = wkmom(k,ia+1,ja+1,ka+1,l) + aaa8*p1

        wkpar(k,ia,  ja,  ka,  l) = wkpar(k,ia,  ja,  ka,  l) + aaa1*p2
        wkpar(k,ia+1,ja,  ka,  l) = wkpar(k,ia+1,ja,  ka,  l) + aaa2*p2
        wkpar(k,ia,  ja+1,ka,  l) = wkpar(k,ia,  ja+1,ka,  l) + aaa3*p2
        wkpar(k,ia+1,ja+1,ka,  l) = wkpar(k,ia+1,ja+1,ka,  l) + aaa4*p2
        wkpar(k,ia,  ja,  ka+1,l) = wkpar(k,ia,  ja,  ka+1,l) + aaa5*p2
        wkpar(k,ia+1,ja,  ka+1,l) = wkpar(k,ia+1,ja,  ka+1,l) + aaa6*p2
        wkpar(k,ia,  ja+1,ka+1,l) = wkpar(k,ia,  ja+1,ka+1,l) + aaa7*p2
        wkpar(k,ia+1,ja+1,ka+1,l) = wkpar(k,ia+1,ja+1,ka+1,l) + aaa8*p2

        wkprp(k,ia,  ja,  ka,  l) = wkprp(k,ia,  ja,  ka,  l) + aaa1*mu1
        wkprp(k,ia+1,ja,  ka,  l) = wkprp(k,ia+1,ja,  ka,  l) + aaa2*mu1
        wkprp(k,ia,  ja+1,ka,  l) = wkprp(k,ia,  ja+1,ka,  l) + aaa3*mu1
        wkprp(k,ia+1,ja+1,ka,  l) = wkprp(k,ia+1,ja+1,ka,  l) + aaa4*mu1
        wkprp(k,ia,  ja,  ka+1,l) = wkprp(k,ia,  ja,  ka+1,l) + aaa5*mu1
        wkprp(k,ia+1,ja,  ka+1,l) = wkprp(k,ia+1,ja,  ka+1,l) + aaa6*mu1
        wkprp(k,ia,  ja+1,ka+1,l) = wkprp(k,ia,  ja+1,ka+1,l) + aaa7*mu1
        wkprp(k,ia+1,ja+1,ka+1,l) = wkprp(k,ia+1,ja+1,ka+1,l) + aaa8*mu1

       end do
      end do
      end do
!$omp end parallel

        do l = 1, lpara
!$omp parallel do
          do i = 1, lrzphi
            do k = 1, nblkd
              dns(i,1,1) = dns(i,1,1) + wkdns(k,i,1,1,l)
              mom(i,1,1) = mom(i,1,1) + wkmom(k,i,1,1,l)
              ppara(i,1,1) = ppara(i,1,1) + wkpar(k,i,1,1,l)
              pperp(i,1,1) = pperp(i,1,1) + wkprp(k,i,1,1,l)
            end do
          end do
        end do

! smoothing
      cwwa = 0.5d0
!      cwwc =-1.d0/6.d0 

!2024-12-21s, correction suggested by Panith Adulsiriswad
      call periodic_particle_mlt4b(dns,mom,ppara,pperp)
      call partsm1(dns,cwwa)
      call partsm1(mom,cwwa)

      call partsm1(ppara,cwwa)
      call partsm1(pperp,cwwa)
!2024-12-21e

! calculate density (per volume)

!      do k = 1, lphi
!      do j = 1, lz
!      do i = 1, lr
!$omp parallel do private(vol)
      do i = 1, lrzphi
        vol = 1.0d0/(grr(i,1,1)*dr*dz*dphi)
        dns(i,1,1) = dns(i,1,1)*vol
        mom(i,1,1) = mom(i,1,1)*vol
        ppara(i,1,1) = ppara(i,1,1)*vol/mass
        pperp(i,1,1) = pperp(i,1,1)*vol*babs(i,1,1)
      end do
!      end do
!      end do


! wall effect

      lr1 = lr - 1
      lz1 = lz - 1

      if(my_rank_r.eq.(mpi_proc_r-1))then
!      do k = 1, lphi
!      do j = 1, lz
!$omp parallel do
       do j = 1, lzphi
          dns(lr1,j,1) = dns(lr1,j,1) + dns(lr,j,1)
          mom(lr1,j,1) = mom(lr1,j,1) + mom(lr,j,1)
          dns(lr,j,1) = 0.0d0
          mom(lr,j,1) = 0.0d0

          ppara(lr1,j,1) = ppara(lr1,j,1) + ppara(lr,j,1)
          pperp(lr1,j,1) = pperp(lr1,j,1) + pperp(lr,j,1)
          ppara(lr,j,1) = 0.0d0
          pperp(lr,j,1) = 0.0d0
       end do
!      end do
      end if

      if(my_rank_r.eq.0)then
!      do k = 1, lphi
!      do j = 1, lz
!$omp parallel do
       do j = 1, lzphi
          dns(2,j,1) = dns(2,j,1) + dns(1,j,1)
          mom(2,j,1) = mom(2,j,1) + mom(1,j,1)
          dns(1,j,1) = 0.0d0
          mom(1,j,1) = 0.0d0

          ppara(2,j,1) = ppara(2,j,1) + ppara(1,j,1)
          pperp(2,j,1) = pperp(2,j,1) + pperp(1,j,1)
          ppara(1,j,1) = 0.0d0
          pperp(1,j,1) = 0.0d0
       end do
!      end do
      end if

      if(my_rank_z.eq.(mpi_proc_z-1))then
!$omp parallel do
       do k = 1, lphi
       do i = 1, lr
          dns(i,lz1,k) = dns(i,lz1,k) + dns(i,lz,k)
          mom(i,lz1,k) = mom(i,lz1,k) + mom(i,lz,k)
          dns(i,lz,k) = 0.0d0
          mom(i,lz,k) = 0.0d0

          ppara(i,lz1,k) = ppara(i,lz1,k) + ppara(i,lz,k)
          pperp(i,lz1,k) = pperp(i,lz1,k) + pperp(i,lz,k)
          ppara(i,lz,k) = 0.0d0
          pperp(i,lz,k) = 0.0d0
       end do
       end do
      end if 

      if(my_rank_z.eq.0)then
!$omp parallel do
       do k = 1, lphi
       do i = 1, lr
          dns(i,2,k) = dns(i,2,k) + dns(i,1,k)
          mom(i,2,k) = mom(i,2,k) + mom(i,1,k)
          dns(i,1,k) = 0.0d0
          mom(i,1,k) = 0.0d0

          ppara(i,2,k) = ppara(i,2,k) + ppara(i,1,k)
          pperp(i,2,k) = pperp(i,2,k) + pperp(i,1,k)
          ppara(i,1,k) = 0.0d0
          pperp(i,1,k) = 0.0d0
       end do
       end do
      end if 

!      if(.not.flag_stored)then
! take only n=1 modes
!      n = 1
!      call n1(n,ppara)
!      call n1(n,pperp)

! take modes n_min <= n <=n_max
!      n_min = 1
!      n_max = 2 
!      call lowpass2_mlt2(n_min,n_max,dns,mom)
!      call lowpass2_mlt2(n_min,n_max,ppara,pperp)
!      end if


!      do k = 1, lphi
!      do j = 1, lz
!      do i = 1, lr
!$omp parallel do
      do i = 1, lrzphi
        dns(i,1,1) = dns(i,1,1) + dns0(i,1,1)
        mom(i,1,1) = mom(i,1,1) + mom0(i,1,1)
        ppara(i,1,1) = ppara(i,1,1) + ppara0(i,1,1)
        pperp(i,1,1) = pperp(i,1,1) + pperp0(i,1,1)/babs0(i,1,1)*babs(i,1,1)
      end do
!      end do
!      end do

!            call wall_clock(d5)
!            t5 = d5 - d4
!            if(my_rank.eq.0)then
!            write(6,*)'density 5, t5=',t5
!            end if

end subroutine density
!--------------------------------------------------------------------
subroutine moments(marker_num,mass,gc &
                  ,dns,mom,ppara,pperp,qpara,qperp &
                  ,dns0,mom0,ppara0,pperp0,qpara0,qperp0)
! calculate pressure
! modified for NEC SX-Aurora TSUBASA 2020-07-10
! gc_x(:,n) = 1:r, 2:z, 3:phi, 4:p, 5:mu, 6:weight, 7:active, 8:p_phi, 9:fff, 10:fnrml
!--------------------------------------------------------------------
      use mpiset
      use grid
      use field, only:babs,babs0
      use gyro, only:marker_each_gyro
      implicit none

      integer::marker_num,nvec
      real(8)::gc(ngc2,marker_num)
      real(8)::mass
      real(8)::dns(lr,lz,lphi),mom(lr,lz,lphi)
      real(8)::ppara(lr,lz,lphi),pperp(lr,lz,lphi)
      real(8)::qpara(lr,lz,lphi),qperp(lr,lz,lphi)
      real(8)::dns0(lr,lz,lphi),mom0(lr,lz,lphi)
      real(8)::ppara0(lr,lz,lphi),pperp0(lr,lz,lphi)
      real(8)::qpara0(lr,lz,lphi),qperp0(lr,lz,lphi)

!NEC      real(8)::wdns(lr,lz,lphi,lpara),wmom(lr,lz,lphi,lpara)
!      real(8)::wpar(lr,lz,lphi,lpara),wprp(lr,lz,lphi,lpara)
!      real(8)::wqar(lr,lz,lphi,lpara),wqrp(lr,lz,lphi,lpara)
!      integer::ijk_a(3,marker_each_gyro)
!NEC      real(8)::aaa(8,marker_each_gyro)

      integer::i,j,k,l,m,n,ia,ia1,ja,ja1,ka,ka1,lr1,lz1
      real(8)::ar,ar1,az,az1,aphi,aphi1
      real(8)::aaa1,aaa2,aaa3,aaa4,aaa5,aaa6,aaa7,aaa8
      real(8)::vol,cwwa,cwwc
!      real(8)::d0,d1,d2,d3,d4,d5,t1,t2,t3,t4,t5
      real(8)::dr1,dz1,dphi1,ma_mi_r
      integer::kr,kz,kphi
      real(8)::p0,p1,p2,mu1,p3,p1mu,mass1 !2015-09-21
      integer::n_min,n_max !2013-05-22
!NEC start
!      integer,parameter::nblkd=255 !20NOV
      integer,parameter::nblkd=15
      real(8)::wkdns(nblkd,lr,lz,lphi,lpara),wkmom(nblkd,lr,lz,lphi,lpara)
      real(8)::wkpar(nblkd,lr,lz,lphi,lpara),wkprp(nblkd,lr,lz,lphi,lpara)
      real(8)::wkqar(nblkd,lr,lz,lphi,lpara),wkqrp(nblkd,lr,lz,lphi,lpara)
      integer::nvec0,nmod

!$omp parallel do
        do l = 1, lpara
          do i = 1, lrzphi
            do k = 1, nblkd
              wkdns(k,i,1,1,l)=0.0d0
              wkmom(k,i,1,1,l)=0.0d0
              wkpar(k,i,1,1,l)=0.0d0
              wkprp(k,i,1,1,l)=0.0d0
              wkqar(k,i,1,1,l)=0.0d0
              wkqrp(k,i,1,1,l)=0.0d0
            end do
          end do
        end do
!NEC end

      mass1 = 1.0d0/mass

      dr1 = 1.0d0/dr
      dz1 = 1.0d0/dz
      dphi1 = 1.0d0/dphi

      kr   = 1 - kr_offset(my_rank)
      kz   = 1 - kz_offset(my_rank)
      kphi = 1 + lphi_shd - kphi_offset(my_rank)
!      kphi = 3 - kphi_offset(my_rank)
      ma_mi_r = major_r-minor_r


!$omp parallel do
        do i = 1, lrzphi
          dns(i,1,1) = 0.0d0
          mom(i,1,1) = 0.0d0
          ppara(i,1,1) = 0.0d0
          pperp(i,1,1) = 0.0d0
          qpara(i,1,1) = 0.0d0
          qperp(i,1,1) = 0.0d0
        end do


      nvec = marker_num/lpara
      nvec0 = marker_num/lpara
      nmod  = mod(marker_num,lpara)

!$omp parallel private(ar1,ar,az1,az,aphi1,aphi &
!$omp& ,nvec,n,m,k,ia,ja,ka,p0,p1,p2,mu1,p3,p1mu &
!$omp& ,aaa1,aaa2,aaa3,aaa4,aaa5,aaa6,aaa7,aaa8)
!$omp do
      do l = 1, lpara
        nvec = nvec0 + min(1,nmod/l)

      do m = 1, nvec, nblkd

       do k = 1, min(nblkd, nvec-m+1)
        n = m + k - 1 + nvec0*(l-1) + min(l-1,nmod)

! gc_x(:,n) = 1:r, 2:z, 3:phi, 4:p, 5:mu, 6:weight, 7:active, 8:p_phi, 9:fff, 10:fnrml
        ia=max(1,min(lr  -1,int((gc(1,n)-ma_mi_r)*dr1  ) + kr  ))
        ja=max(1,min(lz  -1,int(gc(2,n)          *dz1  ) + kz  ))
        ka=max(1,min(lphi-1,int(gc(3,n)        *dphi1) + kphi))

        ar1 = (gc(1,n)-grr(ia,ja,ka ) ) *dr1
        ar  = 1.0d0 - ar1
        az1 = (gc(2,n)-gzz(ia,ja,ka ) ) *dz1
        az  = 1.0d0 - az1
        aphi1 = (gc(3,n)-gphi(ia,ja,ka ) ) *dphi1
        aphi  = 1.0d0 - aphi1

        aaa1 = ar *az *aphi
        aaa2 = ar1*az *aphi
        aaa3 = ar *az1*aphi
        aaa4 = ar1*az1*aphi
        aaa5 = ar *az *aphi1
        aaa6 = ar1*az *aphi1
        aaa7 = ar *az1*aphi1
        aaa8 = ar1*az1*aphi1

        p0 =  gc(6,n)*gc(7,n)
        p1 =  gc(4,n)*p0
        p2 =  gc(4,n)*p1
        mu1 = gc(5,n)*p0
        p3  = gc(4,n)*p2
        p1mu= gc(4,n)*mu1

        wkdns(k,ia,  ja,  ka,  l) = wkdns(k,ia,  ja,  ka,  l) + aaa1*p0
        wkdns(k,ia+1,ja,  ka,  l) = wkdns(k,ia+1,ja,  ka,  l) + aaa2*p0
        wkdns(k,ia,  ja+1,ka,  l) = wkdns(k,ia,  ja+1,ka,  l) + aaa3*p0
        wkdns(k,ia+1,ja+1,ka,  l) = wkdns(k,ia+1,ja+1,ka,  l) + aaa4*p0
        wkdns(k,ia,  ja,  ka+1,l) = wkdns(k,ia,  ja,  ka+1,l) + aaa5*p0
        wkdns(k,ia+1,ja,  ka+1,l) = wkdns(k,ia+1,ja,  ka+1,l) + aaa6*p0
        wkdns(k,ia,  ja+1,ka+1,l) = wkdns(k,ia,  ja+1,ka+1,l) + aaa7*p0
        wkdns(k,ia+1,ja+1,ka+1,l) = wkdns(k,ia+1,ja+1,ka+1,l) + aaa8*p0

        wkmom(k,ia,  ja,  ka,  l) = wkmom(k,ia,  ja,  ka,  l) + aaa1*p1
        wkmom(k,ia+1,ja,  ka,  l) = wkmom(k,ia+1,ja,  ka,  l) + aaa2*p1
        wkmom(k,ia,  ja+1,ka,  l) = wkmom(k,ia,  ja+1,ka,  l) + aaa3*p1
        wkmom(k,ia+1,ja+1,ka,  l) = wkmom(k,ia+1,ja+1,ka,  l) + aaa4*p1
        wkmom(k,ia,  ja,  ka+1,l) = wkmom(k,ia,  ja,  ka+1,l) + aaa5*p1
        wkmom(k,ia+1,ja,  ka+1,l) = wkmom(k,ia+1,ja,  ka+1,l) + aaa6*p1
        wkmom(k,ia,  ja+1,ka+1,l) = wkmom(k,ia,  ja+1,ka+1,l) + aaa7*p1
        wkmom(k,ia+1,ja+1,ka+1,l) = wkmom(k,ia+1,ja+1,ka+1,l) + aaa8*p1

        wkpar(k,ia,  ja,  ka,  l) = wkpar(k,ia,  ja,  ka,  l) + aaa1*p2
        wkpar(k,ia+1,ja,  ka,  l) = wkpar(k,ia+1,ja,  ka,  l) + aaa2*p2
        wkpar(k,ia,  ja+1,ka,  l) = wkpar(k,ia,  ja+1,ka,  l) + aaa3*p2
        wkpar(k,ia+1,ja+1,ka,  l) = wkpar(k,ia+1,ja+1,ka,  l) + aaa4*p2
        wkpar(k,ia,  ja,  ka+1,l) = wkpar(k,ia,  ja,  ka+1,l) + aaa5*p2
        wkpar(k,ia+1,ja,  ka+1,l) = wkpar(k,ia+1,ja,  ka+1,l) + aaa6*p2
        wkpar(k,ia,  ja+1,ka+1,l) = wkpar(k,ia,  ja+1,ka+1,l) + aaa7*p2
        wkpar(k,ia+1,ja+1,ka+1,l) = wkpar(k,ia+1,ja+1,ka+1,l) + aaa8*p2

        wkprp(k,ia,  ja,  ka,  l) = wkprp(k,ia,  ja,  ka,  l) + aaa1*mu1
        wkprp(k,ia+1,ja,  ka,  l) = wkprp(k,ia+1,ja,  ka,  l) + aaa2*mu1
        wkprp(k,ia,  ja+1,ka,  l) = wkprp(k,ia,  ja+1,ka,  l) + aaa3*mu1
        wkprp(k,ia+1,ja+1,ka,  l) = wkprp(k,ia+1,ja+1,ka,  l) + aaa4*mu1
        wkprp(k,ia,  ja,  ka+1,l) = wkprp(k,ia,  ja,  ka+1,l) + aaa5*mu1
        wkprp(k,ia+1,ja,  ka+1,l) = wkprp(k,ia+1,ja,  ka+1,l) + aaa6*mu1
        wkprp(k,ia,  ja+1,ka+1,l) = wkprp(k,ia,  ja+1,ka+1,l) + aaa7*mu1
        wkprp(k,ia+1,ja+1,ka+1,l) = wkprp(k,ia+1,ja+1,ka+1,l) + aaa8*mu1

        wkqar(k,ia,  ja,  ka,  l) = wkqar(k,ia,  ja,  ka,  l) + aaa1*p3
        wkqar(k,ia+1,ja,  ka,  l) = wkqar(k,ia+1,ja,  ka,  l) + aaa2*p3
        wkqar(k,ia,  ja+1,ka,  l) = wkqar(k,ia,  ja+1,ka,  l) + aaa3*p3
        wkqar(k,ia+1,ja+1,ka,  l) = wkqar(k,ia+1,ja+1,ka,  l) + aaa4*p3
        wkqar(k,ia,  ja,  ka+1,l) = wkqar(k,ia,  ja,  ka+1,l) + aaa5*p3
        wkqar(k,ia+1,ja,  ka+1,l) = wkqar(k,ia+1,ja,  ka+1,l) + aaa6*p3
        wkqar(k,ia,  ja+1,ka+1,l) = wkqar(k,ia,  ja+1,ka+1,l) + aaa7*p3
        wkqar(k,ia+1,ja+1,ka+1,l) = wkqar(k,ia+1,ja+1,ka+1,l) + aaa8*p3

        wkqrp(k,ia,  ja,  ka,  l) = wkqrp(k,ia,  ja,  ka,  l) + aaa1*p1mu
        wkqrp(k,ia+1,ja,  ka,  l) = wkqrp(k,ia+1,ja,  ka,  l) + aaa2*p1mu
        wkqrp(k,ia,  ja+1,ka,  l) = wkqrp(k,ia,  ja+1,ka,  l) + aaa3*p1mu
        wkqrp(k,ia+1,ja+1,ka,  l) = wkqrp(k,ia+1,ja+1,ka,  l) + aaa4*p1mu
        wkqrp(k,ia,  ja,  ka+1,l) = wkqrp(k,ia,  ja,  ka+1,l) + aaa5*p1mu
        wkqrp(k,ia+1,ja,  ka+1,l) = wkqrp(k,ia+1,ja,  ka+1,l) + aaa6*p1mu
        wkqrp(k,ia,  ja+1,ka+1,l) = wkqrp(k,ia,  ja+1,ka+1,l) + aaa7*p1mu
        wkqrp(k,ia+1,ja+1,ka+1,l) = wkqrp(k,ia+1,ja+1,ka+1,l) + aaa8*p1mu

       end do
      end do
      end do
!$omp end parallel

        do l = 1, lpara
!$omp parallel do
          do i = 1, lrzphi
            do k = 1, nblkd
              dns(i,1,1) = dns(i,1,1) + wkdns(k,i,1,1,l)
              mom(i,1,1) = mom(i,1,1) + wkmom(k,i,1,1,l)
              ppara(i,1,1) = ppara(i,1,1) + wkpar(k,i,1,1,l)
              pperp(i,1,1) = pperp(i,1,1) + wkprp(k,i,1,1,l)
              qpara(i,1,1) = qpara(i,1,1) + wkqar(k,i,1,1,l)
              qperp(i,1,1) = qperp(i,1,1) + wkqrp(k,i,1,1,l)
            end do
          end do
        end do


! smoothing

       cwwa = 0.5d0
!       cwwc =-1.d0/6.d0 

!2024-12-21s, correction suggested by Panith Adulsiriswad
       call periodic_particle_mlt6b(dns,mom,ppara,pperp,qpara,qperp)
       call partsm1(dns,cwwa)
       call partsm1(mom,cwwa)

       call partsm1(ppara,cwwa)
       call partsm1(pperp,cwwa)

       call partsm1(qpara,cwwa)
       call partsm1(qperp,cwwa)
!2024-12-21e

! calculate density (per volume)

!      do k = 1, lphi
!      do j = 1, lz
!      do i = 1, lr
!$omp parallel do private(vol)
      do i = 1, lrzphi
        vol = 1.0d0/(grr(i,1,1)*dr*dz*dphi)
        dns(i,1,1) = dns(i,1,1)*vol
        mom(i,1,1) = mom(i,1,1)*vol
        ppara(i,1,1) = ppara(i,1,1)*vol*mass1
        pperp(i,1,1) = pperp(i,1,1)*vol*babs(i,1,1)
        qpara(i,1,1) = qpara(i,1,1)*vol*mass1**2
        qperp(i,1,1) = qperp(i,1,1)*vol*babs(i,1,1)*mass1
      end do
!      end do
!      end do

! wall effect

      lr1 = lr - 1
      lz1 = lz - 1

      if(my_rank_r.eq.(mpi_proc_r-1))then
!      do k = 1, lphi
!      do j = 1, lz
!$omp parallel do
       do j = 1, lzphi
          dns(lr1,j,1) = dns(lr1,j,1) + dns(lr,j,1)
          mom(lr1,j,1) = mom(lr1,j,1) + mom(lr,j,1)
          dns(lr,j,1) = 0.0d0
          mom(lr,j,1) = 0.0d0

          ppara(lr1,j,1) = ppara(lr1,j,1) + ppara(lr,j,1)
          pperp(lr1,j,1) = pperp(lr1,j,1) + pperp(lr,j,1)
          ppara(lr,j,1) = 0.0d0
          pperp(lr,j,1) = 0.0d0

          qpara(lr1,j,1) = qpara(lr1,j,1) + qpara(lr,j,1)
          qperp(lr1,j,1) = qperp(lr1,j,1) + qperp(lr,j,1)
          qpara(lr,j,1) = 0.0d0
          qperp(lr,j,1) = 0.0d0
       end do
!      end do
      end if

      if(my_rank_r.eq.0)then
!      do k = 1, lphi
!      do j = 1, lz
!$omp parallel do
       do j = 1, lzphi
          dns(2,j,1) = dns(2,j,1) + dns(1,j,1)
          mom(2,j,1) = mom(2,j,1) + mom(1,j,1)
          dns(1,j,1) = 0.0d0
          mom(1,j,1) = 0.0d0

          ppara(2,j,1) = ppara(2,j,1) + ppara(1,j,1)
          pperp(2,j,1) = pperp(2,j,1) + pperp(1,j,1)
          ppara(1,j,1) = 0.0d0
          pperp(1,j,1) = 0.0d0

          qpara(2,j,1) = qpara(2,j,1) + qpara(1,j,1)
          qperp(2,j,1) = qperp(2,j,1) + qperp(1,j,1)
          qpara(1,j,1) = 0.0d0
          qperp(1,j,1) = 0.0d0
       end do
!      end do
      end if

      if(my_rank_z.eq.(mpi_proc_z-1))then
!$omp parallel do
       do k = 1, lphi
       do i = 1, lr
          dns(i,lz1,k) = dns(i,lz1,k) + dns(i,lz,k)
          mom(i,lz1,k) = mom(i,lz1,k) + mom(i,lz,k)
          dns(i,lz,k) = 0.0d0
          mom(i,lz,k) = 0.0d0

          ppara(i,lz1,k) = ppara(i,lz1,k) + ppara(i,lz,k)
          pperp(i,lz1,k) = pperp(i,lz1,k) + pperp(i,lz,k)
          ppara(i,lz,k) = 0.0d0
          pperp(i,lz,k) = 0.0d0

          qpara(i,lz1,k) = qpara(i,lz1,k) + qpara(i,lz,k)
          qperp(i,lz1,k) = qperp(i,lz1,k) + qperp(i,lz,k)
          qpara(i,lz,k) = 0.0d0
          qperp(i,lz,k) = 0.0d0
       end do
       end do
      end if 

      if(my_rank_z.eq.0)then
!$omp parallel do
       do k = 1, lphi
       do i = 1, lr
          dns(i,2,k) = dns(i,2,k) + dns(i,1,k)
          mom(i,2,k) = mom(i,2,k) + mom(i,1,k)
          dns(i,1,k) = 0.0d0
          mom(i,1,k) = 0.0d0

          ppara(i,2,k) = ppara(i,2,k) + ppara(i,1,k)
          pperp(i,2,k) = pperp(i,2,k) + pperp(i,1,k)
          ppara(i,1,k) = 0.0d0
          pperp(i,1,k) = 0.0d0

          qpara(i,2,k) = qpara(i,2,k) + qpara(i,1,k)
          qperp(i,2,k) = qperp(i,2,k) + qperp(i,1,k)
          qpara(i,1,k) = 0.0d0
          qperp(i,1,k) = 0.0d0
       end do
       end do
      end if 

      if(.not.flag_stored)then
! take modes n_min <= n <=n_max
!      n_min = 1
!      n_max = 2 
!      call lowpass2_mlt2(n_min,n_max,dns,mom)
!      call lowpass2_mlt2(n_min,n_max,ppara,pperp)
!      call lowpass2_mlt2(n_min,n_max,qpara,qperp)
      end if

!2014-08-28s
!$omp parallel do
      do i = 1, lrzphi
        dns(i,1,1) = dns(i,1,1) + dns0(i,1,1)
        mom(i,1,1) = mom(i,1,1) + mom0(i,1,1)
        ppara(i,1,1) = ppara(i,1,1) + ppara0(i,1,1)
        pperp(i,1,1) = pperp(i,1,1) + pperp0(i,1,1)/babs0(i,1,1)*babs(i,1,1)
        qpara(i,1,1) = qpara(i,1,1) + qpara0(i,1,1)
        qperp(i,1,1) = qperp(i,1,1) + qperp0(i,1,1)/babs0(i,1,1)*babs(i,1,1)
      end do
!2014-08-28e

!            call wall_clock(d5)
!            t5 = d5 - d4
!            if(my_rank.eq.0)then
!            write(6,*)'density 5, t5=',t5
!            end if

end subroutine moments
!--------------------------------------------------------------------
subroutine emf_gyro(marker_num,marker_num_gyro,gc,gyro_phys &
                   ,flp,flp_gyro)
! modified for GK simulation on NEC SX-Aurora TSUBASA 2022-01-11
! gc_x(:,n) = 1:r, 2:z, 3:phi, 4:p, 5:mu, 6:weight, 7:active, 8:p_phi, 9:fff, 10:fnrml
! flp(nflp, marker_num); nflp=30, nflp_gyro=7
!--------------------------------------------------------------------
      use mpiset
      use grid
      use field, only:babs,babs0,er,ez,ephi,epara,br,bz,bphi,br0,bz0,bphi0,fld
      use gyro, only:ngyro,marker_each_gyro
      implicit none

      integer::marker_num,nvec,marker_num_gyro
      real(8)::gc(ngc2,marker_num)
      real(8)::gyro_phys(2,marker_num_gyro)
      real(8)::flp(nflp,marker_num) !2015-07-08
      real(8)::flp_gyro(nflp_gyro,marker_num_gyro)
      integer::i,j,k,l,m,n,ia,ja,ka,lr1,lz1
      real(8)::ar,ar1,az,az1,aphi,aphi1
!      real(8)::aaa1,aaa2,aaa3,aaa4,aaa5,aaa6,aaa7,aaa8
      real(8)::dr1,dz1,dphi1,ma_mi_r
      integer::kr,kz,kphi
      integer::nc,in
      real(8)::rngr1,rngr1_ac

!NEC start
      integer::nn,nsta,nend
      integer,parameter::nblkp=1024
!2021-01-11s
      integer::ijk_a(marker_each_gyro,3)
      real(8)::aaa(marker_each_gyro,8)
      integer::nvec0,nmod
!2021-01-11e


      dr1 = 1.0d0/dr
      dz1 = 1.0d0/dz
      dphi1 = 1.0d0/dphi
      rngr1 = 1.0d0/dble(ngyro)

      kr   = 1 - kr_offset(my_rank)
      kz   = 1 - kz_offset(my_rank)
      kphi = 1 + lphi_shd - kphi_offset(my_rank)
!      kphi = 3 - kphi_offset(my_rank)
      ma_mi_r = major_r-minor_r


!$omp parallel
!$omp workshare

! for subroutine extract_em

!      flp_gyro = 0.0d0

      fld(1,:,:,:) = er(:,:,:)
      fld(2,:,:,:) = ez(:,:,:)
      fld(3,:,:,:) = ephi(:,:,:)
      fld(4,:,:,:) = epara(:,:,:)
      fld(5,:,:,:) = br(:,:,:) - br0(:,:,:)
      fld(6,:,:,:) = bz(:,:,:) - bz0(:,:,:)
      fld(7,:,:,:) = bphi(:,:,:) - bphi0(:,:,:)
!$omp end workshare
!$omp end parallel

!$omp parallel private(nsta,nend,n,in,nc &
!$omp&  ,ia,ja,ka &
!$omp&  ,ar1,ar,az1,az,aphi1,aphi &
!$omp&  ,rngr1_ac)
!$omp do schedule(dynamic,1000)
      do nn = 1, marker_num_gyro, nblkp
        nsta = nn
        nend = min(nn+nblkp-1,marker_num_gyro)

      do n = nsta, nend

        nc =(n-1)/ngyro + 1
        ijk_a(n,1)=max(1,min(lr  -1,int((gyro_phys(1,n)-ma_mi_r)*dr1  ) + kr  ))
        ijk_a(n,2)=max(1,min(lz  -1,int(gyro_phys(2,n)          *dz1  ) + kz  ))
        ijk_a(n,3)=max(1,min(lphi-1,int(gc(3,nc)          *dphi1) + kphi))

        ar1  = max(0.0d0, min(1.0d0, (gyro_phys(1,n) - ma_mi_r)*dr1 -  dble(ijk_a(n,1) - kr)  ) )
        ar   = 1.0d0 - ar1
        az1  = max(0.0d0, min(1.0d0, gyro_phys(2,n)*dz1 - dble(ijk_a(n,2) - kz) ) )
        az   = 1.0d0 - az1
        aphi1 = max(0.0d0, min(1.0d0, gc(3,nc)*dphi1 - dble(ijk_a(n,3) - kphi) ) )
        aphi = 1.0d0 - aphi1

        aaa(n,1) = ar *az *aphi
        aaa(n,2) = ar1*az *aphi
        aaa(n,3) = ar *az1*aphi
        aaa(n,4) = ar1*az1*aphi
        aaa(n,5) = ar *az *aphi1
        aaa(n,6) = ar1*az *aphi1
        aaa(n,7) = ar *az1*aphi1
        aaa(n,8) = ar1*az1*aphi1
      end do

!NEC$ outerloop_unroll(4)
      do in = 1, nflp_gyro

        do n = nsta, nend

          ia=ijk_a(n,1)
          ja=ijk_a(n,2)
          ka=ijk_a(n,3)

          flp_gyro(in,n) = fld(in, ia, ja,  ka  )*aaa(n,1) + fld(in, ia+1,ja,  ka  )*aaa(n,2) &
                         + fld(in, ia, ja+1,ka  )*aaa(n,3) + fld(in, ia+1,ja+1,ka  )*aaa(n,4) &
                         + fld(in, ia, ja,  ka+1)*aaa(n,5) + fld(in, ia+1,ja,  ka+1)*aaa(n,6) &
                         + fld(in, ia, ja+1,ka+1)*aaa(n,7) + fld(in, ia+1,ja+1,ka+1)*aaa(n,8)
        end do

      end do

!NEC$ outerloop_unroll(4)
       do in = 1, nflp_gyro

       do nc =(nsta-1)/ngyro+1, nend/ngyro
          rngr1_ac = rngr1*gc(7,nc)
          flp(in,nc) =(flp_gyro(in,ngyro*(nc-1)+1) &
                      +flp_gyro(in,ngyro*(nc-1)+2) &
                      +flp_gyro(in,ngyro*(nc-1)+3) &
                      +flp_gyro(in,ngyro*(nc-1)+4) &
                      )*rngr1_ac
       end do
       end do

      end do !nn
!$omp end parallel

end subroutine emf_gyro
!--------------------------------------------------------------------
subroutine density_gyro(marker_num,marker_num_gyro,mass,gc,gyro_phys &
                       ,dns,mom,ppara,pperp &
                       ,dns0,mom0,ppara0,pperp0)
! 2022-01-11, for FLR case
! modified for NEC SX-Aurora TSUBASA 2020-07-10
! gc_x(:,n) = 1:r, 2:z, 3:phi, 4:p, 5:mu, 6:weight, 7:active, 8:p_phi, 9:fff, 10:fnrml
!--------------------------------------------------------------------
      use mpiset
      use grid
      use field, only:babs,babs0
      use gyro, only:ngyro,marker_each_gyro
      implicit none

      integer::marker_num,nvec,marker_num_gyro
      real(8)::gc(ngc2,marker_num)
      real(8)::gyro_phys(2,marker_num_gyro)
      real(8)::mass
      real(8)::dns(lr,lz,lphi),mom(lr,lz,lphi)
      real(8)::ppara(lr,lz,lphi),pperp(lr,lz,lphi)
      real(8)::dns0(lr,lz,lphi),mom0(lr,lz,lphi)
      real(8)::ppara0(lr,lz,lphi),pperp0(lr,lz,lphi)
      integer::i,j,k,l,m,n,ia,ia1,ja,ja1,ka,ka1,lr1,lz1
      real(8)::ar,ar1,az,az1,aphi,aphi1
      real(8)::aaa1,aaa2,aaa3,aaa4,aaa5,aaa6,aaa7,aaa8
      real(8)::vol,cwwa,cwwc
!      real(8)::d0,d1,d2,d3,d4,d5,t1,t2,t3,t4,t5
      real(8)::dr1,dz1,dphi1,ma_mi_r
      integer::kr,kz,kphi
      real(8)::p0,p1,p2,mu1,mass1 !2015-09-21
      integer::n_min,n_max !2013-05-22
      real(8)::rngr1
      integer::nc
      integer::nvec0,nmod
!NEC start
!      integer,parameter::nblkd=255 !20NOV
!      integer,parameter::nblkd=127 !21JAN
      integer,parameter::nblkd=15
      real(8)::wkdns(nblkd,lr,lz,lphi,lpara),wkmom(nblkd,lr,lz,lphi,lpara)
      real(8)::wkpar(nblkd,lr,lz,lphi,lpara),wkprp(nblkd,lr,lz,lphi,lpara)

        do l = 1, lpara
!$omp parallel do
          do i = 1, lrzphi
            do k = 1, nblkd
              wkdns(k,i,1,1,l)=0.0d0
              wkmom(k,i,1,1,l)=0.0d0
              wkpar(k,i,1,1,l)=0.0d0
              wkprp(k,i,1,1,l)=0.0d0
            end do
          end do
        end do
!NEC end

      mass1 = 1.0d0/mass

      dr1 = 1.0d0/dr
      dz1 = 1.0d0/dz
      dphi1 = 1.0d0/dphi
      rngr1 = 1.0d0/dble(ngyro)

      kr   = 1 - kr_offset(my_rank)
      kz   = 1 - kz_offset(my_rank)
      kphi = 1 + lphi_shd - kphi_offset(my_rank)
!      kphi = 3 - kphi_offset(my_rank)
      ma_mi_r = major_r-minor_r

!$omp parallel do
        do i = 1, lrzphi
          dns(i,1,1) = 0.0d0
          mom(i,1,1) = 0.0d0
          ppara(i,1,1) = 0.0d0
          pperp(i,1,1) = 0.0d0
        end do


      nvec = marker_num_gyro/lpara
      nvec0 = marker_num_gyro/lpara
      nmod  = mod(marker_num_gyro,lpara)

!$omp parallel private(ar1,ar,az1,az,aphi1,aphi &
!$omp& ,nvec,m,k,n,nc,ia,ja,ka,p0,p1,p2,mu1 &
!$omp& ,aaa1,aaa2,aaa3,aaa4,aaa5,aaa6,aaa7,aaa8)
!$omp do
      do l = 1, lpara
        nvec = nvec0 + min(1,nmod/l)

      do m = 1, nvec, nblkd

       do k = 1, min(nblkd, nvec-m+1)
        n = m + k - 1 + nvec0*(l-1) + min(l-1,nmod)
        nc =(n-1)/ngyro + 1

! gc_x(:,n) = 1:r, 2:z, 3:phi, 4:p, 5:mu, 6:weight, 7:active, 8:p_phi, 9:fff, 10:fnrml
        ia=max(1,min(lr  -1,int((gyro_phys(1,n)-ma_mi_r)*dr1  ) + kr  ))
        ja=max(1,min(lz  -1,int(gyro_phys(2,n)          *dz1  ) + kz  ))
        ka=max(1,min(lphi-1,int(gc(3,nc)        *dphi1) + kphi))

        ar1 = (gyro_phys(1,n)-grr(ia,ja,ka ) ) *dr1
        ar  = 1.0d0 - ar1
        az1 = (gyro_phys(2,n)-gzz(ia,ja,ka ) ) *dz1
        az  = 1.0d0 - az1
        aphi1 = (gc(3,nc)-gphi(ia,ja,ka ) ) *dphi1
        aphi  = 1.0d0 - aphi1

        aaa1 = ar *az *aphi
        aaa2 = ar1*az *aphi
        aaa3 = ar *az1*aphi
        aaa4 = ar1*az1*aphi
        aaa5 = ar *az *aphi1
        aaa6 = ar1*az *aphi1
        aaa7 = ar *az1*aphi1
        aaa8 = ar1*az1*aphi1

        p0 =  gc(6,nc)*gc(7,nc)*rngr1 !rngr1 is multiplied
        p1 =  gc(4,nc)*p0
        p2 =  gc(4,nc)*p1
        mu1 = gc(5,nc)*p0

        wkdns(k,ia,  ja,  ka,  l) = wkdns(k,ia,  ja,  ka,  l) + aaa1*p0
        wkdns(k,ia+1,ja,  ka,  l) = wkdns(k,ia+1,ja,  ka,  l) + aaa2*p0
        wkdns(k,ia,  ja+1,ka,  l) = wkdns(k,ia,  ja+1,ka,  l) + aaa3*p0
        wkdns(k,ia+1,ja+1,ka,  l) = wkdns(k,ia+1,ja+1,ka,  l) + aaa4*p0
        wkdns(k,ia,  ja,  ka+1,l) = wkdns(k,ia,  ja,  ka+1,l) + aaa5*p0
        wkdns(k,ia+1,ja,  ka+1,l) = wkdns(k,ia+1,ja,  ka+1,l) + aaa6*p0
        wkdns(k,ia,  ja+1,ka+1,l) = wkdns(k,ia,  ja+1,ka+1,l) + aaa7*p0
        wkdns(k,ia+1,ja+1,ka+1,l) = wkdns(k,ia+1,ja+1,ka+1,l) + aaa8*p0

        wkmom(k,ia,  ja,  ka,  l) = wkmom(k,ia,  ja,  ka,  l) + aaa1*p1
        wkmom(k,ia+1,ja,  ka,  l) = wkmom(k,ia+1,ja,  ka,  l) + aaa2*p1
        wkmom(k,ia,  ja+1,ka,  l) = wkmom(k,ia,  ja+1,ka,  l) + aaa3*p1
        wkmom(k,ia+1,ja+1,ka,  l) = wkmom(k,ia+1,ja+1,ka,  l) + aaa4*p1
        wkmom(k,ia,  ja,  ka+1,l) = wkmom(k,ia,  ja,  ka+1,l) + aaa5*p1
        wkmom(k,ia+1,ja,  ka+1,l) = wkmom(k,ia+1,ja,  ka+1,l) + aaa6*p1
        wkmom(k,ia,  ja+1,ka+1,l) = wkmom(k,ia,  ja+1,ka+1,l) + aaa7*p1
        wkmom(k,ia+1,ja+1,ka+1,l) = wkmom(k,ia+1,ja+1,ka+1,l) + aaa8*p1

        wkpar(k,ia,  ja,  ka,  l) = wkpar(k,ia,  ja,  ka,  l) + aaa1*p2
        wkpar(k,ia+1,ja,  ka,  l) = wkpar(k,ia+1,ja,  ka,  l) + aaa2*p2
        wkpar(k,ia,  ja+1,ka,  l) = wkpar(k,ia,  ja+1,ka,  l) + aaa3*p2
        wkpar(k,ia+1,ja+1,ka,  l) = wkpar(k,ia+1,ja+1,ka,  l) + aaa4*p2
        wkpar(k,ia,  ja,  ka+1,l) = wkpar(k,ia,  ja,  ka+1,l) + aaa5*p2
        wkpar(k,ia+1,ja,  ka+1,l) = wkpar(k,ia+1,ja,  ka+1,l) + aaa6*p2
        wkpar(k,ia,  ja+1,ka+1,l) = wkpar(k,ia,  ja+1,ka+1,l) + aaa7*p2
        wkpar(k,ia+1,ja+1,ka+1,l) = wkpar(k,ia+1,ja+1,ka+1,l) + aaa8*p2

        wkprp(k,ia,  ja,  ka,  l) = wkprp(k,ia,  ja,  ka,  l) + aaa1*mu1
        wkprp(k,ia+1,ja,  ka,  l) = wkprp(k,ia+1,ja,  ka,  l) + aaa2*mu1
        wkprp(k,ia,  ja+1,ka,  l) = wkprp(k,ia,  ja+1,ka,  l) + aaa3*mu1
        wkprp(k,ia+1,ja+1,ka,  l) = wkprp(k,ia+1,ja+1,ka,  l) + aaa4*mu1
        wkprp(k,ia,  ja,  ka+1,l) = wkprp(k,ia,  ja,  ka+1,l) + aaa5*mu1
        wkprp(k,ia+1,ja,  ka+1,l) = wkprp(k,ia+1,ja,  ka+1,l) + aaa6*mu1
        wkprp(k,ia,  ja+1,ka+1,l) = wkprp(k,ia,  ja+1,ka+1,l) + aaa7*mu1
        wkprp(k,ia+1,ja+1,ka+1,l) = wkprp(k,ia+1,ja+1,ka+1,l) + aaa8*mu1
       end do
      end do
      end do
!$omp end parallel

        do l = 1, lpara
!$omp parallel do
          do i = 1, lrzphi
            do k = 1, nblkd
              dns(i,1,1) = dns(i,1,1) + wkdns(k,i,1,1,l)
              mom(i,1,1) = mom(i,1,1) + wkmom(k,i,1,1,l)
              ppara(i,1,1) = ppara(i,1,1) + wkpar(k,i,1,1,l)
              pperp(i,1,1) = pperp(i,1,1) + wkprp(k,i,1,1,l)
            end do
          end do
        end do

! smoothing

       cwwa = 0.5d0
!       cwwc =-1.d0/6.d0 

!2024-12-21s, correction suggested by Panith Adulsiriswad
       call periodic_particle_mlt4b(dns,mom,ppara,pperp)
       call partsm1(dns,cwwa)
       call partsm1(mom,cwwa)

       call partsm1(ppara,cwwa)
       call partsm1(pperp,cwwa)
!2024-12-21e

! calculate density (per volume)

!      do k = 1, lphi
!      do j = 1, lz
!      do i = 1, lr
!$omp parallel do private(vol)
      do i = 1, lrzphi
        vol = 1.0d0/(grr(i,1,1)*dr*dz*dphi)
        dns(i,1,1) = dns(i,1,1)*vol
        mom(i,1,1) = mom(i,1,1)*vol
        ppara(i,1,1) = ppara(i,1,1)*vol/mass
        pperp(i,1,1) = pperp(i,1,1)*vol*babs(i,1,1)
      end do
!      end do
!      end do


! wall effect

      lr1 = lr - 1
      lz1 = lz - 1

      if(my_rank_r.eq.(mpi_proc_r-1))then
!      do k = 1, lphi
!      do j = 1, lz
!$omp parallel do
       do j = 1, lzphi
          dns(lr1,j,1) = dns(lr1,j,1) + dns(lr,j,1)
          mom(lr1,j,1) = mom(lr1,j,1) + mom(lr,j,1)
          dns(lr,j,1) = 0.0d0
          mom(lr,j,1) = 0.0d0

          ppara(lr1,j,1) = ppara(lr1,j,1) + ppara(lr,j,1)
          pperp(lr1,j,1) = pperp(lr1,j,1) + pperp(lr,j,1)
          ppara(lr,j,1) = 0.0d0
          pperp(lr,j,1) = 0.0d0
       end do
!      end do
      end if

      if(my_rank_r.eq.0)then
!      do k = 1, lphi
!      do j = 1, lz
!$omp parallel do
       do j = 1, lzphi
          dns(2,j,1) = dns(2,j,1) + dns(1,j,1)
          mom(2,j,1) = mom(2,j,1) + mom(1,j,1)
          dns(1,j,1) = 0.0d0
          mom(1,j,1) = 0.0d0

          ppara(2,j,1) = ppara(2,j,1) + ppara(1,j,1)
          pperp(2,j,1) = pperp(2,j,1) + pperp(1,j,1)
          ppara(1,j,1) = 0.0d0
          pperp(1,j,1) = 0.0d0
       end do
!      end do
      end if

      if(my_rank_z.eq.(mpi_proc_z-1))then
!$omp parallel do
       do k = 1, lphi
       do i = 1, lr
          dns(i,lz1,k) = dns(i,lz1,k) + dns(i,lz,k)
          mom(i,lz1,k) = mom(i,lz1,k) + mom(i,lz,k)
          dns(i,lz,k) = 0.0d0
          mom(i,lz,k) = 0.0d0

          ppara(i,lz1,k) = ppara(i,lz1,k) + ppara(i,lz,k)
          pperp(i,lz1,k) = pperp(i,lz1,k) + pperp(i,lz,k)
          ppara(i,lz,k) = 0.0d0
          pperp(i,lz,k) = 0.0d0
       end do
       end do
      end if 

      if(my_rank_z.eq.0)then
!$omp parallel do
       do k = 1, lphi
       do i = 1, lr
          dns(i,2,k) = dns(i,2,k) + dns(i,1,k)
          mom(i,2,k) = mom(i,2,k) + mom(i,1,k)
          dns(i,1,k) = 0.0d0
          mom(i,1,k) = 0.0d0

          ppara(i,2,k) = ppara(i,2,k) + ppara(i,1,k)
          pperp(i,2,k) = pperp(i,2,k) + pperp(i,1,k)
          ppara(i,1,k) = 0.0d0
          pperp(i,1,k) = 0.0d0
       end do
       end do
      end if 

! take only n=1 modes
!      n = 1
!      call n1(n,ppara)
!      call n1(n,pperp)

! take modes n_min <= n <=n_max
!      n_min = 1
!      n_max = 2 
!      call lowpass2_mlt2(n_min,n_max,dns,mom)
!      call lowpass2_mlt2(n_min,n_max,ppara,pperp)

!      do k = 1, lphi
!      do j = 1, lz
!      do i = 1, lr
!$omp parallel do
      do i = 1, lrzphi
        dns(i,1,1) = dns(i,1,1) + dns0(i,1,1)
        mom(i,1,1) = mom(i,1,1) + mom0(i,1,1)
        ppara(i,1,1) = ppara(i,1,1) + ppara0(i,1,1)
        pperp(i,1,1) = pperp(i,1,1) + pperp0(i,1,1)/babs0(i,1,1)*babs(i,1,1)
      end do
!      end do
!      end do

end subroutine density_gyro
!--------------------------------------------------------------------
subroutine moments_gyro(marker_num,marker_num_gyro,mass,gc,gyro_phys &
                  ,dns,mom,ppara,pperp,qpara,qperp &
                  ,dns0,mom0,ppara0,pperp0,qpara0,qperp0)
! 2016-02-04, for FLR case
! modified for NEC SX-Aurora TSUBASA 2020-07-10
! gc_x(:,n) = 1:r, 2:z, 3:phi, 4:p, 5:mu, 6:weight, 7:active, 8:p_phi, 9:fff, 10:fnrml
!--------------------------------------------------------------------
      use mpiset
      use grid
      use field, only:babs,babs0
      use gyro, only:ngyro,marker_each_gyro
      implicit none

      integer::marker_num,nvec,marker_num_gyro
      real(8)::gc(ngc2,marker_num)
      real(8)::gyro_phys(2,marker_num_gyro)
      real(8)::mass
      real(8)::dns(lr,lz,lphi),mom(lr,lz,lphi)
      real(8)::ppara(lr,lz,lphi),pperp(lr,lz,lphi)
      real(8)::qpara(lr,lz,lphi),qperp(lr,lz,lphi)
      real(8)::dns0(lr,lz,lphi),mom0(lr,lz,lphi)
      real(8)::ppara0(lr,lz,lphi),pperp0(lr,lz,lphi)
      real(8)::qpara0(lr,lz,lphi),qperp0(lr,lz,lphi)
      integer::i,j,k,l,m,n,ia,ia1,ja,ja1,ka,ka1,lr1,lz1
      real(8)::ar,ar1,az,az1,aphi,aphi1
      real(8)::aaa1,aaa2,aaa3,aaa4,aaa5,aaa6,aaa7,aaa8
      real(8)::vol,cwwa,cwwc
!      real(8)::d0,d1,d2,d3,d4,d5,t1,t2,t3,t4,t5
      real(8)::dr1,dz1,dphi1,ma_mi_r
      integer::kr,kz,kphi
      real(8)::p0,p1,p2,mu1,p3,p1mu,mass1 !2015-09-21
      integer::n_min,n_max !2013-05-22
      real(8)::rngr1
      integer::nc
      integer::nvec0,nmod
!NEC start
!      integer,parameter::nblkd=255 !20NOV
!      integer,parameter::nblkd=127 !21JAN
      integer,parameter::nblkd=15
      real(8)::wkdns(nblkd,lr,lz,lphi,lpara),wkmom(nblkd,lr,lz,lphi,lpara)
      real(8)::wkpar(nblkd,lr,lz,lphi,lpara),wkprp(nblkd,lr,lz,lphi,lpara)
      real(8)::wkqar(nblkd,lr,lz,lphi,lpara),wkqrp(nblkd,lr,lz,lphi,lpara)

        do l = 1, lpara
!$omp parallel do
          do i = 1, lrzphi
            do k = 1, nblkd
              wkdns(k,i,1,1,l)=0.0d0
              wkmom(k,i,1,1,l)=0.0d0
              wkpar(k,i,1,1,l)=0.0d0
              wkprp(k,i,1,1,l)=0.0d0
              wkqar(k,i,1,1,l)=0.0d0
              wkqrp(k,i,1,1,l)=0.0d0
            end do
          end do
        end do
!NEC end

      mass1 = 1.0d0/mass

      dr1 = 1.0d0/dr
      dz1 = 1.0d0/dz
      dphi1 = 1.0d0/dphi
      rngr1 = 1.0d0/dble(ngyro)

      kr   = 1 - kr_offset(my_rank)
      kz   = 1 - kz_offset(my_rank)
      kphi = 1 + lphi_shd - kphi_offset(my_rank)
!      kphi = 3 - kphi_offset(my_rank)
      ma_mi_r = major_r-minor_r

!$omp parallel do
        do i = 1, lrzphi
          dns(i,1,1) = 0.0d0
          mom(i,1,1) = 0.0d0
          ppara(i,1,1) = 0.0d0
          pperp(i,1,1) = 0.0d0
          qpara(i,1,1) = 0.0d0
          qperp(i,1,1) = 0.0d0
        end do


      nvec = marker_num_gyro/lpara
      nvec0 = marker_num_gyro/lpara
      nmod  = mod(marker_num_gyro,lpara)

!$omp parallel private(ar1,ar,az1,az,aphi1,aphi &
!$omp& ,nvec,m,k,n,nc,ia,ja,ka,p0,p1,p2,mu1,p3,p1mu &
!$omp& ,aaa1,aaa2,aaa3,aaa4,aaa5,aaa6,aaa7,aaa8)
!$omp do
      do l = 1, lpara
        nvec = nvec0 + min(1,nmod/l)

      do m = 1, nvec, nblkd

       do k = 1, min(nblkd, nvec-m+1)
        n = m + k - 1 + nvec0*(l-1) + min(l-1,nmod)
        nc =(n-1)/ngyro + 1

! gc_x(:,n) = 1:r, 2:z, 3:phi, 4:p, 5:mu, 6:weight, 7:active, 8:p_phi, 9:fff, 10:fnrml
        ia=max(1,min(lr  -1,int((gyro_phys(1,n)-ma_mi_r)*dr1  ) + kr  ))
        ja=max(1,min(lz  -1,int(gyro_phys(2,n)          *dz1  ) + kz  ))
        ka=max(1,min(lphi-1,int(gc(3,nc)        *dphi1) + kphi))

        ar1 = (gyro_phys(1,n)-grr(ia,ja,ka ) ) *dr1
        ar  = 1.0d0 - ar1
        az1 = (gyro_phys(2,n)-gzz(ia,ja,ka ) ) *dz1
        az  = 1.0d0 - az1
        aphi1 = (gc(3,nc)-gphi(ia,ja,ka ) ) *dphi1
        aphi  = 1.0d0 - aphi1

        aaa1 = ar *az *aphi
        aaa2 = ar1*az *aphi
        aaa3 = ar *az1*aphi
        aaa4 = ar1*az1*aphi
        aaa5 = ar *az *aphi1
        aaa6 = ar1*az *aphi1
        aaa7 = ar *az1*aphi1
        aaa8 = ar1*az1*aphi1

        p0 =  gc(6,nc)*gc(7,nc)*rngr1 !rngr1 is multiplied
        p1 =  gc(4,nc)*p0
        p2 =  gc(4,nc)*p1
        mu1 = gc(5,nc)*p0
        p3  = gc(4,nc)*p2
        p1mu= gc(4,nc)*mu1

        wkdns(k,ia,  ja,  ka,  l) = wkdns(k,ia,  ja,  ka,  l) + aaa1*p0
        wkdns(k,ia+1,ja,  ka,  l) = wkdns(k,ia+1,ja,  ka,  l) + aaa2*p0
        wkdns(k,ia,  ja+1,ka,  l) = wkdns(k,ia,  ja+1,ka,  l) + aaa3*p0
        wkdns(k,ia+1,ja+1,ka,  l) = wkdns(k,ia+1,ja+1,ka,  l) + aaa4*p0
        wkdns(k,ia,  ja,  ka+1,l) = wkdns(k,ia,  ja,  ka+1,l) + aaa5*p0
        wkdns(k,ia+1,ja,  ka+1,l) = wkdns(k,ia+1,ja,  ka+1,l) + aaa6*p0
        wkdns(k,ia,  ja+1,ka+1,l) = wkdns(k,ia,  ja+1,ka+1,l) + aaa7*p0
        wkdns(k,ia+1,ja+1,ka+1,l) = wkdns(k,ia+1,ja+1,ka+1,l) + aaa8*p0

        wkmom(k,ia,  ja,  ka,  l) = wkmom(k,ia,  ja,  ka,  l) + aaa1*p1
        wkmom(k,ia+1,ja,  ka,  l) = wkmom(k,ia+1,ja,  ka,  l) + aaa2*p1
        wkmom(k,ia,  ja+1,ka,  l) = wkmom(k,ia,  ja+1,ka,  l) + aaa3*p1
        wkmom(k,ia+1,ja+1,ka,  l) = wkmom(k,ia+1,ja+1,ka,  l) + aaa4*p1
        wkmom(k,ia,  ja,  ka+1,l) = wkmom(k,ia,  ja,  ka+1,l) + aaa5*p1
        wkmom(k,ia+1,ja,  ka+1,l) = wkmom(k,ia+1,ja,  ka+1,l) + aaa6*p1
        wkmom(k,ia,  ja+1,ka+1,l) = wkmom(k,ia,  ja+1,ka+1,l) + aaa7*p1
        wkmom(k,ia+1,ja+1,ka+1,l) = wkmom(k,ia+1,ja+1,ka+1,l) + aaa8*p1

        wkpar(k,ia,  ja,  ka,  l) = wkpar(k,ia,  ja,  ka,  l) + aaa1*p2
        wkpar(k,ia+1,ja,  ka,  l) = wkpar(k,ia+1,ja,  ka,  l) + aaa2*p2
        wkpar(k,ia,  ja+1,ka,  l) = wkpar(k,ia,  ja+1,ka,  l) + aaa3*p2
        wkpar(k,ia+1,ja+1,ka,  l) = wkpar(k,ia+1,ja+1,ka,  l) + aaa4*p2
        wkpar(k,ia,  ja,  ka+1,l) = wkpar(k,ia,  ja,  ka+1,l) + aaa5*p2
        wkpar(k,ia+1,ja,  ka+1,l) = wkpar(k,ia+1,ja,  ka+1,l) + aaa6*p2
        wkpar(k,ia,  ja+1,ka+1,l) = wkpar(k,ia,  ja+1,ka+1,l) + aaa7*p2
        wkpar(k,ia+1,ja+1,ka+1,l) = wkpar(k,ia+1,ja+1,ka+1,l) + aaa8*p2

        wkprp(k,ia,  ja,  ka,  l) = wkprp(k,ia,  ja,  ka,  l) + aaa1*mu1
        wkprp(k,ia+1,ja,  ka,  l) = wkprp(k,ia+1,ja,  ka,  l) + aaa2*mu1
        wkprp(k,ia,  ja+1,ka,  l) = wkprp(k,ia,  ja+1,ka,  l) + aaa3*mu1
        wkprp(k,ia+1,ja+1,ka,  l) = wkprp(k,ia+1,ja+1,ka,  l) + aaa4*mu1
        wkprp(k,ia,  ja,  ka+1,l) = wkprp(k,ia,  ja,  ka+1,l) + aaa5*mu1
        wkprp(k,ia+1,ja,  ka+1,l) = wkprp(k,ia+1,ja,  ka+1,l) + aaa6*mu1
        wkprp(k,ia,  ja+1,ka+1,l) = wkprp(k,ia,  ja+1,ka+1,l) + aaa7*mu1
        wkprp(k,ia+1,ja+1,ka+1,l) = wkprp(k,ia+1,ja+1,ka+1,l) + aaa8*mu1

        wkqar(k,ia,  ja,  ka,  l) = wkqar(k,ia,  ja,  ka,  l) + aaa1*p3
        wkqar(k,ia+1,ja,  ka,  l) = wkqar(k,ia+1,ja,  ka,  l) + aaa2*p3
        wkqar(k,ia,  ja+1,ka,  l) = wkqar(k,ia,  ja+1,ka,  l) + aaa3*p3
        wkqar(k,ia+1,ja+1,ka,  l) = wkqar(k,ia+1,ja+1,ka,  l) + aaa4*p3
        wkqar(k,ia,  ja,  ka+1,l) = wkqar(k,ia,  ja,  ka+1,l) + aaa5*p3
        wkqar(k,ia+1,ja,  ka+1,l) = wkqar(k,ia+1,ja,  ka+1,l) + aaa6*p3
        wkqar(k,ia,  ja+1,ka+1,l) = wkqar(k,ia,  ja+1,ka+1,l) + aaa7*p3
        wkqar(k,ia+1,ja+1,ka+1,l) = wkqar(k,ia+1,ja+1,ka+1,l) + aaa8*p3

        wkqrp(k,ia,  ja,  ka,  l) = wkqrp(k,ia,  ja,  ka,  l) + aaa1*p1mu
        wkqrp(k,ia+1,ja,  ka,  l) = wkqrp(k,ia+1,ja,  ka,  l) + aaa2*p1mu
        wkqrp(k,ia,  ja+1,ka,  l) = wkqrp(k,ia,  ja+1,ka,  l) + aaa3*p1mu
        wkqrp(k,ia+1,ja+1,ka,  l) = wkqrp(k,ia+1,ja+1,ka,  l) + aaa4*p1mu
        wkqrp(k,ia,  ja,  ka+1,l) = wkqrp(k,ia,  ja,  ka+1,l) + aaa5*p1mu
        wkqrp(k,ia+1,ja,  ka+1,l) = wkqrp(k,ia+1,ja,  ka+1,l) + aaa6*p1mu
        wkqrp(k,ia,  ja+1,ka+1,l) = wkqrp(k,ia,  ja+1,ka+1,l) + aaa7*p1mu
        wkqrp(k,ia+1,ja+1,ka+1,l) = wkqrp(k,ia+1,ja+1,ka+1,l) + aaa8*p1mu

       end do
      end do
      end do
!$omp end parallel

        do l = 1, lpara
!$omp parallel do
          do i = 1, lrzphi
            do k = 1, nblkd
              dns(i,1,1) = dns(i,1,1) + wkdns(k,i,1,1,l)
              mom(i,1,1) = mom(i,1,1) + wkmom(k,i,1,1,l)
              ppara(i,1,1) = ppara(i,1,1) + wkpar(k,i,1,1,l)
              pperp(i,1,1) = pperp(i,1,1) + wkprp(k,i,1,1,l)
              qpara(i,1,1) = qpara(i,1,1) + wkqar(k,i,1,1,l)
              qperp(i,1,1) = qperp(i,1,1) + wkqrp(k,i,1,1,l)
            end do
          end do
        end do

! smoothing

       cwwa = 0.5d0
!       cwwc =-1.d0/6.d0 

!2024-12-21s, correction suggested by Panith Adulsiriswad
       call periodic_particle_mlt6b(dns,mom,ppara,pperp,qpara,qperp)
       call partsm1(dns,cwwa)
       call partsm1(mom,cwwa)

       call partsm1(ppara,cwwa)
       call partsm1(pperp,cwwa)

       call partsm1(qpara,cwwa)
       call partsm1(qperp,cwwa)
!2024-12-21e

! calculate density (per volume)

!      do k = 1, lphi
!      do j = 1, lz
!      do i = 1, lr
!$omp parallel do private(vol)
      do i = 1, lrzphi
        vol = 1.0d0/(grr(i,1,1)*dr*dz*dphi)
        dns(i,1,1) = dns(i,1,1)*vol
        mom(i,1,1) = mom(i,1,1)*vol
        ppara(i,1,1) = ppara(i,1,1)*vol*mass1
        pperp(i,1,1) = pperp(i,1,1)*vol*babs(i,1,1)
        qpara(i,1,1) = qpara(i,1,1)*vol*mass1**2
        qperp(i,1,1) = qperp(i,1,1)*vol*babs(i,1,1)*mass1
      end do
!      end do
!      end do



! wall effect

      lr1 = lr - 1
      lz1 = lz - 1

      if(my_rank_r.eq.(mpi_proc_r-1))then
!      do k = 1, lphi
!      do j = 1, lz
!$omp parallel do
       do j = 1, lzphi
          dns(lr1,j,1) = dns(lr1,j,1) + dns(lr,j,1)
          mom(lr1,j,1) = mom(lr1,j,1) + mom(lr,j,1)
          dns(lr,j,1) = 0.0d0
          mom(lr,j,1) = 0.0d0

          ppara(lr1,j,1) = ppara(lr1,j,1) + ppara(lr,j,1)
          pperp(lr1,j,1) = pperp(lr1,j,1) + pperp(lr,j,1)
          ppara(lr,j,1) = 0.0d0
          pperp(lr,j,1) = 0.0d0

          qpara(lr1,j,1) = qpara(lr1,j,1) + qpara(lr,j,1)
          qperp(lr1,j,1) = qperp(lr1,j,1) + qperp(lr,j,1)
          qpara(lr,j,1) = 0.0d0
          qperp(lr,j,1) = 0.0d0
       end do
!      end do
      end if

      if(my_rank_r.eq.0)then
!      do k = 1, lphi
!      do j = 1, lz
!$omp parallel do
       do j = 1, lzphi
          dns(2,j,1) = dns(2,j,1) + dns(1,j,1)
          mom(2,j,1) = mom(2,j,1) + mom(1,j,1)
          dns(1,j,1) = 0.0d0
          mom(1,j,1) = 0.0d0

          ppara(2,j,1) = ppara(2,j,1) + ppara(1,j,1)
          pperp(2,j,1) = pperp(2,j,1) + pperp(1,j,1)
          ppara(1,j,1) = 0.0d0
          pperp(1,j,1) = 0.0d0

          qpara(2,j,1) = qpara(2,j,1) + qpara(1,j,1)
          qperp(2,j,1) = qperp(2,j,1) + qperp(1,j,1)
          qpara(1,j,1) = 0.0d0
          qperp(1,j,1) = 0.0d0
       end do
!      end do
      end if

      if(my_rank_z.eq.(mpi_proc_z-1))then
!$omp parallel do
       do k = 1, lphi
       do i = 1, lr
          dns(i,lz1,k) = dns(i,lz1,k) + dns(i,lz,k)
          mom(i,lz1,k) = mom(i,lz1,k) + mom(i,lz,k)
          dns(i,lz,k) = 0.0d0
          mom(i,lz,k) = 0.0d0

          ppara(i,lz1,k) = ppara(i,lz1,k) + ppara(i,lz,k)
          pperp(i,lz1,k) = pperp(i,lz1,k) + pperp(i,lz,k)
          ppara(i,lz,k) = 0.0d0
          pperp(i,lz,k) = 0.0d0

          qpara(i,lz1,k) = qpara(i,lz1,k) + qpara(i,lz,k)
          qperp(i,lz1,k) = qperp(i,lz1,k) + qperp(i,lz,k)
          qpara(i,lz,k) = 0.0d0
          qperp(i,lz,k) = 0.0d0
       end do
       end do
      end if 

      if(my_rank_z.eq.0)then
!$omp parallel do
       do k = 1, lphi
       do i = 1, lr
          dns(i,2,k) = dns(i,2,k) + dns(i,1,k)
          mom(i,2,k) = mom(i,2,k) + mom(i,1,k)
          dns(i,1,k) = 0.0d0
          mom(i,1,k) = 0.0d0

          ppara(i,2,k) = ppara(i,2,k) + ppara(i,1,k)
          pperp(i,2,k) = pperp(i,2,k) + pperp(i,1,k)
          ppara(i,1,k) = 0.0d0
          pperp(i,1,k) = 0.0d0

          qpara(i,2,k) = qpara(i,2,k) + qpara(i,1,k)
          qperp(i,2,k) = qperp(i,2,k) + qperp(i,1,k)
          qpara(i,1,k) = 0.0d0
          qperp(i,1,k) = 0.0d0
       end do
       end do
      end if 

!      if(.not.flag_stored)then
! take modes n_min <= n <=n_max
!      n_min = 1
!      n_max = 2
!      call lowpass2_mlt2(n_min,n_max,dns,mom)
!      call lowpass2_mlt2(n_min,n_max,ppara,pperp)
!      call lowpass2_mlt2(n_min,n_max,qpara,qperp)
!      end if

!2014-08-28s
!$omp parallel do
      do i = 1, lrzphi
        dns(i,1,1) = dns(i,1,1) + dns0(i,1,1)
        mom(i,1,1) = mom(i,1,1) + mom0(i,1,1)
        ppara(i,1,1) = ppara(i,1,1) + ppara0(i,1,1)
        pperp(i,1,1) = pperp(i,1,1) + pperp0(i,1,1)/babs0(i,1,1)*babs(i,1,1)
        qpara(i,1,1) = qpara(i,1,1) + qpara0(i,1,1)
        qperp(i,1,1) = qperp(i,1,1) + qperp0(i,1,1)/babs0(i,1,1)*babs(i,1,1)
      end do
!2014-08-28e

!            call wall_clock(d5)
!            t5 = d5 - d4
!            if(my_rank.eq.0)then
!            write(6,*)'density 5, t5=',t5
!            end if

end subroutine moments_gyro

end module aurora
!----------------------------------------------------------------------------


!----------------------------------------------------------------------------
module fx100
contains
!--------------------------------------------------------------------
subroutine push(marker_num,amassp,ep &
               ,type,temp,valpha,deltav,clambda0,dclambda & !2012-06-17
               ,gc,dgc,v &
               ,cf_pphi,pphi_min,pphi_max &
               ,flp,ispecies,igyro)
! type=0: maxwellian, type=1: slowing down, type=2: beam
! type=3: beam with finite pitch angle width !2012-06-17
! igyro=0: w/o FLR, igyro=1: w/ FLR
! modified on 2015-06-23
!   gc_x(:,n) = 1:r, 2:z, 3:phi, 4:p, 5:mu, 6:weight, 7:active, 8:p_phi, 9:fff, 10:fnrml
!--------------------------------------------------------------------
      use mpiset
      use field
      use grid
!      use equi_sol, only:raxis
      use particle, only:nu_krook !2024-04-23
      implicit none

      integer::marker_num,type
      integer::ispecies,igyro
      real(8)::amassp,ep
      real(8)::gc(ngc2,marker_num)
      real(8)::dgc(ngc1,marker_num)
      real(8)::v(marker_num)
      real(8)::cf_pphi(0:lcfpphi),pphi_min,pphi_max !2016-02-04
      real(8)::temp,valpha,deltav
      real(8)::dr1,dz1,dphi1,ep1,amsp1
      integer::n,ia,ia1,ja,ja1,ka,ka1,i
      real(8)::ar,ar1,az,az1,aphi,aphi1
!      real(8)::aaa1,aaa2,aaa3,aaa4,aaa5,aaa6,aaa7,aaa8

!2015-06-04s
      integer::ijk_a(3,marker_each)
      real(8)::aaa(8,marker_each)
!2015-06-04e

      real(8)::flp(nflp,marker_num),detotal
      real(8)::b2,babse,babs0e,b21,bre,bze,bphie
      real(8)::b1,br1a,bz1a,bphi1a,b10,br10,bz10,bphi10,dbpar,dbpr,dbpz,dbpphi
      real(8)::rhopar,orbpr,orbpz,orbpphi,denom1
      real(8)::wvnmlr,wvnmlz,wvnmlphi
      real(8)::dgrdbr,dgrdbz,dgrdbphi,dppar,dppar1,dppar2 !2012-08-31
      real(8)::pvpar,w1nmlr,w1nmlz,w1nmlphi
      real(8)::psinrm,pactive,prof,dprofdpsi,energy0,pphi_n,rminor,bmax
      real(8)::dpphi_n,vpara,vpara0,dvpara,sigma
      real(8)::vlcpart,dvlcpart,vnrm,sd_factor
      real(8)::sqr_pi1
      real(8)::nudt,coef,vrat !2016-01-09
      real(8)::nu_krook_dt !2024-04-23

      integer::kr,kz,kphi
      real(8)::ma_mi_r
      integer::kfl_start,kfl

      real(8)::pt_factor,energy,clambda0,dclambda,clambda(marker_each) !2012-06-17
      real(8)::dwpsi(marker_each),dwenrc(marker_each) !2015-06-23
!      real(8)::psip(marker_each)
!      real(8)::dpsi_dr(marker_each),dpsi_dz(marker_each),dpsi_dphi(marker_each)
      real(8)::bphi1e(marker_each) !2016-01-09

! time derivative of particle position and velocity

      dr1 = 1.0d0/dr
      dz1 = 1.0d0/dz
      dphi1 = 1.0d0/dphi
      ep1 = 1.0d0/ep
      amsp1 = 1.0d0/amassp
      sqr_pi1=1.0d0/sqrt(pi)

      nu_krook_dt = nu_krook * dt !2024-04-23

      kr   = 1 - kr_offset(my_rank)
      kz   = 1 - kz_offset(my_rank)
      kphi = 1 + lphi_shd - kphi_offset(my_rank)
      ma_mi_r = major_r-minor_r

! igyro=0: w/o FLR, igyro/=0: flp(1:7,:) is given in subr. extract_em
      if(igyro.eq.0)then
        kfl_start = 1
      else
        kfl_start = nflp_gyro + 1
      end if

!$omp parallel private(ar1,ar,az1,az,aphi1,aphi)

      if(igyro.eq.0)then !w/o FLR
!$omp workshare
        fld(1,:,:,:) = er(:,:,:)
        fld(2,:,:,:) = ez(:,:,:)
        fld(3,:,:,:) = ephi(:,:,:)
        fld(4,:,:,:) = epara(:,:,:)
        fld(5,:,:,:) = br(:,:,:) - br0(:,:,:)
        fld(6,:,:,:) = bz(:,:,:) - bz0(:,:,:)
        fld(7,:,:,:) = bphi(:,:,:) - bphi0(:,:,:)
!$omp end workshare
      end if

!2025-02-05      
!$omp workshare
      fld(11,:,:,:)= gradbr(:,:,:)
      fld(12,:,:,:)= gradbz(:,:,:)
      fld(13,:,:,:)= gradbphi(:,:,:)
      fld(14,:,:,:)= curvbr(:,:,:)
      fld(15,:,:,:)= curvbz(:,:,:)
      fld(16,:,:,:)= curvbphi(:,:,:)
!$omp end workshare

      if(ispecies.eq.0)then
!$omp workshare
        fld(23,:,:,:)= dns_e0(:,:,:)
        fld(24,:,:,:)= dns_e0_r(:,:,:)
        fld(25,:,:,:)= dns_e0_z(:,:,:)
        fld(26,:,:,:)= dns_e0_phi(:,:,:)
        fld(27,:,:,:)= temp_e0(:,:,:)
        fld(28,:,:,:)= temp_e0_r(:,:,:)
        fld(29,:,:,:)= temp_e0_z(:,:,:)
        fld(30,:,:,:)= temp_e0_phi(:,:,:)
!$omp end workshare
      else if(ispecies.eq.1)then
!$omp workshare
        fld(23,:,:,:)= dns_i0(:,:,:)
        fld(24,:,:,:)= dns_i0_r(:,:,:)
        fld(25,:,:,:)= dns_i0_z(:,:,:)
        fld(26,:,:,:)= dns_i0_phi(:,:,:)
        fld(27,:,:,:)= temp_i0(:,:,:)
        fld(28,:,:,:)= temp_i0_r(:,:,:)
        fld(29,:,:,:)= temp_i0_z(:,:,:)
        fld(30,:,:,:)= temp_i0_phi(:,:,:)
!$omp end workshare
      else if(ispecies.eq.2)then
!$omp workshare
        fld(23,:,:,:)= dns_a0(:,:,:)
        fld(24,:,:,:)= dns_a0_r(:,:,:)
        fld(25,:,:,:)= dns_a0_z(:,:,:)
        fld(26,:,:,:)= dns_a0_phi(:,:,:)
        fld(27,:,:,:)= temp_a0(:,:,:)
        fld(28,:,:,:)= temp_a0_r(:,:,:)
        fld(29,:,:,:)= temp_a0_z(:,:,:)
        fld(30,:,:,:)= temp_a0_phi(:,:,:)
!$omp end workshare
      end if


!poption parallel,schedule(dynamic,1000)
!$omp do schedule(dynamic,1000)
      do n = 1, marker_num
        ijk_a(1,n)=max(1,min(lr  -1,int((gc(1,n)-ma_mi_r)*dr1  ) + kr  ))
        ijk_a(2,n)=max(1,min(lz  -1,int(gc(2,n)          *dz1  ) + kz  ))
        ijk_a(3,n)=max(1,min(lphi-1,int(gc(3,n)          *dphi1) + kphi))

        ar1  = max(0.0d0, min(1.0d0, (gc(1,n) - ma_mi_r)*dr1 -  dble(ijk_a(1,n) - kr)  ) )
        ar   = 1.0d0 - ar1
        az1  = max(0.0d0, min(1.0d0, gc(2,n)*dz1 - dble(ijk_a(2,n) - kz) ) )
        az   = 1.0d0 - az1
        aphi1 = max(0.0d0, min(1.0d0, gc(3,n)*dphi1 - dble(ijk_a(3,n) - kphi) ) )
        aphi = 1.0d0 - aphi1

        aaa(1,n) = ar *az *aphi
        aaa(2,n) = ar1*az *aphi
        aaa(3,n) = ar *az1*aphi
        aaa(4,n) = ar1*az1*aphi
        aaa(5,n) = ar *az *aphi1
        aaa(6,n) = ar1*az *aphi1
        aaa(7,n) = ar *az1*aphi1
        aaa(8,n) = ar1*az1*aphi1
      end do
!$omp end parallel

!poption parallel,schedule(dynamic,1000)
!$omp parallel private(ia,ja,ka)
!$omp do schedule(dynamic,1000)
      do n = 1, marker_num
        ia=ijk_a(1,n)
        ja=ijk_a(2,n)
        ka=ijk_a(3,n)

! fields at each particle position
! igyro=0: w/o FLR, igyro/=0: flp(1:7,:) is given in subr. extract_em
        do kfl = kfl_start, nflp

        flp(kfl,n) = fld(kfl, ia, ja,  ka  )*aaa(1,n) + fld(kfl, ia+1,ja,  ka  )*aaa(2,n) &
                   + fld(kfl, ia, ja+1,ka  )*aaa(3,n) + fld(kfl, ia+1,ja+1,ka  )*aaa(4,n) &
                   + fld(kfl, ia, ja,  ka+1)*aaa(5,n) + fld(kfl, ia+1,ja,  ka+1)*aaa(6,n) &
                   + fld(kfl, ia, ja+1,ka+1)*aaa(7,n) + fld(kfl, ia+1,ja+1,ka+1)*aaa(8,n)
        end do

      end do
!$omp end parallel


!$omp parallel private(bre,bze,bphie,b2,babse,babs0e,b1,b21,br1a,bz1a,bphi1a &
!$omp& ,b10,br10,bz10,bphi10,dbpar,dbpr,dbpz,dbpphi,rhopar,denom1,orbpr,orbpz,orbpphi &
!$omp& ,wvnmlr,wvnmlz,wvnmlphi,dgrdbr,dgrdbz,dgrdbphi,dppar,dppar2,dppar1,pvpar &
!$omp& ,w1nmlr,w1nmlz,w1nmlphi,detotal,energy &
!$omp& ,rminor,bmax,energy0,sigma,vpara0,vpara,dvpara &
!$omp& ,pphi_n,dpphi_n,prof,dprofdpsi)
!$omp do
      do n = 1, marker_num
! flp(5:7,n): delta_br(z,phi)
        bre = flp(5,n) + flp(8,n)
        bze = flp(6,n) + flp(9,n)
        bphie = flp(7,n) + flp(10,n)

        b2  = bre**2 + bze**2 + bphie**2
        babse= max(eps_b, sqrt(b2) )
        babs0e= max(eps_b, sqrt(flp(8,n)**2 + flp(9,n)**2 + flp(10,n)**2) )
        b1 = 1.0d0/babse
        b21= 1.0d0/b2
        br1a = bre*b1
        bz1a = bze*b1
        bphi1a = bphie*b1
        bphi1e(n) = bphi1a !2016-02-4

        b10 = 1.0d0/babs0e
        br10 = flp(8,n)*b10
        bz10 = flp(9,n)*b10
        bphi10 = flp(10,n)*b10

        dbpar = br1a*br10 + bz1a*bz10 + bphi1a*bphi10
        dbpr  = br1a  - br10*dbpar
        dbpz  = bz1a  - bz10*dbpar
        dbpphi= bphi1a- bphi10*dbpar

! guiding center motion

        rhopar = gc(4,n)*ep1*b1

        denom1 = 1.0d0/(1.0d0 + rhopar*(br1a*flp(14,n) + bz1a*flp(15,n) + bphi1a*flp(16,n)))

        orbpr = (br1a + rhopar*flp(14,n))*denom1
        orbpz = (bz1a + rhopar*flp(15,n))*denom1
        orbpphi = (bphi1a + rhopar*flp(16,n))*denom1

! e x b drift

        wvnmlr = (flp(3,n)*bze-flp(2,n)*bphie)*b21*denom1
        wvnmlz = (flp(1,n)*bphie-flp(3,n)*bre)*b21*denom1
        wvnmlphi = (flp(2,n)*bre-flp(1,n)*bze)*b21*denom1

! grad-b drift

        dgrdbr = gc(5,n)*(bphie*flp(12,n) - bze*flp(13,n))*ep1*b21*denom1
        dgrdbz = gc(5,n)*(bre*flp(13,n) - bphie*flp(11,n))*ep1*b21*denom1
        dgrdbphi = gc(5,n)*(bze*flp(11,n) - bre*flp(12,n))*ep1*b21*denom1

! mirror force

        dppar =-gc(5,n)*( flp(11,n)*orbpr &
                      + flp(12,n)*orbpz &
                      + flp(13,n)*orbpphi &
                      )*dt

!2012-08-31
        dppar2=-gc(5,n)*( flp(11,n)*dbpr*denom1 &
                      + flp(12,n)*dbpz*denom1 &
                      + flp(13,n)*orbpphi &
                      )*dt
!2012-08-31 end

! aceeleration due to electric field and curvature drift

!2012-07-07
        dppar1 = ep*((flp(1,n)*flp(14,n) + flp(2,n)*flp(15,n) + flp(3,n)*flp(16,n) &
                      )*rhopar &
                    + flp(4,n) &
                     )*dt*denom1
!2012-07-07 end

! total drift velocity

        pvpar = gc(4,n)*amsp1
        dgc(1,n) = dt*(pvpar*orbpr + wvnmlr + dgrdbr)*gc(7,n)
        dgc(2,n) = dt*(pvpar*orbpz + wvnmlz + dgrdbz)*gc(7,n)
        dgc(3,n) = dt*(pvpar*orbpphi + wvnmlphi + dgrdbphi)/gc(1,n)*gc(7,n)
        dgc(4,n) =(dppar + dppar1)*gc(7,n)

! temporal evolution of weight of high-energy ion particle

        w1nmlr = wvnmlr + pvpar*dbpr*denom1
        w1nmlz = wvnmlz + pvpar*dbpz*denom1
        w1nmlphi = wvnmlphi + pvpar*dbpphi*denom1

        dgc(8,n) =(ep*(w1nmlr*flp(18,n) + w1nmlz*flp(19,n) + w1nmlphi*flp(20,n))*dt &
                  + bphi1a*(w1nmlr*gc(4,n)*dt + gc(1,n)*(dppar1 + dppar2) ) & !2012-08-31
                  )*gc(7,n)

        detotal =(  ep*(flp(1,n)*dgrdbr + flp(2,n)*dgrdbz + flp(3,n)*dgrdbphi)*dt &
                  + dppar1*pvpar )*gc(7,n) &
! the following term considers 'mu* v * grad(B0 - B)', 2025-04-04
                + gc(5,n)*(dgc(1,n)*(flp(31,n)-flp(11,n) ) &
                          +dgc(2,n)*(flp(32,n)-flp(12,n) ) &
                          +dgc(3,n)*(flp(33,n)-flp(13,n) )*gc(1,n) &
                          )
!2025-04-27s
        energy0 = max(1.0d-30, 0.50d0*amassp*pvpar**2 + gc(5,n)*babs0e)
        clambda(n) = gc(5,n)*b0/energy0
        v(n) = sqrt(2.0d0*energy0*amsp1)

!        energy = max(1.0d-30, 0.50d0*amassp*pvpar**2 + gc(5,n)*babse)
!        clambda(n) = gc(5,n)*b0/energy
!        v(n) = sqrt(2.0d0*energy*amsp1)
!2025-04-27e

! weight evolution : weight = f - f0
! d weight /dt = - d f0 /dt
! f0 = f_nrml*prof(psinrm)/(v**3 + flp(21,n)**3)*0.50*erfc((v(n)-valpha)/deltav)

!        rminor = sqrt( (gc(1,n)-raxis)**2 &
!                     + ((gc(2,n)-0.50d0*zleng) )**2 &
!                     )
!        bmax = b0*raxis/(raxis-rminor)
!        energy0 = 0.50d0*amsp1*gc(4,n)**2 + gc(5,n)*babs0e
!        sigma = 0.50d0*(1.0d0 + sign(1.0d0, energy0-gc(5,n)*bmax) ) &
!                      *sign(1.0d0,gc(4,n) )
!        vpara0 = sqrt(2.0d0*(energy0-gc(5,n)*bmin)*amsp1)
!        vpara = vpara0*sigma
!        dvpara = sigma/vpara0*gc(4,n)*dppar1*amsp1**2

!        pphi_n = gc(8,n) - amassp*raxis*vpara
!        dpphi_n= dgc(8,n) - amassp*raxis*dvpara

!        prof = exp(pphi_n/(ep*psimax*0.37d0) )
!        dprofdpsi = prof/(ep*psimax*0.37d0)

!        dwpsi(n) = dpphi_n*dprofdpsi*gc(10,n)
!        dwenrc(n)= detotal/(amassp*v(n) )*prof*gc(10,n)

!2016-08-05s
        dwpsi(n) = ( w1nmlr*flp(24,n) + w1nmlz*flp(25,n) + w1nmlphi*flp(26,n) &
                    +(w1nmlr*flp(28,n) + w1nmlz*flp(29,n) + w1nmlphi*flp(30,n) ) &
                    *0.50d0*(amassp*v(n)**2/flp(27,n) - 3.0d0)*flp(23,n)/flp(27,n) &
                   )*dt*gc(7,n)*gc(10,n)

!        dwpsi(n) = (w1nmlr*flp(24,n) + w1nmlz*flp(25,n) + w1nmlphi*flp(26,n))*dt & !2016-01-09
!                   *gc(7,n)*gc(10,n)
!2016-08-05e

        dwenrc(n)= detotal/(amassp*v(n) )*flp(23,n)*gc(10,n) !2016-01-09

      end do
!$omp end parallel


      if(type.eq.0)then

!$omp parallel do private(vlcpart,dvlcpart)
        do n = 1, marker_num
!2016-08-05s
          vlcpart = exp(-0.50d0*amassp*v(n)**2/flp(27,n) )*flp(27,n)**(-1.5d0) !2013-07-17
          dvlcpart = -amassp*v(n)/flp(27,n)*vlcpart 
!          vlcpart = exp(-0.50d0*amassp*v(n)**2/temp)
!          dvlcpart = -amassp*v(n)/temp*vlcpart 
!2016-08-05e
          dwpsi(n) = dwpsi(n)*vlcpart
          dwenrc(n) = dwenrc(n)*dvlcpart
          dgc(6,n) = - dwpsi(n) - dwenrc(n) &
                     - gc(6,n)*nu_krook_dt & !Krook operator at normalized psi > psi_edge, flp(17)=psi
                     *(0.50d0 + sign(0.50d0, 1.0d0 - flp(17,n)/psimax - psi_edge) ) !2024-04-23
        end do

      else if(type.eq.1.or.type.eq.2)then

!$omp parallel do private(vnrm,sd_factor,vlcpart,dvlcpart)
        do n = 1, marker_num
! for itpa benchmark Oct11, 2009

          vnrm =(v(n)-valpha)/deltav
          sd_factor = 1.0d0/(v(n)**3 + flp(21,n)**3)
          vlcpart = 0.50d0*erfc(vnrm)*sd_factor
          dvlcpart= -( 3.0d0*v(n)**2*vlcpart &
                     + exp(-vnrm**2)*sqr_pi1/deltav &
                     )*sd_factor

          dwpsi(n) = dwpsi(n)*vlcpart
          dwenrc(n) = dwenrc(n)*dvlcpart
          dgc(6,n) = - dwpsi(n) - dwenrc(n) &
                     - gc(6,n)*nu_krook_dt & !Krook operator at normalized psi > psi_edge, flp(17)=psi
                     *(0.50d0 + sign(0.50d0, 1.0d0 - flp(17,n)/psimax - psi_edge) ) !2024-04-23
        end do

!2012-06-17
      else if(type.eq.3)then

!$omp parallel do private(vnrm,sd_factor,pt_factor,vlcpart,dvlcpart)
        do n = 1, marker_num
          vnrm =(v(n)-valpha)/deltav
          sd_factor = 1.0d0/(v(n)**3 + flp(21,n)**3)
          pt_factor = exp(-(clambda(n)-clambda0)**2/dclambda**2)

          vlcpart = 0.50d0*erfc(vnrm)*sd_factor*pt_factor
          dvlcpart= -( 3.0d0*v(n)**2*vlcpart &
                     + exp(-vnrm**2)*sqr_pi1/deltav*pt_factor &
                     )*sd_factor &
                  + 4.0d0*vlcpart*clambda(n)*(clambda(n)-clambda0) &
                         /(v(n)*dclambda**2)

          dwpsi(n) = dwpsi(n)*vlcpart
          dwenrc(n) = dwenrc(n)*dvlcpart
          dgc(6,n) = - dwpsi(n) - dwenrc(n) &
                     - gc(6,n)*nu_krook_dt & !Krook operator at normalized psi > psi_edge, flp(17)=psi
                     *(0.50d0 + sign(0.50d0, 1.0d0 - flp(17,n)/psimax - psi_edge) ) !2024-04-23
        end do
!2012-06-17 end

      end if

! slowing down
      if(type.eq.-5)then
!$omp parallel do private(nudt,vrat,coef)
        do n = 1, marker_num
            dgc(6,n) = 0.0d0 !full-f
            nudt = flp(22,n)*dt*gc(7,n)
            vrat = flp(21,n)/v(n)
            coef =(1.0d0 + vrat**3)*nudt
            dgc(10,n) =-3.0d0*nudt*gc(10,n)
            dgc(5,n)= -2.0d0*gc(5,n)*coef
            dgc(4,n) = dgc(4,n)- gc(4,n)*coef
            dgc(8,n) = dgc(8,n) - gc(4,n)*coef*gc(1,n)*bphi1e(n)
        end do
      end if

end subroutine push
!--------------------------------------------------------------------
subroutine density(marker_num,mass,gc &
                  ,dns,mom,ppara,pperp,dns0,mom0,ppara0,pperp0)
! calculate pressure
! modified on 2015-06-26
! gc_x(:,n) = 1:r, 2:z, 3:phi, 4:p, 5:mu, 6:weight, 7:active, 8:p_phi, 9:fff, 10:fnrml
!--------------------------------------------------------------------
      use mpiset
      use grid
      use field, only:babs,babs0
      implicit none

      integer::marker_num,nvec
      real(8)::gc(ngc2,marker_num)
      real(8)::mass
      real(8)::dns(lr,lz,lphi),mom(lr,lz,lphi)
      real(8)::ppara(lr,lz,lphi),pperp(lr,lz,lphi)
      real(8)::dns0(lr,lz,lphi),mom0(lr,lz,lphi)
      real(8)::ppara0(lr,lz,lphi),pperp0(lr,lz,lphi)
      real(8)::wdns(lr,lz,lphi,lpara),wmom(lr,lz,lphi,lpara)
      real(8)::wpar(lr,lz,lphi,lpara),wprp(lr,lz,lphi,lpara)
      integer::ijk_a(3,marker_each)
      real(8)::aaa(8,marker_each)
      integer::i,j,k,l,m,n,ia,ia1,ja,ja1,ka,ka1,lr1,lz1
      real(8)::ar,ar1,az,az1,aphi,aphi1
!      real(8)::www,aaa1,aaa2,aaa3,aaa4,aaa5,aaa6,aaa7,aaa8
      real(8)::vol,cwwa,cwwc
      real(8)::d0,d1,d2,d3,d4,d5,t1,t2,t3,t4,t5
      real(8)::dr1,dz1,dphi1,ma_mi_r
      integer::kr,kz,kphi
      real(8)::p0,p1,p2,mu1
      integer::n_min,n_max !2013-05-22
      integer::nvec0
      integer::nmod

      dr1 = 1.0d0/dr
      dz1 = 1.0d0/dz
      dphi1 = 1.0d0/dphi

      kr   = 1 - kr_offset(my_rank)
      kz   = 1 - kz_offset(my_rank)
      kphi = 1 + lphi_shd - kphi_offset(my_rank)
!      kphi = 3 - kphi_offset(my_rank)
      ma_mi_r = major_r-minor_r

      if(lpara.ge.2)then

!$omp parallel do
      do k = 1, lpara
        do i = 1, lrzphi
          wdns(i,1,1,k) = 0.0d0
          wmom(i,1,1,k) = 0.0d0
          wpar(i,1,1,k) = 0.0d0
          wprp(i,1,1,k) = 0.0d0
        end do
      end do

      end if

!$omp parallel do
      do i = 1, lrzphi
        dns(i,1,1) = 0.0d0
        mom(i,1,1) = 0.0d0
        ppara(i,1,1) = 0.0d0
        pperp(i,1,1) = 0.0d0
      end do


!poption parallel,schedule(dynamic,1000)
!$omp parallel private(ar1,ar,az1,az,aphi1,aphi)
!$omp do schedule(dynamic,1000)
      do n = 1, marker_num
        ijk_a(1,n)=max(1,min(lr  -1,int((gc(1,n)-ma_mi_r)*dr1  ) + kr  ))
        ijk_a(2,n)=max(1,min(lz  -1,int(gc(2,n)          *dz1  ) + kz  ))
        ijk_a(3,n)=max(1,min(lphi-1,int(gc(3,n)        *dphi1) + kphi))

        ar1 = (gc(1,n)-grr(ijk_a(1,n),ijk_a(2,n),ijk_a(3,n) ) ) *dr1
        ar  = 1.0d0 - ar1
        az1 = (gc(2,n)-gzz(ijk_a(1,n),ijk_a(2,n),ijk_a(3,n) ) ) *dz1
        az  = 1.0d0 - az1
        aphi1 = (gc(3,n)-gphi(ijk_a(1,n),ijk_a(2,n),ijk_a(3,n) ) ) *dphi1
        aphi  = 1.0d0 - aphi1

        aaa(1,n) = ar *az *aphi
        aaa(2,n) = ar1*az *aphi
        aaa(3,n) = ar *az1*aphi
        aaa(4,n) = ar1*az1*aphi
        aaa(5,n) = ar *az *aphi1
        aaa(6,n) = ar1*az *aphi1
        aaa(7,n) = ar *az1*aphi1
        aaa(8,n) = ar1*az1*aphi1
      end do
!$omp end parallel

      nvec = marker_num/lpara
      nvec0 = marker_num/lpara
      nmod  = mod(marker_num,lpara)

      if(lpara.ge.2)then

!$omp parallel private(nvec,m,n,ia,ja,ka,p0,p1,p2,mu1)
!$omp do
      do l = 1, lpara
        nvec = nvec0 + min(1,nmod/l)
      do m = 1, nvec
        n = m +  nvec0*(l-1) + min(l-1,nmod)

        ia=ijk_a(1,n)
        ja=ijk_a(2,n)
        ka=ijk_a(3,n)

! gc_x(:,n) = 1:r, 2:z, 3:phi, 4:p, 5:mu, 6:weight, 7:active, 8:p_phi, 9:fff, 10:fnrml
        p0 =  gc(6,n)*gc(7,n)
        p1 =  gc(4,n)*p0
        p2 =  gc(4,n)*p1
        mu1 = gc(5,n)*p0

        wdns(ia,  ja,  ka,  l) = wdns(ia,  ja,  ka,  l) + aaa(1,n)*p0
        wdns(ia+1,ja,  ka,  l) = wdns(ia+1,ja,  ka,  l) + aaa(2,n)*p0
        wdns(ia,  ja+1,ka,  l) = wdns(ia,  ja+1,ka,  l) + aaa(3,n)*p0
        wdns(ia+1,ja+1,ka,  l) = wdns(ia+1,ja+1,ka,  l) + aaa(4,n)*p0
        wdns(ia,  ja,  ka+1,l) = wdns(ia,  ja,  ka+1,l) + aaa(5,n)*p0
        wdns(ia+1,ja,  ka+1,l) = wdns(ia+1,ja,  ka+1,l) + aaa(6,n)*p0
        wdns(ia,  ja+1,ka+1,l) = wdns(ia,  ja+1,ka+1,l) + aaa(7,n)*p0
        wdns(ia+1,ja+1,ka+1,l) = wdns(ia+1,ja+1,ka+1,l) + aaa(8,n)*p0

        wmom(ia,  ja,  ka,  l) = wmom(ia,  ja,  ka,  l) + aaa(1,n)*p1
        wmom(ia+1,ja,  ka,  l) = wmom(ia+1,ja,  ka,  l) + aaa(2,n)*p1
        wmom(ia,  ja+1,ka,  l) = wmom(ia,  ja+1,ka,  l) + aaa(3,n)*p1
        wmom(ia+1,ja+1,ka,  l) = wmom(ia+1,ja+1,ka,  l) + aaa(4,n)*p1
        wmom(ia,  ja,  ka+1,l) = wmom(ia,  ja,  ka+1,l) + aaa(5,n)*p1
        wmom(ia+1,ja,  ka+1,l) = wmom(ia+1,ja,  ka+1,l) + aaa(6,n)*p1
        wmom(ia,  ja+1,ka+1,l) = wmom(ia,  ja+1,ka+1,l) + aaa(7,n)*p1
        wmom(ia+1,ja+1,ka+1,l) = wmom(ia+1,ja+1,ka+1,l) + aaa(8,n)*p1

        wpar(ia,  ja,  ka,  l) = wpar(ia,  ja,  ka,  l) + aaa(1,n)*p2
        wpar(ia+1,ja,  ka,  l) = wpar(ia+1,ja,  ka,  l) + aaa(2,n)*p2
        wpar(ia,  ja+1,ka,  l) = wpar(ia,  ja+1,ka,  l) + aaa(3,n)*p2
        wpar(ia+1,ja+1,ka,  l) = wpar(ia+1,ja+1,ka,  l) + aaa(4,n)*p2
        wpar(ia,  ja,  ka+1,l) = wpar(ia,  ja,  ka+1,l) + aaa(5,n)*p2
        wpar(ia+1,ja,  ka+1,l) = wpar(ia+1,ja,  ka+1,l) + aaa(6,n)*p2
        wpar(ia,  ja+1,ka+1,l) = wpar(ia,  ja+1,ka+1,l) + aaa(7,n)*p2
        wpar(ia+1,ja+1,ka+1,l) = wpar(ia+1,ja+1,ka+1,l) + aaa(8,n)*p2

        wprp(ia,  ja,  ka,  l) = wprp(ia,  ja,  ka,  l) + aaa(1,n)*mu1
        wprp(ia+1,ja,  ka,  l) = wprp(ia+1,ja,  ka,  l) + aaa(2,n)*mu1
        wprp(ia,  ja+1,ka,  l) = wprp(ia,  ja+1,ka,  l) + aaa(3,n)*mu1
        wprp(ia+1,ja+1,ka,  l) = wprp(ia+1,ja+1,ka,  l) + aaa(4,n)*mu1
        wprp(ia,  ja,  ka+1,l) = wprp(ia,  ja,  ka+1,l) + aaa(5,n)*mu1
        wprp(ia+1,ja,  ka+1,l) = wprp(ia+1,ja,  ka+1,l) + aaa(6,n)*mu1
        wprp(ia,  ja+1,ka+1,l) = wprp(ia,  ja+1,ka+1,l) + aaa(7,n)*mu1
        wprp(ia+1,ja+1,ka+1,l) = wprp(ia+1,ja+1,ka+1,l) + aaa(8,n)*mu1
      end do
      end do
!$omp end parallel


      else

      do n = 1, marker_num

        ia=ijk_a(1,n)
        ja=ijk_a(2,n)
        ka=ijk_a(3,n)

! gc_x(:,n) = 1:r, 2:z, 3:phi, 4:p, 5:mu, 6:weight, 7:active, 8:p_phi, 9:fff, 10:fnrml
        p0 =  gc(6,n)*gc(7,n)
        p1 =  gc(4,n)*p0
        p2 =  gc(4,n)*p1
        mu1 = gc(5,n)*p0

        dns(ia,  ja,  ka  ) = dns(ia,  ja,  ka  ) + aaa(1,n)*p0
        dns(ia+1,ja,  ka  ) = dns(ia+1,ja,  ka  ) + aaa(2,n)*p0
        dns(ia,  ja+1,ka  ) = dns(ia,  ja+1,ka  ) + aaa(3,n)*p0
        dns(ia+1,ja+1,ka  ) = dns(ia+1,ja+1,ka  ) + aaa(4,n)*p0
        dns(ia,  ja,  ka+1) = dns(ia,  ja,  ka+1) + aaa(5,n)*p0
        dns(ia+1,ja,  ka+1) = dns(ia+1,ja,  ka+1) + aaa(6,n)*p0
        dns(ia,  ja+1,ka+1) = dns(ia,  ja+1,ka+1) + aaa(7,n)*p0
        dns(ia+1,ja+1,ka+1) = dns(ia+1,ja+1,ka+1) + aaa(8,n)*p0

        mom(ia,  ja,  ka  ) = mom(ia,  ja,  ka  ) + aaa(1,n)*p1
        mom(ia+1,ja,  ka  ) = mom(ia+1,ja,  ka  ) + aaa(2,n)*p1
        mom(ia,  ja+1,ka  ) = mom(ia,  ja+1,ka  ) + aaa(3,n)*p1
        mom(ia+1,ja+1,ka  ) = mom(ia+1,ja+1,ka  ) + aaa(4,n)*p1
        mom(ia,  ja,  ka+1) = mom(ia,  ja,  ka+1) + aaa(5,n)*p1
        mom(ia+1,ja,  ka+1) = mom(ia+1,ja,  ka+1) + aaa(6,n)*p1
        mom(ia,  ja+1,ka+1) = mom(ia,  ja+1,ka+1) + aaa(7,n)*p1
        mom(ia+1,ja+1,ka+1) = mom(ia+1,ja+1,ka+1) + aaa(8,n)*p1

        ppara(ia,  ja,  ka  ) = ppara(ia,  ja,  ka  ) + aaa(1,n)*p2
        ppara(ia+1,ja,  ka  ) = ppara(ia+1,ja,  ka  ) + aaa(2,n)*p2
        ppara(ia,  ja+1,ka  ) = ppara(ia,  ja+1,ka  ) + aaa(3,n)*p2
        ppara(ia+1,ja+1,ka  ) = ppara(ia+1,ja+1,ka  ) + aaa(4,n)*p2
        ppara(ia,  ja,  ka+1) = ppara(ia,  ja,  ka+1) + aaa(5,n)*p2
        ppara(ia+1,ja,  ka+1) = ppara(ia+1,ja,  ka+1) + aaa(6,n)*p2
        ppara(ia,  ja+1,ka+1) = ppara(ia,  ja+1,ka+1) + aaa(7,n)*p2
        ppara(ia+1,ja+1,ka+1) = ppara(ia+1,ja+1,ka+1) + aaa(8,n)*p2

        pperp(ia,  ja,  ka  ) = pperp(ia,  ja,  ka  ) + aaa(1,n)*mu1
        pperp(ia+1,ja,  ka  ) = pperp(ia+1,ja,  ka  ) + aaa(2,n)*mu1
        pperp(ia,  ja+1,ka  ) = pperp(ia,  ja+1,ka  ) + aaa(3,n)*mu1
        pperp(ia+1,ja+1,ka  ) = pperp(ia+1,ja+1,ka  ) + aaa(4,n)*mu1
        pperp(ia,  ja,  ka+1) = pperp(ia,  ja,  ka+1) + aaa(5,n)*mu1
        pperp(ia+1,ja,  ka+1) = pperp(ia+1,ja,  ka+1) + aaa(6,n)*mu1
        pperp(ia,  ja+1,ka+1) = pperp(ia,  ja+1,ka+1) + aaa(7,n)*mu1
        pperp(ia+1,ja+1,ka+1) = pperp(ia+1,ja+1,ka+1) + aaa(8,n)*mu1
      end do

      end if


      if(lpara.ge.2)then

      do k = 1, lpara
!$omp parallel do
      do i = 1, lrzphi
        dns(i,1,1)   = dns(i,1,1)   + wdns(i,1,1,k)
        mom(i,1,1)   = mom(i,1,1)   + wmom(i,1,1,k)
        ppara(i,1,1) = ppara(i,1,1) + wpar(i,1,1,k)
        pperp(i,1,1) = pperp(i,1,1) + wprp(i,1,1,k)
      end do
      end do

      end if

!2016-12-24s

! smoothing
      cwwa = 0.5d0
!      cwwc =-1.d0/6.d0 

!2024-12-21s, correction suggested by Panith Adulsiriswad
      call periodic_particle_mlt4b(dns,mom,ppara,pperp)
      call partsm1(dns,cwwa)
      call partsm1(mom,cwwa)

      call partsm1(ppara,cwwa)
      call partsm1(pperp,cwwa)
!2024-12-21e

! calculate density (per volume)

!      do k = 1, lphi
!      do j = 1, lz
!      do i = 1, lr
!$omp parallel do private(vol)
      do i = 1, lrzphi
        vol = 1.0d0/(grr(i,1,1)*dr*dz*dphi)
        dns(i,1,1) = dns(i,1,1)*vol
        mom(i,1,1) = mom(i,1,1)*vol
        ppara(i,1,1) = ppara(i,1,1)*vol/mass
        pperp(i,1,1) = pperp(i,1,1)*vol*babs(i,1,1)
      end do
!      end do
!      end do

!2016-12-24e


! wall effect

      lr1 = lr - 1
      lz1 = lz - 1

      if(my_rank_r.eq.(mpi_proc_r-1))then
!      do k = 1, lphi
!      do j = 1, lz
!$omp parallel do
       do j = 1, lzphi
          dns(lr1,j,1) = dns(lr1,j,1) + dns(lr,j,1)
          mom(lr1,j,1) = mom(lr1,j,1) + mom(lr,j,1)
          dns(lr,j,1) = 0.0d0
          mom(lr,j,1) = 0.0d0

          ppara(lr1,j,1) = ppara(lr1,j,1) + ppara(lr,j,1)
          pperp(lr1,j,1) = pperp(lr1,j,1) + pperp(lr,j,1)
          ppara(lr,j,1) = 0.0d0
          pperp(lr,j,1) = 0.0d0
       end do
!      end do
      end if

      if(my_rank_r.eq.0)then
!      do k = 1, lphi
!      do j = 1, lz
!$omp parallel do
       do j = 1, lzphi
          dns(2,j,1) = dns(2,j,1) + dns(1,j,1)
          mom(2,j,1) = mom(2,j,1) + mom(1,j,1)
          dns(1,j,1) = 0.0d0
          mom(1,j,1) = 0.0d0

          ppara(2,j,1) = ppara(2,j,1) + ppara(1,j,1)
          pperp(2,j,1) = pperp(2,j,1) + pperp(1,j,1)
          ppara(1,j,1) = 0.0d0
          pperp(1,j,1) = 0.0d0
       end do
!      end do
      end if

      if(my_rank_z.eq.(mpi_proc_z-1))then
!$omp parallel do
       do k = 1, lphi
       do i = 1, lr
          dns(i,lz1,k) = dns(i,lz1,k) + dns(i,lz,k)
          mom(i,lz1,k) = mom(i,lz1,k) + mom(i,lz,k)
          dns(i,lz,k) = 0.0d0
          mom(i,lz,k) = 0.0d0

          ppara(i,lz1,k) = ppara(i,lz1,k) + ppara(i,lz,k)
          pperp(i,lz1,k) = pperp(i,lz1,k) + pperp(i,lz,k)
          ppara(i,lz,k) = 0.0d0
          pperp(i,lz,k) = 0.0d0
       end do
       end do
      end if 

      if(my_rank_z.eq.0)then
!$omp parallel do
       do k = 1, lphi
       do i = 1, lr
          dns(i,2,k) = dns(i,2,k) + dns(i,1,k)
          mom(i,2,k) = mom(i,2,k) + mom(i,1,k)
          dns(i,1,k) = 0.0d0
          mom(i,1,k) = 0.0d0

          ppara(i,2,k) = ppara(i,2,k) + ppara(i,1,k)
          pperp(i,2,k) = pperp(i,2,k) + pperp(i,1,k)
          ppara(i,1,k) = 0.0d0
          pperp(i,1,k) = 0.0d0
       end do
       end do
      end if 

!      if(.not.flag_stored)then
! take only n=1 modes
!      n = 1
!      call n1(n,ppara)
!      call n1(n,pperp)

! take modes n_min <= n <=n_max
!      n_min = 1
!      n_max = 2
!      call lowpass2_mlt2(n_min,n_max,dns,mom)
!      call lowpass2_mlt2(n_min,n_max,ppara,pperp)
!      end if


!      do k = 1, lphi
!      do j = 1, lz
!      do i = 1, lr
!$omp parallel do
      do i = 1, lrzphi
        dns(i,1,1) = dns(i,1,1) + dns0(i,1,1)
        mom(i,1,1) = mom(i,1,1) + mom0(i,1,1)
        ppara(i,1,1) = ppara(i,1,1) + ppara0(i,1,1)
        pperp(i,1,1) = pperp(i,1,1) + pperp0(i,1,1)/babs0(i,1,1)*babs(i,1,1)
      end do
!      end do
!      end do

!            call wall_clock(d5)
!            t5 = d5 - d4
!            if(my_rank.eq.0)then
!            write(6,*)'density 5, t5=',t5
!            end if

end subroutine density
!--------------------------------------------------------------------
subroutine moments(marker_num,mass,gc &
                  ,dns,mom,ppara,pperp,qpara,qperp &
                  ,dns0,mom0,ppara0,pperp0,qpara0,qperp0)
! calculate pressure
! modified on 2015-06-26
! 2016-01-09, modified similar to subroutine density
! gc_x(:,n) = 1:r, 2:z, 3:phi, 4:p, 5:mu, 6:weight, 7:active, 8:p_phi, 9:fff, 10:fnrml
!--------------------------------------------------------------------
      use mpiset
      use grid
      use field, only:babs,babs0
      use gyro, only:marker_each_gyro
      implicit none

      integer::marker_num,nvec
      real(8)::gc(ngc2,marker_num)
      real(8)::mass
      real(8)::dns(lr,lz,lphi),mom(lr,lz,lphi)
      real(8)::ppara(lr,lz,lphi),pperp(lr,lz,lphi)
      real(8)::qpara(lr,lz,lphi),qperp(lr,lz,lphi)
      real(8)::dns0(lr,lz,lphi),mom0(lr,lz,lphi)
      real(8)::ppara0(lr,lz,lphi),pperp0(lr,lz,lphi)
      real(8)::qpara0(lr,lz,lphi),qperp0(lr,lz,lphi)
      real(8)::wdns(lr,lz,lphi,lpara),wmom(lr,lz,lphi,lpara)
      real(8)::wpar(lr,lz,lphi,lpara),wprp(lr,lz,lphi,lpara)
      real(8)::wqar(lr,lz,lphi,lpara),wqrp(lr,lz,lphi,lpara)
!2015-09-21s
      integer::ijk_a(3,marker_each_gyro)
      real(8)::aaa(8,marker_each_gyro)
!2015-09-21e
      integer::i,j,k,l,m,n,ia,ia1,ja,ja1,ka,ka1,lr1,lz1
      real(8)::ar,ar1,az,az1,aphi,aphi1
!      real(8)::www,aaa1,aaa2,aaa3,aaa4,aaa5,aaa6,aaa7,aaa8
      real(8)::vol,cwwa,cwwc
!      real(8)::d0,d1,d2,d3,d4,d5,t1,t2,t3,t4,t5
      real(8)::dr1,dz1,dphi1,ma_mi_r
      integer::kr,kz,kphi
      real(8)::p0,p1,p2,mu1,p3,p1mu,mass1 !2015-09-21
      integer::n_min,n_max !2013-05-22
      integer::nvec0
      integer::nmod

      mass1 = 1.0d0/mass

      dr1 = 1.0d0/dr
      dz1 = 1.0d0/dz
      dphi1 = 1.0d0/dphi

      kr   = 1 - kr_offset(my_rank)
      kz   = 1 - kz_offset(my_rank)
      kphi = 1 + lphi_shd - kphi_offset(my_rank)
!      kphi = 3 - kphi_offset(my_rank)
      ma_mi_r = major_r-minor_r

      if(lpara.ge.2)then
!$omp parallel do
      do k = 1, lpara
        do i = 1, lrzphi
          wdns(i,1,1,k) = 0.0d0
          wmom(i,1,1,k) = 0.0d0
          wpar(i,1,1,k) = 0.0d0
          wprp(i,1,1,k) = 0.0d0
          wqar(i,1,1,k) = 0.0d0
          wqrp(i,1,1,k) = 0.0d0
        end do
      end do

      end if

      do i = 1, lrzphi
        dns(i,1,1) = 0.0d0
        mom(i,1,1) = 0.0d0
        ppara(i,1,1) = 0.0d0
        pperp(i,1,1) = 0.0d0
        qpara(i,1,1) = 0.0d0
        qperp(i,1,1) = 0.0d0
      end do

!poption parallel,schedule(dynamic,1000)
!$omp parallel private(ar1,ar,az1,az,aphi1,aphi)
!$omp do schedule(dynamic,1000)
      do n = 1, marker_num
        ijk_a(1,n)=max(1,min(lr  -1,int((gc(1,n)-ma_mi_r)*dr1  ) + kr  ))
        ijk_a(2,n)=max(1,min(lz  -1,int(gc(2,n)          *dz1  ) + kz  ))
        ijk_a(3,n)=max(1,min(lphi-1,int(gc(3,n)        *dphi1) + kphi))

        ar1 = (gc(1,n)-grr(ijk_a(1,n),ijk_a(2,n),ijk_a(3,n) ) ) *dr1
        ar  = 1.0d0 - ar1
        az1 = (gc(2,n)-gzz(ijk_a(1,n),ijk_a(2,n),ijk_a(3,n) ) ) *dz1
        az  = 1.0d0 - az1
        aphi1 = (gc(3,n)-gphi(ijk_a(1,n),ijk_a(2,n),ijk_a(3,n) ) ) *dphi1
        aphi  = 1.0d0 - aphi1

        aaa(1,n) = ar *az *aphi
        aaa(2,n) = ar1*az *aphi
        aaa(3,n) = ar *az1*aphi
        aaa(4,n) = ar1*az1*aphi
        aaa(5,n) = ar *az *aphi1
        aaa(6,n) = ar1*az *aphi1
        aaa(7,n) = ar *az1*aphi1
        aaa(8,n) = ar1*az1*aphi1
      end do
!$omp end parallel

      nvec = marker_num/lpara
      nvec0 = marker_num/lpara
      nmod  = mod(marker_num,lpara)

      if(lpara.ge.2)then

!poption parallel
!poption indep(wdns,wmom,wpar,wprp,wqar,wqrp)
!$omp parallel private(nvec,m,n,ia,ja,ka,p0,p1,p2,mu1,p3,p1mu)
!$omp do
      do l = 1, lpara
        nvec = nvec0 + min(1,nmod/l)
      do m = 1, nvec
        n = m +  nvec0*(l-1) + min(l-1,nmod)

        ia=ijk_a(1,n)
        ja=ijk_a(2,n)
        ka=ijk_a(3,n)

! gc_x(:,n) = 1:r, 2:z, 3:phi, 4:p, 5:mu, 6:weight, 7:active, 8:p_phi, 9:fff, 10:fnrml
        p0 =  gc(6,n)*gc(7,n)
        p1 =  gc(4,n)*p0
        p2 =  gc(4,n)*p1
        mu1 = gc(5,n)*p0
        p3  = gc(4,n)*p2
        p1mu= gc(4,n)*mu1

        wdns(ia,  ja,  ka,  l) = wdns(ia,  ja,  ka,  l) + aaa(1,n)*p0
        wdns(ia+1,ja,  ka,  l) = wdns(ia+1,ja,  ka,  l) + aaa(2,n)*p0
        wdns(ia,  ja+1,ka,  l) = wdns(ia,  ja+1,ka,  l) + aaa(3,n)*p0
        wdns(ia+1,ja+1,ka,  l) = wdns(ia+1,ja+1,ka,  l) + aaa(4,n)*p0
        wdns(ia,  ja,  ka+1,l) = wdns(ia,  ja,  ka+1,l) + aaa(5,n)*p0
        wdns(ia+1,ja,  ka+1,l) = wdns(ia+1,ja,  ka+1,l) + aaa(6,n)*p0
        wdns(ia,  ja+1,ka+1,l) = wdns(ia,  ja+1,ka+1,l) + aaa(7,n)*p0
        wdns(ia+1,ja+1,ka+1,l) = wdns(ia+1,ja+1,ka+1,l) + aaa(8,n)*p0

        wmom(ia,  ja,  ka,  l) = wmom(ia,  ja,  ka,  l) + aaa(1,n)*p1
        wmom(ia+1,ja,  ka,  l) = wmom(ia+1,ja,  ka,  l) + aaa(2,n)*p1
        wmom(ia,  ja+1,ka,  l) = wmom(ia,  ja+1,ka,  l) + aaa(3,n)*p1
        wmom(ia+1,ja+1,ka,  l) = wmom(ia+1,ja+1,ka,  l) + aaa(4,n)*p1
        wmom(ia,  ja,  ka+1,l) = wmom(ia,  ja,  ka+1,l) + aaa(5,n)*p1
        wmom(ia+1,ja,  ka+1,l) = wmom(ia+1,ja,  ka+1,l) + aaa(6,n)*p1
        wmom(ia,  ja+1,ka+1,l) = wmom(ia,  ja+1,ka+1,l) + aaa(7,n)*p1
        wmom(ia+1,ja+1,ka+1,l) = wmom(ia+1,ja+1,ka+1,l) + aaa(8,n)*p1

        wpar(ia,  ja,  ka,  l) = wpar(ia,  ja,  ka,  l) + aaa(1,n)*p2
        wpar(ia+1,ja,  ka,  l) = wpar(ia+1,ja,  ka,  l) + aaa(2,n)*p2
        wpar(ia,  ja+1,ka,  l) = wpar(ia,  ja+1,ka,  l) + aaa(3,n)*p2
        wpar(ia+1,ja+1,ka,  l) = wpar(ia+1,ja+1,ka,  l) + aaa(4,n)*p2
        wpar(ia,  ja,  ka+1,l) = wpar(ia,  ja,  ka+1,l) + aaa(5,n)*p2
        wpar(ia+1,ja,  ka+1,l) = wpar(ia+1,ja,  ka+1,l) + aaa(6,n)*p2
        wpar(ia,  ja+1,ka+1,l) = wpar(ia,  ja+1,ka+1,l) + aaa(7,n)*p2
        wpar(ia+1,ja+1,ka+1,l) = wpar(ia+1,ja+1,ka+1,l) + aaa(8,n)*p2

        wprp(ia,  ja,  ka,  l) = wprp(ia,  ja,  ka,  l) + aaa(1,n)*mu1
        wprp(ia+1,ja,  ka,  l) = wprp(ia+1,ja,  ka,  l) + aaa(2,n)*mu1
        wprp(ia,  ja+1,ka,  l) = wprp(ia,  ja+1,ka,  l) + aaa(3,n)*mu1
        wprp(ia+1,ja+1,ka,  l) = wprp(ia+1,ja+1,ka,  l) + aaa(4,n)*mu1
        wprp(ia,  ja,  ka+1,l) = wprp(ia,  ja,  ka+1,l) + aaa(5,n)*mu1
        wprp(ia+1,ja,  ka+1,l) = wprp(ia+1,ja,  ka+1,l) + aaa(6,n)*mu1
        wprp(ia,  ja+1,ka+1,l) = wprp(ia,  ja+1,ka+1,l) + aaa(7,n)*mu1
        wprp(ia+1,ja+1,ka+1,l) = wprp(ia+1,ja+1,ka+1,l) + aaa(8,n)*mu1

        wqar(ia,  ja,  ka,  l) = wqar(ia,  ja,  ka,  l) + aaa(1,n)*p3
        wqar(ia+1,ja,  ka,  l) = wqar(ia+1,ja,  ka,  l) + aaa(2,n)*p3
        wqar(ia,  ja+1,ka,  l) = wqar(ia,  ja+1,ka,  l) + aaa(3,n)*p3
        wqar(ia+1,ja+1,ka,  l) = wqar(ia+1,ja+1,ka,  l) + aaa(4,n)*p3
        wqar(ia,  ja,  ka+1,l) = wqar(ia,  ja,  ka+1,l) + aaa(5,n)*p3
        wqar(ia+1,ja,  ka+1,l) = wqar(ia+1,ja,  ka+1,l) + aaa(6,n)*p3
        wqar(ia,  ja+1,ka+1,l) = wqar(ia,  ja+1,ka+1,l) + aaa(7,n)*p3
        wqar(ia+1,ja+1,ka+1,l) = wqar(ia+1,ja+1,ka+1,l) + aaa(8,n)*p3

        wqrp(ia,  ja,  ka,  l) = wqrp(ia,  ja,  ka,  l) + aaa(1,n)*p1mu
        wqrp(ia+1,ja,  ka,  l) = wqrp(ia+1,ja,  ka,  l) + aaa(2,n)*p1mu
        wqrp(ia,  ja+1,ka,  l) = wqrp(ia,  ja+1,ka,  l) + aaa(3,n)*p1mu
        wqrp(ia+1,ja+1,ka,  l) = wqrp(ia+1,ja+1,ka,  l) + aaa(4,n)*p1mu
        wqrp(ia,  ja,  ka+1,l) = wqrp(ia,  ja,  ka+1,l) + aaa(5,n)*p1mu
        wqrp(ia+1,ja,  ka+1,l) = wqrp(ia+1,ja,  ka+1,l) + aaa(6,n)*p1mu
        wqrp(ia,  ja+1,ka+1,l) = wqrp(ia,  ja+1,ka+1,l) + aaa(7,n)*p1mu
        wqrp(ia+1,ja+1,ka+1,l) = wqrp(ia+1,ja+1,ka+1,l) + aaa(8,n)*p1mu

      end do
      end do
!$omp end parallel

      else

      do n = 1, marker_num

        ia=ijk_a(1,n)
        ja=ijk_a(2,n)
        ka=ijk_a(3,n)

! gc_x(:,n) = 1:r, 2:z, 3:phi, 4:p, 5:mu, 6:weight, 7:active, 8:p_phi, 9:fff, 10:fnrml
        p0 =  gc(6,n)*gc(7,n)
        p1 =  gc(4,n)*p0
        p2 =  gc(4,n)*p1
        mu1 = gc(5,n)*p0
        p3  = gc(4,n)*p2
        p1mu= gc(4,n)*mu1

        dns(ia,  ja,  ka  ) = dns(ia,  ja,  ka  ) + aaa(1,n)*p0
        dns(ia+1,ja,  ka  ) = dns(ia+1,ja,  ka  ) + aaa(2,n)*p0
        dns(ia,  ja+1,ka  ) = dns(ia,  ja+1,ka  ) + aaa(3,n)*p0
        dns(ia+1,ja+1,ka  ) = dns(ia+1,ja+1,ka  ) + aaa(4,n)*p0
        dns(ia,  ja,  ka+1) = dns(ia,  ja,  ka+1) + aaa(5,n)*p0
        dns(ia+1,ja,  ka+1) = dns(ia+1,ja,  ka+1) + aaa(6,n)*p0
        dns(ia,  ja+1,ka+1) = dns(ia,  ja+1,ka+1) + aaa(7,n)*p0
        dns(ia+1,ja+1,ka+1) = dns(ia+1,ja+1,ka+1) + aaa(8,n)*p0

        mom(ia,  ja,  ka  ) = mom(ia,  ja,  ka  ) + aaa(1,n)*p1
        mom(ia+1,ja,  ka  ) = mom(ia+1,ja,  ka  ) + aaa(2,n)*p1
        mom(ia,  ja+1,ka  ) = mom(ia,  ja+1,ka  ) + aaa(3,n)*p1
        mom(ia+1,ja+1,ka  ) = mom(ia+1,ja+1,ka  ) + aaa(4,n)*p1
        mom(ia,  ja,  ka+1) = mom(ia,  ja,  ka+1) + aaa(5,n)*p1
        mom(ia+1,ja,  ka+1) = mom(ia+1,ja,  ka+1) + aaa(6,n)*p1
        mom(ia,  ja+1,ka+1) = mom(ia,  ja+1,ka+1) + aaa(7,n)*p1
        mom(ia+1,ja+1,ka+1) = mom(ia+1,ja+1,ka+1) + aaa(8,n)*p1

        ppara(ia,  ja,  ka  ) = ppara(ia,  ja,  ka  ) + aaa(1,n)*p2
        ppara(ia+1,ja,  ka  ) = ppara(ia+1,ja,  ka  ) + aaa(2,n)*p2
        ppara(ia,  ja+1,ka  ) = ppara(ia,  ja+1,ka  ) + aaa(3,n)*p2
        ppara(ia+1,ja+1,ka  ) = ppara(ia+1,ja+1,ka  ) + aaa(4,n)*p2
        ppara(ia,  ja,  ka+1) = ppara(ia,  ja,  ka+1) + aaa(5,n)*p2
        ppara(ia+1,ja,  ka+1) = ppara(ia+1,ja,  ka+1) + aaa(6,n)*p2
        ppara(ia,  ja+1,ka+1) = ppara(ia,  ja+1,ka+1) + aaa(7,n)*p2
        ppara(ia+1,ja+1,ka+1) = ppara(ia+1,ja+1,ka+1) + aaa(8,n)*p2

        pperp(ia,  ja,  ka  ) = pperp(ia,  ja,  ka  ) + aaa(1,n)*mu1
        pperp(ia+1,ja,  ka  ) = pperp(ia+1,ja,  ka  ) + aaa(2,n)*mu1
        pperp(ia,  ja+1,ka  ) = pperp(ia,  ja+1,ka  ) + aaa(3,n)*mu1
        pperp(ia+1,ja+1,ka  ) = pperp(ia+1,ja+1,ka  ) + aaa(4,n)*mu1
        pperp(ia,  ja,  ka+1) = pperp(ia,  ja,  ka+1) + aaa(5,n)*mu1
        pperp(ia+1,ja,  ka+1) = pperp(ia+1,ja,  ka+1) + aaa(6,n)*mu1
        pperp(ia,  ja+1,ka+1) = pperp(ia,  ja+1,ka+1) + aaa(7,n)*mu1
        pperp(ia+1,ja+1,ka+1) = pperp(ia+1,ja+1,ka+1) + aaa(8,n)*mu1

        qpara(ia,  ja,  ka  ) = qpara(ia,  ja,  ka  ) + aaa(1,n)*p3
        qpara(ia+1,ja,  ka  ) = qpara(ia+1,ja,  ka  ) + aaa(2,n)*p3
        qpara(ia,  ja+1,ka  ) = qpara(ia,  ja+1,ka  ) + aaa(3,n)*p3
        qpara(ia+1,ja+1,ka  ) = qpara(ia+1,ja+1,ka  ) + aaa(4,n)*p3
        qpara(ia,  ja,  ka+1) = qpara(ia,  ja,  ka+1) + aaa(5,n)*p3
        qpara(ia+1,ja,  ka+1) = qpara(ia+1,ja,  ka+1) + aaa(6,n)*p3
        qpara(ia,  ja+1,ka+1) = qpara(ia,  ja+1,ka+1) + aaa(7,n)*p3
        qpara(ia+1,ja+1,ka+1) = qpara(ia+1,ja+1,ka+1) + aaa(8,n)*p3

        qperp(ia,  ja,  ka  ) = qperp(ia,  ja,  ka  ) + aaa(1,n)*p1mu
        qperp(ia+1,ja,  ka  ) = qperp(ia+1,ja,  ka  ) + aaa(2,n)*p1mu
        qperp(ia,  ja+1,ka  ) = qperp(ia,  ja+1,ka  ) + aaa(3,n)*p1mu
        qperp(ia+1,ja+1,ka  ) = qperp(ia+1,ja+1,ka  ) + aaa(4,n)*p1mu
        qperp(ia,  ja,  ka+1) = qperp(ia,  ja,  ka+1) + aaa(5,n)*p1mu
        qperp(ia+1,ja,  ka+1) = qperp(ia+1,ja,  ka+1) + aaa(6,n)*p1mu
        qperp(ia,  ja+1,ka+1) = qperp(ia,  ja+1,ka+1) + aaa(7,n)*p1mu
        qperp(ia+1,ja+1,ka+1) = qperp(ia+1,ja+1,ka+1) + aaa(8,n)*p1mu

      end do

      end if


      if(lpara.ge.2)then

      do k = 1, lpara
!$omp parallel do
      do i = 1, lrzphi
        dns(i,1,1)   = dns(i,1,1)   + wdns(i,1,1,k)
        mom(i,1,1)   = mom(i,1,1)   + wmom(i,1,1,k)
        ppara(i,1,1) = ppara(i,1,1) + wpar(i,1,1,k)
        pperp(i,1,1) = pperp(i,1,1) + wprp(i,1,1,k)
        qpara(i,1,1) = qpara(i,1,1) + wqar(i,1,1,k)
        qperp(i,1,1) = qperp(i,1,1) + wqrp(i,1,1,k)
      end do
      end do

      end if

! smoothing

       cwwa = 0.5d0
!       cwwc =-1.d0/6.d0 

!2024-12-21s, correction suggested by Panith Adulsiriswad
       call periodic_particle_mlt6b(dns,mom,ppara,pperp,qpara,qperp)
       call partsm1(dns,cwwa)
       call partsm1(mom,cwwa)

       call partsm1(ppara,cwwa)
       call partsm1(pperp,cwwa)

       call partsm1(qpara,cwwa)
       call partsm1(qperp,cwwa)
!2024-12-21e

! calculate density (per volume)

!      do k = 1, lphi
!      do j = 1, lz
!      do i = 1, lr
!$omp parallel do private(vol)
      do i = 1, lrzphi
        vol = 1.0d0/(grr(i,1,1)*dr*dz*dphi)
        dns(i,1,1) = dns(i,1,1)*vol
        mom(i,1,1) = mom(i,1,1)*vol
        ppara(i,1,1) = ppara(i,1,1)*vol*mass1
        pperp(i,1,1) = pperp(i,1,1)*vol*babs(i,1,1)
        qpara(i,1,1) = qpara(i,1,1)*vol*mass1**2
        qperp(i,1,1) = qperp(i,1,1)*vol*babs(i,1,1)*mass1
      end do
!      end do
!      end do

! wall effect

      lr1 = lr - 1
      lz1 = lz - 1

      if(my_rank_r.eq.(mpi_proc_r-1))then
!      do k = 1, lphi
!      do j = 1, lz
!$omp parallel do
       do j = 1, lzphi
          dns(lr1,j,1) = dns(lr1,j,1) + dns(lr,j,1)
          mom(lr1,j,1) = mom(lr1,j,1) + mom(lr,j,1)
          dns(lr,j,1) = 0.0d0
          mom(lr,j,1) = 0.0d0

          ppara(lr1,j,1) = ppara(lr1,j,1) + ppara(lr,j,1)
          pperp(lr1,j,1) = pperp(lr1,j,1) + pperp(lr,j,1)
          ppara(lr,j,1) = 0.0d0
          pperp(lr,j,1) = 0.0d0

          qpara(lr1,j,1) = qpara(lr1,j,1) + qpara(lr,j,1)
          qperp(lr1,j,1) = qperp(lr1,j,1) + qperp(lr,j,1)
          qpara(lr,j,1) = 0.0d0
          qperp(lr,j,1) = 0.0d0
       end do
!      end do
      end if

      if(my_rank_r.eq.0)then
!      do k = 1, lphi
!      do j = 1, lz
!$omp parallel do
       do j = 1, lzphi
          dns(2,j,1) = dns(2,j,1) + dns(1,j,1)
          mom(2,j,1) = mom(2,j,1) + mom(1,j,1)
          dns(1,j,1) = 0.0d0
          mom(1,j,1) = 0.0d0

          ppara(2,j,1) = ppara(2,j,1) + ppara(1,j,1)
          pperp(2,j,1) = pperp(2,j,1) + pperp(1,j,1)
          ppara(1,j,1) = 0.0d0
          pperp(1,j,1) = 0.0d0

          qpara(2,j,1) = qpara(2,j,1) + qpara(1,j,1)
          qperp(2,j,1) = qperp(2,j,1) + qperp(1,j,1)
          qpara(1,j,1) = 0.0d0
          qperp(1,j,1) = 0.0d0
       end do
!      end do
      end if

      if(my_rank_z.eq.(mpi_proc_z-1))then
!$omp parallel do
       do k = 1, lphi
       do i = 1, lr
          dns(i,lz1,k) = dns(i,lz1,k) + dns(i,lz,k)
          mom(i,lz1,k) = mom(i,lz1,k) + mom(i,lz,k)
          dns(i,lz,k) = 0.0d0
          mom(i,lz,k) = 0.0d0

          ppara(i,lz1,k) = ppara(i,lz1,k) + ppara(i,lz,k)
          pperp(i,lz1,k) = pperp(i,lz1,k) + pperp(i,lz,k)
          ppara(i,lz,k) = 0.0d0
          pperp(i,lz,k) = 0.0d0

          qpara(i,lz1,k) = qpara(i,lz1,k) + qpara(i,lz,k)
          qperp(i,lz1,k) = qperp(i,lz1,k) + qperp(i,lz,k)
          qpara(i,lz,k) = 0.0d0
          qperp(i,lz,k) = 0.0d0
       end do
       end do
      end if 

      if(my_rank_z.eq.0)then
!$omp parallel do
       do k = 1, lphi
       do i = 1, lr
          dns(i,2,k) = dns(i,2,k) + dns(i,1,k)
          mom(i,2,k) = mom(i,2,k) + mom(i,1,k)
          dns(i,1,k) = 0.0d0
          mom(i,1,k) = 0.0d0

          ppara(i,2,k) = ppara(i,2,k) + ppara(i,1,k)
          pperp(i,2,k) = pperp(i,2,k) + pperp(i,1,k)
          ppara(i,1,k) = 0.0d0
          pperp(i,1,k) = 0.0d0

          qpara(i,2,k) = qpara(i,2,k) + qpara(i,1,k)
          qperp(i,2,k) = qperp(i,2,k) + qperp(i,1,k)
          qpara(i,1,k) = 0.0d0
          qperp(i,1,k) = 0.0d0
       end do
       end do
      end if 

!      if(.not.flag_stored)then
! take modes n_min <= n <=n_max
!      n_min = 1
!      n_max = 2
!      call lowpass2_mlt2(n_min,n_max,dns,mom)
!      call lowpass2_mlt2(n_min,n_max,ppara,pperp)
!      call lowpass2_mlt2(n_min,n_max,qpara,qperp)
!      end if

!2014-08-28s
!$omp parallel do
      do i = 1, lrzphi
        dns(i,1,1) = dns(i,1,1) + dns0(i,1,1)
        mom(i,1,1) = mom(i,1,1) + mom0(i,1,1)
        ppara(i,1,1) = ppara(i,1,1) + ppara0(i,1,1)
        pperp(i,1,1) = pperp(i,1,1) + pperp0(i,1,1)/babs0(i,1,1)*babs(i,1,1)
        qpara(i,1,1) = qpara(i,1,1) + qpara0(i,1,1)
        qperp(i,1,1) = qperp(i,1,1) + qperp0(i,1,1)/babs0(i,1,1)*babs(i,1,1)
      end do
!2014-08-28e

!            call wall_clock(d5)
!            t5 = d5 - d4
!            if(my_rank.eq.0)then
!            write(6,*)'density 5, t5=',t5
!            end if

end subroutine moments
!--------------------------------------------------------------------
subroutine emf_gyro(marker_num,marker_num_gyro,gc,gyro_phys &
                   ,flp,flp_gyro)
! modified for GK simulation on Fugaku 2022-01-11
! gc_x(:,n) = 1:r, 2:z, 3:phi, 4:p, 5:mu, 6:weight, 7:active, 8:p_phi, 9:fff, 10:fnrml
! flp(nflp, marker_num); nflp=30, nflp_gyro=7
!--------------------------------------------------------------------
      use mpiset
      use grid
      use field, only:babs,babs0,er,ez,ephi,epara,br,bz,bphi,br0,bz0,bphi0,fld
      use gyro, only:ngyro,marker_each_gyro
      implicit none

      integer::marker_num,nvec,marker_num_gyro
      real(8)::gc(ngc2,marker_num)
      real(8)::gyro_phys(2,marker_num_gyro)
      real(8)::flp(nflp,marker_num) !2015-07-08
      real(8)::flp_gyro(nflp_gyro,marker_num_gyro)
      integer::i,j,k,l,m,n,ia,ja,ka,lr1,lz1
      real(8)::ar,ar1,az,az1,aphi,aphi1
!      real(8)::aaa1,aaa2,aaa3,aaa4,aaa5,aaa6,aaa7,aaa8
      real(8)::dr1,dz1,dphi1,ma_mi_r
      integer::kr,kz,kphi
      integer::nc,in
      real(8)::rngr1,rngr1_ac
      integer::nvec0
      integer::nmod

!2015-07-08s
      integer::ijk_a(3,marker_each_gyro)
      real(8)::aaa(8,marker_each_gyro)
!2015-07-08e

      dr1 = 1.0d0/dr
      dz1 = 1.0d0/dz
      dphi1 = 1.0d0/dphi
      rngr1 = 1.0d0/dble(ngyro)

      kr   = 1 - kr_offset(my_rank)
      kz   = 1 - kz_offset(my_rank)
      kphi = 1 + lphi_shd - kphi_offset(my_rank)
!      kphi = 3 - kphi_offset(my_rank)
      ma_mi_r = major_r-minor_r


!$omp parallel
!$omp workshare

! for subroutine extract_em

!      flp_gyro = 0.0d0

      fld(1,:,:,:) = er(:,:,:)
      fld(2,:,:,:) = ez(:,:,:)
      fld(3,:,:,:) = ephi(:,:,:)
      fld(4,:,:,:) = epara(:,:,:)
      fld(5,:,:,:) = br(:,:,:) - br0(:,:,:)
      fld(6,:,:,:) = bz(:,:,:) - bz0(:,:,:)
      fld(7,:,:,:) = bphi(:,:,:) - bphi0(:,:,:)
!$omp end workshare
!$omp end parallel


!$omp parallel private(nc,ar1,ar,az1,az,aphi1,aphi)

!poption parallel,schedule(dynamic,1000)
!$omp do schedule(dynamic,1000)
      do n = 1, marker_num_gyro
        nc =(n-1)/ngyro + 1
        ijk_a(1,n)=max(1,min(lr  -1,int((gyro_phys(1,n)-ma_mi_r)*dr1  ) + kr  ))
        ijk_a(2,n)=max(1,min(lz  -1,int(gyro_phys(2,n)          *dz1  ) + kz  ))
        ijk_a(3,n)=max(1,min(lphi-1,int(gc(3,nc)          *dphi1) + kphi))

        ar1  = max(0.0d0, min(1.0d0, (gyro_phys(1,n) - ma_mi_r)*dr1 -  dble(ijk_a(1,n) - kr)  ) )
        ar   = 1.0d0 - ar1
        az1  = max(0.0d0, min(1.0d0, gyro_phys(2,n)*dz1 - dble(ijk_a(2,n) - kz) ) )
        az   = 1.0d0 - az1
        aphi1 = max(0.0d0, min(1.0d0, gc(3,nc)*dphi1 - dble(ijk_a(3,n) - kphi) ) )
        aphi = 1.0d0 - aphi1

        aaa(1,n) = ar *az *aphi
        aaa(2,n) = ar1*az *aphi
        aaa(3,n) = ar *az1*aphi
        aaa(4,n) = ar1*az1*aphi
        aaa(5,n) = ar *az *aphi1
        aaa(6,n) = ar1*az *aphi1
        aaa(7,n) = ar *az1*aphi1
        aaa(8,n) = ar1*az1*aphi1
      end do
!$omp end parallel


!poption parallel,schedule(dynamic,1000)
!$omp parallel private(ia,ja,ka)
!$omp do schedule(dynamic,1000)
      do n = 1, marker_num_gyro
        ia=ijk_a(1,n)
        ja=ijk_a(2,n)
        ka=ijk_a(3,n)

        do in = 1, nflp_gyro
          flp_gyro(in,n) = fld(in, ia, ja,  ka  )*aaa(1,n) + fld(in, ia+1,ja,  ka  )*aaa(2,n) &
                         + fld(in, ia, ja+1,ka  )*aaa(3,n) + fld(in, ia+1,ja+1,ka  )*aaa(4,n) &
                         + fld(in, ia, ja,  ka+1)*aaa(5,n) + fld(in, ia+1,ja,  ka+1)*aaa(6,n) &
                         + fld(in, ia, ja+1,ka+1)*aaa(7,n) + fld(in, ia+1,ja+1,ka+1)*aaa(8,n)
        end do

      end do
!$omp end parallel


!$omp parallel do private(rngr1_ac,in)
      do nc = 1, marker_num
        rngr1_ac = rngr1*gc(7,nc)
        do in = 1, nflp_gyro
          flp(in,nc) =(flp_gyro(in,ngyro*(nc-1)+1) &
                      +flp_gyro(in,ngyro*(nc-1)+2) &
                      +flp_gyro(in,ngyro*(nc-1)+3) &
                      +flp_gyro(in,ngyro*(nc-1)+4) &
                      )*rngr1_ac
        end do
      end do

end subroutine emf_gyro
!--------------------------------------------------------------------
subroutine density_gyro(marker_num,marker_num_gyro,mass,gc,gyro_phys &
                       ,dns,mom,ppara,pperp,dns0,mom0,ppara0,pperp0)
! 2022-01-11, for FLR case
! gc_x(:,n) = 1:r, 2:z, 3:phi, 4:p, 5:mu, 6:weight, 7:active, 8:p_phi, 9:fff, 10:fnrml
!--------------------------------------------------------------------
      use mpiset
      use grid
      use field, only:babs,babs0
      use gyro, only:ngyro,marker_each_gyro
      implicit none

      integer::marker_num,nvec,marker_num_gyro
      real(8)::gc(ngc2,marker_num)
      real(8)::gyro_phys(2,marker_num_gyro)
      real(8)::mass
      real(8)::dns(lr,lz,lphi),mom(lr,lz,lphi)
      real(8)::ppara(lr,lz,lphi),pperp(lr,lz,lphi)
      real(8)::dns0(lr,lz,lphi),mom0(lr,lz,lphi)
      real(8)::ppara0(lr,lz,lphi),pperp0(lr,lz,lphi)
      real(8)::wdns(lr,lz,lphi,lpara),wmom(lr,lz,lphi,lpara)
      real(8)::wpar(lr,lz,lphi,lpara),wprp(lr,lz,lphi,lpara)
!2015-09-21s
      integer::ijk_a(3,marker_each_gyro)
      real(8)::aaa(8,marker_each_gyro)
!2015-09-21e
      integer::i,j,k,l,m,n,ia,ia1,ja,ja1,ka,ka1,lr1,lz1
      real(8)::ar,ar1,az,az1,aphi,aphi1
!      real(8)::www,aaa1,aaa2,aaa3,aaa4,aaa5,aaa6,aaa7,aaa8
      real(8)::vol,cwwa,cwwc
!      real(8)::d0,d1,d2,d3,d4,d5,t1,t2,t3,t4,t5
      real(8)::dr1,dz1,dphi1,ma_mi_r
      integer::kr,kz,kphi
      real(8)::p0,p1,p2,mu1,mass1 !2015-09-21
      integer::n_min,n_max !2013-05-22
      real(8)::rngr1
      integer::nc
      integer::nvec0
      integer::nmod

      mass1 = 1.0d0/mass

      dr1 = 1.0d0/dr
      dz1 = 1.0d0/dz
      dphi1 = 1.0d0/dphi
      rngr1 = 1.0d0/dble(ngyro)

      kr   = 1 - kr_offset(my_rank)
      kz   = 1 - kz_offset(my_rank)
      kphi = 1 + lphi_shd - kphi_offset(my_rank)
!      kphi = 3 - kphi_offset(my_rank)
      ma_mi_r = major_r-minor_r

      if(lpara.ge.2)then
!$omp parallel do
      do k = 1, lpara
        do i = 1, lrzphi
          wdns(i,1,1,k) = 0.0d0
          wmom(i,1,1,k) = 0.0d0
          wpar(i,1,1,k) = 0.0d0
          wprp(i,1,1,k) = 0.0d0
        end do
      end do

      end if

!$omp parallel do
      do i = 1, lrzphi
        dns(i,1,1) = 0.0d0
        mom(i,1,1) = 0.0d0
        ppara(i,1,1) = 0.0d0
        pperp(i,1,1) = 0.0d0
      end do

!$omp parallel private(nc,ar1,ar,az1,az,aphi1,aphi)
!poption parallel,schedule(dynamic,1000)
!$omp do schedule(dynamic,1000)
      do n = 1, marker_num_gyro
        nc =(n-1)/ngyro + 1
        ijk_a(1,n)=max(1,min(lr  -1,int((gyro_phys(1,n)-ma_mi_r)*dr1  ) + kr  ))
        ijk_a(2,n)=max(1,min(lz  -1,int(gyro_phys(2,n)          *dz1  ) + kz  ))
        ijk_a(3,n)=max(1,min(lphi-1,int(gc(3,nc)          *dphi1) + kphi))

        ar1  = max(0.0d0, min(1.0d0, (gyro_phys(1,n) - ma_mi_r)*dr1 -  dble(ijk_a(1,n) - kr)  ) )
        ar   = 1.0d0 - ar1
        az1  = max(0.0d0, min(1.0d0, gyro_phys(2,n)*dz1 - dble(ijk_a(2,n) - kz) ) )
        az   = 1.0d0 - az1
        aphi1 = max(0.0d0, min(1.0d0, gc(3,nc)*dphi1 - dble(ijk_a(3,n) - kphi) ) )
        aphi = 1.0d0 - aphi1

        aaa(1,n) = ar *az *aphi
        aaa(2,n) = ar1*az *aphi
        aaa(3,n) = ar *az1*aphi
        aaa(4,n) = ar1*az1*aphi
        aaa(5,n) = ar *az *aphi1
        aaa(6,n) = ar1*az *aphi1
        aaa(7,n) = ar *az1*aphi1
        aaa(8,n) = ar1*az1*aphi1
      end do
!$omp end parallel

      nvec = marker_num_gyro/lpara
      nvec0 = marker_num_gyro/lpara
      nmod  = mod(marker_num_gyro,lpara)

      if(lpara.ge.2)then

!2017-02-10
!$omp parallel private(nvec,m,n,nc,ia,ja,ka,p0,p1,p2,mu1)
!$omp do
      do l = 1, lpara
        nvec = nvec0 + min(1,nmod/l)
      do m = 1, nvec
        n = m +  nvec0*(l-1) + min(l-1,nmod)
        nc =(n-1)/ngyro + 1

        ia=ijk_a(1,n)
        ja=ijk_a(2,n)
        ka=ijk_a(3,n)

! gc_x(:,n) = 1:r, 2:z, 3:phi, 4:p, 5:mu, 6:weight, 7:active, 8:p_phi, 9:fff, 10:fnrml
        p0 =  gc(6,nc)*gc(7,nc)*rngr1 !rngr1 is multiplied
        p1 =  gc(4,nc)*p0
        p2 =  gc(4,nc)*p1
        mu1 = gc(5,nc)*p0

        wdns(ia,  ja,  ka,  l) = wdns(ia,  ja,  ka,  l) + aaa(1,n)*p0
        wdns(ia+1,ja,  ka,  l) = wdns(ia+1,ja,  ka,  l) + aaa(2,n)*p0
        wdns(ia,  ja+1,ka,  l) = wdns(ia,  ja+1,ka,  l) + aaa(3,n)*p0
        wdns(ia+1,ja+1,ka,  l) = wdns(ia+1,ja+1,ka,  l) + aaa(4,n)*p0
        wdns(ia,  ja,  ka+1,l) = wdns(ia,  ja,  ka+1,l) + aaa(5,n)*p0
        wdns(ia+1,ja,  ka+1,l) = wdns(ia+1,ja,  ka+1,l) + aaa(6,n)*p0
        wdns(ia,  ja+1,ka+1,l) = wdns(ia,  ja+1,ka+1,l) + aaa(7,n)*p0
        wdns(ia+1,ja+1,ka+1,l) = wdns(ia+1,ja+1,ka+1,l) + aaa(8,n)*p0

        wmom(ia,  ja,  ka,  l) = wmom(ia,  ja,  ka,  l) + aaa(1,n)*p1
        wmom(ia+1,ja,  ka,  l) = wmom(ia+1,ja,  ka,  l) + aaa(2,n)*p1
        wmom(ia,  ja+1,ka,  l) = wmom(ia,  ja+1,ka,  l) + aaa(3,n)*p1
        wmom(ia+1,ja+1,ka,  l) = wmom(ia+1,ja+1,ka,  l) + aaa(4,n)*p1
        wmom(ia,  ja,  ka+1,l) = wmom(ia,  ja,  ka+1,l) + aaa(5,n)*p1
        wmom(ia+1,ja,  ka+1,l) = wmom(ia+1,ja,  ka+1,l) + aaa(6,n)*p1
        wmom(ia,  ja+1,ka+1,l) = wmom(ia,  ja+1,ka+1,l) + aaa(7,n)*p1
        wmom(ia+1,ja+1,ka+1,l) = wmom(ia+1,ja+1,ka+1,l) + aaa(8,n)*p1

        wpar(ia,  ja,  ka,  l) = wpar(ia,  ja,  ka,  l) + aaa(1,n)*p2
        wpar(ia+1,ja,  ka,  l) = wpar(ia+1,ja,  ka,  l) + aaa(2,n)*p2
        wpar(ia,  ja+1,ka,  l) = wpar(ia,  ja+1,ka,  l) + aaa(3,n)*p2
        wpar(ia+1,ja+1,ka,  l) = wpar(ia+1,ja+1,ka,  l) + aaa(4,n)*p2
        wpar(ia,  ja,  ka+1,l) = wpar(ia,  ja,  ka+1,l) + aaa(5,n)*p2
        wpar(ia+1,ja,  ka+1,l) = wpar(ia+1,ja,  ka+1,l) + aaa(6,n)*p2
        wpar(ia,  ja+1,ka+1,l) = wpar(ia,  ja+1,ka+1,l) + aaa(7,n)*p2
        wpar(ia+1,ja+1,ka+1,l) = wpar(ia+1,ja+1,ka+1,l) + aaa(8,n)*p2

        wprp(ia,  ja,  ka,  l) = wprp(ia,  ja,  ka,  l) + aaa(1,n)*mu1
        wprp(ia+1,ja,  ka,  l) = wprp(ia+1,ja,  ka,  l) + aaa(2,n)*mu1
        wprp(ia,  ja+1,ka,  l) = wprp(ia,  ja+1,ka,  l) + aaa(3,n)*mu1
        wprp(ia+1,ja+1,ka,  l) = wprp(ia+1,ja+1,ka,  l) + aaa(4,n)*mu1
        wprp(ia,  ja,  ka+1,l) = wprp(ia,  ja,  ka+1,l) + aaa(5,n)*mu1
        wprp(ia+1,ja,  ka+1,l) = wprp(ia+1,ja,  ka+1,l) + aaa(6,n)*mu1
        wprp(ia,  ja+1,ka+1,l) = wprp(ia,  ja+1,ka+1,l) + aaa(7,n)*mu1
        wprp(ia+1,ja+1,ka+1,l) = wprp(ia+1,ja+1,ka+1,l) + aaa(8,n)*mu1
      end do
      end do
!$omp end parallel

      else

      do n = 1, marker_num_gyro
        nc =(n-1)/ngyro + 1

        ia=ijk_a(1,n)
        ja=ijk_a(2,n)
        ka=ijk_a(3,n)

! gc_x(:,n) = 1:r, 2:z, 3:phi, 4:p, 5:mu, 6:weight, 7:active, 8:p_phi, 9:fff, 10:fnrml
        p0 =  gc(6,nc)*gc(7,nc)*rngr1 !rngr1 is multiplied
        p1 =  gc(4,nc)*p0
        p2 =  gc(4,nc)*p1
        mu1 = gc(5,nc)*p0

        dns(ia,  ja,  ka  ) = dns(ia,  ja,  ka  ) + aaa(1,n)*p0
        dns(ia+1,ja,  ka  ) = dns(ia+1,ja,  ka  ) + aaa(2,n)*p0
        dns(ia,  ja+1,ka  ) = dns(ia,  ja+1,ka  ) + aaa(3,n)*p0
        dns(ia+1,ja+1,ka  ) = dns(ia+1,ja+1,ka  ) + aaa(4,n)*p0
        dns(ia,  ja,  ka+1) = dns(ia,  ja,  ka+1) + aaa(5,n)*p0
        dns(ia+1,ja,  ka+1) = dns(ia+1,ja,  ka+1) + aaa(6,n)*p0
        dns(ia,  ja+1,ka+1) = dns(ia,  ja+1,ka+1) + aaa(7,n)*p0
        dns(ia+1,ja+1,ka+1) = dns(ia+1,ja+1,ka+1) + aaa(8,n)*p0

        mom(ia,  ja,  ka  ) = mom(ia,  ja,  ka  ) + aaa(1,n)*p1
        mom(ia+1,ja,  ka  ) = mom(ia+1,ja,  ka  ) + aaa(2,n)*p1
        mom(ia,  ja+1,ka  ) = mom(ia,  ja+1,ka  ) + aaa(3,n)*p1
        mom(ia+1,ja+1,ka  ) = mom(ia+1,ja+1,ka  ) + aaa(4,n)*p1
        mom(ia,  ja,  ka+1) = mom(ia,  ja,  ka+1) + aaa(5,n)*p1
        mom(ia+1,ja,  ka+1) = mom(ia+1,ja,  ka+1) + aaa(6,n)*p1
        mom(ia,  ja+1,ka+1) = mom(ia,  ja+1,ka+1) + aaa(7,n)*p1
        mom(ia+1,ja+1,ka+1) = mom(ia+1,ja+1,ka+1) + aaa(8,n)*p1

        ppara(ia,  ja,  ka  ) = ppara(ia,  ja,  ka  ) + aaa(1,n)*p2
        ppara(ia+1,ja,  ka  ) = ppara(ia+1,ja,  ka  ) + aaa(2,n)*p2
        ppara(ia,  ja+1,ka  ) = ppara(ia,  ja+1,ka  ) + aaa(3,n)*p2
        ppara(ia+1,ja+1,ka  ) = ppara(ia+1,ja+1,ka  ) + aaa(4,n)*p2
        ppara(ia,  ja,  ka+1) = ppara(ia,  ja,  ka+1) + aaa(5,n)*p2
        ppara(ia+1,ja,  ka+1) = ppara(ia+1,ja,  ka+1) + aaa(6,n)*p2
        ppara(ia,  ja+1,ka+1) = ppara(ia,  ja+1,ka+1) + aaa(7,n)*p2
        ppara(ia+1,ja+1,ka+1) = ppara(ia+1,ja+1,ka+1) + aaa(8,n)*p2

        pperp(ia,  ja,  ka  ) = pperp(ia,  ja,  ka  ) + aaa(1,n)*mu1
        pperp(ia+1,ja,  ka  ) = pperp(ia+1,ja,  ka  ) + aaa(2,n)*mu1
        pperp(ia,  ja+1,ka  ) = pperp(ia,  ja+1,ka  ) + aaa(3,n)*mu1
        pperp(ia+1,ja+1,ka  ) = pperp(ia+1,ja+1,ka  ) + aaa(4,n)*mu1
        pperp(ia,  ja,  ka+1) = pperp(ia,  ja,  ka+1) + aaa(5,n)*mu1
        pperp(ia+1,ja,  ka+1) = pperp(ia+1,ja,  ka+1) + aaa(6,n)*mu1
        pperp(ia,  ja+1,ka+1) = pperp(ia,  ja+1,ka+1) + aaa(7,n)*mu1
        pperp(ia+1,ja+1,ka+1) = pperp(ia+1,ja+1,ka+1) + aaa(8,n)*mu1
      end do

      end if


      if(lpara.ge.2)then

      do k = 1, lpara
!$omp parallel do
      do i = 1, lrzphi
        dns(i,1,1)   = dns(i,1,1)   + wdns(i,1,1,k)
        mom(i,1,1)   = mom(i,1,1)   + wmom(i,1,1,k)
        ppara(i,1,1) = ppara(i,1,1) + wpar(i,1,1,k)
        pperp(i,1,1) = pperp(i,1,1) + wprp(i,1,1,k)
      end do
      end do

      end if


! smoothing

       cwwa = 0.5d0
!       cwwc =-1.d0/6.d0 

!2024-12-21s, correction suggested by Panith Adulsiriswad
       call periodic_particle_mlt4b(dns,mom,ppara,pperp)
       call partsm1(dns,cwwa)
       call partsm1(mom,cwwa)

       call partsm1(ppara,cwwa)
       call partsm1(pperp,cwwa)
!2024-12-21e
      
! calculate density (per volume)

!      do k = 1, lphi
!      do j = 1, lz
!      do i = 1, lr
!$omp parallel do private(vol)
      do i = 1, lrzphi
        vol = 1.0d0/(grr(i,1,1)*dr*dz*dphi)
        dns(i,1,1) = dns(i,1,1)*vol
        mom(i,1,1) = mom(i,1,1)*vol
        ppara(i,1,1) = ppara(i,1,1)*vol/mass
        pperp(i,1,1) = pperp(i,1,1)*vol*babs(i,1,1)
      end do
!      end do
!      end do


! wall effect

      lr1 = lr - 1
      lz1 = lz - 1

      if(my_rank_r.eq.(mpi_proc_r-1))then
!      do k = 1, lphi
!      do j = 1, lz
!$omp parallel do
       do j = 1, lzphi
          dns(lr1,j,1) = dns(lr1,j,1) + dns(lr,j,1)
          mom(lr1,j,1) = mom(lr1,j,1) + mom(lr,j,1)
          dns(lr,j,1) = 0.0d0
          mom(lr,j,1) = 0.0d0

          ppara(lr1,j,1) = ppara(lr1,j,1) + ppara(lr,j,1)
          pperp(lr1,j,1) = pperp(lr1,j,1) + pperp(lr,j,1)
          ppara(lr,j,1) = 0.0d0
          pperp(lr,j,1) = 0.0d0
       end do
!      end do
      end if

      if(my_rank_r.eq.0)then
!      do k = 1, lphi
!      do j = 1, lz
!$omp parallel do
       do j = 1, lzphi
          dns(2,j,1) = dns(2,j,1) + dns(1,j,1)
          mom(2,j,1) = mom(2,j,1) + mom(1,j,1)
          dns(1,j,1) = 0.0d0
          mom(1,j,1) = 0.0d0

          ppara(2,j,1) = ppara(2,j,1) + ppara(1,j,1)
          pperp(2,j,1) = pperp(2,j,1) + pperp(1,j,1)
          ppara(1,j,1) = 0.0d0
          pperp(1,j,1) = 0.0d0
       end do
!      end do
      end if

      if(my_rank_z.eq.(mpi_proc_z-1))then
!$omp parallel do
       do k = 1, lphi
       do i = 1, lr
          dns(i,lz1,k) = dns(i,lz1,k) + dns(i,lz,k)
          mom(i,lz1,k) = mom(i,lz1,k) + mom(i,lz,k)
          dns(i,lz,k) = 0.0d0
          mom(i,lz,k) = 0.0d0

          ppara(i,lz1,k) = ppara(i,lz1,k) + ppara(i,lz,k)
          pperp(i,lz1,k) = pperp(i,lz1,k) + pperp(i,lz,k)
          ppara(i,lz,k) = 0.0d0
          pperp(i,lz,k) = 0.0d0
       end do
       end do
      end if 

      if(my_rank_z.eq.0)then
!$omp parallel do
       do k = 1, lphi
       do i = 1, lr
          dns(i,2,k) = dns(i,2,k) + dns(i,1,k)
          mom(i,2,k) = mom(i,2,k) + mom(i,1,k)
          dns(i,1,k) = 0.0d0
          mom(i,1,k) = 0.0d0

          ppara(i,2,k) = ppara(i,2,k) + ppara(i,1,k)
          pperp(i,2,k) = pperp(i,2,k) + pperp(i,1,k)
          ppara(i,1,k) = 0.0d0
          pperp(i,1,k) = 0.0d0
       end do
       end do
      end if 

! take only n=1 modes
!      n = 1
!      call n1(n,ppara)
!      call n1(n,pperp)

! take modes n_min <= n <=n_max
!      n_min = 1
!      n_max = 2
!      call lowpass2_mlt2(n_min,n_max,dns,mom)
!      call lowpass2_mlt2(n_min,n_max,ppara,pperp)

!      do k = 1, lphi
!      do j = 1, lz
!      do i = 1, lr
!$omp parallel do
      do i = 1, lrzphi
        dns(i,1,1) = dns(i,1,1) + dns0(i,1,1)
        mom(i,1,1) = mom(i,1,1) + mom0(i,1,1)
        ppara(i,1,1) = ppara(i,1,1) + ppara0(i,1,1)
        pperp(i,1,1) = pperp(i,1,1) + pperp0(i,1,1)/babs0(i,1,1)*babs(i,1,1)
      end do
!      end do
!      end do

end subroutine density_gyro
!--------------------------------------------------------------------
subroutine moments_gyro(marker_num,marker_num_gyro,mass,gc,gyro_phys &
                  ,dns,mom,ppara,pperp,qpara,qperp &
                  ,dns0,mom0,ppara0,pperp0,qpara0,qperp0)
! 2016-02-04, for FLR case
! gc_x(:,n) = 1:r, 2:z, 3:phi, 4:p, 5:mu, 6:weight, 7:active, 8:p_phi, 9:fff, 10:fnrml
!--------------------------------------------------------------------
      use mpiset
      use grid
      use field, only:babs,babs0
      use gyro, only:ngyro,marker_each_gyro
      implicit none

      integer::marker_num,nvec,marker_num_gyro
      real(8)::gc(ngc2,marker_num)
      real(8)::gyro_phys(2,marker_num_gyro)
      real(8)::mass
      real(8)::dns(lr,lz,lphi),mom(lr,lz,lphi)
      real(8)::ppara(lr,lz,lphi),pperp(lr,lz,lphi)
      real(8)::qpara(lr,lz,lphi),qperp(lr,lz,lphi)
      real(8)::dns0(lr,lz,lphi),mom0(lr,lz,lphi)
      real(8)::ppara0(lr,lz,lphi),pperp0(lr,lz,lphi)
      real(8)::qpara0(lr,lz,lphi),qperp0(lr,lz,lphi)
      real(8)::wdns(lr,lz,lphi,lpara),wmom(lr,lz,lphi,lpara)
      real(8)::wpar(lr,lz,lphi,lpara),wprp(lr,lz,lphi,lpara)
      real(8)::wqar(lr,lz,lphi,lpara),wqrp(lr,lz,lphi,lpara)
!2015-09-21s
      integer::ijk_a(3,marker_each_gyro)
      real(8)::aaa(8,marker_each_gyro)
!2015-09-21e
      integer::i,j,k,l,m,n,ia,ia1,ja,ja1,ka,ka1,lr1,lz1
      real(8)::ar,ar1,az,az1,aphi,aphi1
!      real(8)::www,aaa1,aaa2,aaa3,aaa4,aaa5,aaa6,aaa7,aaa8
      real(8)::vol,cwwa,cwwc
!      real(8)::d0,d1,d2,d3,d4,d5,t1,t2,t3,t4,t5
      real(8)::dr1,dz1,dphi1,ma_mi_r
      integer::kr,kz,kphi
      real(8)::p0,p1,p2,mu1,p3,p1mu,mass1 !2015-09-21
      integer::n_min,n_max !2013-05-22
      real(8)::rngr1
      integer::nc
      integer::nvec0
      integer::nmod

      mass1 = 1.0d0/mass

      dr1 = 1.0d0/dr
      dz1 = 1.0d0/dz
      dphi1 = 1.0d0/dphi
      rngr1 = 1.0d0/dble(ngyro)

      kr   = 1 - kr_offset(my_rank)
      kz   = 1 - kz_offset(my_rank)
      kphi = 1 + lphi_shd - kphi_offset(my_rank)
!      kphi = 3 - kphi_offset(my_rank)
      ma_mi_r = major_r-minor_r

      if(lpara.ge.2)then
!$omp parallel do
      do k = 1, lpara
        do i = 1, lrzphi
          wdns(i,1,1,k) = 0.0d0
          wmom(i,1,1,k) = 0.0d0
          wpar(i,1,1,k) = 0.0d0
          wprp(i,1,1,k) = 0.0d0
          wqar(i,1,1,k) = 0.0d0
          wqrp(i,1,1,k) = 0.0d0
        end do
      end do

      end if

!$omp parallel do
      do i = 1, lrzphi
        dns(i,1,1) = 0.0d0
        mom(i,1,1) = 0.0d0
        ppara(i,1,1) = 0.0d0
        pperp(i,1,1) = 0.0d0
        qpara(i,1,1) = 0.0d0
        qperp(i,1,1) = 0.0d0
      end do

!$omp parallel private(nc,ar1,ar,az1,az,aphi1,aphi)
!poption parallel,schedule(dynamic,1000)
!$omp do schedule(dynamic,1000)
      do n = 1, marker_num_gyro
        nc =(n-1)/ngyro + 1
        ijk_a(1,n)=max(1,min(lr  -1,int((gyro_phys(1,n)-ma_mi_r)*dr1  ) + kr  ))
        ijk_a(2,n)=max(1,min(lz  -1,int(gyro_phys(2,n)          *dz1  ) + kz  ))
        ijk_a(3,n)=max(1,min(lphi-1,int(gc(3,nc)          *dphi1) + kphi))

        ar1  = max(0.0d0, min(1.0d0, (gyro_phys(1,n) - ma_mi_r)*dr1 -  dble(ijk_a(1,n) - kr)  ) )
        ar   = 1.0d0 - ar1
        az1  = max(0.0d0, min(1.0d0, gyro_phys(2,n)*dz1 - dble(ijk_a(2,n) - kz) ) )
        az   = 1.0d0 - az1
        aphi1 = max(0.0d0, min(1.0d0, gc(3,nc)*dphi1 - dble(ijk_a(3,n) - kphi) ) )
        aphi = 1.0d0 - aphi1

        aaa(1,n) = ar *az *aphi
        aaa(2,n) = ar1*az *aphi
        aaa(3,n) = ar *az1*aphi
        aaa(4,n) = ar1*az1*aphi
        aaa(5,n) = ar *az *aphi1
        aaa(6,n) = ar1*az *aphi1
        aaa(7,n) = ar *az1*aphi1
        aaa(8,n) = ar1*az1*aphi1
      end do
!$omp end parallel

      nvec = marker_num_gyro/lpara
      nvec0 = marker_num_gyro/lpara
      nmod  = mod(marker_num_gyro,lpara)

      if(lpara.ge.2)then

!2017-02-10
!$omp parallel private(nvec,m,n,nc,ia,ja,ka,p0,p1,p2,mu1,p3,p1mu)
!$omp do
      do l = 1, lpara
        nvec = nvec0 + min(1,nmod/l)
      do m = 1, nvec
        n = m +  nvec0*(l-1) + min(l-1,nmod)
        nc =(n-1)/ngyro + 1

        ia=ijk_a(1,n)
        ja=ijk_a(2,n)
        ka=ijk_a(3,n)

! gc_x(:,n) = 1:r, 2:z, 3:phi, 4:p, 5:mu, 6:weight, 7:active, 8:p_phi, 9:fff, 10:fnrml
        p0 =  gc(6,nc)*gc(7,nc)*rngr1 !rngr1 is multiplied
        p1 =  gc(4,nc)*p0
        p2 =  gc(4,nc)*p1
        mu1 = gc(5,nc)*p0
        p3  = gc(4,nc)*p2
        p1mu= gc(4,nc)*mu1

        wdns(ia,  ja,  ka,  l) = wdns(ia,  ja,  ka,  l) + aaa(1,n)*p0
        wdns(ia+1,ja,  ka,  l) = wdns(ia+1,ja,  ka,  l) + aaa(2,n)*p0
        wdns(ia,  ja+1,ka,  l) = wdns(ia,  ja+1,ka,  l) + aaa(3,n)*p0
        wdns(ia+1,ja+1,ka,  l) = wdns(ia+1,ja+1,ka,  l) + aaa(4,n)*p0
        wdns(ia,  ja,  ka+1,l) = wdns(ia,  ja,  ka+1,l) + aaa(5,n)*p0
        wdns(ia+1,ja,  ka+1,l) = wdns(ia+1,ja,  ka+1,l) + aaa(6,n)*p0
        wdns(ia,  ja+1,ka+1,l) = wdns(ia,  ja+1,ka+1,l) + aaa(7,n)*p0
        wdns(ia+1,ja+1,ka+1,l) = wdns(ia+1,ja+1,ka+1,l) + aaa(8,n)*p0

        wmom(ia,  ja,  ka,  l) = wmom(ia,  ja,  ka,  l) + aaa(1,n)*p1
        wmom(ia+1,ja,  ka,  l) = wmom(ia+1,ja,  ka,  l) + aaa(2,n)*p1
        wmom(ia,  ja+1,ka,  l) = wmom(ia,  ja+1,ka,  l) + aaa(3,n)*p1
        wmom(ia+1,ja+1,ka,  l) = wmom(ia+1,ja+1,ka,  l) + aaa(4,n)*p1
        wmom(ia,  ja,  ka+1,l) = wmom(ia,  ja,  ka+1,l) + aaa(5,n)*p1
        wmom(ia+1,ja,  ka+1,l) = wmom(ia+1,ja,  ka+1,l) + aaa(6,n)*p1
        wmom(ia,  ja+1,ka+1,l) = wmom(ia,  ja+1,ka+1,l) + aaa(7,n)*p1
        wmom(ia+1,ja+1,ka+1,l) = wmom(ia+1,ja+1,ka+1,l) + aaa(8,n)*p1

        wpar(ia,  ja,  ka,  l) = wpar(ia,  ja,  ka,  l) + aaa(1,n)*p2
        wpar(ia+1,ja,  ka,  l) = wpar(ia+1,ja,  ka,  l) + aaa(2,n)*p2
        wpar(ia,  ja+1,ka,  l) = wpar(ia,  ja+1,ka,  l) + aaa(3,n)*p2
        wpar(ia+1,ja+1,ka,  l) = wpar(ia+1,ja+1,ka,  l) + aaa(4,n)*p2
        wpar(ia,  ja,  ka+1,l) = wpar(ia,  ja,  ka+1,l) + aaa(5,n)*p2
        wpar(ia+1,ja,  ka+1,l) = wpar(ia+1,ja,  ka+1,l) + aaa(6,n)*p2
        wpar(ia,  ja+1,ka+1,l) = wpar(ia,  ja+1,ka+1,l) + aaa(7,n)*p2
        wpar(ia+1,ja+1,ka+1,l) = wpar(ia+1,ja+1,ka+1,l) + aaa(8,n)*p2

        wprp(ia,  ja,  ka,  l) = wprp(ia,  ja,  ka,  l) + aaa(1,n)*mu1
        wprp(ia+1,ja,  ka,  l) = wprp(ia+1,ja,  ka,  l) + aaa(2,n)*mu1
        wprp(ia,  ja+1,ka,  l) = wprp(ia,  ja+1,ka,  l) + aaa(3,n)*mu1
        wprp(ia+1,ja+1,ka,  l) = wprp(ia+1,ja+1,ka,  l) + aaa(4,n)*mu1
        wprp(ia,  ja,  ka+1,l) = wprp(ia,  ja,  ka+1,l) + aaa(5,n)*mu1
        wprp(ia+1,ja,  ka+1,l) = wprp(ia+1,ja,  ka+1,l) + aaa(6,n)*mu1
        wprp(ia,  ja+1,ka+1,l) = wprp(ia,  ja+1,ka+1,l) + aaa(7,n)*mu1
        wprp(ia+1,ja+1,ka+1,l) = wprp(ia+1,ja+1,ka+1,l) + aaa(8,n)*mu1

        wqar(ia,  ja,  ka,  l) = wqar(ia,  ja,  ka,  l) + aaa(1,n)*p3
        wqar(ia+1,ja,  ka,  l) = wqar(ia+1,ja,  ka,  l) + aaa(2,n)*p3
        wqar(ia,  ja+1,ka,  l) = wqar(ia,  ja+1,ka,  l) + aaa(3,n)*p3
        wqar(ia+1,ja+1,ka,  l) = wqar(ia+1,ja+1,ka,  l) + aaa(4,n)*p3
        wqar(ia,  ja,  ka+1,l) = wqar(ia,  ja,  ka+1,l) + aaa(5,n)*p3
        wqar(ia+1,ja,  ka+1,l) = wqar(ia+1,ja,  ka+1,l) + aaa(6,n)*p3
        wqar(ia,  ja+1,ka+1,l) = wqar(ia,  ja+1,ka+1,l) + aaa(7,n)*p3
        wqar(ia+1,ja+1,ka+1,l) = wqar(ia+1,ja+1,ka+1,l) + aaa(8,n)*p3

        wqrp(ia,  ja,  ka,  l) = wqrp(ia,  ja,  ka,  l) + aaa(1,n)*p1mu
        wqrp(ia+1,ja,  ka,  l) = wqrp(ia+1,ja,  ka,  l) + aaa(2,n)*p1mu
        wqrp(ia,  ja+1,ka,  l) = wqrp(ia,  ja+1,ka,  l) + aaa(3,n)*p1mu
        wqrp(ia+1,ja+1,ka,  l) = wqrp(ia+1,ja+1,ka,  l) + aaa(4,n)*p1mu
        wqrp(ia,  ja,  ka+1,l) = wqrp(ia,  ja,  ka+1,l) + aaa(5,n)*p1mu
        wqrp(ia+1,ja,  ka+1,l) = wqrp(ia+1,ja,  ka+1,l) + aaa(6,n)*p1mu
        wqrp(ia,  ja+1,ka+1,l) = wqrp(ia,  ja+1,ka+1,l) + aaa(7,n)*p1mu
        wqrp(ia+1,ja+1,ka+1,l) = wqrp(ia+1,ja+1,ka+1,l) + aaa(8,n)*p1mu

      end do
      end do
!$omp end parallel

      else

      do n = 1, marker_num_gyro
        nc =(n-1)/ngyro + 1

        ia=ijk_a(1,n)
        ja=ijk_a(2,n)
        ka=ijk_a(3,n)

! gc_x(:,n) = 1:r, 2:z, 3:phi, 4:p, 5:mu, 6:weight, 7:active, 8:p_phi, 9:fff, 10:fnrml
        p0 =  gc(6,nc)*gc(7,nc)*rngr1 !rngr1 is multiplied
        p1 =  gc(4,nc)*p0
        p2 =  gc(4,nc)*p1
        mu1 = gc(5,nc)*p0
        p3  = gc(4,nc)*p2
        p1mu= gc(4,nc)*mu1

        dns(ia,  ja,  ka  ) = dns(ia,  ja,  ka  ) + aaa(1,n)*p0
        dns(ia+1,ja,  ka  ) = dns(ia+1,ja,  ka  ) + aaa(2,n)*p0
        dns(ia,  ja+1,ka  ) = dns(ia,  ja+1,ka  ) + aaa(3,n)*p0
        dns(ia+1,ja+1,ka  ) = dns(ia+1,ja+1,ka  ) + aaa(4,n)*p0
        dns(ia,  ja,  ka+1) = dns(ia,  ja,  ka+1) + aaa(5,n)*p0
        dns(ia+1,ja,  ka+1) = dns(ia+1,ja,  ka+1) + aaa(6,n)*p0
        dns(ia,  ja+1,ka+1) = dns(ia,  ja+1,ka+1) + aaa(7,n)*p0
        dns(ia+1,ja+1,ka+1) = dns(ia+1,ja+1,ka+1) + aaa(8,n)*p0

        mom(ia,  ja,  ka  ) = mom(ia,  ja,  ka  ) + aaa(1,n)*p1
        mom(ia+1,ja,  ka  ) = mom(ia+1,ja,  ka  ) + aaa(2,n)*p1
        mom(ia,  ja+1,ka  ) = mom(ia,  ja+1,ka  ) + aaa(3,n)*p1
        mom(ia+1,ja+1,ka  ) = mom(ia+1,ja+1,ka  ) + aaa(4,n)*p1
        mom(ia,  ja,  ka+1) = mom(ia,  ja,  ka+1) + aaa(5,n)*p1
        mom(ia+1,ja,  ka+1) = mom(ia+1,ja,  ka+1) + aaa(6,n)*p1
        mom(ia,  ja+1,ka+1) = mom(ia,  ja+1,ka+1) + aaa(7,n)*p1
        mom(ia+1,ja+1,ka+1) = mom(ia+1,ja+1,ka+1) + aaa(8,n)*p1

        ppara(ia,  ja,  ka  ) = ppara(ia,  ja,  ka  ) + aaa(1,n)*p2
        ppara(ia+1,ja,  ka  ) = ppara(ia+1,ja,  ka  ) + aaa(2,n)*p2
        ppara(ia,  ja+1,ka  ) = ppara(ia,  ja+1,ka  ) + aaa(3,n)*p2
        ppara(ia+1,ja+1,ka  ) = ppara(ia+1,ja+1,ka  ) + aaa(4,n)*p2
        ppara(ia,  ja,  ka+1) = ppara(ia,  ja,  ka+1) + aaa(5,n)*p2
        ppara(ia+1,ja,  ka+1) = ppara(ia+1,ja,  ka+1) + aaa(6,n)*p2
        ppara(ia,  ja+1,ka+1) = ppara(ia,  ja+1,ka+1) + aaa(7,n)*p2
        ppara(ia+1,ja+1,ka+1) = ppara(ia+1,ja+1,ka+1) + aaa(8,n)*p2

        pperp(ia,  ja,  ka  ) = pperp(ia,  ja,  ka  ) + aaa(1,n)*mu1
        pperp(ia+1,ja,  ka  ) = pperp(ia+1,ja,  ka  ) + aaa(2,n)*mu1
        pperp(ia,  ja+1,ka  ) = pperp(ia,  ja+1,ka  ) + aaa(3,n)*mu1
        pperp(ia+1,ja+1,ka  ) = pperp(ia+1,ja+1,ka  ) + aaa(4,n)*mu1
        pperp(ia,  ja,  ka+1) = pperp(ia,  ja,  ka+1) + aaa(5,n)*mu1
        pperp(ia+1,ja,  ka+1) = pperp(ia+1,ja,  ka+1) + aaa(6,n)*mu1
        pperp(ia,  ja+1,ka+1) = pperp(ia,  ja+1,ka+1) + aaa(7,n)*mu1
        pperp(ia+1,ja+1,ka+1) = pperp(ia+1,ja+1,ka+1) + aaa(8,n)*mu1

        qpara(ia,  ja,  ka  ) = qpara(ia,  ja,  ka  ) + aaa(1,n)*p3
        qpara(ia+1,ja,  ka  ) = qpara(ia+1,ja,  ka  ) + aaa(2,n)*p3
        qpara(ia,  ja+1,ka  ) = qpara(ia,  ja+1,ka  ) + aaa(3,n)*p3
        qpara(ia+1,ja+1,ka  ) = qpara(ia+1,ja+1,ka  ) + aaa(4,n)*p3
        qpara(ia,  ja,  ka+1) = qpara(ia,  ja,  ka+1) + aaa(5,n)*p3
        qpara(ia+1,ja,  ka+1) = qpara(ia+1,ja,  ka+1) + aaa(6,n)*p3
        qpara(ia,  ja+1,ka+1) = qpara(ia,  ja+1,ka+1) + aaa(7,n)*p3
        qpara(ia+1,ja+1,ka+1) = qpara(ia+1,ja+1,ka+1) + aaa(8,n)*p3

        qperp(ia,  ja,  ka  ) = qperp(ia,  ja,  ka  ) + aaa(1,n)*p1mu
        qperp(ia+1,ja,  ka  ) = qperp(ia+1,ja,  ka  ) + aaa(2,n)*p1mu
        qperp(ia,  ja+1,ka  ) = qperp(ia,  ja+1,ka  ) + aaa(3,n)*p1mu
        qperp(ia+1,ja+1,ka  ) = qperp(ia+1,ja+1,ka  ) + aaa(4,n)*p1mu
        qperp(ia,  ja,  ka+1) = qperp(ia,  ja,  ka+1) + aaa(5,n)*p1mu
        qperp(ia+1,ja,  ka+1) = qperp(ia+1,ja,  ka+1) + aaa(6,n)*p1mu
        qperp(ia,  ja+1,ka+1) = qperp(ia,  ja+1,ka+1) + aaa(7,n)*p1mu
        qperp(ia+1,ja+1,ka+1) = qperp(ia+1,ja+1,ka+1) + aaa(8,n)*p1mu

      end do

      end if


      if(lpara.ge.2)then

      do k = 1, lpara
!$omp parallel do
      do i = 1, lrzphi
        dns(i,1,1)   = dns(i,1,1)   + wdns(i,1,1,k)
        mom(i,1,1)   = mom(i,1,1)   + wmom(i,1,1,k)
        ppara(i,1,1) = ppara(i,1,1) + wpar(i,1,1,k)
        pperp(i,1,1) = pperp(i,1,1) + wprp(i,1,1,k)
        qpara(i,1,1) = qpara(i,1,1) + wqar(i,1,1,k)
        qperp(i,1,1) = qperp(i,1,1) + wqrp(i,1,1,k)
      end do
      end do

      end if


! smoothing

       cwwa = 0.5d0
!       cwwc =-1.d0/6.d0 

!2024-12-21s, correction suggested by Panith Adulsiriswad
       call periodic_particle_mlt6b(dns,mom,ppara,pperp,qpara,qperp)
       call partsm1(dns,cwwa)
       call partsm1(mom,cwwa)

       call partsm1(ppara,cwwa)
       call partsm1(pperp,cwwa)

       call partsm1(qpara,cwwa)
       call partsm1(qperp,cwwa)
!2024-12-21e

! calculate density (per volume)

!      do k = 1, lphi
!      do j = 1, lz
!      do i = 1, lr
!$omp parallel do private(vol)
      do i = 1, lrzphi
        vol = 1.0d0/(grr(i,1,1)*dr*dz*dphi)
        dns(i,1,1) = dns(i,1,1)*vol
        mom(i,1,1) = mom(i,1,1)*vol
        ppara(i,1,1) = ppara(i,1,1)*vol*mass1
        pperp(i,1,1) = pperp(i,1,1)*vol*babs(i,1,1)
        qpara(i,1,1) = qpara(i,1,1)*vol*mass1**2
        qperp(i,1,1) = qperp(i,1,1)*vol*babs(i,1,1)*mass1
      end do
!      end do
!      end do



! wall effect

      lr1 = lr - 1
      lz1 = lz - 1

      if(my_rank_r.eq.(mpi_proc_r-1))then
!      do k = 1, lphi
!      do j = 1, lz
!$omp parallel do
       do j = 1, lzphi
          dns(lr1,j,1) = dns(lr1,j,1) + dns(lr,j,1)
          mom(lr1,j,1) = mom(lr1,j,1) + mom(lr,j,1)
          dns(lr,j,1) = 0.0d0
          mom(lr,j,1) = 0.0d0

          ppara(lr1,j,1) = ppara(lr1,j,1) + ppara(lr,j,1)
          pperp(lr1,j,1) = pperp(lr1,j,1) + pperp(lr,j,1)
          ppara(lr,j,1) = 0.0d0
          pperp(lr,j,1) = 0.0d0

          qpara(lr1,j,1) = qpara(lr1,j,1) + qpara(lr,j,1)
          qperp(lr1,j,1) = qperp(lr1,j,1) + qperp(lr,j,1)
          qpara(lr,j,1) = 0.0d0
          qperp(lr,j,1) = 0.0d0
       end do
!      end do
      end if

      if(my_rank_r.eq.0)then
!      do k = 1, lphi
!      do j = 1, lz
!$omp parallel do
       do j = 1, lzphi
          dns(2,j,1) = dns(2,j,1) + dns(1,j,1)
          mom(2,j,1) = mom(2,j,1) + mom(1,j,1)
          dns(1,j,1) = 0.0d0
          mom(1,j,1) = 0.0d0

          ppara(2,j,1) = ppara(2,j,1) + ppara(1,j,1)
          pperp(2,j,1) = pperp(2,j,1) + pperp(1,j,1)
          ppara(1,j,1) = 0.0d0
          pperp(1,j,1) = 0.0d0

          qpara(2,j,1) = qpara(2,j,1) + qpara(1,j,1)
          qperp(2,j,1) = qperp(2,j,1) + qperp(1,j,1)
          qpara(1,j,1) = 0.0d0
          qperp(1,j,1) = 0.0d0
       end do
!      end do
      end if

      if(my_rank_z.eq.(mpi_proc_z-1))then
!$omp parallel do
       do k = 1, lphi
       do i = 1, lr
          dns(i,lz1,k) = dns(i,lz1,k) + dns(i,lz,k)
          mom(i,lz1,k) = mom(i,lz1,k) + mom(i,lz,k)
          dns(i,lz,k) = 0.0d0
          mom(i,lz,k) = 0.0d0

          ppara(i,lz1,k) = ppara(i,lz1,k) + ppara(i,lz,k)
          pperp(i,lz1,k) = pperp(i,lz1,k) + pperp(i,lz,k)
          ppara(i,lz,k) = 0.0d0
          pperp(i,lz,k) = 0.0d0

          qpara(i,lz1,k) = qpara(i,lz1,k) + qpara(i,lz,k)
          qperp(i,lz1,k) = qperp(i,lz1,k) + qperp(i,lz,k)
          qpara(i,lz,k) = 0.0d0
          qperp(i,lz,k) = 0.0d0
       end do
       end do
      end if 

      if(my_rank_z.eq.0)then
!$omp parallel do
       do k = 1, lphi
       do i = 1, lr
          dns(i,2,k) = dns(i,2,k) + dns(i,1,k)
          mom(i,2,k) = mom(i,2,k) + mom(i,1,k)
          dns(i,1,k) = 0.0d0
          mom(i,1,k) = 0.0d0

          ppara(i,2,k) = ppara(i,2,k) + ppara(i,1,k)
          pperp(i,2,k) = pperp(i,2,k) + pperp(i,1,k)
          ppara(i,1,k) = 0.0d0
          pperp(i,1,k) = 0.0d0

          qpara(i,2,k) = qpara(i,2,k) + qpara(i,1,k)
          qperp(i,2,k) = qperp(i,2,k) + qperp(i,1,k)
          qpara(i,1,k) = 0.0d0
          qperp(i,1,k) = 0.0d0
       end do
       end do
      end if 

!      if(.not.flag_stored)then
! take modes n_min <= n <=n_max
!      n_min = 1
!      n_max = 2
!      call lowpass2_mlt2(n_min,n_max,dns,mom)
!      call lowpass2_mlt2(n_min,n_max,ppara,pperp)
!      call lowpass2_mlt2(n_min,n_max,qpara,qperp)
!      end if

!2014-08-28s
!$omp parallel do
      do i = 1, lrzphi
        dns(i,1,1) = dns(i,1,1) + dns0(i,1,1)
        mom(i,1,1) = mom(i,1,1) + mom0(i,1,1)
        ppara(i,1,1) = ppara(i,1,1) + ppara0(i,1,1)
        pperp(i,1,1) = pperp(i,1,1) + pperp0(i,1,1)/babs0(i,1,1)*babs(i,1,1)
        qpara(i,1,1) = qpara(i,1,1) + qpara0(i,1,1)
        qperp(i,1,1) = qperp(i,1,1) + qperp0(i,1,1)/babs0(i,1,1)*babs(i,1,1)
      end do
!2014-08-28e

!            call wall_clock(d5)
!            t5 = d5 - d4
!            if(my_rank.eq.0)then
!            write(6,*)'density 5, t5=',t5
!            end if

end subroutine moments_gyro

end module fx100
!----------------------------------------------------------------------------

!----------------------------------------------------------------------------
module fugaku
contains
!--------------------------------------------------------------------
subroutine push(marker_num,amassp,ep &
               ,type,temp,valpha,deltav,clambda0,dclambda & !2012-06-17
               ,gc,dgc,v &
               ,cf_pphi,pphi_min,pphi_max &
               ,flp,ispecies,igyro)
! type=0: maxwellian, type=1: slowing down, type=2: beam
! type=3: beam with finite pitch angle width !2012-06-17
! igyro=0: w/o FLR, igyro=1: w/ FLR
! modified on 2015-06-23
!   gc_x(:,n) = 1:r, 2:z, 3:phi, 4:p, 5:mu, 6:weight, 7:active, 8:p_phi, 9:fff, 10:fnrml
!--------------------------------------------------------------------
      use mpiset
      use field
      use grid
!      use equi_sol, only:raxis
      use particle, only:nu_krook !2024-04-23
      implicit none

      integer::marker_num,type
      integer::ispecies,igyro
      real(8)::amassp,ep
      real(8)::gc(ngc2,marker_num)
      real(8)::dgc(ngc1,marker_num)
      real(8)::v(marker_num)
      real(8)::cf_pphi(0:lcfpphi),pphi_min,pphi_max !2016-02-04
      real(8)::temp,valpha,deltav
      real(8)::dr1,dz1,dphi1,ep1,amsp1
      integer::n,ia,ia1,ja,ja1,ka,ka1,i
      real(8)::ar,ar1,az,az1,aphi,aphi1
!      real(8)::aaa1,aaa2,aaa3,aaa4,aaa5,aaa6,aaa7,aaa8

!2015-06-04s
!      integer::ijk_a(3,marker_each)
!      real(8)::aaa(8,marker_each)
!2015-06-04e

      real(8)::flp(nflp,marker_num),detotal
      real(8)::b2,babse,babs0e,b21,bre,bze,bphie
      real(8)::b1,br1a,bz1a,bphi1a,b10,br10,bz10,bphi10,dbpar,dbpr,dbpz,dbpphi
      real(8)::rhopar,orbpr,orbpz,orbpphi,denom1
      real(8)::wvnmlr,wvnmlz,wvnmlphi
      real(8)::dgrdbr,dgrdbz,dgrdbphi,dppar,dppar1,dppar2 !2012-08-31
      real(8)::pvpar,w1nmlr,w1nmlz,w1nmlphi
      real(8)::psinrm,pactive,prof,dprofdpsi,energy0,pphi_n,rminor,bmax
      real(8)::dpphi_n,vpara,vpara0,dvpara,sigma
      real(8)::vlcpart,dvlcpart,vnrm,sd_factor
      real(8)::sqr_pi1
      real(8)::nudt,coef,vrat !2016-01-09
      real(8)::nu_krook_dt !2024-04-23

      integer::kr,kz,kphi
      real(8)::ma_mi_r
      integer::kfl_start,kfl

      real(8)::pt_factor,energy,clambda0,dclambda,clambda(marker_each) !2012-06-17
      real(8)::dwpsi(marker_each),dwenrc(marker_each) !2015-06-23
!      real(8)::psip(marker_each)
!      real(8)::dpsi_dr(marker_each),dpsi_dz(marker_each),dpsi_dphi(marker_each)
      real(8)::bphi1e(marker_each) !2016-01-09

!2020-11-22s
      integer,parameter:: NB=320
      integer:: in, ib
      integer:: ijk_b(NB,3)
      real(8):: bbb(NB,8)
!2020-11-22e

! time derivative of particle position and velocity

      dr1 = 1.0d0/dr
      dz1 = 1.0d0/dz
      dphi1 = 1.0d0/dphi
      ep1 = 1.0d0/ep
      amsp1 = 1.0d0/amassp
      sqr_pi1=1.0d0/sqrt(pi)

      nu_krook_dt = nu_krook * dt !2024-04-23

      kr   = 1 - kr_offset(my_rank)
      kz   = 1 - kz_offset(my_rank)
      kphi = 1 + lphi_shd - kphi_offset(my_rank)
      ma_mi_r = major_r-minor_r

! igyro=0: w/o FLR, igyro/=0: flp(1:7,:) is given in subr. extract_em
      if(igyro.eq.0)then
        kfl_start = 1
      else
        kfl_start = nflp_gyro + 1
      end if

!$omp parallel

      if(igyro.eq.0)then !w/o FLR
!$omp workshare
        fld(1,:,:,:) = er(:,:,:)
        fld(2,:,:,:) = ez(:,:,:)
        fld(3,:,:,:) = ephi(:,:,:)
        fld(4,:,:,:) = epara(:,:,:)
        fld(5,:,:,:) = br(:,:,:) - br0(:,:,:)
        fld(6,:,:,:) = bz(:,:,:) - bz0(:,:,:)
        fld(7,:,:,:) = bphi(:,:,:) - bphi0(:,:,:)
!$omp end workshare
      end if

!2025-02-05      
!$omp workshare
      fld(11,:,:,:)= gradbr(:,:,:)
      fld(12,:,:,:)= gradbz(:,:,:)
      fld(13,:,:,:)= gradbphi(:,:,:)
      fld(14,:,:,:)= curvbr(:,:,:)
      fld(15,:,:,:)= curvbz(:,:,:)
      fld(16,:,:,:)= curvbphi(:,:,:)
!$omp end workshare

      if(ispecies.eq.0)then
!$omp workshare
        fld(23,:,:,:)= dns_e0(:,:,:)
        fld(24,:,:,:)= dns_e0_r(:,:,:)
        fld(25,:,:,:)= dns_e0_z(:,:,:)
        fld(26,:,:,:)= dns_e0_phi(:,:,:)
        fld(27,:,:,:)= temp_e0(:,:,:)
        fld(28,:,:,:)= temp_e0_r(:,:,:)
        fld(29,:,:,:)= temp_e0_z(:,:,:)
        fld(30,:,:,:)= temp_e0_phi(:,:,:)
!$omp end workshare
      else if(ispecies.eq.1)then
!$omp workshare
        fld(23,:,:,:)= dns_i0(:,:,:)
        fld(24,:,:,:)= dns_i0_r(:,:,:)
        fld(25,:,:,:)= dns_i0_z(:,:,:)
        fld(26,:,:,:)= dns_i0_phi(:,:,:)
        fld(27,:,:,:)= temp_i0(:,:,:)
        fld(28,:,:,:)= temp_i0_r(:,:,:)
        fld(29,:,:,:)= temp_i0_z(:,:,:)
        fld(30,:,:,:)= temp_i0_phi(:,:,:)
!$omp end workshare
      else if(ispecies.eq.2)then
!$omp workshare
        fld(23,:,:,:)= dns_a0(:,:,:)
        fld(24,:,:,:)= dns_a0_r(:,:,:)
        fld(25,:,:,:)= dns_a0_z(:,:,:)
        fld(26,:,:,:)= dns_a0_phi(:,:,:)
        fld(27,:,:,:)= temp_a0(:,:,:)
        fld(28,:,:,:)= temp_a0_r(:,:,:)
        fld(29,:,:,:)= temp_a0_z(:,:,:)
        fld(30,:,:,:)= temp_a0_phi(:,:,:)
!$omp end workshare
      end if
!$omp end parallel


!2020-11-22s
      if(igyro.eq.0)then
!poption parallel,schedule(dynamic,1000)
!$omp parallel private(ia,ja,ka,ar1,ar,az1,az,aphi1,aphi,in,ib,n,ijk_b,bbb)
!$omp do
      do in = 1, marker_num, NB
      do ib=1, min(NB, marker_num-in+1)
        n = in+ib-1

        ia=max(1,min(lr  -1,int((gc(1,n)-ma_mi_r)*dr1  ) + kr  ))
        ja=max(1,min(lz  -1,int(gc(2,n)          *dz1  ) + kz  ))
        ka=max(1,min(lphi-1,int(gc(3,n)          *dphi1) + kphi))

        ar1 = (gc(1,n)-grr(ia,ja,ka ) ) *dr1
        ar  = 1.0d0 - ar1
        az1 = (gc(2,n)-gzz(ia,ja,ka ) ) *dz1
        az  = 1.0d0 - az1
        aphi1 = (gc(3,n)-gphi(ia,ja,ka ) ) *dphi1
        aphi  = 1.0d0 - aphi1

        bbb(ib,1) = ar *az *aphi
        bbb(ib,2) = ar1*az *aphi
        bbb(ib,3) = ar *az1*aphi
        bbb(ib,4) = ar1*az1*aphi
        bbb(ib,5) = ar *az *aphi1
        bbb(ib,6) = ar1*az *aphi1
        bbb(ib,7) = ar *az1*aphi1
        bbb(ib,8) = ar1*az1*aphi1
        ijk_b(ib,1) = ia
        ijk_b(ib,2) = ja
        ijk_b(ib,3) = ka
      end do

! fields at each particle position
! igyro=0: w/o FLR, igyro/=0: flp(1:7,:) is given in subr. extract_em
      do ib=1, min(NB, marker_num-in+1)
        n = in+ib-1
        ia = ijk_b(ib,1)
        ja = ijk_b(ib,2)
        ka = ijk_b(ib,3)
        do kfl = 1, nflp

        flp(kfl,n) = fld(kfl, ia, ja,  ka  )*bbb(ib,1) + fld(kfl, ia+1,ja,    ka)*bbb(ib,2) &
                   + fld(kfl, ia, ja+1,ka  )*bbb(ib,3) + fld(kfl, ia+1,ja+1,  ka)*bbb(ib,4) &
                   + fld(kfl, ia, ja,  ka+1)*bbb(ib,5) + fld(kfl, ia+1,ja,  ka+1)*bbb(ib,6) &
                   + fld(kfl, ia, ja+1,ka+1)*bbb(ib,7) + fld(kfl, ia+1,ja+1,ka+1)*bbb(ib,8)
        end do

      end do
      end do
!$omp end parallel
      else
!poption parallel,schedule(dynamic,1000)
!$omp parallel private(ia,ja,ka,ar1,ar,az1,az,aphi1,aphi,in,ib,n,ijk_b,bbb)
!$omp do
      do in = 1, marker_num, NB
      do ib=1, min(NB, marker_num-in+1)
        n = in+ib-1

        ia=max(1,min(lr  -1,int((gc(1,n)-ma_mi_r)*dr1  ) + kr  ))
        ja=max(1,min(lz  -1,int(gc(2,n)          *dz1  ) + kz  ))
        ka=max(1,min(lphi-1,int(gc(3,n)          *dphi1) + kphi))

        ar1 = (gc(1,n)-grr(ia,ja,ka ) ) *dr1
        ar  = 1.0d0 - ar1
        az1 = (gc(2,n)-gzz(ia,ja,ka ) ) *dz1
        az  = 1.0d0 - az1
        aphi1 = (gc(3,n)-gphi(ia,ja,ka ) ) *dphi1
        aphi  = 1.0d0 - aphi1

        bbb(ib,1) = ar *az *aphi
        bbb(ib,2) = ar1*az *aphi
        bbb(ib,3) = ar *az1*aphi
        bbb(ib,4) = ar1*az1*aphi
        bbb(ib,5) = ar *az *aphi1
        bbb(ib,6) = ar1*az *aphi1
        bbb(ib,7) = ar *az1*aphi1
        bbb(ib,8) = ar1*az1*aphi1
        ijk_b(ib,1) = ia
        ijk_b(ib,2) = ja
        ijk_b(ib,3) = ka
      enddo

! fields at each particle position
! igyro=0: w/o FLR, igyro/=0: flp(1:7,:) is given in subr. extract_em
      do ib=1, min(NB, marker_num-in+1)
        n = in+ib-1
        ia = ijk_b(ib,1)
        ja = ijk_b(ib,2)
        ka = ijk_b(ib,3)
        do kfl = nflp_gyro + 1, nflp

        flp(kfl,n) = fld(kfl, ia, ja,  ka  )*bbb(ib,1) + fld(kfl, ia+1,ja,    ka)*bbb(ib,2) &
                   + fld(kfl, ia, ja+1,ka  )*bbb(ib,3) + fld(kfl, ia+1,ja+1,  ka)*bbb(ib,4) &
                   + fld(kfl, ia, ja,  ka+1)*bbb(ib,5) + fld(kfl, ia+1,ja,  ka+1)*bbb(ib,6) &
                   + fld(kfl, ia, ja+1,ka+1)*bbb(ib,7) + fld(kfl, ia+1,ja+1,ka+1)*bbb(ib,8)
        end do

      end do
      end do
!$omp end parallel
      endif

!2020-11-22e

!$omp parallel private(bre,bze,bphie,b2,babse,babs0e,b1,b21,br1a,bz1a,bphi1a &
!$omp& ,b10,br10,bz10,bphi10,dbpar,dbpr,dbpz,dbpphi,rhopar,denom1,orbpr,orbpz,orbpphi &
!$omp& ,wvnmlr,wvnmlz,wvnmlphi,dgrdbr,dgrdbz,dgrdbphi,dppar,dppar2,dppar1,pvpar &
!$omp& ,w1nmlr,w1nmlz,w1nmlphi,detotal,energy &
!$omp& ,rminor,bmax,energy0,sigma,vpara0,vpara,dvpara &
!$omp& ,pphi_n,dpphi_n,prof,dprofdpsi)
!$omp do
      do n = 1, marker_num
! flp(5:7,n): delta_br(z,phi)
        bre = flp(5,n) + flp(8,n)
        bze = flp(6,n) + flp(9,n)
        bphie = flp(7,n) + flp(10,n)

        b2  = bre**2 + bze**2 + bphie**2
        babse= max(eps_b, sqrt(b2) )
        babs0e= max(eps_b, sqrt(flp(8,n)**2 + flp(9,n)**2 + flp(10,n)**2) )
        b1 = 1.0d0/babse
        b21= 1.0d0/b2
        br1a = bre*b1
        bz1a = bze*b1
        bphi1a = bphie*b1
        bphi1e(n) = bphi1a !2016-02-4

        b10 = 1.0d0/babs0e
        br10 = flp(8,n)*b10
        bz10 = flp(9,n)*b10
        bphi10 = flp(10,n)*b10

        dbpar = br1a*br10 + bz1a*bz10 + bphi1a*bphi10
        dbpr  = br1a  - br10*dbpar
        dbpz  = bz1a  - bz10*dbpar
        dbpphi= bphi1a- bphi10*dbpar

! guiding center motion

        rhopar = gc(4,n)*ep1*b1

        denom1 = 1.0d0/(1.0d0 + rhopar*(br1a*flp(14,n) + bz1a*flp(15,n) + bphi1a*flp(16,n)))

        orbpr = (br1a + rhopar*flp(14,n))*denom1
        orbpz = (bz1a + rhopar*flp(15,n))*denom1
        orbpphi = (bphi1a + rhopar*flp(16,n))*denom1

! e x b drift

        wvnmlr = (flp(3,n)*bze-flp(2,n)*bphie)*b21*denom1
        wvnmlz = (flp(1,n)*bphie-flp(3,n)*bre)*b21*denom1
        wvnmlphi = (flp(2,n)*bre-flp(1,n)*bze)*b21*denom1

! grad-b drift

        dgrdbr = gc(5,n)*(bphie*flp(12,n) - bze*flp(13,n))*ep1*b21*denom1
        dgrdbz = gc(5,n)*(bre*flp(13,n) - bphie*flp(11,n))*ep1*b21*denom1
        dgrdbphi = gc(5,n)*(bze*flp(11,n) - bre*flp(12,n))*ep1*b21*denom1

! mirror force

        dppar =-gc(5,n)*( flp(11,n)*orbpr &
                      + flp(12,n)*orbpz &
                      + flp(13,n)*orbpphi &
                      )*dt

!2012-08-31
        dppar2=-gc(5,n)*( flp(11,n)*dbpr*denom1 &
                      + flp(12,n)*dbpz*denom1 &
                      + flp(13,n)*orbpphi &
                      )*dt
!2012-08-31 end

! aceeleration due to electric field and curvature drift

!2012-07-07
        dppar1 = ep*((flp(1,n)*flp(14,n) + flp(2,n)*flp(15,n) + flp(3,n)*flp(16,n) &
                      )*rhopar &
                    + flp(4,n) &
                     )*dt*denom1
!2012-07-07 end

! total drift velocity

        pvpar = gc(4,n)*amsp1
        dgc(1,n) = dt*(pvpar*orbpr + wvnmlr + dgrdbr)*gc(7,n)
        dgc(2,n) = dt*(pvpar*orbpz + wvnmlz + dgrdbz)*gc(7,n)
        dgc(3,n) = dt*(pvpar*orbpphi + wvnmlphi + dgrdbphi)/gc(1,n)*gc(7,n)
        dgc(4,n) =(dppar + dppar1)*gc(7,n)

! temporal evolution of weight of high-energy ion particle

        w1nmlr = wvnmlr + pvpar*dbpr*denom1
        w1nmlz = wvnmlz + pvpar*dbpz*denom1
        w1nmlphi = wvnmlphi + pvpar*dbpphi*denom1

        dgc(8,n) =(ep*(w1nmlr*flp(18,n) + w1nmlz*flp(19,n) + w1nmlphi*flp(20,n))*dt &
                  + bphi1a*(w1nmlr*gc(4,n)*dt + gc(1,n)*(dppar1 + dppar2) ) & !2012-08-31
                  )*gc(7,n)

        detotal =(  ep*(flp(1,n)*dgrdbr + flp(2,n)*dgrdbz + flp(3,n)*dgrdbphi)*dt &
                + dppar1*pvpar )*gc(7,n) &
! the following term considers 'mu* v * grad(B0 - B)', 2025-04-04
                + gc(5,n)*(dgc(1,n)*(flp(31,n)-flp(11,n) ) &
                          +dgc(2,n)*(flp(32,n)-flp(12,n) ) &
                          +dgc(3,n)*(flp(33,n)-flp(13,n) )*gc(1,n) &
                          )
!2025-04-27s
        energy0 = max(1.0d-30, 0.50d0*amassp*pvpar**2 + gc(5,n)*babs0e)
        clambda(n) = gc(5,n)*b0/energy0
        v(n) = sqrt(2.0d0*energy0*amsp1)

!        energy = max(1.0d-30, 0.50d0*amassp*pvpar**2 + gc(5,n)*babse)
!        clambda(n) = gc(5,n)*b0/energy
!        v(n) = sqrt(2.0d0*energy*amsp1)
!2025-04-27e

! weight evolution : weight = f - f0
! d weight /dt = - d f0 /dt
! f0 = f_nrml*prof(psinrm)/(v**3 + flp(21,n)**3)*0.50*erfc((v(n)-valpha)/deltav)

!        rminor = sqrt( (gc(1,n)-raxis)**2 &
!                     + ((gc(2,n)-0.50d0*zleng) )**2 &
!                     )
!        bmax = b0*raxis/(raxis-rminor)
!        energy0 = 0.50d0*amsp1*gc(4,n)**2 + gc(5,n)*babs0e
!        sigma = 0.50d0*(1.0d0 + sign(1.0d0, energy0-gc(5,n)*bmax) ) &
!                      *sign(1.0d0,gc(4,n) )
!        vpara0 = sqrt(2.0d0*(energy0-gc(5,n)*bmin)*amsp1)
!        vpara = vpara0*sigma
!        dvpara = sigma/vpara0*gc(4,n)*dppar1*amsp1**2

!        pphi_n = gc(8,n) - amassp*raxis*vpara
!        dpphi_n= dgc(8,n) - amassp*raxis*dvpara

!        prof = exp(pphi_n/(ep*psimax*0.37d0) )
!        dprofdpsi = prof/(ep*psimax*0.37d0)

!        dwpsi(n) = dpphi_n*dprofdpsi*gc(10,n)
!        dwenrc(n)= detotal/(amassp*v(n) )*prof*gc(10,n)

!2016-08-05s
        dwpsi(n) = ( w1nmlr*flp(24,n) + w1nmlz*flp(25,n) + w1nmlphi*flp(26,n) &
                    +(w1nmlr*flp(28,n) + w1nmlz*flp(29,n) + w1nmlphi*flp(30,n) ) &
                    *0.50d0*(amassp*v(n)**2/flp(27,n) - 3.0d0)*flp(23,n)/flp(27,n) &
                   )*dt*gc(7,n)*gc(10,n)

!        dwpsi(n) = (w1nmlr*flp(24,n) + w1nmlz*flp(25,n) + w1nmlphi*flp(26,n))*dt & !2016-01-09
!                   *gc(7,n)*gc(10,n)
!2016-08-05e

        dwenrc(n)= detotal/(amassp*v(n) )*flp(23,n)*gc(10,n) !2016-01-09

      end do
!$omp end parallel


      if(type.eq.0)then

!$omp parallel do private(vlcpart,dvlcpart)
        do n = 1, marker_num
!2016-08-05s
          vlcpart = exp(-0.50d0*amassp*v(n)**2/flp(27,n) )*flp(27,n)**(-1.5d0) !2013-07-17
          dvlcpart = -amassp*v(n)/flp(27,n)*vlcpart 
!          vlcpart = exp(-0.50d0*amassp*v(n)**2/temp)
!          dvlcpart = -amassp*v(n)/temp*vlcpart 
!2016-08-05e
          dwpsi(n) = dwpsi(n)*vlcpart
          dwenrc(n) = dwenrc(n)*dvlcpart
          dgc(6,n) = - dwpsi(n) - dwenrc(n) &
                     - gc(6,n)*nu_krook_dt & !Krook operator at normalized psi > psi_edge, flp(17)=psi
                     *(0.50d0 + sign(0.50d0, 1.0d0 - flp(17,n)/psimax - psi_edge) ) !2024-04-23
        end do

      else if(type.eq.1.or.type.eq.2)then

!$omp parallel do private(vnrm,sd_factor,vlcpart,dvlcpart)
        do n = 1, marker_num
! for itpa benchmark Oct11, 2009

          vnrm =(v(n)-valpha)/deltav
          sd_factor = 1.0d0/(v(n)**3 + flp(21,n)**3)
!          vlcpart = 0.50d0*erfc(vnrm)*sd_factor
          vlcpart = 0.50d0*max(0.0d0, min(2.0d0, 1.0d0 - vnrm))*sd_factor !linear function
          dvlcpart= -( 3.0d0*v(n)**2*vlcpart &
!                     + exp(-vnrm**2)*sqr_pi1/deltav &
                     + (0.25d0 - sign(0.25d0, (vnrm+1.0d0)*(vnrm-1.0d0)))/deltav & !linear function
                     )*sd_factor

          dwpsi(n) = dwpsi(n)*vlcpart
          dwenrc(n) = dwenrc(n)*dvlcpart
          dgc(6,n) = - dwpsi(n) - dwenrc(n) &
                     - gc(6,n)*nu_krook_dt & !Krook operator at normalized psi > psi_edge, flp(17)=psi
                     *(0.50d0 + sign(0.50d0, 1.0d0 - flp(17,n)/psimax - psi_edge) ) !2024-04-23
        end do

!2012-06-17
      else if(type.eq.3)then

!$omp parallel do private(vnrm,sd_factor,pt_factor,vlcpart,dvlcpart)
        do n = 1, marker_num
          vnrm =(v(n)-valpha)/deltav
          sd_factor = 1.0d0/(v(n)**3 + flp(21,n)**3)
          pt_factor = exp(-(clambda(n)-clambda0)**2/dclambda**2)

!          vlcpart = 0.50d0*erfc(vnrm)*sd_factor*pt_factor
          vlcpart = 0.50d0*max(0.0d0, min(2.0d0, 1.0d0 - vnrm))*sd_factor*pt_factor !linear function
          dvlcpart= -( 3.0d0*v(n)**2*vlcpart &
!                     + exp(-vnrm**2)*sqr_pi1/deltav*pt_factor &
                     + (0.25d0 - sign(0.25d0, (vnrm+1.0d0)*(vnrm-1.0d0)))/deltav*pt_factor & !linear function
                     )*sd_factor &
                  + 4.0d0*vlcpart*clambda(n)*(clambda(n)-clambda0) &
                         /(v(n)*dclambda**2)

          dwpsi(n) = dwpsi(n)*vlcpart
          dwenrc(n) = dwenrc(n)*dvlcpart
          dgc(6,n) = - dwpsi(n) - dwenrc(n) &
                     - gc(6,n)*nu_krook_dt & !Krook operator at normalized psi > psi_edge, flp(17)=psi
                     *(0.50d0 + sign(0.50d0, 1.0d0 - flp(17,n)/psimax - psi_edge) ) !2024-04-23
        end do
!2012-06-17 end

      end if

! slowing down
      if(type.eq.-5)then
!$omp parallel do private(nudt,vrat,coef)
        do n = 1, marker_num
            dgc(6,n) = 0.0d0 !full-f
            nudt = flp(22,n)*dt*gc(7,n)
            vrat = flp(21,n)/v(n)
            coef =(1.0d0 + vrat**3)*nudt
            dgc(10,n) =-3.0d0*nudt*gc(10,n)
            dgc(5,n)= -2.0d0*gc(5,n)*coef
            dgc(4,n) = dgc(4,n)- gc(4,n)*coef
            dgc(8,n) = dgc(8,n) - gc(4,n)*coef*gc(1,n)*bphi1e(n)
        end do
      end if

end subroutine push
!--------------------------------------------------------------------
subroutine density(marker_num,mass,gc &
                  ,dns,mom,ppara,pperp,dns0,mom0,ppara0,pperp0)
! calculate pressure
! modified for Fugaku 2020-11-22
! gc_x(:,n) = 1:r, 2:z, 3:phi, 4:p, 5:mu, 6:weight, 7:active, 8:p_phi, 9:fff, 10:fnrml
!--------------------------------------------------------------------
      use mpiset
      use grid
      use field, only:babs,babs0
      implicit none

      integer::marker_num,nvec
      real(8)::gc(ngc2,marker_num)
      real(8)::mass
      real(8)::dns(lr,lz,lphi),mom(lr,lz,lphi)
      real(8)::ppara(lr,lz,lphi),pperp(lr,lz,lphi)
      real(8)::dns0(lr,lz,lphi),mom0(lr,lz,lphi)
      real(8)::ppara0(lr,lz,lphi),pperp0(lr,lz,lphi)
      real(8)::wdns(lr,lz,lphi,lpara),wmom(lr,lz,lphi,lpara)
      real(8)::wpar(lr,lz,lphi,lpara),wprp(lr,lz,lphi,lpara)
      integer::ijk_a(3,marker_each)
      real(8)::aaa(8,marker_each)
      integer::i,j,k,l,m,n,ia,ia1,ja,ja1,ka,ka1,lr1,lz1
      real(8)::ar,ar1,az,az1,aphi,aphi1
!      real(8)::www,aaa1,aaa2,aaa3,aaa4,aaa5,aaa6,aaa7,aaa8
      real(8)::vol,cwwa,cwwc
      real(8)::d0,d1,d2,d3,d4,d5,t1,t2,t3,t4,t5
      real(8)::dr1,dz1,dphi1,ma_mi_r
      integer::kr,kz,kphi
      real(8)::p0,p1,p2,mu1
      integer::n_min,n_max !2013-05-22
!2020-11-22s
      integer::nvec0
      integer::nmod
      integer,parameter:: NB=320
      real(8):: bbb(NB,8)
      integer::ijk_b(3,NB)
!2020-11-22e

      dr1 = 1.0d0/dr
      dz1 = 1.0d0/dz
      dphi1 = 1.0d0/dphi

      kr   = 1 - kr_offset(my_rank)
      kz   = 1 - kz_offset(my_rank)
      kphi = 1 + lphi_shd - kphi_offset(my_rank)
!      kphi = 3 - kphi_offset(my_rank)
      ma_mi_r = major_r-minor_r

      if(lpara.ge.2)then

!$omp parallel do
      do k = 1, lpara
        do i = 1, lrzphi
          wdns(i,1,1,k) = 0.0d0
          wmom(i,1,1,k) = 0.0d0
          wpar(i,1,1,k) = 0.0d0
          wprp(i,1,1,k) = 0.0d0
        end do
      end do

      end if

!2020-11-22s

!$omp parallel do
      do i = 1, lrzphi
        dns(i,1,1) = 0.0d0
        mom(i,1,1) = 0.0d0
        ppara(i,1,1) = 0.0d0
        pperp(i,1,1) = 0.0d0
      end do

      nvec = marker_num/lpara
      nvec0 = marker_num/lpara
      nmod  = mod(marker_num,lpara)

      if(lpara.ge.2)then

!$omp parallel private(n,ia,ja,ka,p0,p1,p2,mu1,nvec) &
!$omp& private(ar1,ar,az1,az,aphi1,aphi,ijk_b,bbb)
!$omp do
      do l = 1, lpara

        nvec = nvec0 + min(1,nmod/l)
        do m = 1, nvec, NB

        do i=1, min(NB, nvec-m+1)
        n = m+i-1 +  nvec0*(l-1) + min(l-1,nmod)

        ia=max(1,min(lr  -1,int((gc(1,n)-ma_mi_r)*dr1  ) + kr  ))
        ja=max(1,min(lz  -1,int(gc(2,n)          *dz1  ) + kz  ))
        ka=max(1,min(lphi-1,int(gc(3,n)        *dphi1) + kphi))
        ar1 = (gc(1,n)-grr(ia,ja,ka ) ) *dr1
        ar  = 1.0d0 - ar1
        az1 = (gc(2,n)-gzz(ia,ja,ka ) ) *dz1
        az  = 1.0d0 - az1
        aphi1 = (gc(3,n)-gphi(ia,ja,ka ) ) *dphi1
        aphi  = 1.0d0 - aphi1

        bbb(i,1) = ar *az *aphi
        bbb(i,2) = ar1*az *aphi
        bbb(i,3) = ar *az1*aphi
        bbb(i,4) = ar1*az1*aphi
        bbb(i,5) = ar *az *aphi1
        bbb(i,6) = ar1*az *aphi1
        bbb(i,7) = ar *az1*aphi1
        bbb(i,8) = ar1*az1*aphi1

        ijk_b(1,i)=ia
        ijk_b(2,i)=ja
        ijk_b(3,i)=ka
      enddo

      do i=1, min(NB, nvec-m+1)
        ia=ijk_b(1,i)
        ja=ijk_b(2,i)
        ka=ijk_b(3,i)
        n = m+i-1 +  nvec0*(l-1) + min(l-1,nmod)

! gc_x(:,n) = 1:r, 2:z, 3:phi, 4:p, 5:mu, 6:weight, 7:active, 8:p_phi, 9:fff, 10:fnrml
        p0 =  gc(6,n)*gc(7,n)
        p1 =  gc(4,n)*p0
        p2 =  gc(4,n)*p1
        mu1 = gc(5,n)*p0

        wdns(ia,  ja,  ka,  l) = wdns(ia,  ja,  ka,  l) + bbb(i,1)*p0
        wdns(ia+1,ja,  ka,  l) = wdns(ia+1,ja,  ka,  l) + bbb(i,2)*p0
        wdns(ia,  ja+1,ka,  l) = wdns(ia,  ja+1,ka,  l) + bbb(i,3)*p0
        wdns(ia+1,ja+1,ka,  l) = wdns(ia+1,ja+1,ka,  l) + bbb(i,4)*p0
        wdns(ia,  ja,  ka+1,l) = wdns(ia,  ja,  ka+1,l) + bbb(i,5)*p0
        wdns(ia+1,ja,  ka+1,l) = wdns(ia+1,ja,  ka+1,l) + bbb(i,6)*p0
        wdns(ia,  ja+1,ka+1,l) = wdns(ia,  ja+1,ka+1,l) + bbb(i,7)*p0
        wdns(ia+1,ja+1,ka+1,l) = wdns(ia+1,ja+1,ka+1,l) + bbb(i,8)*p0

        wmom(ia,  ja,  ka,  l) = wmom(ia,  ja,  ka,  l) + bbb(i,1)*p1
        wmom(ia+1,ja,  ka,  l) = wmom(ia+1,ja,  ka,  l) + bbb(i,2)*p1
        wmom(ia,  ja+1,ka,  l) = wmom(ia,  ja+1,ka,  l) + bbb(i,3)*p1
        wmom(ia+1,ja+1,ka,  l) = wmom(ia+1,ja+1,ka,  l) + bbb(i,4)*p1
        wmom(ia,  ja,  ka+1,l) = wmom(ia,  ja,  ka+1,l) + bbb(i,5)*p1
        wmom(ia+1,ja,  ka+1,l) = wmom(ia+1,ja,  ka+1,l) + bbb(i,6)*p1
        wmom(ia,  ja+1,ka+1,l) = wmom(ia,  ja+1,ka+1,l) + bbb(i,7)*p1
        wmom(ia+1,ja+1,ka+1,l) = wmom(ia+1,ja+1,ka+1,l) + bbb(i,8)*p1

        wpar(ia,  ja,  ka,  l) = wpar(ia,  ja,  ka,  l) + bbb(i,1)*p2
        wpar(ia+1,ja,  ka,  l) = wpar(ia+1,ja,  ka,  l) + bbb(i,2)*p2
        wpar(ia,  ja+1,ka,  l) = wpar(ia,  ja+1,ka,  l) + bbb(i,3)*p2
        wpar(ia+1,ja+1,ka,  l) = wpar(ia+1,ja+1,ka,  l) + bbb(i,4)*p2
        wpar(ia,  ja,  ka+1,l) = wpar(ia,  ja,  ka+1,l) + bbb(i,5)*p2
        wpar(ia+1,ja,  ka+1,l) = wpar(ia+1,ja,  ka+1,l) + bbb(i,6)*p2
        wpar(ia,  ja+1,ka+1,l) = wpar(ia,  ja+1,ka+1,l) + bbb(i,7)*p2
        wpar(ia+1,ja+1,ka+1,l) = wpar(ia+1,ja+1,ka+1,l) + bbb(i,8)*p2

        wprp(ia,  ja,  ka,  l) = wprp(ia,  ja,  ka,  l) + bbb(i,1)*mu1
        wprp(ia+1,ja,  ka,  l) = wprp(ia+1,ja,  ka,  l) + bbb(i,2)*mu1
        wprp(ia,  ja+1,ka,  l) = wprp(ia,  ja+1,ka,  l) + bbb(i,3)*mu1
        wprp(ia+1,ja+1,ka,  l) = wprp(ia+1,ja+1,ka,  l) + bbb(i,4)*mu1
        wprp(ia,  ja,  ka+1,l) = wprp(ia,  ja,  ka+1,l) + bbb(i,5)*mu1
        wprp(ia+1,ja,  ka+1,l) = wprp(ia+1,ja,  ka+1,l) + bbb(i,6)*mu1
        wprp(ia,  ja+1,ka+1,l) = wprp(ia,  ja+1,ka+1,l) + bbb(i,7)*mu1
        wprp(ia+1,ja+1,ka+1,l) = wprp(ia+1,ja+1,ka+1,l) + bbb(i,8)*mu1
      enddo

      end do
      end do
!$omp end parallel

      else
!2020-11-22e

!poption parallel,schedule(dynamic,1000)
!$omp parallel private(ar1,ar,az1,az,aphi1,aphi)
!$omp do schedule(dynamic,1000)
      do n = 1, marker_num
        ijk_a(1,n)=max(1,min(lr  -1,int((gc(1,n)-ma_mi_r)*dr1  ) + kr  ))
        ijk_a(2,n)=max(1,min(lz  -1,int(gc(2,n)          *dz1  ) + kz  ))
        ijk_a(3,n)=max(1,min(lphi-1,int(gc(3,n)        *dphi1) + kphi))

        ar1 = (gc(1,n)-grr(ijk_a(1,n),ijk_a(2,n),ijk_a(3,n) ) ) *dr1
        ar  = 1.0d0 - ar1
        az1 = (gc(2,n)-gzz(ijk_a(1,n),ijk_a(2,n),ijk_a(3,n) ) ) *dz1
        az  = 1.0d0 - az1
        aphi1 = (gc(3,n)-gphi(ijk_a(1,n),ijk_a(2,n),ijk_a(3,n) ) ) *dphi1
        aphi  = 1.0d0 - aphi1

        aaa(1,n) = ar *az *aphi
        aaa(2,n) = ar1*az *aphi
        aaa(3,n) = ar *az1*aphi
        aaa(4,n) = ar1*az1*aphi
        aaa(5,n) = ar *az *aphi1
        aaa(6,n) = ar1*az *aphi1
        aaa(7,n) = ar *az1*aphi1
        aaa(8,n) = ar1*az1*aphi1
      end do
!$omp end parallel

      do n = 1, marker_num

        ia=ijk_a(1,n)
        ja=ijk_a(2,n)
        ka=ijk_a(3,n)

! gc_x(:,n) = 1:r, 2:z, 3:phi, 4:p, 5:mu, 6:weight, 7:active, 8:p_phi, 9:fff, 10:fnrml
        p0 =  gc(6,n)*gc(7,n)
        p1 =  gc(4,n)*p0
        p2 =  gc(4,n)*p1
        mu1 = gc(5,n)*p0

        dns(ia,  ja,  ka  ) = dns(ia,  ja,  ka  ) + aaa(1,n)*p0
        dns(ia+1,ja,  ka  ) = dns(ia+1,ja,  ka  ) + aaa(2,n)*p0
        dns(ia,  ja+1,ka  ) = dns(ia,  ja+1,ka  ) + aaa(3,n)*p0
        dns(ia+1,ja+1,ka  ) = dns(ia+1,ja+1,ka  ) + aaa(4,n)*p0
        dns(ia,  ja,  ka+1) = dns(ia,  ja,  ka+1) + aaa(5,n)*p0
        dns(ia+1,ja,  ka+1) = dns(ia+1,ja,  ka+1) + aaa(6,n)*p0
        dns(ia,  ja+1,ka+1) = dns(ia,  ja+1,ka+1) + aaa(7,n)*p0
        dns(ia+1,ja+1,ka+1) = dns(ia+1,ja+1,ka+1) + aaa(8,n)*p0

        mom(ia,  ja,  ka  ) = mom(ia,  ja,  ka  ) + aaa(1,n)*p1
        mom(ia+1,ja,  ka  ) = mom(ia+1,ja,  ka  ) + aaa(2,n)*p1
        mom(ia,  ja+1,ka  ) = mom(ia,  ja+1,ka  ) + aaa(3,n)*p1
        mom(ia+1,ja+1,ka  ) = mom(ia+1,ja+1,ka  ) + aaa(4,n)*p1
        mom(ia,  ja,  ka+1) = mom(ia,  ja,  ka+1) + aaa(5,n)*p1
        mom(ia+1,ja,  ka+1) = mom(ia+1,ja,  ka+1) + aaa(6,n)*p1
        mom(ia,  ja+1,ka+1) = mom(ia,  ja+1,ka+1) + aaa(7,n)*p1
        mom(ia+1,ja+1,ka+1) = mom(ia+1,ja+1,ka+1) + aaa(8,n)*p1

        ppara(ia,  ja,  ka  ) = ppara(ia,  ja,  ka  ) + aaa(1,n)*p2
        ppara(ia+1,ja,  ka  ) = ppara(ia+1,ja,  ka  ) + aaa(2,n)*p2
        ppara(ia,  ja+1,ka  ) = ppara(ia,  ja+1,ka  ) + aaa(3,n)*p2
        ppara(ia+1,ja+1,ka  ) = ppara(ia+1,ja+1,ka  ) + aaa(4,n)*p2
        ppara(ia,  ja,  ka+1) = ppara(ia,  ja,  ka+1) + aaa(5,n)*p2
        ppara(ia+1,ja,  ka+1) = ppara(ia+1,ja,  ka+1) + aaa(6,n)*p2
        ppara(ia,  ja+1,ka+1) = ppara(ia,  ja+1,ka+1) + aaa(7,n)*p2
        ppara(ia+1,ja+1,ka+1) = ppara(ia+1,ja+1,ka+1) + aaa(8,n)*p2

        pperp(ia,  ja,  ka  ) = pperp(ia,  ja,  ka  ) + aaa(1,n)*mu1
        pperp(ia+1,ja,  ka  ) = pperp(ia+1,ja,  ka  ) + aaa(2,n)*mu1
        pperp(ia,  ja+1,ka  ) = pperp(ia,  ja+1,ka  ) + aaa(3,n)*mu1
        pperp(ia+1,ja+1,ka  ) = pperp(ia+1,ja+1,ka  ) + aaa(4,n)*mu1
        pperp(ia,  ja,  ka+1) = pperp(ia,  ja,  ka+1) + aaa(5,n)*mu1
        pperp(ia+1,ja,  ka+1) = pperp(ia+1,ja,  ka+1) + aaa(6,n)*mu1
        pperp(ia,  ja+1,ka+1) = pperp(ia,  ja+1,ka+1) + aaa(7,n)*mu1
        pperp(ia+1,ja+1,ka+1) = pperp(ia+1,ja+1,ka+1) + aaa(8,n)*mu1
      end do

      end if


!2020-11-22s
      if(lpara.ge.2)then

      do k = 1, lpara
!$omp parallel do
      do i = 1, lrzphi
        dns(i,1,1)   = dns(i,1,1)   + wdns(i,1,1,k)
        mom(i,1,1)   = mom(i,1,1)   + wmom(i,1,1,k)
        ppara(i,1,1) = ppara(i,1,1) + wpar(i,1,1,k)
        pperp(i,1,1) = pperp(i,1,1) + wprp(i,1,1,k)
      end do
      end do

      end if
!2020-11-22e


!2016-12-24s

! smoothing
      cwwa = 0.5d0
!      cwwc =-1.d0/6.d0

!2024-12-21s, correction suggested by Panith Adulsiriswad
      call periodic_particle_mlt4b(dns,mom,ppara,pperp)
      call partsm1(dns,cwwa)
      call partsm1(mom,cwwa)

      call partsm1(ppara,cwwa)
      call partsm1(pperp,cwwa)
!2024-12-21e

! calculate density (per volume)

!      do k = 1, lphi
!      do j = 1, lz
!      do i = 1, lr
!$omp parallel do private(vol)
      do i = 1, lrzphi
        vol = 1.0d0/(grr(i,1,1)*dr*dz*dphi)
        dns(i,1,1) = dns(i,1,1)*vol
        mom(i,1,1) = mom(i,1,1)*vol
        ppara(i,1,1) = ppara(i,1,1)*vol/mass
        pperp(i,1,1) = pperp(i,1,1)*vol*babs(i,1,1)
      end do
!      end do
!      end do

!2016-12-24e


! wall effect

      lr1 = lr - 1
      lz1 = lz - 1

      if(my_rank_r.eq.(mpi_proc_r-1))then
!      do k = 1, lphi
!      do j = 1, lz
!$omp parallel do
       do j = 1, lzphi
          dns(lr1,j,1) = dns(lr1,j,1) + dns(lr,j,1)
          mom(lr1,j,1) = mom(lr1,j,1) + mom(lr,j,1)
          dns(lr,j,1) = 0.0d0
          mom(lr,j,1) = 0.0d0

          ppara(lr1,j,1) = ppara(lr1,j,1) + ppara(lr,j,1)
          pperp(lr1,j,1) = pperp(lr1,j,1) + pperp(lr,j,1)
          ppara(lr,j,1) = 0.0d0
          pperp(lr,j,1) = 0.0d0
       end do
!      end do
      end if

      if(my_rank_r.eq.0)then
!      do k = 1, lphi
!      do j = 1, lz
!$omp parallel do
       do j = 1, lzphi
          dns(2,j,1) = dns(2,j,1) + dns(1,j,1)
          mom(2,j,1) = mom(2,j,1) + mom(1,j,1)
          dns(1,j,1) = 0.0d0
          mom(1,j,1) = 0.0d0

          ppara(2,j,1) = ppara(2,j,1) + ppara(1,j,1)
          pperp(2,j,1) = pperp(2,j,1) + pperp(1,j,1)
          ppara(1,j,1) = 0.0d0
          pperp(1,j,1) = 0.0d0
       end do
!      end do
      end if

      if(my_rank_z.eq.(mpi_proc_z-1))then
!$omp parallel do
       do k = 1, lphi
       do i = 1, lr
          dns(i,lz1,k) = dns(i,lz1,k) + dns(i,lz,k)
          mom(i,lz1,k) = mom(i,lz1,k) + mom(i,lz,k)
          dns(i,lz,k) = 0.0d0
          mom(i,lz,k) = 0.0d0

          ppara(i,lz1,k) = ppara(i,lz1,k) + ppara(i,lz,k)
          pperp(i,lz1,k) = pperp(i,lz1,k) + pperp(i,lz,k)
          ppara(i,lz,k) = 0.0d0
          pperp(i,lz,k) = 0.0d0
       end do
       end do
      end if 

      if(my_rank_z.eq.0)then
!$omp parallel do
       do k = 1, lphi
       do i = 1, lr
          dns(i,2,k) = dns(i,2,k) + dns(i,1,k)
          mom(i,2,k) = mom(i,2,k) + mom(i,1,k)
          dns(i,1,k) = 0.0d0
          mom(i,1,k) = 0.0d0

          ppara(i,2,k) = ppara(i,2,k) + ppara(i,1,k)
          pperp(i,2,k) = pperp(i,2,k) + pperp(i,1,k)
          ppara(i,1,k) = 0.0d0
          pperp(i,1,k) = 0.0d0
       end do
       end do
      end if 

!      if(.not.flag_stored)then
! take only n=1 modes
!      n = 1
!      call n1(n,ppara)
!      call n1(n,pperp)

! take modes n_min <= n <=n_max
!      n_min = 1
!      n_max = 2
!      call lowpass2_mlt2(n_min,n_max,dns,mom)
!      call lowpass2_mlt2(n_min,n_max,ppara,pperp)
!      end if


!      do k = 1, lphi
!      do j = 1, lz
!      do i = 1, lr
!$omp parallel do
      do i = 1, lrzphi
        dns(i,1,1) = dns(i,1,1) + dns0(i,1,1)
        mom(i,1,1) = mom(i,1,1) + mom0(i,1,1)
        ppara(i,1,1) = ppara(i,1,1) + ppara0(i,1,1)
        pperp(i,1,1) = pperp(i,1,1) + pperp0(i,1,1)/babs0(i,1,1)*babs(i,1,1)
      end do
!      end do
!      end do

!            call wall_clock(d5)
!            t5 = d5 - d4
!            if(my_rank.eq.0)then
!            write(6,*)'density 5, t5=',t5
!            end if

end subroutine density
!--------------------------------------------------------------------
subroutine moments(marker_num,mass,gc &
                  ,dns,mom,ppara,pperp,qpara,qperp &
                  ,dns0,mom0,ppara0,pperp0,qpara0,qperp0)
! calculate pressure
! modified for Fugaku 2020-11-22
! gc_x(:,n) = 1:r, 2:z, 3:phi, 4:p, 5:mu, 6:weight, 7:active, 8:p_phi, 9:fff, 10:fnrml
!--------------------------------------------------------------------
      use mpiset
      use grid
      use field, only:babs,babs0
      use gyro, only:marker_each_gyro
      implicit none

      integer::marker_num,nvec
      real(8)::gc(ngc2,marker_num)
      real(8)::mass
      real(8)::dns(lr,lz,lphi),mom(lr,lz,lphi)
      real(8)::ppara(lr,lz,lphi),pperp(lr,lz,lphi)
      real(8)::qpara(lr,lz,lphi),qperp(lr,lz,lphi)
      real(8)::dns0(lr,lz,lphi),mom0(lr,lz,lphi)
      real(8)::ppara0(lr,lz,lphi),pperp0(lr,lz,lphi)
      real(8)::qpara0(lr,lz,lphi),qperp0(lr,lz,lphi)
      real(8)::wdns(lr,lz,lphi,lpara),wmom(lr,lz,lphi,lpara)
      real(8)::wpar(lr,lz,lphi,lpara),wprp(lr,lz,lphi,lpara)
      real(8)::wqar(lr,lz,lphi,lpara),wqrp(lr,lz,lphi,lpara)
!2015-09-21s
      integer::ijk_a(3,marker_each_gyro)
      real(8)::aaa(8,marker_each_gyro)
!2015-09-21e
      integer::i,j,k,l,m,n,ia,ia1,ja,ja1,ka,ka1,lr1,lz1
      real(8)::ar,ar1,az,az1,aphi,aphi1
!      real(8)::www,aaa1,aaa2,aaa3,aaa4,aaa5,aaa6,aaa7,aaa8
      real(8)::vol,cwwa,cwwc
!      real(8)::d0,d1,d2,d3,d4,d5,t1,t2,t3,t4,t5
      real(8)::dr1,dz1,dphi1,ma_mi_r
      integer::kr,kz,kphi
      real(8)::p0,p1,p2,mu1,p3,p1mu,mass1 !2015-09-21
      integer::n_min,n_max !2013-05-22
!2020-11-22s
      integer::nvec0
      integer::nmod
      integer,parameter:: NB=320
      real(8):: bbb(NB,8)
      integer::ijk_b(3,NB)
!2020-11-22e

      mass1 = 1.0d0/mass

      dr1 = 1.0d0/dr
      dz1 = 1.0d0/dz
      dphi1 = 1.0d0/dphi

      kr   = 1 - kr_offset(my_rank)
      kz   = 1 - kz_offset(my_rank)
      kphi = 1 + lphi_shd - kphi_offset(my_rank)
!      kphi = 3 - kphi_offset(my_rank)
      ma_mi_r = major_r-minor_r

      if(lpara.ge.2)then
!$omp parallel do
      do k = 1, lpara
        do i = 1, lrzphi
          wdns(i,1,1,k) = 0.0d0
          wmom(i,1,1,k) = 0.0d0
          wpar(i,1,1,k) = 0.0d0
          wprp(i,1,1,k) = 0.0d0
          wqar(i,1,1,k) = 0.0d0
          wqrp(i,1,1,k) = 0.0d0
        end do
      end do

      end if

!2020-11-22s

!$omp parallel do
      do i = 1, lrzphi
        dns(i,1,1) = 0.0d0
        mom(i,1,1) = 0.0d0
        ppara(i,1,1) = 0.0d0
        pperp(i,1,1) = 0.0d0
        qpara(i,1,1) = 0.0d0
        qperp(i,1,1) = 0.0d0
      end do

      nvec = marker_num/lpara
      nvec0 = marker_num/lpara
      nmod  = mod(marker_num,lpara)

      if(lpara.ge.2)then

!$omp parallel private(n,ia,ja,ka,p0,p1,p2,mu1,p3,p1mu,nvec) &
!$omp& private(ar1,ar,az1,az,aphi1,aphi,ijk_b,bbb)
!$omp do
      do l = 1, lpara

        nvec = nvec0 + min(1,nmod/l)
        do m = 1, nvec, NB

        do i=1, min(NB, nvec-m+1)
        n = m+i-1 +  nvec0*(l-1) + min(l-1,nmod)

        ia=max(1,min(lr  -1,int((gc(1,n)-ma_mi_r)*dr1  ) + kr  ))
        ja=max(1,min(lz  -1,int(gc(2,n)          *dz1  ) + kz  ))
        ka=max(1,min(lphi-1,int(gc(3,n)        *dphi1) + kphi))
        ar1 = (gc(1,n)-grr(ia,ja,ka ) ) *dr1
        ar  = 1.0d0 - ar1
        az1 = (gc(2,n)-gzz(ia,ja,ka ) ) *dz1
        az  = 1.0d0 - az1
        aphi1 = (gc(3,n)-gphi(ia,ja,ka ) ) *dphi1
        aphi  = 1.0d0 - aphi1

        bbb(i,1) = ar *az *aphi
        bbb(i,2) = ar1*az *aphi
        bbb(i,3) = ar *az1*aphi
        bbb(i,4) = ar1*az1*aphi
        bbb(i,5) = ar *az *aphi1
        bbb(i,6) = ar1*az *aphi1
        bbb(i,7) = ar *az1*aphi1
        bbb(i,8) = ar1*az1*aphi1

        ijk_b(1,i)=ia
        ijk_b(2,i)=ja
        ijk_b(3,i)=ka
      enddo

      do i=1, min(NB, nvec-m+1)
        ia=ijk_b(1,i)
        ja=ijk_b(2,i)
        ka=ijk_b(3,i)
        n = m+i-1 +  nvec0*(l-1) + min(l-1,nmod)

! gc_x(:,n) = 1:r, 2:z, 3:phi, 4:p, 5:mu, 6:weight, 7:active, 8:p_phi, 9:fff, 10:fnrml
        p0 =  gc(6,n)*gc(7,n)
        p1 =  gc(4,n)*p0
        p2 =  gc(4,n)*p1
        mu1 = gc(5,n)*p0
        p3  = gc(4,n)*p2
        p1mu= gc(4,n)*mu1

        wdns(ia,  ja,  ka,  l) = wdns(ia,  ja,  ka,  l) + bbb(i,1)*p0
        wdns(ia+1,ja,  ka,  l) = wdns(ia+1,ja,  ka,  l) + bbb(i,2)*p0
        wdns(ia,  ja+1,ka,  l) = wdns(ia,  ja+1,ka,  l) + bbb(i,3)*p0
        wdns(ia+1,ja+1,ka,  l) = wdns(ia+1,ja+1,ka,  l) + bbb(i,4)*p0
        wdns(ia,  ja,  ka+1,l) = wdns(ia,  ja,  ka+1,l) + bbb(i,5)*p0
        wdns(ia+1,ja,  ka+1,l) = wdns(ia+1,ja,  ka+1,l) + bbb(i,6)*p0
        wdns(ia,  ja+1,ka+1,l) = wdns(ia,  ja+1,ka+1,l) + bbb(i,7)*p0
        wdns(ia+1,ja+1,ka+1,l) = wdns(ia+1,ja+1,ka+1,l) + bbb(i,8)*p0

        wmom(ia,  ja,  ka,  l) = wmom(ia,  ja,  ka,  l) + bbb(i,1)*p1
        wmom(ia+1,ja,  ka,  l) = wmom(ia+1,ja,  ka,  l) + bbb(i,2)*p1
        wmom(ia,  ja+1,ka,  l) = wmom(ia,  ja+1,ka,  l) + bbb(i,3)*p1
        wmom(ia+1,ja+1,ka,  l) = wmom(ia+1,ja+1,ka,  l) + bbb(i,4)*p1
        wmom(ia,  ja,  ka+1,l) = wmom(ia,  ja,  ka+1,l) + bbb(i,5)*p1
        wmom(ia+1,ja,  ka+1,l) = wmom(ia+1,ja,  ka+1,l) + bbb(i,6)*p1
        wmom(ia,  ja+1,ka+1,l) = wmom(ia,  ja+1,ka+1,l) + bbb(i,7)*p1
        wmom(ia+1,ja+1,ka+1,l) = wmom(ia+1,ja+1,ka+1,l) + bbb(i,8)*p1

        wpar(ia,  ja,  ka,  l) = wpar(ia,  ja,  ka,  l) + bbb(i,1)*p2
        wpar(ia+1,ja,  ka,  l) = wpar(ia+1,ja,  ka,  l) + bbb(i,2)*p2
        wpar(ia,  ja+1,ka,  l) = wpar(ia,  ja+1,ka,  l) + bbb(i,3)*p2
        wpar(ia+1,ja+1,ka,  l) = wpar(ia+1,ja+1,ka,  l) + bbb(i,4)*p2
        wpar(ia,  ja,  ka+1,l) = wpar(ia,  ja,  ka+1,l) + bbb(i,5)*p2
        wpar(ia+1,ja,  ka+1,l) = wpar(ia+1,ja,  ka+1,l) + bbb(i,6)*p2
        wpar(ia,  ja+1,ka+1,l) = wpar(ia,  ja+1,ka+1,l) + bbb(i,7)*p2
        wpar(ia+1,ja+1,ka+1,l) = wpar(ia+1,ja+1,ka+1,l) + bbb(i,8)*p2

        wprp(ia,  ja,  ka,  l) = wprp(ia,  ja,  ka,  l) + bbb(i,1)*mu1
        wprp(ia+1,ja,  ka,  l) = wprp(ia+1,ja,  ka,  l) + bbb(i,2)*mu1
        wprp(ia,  ja+1,ka,  l) = wprp(ia,  ja+1,ka,  l) + bbb(i,3)*mu1
        wprp(ia+1,ja+1,ka,  l) = wprp(ia+1,ja+1,ka,  l) + bbb(i,4)*mu1
        wprp(ia,  ja,  ka+1,l) = wprp(ia,  ja,  ka+1,l) + bbb(i,5)*mu1
        wprp(ia+1,ja,  ka+1,l) = wprp(ia+1,ja,  ka+1,l) + bbb(i,6)*mu1
        wprp(ia,  ja+1,ka+1,l) = wprp(ia,  ja+1,ka+1,l) + bbb(i,7)*mu1
        wprp(ia+1,ja+1,ka+1,l) = wprp(ia+1,ja+1,ka+1,l) + bbb(i,8)*mu1

        wqar(ia,  ja,  ka,  l) = wqar(ia,  ja,  ka,  l) + bbb(i,1)*p3
        wqar(ia+1,ja,  ka,  l) = wqar(ia+1,ja,  ka,  l) + bbb(i,2)*p3
        wqar(ia,  ja+1,ka,  l) = wqar(ia,  ja+1,ka,  l) + bbb(i,3)*p3
        wqar(ia+1,ja+1,ka,  l) = wqar(ia+1,ja+1,ka,  l) + bbb(i,4)*p3
        wqar(ia,  ja,  ka+1,l) = wqar(ia,  ja,  ka+1,l) + bbb(i,5)*p3
        wqar(ia+1,ja,  ka+1,l) = wqar(ia+1,ja,  ka+1,l) + bbb(i,6)*p3
        wqar(ia,  ja+1,ka+1,l) = wqar(ia,  ja+1,ka+1,l) + bbb(i,7)*p3
        wqar(ia+1,ja+1,ka+1,l) = wqar(ia+1,ja+1,ka+1,l) + bbb(i,8)*p3

        wqrp(ia,  ja,  ka,  l) = wqrp(ia,  ja,  ka,  l) + bbb(i,1)*p1mu
        wqrp(ia+1,ja,  ka,  l) = wqrp(ia+1,ja,  ka,  l) + bbb(i,2)*p1mu
        wqrp(ia,  ja+1,ka,  l) = wqrp(ia,  ja+1,ka,  l) + bbb(i,3)*p1mu
        wqrp(ia+1,ja+1,ka,  l) = wqrp(ia+1,ja+1,ka,  l) + bbb(i,4)*p1mu
        wqrp(ia,  ja,  ka+1,l) = wqrp(ia,  ja,  ka+1,l) + bbb(i,5)*p1mu
        wqrp(ia+1,ja,  ka+1,l) = wqrp(ia+1,ja,  ka+1,l) + bbb(i,6)*p1mu
        wqrp(ia,  ja+1,ka+1,l) = wqrp(ia,  ja+1,ka+1,l) + bbb(i,7)*p1mu
        wqrp(ia+1,ja+1,ka+1,l) = wqrp(ia+1,ja+1,ka+1,l) + bbb(i,8)*p1mu
      enddo

      end do
      end do
!$omp end parallel

      else
!2020-11-22e

!poption parallel,schedule(dynamic,1000)
!$omp parallel private(ar1,ar,az1,az,aphi1,aphi)
!$omp do schedule(dynamic,1000)
      do n = 1, marker_num
        ijk_a(1,n)=max(1,min(lr  -1,int((gc(1,n)-ma_mi_r)*dr1  ) + kr  ))
        ijk_a(2,n)=max(1,min(lz  -1,int(gc(2,n)          *dz1  ) + kz  ))
        ijk_a(3,n)=max(1,min(lphi-1,int(gc(3,n)        *dphi1) + kphi))

        ar1 = (gc(1,n)-grr(ijk_a(1,n),ijk_a(2,n),ijk_a(3,n) ) ) *dr1
        ar  = 1.0d0 - ar1
        az1 = (gc(2,n)-gzz(ijk_a(1,n),ijk_a(2,n),ijk_a(3,n) ) ) *dz1
        az  = 1.0d0 - az1
        aphi1 = (gc(3,n)-gphi(ijk_a(1,n),ijk_a(2,n),ijk_a(3,n) ) ) *dphi1
        aphi  = 1.0d0 - aphi1

        aaa(1,n) = ar *az *aphi
        aaa(2,n) = ar1*az *aphi
        aaa(3,n) = ar *az1*aphi
        aaa(4,n) = ar1*az1*aphi
        aaa(5,n) = ar *az *aphi1
        aaa(6,n) = ar1*az *aphi1
        aaa(7,n) = ar *az1*aphi1
        aaa(8,n) = ar1*az1*aphi1
      end do
!$omp end parallel

      do n = 1, marker_num

        ia=ijk_a(1,n)
        ja=ijk_a(2,n)
        ka=ijk_a(3,n)

! gc_x(:,n) = 1:r, 2:z, 3:phi, 4:p, 5:mu, 6:weight, 7:active, 8:p_phi, 9:fff, 10:fnrml
        p0 =  gc(6,n)*gc(7,n)
        p1 =  gc(4,n)*p0
        p2 =  gc(4,n)*p1
        mu1 = gc(5,n)*p0
        p3  = gc(4,n)*p2
        p1mu= gc(4,n)*mu1

        dns(ia,  ja,  ka  ) = dns(ia,  ja,  ka  ) + aaa(1,n)*p0
        dns(ia+1,ja,  ka  ) = dns(ia+1,ja,  ka  ) + aaa(2,n)*p0
        dns(ia,  ja+1,ka  ) = dns(ia,  ja+1,ka  ) + aaa(3,n)*p0
        dns(ia+1,ja+1,ka  ) = dns(ia+1,ja+1,ka  ) + aaa(4,n)*p0
        dns(ia,  ja,  ka+1) = dns(ia,  ja,  ka+1) + aaa(5,n)*p0
        dns(ia+1,ja,  ka+1) = dns(ia+1,ja,  ka+1) + aaa(6,n)*p0
        dns(ia,  ja+1,ka+1) = dns(ia,  ja+1,ka+1) + aaa(7,n)*p0
        dns(ia+1,ja+1,ka+1) = dns(ia+1,ja+1,ka+1) + aaa(8,n)*p0

        mom(ia,  ja,  ka  ) = mom(ia,  ja,  ka  ) + aaa(1,n)*p1
        mom(ia+1,ja,  ka  ) = mom(ia+1,ja,  ka  ) + aaa(2,n)*p1
        mom(ia,  ja+1,ka  ) = mom(ia,  ja+1,ka  ) + aaa(3,n)*p1
        mom(ia+1,ja+1,ka  ) = mom(ia+1,ja+1,ka  ) + aaa(4,n)*p1
        mom(ia,  ja,  ka+1) = mom(ia,  ja,  ka+1) + aaa(5,n)*p1
        mom(ia+1,ja,  ka+1) = mom(ia+1,ja,  ka+1) + aaa(6,n)*p1
        mom(ia,  ja+1,ka+1) = mom(ia,  ja+1,ka+1) + aaa(7,n)*p1
        mom(ia+1,ja+1,ka+1) = mom(ia+1,ja+1,ka+1) + aaa(8,n)*p1

        ppara(ia,  ja,  ka  ) = ppara(ia,  ja,  ka  ) + aaa(1,n)*p2
        ppara(ia+1,ja,  ka  ) = ppara(ia+1,ja,  ka  ) + aaa(2,n)*p2
        ppara(ia,  ja+1,ka  ) = ppara(ia,  ja+1,ka  ) + aaa(3,n)*p2
        ppara(ia+1,ja+1,ka  ) = ppara(ia+1,ja+1,ka  ) + aaa(4,n)*p2
        ppara(ia,  ja,  ka+1) = ppara(ia,  ja,  ka+1) + aaa(5,n)*p2
        ppara(ia+1,ja,  ka+1) = ppara(ia+1,ja,  ka+1) + aaa(6,n)*p2
        ppara(ia,  ja+1,ka+1) = ppara(ia,  ja+1,ka+1) + aaa(7,n)*p2
        ppara(ia+1,ja+1,ka+1) = ppara(ia+1,ja+1,ka+1) + aaa(8,n)*p2

        pperp(ia,  ja,  ka  ) = pperp(ia,  ja,  ka  ) + aaa(1,n)*mu1
        pperp(ia+1,ja,  ka  ) = pperp(ia+1,ja,  ka  ) + aaa(2,n)*mu1
        pperp(ia,  ja+1,ka  ) = pperp(ia,  ja+1,ka  ) + aaa(3,n)*mu1
        pperp(ia+1,ja+1,ka  ) = pperp(ia+1,ja+1,ka  ) + aaa(4,n)*mu1
        pperp(ia,  ja,  ka+1) = pperp(ia,  ja,  ka+1) + aaa(5,n)*mu1
        pperp(ia+1,ja,  ka+1) = pperp(ia+1,ja,  ka+1) + aaa(6,n)*mu1
        pperp(ia,  ja+1,ka+1) = pperp(ia,  ja+1,ka+1) + aaa(7,n)*mu1
        pperp(ia+1,ja+1,ka+1) = pperp(ia+1,ja+1,ka+1) + aaa(8,n)*mu1

        qpara(ia,  ja,  ka  ) = qpara(ia,  ja,  ka  ) + aaa(1,n)*p3
        qpara(ia+1,ja,  ka  ) = qpara(ia+1,ja,  ka  ) + aaa(2,n)*p3
        qpara(ia,  ja+1,ka  ) = qpara(ia,  ja+1,ka  ) + aaa(3,n)*p3
        qpara(ia+1,ja+1,ka  ) = qpara(ia+1,ja+1,ka  ) + aaa(4,n)*p3
        qpara(ia,  ja,  ka+1) = qpara(ia,  ja,  ka+1) + aaa(5,n)*p3
        qpara(ia+1,ja,  ka+1) = qpara(ia+1,ja,  ka+1) + aaa(6,n)*p3
        qpara(ia,  ja+1,ka+1) = qpara(ia,  ja+1,ka+1) + aaa(7,n)*p3
        qpara(ia+1,ja+1,ka+1) = qpara(ia+1,ja+1,ka+1) + aaa(8,n)*p3

        qperp(ia,  ja,  ka  ) = qperp(ia,  ja,  ka  ) + aaa(1,n)*p1mu
        qperp(ia+1,ja,  ka  ) = qperp(ia+1,ja,  ka  ) + aaa(2,n)*p1mu
        qperp(ia,  ja+1,ka  ) = qperp(ia,  ja+1,ka  ) + aaa(3,n)*p1mu
        qperp(ia+1,ja+1,ka  ) = qperp(ia+1,ja+1,ka  ) + aaa(4,n)*p1mu
        qperp(ia,  ja,  ka+1) = qperp(ia,  ja,  ka+1) + aaa(5,n)*p1mu
        qperp(ia+1,ja,  ka+1) = qperp(ia+1,ja,  ka+1) + aaa(6,n)*p1mu
        qperp(ia,  ja+1,ka+1) = qperp(ia,  ja+1,ka+1) + aaa(7,n)*p1mu
        qperp(ia+1,ja+1,ka+1) = qperp(ia+1,ja+1,ka+1) + aaa(8,n)*p1mu

      end do

      end if


!2020-11-22s
      if(lpara.ge.2)then

      do k = 1, lpara
!$omp parallel do
      do i = 1, lrzphi
        dns(i,1,1)   = dns(i,1,1)   + wdns(i,1,1,k)
        mom(i,1,1)   = mom(i,1,1)   + wmom(i,1,1,k)
        ppara(i,1,1) = ppara(i,1,1) + wpar(i,1,1,k)
        pperp(i,1,1) = pperp(i,1,1) + wprp(i,1,1,k)
        qpara(i,1,1) = qpara(i,1,1) + wqar(i,1,1,k)
        qperp(i,1,1) = qperp(i,1,1) + wqrp(i,1,1,k)
      end do
      end do

      end if
!2020-11-22e

! smoothing

       cwwa = 0.5d0
!       cwwc =-1.d0/6.d0 

!2024-12-21s, correction suggested by Panith Adulsiriswad
       call periodic_particle_mlt6b(dns,mom,ppara,pperp,qpara,qperp)
       call partsm1(dns,cwwa)
       call partsm1(mom,cwwa)

       call partsm1(ppara,cwwa)
       call partsm1(pperp,cwwa)

       call partsm1(qpara,cwwa)
       call partsm1(qperp,cwwa)
!2024-12-21e

! calculate density (per volume)

!      do k = 1, lphi
!      do j = 1, lz
!      do i = 1, lr
!$omp parallel do private(vol)
      do i = 1, lrzphi
        vol = 1.0d0/(grr(i,1,1)*dr*dz*dphi)
        dns(i,1,1) = dns(i,1,1)*vol
        mom(i,1,1) = mom(i,1,1)*vol
        ppara(i,1,1) = ppara(i,1,1)*vol*mass1
        pperp(i,1,1) = pperp(i,1,1)*vol*babs(i,1,1)
        qpara(i,1,1) = qpara(i,1,1)*vol*mass1**2
        qperp(i,1,1) = qperp(i,1,1)*vol*babs(i,1,1)*mass1
      end do
!      end do
!      end do

! wall effect

      lr1 = lr - 1
      lz1 = lz - 1

      if(my_rank_r.eq.(mpi_proc_r-1))then
!      do k = 1, lphi
!      do j = 1, lz
!$omp parallel do
       do j = 1, lzphi
          dns(lr1,j,1) = dns(lr1,j,1) + dns(lr,j,1)
          mom(lr1,j,1) = mom(lr1,j,1) + mom(lr,j,1)
          dns(lr,j,1) = 0.0d0
          mom(lr,j,1) = 0.0d0

          ppara(lr1,j,1) = ppara(lr1,j,1) + ppara(lr,j,1)
          pperp(lr1,j,1) = pperp(lr1,j,1) + pperp(lr,j,1)
          ppara(lr,j,1) = 0.0d0
          pperp(lr,j,1) = 0.0d0

          qpara(lr1,j,1) = qpara(lr1,j,1) + qpara(lr,j,1)
          qperp(lr1,j,1) = qperp(lr1,j,1) + qperp(lr,j,1)
          qpara(lr,j,1) = 0.0d0
          qperp(lr,j,1) = 0.0d0
       end do
!      end do
      end if

      if(my_rank_r.eq.0)then
!      do k = 1, lphi
!      do j = 1, lz
!$omp parallel do
       do j = 1, lzphi
          dns(2,j,1) = dns(2,j,1) + dns(1,j,1)
          mom(2,j,1) = mom(2,j,1) + mom(1,j,1)
          dns(1,j,1) = 0.0d0
          mom(1,j,1) = 0.0d0

          ppara(2,j,1) = ppara(2,j,1) + ppara(1,j,1)
          pperp(2,j,1) = pperp(2,j,1) + pperp(1,j,1)
          ppara(1,j,1) = 0.0d0
          pperp(1,j,1) = 0.0d0

          qpara(2,j,1) = qpara(2,j,1) + qpara(1,j,1)
          qperp(2,j,1) = qperp(2,j,1) + qperp(1,j,1)
          qpara(1,j,1) = 0.0d0
          qperp(1,j,1) = 0.0d0
       end do
!      end do
      end if

      if(my_rank_z.eq.(mpi_proc_z-1))then
!$omp parallel do
       do k = 1, lphi
       do i = 1, lr
          dns(i,lz1,k) = dns(i,lz1,k) + dns(i,lz,k)
          mom(i,lz1,k) = mom(i,lz1,k) + mom(i,lz,k)
          dns(i,lz,k) = 0.0d0
          mom(i,lz,k) = 0.0d0

          ppara(i,lz1,k) = ppara(i,lz1,k) + ppara(i,lz,k)
          pperp(i,lz1,k) = pperp(i,lz1,k) + pperp(i,lz,k)
          ppara(i,lz,k) = 0.0d0
          pperp(i,lz,k) = 0.0d0

          qpara(i,lz1,k) = qpara(i,lz1,k) + qpara(i,lz,k)
          qperp(i,lz1,k) = qperp(i,lz1,k) + qperp(i,lz,k)
          qpara(i,lz,k) = 0.0d0
          qperp(i,lz,k) = 0.0d0
       end do
       end do
      end if 

      if(my_rank_z.eq.0)then
!$omp parallel do
       do k = 1, lphi
       do i = 1, lr
          dns(i,2,k) = dns(i,2,k) + dns(i,1,k)
          mom(i,2,k) = mom(i,2,k) + mom(i,1,k)
          dns(i,1,k) = 0.0d0
          mom(i,1,k) = 0.0d0

          ppara(i,2,k) = ppara(i,2,k) + ppara(i,1,k)
          pperp(i,2,k) = pperp(i,2,k) + pperp(i,1,k)
          ppara(i,1,k) = 0.0d0
          pperp(i,1,k) = 0.0d0

          qpara(i,2,k) = qpara(i,2,k) + qpara(i,1,k)
          qperp(i,2,k) = qperp(i,2,k) + qperp(i,1,k)
          qpara(i,1,k) = 0.0d0
          qperp(i,1,k) = 0.0d0
       end do
       end do
      end if 

!      if(.not.flag_stored)then
! take modes n_min <= n <=n_max
!      n_min = 1
!      n_max = 2
!      call lowpass2_mlt2(n_min,n_max,dns,mom)
!      call lowpass2_mlt2(n_min,n_max,ppara,pperp)
!      call lowpass2_mlt2(n_min,n_max,qpara,qperp)
!      end if

!2014-08-28s
!$omp parallel do
      do i = 1, lrzphi
        dns(i,1,1) = dns(i,1,1) + dns0(i,1,1)
        mom(i,1,1) = mom(i,1,1) + mom0(i,1,1)
        ppara(i,1,1) = ppara(i,1,1) + ppara0(i,1,1)
        pperp(i,1,1) = pperp(i,1,1) + pperp0(i,1,1)/babs0(i,1,1)*babs(i,1,1)
        qpara(i,1,1) = qpara(i,1,1) + qpara0(i,1,1)
        qperp(i,1,1) = qperp(i,1,1) + qperp0(i,1,1)/babs0(i,1,1)*babs(i,1,1)
      end do
!2014-08-28e

!            call wall_clock(d5)
!            t5 = d5 - d4
!            if(my_rank.eq.0)then
!            write(6,*)'density 5, t5=',t5
!            end if

end subroutine moments
!--------------------------------------------------------------------
subroutine emf_gyro(marker_num,marker_num_gyro,gc,gyro_phys &
                   ,flp,flp_gyro)
! modified for GK simulation on Fugaku 2022-01-11
! gc_x(:,n) = 1:r, 2:z, 3:phi, 4:p, 5:mu, 6:weight, 7:active, 8:p_phi, 9:fff, 10:fnrml
! flp(nflp, marker_num); nflp=30, nflp_gyro=7
!--------------------------------------------------------------------
      use mpiset
      use grid
      use field, only:babs,babs0,er,ez,ephi,epara,br,bz,bphi,br0,bz0,bphi0,fld
      use gyro, only:ngyro,marker_each_gyro
      implicit none

      integer::marker_num,nvec,marker_num_gyro
      real(8)::gc(ngc2,marker_num)
      real(8)::gyro_phys(2,marker_num_gyro)
      real(8)::flp(nflp,marker_num) !2015-07-08
      real(8)::flp_gyro(nflp_gyro,marker_num_gyro)
      integer::i,j,k,l,m,n,ia,ja,ka,lr1,lz1
      real(8)::ar,ar1,az,az1,aphi,aphi1
!      real(8)::aaa1,aaa2,aaa3,aaa4,aaa5,aaa6,aaa7,aaa8
      real(8)::dr1,dz1,dphi1,ma_mi_r
      integer::kr,kz,kphi
      integer::nc,in
      real(8)::rngr1,rngr1_ac

!2015-07-08s
      integer::ijk_a(3,marker_each_gyro)
      real(8)::aaa(8,marker_each_gyro)
!2015-07-08e

!2020-11-22s
      integer::nvec0
      integer::nmod
      integer,parameter:: NB=320
      real(8):: bbb(NB,8)
      integer::ijk_b(3,NB)
!2020-11-22e

      dr1 = 1.0d0/dr
      dz1 = 1.0d0/dz
      dphi1 = 1.0d0/dphi
      rngr1 = 1.0d0/dble(ngyro)

      kr   = 1 - kr_offset(my_rank)
      kz   = 1 - kz_offset(my_rank)
      kphi = 1 + lphi_shd - kphi_offset(my_rank)
!      kphi = 3 - kphi_offset(my_rank)
      ma_mi_r = major_r-minor_r

!$omp parallel
!$omp workshare

!      flp_gyro = 0.0d0

      fld(1,:,:,:) = er(:,:,:)
      fld(2,:,:,:) = ez(:,:,:)
      fld(3,:,:,:) = ephi(:,:,:)
      fld(4,:,:,:) = epara(:,:,:)
      fld(5,:,:,:) = br(:,:,:) - br0(:,:,:)
      fld(6,:,:,:) = bz(:,:,:) - bz0(:,:,:)
      fld(7,:,:,:) = bphi(:,:,:) - bphi0(:,:,:)
!$omp end workshare
!$omp end parallel


!2020-11-22s
      nvec = marker_num_gyro/lpara
      nvec0 = marker_num_gyro/lpara
      nmod  = mod(marker_num_gyro,lpara)

      if(lpara.ge.2)then

!$omp parallel private(n,nc,ia,ja,ka,nvec) &
!$omp& private(ar1,ar,az1,az,aphi1,aphi,ijk_b,bbb)
!$omp do
      do l = 1, lpara

        nvec = nvec0 + min(1,nmod/l)
        do m = 1, nvec, NB

        do i=1, min(NB, nvec-m+1)
        n = m+i-1 +  nvec0*(l-1) + min(l-1,nmod)
        nc =(n-1)/ngyro + 1

        ia=max(1,min(lr  -1,int((gyro_phys(1,n)-ma_mi_r)*dr1  ) + kr  ))
        ja=max(1,min(lz  -1,int(gyro_phys(2,n)          *dz1  ) + kz  ))
        ka=max(1,min(lphi-1,int(gc(3,nc)        *dphi1) + kphi))
        ar1 = (gyro_phys(1,n)-grr(ia,ja,ka ) ) *dr1
        ar  = 1.0d0 - ar1
        az1 = (gyro_phys(2,n)-gzz(ia,ja,ka ) ) *dz1
        az  = 1.0d0 - az1
        aphi1 = (gc(3,nc)-gphi(ia,ja,ka ) ) *dphi1
        aphi  = 1.0d0 - aphi1

        bbb(i,1) = ar *az *aphi
        bbb(i,2) = ar1*az *aphi
        bbb(i,3) = ar *az1*aphi
        bbb(i,4) = ar1*az1*aphi
        bbb(i,5) = ar *az *aphi1
        bbb(i,6) = ar1*az *aphi1
        bbb(i,7) = ar *az1*aphi1
        bbb(i,8) = ar1*az1*aphi1

        ijk_b(1,i)=ia
        ijk_b(2,i)=ja
        ijk_b(3,i)=ka
        enddo

        do i=1, min(NB, nvec-m+1)
        ia=ijk_b(1,i)
        ja=ijk_b(2,i)
        ka=ijk_b(3,i)
        n = m+i-1 +  nvec0*(l-1) + min(l-1,nmod)

          do in = 1, nflp_gyro
          flp_gyro(in,n) = fld(in, ia, ja,  ka  )*bbb(i,1) + fld(in, ia+1,ja,  ka  )*bbb(i,2) &
                         + fld(in, ia, ja+1,ka  )*bbb(i,3) + fld(in, ia+1,ja+1,ka  )*bbb(i,4) &
                         + fld(in, ia, ja,  ka+1)*bbb(i,5) + fld(in, ia+1,ja,  ka+1)*bbb(i,6) &
                         + fld(in, ia, ja+1,ka+1)*bbb(i,7) + fld(in, ia+1,ja+1,ka+1)*bbb(i,8)
          end do
        end do

     end do
     end do
!$omp end parallel

      else

      do n = 1, marker_num_gyro
        nc =(n-1)/ngyro + 1
        ijk_a(1,n)=max(1,min(lr  -1,int((gyro_phys(1,n)-ma_mi_r)*dr1  ) + kr  ))
        ijk_a(2,n)=max(1,min(lz  -1,int(gyro_phys(2,n)          *dz1  ) + kz  ))
        ijk_a(3,n)=max(1,min(lphi-1,int(gc(3,nc)          *dphi1) + kphi))

        ar1  = max(0.0d0, min(1.0d0, (gyro_phys(1,n) - ma_mi_r)*dr1 -  dble(ijk_a(1,n) - kr)  ) )
        ar   = 1.0d0 - ar1
        az1  = max(0.0d0, min(1.0d0, gyro_phys(2,n)*dz1 - dble(ijk_a(2,n) - kz) ) )
        az   = 1.0d0 - az1
        aphi1 = max(0.0d0, min(1.0d0, gc(3,nc)*dphi1 - dble(ijk_a(3,n) - kphi) ) )
        aphi = 1.0d0 - aphi1

        aaa(1,n) = ar *az *aphi
        aaa(2,n) = ar1*az *aphi
        aaa(3,n) = ar *az1*aphi
        aaa(4,n) = ar1*az1*aphi
        aaa(5,n) = ar *az *aphi1
        aaa(6,n) = ar1*az *aphi1
        aaa(7,n) = ar *az1*aphi1
        aaa(8,n) = ar1*az1*aphi1
      end do

      do n = 1, marker_num_gyro
        ia=ijk_a(1,n)
        ja=ijk_a(2,n)
        ka=ijk_a(3,n)

        do in = 1, nflp_gyro
          flp_gyro(in,n) = fld(in, ia, ja,  ka  )*aaa(1,n) + fld(in, ia+1,ja,  ka  )*aaa(2,n) &
                         + fld(in, ia, ja+1,ka  )*aaa(3,n) + fld(in, ia+1,ja+1,ka  )*aaa(4,n) &
                         + fld(in, ia, ja,  ka+1)*aaa(5,n) + fld(in, ia+1,ja,  ka+1)*aaa(6,n) &
                         + fld(in, ia, ja+1,ka+1)*aaa(7,n) + fld(in, ia+1,ja+1,ka+1)*aaa(8,n)
        end do

      end do

      end if

!$omp parallel do private(rngr1_ac)
      do nc = 1, marker_num
        rngr1_ac = rngr1*gc(7,nc)
        do in = 1, nflp_gyro
          flp(in,nc) =(flp_gyro(in,ngyro*(nc-1)+1) &
                      +flp_gyro(in,ngyro*(nc-1)+2) &
                      +flp_gyro(in,ngyro*(nc-1)+3) &
                      +flp_gyro(in,ngyro*(nc-1)+4) &
                      )*rngr1_ac
        end do
      end do

end subroutine emf_gyro
!--------------------------------------------------------------------
subroutine density_gyro(marker_num,marker_num_gyro,mass,gc,gyro_phys &
                       ,dns,mom,ppara,pperp,dns0,mom0,ppara0,pperp0)
! 2022-01-11, for FLR case
! gc_x(:,n) = 1:r, 2:z, 3:phi, 4:p, 5:mu, 6:weight, 7:active, 8:p_phi, 9:fff, 10:fnrml
!--------------------------------------------------------------------
      use mpiset
      use grid
      use field, only:babs,babs0
      use gyro, only:ngyro,marker_each_gyro
      implicit none

      integer::marker_num,nvec,marker_num_gyro
      real(8)::gc(ngc2,marker_num)
      real(8)::gyro_phys(2,marker_num_gyro)
      real(8)::mass
      real(8)::dns(lr,lz,lphi),mom(lr,lz,lphi)
      real(8)::ppara(lr,lz,lphi),pperp(lr,lz,lphi)
      real(8)::dns0(lr,lz,lphi),mom0(lr,lz,lphi)
      real(8)::ppara0(lr,lz,lphi),pperp0(lr,lz,lphi)
      real(8)::wdns(lr,lz,lphi,lpara),wmom(lr,lz,lphi,lpara)
      real(8)::wpar(lr,lz,lphi,lpara),wprp(lr,lz,lphi,lpara)
!2015-09-21s
      integer::ijk_a(3,marker_each_gyro)
      real(8)::aaa(8,marker_each_gyro)
!2015-09-21e
      integer::i,j,k,l,m,n,ia,ia1,ja,ja1,ka,ka1,lr1,lz1
      real(8)::ar,ar1,az,az1,aphi,aphi1
!      real(8)::www,aaa1,aaa2,aaa3,aaa4,aaa5,aaa6,aaa7,aaa8
      real(8)::vol,cwwa,cwwc
!      real(8)::d0,d1,d2,d3,d4,d5,t1,t2,t3,t4,t5
      real(8)::dr1,dz1,dphi1,ma_mi_r
      integer::kr,kz,kphi
      real(8)::p0,p1,p2,mu1,mass1 !2015-09-21
      integer::n_min,n_max !2013-05-22
      real(8)::rngr1
      integer::nc
!2020-11-22s
      integer::nvec0
      integer::nmod
      integer,parameter:: NB=320
      real(8):: bbb(NB,8)
      integer::ijk_b(3,NB)
!2020-11-22e

      mass1 = 1.0d0/mass

      dr1 = 1.0d0/dr
      dz1 = 1.0d0/dz
      dphi1 = 1.0d0/dphi
      rngr1 = 1.0d0/dble(ngyro)

      kr   = 1 - kr_offset(my_rank)
      kz   = 1 - kz_offset(my_rank)
      kphi = 1 + lphi_shd - kphi_offset(my_rank)
!      kphi = 3 - kphi_offset(my_rank)
      ma_mi_r = major_r-minor_r

      if(lpara.ge.2)then
!$omp parallel do
      do k = 1, lpara
        do i = 1, lrzphi
          wdns(i,1,1,k) = 0.0d0
          wmom(i,1,1,k) = 0.0d0
          wpar(i,1,1,k) = 0.0d0
          wprp(i,1,1,k) = 0.0d0
        end do
      end do

      end if

!2020-11-22s

!$omp parallel do
      do i = 1, lrzphi
        dns(i,1,1) = 0.0d0
        mom(i,1,1) = 0.0d0
        ppara(i,1,1) = 0.0d0
        pperp(i,1,1) = 0.0d0
      end do

      nvec = marker_num_gyro/lpara
      nvec0 = marker_num_gyro/lpara
      nmod  = mod(marker_num_gyro,lpara)

      if(lpara.ge.2)then

!$omp parallel private(n,nc,ia,ja,ka,p0,p1,p2,mu1,nvec) &
!$omp& private(ar1,ar,az1,az,aphi1,aphi,ijk_b,bbb)
!$omp do
      do l = 1, lpara

        nvec = nvec0 + min(1,nmod/l)
        do m = 1, nvec, NB

        do i=1, min(NB, nvec-m+1)
        n = m+i-1 +  nvec0*(l-1) + min(l-1,nmod)
        nc =(n-1)/ngyro + 1

        ia=max(1,min(lr  -1,int((gyro_phys(1,n)-ma_mi_r)*dr1  ) + kr  ))
        ja=max(1,min(lz  -1,int(gyro_phys(2,n)          *dz1  ) + kz  ))
        ka=max(1,min(lphi-1,int(gc(3,nc)        *dphi1) + kphi))
        ar1 = (gyro_phys(1,n)-grr(ia,ja,ka ) ) *dr1
        ar  = 1.0d0 - ar1
        az1 = (gyro_phys(2,n)-gzz(ia,ja,ka ) ) *dz1
        az  = 1.0d0 - az1
        aphi1 = (gc(3,nc)-gphi(ia,ja,ka ) ) *dphi1
        aphi  = 1.0d0 - aphi1

        bbb(i,1) = ar *az *aphi
        bbb(i,2) = ar1*az *aphi
        bbb(i,3) = ar *az1*aphi
        bbb(i,4) = ar1*az1*aphi
        bbb(i,5) = ar *az *aphi1
        bbb(i,6) = ar1*az *aphi1
        bbb(i,7) = ar *az1*aphi1
        bbb(i,8) = ar1*az1*aphi1

        ijk_b(1,i)=ia
        ijk_b(2,i)=ja
        ijk_b(3,i)=ka
      enddo

      do i=1, min(NB, nvec-m+1)
        ia=ijk_b(1,i)
        ja=ijk_b(2,i)
        ka=ijk_b(3,i)
        n = m+i-1 +  nvec0*(l-1) + min(l-1,nmod)
        nc =(n-1)/ngyro + 1

! gc_x(:,n) = 1:r, 2:z, 3:phi, 4:p, 5:mu, 6:weight, 7:active, 8:p_phi, 9:fff, 10:fnrml
        p0 =  gc(6,nc)*gc(7,nc)*rngr1 !rngr1 is multiplied
        p1 =  gc(4,nc)*p0
        p2 =  gc(4,nc)*p1
        mu1 = gc(5,nc)*p0

        wdns(ia,  ja,  ka,  l) = wdns(ia,  ja,  ka,  l) + bbb(i,1)*p0
        wdns(ia+1,ja,  ka,  l) = wdns(ia+1,ja,  ka,  l) + bbb(i,2)*p0
        wdns(ia,  ja+1,ka,  l) = wdns(ia,  ja+1,ka,  l) + bbb(i,3)*p0
        wdns(ia+1,ja+1,ka,  l) = wdns(ia+1,ja+1,ka,  l) + bbb(i,4)*p0
        wdns(ia,  ja,  ka+1,l) = wdns(ia,  ja,  ka+1,l) + bbb(i,5)*p0
        wdns(ia+1,ja,  ka+1,l) = wdns(ia+1,ja,  ka+1,l) + bbb(i,6)*p0
        wdns(ia,  ja+1,ka+1,l) = wdns(ia,  ja+1,ka+1,l) + bbb(i,7)*p0
        wdns(ia+1,ja+1,ka+1,l) = wdns(ia+1,ja+1,ka+1,l) + bbb(i,8)*p0

        wmom(ia,  ja,  ka,  l) = wmom(ia,  ja,  ka,  l) + bbb(i,1)*p1
        wmom(ia+1,ja,  ka,  l) = wmom(ia+1,ja,  ka,  l) + bbb(i,2)*p1
        wmom(ia,  ja+1,ka,  l) = wmom(ia,  ja+1,ka,  l) + bbb(i,3)*p1
        wmom(ia+1,ja+1,ka,  l) = wmom(ia+1,ja+1,ka,  l) + bbb(i,4)*p1
        wmom(ia,  ja,  ka+1,l) = wmom(ia,  ja,  ka+1,l) + bbb(i,5)*p1
        wmom(ia+1,ja,  ka+1,l) = wmom(ia+1,ja,  ka+1,l) + bbb(i,6)*p1
        wmom(ia,  ja+1,ka+1,l) = wmom(ia,  ja+1,ka+1,l) + bbb(i,7)*p1
        wmom(ia+1,ja+1,ka+1,l) = wmom(ia+1,ja+1,ka+1,l) + bbb(i,8)*p1

        wpar(ia,  ja,  ka,  l) = wpar(ia,  ja,  ka,  l) + bbb(i,1)*p2
        wpar(ia+1,ja,  ka,  l) = wpar(ia+1,ja,  ka,  l) + bbb(i,2)*p2
        wpar(ia,  ja+1,ka,  l) = wpar(ia,  ja+1,ka,  l) + bbb(i,3)*p2
        wpar(ia+1,ja+1,ka,  l) = wpar(ia+1,ja+1,ka,  l) + bbb(i,4)*p2
        wpar(ia,  ja,  ka+1,l) = wpar(ia,  ja,  ka+1,l) + bbb(i,5)*p2
        wpar(ia+1,ja,  ka+1,l) = wpar(ia+1,ja,  ka+1,l) + bbb(i,6)*p2
        wpar(ia,  ja+1,ka+1,l) = wpar(ia,  ja+1,ka+1,l) + bbb(i,7)*p2
        wpar(ia+1,ja+1,ka+1,l) = wpar(ia+1,ja+1,ka+1,l) + bbb(i,8)*p2

        wprp(ia,  ja,  ka,  l) = wprp(ia,  ja,  ka,  l) + bbb(i,1)*mu1
        wprp(ia+1,ja,  ka,  l) = wprp(ia+1,ja,  ka,  l) + bbb(i,2)*mu1
        wprp(ia,  ja+1,ka,  l) = wprp(ia,  ja+1,ka,  l) + bbb(i,3)*mu1
        wprp(ia+1,ja+1,ka,  l) = wprp(ia+1,ja+1,ka,  l) + bbb(i,4)*mu1
        wprp(ia,  ja,  ka+1,l) = wprp(ia,  ja,  ka+1,l) + bbb(i,5)*mu1
        wprp(ia+1,ja,  ka+1,l) = wprp(ia+1,ja,  ka+1,l) + bbb(i,6)*mu1
        wprp(ia,  ja+1,ka+1,l) = wprp(ia,  ja+1,ka+1,l) + bbb(i,7)*mu1
        wprp(ia+1,ja+1,ka+1,l) = wprp(ia+1,ja+1,ka+1,l) + bbb(i,8)*mu1
      enddo

      end do
      end do
!$omp end parallel

      else
!2020-11-22e


!$omp parallel private(nc,ar1,ar,az1,az,aphi1,aphi)
!poption parallel,schedule(dynamic,1000)
!$omp do schedule(dynamic,1000)
      do n = 1, marker_num_gyro
        nc =(n-1)/ngyro + 1
        ijk_a(1,n)=max(1,min(lr  -1,int((gyro_phys(1,n)-ma_mi_r)*dr1  ) + kr  ))
        ijk_a(2,n)=max(1,min(lz  -1,int(gyro_phys(2,n)          *dz1  ) + kz  ))
        ijk_a(3,n)=max(1,min(lphi-1,int(gc(3,nc)          *dphi1) + kphi))

        ar1  = max(0.0d0, min(1.0d0, (gyro_phys(1,n) - ma_mi_r)*dr1 -  dble(ijk_a(1,n) - kr)  ) )
        ar   = 1.0d0 - ar1
        az1  = max(0.0d0, min(1.0d0, gyro_phys(2,n)*dz1 - dble(ijk_a(2,n) - kz) ) )
        az   = 1.0d0 - az1
        aphi1 = max(0.0d0, min(1.0d0, gc(3,nc)*dphi1 - dble(ijk_a(3,n) - kphi) ) )
        aphi = 1.0d0 - aphi1

        aaa(1,n) = ar *az *aphi
        aaa(2,n) = ar1*az *aphi
        aaa(3,n) = ar *az1*aphi
        aaa(4,n) = ar1*az1*aphi
        aaa(5,n) = ar *az *aphi1
        aaa(6,n) = ar1*az *aphi1
        aaa(7,n) = ar *az1*aphi1
        aaa(8,n) = ar1*az1*aphi1
      end do
!$omp end parallel

      do n = 1, marker_num_gyro

        nc =(n-1)/ngyro + 1

        ia=ijk_a(1,n)
        ja=ijk_a(2,n)
        ka=ijk_a(3,n)

! gc_x(:,n) = 1:r, 2:z, 3:phi, 4:p, 5:mu, 6:weight, 7:active, 8:p_phi, 9:fff, 10:fnrml
        p0 =  gc(6,nc)*gc(7,nc)*rngr1 !rngr1 is multiplied
        p1 =  gc(4,nc)*p0
        p2 =  gc(4,nc)*p1
        mu1 = gc(5,nc)*p0

        dns(ia,  ja,  ka  ) = dns(ia,  ja,  ka  ) + aaa(1,n)*p0
        dns(ia+1,ja,  ka  ) = dns(ia+1,ja,  ka  ) + aaa(2,n)*p0
        dns(ia,  ja+1,ka  ) = dns(ia,  ja+1,ka  ) + aaa(3,n)*p0
        dns(ia+1,ja+1,ka  ) = dns(ia+1,ja+1,ka  ) + aaa(4,n)*p0
        dns(ia,  ja,  ka+1) = dns(ia,  ja,  ka+1) + aaa(5,n)*p0
        dns(ia+1,ja,  ka+1) = dns(ia+1,ja,  ka+1) + aaa(6,n)*p0
        dns(ia,  ja+1,ka+1) = dns(ia,  ja+1,ka+1) + aaa(7,n)*p0
        dns(ia+1,ja+1,ka+1) = dns(ia+1,ja+1,ka+1) + aaa(8,n)*p0

        mom(ia,  ja,  ka  ) = mom(ia,  ja,  ka  ) + aaa(1,n)*p1
        mom(ia+1,ja,  ka  ) = mom(ia+1,ja,  ka  ) + aaa(2,n)*p1
        mom(ia,  ja+1,ka  ) = mom(ia,  ja+1,ka  ) + aaa(3,n)*p1
        mom(ia+1,ja+1,ka  ) = mom(ia+1,ja+1,ka  ) + aaa(4,n)*p1
        mom(ia,  ja,  ka+1) = mom(ia,  ja,  ka+1) + aaa(5,n)*p1
        mom(ia+1,ja,  ka+1) = mom(ia+1,ja,  ka+1) + aaa(6,n)*p1
        mom(ia,  ja+1,ka+1) = mom(ia,  ja+1,ka+1) + aaa(7,n)*p1
        mom(ia+1,ja+1,ka+1) = mom(ia+1,ja+1,ka+1) + aaa(8,n)*p1

        ppara(ia,  ja,  ka  ) = ppara(ia,  ja,  ka  ) + aaa(1,n)*p2
        ppara(ia+1,ja,  ka  ) = ppara(ia+1,ja,  ka  ) + aaa(2,n)*p2
        ppara(ia,  ja+1,ka  ) = ppara(ia,  ja+1,ka  ) + aaa(3,n)*p2
        ppara(ia+1,ja+1,ka  ) = ppara(ia+1,ja+1,ka  ) + aaa(4,n)*p2
        ppara(ia,  ja,  ka+1) = ppara(ia,  ja,  ka+1) + aaa(5,n)*p2
        ppara(ia+1,ja,  ka+1) = ppara(ia+1,ja,  ka+1) + aaa(6,n)*p2
        ppara(ia,  ja+1,ka+1) = ppara(ia,  ja+1,ka+1) + aaa(7,n)*p2
        ppara(ia+1,ja+1,ka+1) = ppara(ia+1,ja+1,ka+1) + aaa(8,n)*p2

        pperp(ia,  ja,  ka  ) = pperp(ia,  ja,  ka  ) + aaa(1,n)*mu1
        pperp(ia+1,ja,  ka  ) = pperp(ia+1,ja,  ka  ) + aaa(2,n)*mu1
        pperp(ia,  ja+1,ka  ) = pperp(ia,  ja+1,ka  ) + aaa(3,n)*mu1
        pperp(ia+1,ja+1,ka  ) = pperp(ia+1,ja+1,ka  ) + aaa(4,n)*mu1
        pperp(ia,  ja,  ka+1) = pperp(ia,  ja,  ka+1) + aaa(5,n)*mu1
        pperp(ia+1,ja,  ka+1) = pperp(ia+1,ja,  ka+1) + aaa(6,n)*mu1
        pperp(ia,  ja+1,ka+1) = pperp(ia,  ja+1,ka+1) + aaa(7,n)*mu1
        pperp(ia+1,ja+1,ka+1) = pperp(ia+1,ja+1,ka+1) + aaa(8,n)*mu1
      end do

      end if


      if(lpara.ge.2)then

      do k = 1, lpara
!$omp parallel do
      do i = 1, lrzphi
        dns(i,1,1)   = dns(i,1,1)   + wdns(i,1,1,k)
        mom(i,1,1)   = mom(i,1,1)   + wmom(i,1,1,k)
        ppara(i,1,1) = ppara(i,1,1) + wpar(i,1,1,k)
        pperp(i,1,1) = pperp(i,1,1) + wprp(i,1,1,k)
      end do
      end do

      end if


! smoothing

       cwwa = 0.5d0
!       cwwc =-1.d0/6.d0 

!2024-12-21s, correction suggested by Panith Adulsiriswad
       call periodic_particle_mlt4b(dns,mom,ppara,pperp)
       call partsm1(dns,cwwa)
       call partsm1(mom,cwwa)

       call partsm1(ppara,cwwa)
       call partsm1(pperp,cwwa)
!2024-12-21e

! calculate density (per volume)

!      do k = 1, lphi
!      do j = 1, lz
!      do i = 1, lr
!$omp parallel do private(vol)
      do i = 1, lrzphi
        vol = 1.0d0/(grr(i,1,1)*dr*dz*dphi)
        dns(i,1,1) = dns(i,1,1)*vol
        mom(i,1,1) = mom(i,1,1)*vol
        ppara(i,1,1) = ppara(i,1,1)*vol/mass
        pperp(i,1,1) = pperp(i,1,1)*vol*babs(i,1,1)
      end do
!      end do
!      end do


! wall effect

      lr1 = lr - 1
      lz1 = lz - 1

      if(my_rank_r.eq.(mpi_proc_r-1))then
!      do k = 1, lphi
!      do j = 1, lz
!$omp parallel do
       do j = 1, lzphi
          dns(lr1,j,1) = dns(lr1,j,1) + dns(lr,j,1)
          mom(lr1,j,1) = mom(lr1,j,1) + mom(lr,j,1)
          dns(lr,j,1) = 0.0d0
          mom(lr,j,1) = 0.0d0

          ppara(lr1,j,1) = ppara(lr1,j,1) + ppara(lr,j,1)
          pperp(lr1,j,1) = pperp(lr1,j,1) + pperp(lr,j,1)
          ppara(lr,j,1) = 0.0d0
          pperp(lr,j,1) = 0.0d0
       end do
!      end do
      end if

      if(my_rank_r.eq.0)then
!      do k = 1, lphi
!      do j = 1, lz
!$omp parallel do
       do j = 1, lzphi
          dns(2,j,1) = dns(2,j,1) + dns(1,j,1)
          mom(2,j,1) = mom(2,j,1) + mom(1,j,1)
          dns(1,j,1) = 0.0d0
          mom(1,j,1) = 0.0d0

          ppara(2,j,1) = ppara(2,j,1) + ppara(1,j,1)
          pperp(2,j,1) = pperp(2,j,1) + pperp(1,j,1)
          ppara(1,j,1) = 0.0d0
          pperp(1,j,1) = 0.0d0
       end do
!      end do
      end if

      if(my_rank_z.eq.(mpi_proc_z-1))then
!$omp parallel do
       do k = 1, lphi
       do i = 1, lr
          dns(i,lz1,k) = dns(i,lz1,k) + dns(i,lz,k)
          mom(i,lz1,k) = mom(i,lz1,k) + mom(i,lz,k)
          dns(i,lz,k) = 0.0d0
          mom(i,lz,k) = 0.0d0

          ppara(i,lz1,k) = ppara(i,lz1,k) + ppara(i,lz,k)
          pperp(i,lz1,k) = pperp(i,lz1,k) + pperp(i,lz,k)
          ppara(i,lz,k) = 0.0d0
          pperp(i,lz,k) = 0.0d0
       end do
       end do
      end if 

      if(my_rank_z.eq.0)then
!$omp parallel do
       do k = 1, lphi
       do i = 1, lr
          dns(i,2,k) = dns(i,2,k) + dns(i,1,k)
          mom(i,2,k) = mom(i,2,k) + mom(i,1,k)
          dns(i,1,k) = 0.0d0
          mom(i,1,k) = 0.0d0

          ppara(i,2,k) = ppara(i,2,k) + ppara(i,1,k)
          pperp(i,2,k) = pperp(i,2,k) + pperp(i,1,k)
          ppara(i,1,k) = 0.0d0
          pperp(i,1,k) = 0.0d0
       end do
       end do
      end if 

! take only n=1 modes
!      n = 1
!      call n1(n,ppara)
!      call n1(n,pperp)

! take modes n_min <= n <=n_max
!      n_min = 1
!      n_max = 2
!      call lowpass2_mlt2(n_min,n_max,dns,mom)
!      call lowpass2_mlt2(n_min,n_max,ppara,pperp)

!      do k = 1, lphi
!      do j = 1, lz
!      do i = 1, lr
!$omp parallel do
      do i = 1, lrzphi
        dns(i,1,1) = dns(i,1,1) + dns0(i,1,1)
        mom(i,1,1) = mom(i,1,1) + mom0(i,1,1)
        ppara(i,1,1) = ppara(i,1,1) + ppara0(i,1,1)
        pperp(i,1,1) = pperp(i,1,1) + pperp0(i,1,1)/babs0(i,1,1)*babs(i,1,1)
      end do
!      end do
!      end do

end subroutine density_gyro
!--------------------------------------------------------------------
subroutine moments_gyro(marker_num,marker_num_gyro,mass,gc,gyro_phys &
                  ,dns,mom,ppara,pperp,qpara,qperp &
                  ,dns0,mom0,ppara0,pperp0,qpara0,qperp0)
! 2016-02-04, for FLR case
! gc_x(:,n) = 1:r, 2:z, 3:phi, 4:p, 5:mu, 6:weight, 7:active, 8:p_phi, 9:fff, 10:fnrml
!--------------------------------------------------------------------
      use mpiset
      use grid
      use field, only:babs,babs0
      use gyro, only:ngyro,marker_each_gyro
      implicit none

      integer::marker_num,nvec,marker_num_gyro
      real(8)::gc(ngc2,marker_num)
      real(8)::gyro_phys(2,marker_num_gyro)
      real(8)::mass
      real(8)::dns(lr,lz,lphi),mom(lr,lz,lphi)
      real(8)::ppara(lr,lz,lphi),pperp(lr,lz,lphi)
      real(8)::qpara(lr,lz,lphi),qperp(lr,lz,lphi)
      real(8)::dns0(lr,lz,lphi),mom0(lr,lz,lphi)
      real(8)::ppara0(lr,lz,lphi),pperp0(lr,lz,lphi)
      real(8)::qpara0(lr,lz,lphi),qperp0(lr,lz,lphi)
      real(8)::wdns(lr,lz,lphi,lpara),wmom(lr,lz,lphi,lpara)
      real(8)::wpar(lr,lz,lphi,lpara),wprp(lr,lz,lphi,lpara)
      real(8)::wqar(lr,lz,lphi,lpara),wqrp(lr,lz,lphi,lpara)
!2015-09-21s
      integer::ijk_a(3,marker_each_gyro)
      real(8)::aaa(8,marker_each_gyro)
!2015-09-21e
      integer::i,j,k,l,m,n,ia,ia1,ja,ja1,ka,ka1,lr1,lz1
      real(8)::ar,ar1,az,az1,aphi,aphi1
!      real(8)::www,aaa1,aaa2,aaa3,aaa4,aaa5,aaa6,aaa7,aaa8
      real(8)::vol,cwwa,cwwc
!      real(8)::d0,d1,d2,d3,d4,d5,t1,t2,t3,t4,t5
      real(8)::dr1,dz1,dphi1,ma_mi_r
      integer::kr,kz,kphi
      real(8)::p0,p1,p2,mu1,p3,p1mu,mass1 !2015-09-21
      integer::n_min,n_max !2013-05-22
      real(8)::rngr1
      integer::nc
!2020-11-22s
      integer::nvec0
      integer::nmod
      integer,parameter:: NB=320
      real(8):: bbb(NB,8)
      integer::ijk_b(3,NB)
!2020-11-22e

      mass1 = 1.0d0/mass

      dr1 = 1.0d0/dr
      dz1 = 1.0d0/dz
      dphi1 = 1.0d0/dphi
      rngr1 = 1.0d0/dble(ngyro)

      kr   = 1 - kr_offset(my_rank)
      kz   = 1 - kz_offset(my_rank)
      kphi = 1 + lphi_shd - kphi_offset(my_rank)
!      kphi = 3 - kphi_offset(my_rank)
      ma_mi_r = major_r-minor_r

      if(lpara.ge.2)then
!$omp parallel do
      do k = 1, lpara
        do i = 1, lrzphi
          wdns(i,1,1,k) = 0.0d0
          wmom(i,1,1,k) = 0.0d0
          wpar(i,1,1,k) = 0.0d0
          wprp(i,1,1,k) = 0.0d0
          wqar(i,1,1,k) = 0.0d0
          wqrp(i,1,1,k) = 0.0d0
        end do
      end do

      end if

!2020-11-22s

!$omp parallel do
      do i = 1, lrzphi
        dns(i,1,1) = 0.0d0
        mom(i,1,1) = 0.0d0
        ppara(i,1,1) = 0.0d0
        pperp(i,1,1) = 0.0d0
        qpara(i,1,1) = 0.0d0
        qperp(i,1,1) = 0.0d0
      end do

      nvec = marker_num_gyro/lpara
      nvec0 = marker_num_gyro/lpara
      nmod  = mod(marker_num_gyro,lpara)

      if(lpara.ge.2)then

!$omp parallel private(n,nc,ia,ja,ka,p0,p1,p2,mu1,p3,p1mu,nvec) &
!$omp& private(ar1,ar,az1,az,aphi1,aphi,ijk_b,bbb)
!$omp do
      do l = 1, lpara

        nvec = nvec0 + min(1,nmod/l)
        do m = 1, nvec, NB

        do i=1, min(NB, nvec-m+1)
        n = m+i-1 +  nvec0*(l-1) + min(l-1,nmod)
        nc =(n-1)/ngyro + 1

        ia=max(1,min(lr  -1,int((gyro_phys(1,n)-ma_mi_r)*dr1  ) + kr  ))
        ja=max(1,min(lz  -1,int(gyro_phys(2,n)          *dz1  ) + kz  ))
        ka=max(1,min(lphi-1,int(gc(3,nc)        *dphi1) + kphi))
        ar1 = (gyro_phys(1,n)-grr(ia,ja,ka ) ) *dr1
        ar  = 1.0d0 - ar1
        az1 = (gyro_phys(2,n)-gzz(ia,ja,ka ) ) *dz1
        az  = 1.0d0 - az1
        aphi1 = (gc(3,nc)-gphi(ia,ja,ka ) ) *dphi1
        aphi  = 1.0d0 - aphi1

        bbb(i,1) = ar *az *aphi
        bbb(i,2) = ar1*az *aphi
        bbb(i,3) = ar *az1*aphi
        bbb(i,4) = ar1*az1*aphi
        bbb(i,5) = ar *az *aphi1
        bbb(i,6) = ar1*az *aphi1
        bbb(i,7) = ar *az1*aphi1
        bbb(i,8) = ar1*az1*aphi1

        ijk_b(1,i)=ia
        ijk_b(2,i)=ja
        ijk_b(3,i)=ka
      enddo

      do i=1, min(NB, nvec-m+1)
        ia=ijk_b(1,i)
        ja=ijk_b(2,i)
        ka=ijk_b(3,i)
        n = m+i-1 +  nvec0*(l-1) + min(l-1,nmod)
        nc =(n-1)/ngyro + 1

! gc_x(:,n) = 1:r, 2:z, 3:phi, 4:p, 5:mu, 6:weight, 7:active, 8:p_phi, 9:fff, 10:fnrml
        p0 =  gc(6,nc)*gc(7,nc)*rngr1 !rngr1 is multiplied
        p1 =  gc(4,nc)*p0
        p2 =  gc(4,nc)*p1
        mu1 = gc(5,nc)*p0
        p3  = gc(4,nc)*p2
        p1mu= gc(4,nc)*mu1

        wdns(ia,  ja,  ka,  l) = wdns(ia,  ja,  ka,  l) + bbb(i,1)*p0
        wdns(ia+1,ja,  ka,  l) = wdns(ia+1,ja,  ka,  l) + bbb(i,2)*p0
        wdns(ia,  ja+1,ka,  l) = wdns(ia,  ja+1,ka,  l) + bbb(i,3)*p0
        wdns(ia+1,ja+1,ka,  l) = wdns(ia+1,ja+1,ka,  l) + bbb(i,4)*p0
        wdns(ia,  ja,  ka+1,l) = wdns(ia,  ja,  ka+1,l) + bbb(i,5)*p0
        wdns(ia+1,ja,  ka+1,l) = wdns(ia+1,ja,  ka+1,l) + bbb(i,6)*p0
        wdns(ia,  ja+1,ka+1,l) = wdns(ia,  ja+1,ka+1,l) + bbb(i,7)*p0
        wdns(ia+1,ja+1,ka+1,l) = wdns(ia+1,ja+1,ka+1,l) + bbb(i,8)*p0

        wmom(ia,  ja,  ka,  l) = wmom(ia,  ja,  ka,  l) + bbb(i,1)*p1
        wmom(ia+1,ja,  ka,  l) = wmom(ia+1,ja,  ka,  l) + bbb(i,2)*p1
        wmom(ia,  ja+1,ka,  l) = wmom(ia,  ja+1,ka,  l) + bbb(i,3)*p1
        wmom(ia+1,ja+1,ka,  l) = wmom(ia+1,ja+1,ka,  l) + bbb(i,4)*p1
        wmom(ia,  ja,  ka+1,l) = wmom(ia,  ja,  ka+1,l) + bbb(i,5)*p1
        wmom(ia+1,ja,  ka+1,l) = wmom(ia+1,ja,  ka+1,l) + bbb(i,6)*p1
        wmom(ia,  ja+1,ka+1,l) = wmom(ia,  ja+1,ka+1,l) + bbb(i,7)*p1
        wmom(ia+1,ja+1,ka+1,l) = wmom(ia+1,ja+1,ka+1,l) + bbb(i,8)*p1

        wpar(ia,  ja,  ka,  l) = wpar(ia,  ja,  ka,  l) + bbb(i,1)*p2
        wpar(ia+1,ja,  ka,  l) = wpar(ia+1,ja,  ka,  l) + bbb(i,2)*p2
        wpar(ia,  ja+1,ka,  l) = wpar(ia,  ja+1,ka,  l) + bbb(i,3)*p2
        wpar(ia+1,ja+1,ka,  l) = wpar(ia+1,ja+1,ka,  l) + bbb(i,4)*p2
        wpar(ia,  ja,  ka+1,l) = wpar(ia,  ja,  ka+1,l) + bbb(i,5)*p2
        wpar(ia+1,ja,  ka+1,l) = wpar(ia+1,ja,  ka+1,l) + bbb(i,6)*p2
        wpar(ia,  ja+1,ka+1,l) = wpar(ia,  ja+1,ka+1,l) + bbb(i,7)*p2
        wpar(ia+1,ja+1,ka+1,l) = wpar(ia+1,ja+1,ka+1,l) + bbb(i,8)*p2

        wprp(ia,  ja,  ka,  l) = wprp(ia,  ja,  ka,  l) + bbb(i,1)*mu1
        wprp(ia+1,ja,  ka,  l) = wprp(ia+1,ja,  ka,  l) + bbb(i,2)*mu1
        wprp(ia,  ja+1,ka,  l) = wprp(ia,  ja+1,ka,  l) + bbb(i,3)*mu1
        wprp(ia+1,ja+1,ka,  l) = wprp(ia+1,ja+1,ka,  l) + bbb(i,4)*mu1
        wprp(ia,  ja,  ka+1,l) = wprp(ia,  ja,  ka+1,l) + bbb(i,5)*mu1
        wprp(ia+1,ja,  ka+1,l) = wprp(ia+1,ja,  ka+1,l) + bbb(i,6)*mu1
        wprp(ia,  ja+1,ka+1,l) = wprp(ia,  ja+1,ka+1,l) + bbb(i,7)*mu1
        wprp(ia+1,ja+1,ka+1,l) = wprp(ia+1,ja+1,ka+1,l) + bbb(i,8)*mu1

        wqar(ia,  ja,  ka,  l) = wqar(ia,  ja,  ka,  l) + bbb(i,1)*p3
        wqar(ia+1,ja,  ka,  l) = wqar(ia+1,ja,  ka,  l) + bbb(i,2)*p3
        wqar(ia,  ja+1,ka,  l) = wqar(ia,  ja+1,ka,  l) + bbb(i,3)*p3
        wqar(ia+1,ja+1,ka,  l) = wqar(ia+1,ja+1,ka,  l) + bbb(i,4)*p3
        wqar(ia,  ja,  ka+1,l) = wqar(ia,  ja,  ka+1,l) + bbb(i,5)*p3
        wqar(ia+1,ja,  ka+1,l) = wqar(ia+1,ja,  ka+1,l) + bbb(i,6)*p3
        wqar(ia,  ja+1,ka+1,l) = wqar(ia,  ja+1,ka+1,l) + bbb(i,7)*p3
        wqar(ia+1,ja+1,ka+1,l) = wqar(ia+1,ja+1,ka+1,l) + bbb(i,8)*p3

        wqrp(ia,  ja,  ka,  l) = wqrp(ia,  ja,  ka,  l) + bbb(i,1)*p1mu
        wqrp(ia+1,ja,  ka,  l) = wqrp(ia+1,ja,  ka,  l) + bbb(i,2)*p1mu
        wqrp(ia,  ja+1,ka,  l) = wqrp(ia,  ja+1,ka,  l) + bbb(i,3)*p1mu
        wqrp(ia+1,ja+1,ka,  l) = wqrp(ia+1,ja+1,ka,  l) + bbb(i,4)*p1mu
        wqrp(ia,  ja,  ka+1,l) = wqrp(ia,  ja,  ka+1,l) + bbb(i,5)*p1mu
        wqrp(ia+1,ja,  ka+1,l) = wqrp(ia+1,ja,  ka+1,l) + bbb(i,6)*p1mu
        wqrp(ia,  ja+1,ka+1,l) = wqrp(ia,  ja+1,ka+1,l) + bbb(i,7)*p1mu
        wqrp(ia+1,ja+1,ka+1,l) = wqrp(ia+1,ja+1,ka+1,l) + bbb(i,8)*p1mu
      enddo

      end do
      end do
!$omp end parallel

      else
!2020-11-22e


!$omp parallel private(nc,ar1,ar,az1,az,aphi1,aphi)
!poption parallel,schedule(dynamic,1000)
!$omp do schedule(dynamic,1000)
      do n = 1, marker_num_gyro
        nc =(n-1)/ngyro + 1
        ijk_a(1,n)=max(1,min(lr  -1,int((gyro_phys(1,n)-ma_mi_r)*dr1  ) + kr  ))
        ijk_a(2,n)=max(1,min(lz  -1,int(gyro_phys(2,n)          *dz1  ) + kz  ))
        ijk_a(3,n)=max(1,min(lphi-1,int(gc(3,nc)          *dphi1) + kphi))

        ar1  = max(0.0d0, min(1.0d0, (gyro_phys(1,n) - ma_mi_r)*dr1 -  dble(ijk_a(1,n) - kr)  ) )
        ar   = 1.0d0 - ar1
        az1  = max(0.0d0, min(1.0d0, gyro_phys(2,n)*dz1 - dble(ijk_a(2,n) - kz) ) )
        az   = 1.0d0 - az1
        aphi1 = max(0.0d0, min(1.0d0, gc(3,nc)*dphi1 - dble(ijk_a(3,n) - kphi) ) )
        aphi = 1.0d0 - aphi1

        aaa(1,n) = ar *az *aphi
        aaa(2,n) = ar1*az *aphi
        aaa(3,n) = ar *az1*aphi
        aaa(4,n) = ar1*az1*aphi
        aaa(5,n) = ar *az *aphi1
        aaa(6,n) = ar1*az *aphi1
        aaa(7,n) = ar *az1*aphi1
        aaa(8,n) = ar1*az1*aphi1
      end do
!$omp end parallel

      do n = 1, marker_num_gyro

        nc =(n-1)/ngyro + 1

        ia=ijk_a(1,n)
        ja=ijk_a(2,n)
        ka=ijk_a(3,n)

! gc_x(:,n) = 1:r, 2:z, 3:phi, 4:p, 5:mu, 6:weight, 7:active, 8:p_phi, 9:fff, 10:fnrml
        p0 =  gc(6,nc)*gc(7,nc)*rngr1 !rngr1 is multiplied
        p1 =  gc(4,nc)*p0
        p2 =  gc(4,nc)*p1
        mu1 = gc(5,nc)*p0
        p3  = gc(4,nc)*p2
        p1mu= gc(4,nc)*mu1

        dns(ia,  ja,  ka  ) = dns(ia,  ja,  ka  ) + aaa(1,n)*p0
        dns(ia+1,ja,  ka  ) = dns(ia+1,ja,  ka  ) + aaa(2,n)*p0
        dns(ia,  ja+1,ka  ) = dns(ia,  ja+1,ka  ) + aaa(3,n)*p0
        dns(ia+1,ja+1,ka  ) = dns(ia+1,ja+1,ka  ) + aaa(4,n)*p0
        dns(ia,  ja,  ka+1) = dns(ia,  ja,  ka+1) + aaa(5,n)*p0
        dns(ia+1,ja,  ka+1) = dns(ia+1,ja,  ka+1) + aaa(6,n)*p0
        dns(ia,  ja+1,ka+1) = dns(ia,  ja+1,ka+1) + aaa(7,n)*p0
        dns(ia+1,ja+1,ka+1) = dns(ia+1,ja+1,ka+1) + aaa(8,n)*p0

        mom(ia,  ja,  ka  ) = mom(ia,  ja,  ka  ) + aaa(1,n)*p1
        mom(ia+1,ja,  ka  ) = mom(ia+1,ja,  ka  ) + aaa(2,n)*p1
        mom(ia,  ja+1,ka  ) = mom(ia,  ja+1,ka  ) + aaa(3,n)*p1
        mom(ia+1,ja+1,ka  ) = mom(ia+1,ja+1,ka  ) + aaa(4,n)*p1
        mom(ia,  ja,  ka+1) = mom(ia,  ja,  ka+1) + aaa(5,n)*p1
        mom(ia+1,ja,  ka+1) = mom(ia+1,ja,  ka+1) + aaa(6,n)*p1
        mom(ia,  ja+1,ka+1) = mom(ia,  ja+1,ka+1) + aaa(7,n)*p1
        mom(ia+1,ja+1,ka+1) = mom(ia+1,ja+1,ka+1) + aaa(8,n)*p1

        ppara(ia,  ja,  ka  ) = ppara(ia,  ja,  ka  ) + aaa(1,n)*p2
        ppara(ia+1,ja,  ka  ) = ppara(ia+1,ja,  ka  ) + aaa(2,n)*p2
        ppara(ia,  ja+1,ka  ) = ppara(ia,  ja+1,ka  ) + aaa(3,n)*p2
        ppara(ia+1,ja+1,ka  ) = ppara(ia+1,ja+1,ka  ) + aaa(4,n)*p2
        ppara(ia,  ja,  ka+1) = ppara(ia,  ja,  ka+1) + aaa(5,n)*p2
        ppara(ia+1,ja,  ka+1) = ppara(ia+1,ja,  ka+1) + aaa(6,n)*p2
        ppara(ia,  ja+1,ka+1) = ppara(ia,  ja+1,ka+1) + aaa(7,n)*p2
        ppara(ia+1,ja+1,ka+1) = ppara(ia+1,ja+1,ka+1) + aaa(8,n)*p2

        pperp(ia,  ja,  ka  ) = pperp(ia,  ja,  ka  ) + aaa(1,n)*mu1
        pperp(ia+1,ja,  ka  ) = pperp(ia+1,ja,  ka  ) + aaa(2,n)*mu1
        pperp(ia,  ja+1,ka  ) = pperp(ia,  ja+1,ka  ) + aaa(3,n)*mu1
        pperp(ia+1,ja+1,ka  ) = pperp(ia+1,ja+1,ka  ) + aaa(4,n)*mu1
        pperp(ia,  ja,  ka+1) = pperp(ia,  ja,  ka+1) + aaa(5,n)*mu1
        pperp(ia+1,ja,  ka+1) = pperp(ia+1,ja,  ka+1) + aaa(6,n)*mu1
        pperp(ia,  ja+1,ka+1) = pperp(ia,  ja+1,ka+1) + aaa(7,n)*mu1
        pperp(ia+1,ja+1,ka+1) = pperp(ia+1,ja+1,ka+1) + aaa(8,n)*mu1

        qpara(ia,  ja,  ka  ) = qpara(ia,  ja,  ka  ) + aaa(1,n)*p3
        qpara(ia+1,ja,  ka  ) = qpara(ia+1,ja,  ka  ) + aaa(2,n)*p3
        qpara(ia,  ja+1,ka  ) = qpara(ia,  ja+1,ka  ) + aaa(3,n)*p3
        qpara(ia+1,ja+1,ka  ) = qpara(ia+1,ja+1,ka  ) + aaa(4,n)*p3
        qpara(ia,  ja,  ka+1) = qpara(ia,  ja,  ka+1) + aaa(5,n)*p3
        qpara(ia+1,ja,  ka+1) = qpara(ia+1,ja,  ka+1) + aaa(6,n)*p3
        qpara(ia,  ja+1,ka+1) = qpara(ia,  ja+1,ka+1) + aaa(7,n)*p3
        qpara(ia+1,ja+1,ka+1) = qpara(ia+1,ja+1,ka+1) + aaa(8,n)*p3

        qperp(ia,  ja,  ka  ) = qperp(ia,  ja,  ka  ) + aaa(1,n)*p1mu
        qperp(ia+1,ja,  ka  ) = qperp(ia+1,ja,  ka  ) + aaa(2,n)*p1mu
        qperp(ia,  ja+1,ka  ) = qperp(ia,  ja+1,ka  ) + aaa(3,n)*p1mu
        qperp(ia+1,ja+1,ka  ) = qperp(ia+1,ja+1,ka  ) + aaa(4,n)*p1mu
        qperp(ia,  ja,  ka+1) = qperp(ia,  ja,  ka+1) + aaa(5,n)*p1mu
        qperp(ia+1,ja,  ka+1) = qperp(ia+1,ja,  ka+1) + aaa(6,n)*p1mu
        qperp(ia,  ja+1,ka+1) = qperp(ia,  ja+1,ka+1) + aaa(7,n)*p1mu
        qperp(ia+1,ja+1,ka+1) = qperp(ia+1,ja+1,ka+1) + aaa(8,n)*p1mu

      end do

      end if


      if(lpara.ge.2)then

      do k = 1, lpara
!$omp parallel do
      do i = 1, lrzphi
        dns(i,1,1)   = dns(i,1,1)   + wdns(i,1,1,k)
        mom(i,1,1)   = mom(i,1,1)   + wmom(i,1,1,k)
        ppara(i,1,1) = ppara(i,1,1) + wpar(i,1,1,k)
        pperp(i,1,1) = pperp(i,1,1) + wprp(i,1,1,k)
        qpara(i,1,1) = qpara(i,1,1) + wqar(i,1,1,k)
        qperp(i,1,1) = qperp(i,1,1) + wqrp(i,1,1,k)
      end do
      end do

      end if


! smoothing

       cwwa = 0.5d0
!       cwwc =-1.d0/6.d0 

!2024-12-21s, correction suggested by Panith Adulsiriswad
       call periodic_particle_mlt6b(dns,mom,ppara,pperp,qpara,qperp)
       call partsm1(dns,cwwa)
       call partsm1(mom,cwwa)

       call partsm1(ppara,cwwa)
       call partsm1(pperp,cwwa)

       call partsm1(qpara,cwwa)
       call partsm1(qperp,cwwa)
!2024-12-21e

! calculate density (per volume)

!      do k = 1, lphi
!      do j = 1, lz
!      do i = 1, lr
!$omp parallel do private(vol)
      do i = 1, lrzphi
        vol = 1.0d0/(grr(i,1,1)*dr*dz*dphi)
        dns(i,1,1) = dns(i,1,1)*vol
        mom(i,1,1) = mom(i,1,1)*vol
        ppara(i,1,1) = ppara(i,1,1)*vol*mass1
        pperp(i,1,1) = pperp(i,1,1)*vol*babs(i,1,1)
        qpara(i,1,1) = qpara(i,1,1)*vol*mass1**2
        qperp(i,1,1) = qperp(i,1,1)*vol*babs(i,1,1)*mass1
      end do
!      end do
!      end do



! wall effect

      lr1 = lr - 1
      lz1 = lz - 1

      if(my_rank_r.eq.(mpi_proc_r-1))then
!      do k = 1, lphi
!      do j = 1, lz
!$omp parallel do
       do j = 1, lzphi
          dns(lr1,j,1) = dns(lr1,j,1) + dns(lr,j,1)
          mom(lr1,j,1) = mom(lr1,j,1) + mom(lr,j,1)
          dns(lr,j,1) = 0.0d0
          mom(lr,j,1) = 0.0d0

          ppara(lr1,j,1) = ppara(lr1,j,1) + ppara(lr,j,1)
          pperp(lr1,j,1) = pperp(lr1,j,1) + pperp(lr,j,1)
          ppara(lr,j,1) = 0.0d0
          pperp(lr,j,1) = 0.0d0

          qpara(lr1,j,1) = qpara(lr1,j,1) + qpara(lr,j,1)
          qperp(lr1,j,1) = qperp(lr1,j,1) + qperp(lr,j,1)
          qpara(lr,j,1) = 0.0d0
          qperp(lr,j,1) = 0.0d0
       end do
!      end do
      end if

      if(my_rank_r.eq.0)then
!      do k = 1, lphi
!      do j = 1, lz
!$omp parallel do
       do j = 1, lzphi
          dns(2,j,1) = dns(2,j,1) + dns(1,j,1)
          mom(2,j,1) = mom(2,j,1) + mom(1,j,1)
          dns(1,j,1) = 0.0d0
          mom(1,j,1) = 0.0d0

          ppara(2,j,1) = ppara(2,j,1) + ppara(1,j,1)
          pperp(2,j,1) = pperp(2,j,1) + pperp(1,j,1)
          ppara(1,j,1) = 0.0d0
          pperp(1,j,1) = 0.0d0

          qpara(2,j,1) = qpara(2,j,1) + qpara(1,j,1)
          qperp(2,j,1) = qperp(2,j,1) + qperp(1,j,1)
          qpara(1,j,1) = 0.0d0
          qperp(1,j,1) = 0.0d0
       end do
!      end do
      end if

      if(my_rank_z.eq.(mpi_proc_z-1))then
!$omp parallel do
       do k = 1, lphi
       do i = 1, lr
          dns(i,lz1,k) = dns(i,lz1,k) + dns(i,lz,k)
          mom(i,lz1,k) = mom(i,lz1,k) + mom(i,lz,k)
          dns(i,lz,k) = 0.0d0
          mom(i,lz,k) = 0.0d0

          ppara(i,lz1,k) = ppara(i,lz1,k) + ppara(i,lz,k)
          pperp(i,lz1,k) = pperp(i,lz1,k) + pperp(i,lz,k)
          ppara(i,lz,k) = 0.0d0
          pperp(i,lz,k) = 0.0d0

          qpara(i,lz1,k) = qpara(i,lz1,k) + qpara(i,lz,k)
          qperp(i,lz1,k) = qperp(i,lz1,k) + qperp(i,lz,k)
          qpara(i,lz,k) = 0.0d0
          qperp(i,lz,k) = 0.0d0
       end do
       end do
      end if 

      if(my_rank_z.eq.0)then
!$omp parallel do
       do k = 1, lphi
       do i = 1, lr
          dns(i,2,k) = dns(i,2,k) + dns(i,1,k)
          mom(i,2,k) = mom(i,2,k) + mom(i,1,k)
          dns(i,1,k) = 0.0d0
          mom(i,1,k) = 0.0d0

          ppara(i,2,k) = ppara(i,2,k) + ppara(i,1,k)
          pperp(i,2,k) = pperp(i,2,k) + pperp(i,1,k)
          ppara(i,1,k) = 0.0d0
          pperp(i,1,k) = 0.0d0

          qpara(i,2,k) = qpara(i,2,k) + qpara(i,1,k)
          qperp(i,2,k) = qperp(i,2,k) + qperp(i,1,k)
          qpara(i,1,k) = 0.0d0
          qperp(i,1,k) = 0.0d0
       end do
       end do
      end if 

!      if(.not.flag_stored)then
! take modes n_min <= n <=n_max
!      n_min = 1
!      n_max = 2
!      call lowpass2_mlt2(n_min,n_max,dns,mom)
!      call lowpass2_mlt2(n_min,n_max,ppara,pperp)
!      call lowpass2_mlt2(n_min,n_max,qpara,qperp)
!      end if

!2014-08-28s
!$omp parallel do
      do i = 1, lrzphi
        dns(i,1,1) = dns(i,1,1) + dns0(i,1,1)
        mom(i,1,1) = mom(i,1,1) + mom0(i,1,1)
        ppara(i,1,1) = ppara(i,1,1) + ppara0(i,1,1)
        pperp(i,1,1) = pperp(i,1,1) + pperp0(i,1,1)/babs0(i,1,1)*babs(i,1,1)
        qpara(i,1,1) = qpara(i,1,1) + qpara0(i,1,1)
        qperp(i,1,1) = qperp(i,1,1) + qperp0(i,1,1)/babs0(i,1,1)*babs(i,1,1)
      end do
!2014-08-28e

!            call wall_clock(d5)
!            t5 = d5 - d4
!            if(my_rank.eq.0)then
!            write(6,*)'density 5, t5=',t5
!            end if

end subroutine moments_gyro

end module fugaku
!----------------------------------------------------------------------------

!----------------------------------------------------------------------------
module xeon6900p
contains
!--------------------------------------------------------------------
subroutine push(marker_num,amassp,ep &
               ,type,temp,valpha,deltav,clambda0,dclambda & !2012-06-17
               ,gc,dgc,v &
               ,cf_pphi,pphi_min,pphi_max &
               ,flp,ispecies,igyro)
! type=0: maxwellian, type=1: slowing down, type=2: beam
! type=3: beam with finite pitch angle width !2012-06-17
! igyro=0: w/o FLR, igyro=1: w/ FLR
! modified for Subsystem A of Plasma Simulator with Xeon 6900P, 2025-07-11
!   gc_x(:,n) = 1:r, 2:z, 3:phi, 4:p, 5:mu, 6:weight, 7:active, 8:p_phi, 9:fff, 10:fnrml
!--------------------------------------------------------------------
      use mpiset
      use field
      use grid
!      use equi_sol, only:raxis
      use particle, only:nu_krook !2024-04-23
      implicit none

      integer::marker_num,type
      integer::ispecies,igyro
      real(8)::amassp,ep
      real(8)::gc(ngc2,marker_num)
      real(8)::dgc(ngc1,marker_num)
      real(8)::v(marker_num)
      real(8)::cf_pphi(0:lcfpphi),pphi_min,pphi_max !2016-02-04
      real(8)::temp,valpha,deltav
      real(8)::dr1,dz1,dphi1,ep1,amsp1
      integer::n,ia,ia1,ja,ja1,ka,ka1,i
      real(8)::ar,ar1,az,az1,aphi,aphi1
      real(8)::aaa1,aaa2,aaa3,aaa4,aaa5,aaa6,aaa7,aaa8 !2025-07-11

!2015-06-04s
!      integer::ijk_a(3,marker_each)
!      real(8)::aaa(8,marker_each)
!2015-06-04e

      real(8)::flp(nflp,marker_num),detotal
      real(8)::b2,babse,babs0e,b21,bre,bze,bphie
      real(8)::b1,br1a,bz1a,bphi1a,b10,br10,bz10,bphi10,dbpar,dbpr,dbpz,dbpphi
      real(8)::rhopar,orbpr,orbpz,orbpphi,denom1
      real(8)::wvnmlr,wvnmlz,wvnmlphi
      real(8)::dgrdbr,dgrdbz,dgrdbphi,dppar,dppar1,dppar2 !2012-08-31
      real(8)::pvpar,w1nmlr,w1nmlz,w1nmlphi
      real(8)::psinrm,pactive,prof,dprofdpsi,energy0,pphi_n,rminor,bmax
      real(8)::dpphi_n,vpara,vpara0,dvpara,sigma
      real(8)::vlcpart,dvlcpart,vnrm,sd_factor
      real(8)::sqr_pi1
      real(8)::nudt,coef,vrat !2016-01-09
      real(8)::nu_krook_dt !2024-04-23

      integer::kr,kz,kphi
      real(8)::ma_mi_r
      integer::kfl_start,kfl

      real(8)::pt_factor,energy,clambda0,dclambda,clambda(marker_each) !2012-06-17
      real(8)::dwpsi(marker_each),dwenrc(marker_each) !2015-06-23
!      real(8)::psip(marker_each)
!      real(8)::dpsi_dr(marker_each),dpsi_dz(marker_each),dpsi_dphi(marker_each)
      real(8)::bphi1e(marker_each) !2016-01-09

      integer,parameter::nblkp=1024 !NEC
!2021-01-11s
!      integer::ijk_a(nblkp,3)
!      real(8)::aaa(nblkp,8)
!      integer::nn,nsta,nend,ivect
      integer::nn,nsta,nend !2025-07-11
!2021-01-11e


! time derivative of particle position and velocity

      dr1 = 1.0d0/dr
      dz1 = 1.0d0/dz
      dphi1 = 1.0d0/dphi
      ep1 = 1.0d0/ep
      amsp1 = 1.0d0/amassp
      sqr_pi1=1.0d0/sqrt(pi)

      nu_krook_dt = nu_krook * dt !2024-04-23

      kr   = 1 - kr_offset(my_rank)
      kz   = 1 - kz_offset(my_rank)
      kphi = 1 + lphi_shd - kphi_offset(my_rank)
      ma_mi_r = major_r-minor_r

! igyro=0: w/o FLR, igyro/=0: flp(1:7,:) is given in subr. extract_em
      if(igyro.eq.0)then
        kfl_start = 1
      else
        kfl_start = nflp_gyro + 1
      end if

!      call ftrace_region_begin('push0')

!$omp parallel
      if(igyro.eq.0)then !w/o FLR
!$omp workshare
        fld(1,:,:,:) = er(:,:,:)
        fld(2,:,:,:) = ez(:,:,:)
        fld(3,:,:,:) = ephi(:,:,:)
        fld(4,:,:,:) = epara(:,:,:)
        fld(5,:,:,:) = br(:,:,:) - br0(:,:,:)
        fld(6,:,:,:) = bz(:,:,:) - bz0(:,:,:)
        fld(7,:,:,:) = bphi(:,:,:) - bphi0(:,:,:)
!$omp end workshare
      end if

!2025-02-05      
!$omp workshare
      fld(11,:,:,:)= gradbr(:,:,:)
      fld(12,:,:,:)= gradbz(:,:,:)
      fld(13,:,:,:)= gradbphi(:,:,:)
      fld(14,:,:,:)= curvbr(:,:,:)
      fld(15,:,:,:)= curvbz(:,:,:)
      fld(16,:,:,:)= curvbphi(:,:,:)
!$omp end workshare

      if(ispecies.eq.0)then
!$omp workshare
        fld(23,:,:,:)= dns_e0(:,:,:)
        fld(24,:,:,:)= dns_e0_r(:,:,:)
        fld(25,:,:,:)= dns_e0_z(:,:,:)
        fld(26,:,:,:)= dns_e0_phi(:,:,:)
        fld(27,:,:,:)= temp_e0(:,:,:)
        fld(28,:,:,:)= temp_e0_r(:,:,:)
        fld(29,:,:,:)= temp_e0_z(:,:,:)
        fld(30,:,:,:)= temp_e0_phi(:,:,:)
!$omp end workshare
      else if(ispecies.eq.1)then
!$omp workshare
        fld(23,:,:,:)= dns_i0(:,:,:)
        fld(24,:,:,:)= dns_i0_r(:,:,:)
        fld(25,:,:,:)= dns_i0_z(:,:,:)
        fld(26,:,:,:)= dns_i0_phi(:,:,:)
        fld(27,:,:,:)= temp_i0(:,:,:)
        fld(28,:,:,:)= temp_i0_r(:,:,:)
        fld(29,:,:,:)= temp_i0_z(:,:,:)
        fld(30,:,:,:)= temp_i0_phi(:,:,:)
!$omp end workshare
      else if(ispecies.eq.2)then
!$omp workshare
        fld(23,:,:,:)= dns_a0(:,:,:)
        fld(24,:,:,:)= dns_a0_r(:,:,:)
        fld(25,:,:,:)= dns_a0_z(:,:,:)
        fld(26,:,:,:)= dns_a0_phi(:,:,:)
        fld(27,:,:,:)= temp_a0(:,:,:)
        fld(28,:,:,:)= temp_a0_r(:,:,:)
        fld(29,:,:,:)= temp_a0_z(:,:,:)
        fld(30,:,:,:)= temp_a0_phi(:,:,:)
!$omp end workshare
      end if
!$omp end parallel

!      call ftrace_region_end('push0')

! !$omp parallel private(nsta,nend,n,ivect,kfl,ijk_a,ia,ja,ka &
!$omp parallel private(nsta,nend,n,kfl,ia,ja,ka & !2025-07-11
! !$omp& ,aaa,ar1,ar,az1,az,aphi1,aphi &
!$omp& ,ar1,ar,az1,az,aphi1,aphi &
!$omp& ,aaa1,aaa2,aaa3,aaa4,aaa5,aaa6,aaa7,aaa8 & !2025-07-11
!$omp& ,bre,bze,bphie,b2,babse,babs0e,b1,b21,br1a,bz1a,bphi1a &
!$omp& ,b10,br10,bz10,bphi10,dbpar,dbpr,dbpz,dbpphi,rhopar,denom1,orbpr,orbpz,orbpphi &
!$omp& ,wvnmlr,wvnmlz,wvnmlphi,dgrdbr,dgrdbz,dgrdbphi,dppar,dppar2,dppar1,pvpar &
!$omp& ,w1nmlr,w1nmlz,w1nmlphi,detotal,energy &
!$omp& ,rminor,bmax,energy0,sigma,vpara0,vpara,dvpara &
!$omp& ,pphi_n,dpphi_n,prof,dprofdpsi &
!$omp& ,vlcpart,dvlcpart,vnrm,sd_factor,pt_factor,nudt,vrat,coef)
!$omp do schedule(dynamic,1000)
      do nn = 1, marker_num, nblkp
        nsta = nn
        nend = min(nn+nblkp-1,marker_num)

!      call ftrace_region_begin('push1')

!2025-07-11s
      do n = nsta, nend
        ia = max(1,min(lr  -1,int((gc(1,n)-ma_mi_r)*dr1  ) + kr  ))
        ja = max(1,min(lz  -1,int(gc(2,n)          *dz1  ) + kz  ))
        ka = max(1,min(lphi-1,int(gc(3,n)          *dphi1) + kphi))

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

! fields at each particle position
! igyro=0: w/o FLR, igyro/=0: flp(1:7,:) is given in subr. extract_em

!      call ftrace_region_end('push1')
!      call ftrace_region_begin('push2')

        do kfl = kfl_start, nflp
        flp(kfl,n) = fld(kfl, ia, ja,  ka  )*aaa1 + fld(kfl, ia+1,ja,  ka  )*aaa2 &
                   + fld(kfl, ia, ja+1,ka  )*aaa3 + fld(kfl, ia+1,ja+1,ka  )*aaa4 &
                   + fld(kfl, ia, ja,  ka+1)*aaa5 + fld(kfl, ia+1,ja,  ka+1)*aaa6 &
                   + fld(kfl, ia, ja+1,ka+1)*aaa7 + fld(kfl, ia+1,ja+1,ka+1)*aaa8
        end do

      end do
!2025-07-11e

!      call ftrace_region_end('push2')
!      call ftrace_region_begin('push3')

      do n = nsta, nend

! flp(5:7,n): delta_br(z,phi)
        bre = flp(5,n) + flp(8,n)
        bze = flp(6,n) + flp(9,n)
        bphie = flp(7,n) + flp(10,n)

        b2  = bre**2 + bze**2 + bphie**2
        babse= max(eps_b, sqrt(b2) )
        babs0e= max(eps_b, sqrt(flp(8,n)**2 + flp(9,n)**2 + flp(10,n)**2) )
        b1 = 1.0d0/babse
        b21= 1.0d0/b2
        br1a = bre*b1
        bz1a = bze*b1
        bphi1a = bphie*b1
        bphi1e(n) = bphi1a !2016-02-4

        b10 = 1.0d0/babs0e
        br10 = flp(8,n)*b10
        bz10 = flp(9,n)*b10
        bphi10 = flp(10,n)*b10

        dbpar = br1a*br10 + bz1a*bz10 + bphi1a*bphi10
        dbpr  = br1a  - br10*dbpar
        dbpz  = bz1a  - bz10*dbpar
        dbpphi= bphi1a- bphi10*dbpar

! guiding center motion

        rhopar = gc(4,n)*ep1*b1

        denom1 = 1.0d0/(1.0d0 + rhopar*(br1a*flp(14,n) + bz1a*flp(15,n) + bphi1a*flp(16,n)))

        orbpr = (br1a + rhopar*flp(14,n))*denom1
        orbpz = (bz1a + rhopar*flp(15,n))*denom1
        orbpphi = (bphi1a + rhopar*flp(16,n))*denom1

! e x b drift

        wvnmlr = (flp(3,n)*bze-flp(2,n)*bphie)*b21*denom1
        wvnmlz = (flp(1,n)*bphie-flp(3,n)*bre)*b21*denom1
        wvnmlphi = (flp(2,n)*bre-flp(1,n)*bze)*b21*denom1

! grad-b drift

        dgrdbr = gc(5,n)*(bphie*flp(12,n) - bze*flp(13,n))*ep1*b21*denom1
        dgrdbz = gc(5,n)*(bre*flp(13,n) - bphie*flp(11,n))*ep1*b21*denom1
        dgrdbphi = gc(5,n)*(bze*flp(11,n) - bre*flp(12,n))*ep1*b21*denom1

! mirror force

        dppar =-gc(5,n)*( flp(11,n)*orbpr &
                      + flp(12,n)*orbpz &
                      + flp(13,n)*orbpphi &
                      )*dt

!2012-08-31
        dppar2=-gc(5,n)*( flp(11,n)*dbpr*denom1 &
                      + flp(12,n)*dbpz*denom1 &
                      + flp(13,n)*orbpphi &
                      )*dt
!2012-08-31 end

! aceeleration due to electric field and curvature drift

!2012-07-07
        dppar1 = ep*((flp(1,n)*flp(14,n) + flp(2,n)*flp(15,n) + flp(3,n)*flp(16,n) &
                      )*rhopar &
                    + flp(4,n) &
                     )*dt*denom1
!2012-07-07 end

! total drift velocity

        pvpar = gc(4,n)*amsp1
        dgc(1,n) = dt*(pvpar*orbpr + wvnmlr + dgrdbr)*gc(7,n)
        dgc(2,n) = dt*(pvpar*orbpz + wvnmlz + dgrdbz)*gc(7,n)
        dgc(3,n) = dt*(pvpar*orbpphi + wvnmlphi + dgrdbphi)/gc(1,n)*gc(7,n)
        dgc(4,n) =(dppar + dppar1)*gc(7,n)

! temporal evolution of weight of high-energy ion particle

        w1nmlr = wvnmlr + pvpar*dbpr*denom1
        w1nmlz = wvnmlz + pvpar*dbpz*denom1
        w1nmlphi = wvnmlphi + pvpar*dbpphi*denom1

        dgc(8,n) =(ep*(w1nmlr*flp(18,n) + w1nmlz*flp(19,n) + w1nmlphi*flp(20,n))*dt &
                  + bphi1a*(w1nmlr*gc(4,n)*dt + gc(1,n)*(dppar1 + dppar2) ) & !2012-08-31
                  )*gc(7,n)

        detotal =(  ep*(flp(1,n)*dgrdbr + flp(2,n)*dgrdbz + flp(3,n)*dgrdbphi)*dt &
                + dppar1*pvpar )*gc(7,n) &
! the following term considers 'mu* v * grad(B0 - B)', 2025-04-04
                + gc(5,n)*(dgc(1,n)*(flp(31,n)-flp(11,n) ) &
                          +dgc(2,n)*(flp(32,n)-flp(12,n) ) &
                          +dgc(3,n)*(flp(33,n)-flp(13,n) )*gc(1,n) &
                          )
!2025-04-27s
        energy0 = max(1.0d-30, 0.50d0*amassp*pvpar**2 + gc(5,n)*babs0e)
        clambda(n) = gc(5,n)*b0/energy0
        v(n) = sqrt(2.0d0*energy0*amsp1)

!        energy = max(1.0d-30, 0.50d0*amassp*pvpar**2 + gc(5,n)*babse)
!        clambda(n) = gc(5,n)*b0/energy
!        v(n) = sqrt(2.0d0*energy*amsp1)
!2025-04-27e

! weight evolution : weight = f - f0
! d weight /dt = - d f0 /dt
! f0 = f_nrml*prof(psinrm)/(v**3 + flp(21,n)**3)*0.50*erfc((v(n)-valpha)/deltav)

!        rminor = sqrt( (gc(1,n)-raxis)**2 &
!                     + ((gc(2,n)-0.50d0*zleng) )**2 &
!                     )
!        bmax = b0*raxis/(raxis-rminor)
!        energy0 = 0.50d0*amsp1*gc(4,n)**2 + gc(5,n)*babs0e
!        sigma = 0.50d0*(1.0d0 + sign(1.0d0, energy0-gc(5,n)*bmax) ) &
!                      *sign(1.0d0,gc(4,n) )
!        vpara0 = sqrt(2.0d0*(energy0-gc(5,n)*bmin)*amsp1)
!        vpara = vpara0*sigma
!        dvpara = sigma/vpara0*gc(4,n)*dppar1*amsp1**2

!        pphi_n = gc(8,n) - amassp*raxis*vpara
!        dpphi_n= dgc(8,n) - amassp*raxis*dvpara

!        prof = exp(pphi_n/(ep*psimax*0.37d0) )
!        dprofdpsi = prof/(ep*psimax*0.37d0)

!        dwpsi(n) = dpphi_n*dprofdpsi*gc(10,n)
!        dwenrc(n)= detotal/(amassp*v(n) )*prof*gc(10,n)

!2016-08-05s
        dwpsi(n) = ( w1nmlr*flp(24,n) + w1nmlz*flp(25,n) + w1nmlphi*flp(26,n) &
                    +(w1nmlr*flp(28,n) + w1nmlz*flp(29,n) + w1nmlphi*flp(30,n) ) &
                    *0.50d0*(amassp*v(n)**2/flp(27,n) - 3.0d0)*flp(23,n)/flp(27,n) &
                   )*dt*gc(7,n)*gc(10,n)

!        dwpsi(n) = (w1nmlr*flp(24,n) + w1nmlz*flp(25,n) + w1nmlphi*flp(26,n))*dt & !2016-01-09
!                   *gc(7,n)*gc(10,n)
!2016-08-05e

        dwenrc(n)= detotal/(amassp*v(n) )*flp(23,n)*gc(10,n) !2016-01-09

      end do

!      call ftrace_region_end('push3')
!      call ftrace_region_begin('push4')

      if(type.eq.0)then

        do n = nsta, nend
!2016-08-05s
          vlcpart = exp(-0.50d0*amassp*v(n)**2/flp(27,n) )*flp(27,n)**(-1.5d0) !2013-07-17
          dvlcpart = -amassp*v(n)/flp(27,n)*vlcpart 
!          vlcpart = exp(-0.50d0*amassp*v(n)**2/temp)
!          dvlcpart = -amassp*v(n)/temp*vlcpart 
!2016-08-05e
          dwpsi(n) = dwpsi(n)*vlcpart
          dwenrc(n) = dwenrc(n)*dvlcpart
          dgc(6,n) = - dwpsi(n) - dwenrc(n) &
                     - gc(6,n)*nu_krook_dt & !Krook operator at normalized psi > psi_edge, flp(17)=psi
                     *(0.50d0 + sign(0.50d0, 1.0d0 - flp(17,n)/psimax - psi_edge) ) !2024-04-23
        end do

      else if(type.eq.1.or.type.eq.2)then

        do n = nsta, nend
          vnrm =(v(n)-valpha)/deltav
          sd_factor = 1.0d0/(v(n)**3 + flp(21,n)**3)
          vlcpart = 0.50d0*erfc(vnrm)*sd_factor
          dvlcpart= -( 3.0d0*v(n)**2*vlcpart &
                     + exp(-vnrm**2)*sqr_pi1/deltav &
                     )*sd_factor

          dwpsi(n) = dwpsi(n)*vlcpart
          dwenrc(n) = dwenrc(n)*dvlcpart
          dgc(6,n) = - dwpsi(n) - dwenrc(n) &
                     - gc(6,n)*nu_krook_dt & !Krook operator at normalized psi > psi_edge, flp(17)=psi
                     *(0.50d0 + sign(0.50d0, 1.0d0 - flp(17,n)/psimax - psi_edge) ) !2024-04-23
        end do

!2012-06-17
      else if(type.eq.3)then

        do n = nsta, nend
          vnrm =(v(n)-valpha)/deltav
          sd_factor = 1.0d0/(v(n)**3 + flp(21,n)**3)
          pt_factor = exp(-(clambda(n)-clambda0)**2/dclambda**2)

          vlcpart = 0.50d0*erfc(vnrm)*sd_factor*pt_factor
          dvlcpart= -( 3.0d0*v(n)**2*vlcpart &
                     + exp(-vnrm**2)*sqr_pi1/deltav*pt_factor &
                     )*sd_factor &
                  + 4.0d0*vlcpart*clambda(n)*(clambda(n)-clambda0) &
                         /(v(n)*dclambda**2)

          dwpsi(n) = dwpsi(n)*vlcpart
          dwenrc(n) = dwenrc(n)*dvlcpart
          dgc(6,n) = - dwpsi(n) - dwenrc(n) &
                     - gc(6,n)*nu_krook_dt & !Krook operator at normalized psi > psi_edge, flp(17)=psi
                     *(0.50d0 + sign(0.50d0, 1.0d0 - flp(17,n)/psimax - psi_edge) ) !2024-04-23
        end do
!2012-06-17 end

      end if

! slowing down
      if(type.eq.-5)then
        do n = nsta, nend
            dgc(6,n) = 0.0d0 !full-f
            nudt = flp(22,n)*dt*gc(7,n)
            vrat = flp(21,n)/v(n)
            coef =(1.0d0 + vrat**3)*nudt
            dgc(10,n) =-3.0d0*nudt*gc(10,n)
            dgc(5,n)= -2.0d0*gc(5,n)*coef
            dgc(4,n) = dgc(4,n)- gc(4,n)*coef
            dgc(8,n) = dgc(8,n) - gc(4,n)*coef*gc(1,n)*bphi1e(n)
        end do
      end if

!      call ftrace_region_end('push4')

      end do !nn
!$omp end parallel

end subroutine push
!--------------------------------------------------------------------
subroutine density(marker_num,mass,gc &
                  ,dns,mom,ppara,pperp,dns0,mom0,ppara0,pperp0)
! calculate pressure
! modified for Subsystem A of Plasma Simulator with Xeon 6900P, 2025-07-11
! gc_x(:,n) = 1:r, 2:z, 3:phi, 4:p, 5:mu, 6:weight, 7:active, 8:p_phi, 9:fff, 10:fnrml
!--------------------------------------------------------------------
      use mpiset
      use grid
      use field, only:babs,babs0
      implicit none

      integer::marker_num,nvec
      real(8)::gc(ngc2,marker_num)
      real(8)::mass
      real(8)::dns(lr,lz,lphi),mom(lr,lz,lphi)
      real(8)::ppara(lr,lz,lphi),pperp(lr,lz,lphi)
      real(8)::dns0(lr,lz,lphi),mom0(lr,lz,lphi)
      real(8)::ppara0(lr,lz,lphi),pperp0(lr,lz,lphi)
      real(8)::wdns(lr,lz,lphi,lpara),wmom(lr,lz,lphi,lpara)
      real(8)::wpar(lr,lz,lphi,lpara),wprp(lr,lz,lphi,lpara)
!      integer::ijk_a(3,marker_each)
!      real(8)::aaa(8,marker_each)
      integer::i,j,k,l,m,n,ia,ia1,ja,ja1,ka,ka1,lr1,lz1
      real(8)::ar,ar1,az,az1,aphi,aphi1
      real(8)::aaa1,aaa2,aaa3,aaa4,aaa5,aaa6,aaa7,aaa8 !2025-07-11
      real(8)::vol,cwwa,cwwc
      real(8)::d0,d1,d2,d3,d4,d5,t1,t2,t3,t4,t5
      real(8)::dr1,dz1,dphi1,ma_mi_r
      integer::kr,kz,kphi
      real(8)::p0,p1,p2,mu1
      integer::n_min,n_max !2013-05-22
      integer::nvec0
      integer::nmod

      dr1 = 1.0d0/dr
      dz1 = 1.0d0/dz
      dphi1 = 1.0d0/dphi

      kr   = 1 - kr_offset(my_rank)
      kz   = 1 - kz_offset(my_rank)
      kphi = 1 + lphi_shd - kphi_offset(my_rank)
!      kphi = 3 - kphi_offset(my_rank)
      ma_mi_r = major_r-minor_r

      if(lpara.ge.2)then

!$omp parallel do
      do k = 1, lpara
        do i = 1, lrzphi
          wdns(i,1,1,k) = 0.0d0
          wmom(i,1,1,k) = 0.0d0
          wpar(i,1,1,k) = 0.0d0
          wprp(i,1,1,k) = 0.0d0
        end do
      end do

      end if

!$omp parallel do
      do i = 1, lrzphi
        dns(i,1,1) = 0.0d0
        mom(i,1,1) = 0.0d0
        ppara(i,1,1) = 0.0d0
        pperp(i,1,1) = 0.0d0
      end do

!2025-07-11s      
      if(lpara.ge.2)then

        nvec = marker_num/lpara
        nvec0 = marker_num/lpara
        nmod  = mod(marker_num,lpara)

!$omp parallel private(nvec,l,m,n,p0,p1,p2,mu1 &
!$omp& ,ia,ja,ka,ar1,ar,az1,az,aphi1,aphi &
!$omp& ,aaa1,aaa2,aaa3,aaa4,aaa5,aaa6,aaa7,aaa8)
!$omp do
      do l = 1, lpara
        nvec = nvec0 + min(1,nmod/l)
      do m = 1, nvec
        n = m +  nvec0*(l-1) + min(l-1,nmod)

        ia = max(1,min(lr  -1,int((gc(1,n)-ma_mi_r)*dr1  ) + kr  ))
        ja = max(1,min(lz  -1,int(gc(2,n)          *dz1  ) + kz  ))
        ka = max(1,min(lphi-1,int(gc(3,n)          *dphi1) + kphi))

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

! gc_x(:,n) = 1:r, 2:z, 3:phi, 4:p, 5:mu, 6:weight, 7:active, 8:p_phi, 9:fff, 10:fnrml
        p0 =  gc(6,n)*gc(7,n)
        p1 =  gc(4,n)*p0
        p2 =  gc(4,n)*p1
        mu1 = gc(5,n)*p0

        wdns(ia,  ja,  ka,  l) = wdns(ia,  ja,  ka,  l) + aaa1*p0
        wdns(ia+1,ja,  ka,  l) = wdns(ia+1,ja,  ka,  l) + aaa2*p0
        wdns(ia,  ja+1,ka,  l) = wdns(ia,  ja+1,ka,  l) + aaa3*p0
        wdns(ia+1,ja+1,ka,  l) = wdns(ia+1,ja+1,ka,  l) + aaa4*p0
        wdns(ia,  ja,  ka+1,l) = wdns(ia,  ja,  ka+1,l) + aaa5*p0
        wdns(ia+1,ja,  ka+1,l) = wdns(ia+1,ja,  ka+1,l) + aaa6*p0
        wdns(ia,  ja+1,ka+1,l) = wdns(ia,  ja+1,ka+1,l) + aaa7*p0
        wdns(ia+1,ja+1,ka+1,l) = wdns(ia+1,ja+1,ka+1,l) + aaa8*p0

        wmom(ia,  ja,  ka,  l) = wmom(ia,  ja,  ka,  l) + aaa1*p1
        wmom(ia+1,ja,  ka,  l) = wmom(ia+1,ja,  ka,  l) + aaa2*p1
        wmom(ia,  ja+1,ka,  l) = wmom(ia,  ja+1,ka,  l) + aaa3*p1
        wmom(ia+1,ja+1,ka,  l) = wmom(ia+1,ja+1,ka,  l) + aaa4*p1
        wmom(ia,  ja,  ka+1,l) = wmom(ia,  ja,  ka+1,l) + aaa5*p1
        wmom(ia+1,ja,  ka+1,l) = wmom(ia+1,ja,  ka+1,l) + aaa6*p1
        wmom(ia,  ja+1,ka+1,l) = wmom(ia,  ja+1,ka+1,l) + aaa7*p1
        wmom(ia+1,ja+1,ka+1,l) = wmom(ia+1,ja+1,ka+1,l) + aaa8*p1

        wpar(ia,  ja,  ka,  l) = wpar(ia,  ja,  ka,  l) + aaa1*p2
        wpar(ia+1,ja,  ka,  l) = wpar(ia+1,ja,  ka,  l) + aaa2*p2
        wpar(ia,  ja+1,ka,  l) = wpar(ia,  ja+1,ka,  l) + aaa3*p2
        wpar(ia+1,ja+1,ka,  l) = wpar(ia+1,ja+1,ka,  l) + aaa4*p2
        wpar(ia,  ja,  ka+1,l) = wpar(ia,  ja,  ka+1,l) + aaa5*p2
        wpar(ia+1,ja,  ka+1,l) = wpar(ia+1,ja,  ka+1,l) + aaa6*p2
        wpar(ia,  ja+1,ka+1,l) = wpar(ia,  ja+1,ka+1,l) + aaa7*p2
        wpar(ia+1,ja+1,ka+1,l) = wpar(ia+1,ja+1,ka+1,l) + aaa8*p2

        wprp(ia,  ja,  ka,  l) = wprp(ia,  ja,  ka,  l) + aaa1*mu1
        wprp(ia+1,ja,  ka,  l) = wprp(ia+1,ja,  ka,  l) + aaa2*mu1
        wprp(ia,  ja+1,ka,  l) = wprp(ia,  ja+1,ka,  l) + aaa3*mu1
        wprp(ia+1,ja+1,ka,  l) = wprp(ia+1,ja+1,ka,  l) + aaa4*mu1
        wprp(ia,  ja,  ka+1,l) = wprp(ia,  ja,  ka+1,l) + aaa5*mu1
        wprp(ia+1,ja,  ka+1,l) = wprp(ia+1,ja,  ka+1,l) + aaa6*mu1
        wprp(ia,  ja+1,ka+1,l) = wprp(ia,  ja+1,ka+1,l) + aaa7*mu1
        wprp(ia+1,ja+1,ka+1,l) = wprp(ia+1,ja+1,ka+1,l) + aaa8*mu1
      end do
      end do
!$omp end parallel


      else

      do n = 1, marker_num
        ia = max(1,min(lr  -1,int((gc(1,n)-ma_mi_r)*dr1  ) + kr  ))
        ja = max(1,min(lz  -1,int(gc(2,n)          *dz1  ) + kz  ))
        ka = max(1,min(lphi-1,int(gc(3,n)          *dphi1) + kphi))

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

! gc_x(:,n) = 1:r, 2:z, 3:phi, 4:p, 5:mu, 6:weight, 7:active, 8:p_phi, 9:fff, 10:fnrml
        p0 =  gc(6,n)*gc(7,n)
        p1 =  gc(4,n)*p0
        p2 =  gc(4,n)*p1
        mu1 = gc(5,n)*p0

        dns(ia,  ja,  ka  ) = dns(ia,  ja,  ka  ) + aaa1*p0
        dns(ia+1,ja,  ka  ) = dns(ia+1,ja,  ka  ) + aaa2*p0
        dns(ia,  ja+1,ka  ) = dns(ia,  ja+1,ka  ) + aaa3*p0
        dns(ia+1,ja+1,ka  ) = dns(ia+1,ja+1,ka  ) + aaa4*p0
        dns(ia,  ja,  ka+1) = dns(ia,  ja,  ka+1) + aaa5*p0
        dns(ia+1,ja,  ka+1) = dns(ia+1,ja,  ka+1) + aaa6*p0
        dns(ia,  ja+1,ka+1) = dns(ia,  ja+1,ka+1) + aaa7*p0
        dns(ia+1,ja+1,ka+1) = dns(ia+1,ja+1,ka+1) + aaa8*p0

        mom(ia,  ja,  ka  ) = mom(ia,  ja,  ka  ) + aaa1*p1
        mom(ia+1,ja,  ka  ) = mom(ia+1,ja,  ka  ) + aaa2*p1
        mom(ia,  ja+1,ka  ) = mom(ia,  ja+1,ka  ) + aaa3*p1
        mom(ia+1,ja+1,ka  ) = mom(ia+1,ja+1,ka  ) + aaa4*p1
        mom(ia,  ja,  ka+1) = mom(ia,  ja,  ka+1) + aaa5*p1
        mom(ia+1,ja,  ka+1) = mom(ia+1,ja,  ka+1) + aaa6*p1
        mom(ia,  ja+1,ka+1) = mom(ia,  ja+1,ka+1) + aaa7*p1
        mom(ia+1,ja+1,ka+1) = mom(ia+1,ja+1,ka+1) + aaa8*p1

        ppara(ia,  ja,  ka  ) = ppara(ia,  ja,  ka  ) + aaa1*p2
        ppara(ia+1,ja,  ka  ) = ppara(ia+1,ja,  ka  ) + aaa2*p2
        ppara(ia,  ja+1,ka  ) = ppara(ia,  ja+1,ka  ) + aaa3*p2
        ppara(ia+1,ja+1,ka  ) = ppara(ia+1,ja+1,ka  ) + aaa4*p2
        ppara(ia,  ja,  ka+1) = ppara(ia,  ja,  ka+1) + aaa5*p2
        ppara(ia+1,ja,  ka+1) = ppara(ia+1,ja,  ka+1) + aaa6*p2
        ppara(ia,  ja+1,ka+1) = ppara(ia,  ja+1,ka+1) + aaa7*p2
        ppara(ia+1,ja+1,ka+1) = ppara(ia+1,ja+1,ka+1) + aaa8*p2

        pperp(ia,  ja,  ka  ) = pperp(ia,  ja,  ka  ) + aaa1*mu1
        pperp(ia+1,ja,  ka  ) = pperp(ia+1,ja,  ka  ) + aaa2*mu1
        pperp(ia,  ja+1,ka  ) = pperp(ia,  ja+1,ka  ) + aaa3*mu1
        pperp(ia+1,ja+1,ka  ) = pperp(ia+1,ja+1,ka  ) + aaa4*mu1
        pperp(ia,  ja,  ka+1) = pperp(ia,  ja,  ka+1) + aaa5*mu1
        pperp(ia+1,ja,  ka+1) = pperp(ia+1,ja,  ka+1) + aaa6*mu1
        pperp(ia,  ja+1,ka+1) = pperp(ia,  ja+1,ka+1) + aaa7*mu1
        pperp(ia+1,ja+1,ka+1) = pperp(ia+1,ja+1,ka+1) + aaa8*mu1
      end do

      end if
!2025-07-11e

      if(lpara.ge.2)then

      do k = 1, lpara
!$omp parallel do
      do i = 1, lrzphi
        dns(i,1,1)   = dns(i,1,1)   + wdns(i,1,1,k)
        mom(i,1,1)   = mom(i,1,1)   + wmom(i,1,1,k)
        ppara(i,1,1) = ppara(i,1,1) + wpar(i,1,1,k)
        pperp(i,1,1) = pperp(i,1,1) + wprp(i,1,1,k)
      end do
      end do

      end if

!2016-12-24s

! smoothing
      cwwa = 0.5d0
!      cwwc =-1.d0/6.d0 

!2024-12-21s, correction suggested by Panith Adulsiriswad
      call periodic_particle_mlt4b(dns,mom,ppara,pperp)
      call partsm1(dns,cwwa)
      call partsm1(mom,cwwa)

      call partsm1(ppara,cwwa)
      call partsm1(pperp,cwwa)
!2024-12-21e

! calculate density (per volume)

!      do k = 1, lphi
!      do j = 1, lz
!      do i = 1, lr
!$omp parallel do private(vol)
      do i = 1, lrzphi
        vol = 1.0d0/(grr(i,1,1)*dr*dz*dphi)
        dns(i,1,1) = dns(i,1,1)*vol
        mom(i,1,1) = mom(i,1,1)*vol
        ppara(i,1,1) = ppara(i,1,1)*vol/mass
        pperp(i,1,1) = pperp(i,1,1)*vol*babs(i,1,1)
      end do
!      end do
!      end do

!2016-12-24e


! wall effect

      lr1 = lr - 1
      lz1 = lz - 1

      if(my_rank_r.eq.(mpi_proc_r-1))then
!      do k = 1, lphi
!      do j = 1, lz
!$omp parallel do
       do j = 1, lzphi
          dns(lr1,j,1) = dns(lr1,j,1) + dns(lr,j,1)
          mom(lr1,j,1) = mom(lr1,j,1) + mom(lr,j,1)
          dns(lr,j,1) = 0.0d0
          mom(lr,j,1) = 0.0d0

          ppara(lr1,j,1) = ppara(lr1,j,1) + ppara(lr,j,1)
          pperp(lr1,j,1) = pperp(lr1,j,1) + pperp(lr,j,1)
          ppara(lr,j,1) = 0.0d0
          pperp(lr,j,1) = 0.0d0
       end do
!      end do
      end if

      if(my_rank_r.eq.0)then
!      do k = 1, lphi
!      do j = 1, lz
!$omp parallel do
       do j = 1, lzphi
          dns(2,j,1) = dns(2,j,1) + dns(1,j,1)
          mom(2,j,1) = mom(2,j,1) + mom(1,j,1)
          dns(1,j,1) = 0.0d0
          mom(1,j,1) = 0.0d0

          ppara(2,j,1) = ppara(2,j,1) + ppara(1,j,1)
          pperp(2,j,1) = pperp(2,j,1) + pperp(1,j,1)
          ppara(1,j,1) = 0.0d0
          pperp(1,j,1) = 0.0d0
       end do
!      end do
      end if

      if(my_rank_z.eq.(mpi_proc_z-1))then
!$omp parallel do
       do k = 1, lphi
       do i = 1, lr
          dns(i,lz1,k) = dns(i,lz1,k) + dns(i,lz,k)
          mom(i,lz1,k) = mom(i,lz1,k) + mom(i,lz,k)
          dns(i,lz,k) = 0.0d0
          mom(i,lz,k) = 0.0d0

          ppara(i,lz1,k) = ppara(i,lz1,k) + ppara(i,lz,k)
          pperp(i,lz1,k) = pperp(i,lz1,k) + pperp(i,lz,k)
          ppara(i,lz,k) = 0.0d0
          pperp(i,lz,k) = 0.0d0
       end do
       end do
      end if 

      if(my_rank_z.eq.0)then
!$omp parallel do
       do k = 1, lphi
       do i = 1, lr
          dns(i,2,k) = dns(i,2,k) + dns(i,1,k)
          mom(i,2,k) = mom(i,2,k) + mom(i,1,k)
          dns(i,1,k) = 0.0d0
          mom(i,1,k) = 0.0d0

          ppara(i,2,k) = ppara(i,2,k) + ppara(i,1,k)
          pperp(i,2,k) = pperp(i,2,k) + pperp(i,1,k)
          ppara(i,1,k) = 0.0d0
          pperp(i,1,k) = 0.0d0
       end do
       end do
      end if 

!      if(.not.flag_stored)then
! take only n=1 modes
!      n = 1
!      call n1(n,ppara)
!      call n1(n,pperp)

! take modes n_min <= n <=n_max
!      n_min = 1
!      n_max = 2
!      call lowpass2_mlt2(n_min,n_max,dns,mom)
!      call lowpass2_mlt2(n_min,n_max,ppara,pperp)
!      end if


!      do k = 1, lphi
!      do j = 1, lz
!      do i = 1, lr
!$omp parallel do
      do i = 1, lrzphi
        dns(i,1,1) = dns(i,1,1) + dns0(i,1,1)
        mom(i,1,1) = mom(i,1,1) + mom0(i,1,1)
        ppara(i,1,1) = ppara(i,1,1) + ppara0(i,1,1)
        pperp(i,1,1) = pperp(i,1,1) + pperp0(i,1,1)/babs0(i,1,1)*babs(i,1,1)
      end do
!      end do
!      end do

!            call wall_clock(d5)
!            t5 = d5 - d4
!            if(my_rank.eq.0)then
!            write(6,*)'density 5, t5=',t5
!            end if

end subroutine density
!--------------------------------------------------------------------
subroutine moments(marker_num,mass,gc &
                  ,dns,mom,ppara,pperp,qpara,qperp &
                  ,dns0,mom0,ppara0,pperp0,qpara0,qperp0)
! calculate pressure
! modified for Subsystem A of Plasma Simulator with Xeon 6900P, 2025-07-11
! 2016-01-09, modified similar to subroutine density
! gc_x(:,n) = 1:r, 2:z, 3:phi, 4:p, 5:mu, 6:weight, 7:active, 8:p_phi, 9:fff, 10:fnrml
!--------------------------------------------------------------------
      use mpiset
      use grid
      use field, only:babs,babs0
      use gyro, only:marker_each_gyro
      implicit none

      integer::marker_num,nvec
      real(8)::gc(ngc2,marker_num)
      real(8)::mass
      real(8)::dns(lr,lz,lphi),mom(lr,lz,lphi)
      real(8)::ppara(lr,lz,lphi),pperp(lr,lz,lphi)
      real(8)::qpara(lr,lz,lphi),qperp(lr,lz,lphi)
      real(8)::dns0(lr,lz,lphi),mom0(lr,lz,lphi)
      real(8)::ppara0(lr,lz,lphi),pperp0(lr,lz,lphi)
      real(8)::qpara0(lr,lz,lphi),qperp0(lr,lz,lphi)
      real(8)::wdns(lr,lz,lphi,lpara),wmom(lr,lz,lphi,lpara)
      real(8)::wpar(lr,lz,lphi,lpara),wprp(lr,lz,lphi,lpara)
      real(8)::wqar(lr,lz,lphi,lpara),wqrp(lr,lz,lphi,lpara)
!2015-09-21s
!      integer::ijk_a(3,marker_each_gyro)
!      real(8)::aaa(8,marker_each_gyro)
!2015-09-21e
      integer::i,j,k,l,m,n,ia,ia1,ja,ja1,ka,ka1,lr1,lz1
      real(8)::ar,ar1,az,az1,aphi,aphi1
      real(8)::aaa1,aaa2,aaa3,aaa4,aaa5,aaa6,aaa7,aaa8 !2025-07-11
      real(8)::vol,cwwa,cwwc
!      real(8)::d0,d1,d2,d3,d4,d5,t1,t2,t3,t4,t5
      real(8)::dr1,dz1,dphi1,ma_mi_r
      integer::kr,kz,kphi
      real(8)::p0,p1,p2,mu1,p3,p1mu,mass1 !2015-09-21
      integer::n_min,n_max !2013-05-22
      integer::nvec0
      integer::nmod

      mass1 = 1.0d0/mass

      dr1 = 1.0d0/dr
      dz1 = 1.0d0/dz
      dphi1 = 1.0d0/dphi

      kr   = 1 - kr_offset(my_rank)
      kz   = 1 - kz_offset(my_rank)
      kphi = 1 + lphi_shd - kphi_offset(my_rank)
!      kphi = 3 - kphi_offset(my_rank)
      ma_mi_r = major_r-minor_r

      if(lpara.ge.2)then
!$omp parallel do
      do k = 1, lpara
        do i = 1, lrzphi
          wdns(i,1,1,k) = 0.0d0
          wmom(i,1,1,k) = 0.0d0
          wpar(i,1,1,k) = 0.0d0
          wprp(i,1,1,k) = 0.0d0
          wqar(i,1,1,k) = 0.0d0
          wqrp(i,1,1,k) = 0.0d0
        end do
      end do

      end if

      do i = 1, lrzphi
        dns(i,1,1) = 0.0d0
        mom(i,1,1) = 0.0d0
        ppara(i,1,1) = 0.0d0
        pperp(i,1,1) = 0.0d0
        qpara(i,1,1) = 0.0d0
        qperp(i,1,1) = 0.0d0
      end do

!2025-07-11s      
      if(lpara.ge.2)then

        nvec = marker_num/lpara
        nvec0 = marker_num/lpara
        nmod  = mod(marker_num,lpara)

!$omp parallel private(nvec,l,m,n,p0,p1,p2,mu1,p3,p1mu &
!$omp& ,ia,ja,ka,ar1,ar,az1,az,aphi1,aphi &
!$omp& ,aaa1,aaa2,aaa3,aaa4,aaa5,aaa6,aaa7,aaa8)
!$omp do
      do l = 1, lpara
        nvec = nvec0 + min(1,nmod/l)
      do m = 1, nvec
        n = m +  nvec0*(l-1) + min(l-1,nmod)

        ia = max(1,min(lr  -1,int((gc(1,n)-ma_mi_r)*dr1  ) + kr  ))
        ja = max(1,min(lz  -1,int(gc(2,n)          *dz1  ) + kz  ))
        ka = max(1,min(lphi-1,int(gc(3,n)          *dphi1) + kphi))

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

! gc_x(:,n) = 1:r, 2:z, 3:phi, 4:p, 5:mu, 6:weight, 7:active, 8:p_phi, 9:fff, 10:fnrml
        p0 =  gc(6,n)*gc(7,n)
        p1 =  gc(4,n)*p0
        p2 =  gc(4,n)*p1
        mu1 = gc(5,n)*p0
        p3  = gc(4,n)*p2
        p1mu= gc(4,n)*mu1

        wdns(ia,  ja,  ka,  l) = wdns(ia,  ja,  ka,  l) + aaa1*p0
        wdns(ia+1,ja,  ka,  l) = wdns(ia+1,ja,  ka,  l) + aaa2*p0
        wdns(ia,  ja+1,ka,  l) = wdns(ia,  ja+1,ka,  l) + aaa3*p0
        wdns(ia+1,ja+1,ka,  l) = wdns(ia+1,ja+1,ka,  l) + aaa4*p0
        wdns(ia,  ja,  ka+1,l) = wdns(ia,  ja,  ka+1,l) + aaa5*p0
        wdns(ia+1,ja,  ka+1,l) = wdns(ia+1,ja,  ka+1,l) + aaa6*p0
        wdns(ia,  ja+1,ka+1,l) = wdns(ia,  ja+1,ka+1,l) + aaa7*p0
        wdns(ia+1,ja+1,ka+1,l) = wdns(ia+1,ja+1,ka+1,l) + aaa8*p0

        wmom(ia,  ja,  ka,  l) = wmom(ia,  ja,  ka,  l) + aaa1*p1
        wmom(ia+1,ja,  ka,  l) = wmom(ia+1,ja,  ka,  l) + aaa2*p1
        wmom(ia,  ja+1,ka,  l) = wmom(ia,  ja+1,ka,  l) + aaa3*p1
        wmom(ia+1,ja+1,ka,  l) = wmom(ia+1,ja+1,ka,  l) + aaa4*p1
        wmom(ia,  ja,  ka+1,l) = wmom(ia,  ja,  ka+1,l) + aaa5*p1
        wmom(ia+1,ja,  ka+1,l) = wmom(ia+1,ja,  ka+1,l) + aaa6*p1
        wmom(ia,  ja+1,ka+1,l) = wmom(ia,  ja+1,ka+1,l) + aaa7*p1
        wmom(ia+1,ja+1,ka+1,l) = wmom(ia+1,ja+1,ka+1,l) + aaa8*p1

        wpar(ia,  ja,  ka,  l) = wpar(ia,  ja,  ka,  l) + aaa1*p2
        wpar(ia+1,ja,  ka,  l) = wpar(ia+1,ja,  ka,  l) + aaa2*p2
        wpar(ia,  ja+1,ka,  l) = wpar(ia,  ja+1,ka,  l) + aaa3*p2
        wpar(ia+1,ja+1,ka,  l) = wpar(ia+1,ja+1,ka,  l) + aaa4*p2
        wpar(ia,  ja,  ka+1,l) = wpar(ia,  ja,  ka+1,l) + aaa5*p2
        wpar(ia+1,ja,  ka+1,l) = wpar(ia+1,ja,  ka+1,l) + aaa6*p2
        wpar(ia,  ja+1,ka+1,l) = wpar(ia,  ja+1,ka+1,l) + aaa7*p2
        wpar(ia+1,ja+1,ka+1,l) = wpar(ia+1,ja+1,ka+1,l) + aaa8*p2

        wprp(ia,  ja,  ka,  l) = wprp(ia,  ja,  ka,  l) + aaa1*mu1
        wprp(ia+1,ja,  ka,  l) = wprp(ia+1,ja,  ka,  l) + aaa2*mu1
        wprp(ia,  ja+1,ka,  l) = wprp(ia,  ja+1,ka,  l) + aaa3*mu1
        wprp(ia+1,ja+1,ka,  l) = wprp(ia+1,ja+1,ka,  l) + aaa4*mu1
        wprp(ia,  ja,  ka+1,l) = wprp(ia,  ja,  ka+1,l) + aaa5*mu1
        wprp(ia+1,ja,  ka+1,l) = wprp(ia+1,ja,  ka+1,l) + aaa6*mu1
        wprp(ia,  ja+1,ka+1,l) = wprp(ia,  ja+1,ka+1,l) + aaa7*mu1
        wprp(ia+1,ja+1,ka+1,l) = wprp(ia+1,ja+1,ka+1,l) + aaa8*mu1

        wqar(ia,  ja,  ka,  l) = wqar(ia,  ja,  ka,  l) + aaa1*p3
        wqar(ia+1,ja,  ka,  l) = wqar(ia+1,ja,  ka,  l) + aaa2*p3
        wqar(ia,  ja+1,ka,  l) = wqar(ia,  ja+1,ka,  l) + aaa3*p3
        wqar(ia+1,ja+1,ka,  l) = wqar(ia+1,ja+1,ka,  l) + aaa4*p3
        wqar(ia,  ja,  ka+1,l) = wqar(ia,  ja,  ka+1,l) + aaa5*p3
        wqar(ia+1,ja,  ka+1,l) = wqar(ia+1,ja,  ka+1,l) + aaa6*p3
        wqar(ia,  ja+1,ka+1,l) = wqar(ia,  ja+1,ka+1,l) + aaa7*p3
        wqar(ia+1,ja+1,ka+1,l) = wqar(ia+1,ja+1,ka+1,l) + aaa8*p3

        wqrp(ia,  ja,  ka,  l) = wqrp(ia,  ja,  ka,  l) + aaa1*p1mu
        wqrp(ia+1,ja,  ka,  l) = wqrp(ia+1,ja,  ka,  l) + aaa2*p1mu
        wqrp(ia,  ja+1,ka,  l) = wqrp(ia,  ja+1,ka,  l) + aaa3*p1mu
        wqrp(ia+1,ja+1,ka,  l) = wqrp(ia+1,ja+1,ka,  l) + aaa4*p1mu
        wqrp(ia,  ja,  ka+1,l) = wqrp(ia,  ja,  ka+1,l) + aaa5*p1mu
        wqrp(ia+1,ja,  ka+1,l) = wqrp(ia+1,ja,  ka+1,l) + aaa6*p1mu
        wqrp(ia,  ja+1,ka+1,l) = wqrp(ia,  ja+1,ka+1,l) + aaa7*p1mu
        wqrp(ia+1,ja+1,ka+1,l) = wqrp(ia+1,ja+1,ka+1,l) + aaa8*p1mu
      end do
      end do
!$omp end parallel

      else

      do n = 1, marker_num
        ia = max(1,min(lr  -1,int((gc(1,n)-ma_mi_r)*dr1  ) + kr  ))
        ja = max(1,min(lz  -1,int(gc(2,n)          *dz1  ) + kz  ))
        ka = max(1,min(lphi-1,int(gc(3,n)          *dphi1) + kphi))

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

! gc_x(:,n) = 1:r, 2:z, 3:phi, 4:p, 5:mu, 6:weight, 7:active, 8:p_phi, 9:fff, 10:fnrml
        p0 =  gc(6,n)*gc(7,n)
        p1 =  gc(4,n)*p0
        p2 =  gc(4,n)*p1
        mu1 = gc(5,n)*p0
        p3  = gc(4,n)*p2
        p1mu= gc(4,n)*mu1

        dns(ia,  ja,  ka  ) = dns(ia,  ja,  ka  ) + aaa1*p0
        dns(ia+1,ja,  ka  ) = dns(ia+1,ja,  ka  ) + aaa2*p0
        dns(ia,  ja+1,ka  ) = dns(ia,  ja+1,ka  ) + aaa3*p0
        dns(ia+1,ja+1,ka  ) = dns(ia+1,ja+1,ka  ) + aaa4*p0
        dns(ia,  ja,  ka+1) = dns(ia,  ja,  ka+1) + aaa5*p0
        dns(ia+1,ja,  ka+1) = dns(ia+1,ja,  ka+1) + aaa6*p0
        dns(ia,  ja+1,ka+1) = dns(ia,  ja+1,ka+1) + aaa7*p0
        dns(ia+1,ja+1,ka+1) = dns(ia+1,ja+1,ka+1) + aaa8*p0

        mom(ia,  ja,  ka  ) = mom(ia,  ja,  ka  ) + aaa1*p1
        mom(ia+1,ja,  ka  ) = mom(ia+1,ja,  ka  ) + aaa2*p1
        mom(ia,  ja+1,ka  ) = mom(ia,  ja+1,ka  ) + aaa3*p1
        mom(ia+1,ja+1,ka  ) = mom(ia+1,ja+1,ka  ) + aaa4*p1
        mom(ia,  ja,  ka+1) = mom(ia,  ja,  ka+1) + aaa5*p1
        mom(ia+1,ja,  ka+1) = mom(ia+1,ja,  ka+1) + aaa6*p1
        mom(ia,  ja+1,ka+1) = mom(ia,  ja+1,ka+1) + aaa7*p1
        mom(ia+1,ja+1,ka+1) = mom(ia+1,ja+1,ka+1) + aaa8*p1

        ppara(ia,  ja,  ka  ) = ppara(ia,  ja,  ka  ) + aaa1*p2
        ppara(ia+1,ja,  ka  ) = ppara(ia+1,ja,  ka  ) + aaa2*p2
        ppara(ia,  ja+1,ka  ) = ppara(ia,  ja+1,ka  ) + aaa3*p2
        ppara(ia+1,ja+1,ka  ) = ppara(ia+1,ja+1,ka  ) + aaa4*p2
        ppara(ia,  ja,  ka+1) = ppara(ia,  ja,  ka+1) + aaa5*p2
        ppara(ia+1,ja,  ka+1) = ppara(ia+1,ja,  ka+1) + aaa6*p2
        ppara(ia,  ja+1,ka+1) = ppara(ia,  ja+1,ka+1) + aaa7*p2
        ppara(ia+1,ja+1,ka+1) = ppara(ia+1,ja+1,ka+1) + aaa8*p2

        pperp(ia,  ja,  ka  ) = pperp(ia,  ja,  ka  ) + aaa1*mu1
        pperp(ia+1,ja,  ka  ) = pperp(ia+1,ja,  ka  ) + aaa2*mu1
        pperp(ia,  ja+1,ka  ) = pperp(ia,  ja+1,ka  ) + aaa3*mu1
        pperp(ia+1,ja+1,ka  ) = pperp(ia+1,ja+1,ka  ) + aaa4*mu1
        pperp(ia,  ja,  ka+1) = pperp(ia,  ja,  ka+1) + aaa5*mu1
        pperp(ia+1,ja,  ka+1) = pperp(ia+1,ja,  ka+1) + aaa6*mu1
        pperp(ia,  ja+1,ka+1) = pperp(ia,  ja+1,ka+1) + aaa7*mu1
        pperp(ia+1,ja+1,ka+1) = pperp(ia+1,ja+1,ka+1) + aaa8*mu1

        qpara(ia,  ja,  ka  ) = qpara(ia,  ja,  ka  ) + aaa1*p3
        qpara(ia+1,ja,  ka  ) = qpara(ia+1,ja,  ka  ) + aaa2*p3
        qpara(ia,  ja+1,ka  ) = qpara(ia,  ja+1,ka  ) + aaa3*p3
        qpara(ia+1,ja+1,ka  ) = qpara(ia+1,ja+1,ka  ) + aaa4*p3
        qpara(ia,  ja,  ka+1) = qpara(ia,  ja,  ka+1) + aaa5*p3
        qpara(ia+1,ja,  ka+1) = qpara(ia+1,ja,  ka+1) + aaa6*p3
        qpara(ia,  ja+1,ka+1) = qpara(ia,  ja+1,ka+1) + aaa7*p3
        qpara(ia+1,ja+1,ka+1) = qpara(ia+1,ja+1,ka+1) + aaa8*p3

        qperp(ia,  ja,  ka  ) = qperp(ia,  ja,  ka  ) + aaa1*p1mu
        qperp(ia+1,ja,  ka  ) = qperp(ia+1,ja,  ka  ) + aaa2*p1mu
        qperp(ia,  ja+1,ka  ) = qperp(ia,  ja+1,ka  ) + aaa3*p1mu
        qperp(ia+1,ja+1,ka  ) = qperp(ia+1,ja+1,ka  ) + aaa4*p1mu
        qperp(ia,  ja,  ka+1) = qperp(ia,  ja,  ka+1) + aaa5*p1mu
        qperp(ia+1,ja,  ka+1) = qperp(ia+1,ja,  ka+1) + aaa6*p1mu
        qperp(ia,  ja+1,ka+1) = qperp(ia,  ja+1,ka+1) + aaa7*p1mu
        qperp(ia+1,ja+1,ka+1) = qperp(ia+1,ja+1,ka+1) + aaa8*p1mu
      end do

      end if
!2025-07-11e

      if(lpara.ge.2)then

      do k = 1, lpara
!$omp parallel do
      do i = 1, lrzphi
        dns(i,1,1)   = dns(i,1,1)   + wdns(i,1,1,k)
        mom(i,1,1)   = mom(i,1,1)   + wmom(i,1,1,k)
        ppara(i,1,1) = ppara(i,1,1) + wpar(i,1,1,k)
        pperp(i,1,1) = pperp(i,1,1) + wprp(i,1,1,k)
        qpara(i,1,1) = qpara(i,1,1) + wqar(i,1,1,k)
        qperp(i,1,1) = qperp(i,1,1) + wqrp(i,1,1,k)
      end do
      end do

      end if

! smoothing

       cwwa = 0.5d0
!       cwwc =-1.d0/6.d0 

!2024-12-21s, correction suggested by Panith Adulsiriswad
       call periodic_particle_mlt6b(dns,mom,ppara,pperp,qpara,qperp)
       call partsm1(dns,cwwa)
       call partsm1(mom,cwwa)

       call partsm1(ppara,cwwa)
       call partsm1(pperp,cwwa)

       call partsm1(qpara,cwwa)
       call partsm1(qperp,cwwa)
!2024-12-21e

! calculate density (per volume)

!      do k = 1, lphi
!      do j = 1, lz
!      do i = 1, lr
!$omp parallel do private(vol)
      do i = 1, lrzphi
        vol = 1.0d0/(grr(i,1,1)*dr*dz*dphi)
        dns(i,1,1) = dns(i,1,1)*vol
        mom(i,1,1) = mom(i,1,1)*vol
        ppara(i,1,1) = ppara(i,1,1)*vol*mass1
        pperp(i,1,1) = pperp(i,1,1)*vol*babs(i,1,1)
        qpara(i,1,1) = qpara(i,1,1)*vol*mass1**2
        qperp(i,1,1) = qperp(i,1,1)*vol*babs(i,1,1)*mass1
      end do
!      end do
!      end do

! wall effect

      lr1 = lr - 1
      lz1 = lz - 1

      if(my_rank_r.eq.(mpi_proc_r-1))then
!      do k = 1, lphi
!      do j = 1, lz
!$omp parallel do
       do j = 1, lzphi
          dns(lr1,j,1) = dns(lr1,j,1) + dns(lr,j,1)
          mom(lr1,j,1) = mom(lr1,j,1) + mom(lr,j,1)
          dns(lr,j,1) = 0.0d0
          mom(lr,j,1) = 0.0d0

          ppara(lr1,j,1) = ppara(lr1,j,1) + ppara(lr,j,1)
          pperp(lr1,j,1) = pperp(lr1,j,1) + pperp(lr,j,1)
          ppara(lr,j,1) = 0.0d0
          pperp(lr,j,1) = 0.0d0

          qpara(lr1,j,1) = qpara(lr1,j,1) + qpara(lr,j,1)
          qperp(lr1,j,1) = qperp(lr1,j,1) + qperp(lr,j,1)
          qpara(lr,j,1) = 0.0d0
          qperp(lr,j,1) = 0.0d0
       end do
!      end do
      end if

      if(my_rank_r.eq.0)then
!      do k = 1, lphi
!      do j = 1, lz
!$omp parallel do
       do j = 1, lzphi
          dns(2,j,1) = dns(2,j,1) + dns(1,j,1)
          mom(2,j,1) = mom(2,j,1) + mom(1,j,1)
          dns(1,j,1) = 0.0d0
          mom(1,j,1) = 0.0d0

          ppara(2,j,1) = ppara(2,j,1) + ppara(1,j,1)
          pperp(2,j,1) = pperp(2,j,1) + pperp(1,j,1)
          ppara(1,j,1) = 0.0d0
          pperp(1,j,1) = 0.0d0

          qpara(2,j,1) = qpara(2,j,1) + qpara(1,j,1)
          qperp(2,j,1) = qperp(2,j,1) + qperp(1,j,1)
          qpara(1,j,1) = 0.0d0
          qperp(1,j,1) = 0.0d0
       end do
!      end do
      end if

      if(my_rank_z.eq.(mpi_proc_z-1))then
!$omp parallel do
       do k = 1, lphi
       do i = 1, lr
          dns(i,lz1,k) = dns(i,lz1,k) + dns(i,lz,k)
          mom(i,lz1,k) = mom(i,lz1,k) + mom(i,lz,k)
          dns(i,lz,k) = 0.0d0
          mom(i,lz,k) = 0.0d0

          ppara(i,lz1,k) = ppara(i,lz1,k) + ppara(i,lz,k)
          pperp(i,lz1,k) = pperp(i,lz1,k) + pperp(i,lz,k)
          ppara(i,lz,k) = 0.0d0
          pperp(i,lz,k) = 0.0d0

          qpara(i,lz1,k) = qpara(i,lz1,k) + qpara(i,lz,k)
          qperp(i,lz1,k) = qperp(i,lz1,k) + qperp(i,lz,k)
          qpara(i,lz,k) = 0.0d0
          qperp(i,lz,k) = 0.0d0
       end do
       end do
      end if 

      if(my_rank_z.eq.0)then
!$omp parallel do
       do k = 1, lphi
       do i = 1, lr
          dns(i,2,k) = dns(i,2,k) + dns(i,1,k)
          mom(i,2,k) = mom(i,2,k) + mom(i,1,k)
          dns(i,1,k) = 0.0d0
          mom(i,1,k) = 0.0d0

          ppara(i,2,k) = ppara(i,2,k) + ppara(i,1,k)
          pperp(i,2,k) = pperp(i,2,k) + pperp(i,1,k)
          ppara(i,1,k) = 0.0d0
          pperp(i,1,k) = 0.0d0

          qpara(i,2,k) = qpara(i,2,k) + qpara(i,1,k)
          qperp(i,2,k) = qperp(i,2,k) + qperp(i,1,k)
          qpara(i,1,k) = 0.0d0
          qperp(i,1,k) = 0.0d0
       end do
       end do
      end if 

!      if(.not.flag_stored)then
! take modes n_min <= n <=n_max
!      n_min = 1
!      n_max = 2
!      call lowpass2_mlt2(n_min,n_max,dns,mom)
!      call lowpass2_mlt2(n_min,n_max,ppara,pperp)
!      call lowpass2_mlt2(n_min,n_max,qpara,qperp)
!      end if

!2014-08-28s
!$omp parallel do
      do i = 1, lrzphi
        dns(i,1,1) = dns(i,1,1) + dns0(i,1,1)
        mom(i,1,1) = mom(i,1,1) + mom0(i,1,1)
        ppara(i,1,1) = ppara(i,1,1) + ppara0(i,1,1)
        pperp(i,1,1) = pperp(i,1,1) + pperp0(i,1,1)/babs0(i,1,1)*babs(i,1,1)
        qpara(i,1,1) = qpara(i,1,1) + qpara0(i,1,1)
        qperp(i,1,1) = qperp(i,1,1) + qperp0(i,1,1)/babs0(i,1,1)*babs(i,1,1)
      end do
!2014-08-28e

!            call wall_clock(d5)
!            t5 = d5 - d4
!            if(my_rank.eq.0)then
!            write(6,*)'density 5, t5=',t5
!            end if

end subroutine moments
!--------------------------------------------------------------------
subroutine emf_gyro(marker_num,marker_num_gyro,gc,gyro_phys &
                   ,flp,flp_gyro)
! modified for GK simulation on Fugaku 2022-01-11
! modified for Subsystem A of Plasma Simulator with Xeon 6900P, 2025-07-11
! gc_x(:,n) = 1:r, 2:z, 3:phi, 4:p, 5:mu, 6:weight, 7:active, 8:p_phi, 9:fff, 10:fnrml
! flp(nflp, marker_num); nflp=30, nflp_gyro=7
!--------------------------------------------------------------------
      use mpiset
      use grid
      use field, only:babs,babs0,er,ez,ephi,epara,br,bz,bphi,br0,bz0,bphi0,fld
      use gyro, only:ngyro,marker_each_gyro
      implicit none

      integer::marker_num,nvec,marker_num_gyro
      real(8)::gc(ngc2,marker_num)
      real(8)::gyro_phys(2,marker_num_gyro)
      real(8)::flp(nflp,marker_num) !2015-07-08
      real(8)::flp_gyro(nflp_gyro,marker_num_gyro)
      integer::i,j,k,l,m,n,ia,ja,ka,lr1,lz1
      real(8)::ar,ar1,az,az1,aphi,aphi1
      real(8)::aaa1,aaa2,aaa3,aaa4,aaa5,aaa6,aaa7,aaa8 !2025-07-11
      real(8)::dr1,dz1,dphi1,ma_mi_r
      integer::kr,kz,kphi
      integer::nc,in
      real(8)::rngr1,rngr1_ac
      integer::nvec0
      integer::nmod

!2015-07-08s
!      integer::ijk_a(3,marker_each_gyro)
!      real(8)::aaa(8,marker_each_gyro)
!2015-07-08e

      dr1 = 1.0d0/dr
      dz1 = 1.0d0/dz
      dphi1 = 1.0d0/dphi
      rngr1 = 1.0d0/dble(ngyro)

      kr   = 1 - kr_offset(my_rank)
      kz   = 1 - kz_offset(my_rank)
      kphi = 1 + lphi_shd - kphi_offset(my_rank)
!      kphi = 3 - kphi_offset(my_rank)
      ma_mi_r = major_r-minor_r


!$omp parallel
!$omp workshare

! for subroutine extract_em

!      flp_gyro = 0.0d0

      fld(1,:,:,:) = er(:,:,:)
      fld(2,:,:,:) = ez(:,:,:)
      fld(3,:,:,:) = ephi(:,:,:)
      fld(4,:,:,:) = epara(:,:,:)
      fld(5,:,:,:) = br(:,:,:) - br0(:,:,:)
      fld(6,:,:,:) = bz(:,:,:) - bz0(:,:,:)
      fld(7,:,:,:) = bphi(:,:,:) - bphi0(:,:,:)
!$omp end workshare
!$omp end parallel

!2025-07-11s
!$omp parallel private(n,nc,in,ia,ja,ka,ar1,ar,az1,az,aphi1,aphi &
!$omp& ,aaa1,aaa2,aaa3,aaa4,aaa5,aaa6,aaa7,aaa8)
!$omp do schedule(dynamic,1000)
      do n = 1, marker_num_gyro
        nc =(n-1)/ngyro + 1

        ia = max(1,min(lr  -1,int((gyro_phys(1,n)-ma_mi_r)*dr1  ) + kr  ))
        ja = max(1,min(lz  -1,int(gyro_phys(2,n)          *dz1  ) + kz  ))
        ka = max(1,min(lphi-1,int(gc(3,nc)          *dphi1) + kphi))

        ar1  = max(0.0d0, min(1.0d0, (gyro_phys(1,n) - ma_mi_r)*dr1 -  dble(ia - kr)  ) )
        ar   = 1.0d0 - ar1
        az1  = max(0.0d0, min(1.0d0, gyro_phys(2,n)*dz1 - dble(ja - kz) ) )
        az   = 1.0d0 - az1
        aphi1 = max(0.0d0, min(1.0d0, gc(3,nc)*dphi1 - dble(ka - kphi) ) )
        aphi = 1.0d0 - aphi1

        aaa1 = ar *az *aphi
        aaa2 = ar1*az *aphi
        aaa3 = ar *az1*aphi
        aaa4 = ar1*az1*aphi
        aaa5 = ar *az *aphi1
        aaa6 = ar1*az *aphi1
        aaa7 = ar *az1*aphi1
        aaa8 = ar1*az1*aphi1

        do in = 1, nflp_gyro
          flp_gyro(in,n) = fld(in, ia, ja,  ka  )*aaa1 + fld(in, ia+1,ja,  ka  )*aaa2 &
                         + fld(in, ia, ja+1,ka  )*aaa3 + fld(in, ia+1,ja+1,ka  )*aaa4 &
                         + fld(in, ia, ja,  ka+1)*aaa5 + fld(in, ia+1,ja,  ka+1)*aaa6 &
                         + fld(in, ia, ja+1,ka+1)*aaa7 + fld(in, ia+1,ja+1,ka+1)*aaa8
        end do

      end do
!$omp end parallel
!2025-07-11e

!$omp parallel do private(rngr1_ac,in)
      do nc = 1, marker_num
        rngr1_ac = rngr1*gc(7,nc)
        do in = 1, nflp_gyro
          flp(in,nc) =(flp_gyro(in,ngyro*(nc-1)+1) &
                      +flp_gyro(in,ngyro*(nc-1)+2) &
                      +flp_gyro(in,ngyro*(nc-1)+3) &
                      +flp_gyro(in,ngyro*(nc-1)+4) &
                      )*rngr1_ac
        end do
      end do

end subroutine emf_gyro
!--------------------------------------------------------------------
subroutine density_gyro(marker_num,marker_num_gyro,mass,gc,gyro_phys &
                       ,dns,mom,ppara,pperp,dns0,mom0,ppara0,pperp0)
! 2022-01-11, for FLR case
! modified for Subsystem A of Plasma Simulator with Xeon 6900P, 2025-07-11
! gc_x(:,n) = 1:r, 2:z, 3:phi, 4:p, 5:mu, 6:weight, 7:active, 8:p_phi, 9:fff, 10:fnrml
!--------------------------------------------------------------------
      use mpiset
      use grid
      use field, only:babs,babs0
      use gyro, only:ngyro,marker_each_gyro
      implicit none

      integer::marker_num,nvec,marker_num_gyro
      real(8)::gc(ngc2,marker_num)
      real(8)::gyro_phys(2,marker_num_gyro)
      real(8)::mass
      real(8)::dns(lr,lz,lphi),mom(lr,lz,lphi)
      real(8)::ppara(lr,lz,lphi),pperp(lr,lz,lphi)
      real(8)::dns0(lr,lz,lphi),mom0(lr,lz,lphi)
      real(8)::ppara0(lr,lz,lphi),pperp0(lr,lz,lphi)
      real(8)::wdns(lr,lz,lphi,lpara),wmom(lr,lz,lphi,lpara)
      real(8)::wpar(lr,lz,lphi,lpara),wprp(lr,lz,lphi,lpara)
!2015-09-21s
!      integer::ijk_a(3,marker_each_gyro)
!      real(8)::aaa(8,marker_each_gyro)
!2015-09-21e
      integer::i,j,k,l,m,n,ia,ia1,ja,ja1,ka,ka1,lr1,lz1
      real(8)::ar,ar1,az,az1,aphi,aphi1
      real(8)::aaa1,aaa2,aaa3,aaa4,aaa5,aaa6,aaa7,aaa8 !2025-07-11
      real(8)::vol,cwwa,cwwc
!      real(8)::d0,d1,d2,d3,d4,d5,t1,t2,t3,t4,t5
      real(8)::dr1,dz1,dphi1,ma_mi_r
      integer::kr,kz,kphi
      real(8)::p0,p1,p2,mu1,mass1 !2015-09-21
      integer::n_min,n_max !2013-05-22
      real(8)::rngr1
      integer::nc
      integer::nvec0
      integer::nmod

      mass1 = 1.0d0/mass

      dr1 = 1.0d0/dr
      dz1 = 1.0d0/dz
      dphi1 = 1.0d0/dphi
      rngr1 = 1.0d0/dble(ngyro)

      kr   = 1 - kr_offset(my_rank)
      kz   = 1 - kz_offset(my_rank)
      kphi = 1 + lphi_shd - kphi_offset(my_rank)
!      kphi = 3 - kphi_offset(my_rank)
      ma_mi_r = major_r-minor_r

      if(lpara.ge.2)then
!$omp parallel do
      do k = 1, lpara
        do i = 1, lrzphi
          wdns(i,1,1,k) = 0.0d0
          wmom(i,1,1,k) = 0.0d0
          wpar(i,1,1,k) = 0.0d0
          wprp(i,1,1,k) = 0.0d0
        end do
      end do

      end if

!$omp parallel do
      do i = 1, lrzphi
        dns(i,1,1) = 0.0d0
        mom(i,1,1) = 0.0d0
        ppara(i,1,1) = 0.0d0
        pperp(i,1,1) = 0.0d0
      end do

!2025-07-11s
      if(lpara.ge.2)then
        nvec = marker_num_gyro/lpara
        nvec0 = marker_num_gyro/lpara
        nmod  = mod(marker_num_gyro,lpara)

!$omp parallel private(nvec,l,m,n,nc,p0,p1,p2,mu1 &
!$omp& ,ia,ja,ka,ar1,ar,az1,az,aphi1,aphi &
!$omp& ,aaa1,aaa2,aaa3,aaa4,aaa5,aaa6,aaa7,aaa8)
!$omp do
      do l = 1, lpara
        nvec = nvec0 + min(1,nmod/l)
      do m = 1, nvec
        n = m +  nvec0*(l-1) + min(l-1,nmod)
        nc =(n-1)/ngyro + 1

        ia = max(1,min(lr  -1,int((gyro_phys(1,n)-ma_mi_r)*dr1  ) + kr  ))
        ja = max(1,min(lz  -1,int(gyro_phys(2,n)          *dz1  ) + kz  ))
        ka = max(1,min(lphi-1,int(gc(3,nc)          *dphi1) + kphi))

        ar1  = max(0.0d0, min(1.0d0, (gyro_phys(1,n) - ma_mi_r)*dr1 -  dble(ia - kr)  ) )
        ar   = 1.0d0 - ar1
        az1  = max(0.0d0, min(1.0d0, gyro_phys(2,n)*dz1 - dble(ja - kz) ) )
        az   = 1.0d0 - az1
        aphi1 = max(0.0d0, min(1.0d0, gc(3,nc)*dphi1 - dble(ka - kphi) ) )
        aphi = 1.0d0 - aphi1

        aaa1 = ar *az *aphi
        aaa2 = ar1*az *aphi
        aaa3 = ar *az1*aphi
        aaa4 = ar1*az1*aphi
        aaa5 = ar *az *aphi1
        aaa6 = ar1*az *aphi1
        aaa7 = ar *az1*aphi1
        aaa8 = ar1*az1*aphi1

! gc_x(:,n) = 1:r, 2:z, 3:phi, 4:p, 5:mu, 6:weight, 7:active, 8:p_phi, 9:fff, 10:fnrml
        p0 =  gc(6,nc)*gc(7,nc)*rngr1 !rngr1 is multiplied
        p1 =  gc(4,nc)*p0
        p2 =  gc(4,nc)*p1
        mu1 = gc(5,nc)*p0

        wdns(ia,  ja,  ka,  l) = wdns(ia,  ja,  ka,  l) + aaa1*p0
        wdns(ia+1,ja,  ka,  l) = wdns(ia+1,ja,  ka,  l) + aaa2*p0
        wdns(ia,  ja+1,ka,  l) = wdns(ia,  ja+1,ka,  l) + aaa3*p0
        wdns(ia+1,ja+1,ka,  l) = wdns(ia+1,ja+1,ka,  l) + aaa4*p0
        wdns(ia,  ja,  ka+1,l) = wdns(ia,  ja,  ka+1,l) + aaa5*p0
        wdns(ia+1,ja,  ka+1,l) = wdns(ia+1,ja,  ka+1,l) + aaa6*p0
        wdns(ia,  ja+1,ka+1,l) = wdns(ia,  ja+1,ka+1,l) + aaa7*p0
        wdns(ia+1,ja+1,ka+1,l) = wdns(ia+1,ja+1,ka+1,l) + aaa8*p0

        wmom(ia,  ja,  ka,  l) = wmom(ia,  ja,  ka,  l) + aaa1*p1
        wmom(ia+1,ja,  ka,  l) = wmom(ia+1,ja,  ka,  l) + aaa2*p1
        wmom(ia,  ja+1,ka,  l) = wmom(ia,  ja+1,ka,  l) + aaa3*p1
        wmom(ia+1,ja+1,ka,  l) = wmom(ia+1,ja+1,ka,  l) + aaa4*p1
        wmom(ia,  ja,  ka+1,l) = wmom(ia,  ja,  ka+1,l) + aaa5*p1
        wmom(ia+1,ja,  ka+1,l) = wmom(ia+1,ja,  ka+1,l) + aaa6*p1
        wmom(ia,  ja+1,ka+1,l) = wmom(ia,  ja+1,ka+1,l) + aaa7*p1
        wmom(ia+1,ja+1,ka+1,l) = wmom(ia+1,ja+1,ka+1,l) + aaa8*p1

        wpar(ia,  ja,  ka,  l) = wpar(ia,  ja,  ka,  l) + aaa1*p2
        wpar(ia+1,ja,  ka,  l) = wpar(ia+1,ja,  ka,  l) + aaa2*p2
        wpar(ia,  ja+1,ka,  l) = wpar(ia,  ja+1,ka,  l) + aaa3*p2
        wpar(ia+1,ja+1,ka,  l) = wpar(ia+1,ja+1,ka,  l) + aaa4*p2
        wpar(ia,  ja,  ka+1,l) = wpar(ia,  ja,  ka+1,l) + aaa5*p2
        wpar(ia+1,ja,  ka+1,l) = wpar(ia+1,ja,  ka+1,l) + aaa6*p2
        wpar(ia,  ja+1,ka+1,l) = wpar(ia,  ja+1,ka+1,l) + aaa7*p2
        wpar(ia+1,ja+1,ka+1,l) = wpar(ia+1,ja+1,ka+1,l) + aaa8*p2

        wprp(ia,  ja,  ka,  l) = wprp(ia,  ja,  ka,  l) + aaa1*mu1
        wprp(ia+1,ja,  ka,  l) = wprp(ia+1,ja,  ka,  l) + aaa2*mu1
        wprp(ia,  ja+1,ka,  l) = wprp(ia,  ja+1,ka,  l) + aaa3*mu1
        wprp(ia+1,ja+1,ka,  l) = wprp(ia+1,ja+1,ka,  l) + aaa4*mu1
        wprp(ia,  ja,  ka+1,l) = wprp(ia,  ja,  ka+1,l) + aaa5*mu1
        wprp(ia+1,ja,  ka+1,l) = wprp(ia+1,ja,  ka+1,l) + aaa6*mu1
        wprp(ia,  ja+1,ka+1,l) = wprp(ia,  ja+1,ka+1,l) + aaa7*mu1
        wprp(ia+1,ja+1,ka+1,l) = wprp(ia+1,ja+1,ka+1,l) + aaa8*mu1
      end do
      end do
!$omp end parallel

      else

      do n = 1, marker_num_gyro
        nc =(n-1)/ngyro + 1

        ia = max(1,min(lr  -1,int((gyro_phys(1,n)-ma_mi_r)*dr1  ) + kr  ))
        ja = max(1,min(lz  -1,int(gyro_phys(2,n)          *dz1  ) + kz  ))
        ka = max(1,min(lphi-1,int(gc(3,nc)          *dphi1) + kphi))

        ar1  = max(0.0d0, min(1.0d0, (gyro_phys(1,n) - ma_mi_r)*dr1 -  dble(ia - kr)  ) )
        ar   = 1.0d0 - ar1
        az1  = max(0.0d0, min(1.0d0, gyro_phys(2,n)*dz1 - dble(ja - kz) ) )
        az   = 1.0d0 - az1
        aphi1 = max(0.0d0, min(1.0d0, gc(3,nc)*dphi1 - dble(ka - kphi) ) )
        aphi = 1.0d0 - aphi1

        aaa1 = ar *az *aphi
        aaa2 = ar1*az *aphi
        aaa3 = ar *az1*aphi
        aaa4 = ar1*az1*aphi
        aaa5 = ar *az *aphi1
        aaa6 = ar1*az *aphi1
        aaa7 = ar *az1*aphi1
        aaa8 = ar1*az1*aphi1

! gc_x(:,n) = 1:r, 2:z, 3:phi, 4:p, 5:mu, 6:weight, 7:active, 8:p_phi, 9:fff, 10:fnrml
        p0 =  gc(6,nc)*gc(7,nc)*rngr1 !rngr1 is multiplied
        p1 =  gc(4,nc)*p0
        p2 =  gc(4,nc)*p1
        mu1 = gc(5,nc)*p0

        dns(ia,  ja,  ka  ) = dns(ia,  ja,  ka  ) + aaa1*p0
        dns(ia+1,ja,  ka  ) = dns(ia+1,ja,  ka  ) + aaa2*p0
        dns(ia,  ja+1,ka  ) = dns(ia,  ja+1,ka  ) + aaa3*p0
        dns(ia+1,ja+1,ka  ) = dns(ia+1,ja+1,ka  ) + aaa4*p0
        dns(ia,  ja,  ka+1) = dns(ia,  ja,  ka+1) + aaa5*p0
        dns(ia+1,ja,  ka+1) = dns(ia+1,ja,  ka+1) + aaa6*p0
        dns(ia,  ja+1,ka+1) = dns(ia,  ja+1,ka+1) + aaa7*p0
        dns(ia+1,ja+1,ka+1) = dns(ia+1,ja+1,ka+1) + aaa8*p0

        mom(ia,  ja,  ka  ) = mom(ia,  ja,  ka  ) + aaa1*p1
        mom(ia+1,ja,  ka  ) = mom(ia+1,ja,  ka  ) + aaa2*p1
        mom(ia,  ja+1,ka  ) = mom(ia,  ja+1,ka  ) + aaa3*p1
        mom(ia+1,ja+1,ka  ) = mom(ia+1,ja+1,ka  ) + aaa4*p1
        mom(ia,  ja,  ka+1) = mom(ia,  ja,  ka+1) + aaa5*p1
        mom(ia+1,ja,  ka+1) = mom(ia+1,ja,  ka+1) + aaa6*p1
        mom(ia,  ja+1,ka+1) = mom(ia,  ja+1,ka+1) + aaa7*p1
        mom(ia+1,ja+1,ka+1) = mom(ia+1,ja+1,ka+1) + aaa8*p1

        ppara(ia,  ja,  ka  ) = ppara(ia,  ja,  ka  ) + aaa1*p2
        ppara(ia+1,ja,  ka  ) = ppara(ia+1,ja,  ka  ) + aaa2*p2
        ppara(ia,  ja+1,ka  ) = ppara(ia,  ja+1,ka  ) + aaa3*p2
        ppara(ia+1,ja+1,ka  ) = ppara(ia+1,ja+1,ka  ) + aaa4*p2
        ppara(ia,  ja,  ka+1) = ppara(ia,  ja,  ka+1) + aaa5*p2
        ppara(ia+1,ja,  ka+1) = ppara(ia+1,ja,  ka+1) + aaa6*p2
        ppara(ia,  ja+1,ka+1) = ppara(ia,  ja+1,ka+1) + aaa7*p2
        ppara(ia+1,ja+1,ka+1) = ppara(ia+1,ja+1,ka+1) + aaa8*p2

        pperp(ia,  ja,  ka  ) = pperp(ia,  ja,  ka  ) + aaa1*mu1
        pperp(ia+1,ja,  ka  ) = pperp(ia+1,ja,  ka  ) + aaa2*mu1
        pperp(ia,  ja+1,ka  ) = pperp(ia,  ja+1,ka  ) + aaa3*mu1
        pperp(ia+1,ja+1,ka  ) = pperp(ia+1,ja+1,ka  ) + aaa4*mu1
        pperp(ia,  ja,  ka+1) = pperp(ia,  ja,  ka+1) + aaa5*mu1
        pperp(ia+1,ja,  ka+1) = pperp(ia+1,ja,  ka+1) + aaa6*mu1
        pperp(ia,  ja+1,ka+1) = pperp(ia,  ja+1,ka+1) + aaa7*mu1
        pperp(ia+1,ja+1,ka+1) = pperp(ia+1,ja+1,ka+1) + aaa8*mu1
      end do

      end if
!2025-07-11e


      if(lpara.ge.2)then

      do k = 1, lpara
!$omp parallel do
      do i = 1, lrzphi
        dns(i,1,1)   = dns(i,1,1)   + wdns(i,1,1,k)
        mom(i,1,1)   = mom(i,1,1)   + wmom(i,1,1,k)
        ppara(i,1,1) = ppara(i,1,1) + wpar(i,1,1,k)
        pperp(i,1,1) = pperp(i,1,1) + wprp(i,1,1,k)
      end do
      end do

      end if


! smoothing

       cwwa = 0.5d0
!       cwwc =-1.d0/6.d0 

!2024-12-21s, correction suggested by Panith Adulsiriswad
       call periodic_particle_mlt4b(dns,mom,ppara,pperp)
       call partsm1(dns,cwwa)
       call partsm1(mom,cwwa)

       call partsm1(ppara,cwwa)
       call partsm1(pperp,cwwa)
!2024-12-21e
      
! calculate density (per volume)

!      do k = 1, lphi
!      do j = 1, lz
!      do i = 1, lr
!$omp parallel do private(vol)
      do i = 1, lrzphi
        vol = 1.0d0/(grr(i,1,1)*dr*dz*dphi)
        dns(i,1,1) = dns(i,1,1)*vol
        mom(i,1,1) = mom(i,1,1)*vol
        ppara(i,1,1) = ppara(i,1,1)*vol/mass
        pperp(i,1,1) = pperp(i,1,1)*vol*babs(i,1,1)
      end do
!      end do
!      end do


! wall effect

      lr1 = lr - 1
      lz1 = lz - 1

      if(my_rank_r.eq.(mpi_proc_r-1))then
!      do k = 1, lphi
!      do j = 1, lz
!$omp parallel do
       do j = 1, lzphi
          dns(lr1,j,1) = dns(lr1,j,1) + dns(lr,j,1)
          mom(lr1,j,1) = mom(lr1,j,1) + mom(lr,j,1)
          dns(lr,j,1) = 0.0d0
          mom(lr,j,1) = 0.0d0

          ppara(lr1,j,1) = ppara(lr1,j,1) + ppara(lr,j,1)
          pperp(lr1,j,1) = pperp(lr1,j,1) + pperp(lr,j,1)
          ppara(lr,j,1) = 0.0d0
          pperp(lr,j,1) = 0.0d0
       end do
!      end do
      end if

      if(my_rank_r.eq.0)then
!      do k = 1, lphi
!      do j = 1, lz
!$omp parallel do
       do j = 1, lzphi
          dns(2,j,1) = dns(2,j,1) + dns(1,j,1)
          mom(2,j,1) = mom(2,j,1) + mom(1,j,1)
          dns(1,j,1) = 0.0d0
          mom(1,j,1) = 0.0d0

          ppara(2,j,1) = ppara(2,j,1) + ppara(1,j,1)
          pperp(2,j,1) = pperp(2,j,1) + pperp(1,j,1)
          ppara(1,j,1) = 0.0d0
          pperp(1,j,1) = 0.0d0
       end do
!      end do
      end if

      if(my_rank_z.eq.(mpi_proc_z-1))then
!$omp parallel do
       do k = 1, lphi
       do i = 1, lr
          dns(i,lz1,k) = dns(i,lz1,k) + dns(i,lz,k)
          mom(i,lz1,k) = mom(i,lz1,k) + mom(i,lz,k)
          dns(i,lz,k) = 0.0d0
          mom(i,lz,k) = 0.0d0

          ppara(i,lz1,k) = ppara(i,lz1,k) + ppara(i,lz,k)
          pperp(i,lz1,k) = pperp(i,lz1,k) + pperp(i,lz,k)
          ppara(i,lz,k) = 0.0d0
          pperp(i,lz,k) = 0.0d0
       end do
       end do
      end if 

      if(my_rank_z.eq.0)then
!$omp parallel do
       do k = 1, lphi
       do i = 1, lr
          dns(i,2,k) = dns(i,2,k) + dns(i,1,k)
          mom(i,2,k) = mom(i,2,k) + mom(i,1,k)
          dns(i,1,k) = 0.0d0
          mom(i,1,k) = 0.0d0

          ppara(i,2,k) = ppara(i,2,k) + ppara(i,1,k)
          pperp(i,2,k) = pperp(i,2,k) + pperp(i,1,k)
          ppara(i,1,k) = 0.0d0
          pperp(i,1,k) = 0.0d0
       end do
       end do
      end if 

! take only n=1 modes
!      n = 1
!      call n1(n,ppara)
!      call n1(n,pperp)

! take modes n_min <= n <=n_max
!      n_min = 1
!      n_max = 2
!      call lowpass2_mlt2(n_min,n_max,dns,mom)
!      call lowpass2_mlt2(n_min,n_max,ppara,pperp)

!      do k = 1, lphi
!      do j = 1, lz
!      do i = 1, lr
!$omp parallel do
      do i = 1, lrzphi
        dns(i,1,1) = dns(i,1,1) + dns0(i,1,1)
        mom(i,1,1) = mom(i,1,1) + mom0(i,1,1)
        ppara(i,1,1) = ppara(i,1,1) + ppara0(i,1,1)
        pperp(i,1,1) = pperp(i,1,1) + pperp0(i,1,1)/babs0(i,1,1)*babs(i,1,1)
      end do
!      end do
!      end do

end subroutine density_gyro
!--------------------------------------------------------------------
subroutine moments_gyro(marker_num,marker_num_gyro,mass,gc,gyro_phys &
                  ,dns,mom,ppara,pperp,qpara,qperp &
                  ,dns0,mom0,ppara0,pperp0,qpara0,qperp0)
! 2016-02-04, for FLR case
! modified for Subsystem A of Plasma Simulator with Xeon 6900P, 2025-07-11
! gc_x(:,n) = 1:r, 2:z, 3:phi, 4:p, 5:mu, 6:weight, 7:active, 8:p_phi, 9:fff, 10:fnrml
!--------------------------------------------------------------------
      use mpiset
      use grid
      use field, only:babs,babs0
      use gyro, only:ngyro,marker_each_gyro
      implicit none

      integer::marker_num,nvec,marker_num_gyro
      real(8)::gc(ngc2,marker_num)
      real(8)::gyro_phys(2,marker_num_gyro)
      real(8)::mass
      real(8)::dns(lr,lz,lphi),mom(lr,lz,lphi)
      real(8)::ppara(lr,lz,lphi),pperp(lr,lz,lphi)
      real(8)::qpara(lr,lz,lphi),qperp(lr,lz,lphi)
      real(8)::dns0(lr,lz,lphi),mom0(lr,lz,lphi)
      real(8)::ppara0(lr,lz,lphi),pperp0(lr,lz,lphi)
      real(8)::qpara0(lr,lz,lphi),qperp0(lr,lz,lphi)
      real(8)::wdns(lr,lz,lphi,lpara),wmom(lr,lz,lphi,lpara)
      real(8)::wpar(lr,lz,lphi,lpara),wprp(lr,lz,lphi,lpara)
      real(8)::wqar(lr,lz,lphi,lpara),wqrp(lr,lz,lphi,lpara)
!2015-09-21s
!      integer::ijk_a(3,marker_each_gyro)
!      real(8)::aaa(8,marker_each_gyro)
!2015-09-21e
      integer::i,j,k,l,m,n,ia,ia1,ja,ja1,ka,ka1,lr1,lz1
      real(8)::ar,ar1,az,az1,aphi,aphi1
      real(8)::aaa1,aaa2,aaa3,aaa4,aaa5,aaa6,aaa7,aaa8 !2025-07-11e
      real(8)::vol,cwwa,cwwc
!      real(8)::d0,d1,d2,d3,d4,d5,t1,t2,t3,t4,t5
      real(8)::dr1,dz1,dphi1,ma_mi_r
      integer::kr,kz,kphi
      real(8)::p0,p1,p2,mu1,p3,p1mu,mass1 !2015-09-21
      integer::n_min,n_max !2013-05-22
      real(8)::rngr1
      integer::nc
      integer::nvec0
      integer::nmod

      mass1 = 1.0d0/mass

      dr1 = 1.0d0/dr
      dz1 = 1.0d0/dz
      dphi1 = 1.0d0/dphi
      rngr1 = 1.0d0/dble(ngyro)

      kr   = 1 - kr_offset(my_rank)
      kz   = 1 - kz_offset(my_rank)
      kphi = 1 + lphi_shd - kphi_offset(my_rank)
!      kphi = 3 - kphi_offset(my_rank)
      ma_mi_r = major_r-minor_r

      if(lpara.ge.2)then
!$omp parallel do
      do k = 1, lpara
        do i = 1, lrzphi
          wdns(i,1,1,k) = 0.0d0
          wmom(i,1,1,k) = 0.0d0
          wpar(i,1,1,k) = 0.0d0
          wprp(i,1,1,k) = 0.0d0
          wqar(i,1,1,k) = 0.0d0
          wqrp(i,1,1,k) = 0.0d0
        end do
      end do

      end if

!$omp parallel do
      do i = 1, lrzphi
        dns(i,1,1) = 0.0d0
        mom(i,1,1) = 0.0d0
        ppara(i,1,1) = 0.0d0
        pperp(i,1,1) = 0.0d0
        qpara(i,1,1) = 0.0d0
        qperp(i,1,1) = 0.0d0
      end do

!2025-07-11s
      if(lpara.ge.2)then
        nvec = marker_num_gyro/lpara
        nvec0 = marker_num_gyro/lpara
        nmod  = mod(marker_num_gyro,lpara)

!$omp parallel private(nvec,l,m,n,nc,p0,p1,p2,mu1,p3,p1mu &
!$omp& ,ia,ja,ka,ar1,ar,az1,az,aphi1,aphi &
!$omp& ,aaa1,aaa2,aaa3,aaa4,aaa5,aaa6,aaa7,aaa8)
!$omp do
      do l = 1, lpara
        nvec = nvec0 + min(1,nmod/l)
      do m = 1, nvec
        n = m +  nvec0*(l-1) + min(l-1,nmod)
        nc =(n-1)/ngyro + 1

        ia = max(1,min(lr  -1,int((gyro_phys(1,n)-ma_mi_r)*dr1  ) + kr  ))
        ja = max(1,min(lz  -1,int(gyro_phys(2,n)          *dz1  ) + kz  ))
        ka = max(1,min(lphi-1,int(gc(3,nc)          *dphi1) + kphi))

        ar1  = max(0.0d0, min(1.0d0, (gyro_phys(1,n) - ma_mi_r)*dr1 -  dble(ia - kr)  ) )
        ar   = 1.0d0 - ar1
        az1  = max(0.0d0, min(1.0d0, gyro_phys(2,n)*dz1 - dble(ja - kz) ) )
        az   = 1.0d0 - az1
        aphi1 = max(0.0d0, min(1.0d0, gc(3,nc)*dphi1 - dble(ka - kphi) ) )
        aphi = 1.0d0 - aphi1

        aaa1 = ar *az *aphi
        aaa2 = ar1*az *aphi
        aaa3 = ar *az1*aphi
        aaa4 = ar1*az1*aphi
        aaa5 = ar *az *aphi1
        aaa6 = ar1*az *aphi1
        aaa7 = ar *az1*aphi1
        aaa8 = ar1*az1*aphi1

! gc_x(:,n) = 1:r, 2:z, 3:phi, 4:p, 5:mu, 6:weight, 7:active, 8:p_phi, 9:fff, 10:fnrml
        p0 =  gc(6,nc)*gc(7,nc)*rngr1 !rngr1 is multiplied
        p1 =  gc(4,nc)*p0
        p2 =  gc(4,nc)*p1
        mu1 = gc(5,nc)*p0
        p3  = gc(4,nc)*p2
        p1mu= gc(4,nc)*mu1

        wdns(ia,  ja,  ka,  l) = wdns(ia,  ja,  ka,  l) + aaa1*p0
        wdns(ia+1,ja,  ka,  l) = wdns(ia+1,ja,  ka,  l) + aaa2*p0
        wdns(ia,  ja+1,ka,  l) = wdns(ia,  ja+1,ka,  l) + aaa3*p0
        wdns(ia+1,ja+1,ka,  l) = wdns(ia+1,ja+1,ka,  l) + aaa4*p0
        wdns(ia,  ja,  ka+1,l) = wdns(ia,  ja,  ka+1,l) + aaa5*p0
        wdns(ia+1,ja,  ka+1,l) = wdns(ia+1,ja,  ka+1,l) + aaa6*p0
        wdns(ia,  ja+1,ka+1,l) = wdns(ia,  ja+1,ka+1,l) + aaa7*p0
        wdns(ia+1,ja+1,ka+1,l) = wdns(ia+1,ja+1,ka+1,l) + aaa8*p0

        wmom(ia,  ja,  ka,  l) = wmom(ia,  ja,  ka,  l) + aaa1*p1
        wmom(ia+1,ja,  ka,  l) = wmom(ia+1,ja,  ka,  l) + aaa2*p1
        wmom(ia,  ja+1,ka,  l) = wmom(ia,  ja+1,ka,  l) + aaa3*p1
        wmom(ia+1,ja+1,ka,  l) = wmom(ia+1,ja+1,ka,  l) + aaa4*p1
        wmom(ia,  ja,  ka+1,l) = wmom(ia,  ja,  ka+1,l) + aaa5*p1
        wmom(ia+1,ja,  ka+1,l) = wmom(ia+1,ja,  ka+1,l) + aaa6*p1
        wmom(ia,  ja+1,ka+1,l) = wmom(ia,  ja+1,ka+1,l) + aaa7*p1
        wmom(ia+1,ja+1,ka+1,l) = wmom(ia+1,ja+1,ka+1,l) + aaa8*p1

        wpar(ia,  ja,  ka,  l) = wpar(ia,  ja,  ka,  l) + aaa1*p2
        wpar(ia+1,ja,  ka,  l) = wpar(ia+1,ja,  ka,  l) + aaa2*p2
        wpar(ia,  ja+1,ka,  l) = wpar(ia,  ja+1,ka,  l) + aaa3*p2
        wpar(ia+1,ja+1,ka,  l) = wpar(ia+1,ja+1,ka,  l) + aaa4*p2
        wpar(ia,  ja,  ka+1,l) = wpar(ia,  ja,  ka+1,l) + aaa5*p2
        wpar(ia+1,ja,  ka+1,l) = wpar(ia+1,ja,  ka+1,l) + aaa6*p2
        wpar(ia,  ja+1,ka+1,l) = wpar(ia,  ja+1,ka+1,l) + aaa7*p2
        wpar(ia+1,ja+1,ka+1,l) = wpar(ia+1,ja+1,ka+1,l) + aaa8*p2

        wprp(ia,  ja,  ka,  l) = wprp(ia,  ja,  ka,  l) + aaa1*mu1
        wprp(ia+1,ja,  ka,  l) = wprp(ia+1,ja,  ka,  l) + aaa2*mu1
        wprp(ia,  ja+1,ka,  l) = wprp(ia,  ja+1,ka,  l) + aaa3*mu1
        wprp(ia+1,ja+1,ka,  l) = wprp(ia+1,ja+1,ka,  l) + aaa4*mu1
        wprp(ia,  ja,  ka+1,l) = wprp(ia,  ja,  ka+1,l) + aaa5*mu1
        wprp(ia+1,ja,  ka+1,l) = wprp(ia+1,ja,  ka+1,l) + aaa6*mu1
        wprp(ia,  ja+1,ka+1,l) = wprp(ia,  ja+1,ka+1,l) + aaa7*mu1
        wprp(ia+1,ja+1,ka+1,l) = wprp(ia+1,ja+1,ka+1,l) + aaa8*mu1

        wqar(ia,  ja,  ka,  l) = wqar(ia,  ja,  ka,  l) + aaa1*p3
        wqar(ia+1,ja,  ka,  l) = wqar(ia+1,ja,  ka,  l) + aaa2*p3
        wqar(ia,  ja+1,ka,  l) = wqar(ia,  ja+1,ka,  l) + aaa3*p3
        wqar(ia+1,ja+1,ka,  l) = wqar(ia+1,ja+1,ka,  l) + aaa4*p3
        wqar(ia,  ja,  ka+1,l) = wqar(ia,  ja,  ka+1,l) + aaa5*p3
        wqar(ia+1,ja,  ka+1,l) = wqar(ia+1,ja,  ka+1,l) + aaa6*p3
        wqar(ia,  ja+1,ka+1,l) = wqar(ia,  ja+1,ka+1,l) + aaa7*p3
        wqar(ia+1,ja+1,ka+1,l) = wqar(ia+1,ja+1,ka+1,l) + aaa8*p3

        wqrp(ia,  ja,  ka,  l) = wqrp(ia,  ja,  ka,  l) + aaa1*p1mu
        wqrp(ia+1,ja,  ka,  l) = wqrp(ia+1,ja,  ka,  l) + aaa2*p1mu
        wqrp(ia,  ja+1,ka,  l) = wqrp(ia,  ja+1,ka,  l) + aaa3*p1mu
        wqrp(ia+1,ja+1,ka,  l) = wqrp(ia+1,ja+1,ka,  l) + aaa4*p1mu
        wqrp(ia,  ja,  ka+1,l) = wqrp(ia,  ja,  ka+1,l) + aaa5*p1mu
        wqrp(ia+1,ja,  ka+1,l) = wqrp(ia+1,ja,  ka+1,l) + aaa6*p1mu
        wqrp(ia,  ja+1,ka+1,l) = wqrp(ia,  ja+1,ka+1,l) + aaa7*p1mu
        wqrp(ia+1,ja+1,ka+1,l) = wqrp(ia+1,ja+1,ka+1,l) + aaa8*p1mu

      end do
      end do
!$omp end parallel

      else

      do n = 1, marker_num_gyro
        nc =(n-1)/ngyro + 1

        ia = max(1,min(lr  -1,int((gyro_phys(1,n)-ma_mi_r)*dr1  ) + kr  ))
        ja = max(1,min(lz  -1,int(gyro_phys(2,n)          *dz1  ) + kz  ))
        ka = max(1,min(lphi-1,int(gc(3,nc)          *dphi1) + kphi))

        ar1  = max(0.0d0, min(1.0d0, (gyro_phys(1,n) - ma_mi_r)*dr1 -  dble(ia - kr)  ) )
        ar   = 1.0d0 - ar1
        az1  = max(0.0d0, min(1.0d0, gyro_phys(2,n)*dz1 - dble(ja - kz) ) )
        az   = 1.0d0 - az1
        aphi1 = max(0.0d0, min(1.0d0, gc(3,nc)*dphi1 - dble(ka - kphi) ) )
        aphi = 1.0d0 - aphi1

        aaa1 = ar *az *aphi
        aaa2 = ar1*az *aphi
        aaa3 = ar *az1*aphi
        aaa4 = ar1*az1*aphi
        aaa5 = ar *az *aphi1
        aaa6 = ar1*az *aphi1
        aaa7 = ar *az1*aphi1
        aaa8 = ar1*az1*aphi1

! gc_x(:,n) = 1:r, 2:z, 3:phi, 4:p, 5:mu, 6:weight, 7:active, 8:p_phi, 9:fff, 10:fnrml
        p0 =  gc(6,nc)*gc(7,nc)*rngr1 !rngr1 is multiplied
        p1 =  gc(4,nc)*p0
        p2 =  gc(4,nc)*p1
        mu1 = gc(5,nc)*p0
        p3  = gc(4,nc)*p2
        p1mu= gc(4,nc)*mu1

        dns(ia,  ja,  ka  ) = dns(ia,  ja,  ka  ) + aaa1*p0
        dns(ia+1,ja,  ka  ) = dns(ia+1,ja,  ka  ) + aaa2*p0
        dns(ia,  ja+1,ka  ) = dns(ia,  ja+1,ka  ) + aaa3*p0
        dns(ia+1,ja+1,ka  ) = dns(ia+1,ja+1,ka  ) + aaa4*p0
        dns(ia,  ja,  ka+1) = dns(ia,  ja,  ka+1) + aaa5*p0
        dns(ia+1,ja,  ka+1) = dns(ia+1,ja,  ka+1) + aaa6*p0
        dns(ia,  ja+1,ka+1) = dns(ia,  ja+1,ka+1) + aaa7*p0
        dns(ia+1,ja+1,ka+1) = dns(ia+1,ja+1,ka+1) + aaa8*p0

        mom(ia,  ja,  ka  ) = mom(ia,  ja,  ka  ) + aaa1*p1
        mom(ia+1,ja,  ka  ) = mom(ia+1,ja,  ka  ) + aaa2*p1
        mom(ia,  ja+1,ka  ) = mom(ia,  ja+1,ka  ) + aaa3*p1
        mom(ia+1,ja+1,ka  ) = mom(ia+1,ja+1,ka  ) + aaa4*p1
        mom(ia,  ja,  ka+1) = mom(ia,  ja,  ka+1) + aaa5*p1
        mom(ia+1,ja,  ka+1) = mom(ia+1,ja,  ka+1) + aaa6*p1
        mom(ia,  ja+1,ka+1) = mom(ia,  ja+1,ka+1) + aaa7*p1
        mom(ia+1,ja+1,ka+1) = mom(ia+1,ja+1,ka+1) + aaa8*p1

        ppara(ia,  ja,  ka  ) = ppara(ia,  ja,  ka  ) + aaa1*p2
        ppara(ia+1,ja,  ka  ) = ppara(ia+1,ja,  ka  ) + aaa2*p2
        ppara(ia,  ja+1,ka  ) = ppara(ia,  ja+1,ka  ) + aaa3*p2
        ppara(ia+1,ja+1,ka  ) = ppara(ia+1,ja+1,ka  ) + aaa4*p2
        ppara(ia,  ja,  ka+1) = ppara(ia,  ja,  ka+1) + aaa5*p2
        ppara(ia+1,ja,  ka+1) = ppara(ia+1,ja,  ka+1) + aaa6*p2
        ppara(ia,  ja+1,ka+1) = ppara(ia,  ja+1,ka+1) + aaa7*p2
        ppara(ia+1,ja+1,ka+1) = ppara(ia+1,ja+1,ka+1) + aaa8*p2

        pperp(ia,  ja,  ka  ) = pperp(ia,  ja,  ka  ) + aaa1*mu1
        pperp(ia+1,ja,  ka  ) = pperp(ia+1,ja,  ka  ) + aaa2*mu1
        pperp(ia,  ja+1,ka  ) = pperp(ia,  ja+1,ka  ) + aaa3*mu1
        pperp(ia+1,ja+1,ka  ) = pperp(ia+1,ja+1,ka  ) + aaa4*mu1
        pperp(ia,  ja,  ka+1) = pperp(ia,  ja,  ka+1) + aaa5*mu1
        pperp(ia+1,ja,  ka+1) = pperp(ia+1,ja,  ka+1) + aaa6*mu1
        pperp(ia,  ja+1,ka+1) = pperp(ia,  ja+1,ka+1) + aaa7*mu1
        pperp(ia+1,ja+1,ka+1) = pperp(ia+1,ja+1,ka+1) + aaa8*mu1

        qpara(ia,  ja,  ka  ) = qpara(ia,  ja,  ka  ) + aaa1*p3
        qpara(ia+1,ja,  ka  ) = qpara(ia+1,ja,  ka  ) + aaa2*p3
        qpara(ia,  ja+1,ka  ) = qpara(ia,  ja+1,ka  ) + aaa3*p3
        qpara(ia+1,ja+1,ka  ) = qpara(ia+1,ja+1,ka  ) + aaa4*p3
        qpara(ia,  ja,  ka+1) = qpara(ia,  ja,  ka+1) + aaa5*p3
        qpara(ia+1,ja,  ka+1) = qpara(ia+1,ja,  ka+1) + aaa6*p3
        qpara(ia,  ja+1,ka+1) = qpara(ia,  ja+1,ka+1) + aaa7*p3
        qpara(ia+1,ja+1,ka+1) = qpara(ia+1,ja+1,ka+1) + aaa8*p3

        qperp(ia,  ja,  ka  ) = qperp(ia,  ja,  ka  ) + aaa1*p1mu
        qperp(ia+1,ja,  ka  ) = qperp(ia+1,ja,  ka  ) + aaa2*p1mu
        qperp(ia,  ja+1,ka  ) = qperp(ia,  ja+1,ka  ) + aaa3*p1mu
        qperp(ia+1,ja+1,ka  ) = qperp(ia+1,ja+1,ka  ) + aaa4*p1mu
        qperp(ia,  ja,  ka+1) = qperp(ia,  ja,  ka+1) + aaa5*p1mu
        qperp(ia+1,ja,  ka+1) = qperp(ia+1,ja,  ka+1) + aaa6*p1mu
        qperp(ia,  ja+1,ka+1) = qperp(ia,  ja+1,ka+1) + aaa7*p1mu
        qperp(ia+1,ja+1,ka+1) = qperp(ia+1,ja+1,ka+1) + aaa8*p1mu

      end do

      end if
!2025-07-11e

      if(lpara.ge.2)then

      do k = 1, lpara
!$omp parallel do
      do i = 1, lrzphi
        dns(i,1,1)   = dns(i,1,1)   + wdns(i,1,1,k)
        mom(i,1,1)   = mom(i,1,1)   + wmom(i,1,1,k)
        ppara(i,1,1) = ppara(i,1,1) + wpar(i,1,1,k)
        pperp(i,1,1) = pperp(i,1,1) + wprp(i,1,1,k)
        qpara(i,1,1) = qpara(i,1,1) + wqar(i,1,1,k)
        qperp(i,1,1) = qperp(i,1,1) + wqrp(i,1,1,k)
      end do
      end do

      end if


! smoothing

       cwwa = 0.5d0
!       cwwc =-1.d0/6.d0 

!2024-12-21s, correction suggested by Panith Adulsiriswad
       call periodic_particle_mlt6b(dns,mom,ppara,pperp,qpara,qperp)
       call partsm1(dns,cwwa)
       call partsm1(mom,cwwa)

       call partsm1(ppara,cwwa)
       call partsm1(pperp,cwwa)

       call partsm1(qpara,cwwa)
       call partsm1(qperp,cwwa)
!2024-12-21e

! calculate density (per volume)

!      do k = 1, lphi
!      do j = 1, lz
!      do i = 1, lr
!$omp parallel do private(vol)
      do i = 1, lrzphi
        vol = 1.0d0/(grr(i,1,1)*dr*dz*dphi)
        dns(i,1,1) = dns(i,1,1)*vol
        mom(i,1,1) = mom(i,1,1)*vol
        ppara(i,1,1) = ppara(i,1,1)*vol*mass1
        pperp(i,1,1) = pperp(i,1,1)*vol*babs(i,1,1)
        qpara(i,1,1) = qpara(i,1,1)*vol*mass1**2
        qperp(i,1,1) = qperp(i,1,1)*vol*babs(i,1,1)*mass1
      end do
!      end do
!      end do



! wall effect

      lr1 = lr - 1
      lz1 = lz - 1

      if(my_rank_r.eq.(mpi_proc_r-1))then
!      do k = 1, lphi
!      do j = 1, lz
!$omp parallel do
       do j = 1, lzphi
          dns(lr1,j,1) = dns(lr1,j,1) + dns(lr,j,1)
          mom(lr1,j,1) = mom(lr1,j,1) + mom(lr,j,1)
          dns(lr,j,1) = 0.0d0
          mom(lr,j,1) = 0.0d0

          ppara(lr1,j,1) = ppara(lr1,j,1) + ppara(lr,j,1)
          pperp(lr1,j,1) = pperp(lr1,j,1) + pperp(lr,j,1)
          ppara(lr,j,1) = 0.0d0
          pperp(lr,j,1) = 0.0d0

          qpara(lr1,j,1) = qpara(lr1,j,1) + qpara(lr,j,1)
          qperp(lr1,j,1) = qperp(lr1,j,1) + qperp(lr,j,1)
          qpara(lr,j,1) = 0.0d0
          qperp(lr,j,1) = 0.0d0
       end do
!      end do
      end if

      if(my_rank_r.eq.0)then
!      do k = 1, lphi
!      do j = 1, lz
!$omp parallel do
       do j = 1, lzphi
          dns(2,j,1) = dns(2,j,1) + dns(1,j,1)
          mom(2,j,1) = mom(2,j,1) + mom(1,j,1)
          dns(1,j,1) = 0.0d0
          mom(1,j,1) = 0.0d0

          ppara(2,j,1) = ppara(2,j,1) + ppara(1,j,1)
          pperp(2,j,1) = pperp(2,j,1) + pperp(1,j,1)
          ppara(1,j,1) = 0.0d0
          pperp(1,j,1) = 0.0d0

          qpara(2,j,1) = qpara(2,j,1) + qpara(1,j,1)
          qperp(2,j,1) = qperp(2,j,1) + qperp(1,j,1)
          qpara(1,j,1) = 0.0d0
          qperp(1,j,1) = 0.0d0
       end do
!      end do
      end if

      if(my_rank_z.eq.(mpi_proc_z-1))then
!$omp parallel do
       do k = 1, lphi
       do i = 1, lr
          dns(i,lz1,k) = dns(i,lz1,k) + dns(i,lz,k)
          mom(i,lz1,k) = mom(i,lz1,k) + mom(i,lz,k)
          dns(i,lz,k) = 0.0d0
          mom(i,lz,k) = 0.0d0

          ppara(i,lz1,k) = ppara(i,lz1,k) + ppara(i,lz,k)
          pperp(i,lz1,k) = pperp(i,lz1,k) + pperp(i,lz,k)
          ppara(i,lz,k) = 0.0d0
          pperp(i,lz,k) = 0.0d0

          qpara(i,lz1,k) = qpara(i,lz1,k) + qpara(i,lz,k)
          qperp(i,lz1,k) = qperp(i,lz1,k) + qperp(i,lz,k)
          qpara(i,lz,k) = 0.0d0
          qperp(i,lz,k) = 0.0d0
       end do
       end do
      end if 

      if(my_rank_z.eq.0)then
!$omp parallel do
       do k = 1, lphi
       do i = 1, lr
          dns(i,2,k) = dns(i,2,k) + dns(i,1,k)
          mom(i,2,k) = mom(i,2,k) + mom(i,1,k)
          dns(i,1,k) = 0.0d0
          mom(i,1,k) = 0.0d0

          ppara(i,2,k) = ppara(i,2,k) + ppara(i,1,k)
          pperp(i,2,k) = pperp(i,2,k) + pperp(i,1,k)
          ppara(i,1,k) = 0.0d0
          pperp(i,1,k) = 0.0d0

          qpara(i,2,k) = qpara(i,2,k) + qpara(i,1,k)
          qperp(i,2,k) = qperp(i,2,k) + qperp(i,1,k)
          qpara(i,1,k) = 0.0d0
          qperp(i,1,k) = 0.0d0
       end do
       end do
      end if 

!      if(.not.flag_stored)then
! take modes n_min <= n <=n_max
!      n_min = 1
!      n_max = 2
!      call lowpass2_mlt2(n_min,n_max,dns,mom)
!      call lowpass2_mlt2(n_min,n_max,ppara,pperp)
!      call lowpass2_mlt2(n_min,n_max,qpara,qperp)
!      end if

!2014-08-28s
!$omp parallel do
      do i = 1, lrzphi
        dns(i,1,1) = dns(i,1,1) + dns0(i,1,1)
        mom(i,1,1) = mom(i,1,1) + mom0(i,1,1)
        ppara(i,1,1) = ppara(i,1,1) + ppara0(i,1,1)
        pperp(i,1,1) = pperp(i,1,1) + pperp0(i,1,1)/babs0(i,1,1)*babs(i,1,1)
        qpara(i,1,1) = qpara(i,1,1) + qpara0(i,1,1)
        qperp(i,1,1) = qperp(i,1,1) + qperp0(i,1,1)/babs0(i,1,1)*babs(i,1,1)
      end do
!2014-08-28e

!            call wall_clock(d5)
!            t5 = d5 - d4
!            if(my_rank.eq.0)then
!            write(6,*)'density 5, t5=',t5
!            end if

end subroutine moments_gyro

end module xeon6900p
!----------------------------------------------------------------------------

!----------------------------------------------------------------------------
module amdgpu
contains
!--------------------------------------------------------------------
subroutine push(marker_num,amassp,ep &
               ,type,temp,valpha,deltav,clambda0,dclambda & !2012-06-17
               ,gc,dgc,v &
               ,cf_pphi,pphi_min,pphi_max &
               ,flp,ispecies,igyro)
! type=0: maxwellian, type=1: slowing down, type=2: beam
! type=3: beam with finite pitch angle width !2012-06-17
! igyro=0: w/o FLR, igyro=1: w/ FLR
! modified for NEC SX-Aurora TSUBASA 2020-07-10
!   gc_x(:,n) = 1:r, 2:z, 3:phi, 4:p, 5:mu, 6:weight, 7:active, 8:p_phi, 9:fff, 10:fnrml
!--------------------------------------------------------------------
      use mpiset
      use field
      use grid
!      use equi_sol, only:raxis
      use particle, only:nu_krook !2024-04-23
      implicit none

      integer::marker_num,type
      integer::ispecies,igyro
      real(8)::amassp,ep
      real(8)::gc(ngc2,marker_num)
      real(8)::dgc(ngc1,marker_num)
      real(8)::v(marker_num)
      real(8)::cf_pphi(0:lcfpphi),pphi_min,pphi_max !2016-02-04
      real(8)::temp,valpha,deltav
      real(8)::dr1,dz1,dphi1,ep1,amsp1
      integer::n,ia,ia1,ja,ja1,ka,ka1,i
      real(8)::ar,ar1,az,az1,aphi,aphi1
!      real(8)::aaa1,aaa2,aaa3,aaa4,aaa5,aaa6,aaa7,aaa8

!2015-06-04s
!      integer::ijk_a(3,marker_each)
!      real(8)::aaa(8,marker_each)
!2015-06-04e

      real(8)::flp(nflp,marker_num),detotal
      real(8)::b2,babse,babs0e,b21,bre,bze,bphie
      real(8)::b1,br1a,bz1a,bphi1a,b10,br10,bz10,bphi10,dbpar,dbpr,dbpz,dbpphi
      real(8)::rhopar,orbpr,orbpz,orbpphi,denom1
      real(8)::wvnmlr,wvnmlz,wvnmlphi
      real(8)::dgrdbr,dgrdbz,dgrdbphi,dppar,dppar1,dppar2 !2012-08-31
      real(8)::pvpar,w1nmlr,w1nmlz,w1nmlphi
      real(8)::psinrm,pactive,prof,dprofdpsi,energy0,pphi_n,rminor,bmax
      real(8)::dpphi_n,vpara,vpara0,dvpara,sigma
      real(8)::vlcpart,dvlcpart,vnrm,sd_factor
      real(8)::sqr_pi1
      real(8)::nudt,coef,vrat !2016-01-09
      real(8)::nu_krook_dt !2024-04-23

      integer::kr,kz,kphi
      real(8)::ma_mi_r
      integer::kfl_start,kfl

      real(8)::pt_factor,energy,clambda0,dclambda,clambda(marker_each) !2012-06-17
      real(8)::dwpsi(marker_each),dwenrc(marker_each) !2015-06-23
!      real(8)::psip(marker_each)
!      real(8)::dpsi_dr(marker_each),dpsi_dz(marker_each),dpsi_dphi(marker_each)
      real(8)::bphi1e(marker_each) !2016-01-09

      integer,parameter::nblkp=1024*1024 !NEC
!2021-01-11s
      integer::ijk_a(nblkp,3)
      real(8)::aaa(nblkp,8)
      integer::nn,nsta,nend,ivect
!2021-01-11e
      integer::j,k


! time derivative of particle position and velocity

      dr1 = 1.0d0/dr
      dz1 = 1.0d0/dz
      dphi1 = 1.0d0/dphi
      ep1 = 1.0d0/ep
      amsp1 = 1.0d0/amassp
      sqr_pi1=1.0d0/sqrt(pi)

      nu_krook_dt = nu_krook * dt !2024-04-23

      kr   = 1 - kr_offset(my_rank)
      kz   = 1 - kz_offset(my_rank)
      kphi = 1 + lphi_shd - kphi_offset(my_rank)
      ma_mi_r = major_r-minor_r

! igyro=0: w/o FLR, igyro/=0: flp(1:7,:) is given in subr. extract_em
      if(igyro.eq.0)then
        kfl_start = 1
      else
        kfl_start = nflp_gyro + 1
      end if

!      call ftrace_region_begin('push0')

      if(igyro.eq.0)then !w/o FLR
!$omp target teams distribute parallel do collapse(3) private(k,j,i)
        do k=1,lphi
        do j=1,lz
        do i=1,lr
        fld(1,i,j,k) = er(i,j,k)
        fld(2,i,j,k) = ez(i,j,k)
        fld(3,i,j,k) = ephi(i,j,k)
        fld(4,i,j,k) = epara(i,j,k)
        fld(5,i,j,k) = br(i,j,k) - br0(i,j,k)
        fld(6,i,j,k) = bz(i,j,k) - bz0(i,j,k)
        fld(7,i,j,k) = bphi(i,j,k) - bphi0(i,j,k)
        end do
        end do
        end do
      end if

!2025-09-03
!$omp target teams distribute parallel do collapse(3) private(k,j,i)
        do k=1,lphi
        do j=1,lz
        do i=1,lr
        fld(11,i,j,k) = gradbr(i,j,k)
        fld(12,i,j,k) = gradbz(i,j,k)
        fld(13,i,j,k) = gradbphi(i,j,k)
        fld(14,i,j,k) = curvbr(i,j,k)
        fld(15,i,j,k) = curvbz(i,j,k)
        fld(16,i,j,k) = curvbphi(i,j,k)
        end do
        end do
        end do

      if(ispecies.eq.0)then
!$omp target teams distribute parallel do collapse(3) private(k,j,i)
        do k=1,lphi
        do j=1,lz
        do i=1,lr
        fld(23,i,j,k)= dns_e0(i,j,k)
        fld(24,i,j,k)= dns_e0_r(i,j,k)
        fld(25,i,j,k)= dns_e0_z(i,j,k)
        fld(26,i,j,k)= dns_e0_phi(i,j,k)
        fld(27,i,j,k)= temp_e0(i,j,k)
        fld(28,i,j,k)= temp_e0_r(i,j,k)
        fld(29,i,j,k)= temp_e0_z(i,j,k)
        fld(30,i,j,k)= temp_e0_phi(i,j,k)
        end do
        end do
        end do
      else if(ispecies.eq.1)then
!$omp target teams distribute parallel do collapse(3) private(k,j,i)
        do k=1,lphi
        do j=1,lz
        do i=1,lr
        fld(23,i,j,k)= dns_i0(i,j,k)
        fld(24,i,j,k)= dns_i0_r(i,j,k)
        fld(25,i,j,k)= dns_i0_z(i,j,k)
        fld(26,i,j,k)= dns_i0_phi(i,j,k)
        fld(27,i,j,k)= temp_i0(i,j,k)
        fld(28,i,j,k)= temp_i0_r(i,j,k)
        fld(29,i,j,k)= temp_i0_z(i,j,k)
        fld(30,i,j,k)= temp_i0_phi(i,j,k)
        end do
        end do
        end do
      else if(ispecies.eq.2)then
!$omp target teams distribute parallel do collapse(3) private(k,j,i)
        do k=1,lphi
        do j=1,lz
        do i=1,lr
        fld(23,i,j,k)= dns_a0(i,j,k)
        fld(24,i,j,k)= dns_a0_r(i,j,k)
        fld(25,i,j,k)= dns_a0_z(i,j,k)
        fld(26,i,j,k)= dns_a0_phi(i,j,k)
        fld(27,i,j,k)= temp_a0(i,j,k)
        fld(28,i,j,k)= temp_a0_r(i,j,k)
        fld(29,i,j,k)= temp_a0_z(i,j,k)
        fld(30,i,j,k)= temp_a0_phi(i,j,k)
        end do
        end do
        end do
      end if

!      call ftrace_region_end('push0')

      do nn = 1, marker_num, nblkp
        nsta = nn
        nend = min(nn+nblkp-1,marker_num)

!      call ftrace_region_begin('push1')

!$omp target teams distribute parallel do private(n,ivect,ar1,ar,az1,az,aphi1,aphi)
      do n = nsta, nend
        ivect = n - nsta + 1

        ijk_a(ivect,1)=max(1,min(lr  -1,int((gc(1,n)-ma_mi_r)*dr1  ) + kr  ))
        ijk_a(ivect,2)=max(1,min(lz  -1,int(gc(2,n)          *dz1  ) + kz  ))
        ijk_a(ivect,3)=max(1,min(lphi-1,int(gc(3,n)          *dphi1) + kphi))

        ar1  = max(0.0d0, min(1.0d0, (gc(1,n) - ma_mi_r)*dr1 -  dble(ijk_a(ivect,1) - kr)  ) )
        ar   = 1.0d0 - ar1
        az1  = max(0.0d0, min(1.0d0, gc(2,n)*dz1 - dble(ijk_a(ivect,2) - kz) ) )
        az   = 1.0d0 - az1
        aphi1 = max(0.0d0, min(1.0d0, gc(3,n)*dphi1 - dble(ijk_a(ivect,3) - kphi) ) )
        aphi = 1.0d0 - aphi1

        aaa(ivect,1) = ar *az *aphi
        aaa(ivect,2) = ar1*az *aphi
        aaa(ivect,3) = ar *az1*aphi
        aaa(ivect,4) = ar1*az1*aphi
        aaa(ivect,5) = ar *az *aphi1
        aaa(ivect,6) = ar1*az *aphi1
        aaa(ivect,7) = ar *az1*aphi1
        aaa(ivect,8) = ar1*az1*aphi1
      end do


! fields at each particle position
! igyro=0: w/o FLR, igyro/=0: flp(1:7,:) is given in subr. extract_em

!      call ftrace_region_end('push1')
!      call ftrace_region_begin('push2')

!NEC$ outerloop_unroll(4)
!$omp target teams distribute parallel do collapse(2) private(kfl,n,ivect,ia,ja,ka)
      do n = nsta, nend

      do kfl = kfl_start, nflp
        ivect = n - nsta + 1

        ia=ijk_a(ivect,1)
        ja=ijk_a(ivect,2)
        ka=ijk_a(ivect,3)

        flp(kfl,n) = fld(kfl, ia, ja,  ka  )*aaa(ivect,1) + fld(kfl, ia+1,ja,  ka  )*aaa(ivect,2) &
                   + fld(kfl, ia, ja+1,ka  )*aaa(ivect,3) + fld(kfl, ia+1,ja+1,ka  )*aaa(ivect,4) &
                   + fld(kfl, ia, ja,  ka+1)*aaa(ivect,5) + fld(kfl, ia+1,ja,  ka+1)*aaa(ivect,6) &
                   + fld(kfl, ia, ja+1,ka+1)*aaa(ivect,7) + fld(kfl, ia+1,ja+1,ka+1)*aaa(ivect,8)
        end do

      end do

!      call ftrace_region_end('push2')
!      call ftrace_region_begin('push3')

!$omp  target teams distribute parallel do
!$omp& private(n,bre,bze,bphie,b2,babse,babs0e,b1,b21,br1a,bz1a,bphi1a,b10,br10,bz10,bphi10,dbpar,dbpr,dbpz,dbpphi)
!$omp& private(rhopar,denom1,orbpr,orbpz,orbpphi,wvnmlr,wvnmlz,wvnmlphi,dgrdbr,dgrdbz,dgrdbphi,dppar,dppar2)
!$omp& private(dppar1,pvpar,w1nmlr,w1nmlz,w1nmlphi,detotal,energy0)
      do n = nsta, nend

! flp(5:7,n): delta_br(z,phi)
        bre = flp(5,n) + flp(8,n)
        bze = flp(6,n) + flp(9,n)
        bphie = flp(7,n) + flp(10,n)

        b2  = bre**2 + bze**2 + bphie**2
        babse= max(eps_b, sqrt(b2) )
        babs0e= max(eps_b, sqrt(flp(8,n)**2 + flp(9,n)**2 + flp(10,n)**2) )
        b1 = 1.0d0/babse
        b21= 1.0d0/b2
        br1a = bre*b1
        bz1a = bze*b1
        bphi1a = bphie*b1
        bphi1e(n) = bphi1a !2016-02-4

        b10 = 1.0d0/babs0e
        br10 = flp(8,n)*b10
        bz10 = flp(9,n)*b10
        bphi10 = flp(10,n)*b10

        dbpar = br1a*br10 + bz1a*bz10 + bphi1a*bphi10
        dbpr  = br1a  - br10*dbpar
        dbpz  = bz1a  - bz10*dbpar
        dbpphi= bphi1a- bphi10*dbpar

! guiding center motion

        rhopar = gc(4,n)*ep1*b1

        denom1 = 1.0d0/(1.0d0 + rhopar*(br1a*flp(14,n) + bz1a*flp(15,n) + bphi1a*flp(16,n)))

        orbpr = (br1a + rhopar*flp(14,n))*denom1
        orbpz = (bz1a + rhopar*flp(15,n))*denom1
        orbpphi = (bphi1a + rhopar*flp(16,n))*denom1

! e x b drift

        wvnmlr = (flp(3,n)*bze-flp(2,n)*bphie)*b21*denom1
        wvnmlz = (flp(1,n)*bphie-flp(3,n)*bre)*b21*denom1
        wvnmlphi = (flp(2,n)*bre-flp(1,n)*bze)*b21*denom1

! grad-b drift

        dgrdbr = gc(5,n)*(bphie*flp(12,n) - bze*flp(13,n))*ep1*b21*denom1
        dgrdbz = gc(5,n)*(bre*flp(13,n) - bphie*flp(11,n))*ep1*b21*denom1
        dgrdbphi = gc(5,n)*(bze*flp(11,n) - bre*flp(12,n))*ep1*b21*denom1

! mirror force

        dppar =-gc(5,n)*( flp(11,n)*orbpr &
                      + flp(12,n)*orbpz &
                      + flp(13,n)*orbpphi &
                      )*dt

!2012-08-31
        dppar2=-gc(5,n)*( flp(11,n)*dbpr*denom1 &
                      + flp(12,n)*dbpz*denom1 &
                      + flp(13,n)*orbpphi &
                      )*dt
!2012-08-31 end

! aceeleration due to electric field and curvature drift

!2012-07-07
        dppar1 = ep*((flp(1,n)*flp(14,n) + flp(2,n)*flp(15,n) + flp(3,n)*flp(16,n) &
                      )*rhopar &
                    + flp(4,n) &
                     )*dt*denom1
!2012-07-07 end

! total drift velocity

        pvpar = gc(4,n)*amsp1
        dgc(1,n) = dt*(pvpar*orbpr + wvnmlr + dgrdbr)*gc(7,n)
        dgc(2,n) = dt*(pvpar*orbpz + wvnmlz + dgrdbz)*gc(7,n)
        dgc(3,n) = dt*(pvpar*orbpphi + wvnmlphi + dgrdbphi)/gc(1,n)*gc(7,n)
        dgc(4,n) =(dppar + dppar1)*gc(7,n)

! temporal evolution of weight of high-energy ion particle

        w1nmlr = wvnmlr + pvpar*dbpr*denom1
        w1nmlz = wvnmlz + pvpar*dbpz*denom1
        w1nmlphi = wvnmlphi + pvpar*dbpphi*denom1

        dgc(8,n) =(ep*(w1nmlr*flp(18,n) + w1nmlz*flp(19,n) + w1nmlphi*flp(20,n))*dt &
                  + bphi1a*(w1nmlr*gc(4,n)*dt + gc(1,n)*(dppar1 + dppar2) ) & !2012-08-31
                  )*gc(7,n)

        detotal =(  ep*(flp(1,n)*dgrdbr + flp(2,n)*dgrdbz + flp(3,n)*dgrdbphi)*dt &
                + dppar1*pvpar )*gc(7,n) &
! the following term considers 'mu* v * grad(B0 - B)', 2025-04-04
                + gc(5,n)*(dgc(1,n)*(flp(31,n)-flp(11,n) ) &
                          +dgc(2,n)*(flp(32,n)-flp(12,n) ) &
                          +dgc(3,n)*(flp(33,n)-flp(13,n) )*gc(1,n) &
                          )
!2025-04-27s
        energy0 = max(1.0d-30, 0.50d0*amassp*pvpar**2 + gc(5,n)*babs0e)
        clambda(n) = gc(5,n)*b0/energy0
        v(n) = sqrt(2.0d0*energy0*amsp1)

!        energy = max(1.0d-30, 0.50d0*amassp*pvpar**2 + gc(5,n)*babse)
!        clambda(n) = gc(5,n)*b0/energy
!        v(n) = sqrt(2.0d0*energy*amsp1)
!2025-04-27e

! weight evolution : weight = f - f0
! d weight /dt = - d f0 /dt
! f0 = f_nrml*prof(psinrm)/(v**3 + flp(21,n)**3)*0.50*erfc((v(n)-valpha)/deltav)

!        rminor = sqrt( (gc(1,n)-raxis)**2 &
!                     + ((gc(2,n)-0.50d0*zleng) )**2 &
!                     )
!        bmax = b0*raxis/(raxis-rminor)
!        energy0 = 0.50d0*amsp1*gc(4,n)**2 + gc(5,n)*babs0e
!        sigma = 0.50d0*(1.0d0 + sign(1.0d0, energy0-gc(5,n)*bmax) ) &
!                      *sign(1.0d0,gc(4,n) )
!        vpara0 = sqrt(2.0d0*(energy0-gc(5,n)*bmin)*amsp1)
!        vpara = vpara0*sigma
!        dvpara = sigma/vpara0*gc(4,n)*dppar1*amsp1**2

!        pphi_n = gc(8,n) - amassp*raxis*vpara
!        dpphi_n= dgc(8,n) - amassp*raxis*dvpara

!        prof = exp(pphi_n/(ep*psimax*0.37d0) )
!        dprofdpsi = prof/(ep*psimax*0.37d0)

!        dwpsi(n) = dpphi_n*dprofdpsi*gc(10,n)
!        dwenrc(n)= detotal/(amassp*v(n) )*prof*gc(10,n)

!2016-08-05s
        dwpsi(n) = ( w1nmlr*flp(24,n) + w1nmlz*flp(25,n) + w1nmlphi*flp(26,n) &
                    +(w1nmlr*flp(28,n) + w1nmlz*flp(29,n) + w1nmlphi*flp(30,n) ) &
                    *0.50d0*(amassp*v(n)**2/flp(27,n) - 3.0d0)*flp(23,n)/flp(27,n) &
                   )*dt*gc(7,n)*gc(10,n)

!        dwpsi(n) = (w1nmlr*flp(24,n) + w1nmlz*flp(25,n) + w1nmlphi*flp(26,n))*dt & !2016-01-09
!                   *gc(7,n)*gc(10,n)
!2016-08-05e

        dwenrc(n)= detotal/(amassp*v(n) )*flp(23,n)*gc(10,n) !2016-01-09

      end do

!      call ftrace_region_end('push3')
!      call ftrace_region_begin('push4')

      if(type.eq.0)then

!$omp target teams distribute parallel do private(n,vlcpart,dvlcpart)
        do n = nsta, nend
!2016-08-05s
          vlcpart = exp(-0.50d0*amassp*v(n)**2/flp(27,n) )*flp(27,n)**(-1.5d0) !2013-07-17
          dvlcpart = -amassp*v(n)/flp(27,n)*vlcpart 
!          vlcpart = exp(-0.50d0*amassp*v(n)**2/temp)
!          dvlcpart = -amassp*v(n)/temp*vlcpart 
!2016-08-05e
          dwpsi(n) = dwpsi(n)*vlcpart
          dwenrc(n) = dwenrc(n)*dvlcpart
          dgc(6,n) = - dwpsi(n) - dwenrc(n) &
                     - gc(6,n)*nu_krook_dt & !Krook operator at normalized psi > psi_edge, flp(17)=psi
                     *(0.50d0 + sign(0.50d0, 1.0d0 - flp(17,n)/psimax - psi_edge) ) !2024-04-23
        end do

      else if(type.eq.1.or.type.eq.2)then

!$omp target teams distribute parallel do private(n,vnrm,sd_factor,vlcpart,dvlcpart)
        do n = nsta, nend
          vnrm =(v(n)-valpha)/deltav
          sd_factor = 1.0d0/(v(n)**3 + flp(21,n)**3)
          vlcpart = 0.50d0*derfc(vnrm)*sd_factor
          dvlcpart= -( 3.0d0*v(n)**2*vlcpart &
                     + exp(-vnrm**2)*sqr_pi1/deltav &
                     )*sd_factor

          dwpsi(n) = dwpsi(n)*vlcpart
          dwenrc(n) = dwenrc(n)*dvlcpart
          dgc(6,n) = - dwpsi(n) - dwenrc(n) &
                     - gc(6,n)*nu_krook_dt & !Krook operator at normalized psi > psi_edge, flp(17)=psi
                     *(0.50d0 + sign(0.50d0, 1.0d0 - flp(17,n)/psimax - psi_edge) ) !2024-04-23
        end do

!2012-06-17
      else if(type.eq.3)then

!$omp target teams distribute parallel do private(n,vnrm,sd_factor,pt_factor,vlcpart,dvlcpart)
        do n = nsta, nend
          vnrm =(v(n)-valpha)/deltav
          sd_factor = 1.0d0/(v(n)**3 + flp(21,n)**3)
          pt_factor = exp(-(clambda(n)-clambda0)**2/dclambda**2)

          vlcpart = 0.50d0*derfc(vnrm)*sd_factor*pt_factor
          dvlcpart= -( 3.0d0*v(n)**2*vlcpart &
                     + exp(-vnrm**2)*sqr_pi1/deltav*pt_factor &
                     )*sd_factor &
                  + 4.0d0*vlcpart*clambda(n)*(clambda(n)-clambda0) &
                         /(v(n)*dclambda**2)

          dwpsi(n) = dwpsi(n)*vlcpart
          dwenrc(n) = dwenrc(n)*dvlcpart
          dgc(6,n) = - dwpsi(n) - dwenrc(n) &
                     - gc(6,n)*nu_krook_dt & !Krook operator at normalized psi > psi_edge, flp(17)=psi
                     *(0.50d0 + sign(0.50d0, 1.0d0 - flp(17,n)/psimax - psi_edge) ) !2024-04-23
        end do
!2012-06-17 end

      end if

! slowing down
      if(type.eq.-5)then
!$omp target teams distribute parallel do private(n,nudt,vrat,coef)
        do n = nsta, nend
            dgc(6,n) = 0.0d0 !full-f
            nudt = flp(22,n)*dt*gc(7,n)
            vrat = flp(21,n)/v(n)
            coef =(1.0d0 + vrat**3)*nudt
            dgc(10,n) =-3.0d0*nudt*gc(10,n)
            dgc(5,n)= -2.0d0*gc(5,n)*coef
            dgc(4,n) = dgc(4,n)- gc(4,n)*coef
            dgc(8,n) = dgc(8,n) - gc(4,n)*coef*gc(1,n)*bphi1e(n)
        end do
      end if

!      call ftrace_region_end('push4')

      end do !nn

end subroutine push
!--------------------------------------------------------------------
subroutine density(marker_num,mass,gc &
                  ,dns,mom,ppara,pperp &
                  ,dns0,mom0,ppara0,pperp0)
! gc_x(:,n) = 1:r, 2:z, 3:phi, 4:p, 5:mu, 6:weight, 7:active, 8:p_phi, 9:fff, 10:fnrml
!--------------------------------------------------------------------
      use mpiset
      use grid
      use field, only:babs,babs0
!      use mod_timer ! NEC TIMER
      implicit none

      integer::marker_num
      real(8)::gc(ngc2,marker_num)
      real(8)::mass
      real(8)::dns(lr,lz,lphi),mom(lr,lz,lphi)
      real(8)::ppara(lr,lz,lphi),pperp(lr,lz,lphi)
      real(8)::dns0(lr,lz,lphi),mom0(lr,lz,lphi)
      real(8)::ppara0(lr,lz,lphi),pperp0(lr,lz,lphi)
      integer::i,j,k,l,n,ia,ia1,ja,ja1,ka,ka1,lr1,lz1
      real(8)::ar,ar1,az,az1,aphi,aphi1
      real(8)::aaa1,aaa2,aaa3,aaa4,aaa5,aaa6,aaa7,aaa8
      real(8)::vol,cwwa,cwwc
      real(8)::dr1,dz1,dphi1,ma_mi_r
      integer::kr,kz,kphi
      real(8)::p0,p1,p2,mu1,mass1
      integer::n_min,n_max
      integer::ib
      integer,parameter::nwork=32
      !real(8),dimension(lr,lz,lphi,nwork)::dns_,mom_,ppara_,pperp_
      real(8),dimension(nwork,lr,lz,lphi)::dns_,mom_,ppara_,pperp_
!      real(8)::d0,d1,d2,d3,d4,d5,t1,t2,t3,t4,t5

      mass1 = 1.0d0/mass

      dr1 = 1.0d0/dr
      dz1 = 1.0d0/dz
      dphi1 = 1.0d0/dphi

      kr   = 1 - kr_offset(my_rank)
      kz   = 1 - kz_offset(my_rank)
      kphi = 1 + lphi_shd - kphi_offset(my_rank)
!      kphi = 3 - kphi_offset(my_rank)
      ma_mi_r = major_r-minor_r

!      call start_timer('  density:loop02') ! NEC TIMER
!$omp target teams distribute parallel do private(i)
        do i = 1, lrzphi
          dns(i,1,1) = 0.0d0
          mom(i,1,1) = 0.0d0
          ppara(i,1,1) = 0.0d0
          pperp(i,1,1) = 0.0d0
        end do
!      call stop_timer('  density:loop02') ! NEC TIMER

!      call start_timer('  density:loop03before') ! NEC TIMER
!$omp target teams distribute parallel do private(i)
        do i = 1, lrzphi*nwork
          dns_(i,1,1,1) = 0.0d0
          mom_(i,1,1,1) = 0.0d0
          ppara_(i,1,1,1) = 0.0d0
          pperp_(i,1,1,1) = 0.0d0
        end do
!      call stop_timer('  density:loop03before') ! NEC TIMER

!      call start_timer('  density:loop03') ! NEC TIMER
!$omp  target teams distribute parallel do
!$omp& private(n,ia,ja,ka,ar1,ar,az1,az,aphi1,aphi,aaa1,aaa2,aaa3,aaa4,aaa5,aaa6,aaa7,aaa8)
!$omp& private(p0,p1,p2,mu1,ib)
       do n = 1, marker_num
        ib = mod(n-1,nwork)+1

! gc_x(:,n) = 1:r, 2:z, 3:phi, 4:p, 5:mu, 6:weight, 7:active, 8:p_phi, 9:fff, 10:fnrml
        ia=max(1,min(lr  -1,int((gc(1,n)-ma_mi_r)*dr1  ) + kr  ))
        ja=max(1,min(lz  -1,int(gc(2,n)          *dz1  ) + kz  ))
        ka=max(1,min(lphi-1,int(gc(3,n)        *dphi1) + kphi))

        ar1 = (gc(1,n)-grr(ia,ja,ka ) ) *dr1
        ar  = 1.0d0 - ar1
        az1 = (gc(2,n)-gzz(ia,ja,ka ) ) *dz1
        az  = 1.0d0 - az1
        aphi1 = (gc(3,n)-gphi(ia,ja,ka ) ) *dphi1
        aphi  = 1.0d0 - aphi1

        aaa1 = ar *az *aphi
        aaa2 = ar1*az *aphi
        aaa3 = ar *az1*aphi
        aaa4 = ar1*az1*aphi
        aaa5 = ar *az *aphi1
        aaa6 = ar1*az *aphi1
        aaa7 = ar *az1*aphi1
        aaa8 = ar1*az1*aphi1

        p0 =  gc(6,n)*gc(7,n)
        p1 =  gc(4,n)*p0
        p2 =  gc(4,n)*p1
        mu1 = gc(5,n)*p0

!$omp atomic update
        dns_(ib,ia,  ja,  ka  ) = dns_(ib,ia,  ja,  ka  ) + aaa1*p0
!$omp atomic update
        dns_(ib,ia+1,ja,  ka  ) = dns_(ib,ia+1,ja,  ka  ) + aaa2*p0
!$omp atomic update
        dns_(ib,ia,  ja+1,ka  ) = dns_(ib,ia,  ja+1,ka  ) + aaa3*p0
!$omp atomic update
        dns_(ib,ia+1,ja+1,ka  ) = dns_(ib,ia+1,ja+1,ka  ) + aaa4*p0
!$omp atomic update
        dns_(ib,ia,  ja,  ka+1) = dns_(ib,ia,  ja,  ka+1) + aaa5*p0
!$omp atomic update
        dns_(ib,ia+1,ja,  ka+1) = dns_(ib,ia+1,ja,  ka+1) + aaa6*p0
!$omp atomic update
        dns_(ib,ia,  ja+1,ka+1) = dns_(ib,ia,  ja+1,ka+1) + aaa7*p0
!$omp atomic update
        dns_(ib,ia+1,ja+1,ka+1) = dns_(ib,ia+1,ja+1,ka+1) + aaa8*p0

!$omp atomic update
        mom_(ib,ia,  ja,  ka  ) = mom_(ib,ia,  ja,  ka  ) + aaa1*p1
!$omp atomic update
        mom_(ib,ia+1,ja,  ka  ) = mom_(ib,ia+1,ja,  ka  ) + aaa2*p1
!$omp atomic update
        mom_(ib,ia,  ja+1,ka  ) = mom_(ib,ia,  ja+1,ka  ) + aaa3*p1
!$omp atomic update
        mom_(ib,ia+1,ja+1,ka  ) = mom_(ib,ia+1,ja+1,ka  ) + aaa4*p1
!$omp atomic update
        mom_(ib,ia,  ja,  ka+1) = mom_(ib,ia,  ja,  ka+1) + aaa5*p1
!$omp atomic update
        mom_(ib,ia+1,ja,  ka+1) = mom_(ib,ia+1,ja,  ka+1) + aaa6*p1
!$omp atomic update
        mom_(ib,ia,  ja+1,ka+1) = mom_(ib,ia,  ja+1,ka+1) + aaa7*p1
!$omp atomic update
        mom_(ib,ia+1,ja+1,ka+1) = mom_(ib,ia+1,ja+1,ka+1) + aaa8*p1

!$omp atomic update
        ppara_(ib,ia,  ja,  ka  ) = ppara_(ib,ia,  ja,  ka  ) + aaa1*p2
!$omp atomic update
        ppara_(ib,ia+1,ja,  ka  ) = ppara_(ib,ia+1,ja,  ka  ) + aaa2*p2
!$omp atomic update
        ppara_(ib,ia,  ja+1,ka  ) = ppara_(ib,ia,  ja+1,ka  ) + aaa3*p2
!$omp atomic update
        ppara_(ib,ia+1,ja+1,ka  ) = ppara_(ib,ia+1,ja+1,ka  ) + aaa4*p2
!$omp atomic update
        ppara_(ib,ia,  ja,  ka+1) = ppara_(ib,ia,  ja,  ka+1) + aaa5*p2
!$omp atomic update
        ppara_(ib,ia+1,ja,  ka+1) = ppara_(ib,ia+1,ja,  ka+1) + aaa6*p2
!$omp atomic update
        ppara_(ib,ia,  ja+1,ka+1) = ppara_(ib,ia,  ja+1,ka+1) + aaa7*p2
!$omp atomic update
        ppara_(ib,ia+1,ja+1,ka+1) = ppara_(ib,ia+1,ja+1,ka+1) + aaa8*p2

!$omp atomic update
        pperp_(ib,ia,  ja,  ka  ) = pperp_(ib,ia,  ja,  ka  ) + aaa1*mu1
!$omp atomic update
        pperp_(ib,ia+1,ja,  ka  ) = pperp_(ib,ia+1,ja,  ka  ) + aaa2*mu1
!$omp atomic update
        pperp_(ib,ia,  ja+1,ka  ) = pperp_(ib,ia,  ja+1,ka  ) + aaa3*mu1
!$omp atomic update
        pperp_(ib,ia+1,ja+1,ka  ) = pperp_(ib,ia+1,ja+1,ka  ) + aaa4*mu1
!$omp atomic update
        pperp_(ib,ia,  ja,  ka+1) = pperp_(ib,ia,  ja,  ka+1) + aaa5*mu1
!$omp atomic update
        pperp_(ib,ia+1,ja,  ka+1) = pperp_(ib,ia+1,ja,  ka+1) + aaa6*mu1
!$omp atomic update
        pperp_(ib,ia,  ja+1,ka+1) = pperp_(ib,ia,  ja+1,ka+1) + aaa7*mu1
!$omp atomic update
        pperp_(ib,ia+1,ja+1,ka+1) = pperp_(ib,ia+1,ja+1,ka+1) + aaa8*mu1
       end do
!      call stop_timer('  density:loop03') ! NEC TIMER

!      call start_timer('  density:loop03after') ! NEC TIMER
!$omp target teams distribute parallel do private(i,j)
        do i = 1, lrzphi
        do j = 1, nwork
          dns(i,1,1) = dns(i,1,1) + dns_(j,i,1,1)
          mom(i,1,1) = mom(i,1,1) + mom_(j,i,1,1)
          ppara(i,1,1) = ppara(i,1,1) + ppara_(j,i,1,1)
          pperp(i,1,1) = pperp(i,1,1) + pperp_(j,i,1,1)
        end do
        end do
!      call stop_timer('  density:loop03after') ! NEC TIMER
! smoothing

       cwwa = 0.5d0
!       cwwc =-1.d0/6.d0 

!2024-12-21s, correction suggested by Panith Adulsiriswad
!      call start_timer('  density:com01') ! NEC TIMER
       call periodic_particle_mlt4b(dns,mom,ppara,pperp)
       call partsm1(dns,cwwa)
       call partsm1(mom,cwwa)

       call partsm1(ppara,cwwa)
       call partsm1(pperp,cwwa)
!      call stop_timer('  density:com01') ! NEC TIMER
!2024-12-21e

! calculate density (per volume)

!      do k = 1, lphi
!      do j = 1, lz
!      do i = 1, lr
!      call start_timer('  density:loop05') ! NEC TIMER
!$omp target teams distribute parallel do private(i,vol) !corrected by R. Seki
      do i = 1, lrzphi
        vol = 1.0d0/(grr(i,1,1)*dr*dz*dphi)
        dns(i,1,1) = dns(i,1,1)*vol
        mom(i,1,1) = mom(i,1,1)*vol
        ppara(i,1,1) = ppara(i,1,1)*vol/mass
        pperp(i,1,1) = pperp(i,1,1)*vol*babs(i,1,1)
      end do
!      call stop_timer('  density:loop05') ! NEC TIMER
!      end do
!      end do


! wall effect

      lr1 = lr - 1
      lz1 = lz - 1

      if(my_rank_r.eq.(mpi_proc_r-1))then
!      do k = 1, lphi
!      do j = 1, lz
!      call start_timer('  density:loop06') ! NEC TIMER
!$omp target teams distribute parallel do private(j)
       do j = 1, lzphi
          dns(lr1,j,1) = dns(lr1,j,1) + dns(lr,j,1)
          mom(lr1,j,1) = mom(lr1,j,1) + mom(lr,j,1)
          dns(lr,j,1) = 0.0d0
          mom(lr,j,1) = 0.0d0

          ppara(lr1,j,1) = ppara(lr1,j,1) + ppara(lr,j,1)
          pperp(lr1,j,1) = pperp(lr1,j,1) + pperp(lr,j,1)
          ppara(lr,j,1) = 0.0d0
          pperp(lr,j,1) = 0.0d0
       end do
!      call stop_timer('  density:loop06') ! NEC TIMER
!      end do
      end if

      if(my_rank_r.eq.0)then
!      do k = 1, lphi
!      do j = 1, lz
!      call start_timer('  density:loop07') ! NEC TIMER
!$omp target teams distribute parallel do private(j)
       do j = 1, lzphi
          dns(2,j,1) = dns(2,j,1) + dns(1,j,1)
          mom(2,j,1) = mom(2,j,1) + mom(1,j,1)
          dns(1,j,1) = 0.0d0
          mom(1,j,1) = 0.0d0

          ppara(2,j,1) = ppara(2,j,1) + ppara(1,j,1)
          pperp(2,j,1) = pperp(2,j,1) + pperp(1,j,1)
          ppara(1,j,1) = 0.0d0
          pperp(1,j,1) = 0.0d0
       end do
!      call stop_timer('  density:loop07') ! NEC TIMER
!      end do
      end if

      if(my_rank_z.eq.(mpi_proc_z-1))then
!      call start_timer('  density:loop08') ! NEC TIMER
!$omp target teams distribute parallel do private(k,i) collapse(2)
       do k = 1, lphi
       do i = 1, lr
          dns(i,lz1,k) = dns(i,lz1,k) + dns(i,lz,k)
          mom(i,lz1,k) = mom(i,lz1,k) + mom(i,lz,k)
          dns(i,lz,k) = 0.0d0
          mom(i,lz,k) = 0.0d0

          ppara(i,lz1,k) = ppara(i,lz1,k) + ppara(i,lz,k)
          pperp(i,lz1,k) = pperp(i,lz1,k) + pperp(i,lz,k)
          ppara(i,lz,k) = 0.0d0
          pperp(i,lz,k) = 0.0d0
       end do
       end do
!      call stop_timer('  density:loop08') ! NEC TIMER
      end if 

      if(my_rank_z.eq.0)then
!      call start_timer('  density:loop09') ! NEC TIMER
!$omp target teams distribute parallel do private(k,i) collapse(2)
       do k = 1, lphi
       do i = 1, lr
          dns(i,2,k) = dns(i,2,k) + dns(i,1,k)
          mom(i,2,k) = mom(i,2,k) + mom(i,1,k)
          dns(i,1,k) = 0.0d0
          mom(i,1,k) = 0.0d0

          ppara(i,2,k) = ppara(i,2,k) + ppara(i,1,k)
          pperp(i,2,k) = pperp(i,2,k) + pperp(i,1,k)
          ppara(i,1,k) = 0.0d0
          pperp(i,1,k) = 0.0d0
       end do
       end do
!      call stop_timer('  density:loop09') ! NEC TIMER
      end if 

! take only n=1 modes
!      n = 1
!      call n1(n,ppara)
!      call n1(n,pperp)

! take modes n_min <= n <=n_max
!      n_min = 1
!      n_max = 2 
!      call start_timer('  density:com02') ! NEC TIMER
!      call lowpass2_mlt2(n_min,n_max,dns,mom)
!      call lowpass2_mlt2(n_min,n_max,ppara,pperp)
!      call stop_timer('  density:com02') ! NEC TIMER

!      do k = 1, lphi
!      do j = 1, lz
!      do i = 1, lr
!      call start_timer('  density:loop10') ! NEC TIMER
!$omp target teams distribute parallel do private(i)
      do i = 1, lrzphi
        dns(i,1,1) = dns(i,1,1) + dns0(i,1,1)
        mom(i,1,1) = mom(i,1,1) + mom0(i,1,1)
        ppara(i,1,1) = ppara(i,1,1) + ppara0(i,1,1)
        pperp(i,1,1) = pperp(i,1,1) + pperp0(i,1,1)/babs0(i,1,1)*babs(i,1,1)
      end do
!      call stop_timer('  density:loop10') ! NEC TIMER
!      end do
!      end do

end subroutine density
!--------------------------------------------------------------------
subroutine moments(marker_num,mass,gc &
                  ,dns,mom,ppara,pperp,qpara,qperp &
                  ,dns0,mom0,ppara0,pperp0,qpara0,qperp0)
! gc_x(:,n) = 1:r, 2:z, 3:phi, 4:p, 5:mu, 6:weight, 7:active, 8:p_phi, 9:fff, 10:fnrml
!--------------------------------------------------------------------
      use mpiset
      use grid
      use field, only:babs,babs0
!      use mod_timer ! NEC TIMER
      implicit none

      integer::marker_num,nvec,marker_num_gyro
      real(8)::gc(ngc2,marker_num)
      real(8)::mass
      real(8)::dns(lr,lz,lphi),mom(lr,lz,lphi)
      real(8)::ppara(lr,lz,lphi),pperp(lr,lz,lphi)
      real(8)::qpara(lr,lz,lphi),qperp(lr,lz,lphi)
      real(8)::dns0(lr,lz,lphi),mom0(lr,lz,lphi)
      real(8)::ppara0(lr,lz,lphi),pperp0(lr,lz,lphi)
      real(8)::qpara0(lr,lz,lphi),qperp0(lr,lz,lphi)
      integer::i,j,k,l,m,n,ia,ia1,ja,ja1,ka,ka1,lr1,lz1
      real(8)::ar,ar1,az,az1,aphi,aphi1
      real(8)::aaa1,aaa2,aaa3,aaa4,aaa5,aaa6,aaa7,aaa8
      real(8)::vol,cwwa,cwwc
      real(8)::dr1,dz1,dphi1,ma_mi_r
      integer::kr,kz,kphi
      real(8)::p0,p1,p2,mu1,p3,p1mu,mass1 !2015-09-21
      integer::n_min,n_max !2013-05-22
      integer::ib
      integer,parameter::nwork=32
      real(8),dimension(nwork,lr,lz,lphi)::dns_,mom_,ppara_,pperp_,qpara_,qperp_
!      real(8)::d0,d1,d2,d3,d4,d5,t1,t2,t3,t4,t5

      mass1 = 1.0d0/mass

      dr1 = 1.0d0/dr
      dz1 = 1.0d0/dz
      dphi1 = 1.0d0/dphi

      kr   = 1 - kr_offset(my_rank)
      kz   = 1 - kz_offset(my_rank)
      kphi = 1 + lphi_shd - kphi_offset(my_rank)
!      kphi = 3 - kphi_offset(my_rank)
      ma_mi_r = major_r-minor_r

!      call start_timer('  moments:loop02') ! NEC TIMER
!$omp target teams distribute parallel do private(i)
        do i = 1, lrzphi
          dns(i,1,1) = 0.0d0
          mom(i,1,1) = 0.0d0
          ppara(i,1,1) = 0.0d0
          pperp(i,1,1) = 0.0d0
          qpara(i,1,1) = 0.0d0
          qperp(i,1,1) = 0.0d0
        end do
!      call stop_timer('  moments:loop02') ! NEC TIMER

!      call start_timer('  moments:loop03before') ! NEC TIMER
!$omp target teams distribute parallel do private(i)
        do i = 1, lrzphi*nwork
          dns_(i,1,1,1) = 0.0d0
          mom_(i,1,1,1) = 0.0d0
          ppara_(i,1,1,1) = 0.0d0
          pperp_(i,1,1,1) = 0.0d0
          qpara_(i,1,1,1) = 0.0d0
          qperp_(i,1,1,1) = 0.0d0
        end do
!      call stop_timer('  moments:loop03before') ! NEC TIMER

!      call start_timer('  moments:loop03') ! NEC TIMER
!$omp  target teams distribute parallel do
!$omp& private(n,ia,ja,ka,ar1,ar,az1,az,aphi1,aphi,aaa1,aaa2,aaa3,aaa4,aaa5,aaa6,aaa7,aaa8)
!$omp& private(p0,p1,p2,mu1,p3,p1mu,ib)
       do n = 1, marker_num_gyro
        ib = mod(n-1,nwork)+1

! gc_x(:,n) = 1:r, 2:z, 3:phi, 4:p, 5:mu, 6:weight, 7:active, 8:p_phi, 9:fff, 10:fnrml
        ia=max(1,min(lr  -1,int((gc(1,n)-ma_mi_r)*dr1  ) + kr  ))
        ja=max(1,min(lz  -1,int(gc(2,n)          *dz1  ) + kz  ))
        ka=max(1,min(lphi-1,int(gc(3,n)        *dphi1) + kphi))

        ar1 = (gc(1,n)-grr(ia,ja,ka ) ) *dr1
        ar  = 1.0d0 - ar1
        az1 = (gc(2,n)-gzz(ia,ja,ka ) ) *dz1
        az  = 1.0d0 - az1
        aphi1 = (gc(3,n)-gphi(ia,ja,ka ) ) *dphi1
        aphi  = 1.0d0 - aphi1

        aaa1 = ar *az *aphi
        aaa2 = ar1*az *aphi
        aaa3 = ar *az1*aphi
        aaa4 = ar1*az1*aphi
        aaa5 = ar *az *aphi1
        aaa6 = ar1*az *aphi1
        aaa7 = ar *az1*aphi1
        aaa8 = ar1*az1*aphi1

        p0 =  gc(6,n)*gc(7,n)
        p1 =  gc(4,n)*p0
        p2 =  gc(4,n)*p1
        mu1 = gc(5,n)*p0
        p3  = gc(4,n)*p2
        p1mu= gc(4,n)*mu1

!$omp atomic update
        dns_(ib,ia,  ja,  ka  ) = dns_(ib,ia,  ja,  ka  ) + aaa1*p0
!$omp atomic update
        dns_(ib,ia+1,ja,  ka  ) = dns_(ib,ia+1,ja,  ka  ) + aaa2*p0
!$omp atomic update
        dns_(ib,ia,  ja+1,ka  ) = dns_(ib,ia,  ja+1,ka  ) + aaa3*p0
!$omp atomic update
        dns_(ib,ia+1,ja+1,ka  ) = dns_(ib,ia+1,ja+1,ka  ) + aaa4*p0
!$omp atomic update
        dns_(ib,ia,  ja,  ka+1) = dns_(ib,ia,  ja,  ka+1) + aaa5*p0
!$omp atomic update
        dns_(ib,ia+1,ja,  ka+1) = dns_(ib,ia+1,ja,  ka+1) + aaa6*p0
!$omp atomic update
        dns_(ib,ia,  ja+1,ka+1) = dns_(ib,ia,  ja+1,ka+1) + aaa7*p0
!$omp atomic update
        dns_(ib,ia+1,ja+1,ka+1) = dns_(ib,ia+1,ja+1,ka+1) + aaa8*p0

!$omp atomic update
        mom_(ib,ia,  ja,  ka  ) = mom_(ib,ia,  ja,  ka  ) + aaa1*p1
!$omp atomic update
        mom_(ib,ia+1,ja,  ka  ) = mom_(ib,ia+1,ja,  ka  ) + aaa2*p1
!$omp atomic update
        mom_(ib,ia,  ja+1,ka  ) = mom_(ib,ia,  ja+1,ka  ) + aaa3*p1
!$omp atomic update
        mom_(ib,ia+1,ja+1,ka  ) = mom_(ib,ia+1,ja+1,ka  ) + aaa4*p1
!$omp atomic update
        mom_(ib,ia,  ja,  ka+1) = mom_(ib,ia,  ja,  ka+1) + aaa5*p1
!$omp atomic update
        mom_(ib,ia+1,ja,  ka+1) = mom_(ib,ia+1,ja,  ka+1) + aaa6*p1
!$omp atomic update
        mom_(ib,ia,  ja+1,ka+1) = mom_(ib,ia,  ja+1,ka+1) + aaa7*p1
!$omp atomic update
        mom_(ib,ia+1,ja+1,ka+1) = mom_(ib,ia+1,ja+1,ka+1) + aaa8*p1

!$omp atomic update
        ppara_(ib,ia,  ja,  ka  ) = ppara_(ib,ia,  ja,  ka  ) + aaa1*p2
!$omp atomic update
        ppara_(ib,ia+1,ja,  ka  ) = ppara_(ib,ia+1,ja,  ka  ) + aaa2*p2
!$omp atomic update
        ppara_(ib,ia,  ja+1,ka  ) = ppara_(ib,ia,  ja+1,ka  ) + aaa3*p2
!$omp atomic update
        ppara_(ib,ia+1,ja+1,ka  ) = ppara_(ib,ia+1,ja+1,ka  ) + aaa4*p2
!$omp atomic update
        ppara_(ib,ia,  ja,  ka+1) = ppara_(ib,ia,  ja,  ka+1) + aaa5*p2
!$omp atomic update
        ppara_(ib,ia+1,ja,  ka+1) = ppara_(ib,ia+1,ja,  ka+1) + aaa6*p2
!$omp atomic update
        ppara_(ib,ia,  ja+1,ka+1) = ppara_(ib,ia,  ja+1,ka+1) + aaa7*p2
!$omp atomic update
        ppara_(ib,ia+1,ja+1,ka+1) = ppara_(ib,ia+1,ja+1,ka+1) + aaa8*p2

!$omp atomic update
        pperp_(ib,ia,  ja,  ka  ) = pperp_(ib,ia,  ja,  ka  ) + aaa1*mu1
!$omp atomic update
        pperp_(ib,ia+1,ja,  ka  ) = pperp_(ib,ia+1,ja,  ka  ) + aaa2*mu1
!$omp atomic update
        pperp_(ib,ia,  ja+1,ka  ) = pperp_(ib,ia,  ja+1,ka  ) + aaa3*mu1
!$omp atomic update
        pperp_(ib,ia+1,ja+1,ka  ) = pperp_(ib,ia+1,ja+1,ka  ) + aaa4*mu1
!$omp atomic update
        pperp_(ib,ia,  ja,  ka+1) = pperp_(ib,ia,  ja,  ka+1) + aaa5*mu1
!$omp atomic update
        pperp_(ib,ia+1,ja,  ka+1) = pperp_(ib,ia+1,ja,  ka+1) + aaa6*mu1
!$omp atomic update
        pperp_(ib,ia,  ja+1,ka+1) = pperp_(ib,ia,  ja+1,ka+1) + aaa7*mu1
!$omp atomic update
        pperp_(ib,ia+1,ja+1,ka+1) = pperp_(ib,ia+1,ja+1,ka+1) + aaa8*mu1

!$omp atomic update
        qpara_(ib,ia,  ja,  ka  ) = qpara_(ib,ia,  ja,  ka  ) + aaa1*p3
!$omp atomic update
        qpara_(ib,ia+1,ja,  ka  ) = qpara_(ib,ia+1,ja,  ka  ) + aaa2*p3
!$omp atomic update
        qpara_(ib,ia,  ja+1,ka  ) = qpara_(ib,ia,  ja+1,ka  ) + aaa3*p3
!$omp atomic update
        qpara_(ib,ia+1,ja+1,ka  ) = qpara_(ib,ia+1,ja+1,ka  ) + aaa4*p3
!$omp atomic update
        qpara_(ib,ia,  ja,  ka+1) = qpara_(ib,ia,  ja,  ka+1) + aaa5*p3
!$omp atomic update
        qpara_(ib,ia+1,ja,  ka+1) = qpara_(ib,ia+1,ja,  ka+1) + aaa6*p3
!$omp atomic update
        qpara_(ib,ia,  ja+1,ka+1) = qpara_(ib,ia,  ja+1,ka+1) + aaa7*p3
!$omp atomic update
        qpara_(ib,ia+1,ja+1,ka+1) = qpara_(ib,ia+1,ja+1,ka+1) + aaa8*p3

!$omp atomic update
        qperp_(ib,ia,  ja,  ka  ) = qperp_(ib,ia,  ja,  ka  ) + aaa1*p1mu
!$omp atomic update
        qperp_(ib,ia+1,ja,  ka  ) = qperp_(ib,ia+1,ja,  ka  ) + aaa2*p1mu
!$omp atomic update
        qperp_(ib,ia,  ja+1,ka  ) = qperp_(ib,ia,  ja+1,ka  ) + aaa3*p1mu
!$omp atomic update
        qperp_(ib,ia+1,ja+1,ka  ) = qperp_(ib,ia+1,ja+1,ka  ) + aaa4*p1mu
!$omp atomic update
        qperp_(ib,ia,  ja,  ka+1) = qperp_(ib,ia,  ja,  ka+1) + aaa5*p1mu
!$omp atomic update
        qperp_(ib,ia+1,ja,  ka+1) = qperp_(ib,ia+1,ja,  ka+1) + aaa6*p1mu
!$omp atomic update
        qperp_(ib,ia,  ja+1,ka+1) = qperp_(ib,ia,  ja+1,ka+1) + aaa7*p1mu
!$omp atomic update
        qperp_(ib,ia+1,ja+1,ka+1) = qperp_(ib,ia+1,ja+1,ka+1) + aaa8*p1mu
       end do
!      call stop_timer('  moments:loop03') ! NEC TIMER

!      call start_timer('  moments:loop03after') ! NEC TIMER

!$omp target teams distribute parallel do private(i,j)
        do i = 1, lrzphi
        do j = 1, nwork
          dns(i,1,1) = dns(i,1,1) + dns_(j,i,1,1)
          mom(i,1,1) = mom(i,1,1) + mom_(j,i,1,1)
          ppara(i,1,1) = ppara(i,1,1) + ppara_(j,i,1,1)
          pperp(i,1,1) = pperp(i,1,1) + pperp_(j,i,1,1)
          qpara(i,1,1) = qpara(i,1,1) + qpara_(j,i,1,1)
          qperp(i,1,1) = qperp(i,1,1) + qperp_(j,i,1,1)
        end do
        end do

!      call stop_timer('  moments:loop03after') ! NEC TIMER
! smoothing

       cwwa = 0.5d0
!       cwwc =-1.d0/6.d0 

!2024-12-21s, correction suggested by Panith Adulsiriswad
!      call start_timer('  moments:com01') ! NEC TIMER
       call periodic_particle_mlt6b(dns,mom,ppara,pperp,qpara,qperp)
       call partsm1(dns,cwwa)
       call partsm1(mom,cwwa)

       call partsm1(ppara,cwwa)
       call partsm1(pperp,cwwa)

       call partsm1(qpara,cwwa)
       call partsm1(qperp,cwwa)
!      call stop_timer('  moments:com01') ! NEC TIMER
!2024-12-21e

! calculate density (per volume)

!      do k = 1, lphi
!      do j = 1, lz
!      do i = 1, lr
!      call start_timer('  moments:loop05') ! NEC TIMER
!$omp target teams distribute parallel do private(i,vol) !corrected by R. Seki
      do i = 1, lrzphi
        vol = 1.0d0/(grr(i,1,1)*dr*dz*dphi)
        dns(i,1,1) = dns(i,1,1)*vol
        mom(i,1,1) = mom(i,1,1)*vol
        ppara(i,1,1) = ppara(i,1,1)*vol*mass1
        pperp(i,1,1) = pperp(i,1,1)*vol*babs(i,1,1)
        qpara(i,1,1) = qpara(i,1,1)*vol*mass1**2
        qperp(i,1,1) = qperp(i,1,1)*vol*babs(i,1,1)*mass1
      end do
!      call stop_timer('  moments:loop05') ! NEC TIMER
!      end do
!      end do



! wall effect

      lr1 = lr - 1
      lz1 = lz - 1

      if(my_rank_r.eq.(mpi_proc_r-1))then
!      do k = 1, lphi
!      do j = 1, lz
!      call start_timer('  moments:loop06') ! NEC TIMER
!$omp target teams distribute parallel do private(j)
       do j = 1, lzphi
          dns(lr1,j,1) = dns(lr1,j,1) + dns(lr,j,1)
          mom(lr1,j,1) = mom(lr1,j,1) + mom(lr,j,1)
          dns(lr,j,1) = 0.0d0
          mom(lr,j,1) = 0.0d0

          ppara(lr1,j,1) = ppara(lr1,j,1) + ppara(lr,j,1)
          pperp(lr1,j,1) = pperp(lr1,j,1) + pperp(lr,j,1)
          ppara(lr,j,1) = 0.0d0
          pperp(lr,j,1) = 0.0d0

          qpara(lr1,j,1) = qpara(lr1,j,1) + qpara(lr,j,1)
          qperp(lr1,j,1) = qperp(lr1,j,1) + qperp(lr,j,1)
          qpara(lr,j,1) = 0.0d0
          qperp(lr,j,1) = 0.0d0
       end do
!      call stop_timer('  moments:loop06') ! NEC TIMER
!      end do
      end if

      if(my_rank_r.eq.0)then
!      do k = 1, lphi
!      do j = 1, lz
!      call start_timer('  moments:loop07') ! NEC TIMER
!$omp target teams distribute parallel do private(j)
       do j = 1, lzphi
          dns(2,j,1) = dns(2,j,1) + dns(1,j,1)
          mom(2,j,1) = mom(2,j,1) + mom(1,j,1)
          dns(1,j,1) = 0.0d0
          mom(1,j,1) = 0.0d0

          ppara(2,j,1) = ppara(2,j,1) + ppara(1,j,1)
          pperp(2,j,1) = pperp(2,j,1) + pperp(1,j,1)
          ppara(1,j,1) = 0.0d0
          pperp(1,j,1) = 0.0d0

          qpara(2,j,1) = qpara(2,j,1) + qpara(1,j,1)
          qperp(2,j,1) = qperp(2,j,1) + qperp(1,j,1)
          qpara(1,j,1) = 0.0d0
          qperp(1,j,1) = 0.0d0
       end do
!      call stop_timer('  moments:loop07') ! NEC TIMER
!      end do
      end if

      if(my_rank_z.eq.(mpi_proc_z-1))then
!      call start_timer('  moments:loop08') ! NEC TIMER
!$omp target teams distribute parallel do private(k,i) collapse(2)
       do k = 1, lphi
       do i = 1, lr
          dns(i,lz1,k) = dns(i,lz1,k) + dns(i,lz,k)
          mom(i,lz1,k) = mom(i,lz1,k) + mom(i,lz,k)
          dns(i,lz,k) = 0.0d0
          mom(i,lz,k) = 0.0d0

          ppara(i,lz1,k) = ppara(i,lz1,k) + ppara(i,lz,k)
          pperp(i,lz1,k) = pperp(i,lz1,k) + pperp(i,lz,k)
          ppara(i,lz,k) = 0.0d0
          pperp(i,lz,k) = 0.0d0

          qpara(i,lz1,k) = qpara(i,lz1,k) + qpara(i,lz,k)
          qperp(i,lz1,k) = qperp(i,lz1,k) + qperp(i,lz,k)
          qpara(i,lz,k) = 0.0d0
          qperp(i,lz,k) = 0.0d0
       end do
       end do
!      call stop_timer('  moments:loop08') ! NEC TIMER
      end if 

      if(my_rank_z.eq.0)then
!      call start_timer('  moments:loop09') ! NEC TIMER
!$omp target teams distribute parallel do private(k,i) collapse(2)
       do k = 1, lphi
       do i = 1, lr
          dns(i,2,k) = dns(i,2,k) + dns(i,1,k)
          mom(i,2,k) = mom(i,2,k) + mom(i,1,k)
          dns(i,1,k) = 0.0d0
          mom(i,1,k) = 0.0d0

          ppara(i,2,k) = ppara(i,2,k) + ppara(i,1,k)
          pperp(i,2,k) = pperp(i,2,k) + pperp(i,1,k)
          ppara(i,1,k) = 0.0d0
          pperp(i,1,k) = 0.0d0

          qpara(i,2,k) = qpara(i,2,k) + qpara(i,1,k)
          qperp(i,2,k) = qperp(i,2,k) + qperp(i,1,k)
          qpara(i,1,k) = 0.0d0
          qperp(i,1,k) = 0.0d0
       end do
       end do
!      call stop_timer('  moments:loop09') ! NEC TIMER
      end if 

      if(.not.flag_stored)then !2012-05-02
! take modes n_min <= n <=n_max
!      n_min = 1
!      n_max = 2
!      call start_timer('  moments:com02') ! NEC TIMER
!      call lowpass2_mlt2(n_min,n_max,dns,mom)
!      call lowpass2_mlt2(n_min,n_max,ppara,pperp)
!      call lowpass2_mlt2(n_min,n_max,qpara,qperp)
!      call stop_timer('  moments:com02') ! NEC TIMER
!2013-05-22 end

      end if !2012-05-02

!2014-08-28s
!      call start_timer('  moments:loop10') ! NEC TIMER
!$omp target teams distribute parallel do private(i)
      do i = 1, lrzphi
        dns(i,1,1) = dns(i,1,1) + dns0(i,1,1)
        mom(i,1,1) = mom(i,1,1) + mom0(i,1,1)
        ppara(i,1,1) = ppara(i,1,1) + ppara0(i,1,1)
        pperp(i,1,1) = pperp(i,1,1) + pperp0(i,1,1)/babs0(i,1,1)*babs(i,1,1)
        qpara(i,1,1) = qpara(i,1,1) + qpara0(i,1,1)
        qperp(i,1,1) = qperp(i,1,1) + qperp0(i,1,1)/babs0(i,1,1)*babs(i,1,1)
      end do
!      call stop_timer('  moments:loop10') ! NEC TIMER
!2014-08-28e

!            call wall_clock(d5)
!            t5 = d5 - d4
!            if(my_rank.eq.0)then
!            write(6,*)'density 5, t5=',t5
!            end if

end subroutine moments
!--------------------------------------------------------------------
subroutine emf_gyro(marker_num,marker_num_gyro,gc,gyro_phys &
                   ,flp,flp_gyro)
! modified for GK simulation on NEC SX-Aurora TSUBASA 2022-01-11
! gc_x(:,n) = 1:r, 2:z, 3:phi, 4:p, 5:mu, 6:weight, 7:active, 8:p_phi, 9:fff, 10:fnrml
! flp(nflp, marker_num); nflp=30, nflp_gyro=7
!--------------------------------------------------------------------
      use mpiset
      use grid
      use field, only:babs,babs0,er,ez,ephi,epara,br,bz,bphi,br0,bz0,bphi0,fld
      use gyro, only:ngyro,marker_each_gyro
      implicit none

      integer::marker_num,nvec,marker_num_gyro
      real(8)::gc(ngc2,marker_num)
      real(8)::gyro_phys(2,marker_num_gyro)
      real(8)::flp(nflp,marker_num) !2015-07-08
      real(8)::flp_gyro(nflp_gyro,marker_num_gyro)
      integer::i,j,k,l,m,n,ia,ja,ka,lr1,lz1
      real(8)::ar,ar1,az,az1,aphi,aphi1
!      real(8)::aaa1,aaa2,aaa3,aaa4,aaa5,aaa6,aaa7,aaa8
      real(8)::dr1,dz1,dphi1,ma_mi_r
      integer::kr,kz,kphi
      integer::nc,in
      real(8)::rngr1,rngr1_ac

!NEC start
      integer::nn,nsta,nend
      integer,parameter::nblkp=1024*1024
!2021-01-11s
      integer::ijk_a(marker_each_gyro,3)
      real(8)::aaa(marker_each_gyro,8)
      integer::nvec0,nmod
!2021-01-11e


      dr1 = 1.0d0/dr
      dz1 = 1.0d0/dz
      dphi1 = 1.0d0/dphi
      rngr1 = 1.0d0/dble(ngyro)

      kr   = 1 - kr_offset(my_rank)
      kz   = 1 - kz_offset(my_rank)
      kphi = 1 + lphi_shd - kphi_offset(my_rank)
!      kphi = 3 - kphi_offset(my_rank)
      ma_mi_r = major_r-minor_r



! for subroutine extract_em

!      flp_gyro = 0.0d0

!$omp target teams distribute parallel do collapse(3) private(k,j,i)
      do k=1,lphi
      do j=1,lz
      do i=1,lr
      fld(1,i,j,k) = er(i,j,k)
      fld(2,i,j,k) = ez(i,j,k)
      fld(3,i,j,k) = ephi(i,j,k)
      fld(4,i,j,k) = epara(i,j,k)
      fld(5,i,j,k) = br(i,j,k) - br0(i,j,k)
      fld(6,i,j,k) = bz(i,j,k) - bz0(i,j,k)
      fld(7,i,j,k) = bphi(i,j,k) - bphi0(i,j,k)
      end do
      end do
      end do

      do nn = 1, marker_num_gyro, nblkp
        nsta = nn
        nend = min(nn+nblkp-1,marker_num_gyro)

!$omp target teams distribute parallel do private(n,nc,ar1,ar,az1,aphi1,aphi)
      do n = nsta, nend

        nc =(n-1)/ngyro + 1
        ijk_a(n,1)=max(1,min(lr  -1,int((gyro_phys(1,n)-ma_mi_r)*dr1  ) + kr  ))
        ijk_a(n,2)=max(1,min(lz  -1,int(gyro_phys(2,n)          *dz1  ) + kz  ))
        ijk_a(n,3)=max(1,min(lphi-1,int(gc(3,nc)          *dphi1) + kphi))

        ar1  = max(0.0d0, min(1.0d0, (gyro_phys(1,n) - ma_mi_r)*dr1 -  dble(ijk_a(n,1) - kr)  ) )
        ar   = 1.0d0 - ar1
        az1  = max(0.0d0, min(1.0d0, gyro_phys(2,n)*dz1 - dble(ijk_a(n,2) - kz) ) )
        az   = 1.0d0 - az1
        aphi1 = max(0.0d0, min(1.0d0, gc(3,nc)*dphi1 - dble(ijk_a(n,3) - kphi) ) )
        aphi = 1.0d0 - aphi1

        aaa(n,1) = ar *az *aphi
        aaa(n,2) = ar1*az *aphi
        aaa(n,3) = ar *az1*aphi
        aaa(n,4) = ar1*az1*aphi
        aaa(n,5) = ar *az *aphi1
        aaa(n,6) = ar1*az *aphi1
        aaa(n,7) = ar *az1*aphi1
        aaa(n,8) = ar1*az1*aphi1
      end do

!NEC$ outerloop_unroll(4)
!$omp target teams distribute parallel do collapse(2) private(in,n,ia,ja,ka)
      do n = nsta, nend

        do in = 1, nflp_gyro

          ia=ijk_a(n,1)
          ja=ijk_a(n,2)
          ka=ijk_a(n,3)

          flp_gyro(in,n) = fld(in, ia, ja,  ka  )*aaa(n,1) + fld(in, ia+1,ja,  ka  )*aaa(n,2) &
                         + fld(in, ia, ja+1,ka  )*aaa(n,3) + fld(in, ia+1,ja+1,ka  )*aaa(n,4) &
                         + fld(in, ia, ja,  ka+1)*aaa(n,5) + fld(in, ia+1,ja,  ka+1)*aaa(n,6) &
                         + fld(in, ia, ja+1,ka+1)*aaa(n,7) + fld(in, ia+1,ja+1,ka+1)*aaa(n,8)
        end do

      end do

!NEC$ outerloop_unroll(4)
!$omp target teams distribute parallel do collapse(2) private(in,nc,rngr1_ac)
       do nc =(nsta-1)/ngyro+1, nend/ngyro
       do in = 1, nflp_gyro
          rngr1_ac = rngr1*gc(7,nc)
          flp(in,nc) =(flp_gyro(in,ngyro*(nc-1)+1) &
                      +flp_gyro(in,ngyro*(nc-1)+2) &
                      +flp_gyro(in,ngyro*(nc-1)+3) &
                      +flp_gyro(in,ngyro*(nc-1)+4) &
                      )*rngr1_ac
       end do
       end do

      end do !nn

end subroutine emf_gyro
!--------------------------------------------------------------------
subroutine density_gyro(marker_num,marker_num_gyro,mass,gc,gyro_phys &
                       ,dns,mom,ppara,pperp &
                       ,dns0,mom0,ppara0,pperp0)
! 2022-01-11, for FLR case
! modified for NEC SX-Aurora TSUBASA 2020-07-10
! gc_x(:,n) = 1:r, 2:z, 3:phi, 4:p, 5:mu, 6:weight, 7:active, 8:p_phi, 9:fff, 10:fnrml
!--------------------------------------------------------------------
      use mpiset
      use grid
      use field, only:babs,babs0
      use gyro, only:ngyro,marker_each_gyro
!      use mod_timer ! NEC TIMER
      implicit none

      integer::marker_num,marker_num_gyro
      real(8)::gc(ngc2,marker_num)
      real(8)::gyro_phys(2,marker_num_gyro)
      real(8)::mass
      real(8)::dns(lr,lz,lphi),mom(lr,lz,lphi)
      real(8)::ppara(lr,lz,lphi),pperp(lr,lz,lphi)
      real(8)::dns0(lr,lz,lphi),mom0(lr,lz,lphi)
      real(8)::ppara0(lr,lz,lphi),pperp0(lr,lz,lphi)
      integer::i,j,k,l,n,ia,ia1,ja,ja1,ka,ka1,lr1,lz1
      real(8)::ar,ar1,az,az1,aphi,aphi1
      real(8)::aaa1,aaa2,aaa3,aaa4,aaa5,aaa6,aaa7,aaa8
      real(8)::vol,cwwa,cwwc
!      real(8)::d0,d1,d2,d3,d4,d5,t1,t2,t3,t4,t5
      real(8)::dr1,dz1,dphi1,ma_mi_r
      integer::kr,kz,kphi
      real(8)::p0,p1,p2,mu1,mass1 !2015-09-21
      integer::n_min,n_max !2013-05-22
      real(8)::rngr1
      integer::nc
      integer::ib
      integer,parameter::nwork=32
      !real(8),dimension(lr,lz,lphi,nwork)::dns_,mom_,ppara_,pperp_
      real(8),dimension(nwork,lr,lz,lphi)::dns_,mom_,ppara_,pperp_

      mass1 = 1.0d0/mass

      dr1 = 1.0d0/dr
      dz1 = 1.0d0/dz
      dphi1 = 1.0d0/dphi
      rngr1 = 1.0d0/dble(ngyro)

      kr   = 1 - kr_offset(my_rank)
      kz   = 1 - kz_offset(my_rank)
      kphi = 1 + lphi_shd - kphi_offset(my_rank)
!      kphi = 3 - kphi_offset(my_rank)
      ma_mi_r = major_r-minor_r

!      call start_timer('  density_gyro:loop02') ! NEC TIMER
!$omp target teams distribute parallel do private(i)
        do i = 1, lrzphi
          dns(i,1,1) = 0.0d0
          mom(i,1,1) = 0.0d0
          ppara(i,1,1) = 0.0d0
          pperp(i,1,1) = 0.0d0
        end do
!      call stop_timer('  density_gyro:loop02') ! NEC TIMER

!      call start_timer('  density_gyro:loop03before') ! NEC TIMER
!$omp target teams distribute parallel do private(i)
        do i = 1, lrzphi*nwork
          dns_(i,1,1,1) = 0.0d0
          mom_(i,1,1,1) = 0.0d0
          ppara_(i,1,1,1) = 0.0d0
          pperp_(i,1,1,1) = 0.0d0
        end do
!      call stop_timer('  density_gyro:loop03before') ! NEC TIMER

!      call start_timer('  density_gyro:loop03') ! NEC TIMER
!$omp  target teams distribute parallel do
!$omp& private(n,nc,ia,ja,ka,ar1,ar,az1,az,aphi1,aphi,aaa1,aaa2,aaa3,aaa4,aaa5,aaa6,aaa7,aaa8)
!$omp& private(p0,p1,p2,mu1,ib)
       do n = 1, marker_num_gyro
        ib = mod(n-1,nwork)+1
        nc =(n-1)/ngyro + 1

! gc_x(:,n) = 1:r, 2:z, 3:phi, 4:p, 5:mu, 6:weight, 7:active, 8:p_phi, 9:fff, 10:fnrml
        ia=max(1,min(lr  -1,int((gyro_phys(1,n)-ma_mi_r)*dr1  ) + kr  ))
        ja=max(1,min(lz  -1,int(gyro_phys(2,n)          *dz1  ) + kz  ))
        ka=max(1,min(lphi-1,int(gc(3,nc)        *dphi1) + kphi))

        ar1 = (gyro_phys(1,n)-grr(ia,ja,ka ) ) *dr1
        ar  = 1.0d0 - ar1
        az1 = (gyro_phys(2,n)-gzz(ia,ja,ka ) ) *dz1
        az  = 1.0d0 - az1
        aphi1 = (gc(3,nc)-gphi(ia,ja,ka ) ) *dphi1
        aphi  = 1.0d0 - aphi1

        aaa1 = ar *az *aphi
        aaa2 = ar1*az *aphi
        aaa3 = ar *az1*aphi
        aaa4 = ar1*az1*aphi
        aaa5 = ar *az *aphi1
        aaa6 = ar1*az *aphi1
        aaa7 = ar *az1*aphi1
        aaa8 = ar1*az1*aphi1

        p0 =  gc(6,nc)*gc(7,nc)*rngr1 !rngr1 is multiplied
        p1 =  gc(4,nc)*p0
        p2 =  gc(4,nc)*p1
        mu1 = gc(5,nc)*p0

!$omp atomic update
        dns_(ib,ia,  ja,  ka  ) = dns_(ib,ia,  ja,  ka  ) + aaa1*p0
!$omp atomic update
        dns_(ib,ia+1,ja,  ka  ) = dns_(ib,ia+1,ja,  ka  ) + aaa2*p0
!$omp atomic update
        dns_(ib,ia,  ja+1,ka  ) = dns_(ib,ia,  ja+1,ka  ) + aaa3*p0
!$omp atomic update
        dns_(ib,ia+1,ja+1,ka  ) = dns_(ib,ia+1,ja+1,ka  ) + aaa4*p0
!$omp atomic update
        dns_(ib,ia,  ja,  ka+1) = dns_(ib,ia,  ja,  ka+1) + aaa5*p0
!$omp atomic update
        dns_(ib,ia+1,ja,  ka+1) = dns_(ib,ia+1,ja,  ka+1) + aaa6*p0
!$omp atomic update
        dns_(ib,ia,  ja+1,ka+1) = dns_(ib,ia,  ja+1,ka+1) + aaa7*p0
!$omp atomic update
        dns_(ib,ia+1,ja+1,ka+1) = dns_(ib,ia+1,ja+1,ka+1) + aaa8*p0

!$omp atomic update
        mom_(ib,ia,  ja,  ka  ) = mom_(ib,ia,  ja,  ka  ) + aaa1*p1
!$omp atomic update
        mom_(ib,ia+1,ja,  ka  ) = mom_(ib,ia+1,ja,  ka  ) + aaa2*p1
!$omp atomic update
        mom_(ib,ia,  ja+1,ka  ) = mom_(ib,ia,  ja+1,ka  ) + aaa3*p1
!$omp atomic update
        mom_(ib,ia+1,ja+1,ka  ) = mom_(ib,ia+1,ja+1,ka  ) + aaa4*p1
!$omp atomic update
        mom_(ib,ia,  ja,  ka+1) = mom_(ib,ia,  ja,  ka+1) + aaa5*p1
!$omp atomic update
        mom_(ib,ia+1,ja,  ka+1) = mom_(ib,ia+1,ja,  ka+1) + aaa6*p1
!$omp atomic update
        mom_(ib,ia,  ja+1,ka+1) = mom_(ib,ia,  ja+1,ka+1) + aaa7*p1
!$omp atomic update
        mom_(ib,ia+1,ja+1,ka+1) = mom_(ib,ia+1,ja+1,ka+1) + aaa8*p1

!$omp atomic update
        ppara_(ib,ia,  ja,  ka  ) = ppara_(ib,ia,  ja,  ka  ) + aaa1*p2
!$omp atomic update
        ppara_(ib,ia+1,ja,  ka  ) = ppara_(ib,ia+1,ja,  ka  ) + aaa2*p2
!$omp atomic update
        ppara_(ib,ia,  ja+1,ka  ) = ppara_(ib,ia,  ja+1,ka  ) + aaa3*p2
!$omp atomic update
        ppara_(ib,ia+1,ja+1,ka  ) = ppara_(ib,ia+1,ja+1,ka  ) + aaa4*p2
!$omp atomic update
        ppara_(ib,ia,  ja,  ka+1) = ppara_(ib,ia,  ja,  ka+1) + aaa5*p2
!$omp atomic update
        ppara_(ib,ia+1,ja,  ka+1) = ppara_(ib,ia+1,ja,  ka+1) + aaa6*p2
!$omp atomic update
        ppara_(ib,ia,  ja+1,ka+1) = ppara_(ib,ia,  ja+1,ka+1) + aaa7*p2
!$omp atomic update
        ppara_(ib,ia+1,ja+1,ka+1) = ppara_(ib,ia+1,ja+1,ka+1) + aaa8*p2

!$omp atomic update
        pperp_(ib,ia,  ja,  ka  ) = pperp_(ib,ia,  ja,  ka  ) + aaa1*mu1
!$omp atomic update
        pperp_(ib,ia+1,ja,  ka  ) = pperp_(ib,ia+1,ja,  ka  ) + aaa2*mu1
!$omp atomic update
        pperp_(ib,ia,  ja+1,ka  ) = pperp_(ib,ia,  ja+1,ka  ) + aaa3*mu1
!$omp atomic update
        pperp_(ib,ia+1,ja+1,ka  ) = pperp_(ib,ia+1,ja+1,ka  ) + aaa4*mu1
!$omp atomic update
        pperp_(ib,ia,  ja,  ka+1) = pperp_(ib,ia,  ja,  ka+1) + aaa5*mu1
!$omp atomic update
        pperp_(ib,ia+1,ja,  ka+1) = pperp_(ib,ia+1,ja,  ka+1) + aaa6*mu1
!$omp atomic update
        pperp_(ib,ia,  ja+1,ka+1) = pperp_(ib,ia,  ja+1,ka+1) + aaa7*mu1
!$omp atomic update
        pperp_(ib,ia+1,ja+1,ka+1) = pperp_(ib,ia+1,ja+1,ka+1) + aaa8*mu1
       end do
!      call stop_timer('  density_gyro:loop03') ! NEC TIMER

!      call start_timer('  density_gyro:loop03after') ! NEC TIMER
!$omp target teams distribute parallel do private(i,j)
        do i = 1, lrzphi
        do j = 1, nwork
          dns(i,1,1) = dns(i,1,1) + dns_(j,i,1,1)
          mom(i,1,1) = mom(i,1,1) + mom_(j,i,1,1)
          ppara(i,1,1) = ppara(i,1,1) + ppara_(j,i,1,1)
          pperp(i,1,1) = pperp(i,1,1) + pperp_(j,i,1,1)
        end do
        end do
!      call stop_timer('  density_gyro:loop03after') ! NEC TIMER
! smoothing

       cwwa = 0.5d0
!       cwwc =-1.d0/6.d0 

!2024-12-21s, correction suggested by Panith Adulsiriswad
!      call start_timer('  density_gyro:com01') ! NEC TIMER
       call periodic_particle_mlt4b(dns,mom,ppara,pperp)
       call partsm1(dns,cwwa)
       call partsm1(mom,cwwa)

       call partsm1(ppara,cwwa)
       call partsm1(pperp,cwwa)
!      call stop_timer('  density_gyro:com01') ! NEC TIMER
!2024-12-21e

! calculate density (per volume)

!      do k = 1, lphi
!      do j = 1, lz
!      do i = 1, lr
!      call start_timer('  density_gyro:loop05') ! NEC TIMER
!$omp target teams distribute parallel do private(i,vol) !corrected by R. Seki
      do i = 1, lrzphi
        vol = 1.0d0/(grr(i,1,1)*dr*dz*dphi)
        dns(i,1,1) = dns(i,1,1)*vol
        mom(i,1,1) = mom(i,1,1)*vol
        ppara(i,1,1) = ppara(i,1,1)*vol/mass
        pperp(i,1,1) = pperp(i,1,1)*vol*babs(i,1,1)
      end do
!      call stop_timer('  density_gyro:loop05') ! NEC TIMER
!      end do
!      end do


! wall effect

      lr1 = lr - 1
      lz1 = lz - 1

      if(my_rank_r.eq.(mpi_proc_r-1))then
!      do k = 1, lphi
!      do j = 1, lz
!      call start_timer('  density_gyro:loop06') ! NEC TIMER
!$omp target teams distribute parallel do private(j)
       do j = 1, lzphi
          dns(lr1,j,1) = dns(lr1,j,1) + dns(lr,j,1)
          mom(lr1,j,1) = mom(lr1,j,1) + mom(lr,j,1)
          dns(lr,j,1) = 0.0d0
          mom(lr,j,1) = 0.0d0

          ppara(lr1,j,1) = ppara(lr1,j,1) + ppara(lr,j,1)
          pperp(lr1,j,1) = pperp(lr1,j,1) + pperp(lr,j,1)
          ppara(lr,j,1) = 0.0d0
          pperp(lr,j,1) = 0.0d0
       end do
!      call stop_timer('  density_gyro:loop06') ! NEC TIMER
!      end do
      end if

      if(my_rank_r.eq.0)then
!      do k = 1, lphi
!      do j = 1, lz
!      call start_timer('  density_gyro:loop07') ! NEC TIMER
!$omp target teams distribute parallel do private(j)
       do j = 1, lzphi
          dns(2,j,1) = dns(2,j,1) + dns(1,j,1)
          mom(2,j,1) = mom(2,j,1) + mom(1,j,1)
          dns(1,j,1) = 0.0d0
          mom(1,j,1) = 0.0d0

          ppara(2,j,1) = ppara(2,j,1) + ppara(1,j,1)
          pperp(2,j,1) = pperp(2,j,1) + pperp(1,j,1)
          ppara(1,j,1) = 0.0d0
          pperp(1,j,1) = 0.0d0
       end do
!      call stop_timer('  density_gyro:loop07') ! NEC TIMER
!      end do
      end if

      if(my_rank_z.eq.(mpi_proc_z-1))then
!      call start_timer('  density_gyro:loop08') ! NEC TIMER
!$omp target teams distribute parallel do private(k,i) collapse(2)
       do k = 1, lphi
       do i = 1, lr
          dns(i,lz1,k) = dns(i,lz1,k) + dns(i,lz,k)
          mom(i,lz1,k) = mom(i,lz1,k) + mom(i,lz,k)
          dns(i,lz,k) = 0.0d0
          mom(i,lz,k) = 0.0d0

          ppara(i,lz1,k) = ppara(i,lz1,k) + ppara(i,lz,k)
          pperp(i,lz1,k) = pperp(i,lz1,k) + pperp(i,lz,k)
          ppara(i,lz,k) = 0.0d0
          pperp(i,lz,k) = 0.0d0
       end do
       end do
!      call stop_timer('  density_gyro:loop08') ! NEC TIMER
      end if 

      if(my_rank_z.eq.0)then
!      call start_timer('  density_gyro:loop09') ! NEC TIMER
!$omp target teams distribute parallel do private(k,i) collapse(2)
       do k = 1, lphi
       do i = 1, lr
          dns(i,2,k) = dns(i,2,k) + dns(i,1,k)
          mom(i,2,k) = mom(i,2,k) + mom(i,1,k)
          dns(i,1,k) = 0.0d0
          mom(i,1,k) = 0.0d0

          ppara(i,2,k) = ppara(i,2,k) + ppara(i,1,k)
          pperp(i,2,k) = pperp(i,2,k) + pperp(i,1,k)
          ppara(i,1,k) = 0.0d0
          pperp(i,1,k) = 0.0d0
       end do
       end do
!      call stop_timer('  density_gyro:loop09') ! NEC TIMER
      end if 

! take only n=1 modes
!      n = 1
!      call n1(n,ppara)
!      call n1(n,pperp)

! take modes n_min <= n <=n_max
!      n_min = 1
!      n_max = 2 
!      call start_timer('  density_gyro:com02') ! NEC TIMER
!      call lowpass2_mlt2(n_min,n_max,dns,mom)
!      call lowpass2_mlt2(n_min,n_max,ppara,pperp)
!      call stop_timer('  density_gyro:com02') ! NEC TIMER

!      do k = 1, lphi
!      do j = 1, lz
!      do i = 1, lr
!      call start_timer('  density_gyro:loop10') ! NEC TIMER
!$omp target teams distribute parallel do private(i)
      do i = 1, lrzphi
        dns(i,1,1) = dns(i,1,1) + dns0(i,1,1)
        mom(i,1,1) = mom(i,1,1) + mom0(i,1,1)
        ppara(i,1,1) = ppara(i,1,1) + ppara0(i,1,1)
        pperp(i,1,1) = pperp(i,1,1) + pperp0(i,1,1)/babs0(i,1,1)*babs(i,1,1)
      end do
!      call stop_timer('  density_gyro:loop10') ! NEC TIMER
!      end do
!      end do

end subroutine density_gyro
!--------------------------------------------------------------------
subroutine moments_gyro(marker_num,marker_num_gyro,mass,gc,gyro_phys &
                  ,dns,mom,ppara,pperp,qpara,qperp &
                  ,dns0,mom0,ppara0,pperp0,qpara0,qperp0)
! 2016-02-04, for FLR case
! modified for NEC SX-Aurora TSUBASA 2020-07-10
! gc_x(:,n) = 1:r, 2:z, 3:phi, 4:p, 5:mu, 6:weight, 7:active, 8:p_phi, 9:fff, 10:fnrml
!--------------------------------------------------------------------
      use mpiset
      use grid
      use field, only:babs,babs0
      use gyro, only:ngyro,marker_each_gyro
!      use mod_timer ! NEC TIMER
      implicit none

      integer::marker_num,nvec,marker_num_gyro
      real(8)::gc(ngc2,marker_num)
      real(8)::gyro_phys(2,marker_num_gyro)
      real(8)::mass
      real(8)::dns(lr,lz,lphi),mom(lr,lz,lphi)
      real(8)::ppara(lr,lz,lphi),pperp(lr,lz,lphi)
      real(8)::qpara(lr,lz,lphi),qperp(lr,lz,lphi)
      real(8)::dns0(lr,lz,lphi),mom0(lr,lz,lphi)
      real(8)::ppara0(lr,lz,lphi),pperp0(lr,lz,lphi)
      real(8)::qpara0(lr,lz,lphi),qperp0(lr,lz,lphi)
      integer::i,j,k,l,m,n,ia,ia1,ja,ja1,ka,ka1,lr1,lz1
      real(8)::ar,ar1,az,az1,aphi,aphi1
      real(8)::aaa1,aaa2,aaa3,aaa4,aaa5,aaa6,aaa7,aaa8
      real(8)::vol,cwwa,cwwc
!      real(8)::d0,d1,d2,d3,d4,d5,t1,t2,t3,t4,t5
      real(8)::dr1,dz1,dphi1,ma_mi_r
      integer::kr,kz,kphi
      real(8)::p0,p1,p2,mu1,p3,p1mu,mass1 !2015-09-21
      integer::n_min,n_max !2013-05-22
      real(8)::rngr1
      integer::nc
      integer::ib
      integer,parameter::nwork=32
      real(8),dimension(nwork,lr,lz,lphi)::dns_,mom_,ppara_,pperp_,qpara_,qperp_

      mass1 = 1.0d0/mass

      dr1 = 1.0d0/dr
      dz1 = 1.0d0/dz
      dphi1 = 1.0d0/dphi
      rngr1 = 1.0d0/dble(ngyro)

      kr   = 1 - kr_offset(my_rank)
      kz   = 1 - kz_offset(my_rank)
      kphi = 1 + lphi_shd - kphi_offset(my_rank)
!      kphi = 3 - kphi_offset(my_rank)
      ma_mi_r = major_r-minor_r

!      call start_timer('  moments_gyro:loop02') ! NEC TIMER
!$omp target teams distribute parallel do private(i)
        do i = 1, lrzphi
          dns(i,1,1) = 0.0d0
          mom(i,1,1) = 0.0d0
          ppara(i,1,1) = 0.0d0
          pperp(i,1,1) = 0.0d0
          qpara(i,1,1) = 0.0d0
          qperp(i,1,1) = 0.0d0
        end do
!      call stop_timer('  moments_gyro:loop02') ! NEC TIMER

!      call start_timer('  moments_gyro:loop03before') ! NEC TIMER
!$omp target teams distribute parallel do private(i)
        do i = 1, lrzphi*nwork
          dns_(i,1,1,1) = 0.0d0
          mom_(i,1,1,1) = 0.0d0
          ppara_(i,1,1,1) = 0.0d0
          pperp_(i,1,1,1) = 0.0d0
          qpara_(i,1,1,1) = 0.0d0
          qperp_(i,1,1,1) = 0.0d0
        end do
!      call stop_timer('  moments_gyro:loop03before') ! NEC TIMER

!      call start_timer('  moments_gyro:loop03') ! NEC TIMER
!$omp  target teams distribute parallel do
!$omp& private(n,nc,ia,ja,ka,ar1,ar,az1,az,aphi1,aphi,aaa1,aaa2,aaa3,aaa4,aaa5,aaa6,aaa7,aaa8)
!$omp& private(p0,p1,p2,mu1,p3,p1mu,ib)
       do n = 1, marker_num_gyro
        ib = mod(n-1,nwork)+1
        nc =(n-1)/ngyro + 1

! gc_x(:,n) = 1:r, 2:z, 3:phi, 4:p, 5:mu, 6:weight, 7:active, 8:p_phi, 9:fff, 10:fnrml
        ia=max(1,min(lr  -1,int((gyro_phys(1,n)-ma_mi_r)*dr1  ) + kr  ))
        ja=max(1,min(lz  -1,int(gyro_phys(2,n)          *dz1  ) + kz  ))
        ka=max(1,min(lphi-1,int(gc(3,nc)        *dphi1) + kphi))

        ar1 = (gyro_phys(1,n)-grr(ia,ja,ka ) ) *dr1
        ar  = 1.0d0 - ar1
        az1 = (gyro_phys(2,n)-gzz(ia,ja,ka ) ) *dz1
        az  = 1.0d0 - az1
        aphi1 = (gc(3,nc)-gphi(ia,ja,ka ) ) *dphi1
        aphi  = 1.0d0 - aphi1

        aaa1 = ar *az *aphi
        aaa2 = ar1*az *aphi
        aaa3 = ar *az1*aphi
        aaa4 = ar1*az1*aphi
        aaa5 = ar *az *aphi1
        aaa6 = ar1*az *aphi1
        aaa7 = ar *az1*aphi1
        aaa8 = ar1*az1*aphi1

        p0 =  gc(6,nc)*gc(7,nc)*rngr1 !rngr1 is multiplied
        p1 =  gc(4,nc)*p0
        p2 =  gc(4,nc)*p1
        mu1 = gc(5,nc)*p0
        p3  = gc(4,nc)*p2
        p1mu= gc(4,nc)*mu1

!$omp atomic update
        dns_(ib,ia,  ja,  ka  ) = dns_(ib,ia,  ja,  ka  ) + aaa1*p0
!$omp atomic update
        dns_(ib,ia+1,ja,  ka  ) = dns_(ib,ia+1,ja,  ka  ) + aaa2*p0
!$omp atomic update
        dns_(ib,ia,  ja+1,ka  ) = dns_(ib,ia,  ja+1,ka  ) + aaa3*p0
!$omp atomic update
        dns_(ib,ia+1,ja+1,ka  ) = dns_(ib,ia+1,ja+1,ka  ) + aaa4*p0
!$omp atomic update
        dns_(ib,ia,  ja,  ka+1) = dns_(ib,ia,  ja,  ka+1) + aaa5*p0
!$omp atomic update
        dns_(ib,ia+1,ja,  ka+1) = dns_(ib,ia+1,ja,  ka+1) + aaa6*p0
!$omp atomic update
        dns_(ib,ia,  ja+1,ka+1) = dns_(ib,ia,  ja+1,ka+1) + aaa7*p0
!$omp atomic update
        dns_(ib,ia+1,ja+1,ka+1) = dns_(ib,ia+1,ja+1,ka+1) + aaa8*p0

!$omp atomic update
        mom_(ib,ia,  ja,  ka  ) = mom_(ib,ia,  ja,  ka  ) + aaa1*p1
!$omp atomic update
        mom_(ib,ia+1,ja,  ka  ) = mom_(ib,ia+1,ja,  ka  ) + aaa2*p1
!$omp atomic update
        mom_(ib,ia,  ja+1,ka  ) = mom_(ib,ia,  ja+1,ka  ) + aaa3*p1
!$omp atomic update
        mom_(ib,ia+1,ja+1,ka  ) = mom_(ib,ia+1,ja+1,ka  ) + aaa4*p1
!$omp atomic update
        mom_(ib,ia,  ja,  ka+1) = mom_(ib,ia,  ja,  ka+1) + aaa5*p1
!$omp atomic update
        mom_(ib,ia+1,ja,  ka+1) = mom_(ib,ia+1,ja,  ka+1) + aaa6*p1
!$omp atomic update
        mom_(ib,ia,  ja+1,ka+1) = mom_(ib,ia,  ja+1,ka+1) + aaa7*p1
!$omp atomic update
        mom_(ib,ia+1,ja+1,ka+1) = mom_(ib,ia+1,ja+1,ka+1) + aaa8*p1

!$omp atomic update
        ppara_(ib,ia,  ja,  ka  ) = ppara_(ib,ia,  ja,  ka  ) + aaa1*p2
!$omp atomic update
        ppara_(ib,ia+1,ja,  ka  ) = ppara_(ib,ia+1,ja,  ka  ) + aaa2*p2
!$omp atomic update
        ppara_(ib,ia,  ja+1,ka  ) = ppara_(ib,ia,  ja+1,ka  ) + aaa3*p2
!$omp atomic update
        ppara_(ib,ia+1,ja+1,ka  ) = ppara_(ib,ia+1,ja+1,ka  ) + aaa4*p2
!$omp atomic update
        ppara_(ib,ia,  ja,  ka+1) = ppara_(ib,ia,  ja,  ka+1) + aaa5*p2
!$omp atomic update
        ppara_(ib,ia+1,ja,  ka+1) = ppara_(ib,ia+1,ja,  ka+1) + aaa6*p2
!$omp atomic update
        ppara_(ib,ia,  ja+1,ka+1) = ppara_(ib,ia,  ja+1,ka+1) + aaa7*p2
!$omp atomic update
        ppara_(ib,ia+1,ja+1,ka+1) = ppara_(ib,ia+1,ja+1,ka+1) + aaa8*p2

!$omp atomic update
        pperp_(ib,ia,  ja,  ka  ) = pperp_(ib,ia,  ja,  ka  ) + aaa1*mu1
!$omp atomic update
        pperp_(ib,ia+1,ja,  ka  ) = pperp_(ib,ia+1,ja,  ka  ) + aaa2*mu1
!$omp atomic update
        pperp_(ib,ia,  ja+1,ka  ) = pperp_(ib,ia,  ja+1,ka  ) + aaa3*mu1
!$omp atomic update
        pperp_(ib,ia+1,ja+1,ka  ) = pperp_(ib,ia+1,ja+1,ka  ) + aaa4*mu1
!$omp atomic update
        pperp_(ib,ia,  ja,  ka+1) = pperp_(ib,ia,  ja,  ka+1) + aaa5*mu1
!$omp atomic update
        pperp_(ib,ia+1,ja,  ka+1) = pperp_(ib,ia+1,ja,  ka+1) + aaa6*mu1
!$omp atomic update
        pperp_(ib,ia,  ja+1,ka+1) = pperp_(ib,ia,  ja+1,ka+1) + aaa7*mu1
!$omp atomic update
        pperp_(ib,ia+1,ja+1,ka+1) = pperp_(ib,ia+1,ja+1,ka+1) + aaa8*mu1

!$omp atomic update
        qpara_(ib,ia,  ja,  ka  ) = qpara_(ib,ia,  ja,  ka  ) + aaa1*p3
!$omp atomic update
        qpara_(ib,ia+1,ja,  ka  ) = qpara_(ib,ia+1,ja,  ka  ) + aaa2*p3
!$omp atomic update
        qpara_(ib,ia,  ja+1,ka  ) = qpara_(ib,ia,  ja+1,ka  ) + aaa3*p3
!$omp atomic update
        qpara_(ib,ia+1,ja+1,ka  ) = qpara_(ib,ia+1,ja+1,ka  ) + aaa4*p3
!$omp atomic update
        qpara_(ib,ia,  ja,  ka+1) = qpara_(ib,ia,  ja,  ka+1) + aaa5*p3
!$omp atomic update
        qpara_(ib,ia+1,ja,  ka+1) = qpara_(ib,ia+1,ja,  ka+1) + aaa6*p3
!$omp atomic update
        qpara_(ib,ia,  ja+1,ka+1) = qpara_(ib,ia,  ja+1,ka+1) + aaa7*p3
!$omp atomic update
        qpara_(ib,ia+1,ja+1,ka+1) = qpara_(ib,ia+1,ja+1,ka+1) + aaa8*p3

!$omp atomic update
        qperp_(ib,ia,  ja,  ka  ) = qperp_(ib,ia,  ja,  ka  ) + aaa1*p1mu
!$omp atomic update
        qperp_(ib,ia+1,ja,  ka  ) = qperp_(ib,ia+1,ja,  ka  ) + aaa2*p1mu
!$omp atomic update
        qperp_(ib,ia,  ja+1,ka  ) = qperp_(ib,ia,  ja+1,ka  ) + aaa3*p1mu
!$omp atomic update
        qperp_(ib,ia+1,ja+1,ka  ) = qperp_(ib,ia+1,ja+1,ka  ) + aaa4*p1mu
!$omp atomic update
        qperp_(ib,ia,  ja,  ka+1) = qperp_(ib,ia,  ja,  ka+1) + aaa5*p1mu
!$omp atomic update
        qperp_(ib,ia+1,ja,  ka+1) = qperp_(ib,ia+1,ja,  ka+1) + aaa6*p1mu
!$omp atomic update
        qperp_(ib,ia,  ja+1,ka+1) = qperp_(ib,ia,  ja+1,ka+1) + aaa7*p1mu
!$omp atomic update
        qperp_(ib,ia+1,ja+1,ka+1) = qperp_(ib,ia+1,ja+1,ka+1) + aaa8*p1mu
       end do
!      call stop_timer('  moments_gyro:loop03') ! NEC TIMER

!      call start_timer('  moments_gyro:loop03after') ! NEC TIMER

!$omp target teams distribute parallel do private(i,j)
        do i = 1, lrzphi
        do j = 1, nwork
          dns(i,1,1) = dns(i,1,1) + dns_(j,i,1,1)
          mom(i,1,1) = mom(i,1,1) + mom_(j,i,1,1)
          ppara(i,1,1) = ppara(i,1,1) + ppara_(j,i,1,1)
          pperp(i,1,1) = pperp(i,1,1) + pperp_(j,i,1,1)
          qpara(i,1,1) = qpara(i,1,1) + qpara_(j,i,1,1)
          qperp(i,1,1) = qperp(i,1,1) + qperp_(j,i,1,1)
        end do
        end do

!      call stop_timer('  moments_gyro:loop03after') ! NEC TIMER
! smoothing

       cwwa = 0.5d0
!       cwwc =-1.d0/6.d0 

!2024-12-21s, correction suggested by Panith Adulsiriswad
!      call start_timer('  moments_gyro:com01') ! NEC TIMER
       call periodic_particle_mlt6b(dns,mom,ppara,pperp,qpara,qperp)
       call partsm1(dns,cwwa)
       call partsm1(mom,cwwa)

       call partsm1(ppara,cwwa)
       call partsm1(pperp,cwwa)

       call partsm1(qpara,cwwa)
       call partsm1(qperp,cwwa)
!      call stop_timer('  moments_gyro:com01') ! NEC TIMER
!2024-12-21e

! calculate density (per volume)

!      do k = 1, lphi
!      do j = 1, lz
!      do i = 1, lr
!      call start_timer('  moments_gyro:loop05') ! NEC TIMER
!$omp target teams distribute parallel do private(i,vol) !corrected by R. Seki
      do i = 1, lrzphi
        vol = 1.0d0/(grr(i,1,1)*dr*dz*dphi)
        dns(i,1,1) = dns(i,1,1)*vol
        mom(i,1,1) = mom(i,1,1)*vol
        ppara(i,1,1) = ppara(i,1,1)*vol*mass1
        pperp(i,1,1) = pperp(i,1,1)*vol*babs(i,1,1)
        qpara(i,1,1) = qpara(i,1,1)*vol*mass1**2
        qperp(i,1,1) = qperp(i,1,1)*vol*babs(i,1,1)*mass1
      end do
!      call stop_timer('  moments_gyro:loop05') ! NEC TIMER
!      end do
!      end do



! wall effect

      lr1 = lr - 1
      lz1 = lz - 1

      if(my_rank_r.eq.(mpi_proc_r-1))then
!      do k = 1, lphi
!      do j = 1, lz
!      call start_timer('  moments_gyro:loop06') ! NEC TIMER
!$omp target teams distribute parallel do private(j)
       do j = 1, lzphi
          dns(lr1,j,1) = dns(lr1,j,1) + dns(lr,j,1)
          mom(lr1,j,1) = mom(lr1,j,1) + mom(lr,j,1)
          dns(lr,j,1) = 0.0d0
          mom(lr,j,1) = 0.0d0

          ppara(lr1,j,1) = ppara(lr1,j,1) + ppara(lr,j,1)
          pperp(lr1,j,1) = pperp(lr1,j,1) + pperp(lr,j,1)
          ppara(lr,j,1) = 0.0d0
          pperp(lr,j,1) = 0.0d0

          qpara(lr1,j,1) = qpara(lr1,j,1) + qpara(lr,j,1)
          qperp(lr1,j,1) = qperp(lr1,j,1) + qperp(lr,j,1)
          qpara(lr,j,1) = 0.0d0
          qperp(lr,j,1) = 0.0d0
       end do
!      call stop_timer('  moments_gyro:loop06') ! NEC TIMER
!      end do
      end if

      if(my_rank_r.eq.0)then
!      do k = 1, lphi
!      do j = 1, lz
!      call start_timer('  moments_gyro:loop07') ! NEC TIMER
!$omp target teams distribute parallel do private(j)
       do j = 1, lzphi
          dns(2,j,1) = dns(2,j,1) + dns(1,j,1)
          mom(2,j,1) = mom(2,j,1) + mom(1,j,1)
          dns(1,j,1) = 0.0d0
          mom(1,j,1) = 0.0d0

          ppara(2,j,1) = ppara(2,j,1) + ppara(1,j,1)
          pperp(2,j,1) = pperp(2,j,1) + pperp(1,j,1)
          ppara(1,j,1) = 0.0d0
          pperp(1,j,1) = 0.0d0

          qpara(2,j,1) = qpara(2,j,1) + qpara(1,j,1)
          qperp(2,j,1) = qperp(2,j,1) + qperp(1,j,1)
          qpara(1,j,1) = 0.0d0
          qperp(1,j,1) = 0.0d0
       end do
!      call stop_timer('  moments_gyro:loop07') ! NEC TIMER
!      end do
      end if

      if(my_rank_z.eq.(mpi_proc_z-1))then
!      call start_timer('  moments_gyro:loop08') ! NEC TIMER
!$omp target teams distribute parallel do private(k,i) collapse(2)
       do k = 1, lphi
       do i = 1, lr
          dns(i,lz1,k) = dns(i,lz1,k) + dns(i,lz,k)
          mom(i,lz1,k) = mom(i,lz1,k) + mom(i,lz,k)
          dns(i,lz,k) = 0.0d0
          mom(i,lz,k) = 0.0d0

          ppara(i,lz1,k) = ppara(i,lz1,k) + ppara(i,lz,k)
          pperp(i,lz1,k) = pperp(i,lz1,k) + pperp(i,lz,k)
          ppara(i,lz,k) = 0.0d0
          pperp(i,lz,k) = 0.0d0

          qpara(i,lz1,k) = qpara(i,lz1,k) + qpara(i,lz,k)
          qperp(i,lz1,k) = qperp(i,lz1,k) + qperp(i,lz,k)
          qpara(i,lz,k) = 0.0d0
          qperp(i,lz,k) = 0.0d0
       end do
       end do
!      call stop_timer('  moments_gyro:loop08') ! NEC TIMER
      end if 

      if(my_rank_z.eq.0)then
!      call start_timer('  moments_gyro:loop09') ! NEC TIMER
!$omp target teams distribute parallel do private(k,i) collapse(2)
       do k = 1, lphi
       do i = 1, lr
          dns(i,2,k) = dns(i,2,k) + dns(i,1,k)
          mom(i,2,k) = mom(i,2,k) + mom(i,1,k)
          dns(i,1,k) = 0.0d0
          mom(i,1,k) = 0.0d0

          ppara(i,2,k) = ppara(i,2,k) + ppara(i,1,k)
          pperp(i,2,k) = pperp(i,2,k) + pperp(i,1,k)
          ppara(i,1,k) = 0.0d0
          pperp(i,1,k) = 0.0d0

          qpara(i,2,k) = qpara(i,2,k) + qpara(i,1,k)
          qperp(i,2,k) = qperp(i,2,k) + qperp(i,1,k)
          qpara(i,1,k) = 0.0d0
          qperp(i,1,k) = 0.0d0
       end do
       end do
!      call stop_timer('  moments_gyro:loop09') ! NEC TIMER
      end if 

      if(.not.flag_stored)then !2012-05-02
! take modes n_min <= n <=n_max
!      n_min = 1
!      n_max = 2
!      call start_timer('  moments_gyro:com02') ! NEC TIMER
!      call lowpass2_mlt2(n_min,n_max,dns,mom)
!      call lowpass2_mlt2(n_min,n_max,ppara,pperp)
!      call lowpass2_mlt2(n_min,n_max,qpara,qperp)
!      call stop_timer('  moments_gyro:com02') ! NEC TIMER
!2013-05-22 end

      end if !2012-05-02

!2014-08-28s
!      call start_timer('  moments_gyro:loop10') ! NEC TIMER
!$omp target teams distribute parallel do private(i)
      do i = 1, lrzphi
        dns(i,1,1) = dns(i,1,1) + dns0(i,1,1)
        mom(i,1,1) = mom(i,1,1) + mom0(i,1,1)
        ppara(i,1,1) = ppara(i,1,1) + ppara0(i,1,1)
        pperp(i,1,1) = pperp(i,1,1) + pperp0(i,1,1)/babs0(i,1,1)*babs(i,1,1)
        qpara(i,1,1) = qpara(i,1,1) + qpara0(i,1,1)
        qperp(i,1,1) = qperp(i,1,1) + qperp0(i,1,1)/babs0(i,1,1)*babs(i,1,1)
      end do
!      call stop_timer('  moments_gyro:loop10') ! NEC TIMER
!2014-08-28e

!            call wall_clock(d5)
!            t5 = d5 - d4
!            if(my_rank.eq.0)then
!            write(6,*)'density 5, t5=',t5
!            end if

end subroutine moments_gyro

end module amdgpu
