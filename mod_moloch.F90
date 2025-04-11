!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!    This file is part of ICTP RegCM.
!
!    ICTP RegCM is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    ICTP RegCM is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with ICTP RegCM.  If not, see <http://www.gnu.org/licenses/>
!
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

module mod_moloch

  use, intrinsic :: iso_fortran_env
  use mod_intkinds
  use mod_realkinds
  use mod_atmosphere
  use mod_constants
  use mod_dynparam
  use mod_runparams
  use mod_memutil
  use mod_mppparam
  use mod_zita

  implicit none

  private

  ! generalized vertical velocity
  real(rkx) , pointer , dimension(:,:,:) :: s => null( )
  ! nonhydrostatic term in pressure gradient force
  ! tridiagonal inversion
  real(rkx) , pointer , dimension(:,:,:) :: wx => null( )
  real(rkx) , pointer , dimension(:,:,:) :: deltaw => null( )
  real(rkx) , pointer , dimension(:,:,:) :: wz => null( )
  real(rkx) , pointer , dimension(:,:) :: mx2 => null( )
  real(rkx) , pointer , dimension(:,:) :: rmu => null( )
  real(rkx) , pointer , dimension(:,:) :: rmv => null( )
  real(rkx) , pointer , dimension(:,:,:) :: p0 => null( )
  real(rkx) , pointer , dimension(:,:,:) :: zdiv2 => null( )

  real(rkx) , dimension(:) , pointer :: gzitak => null( )
  real(rkx) , dimension(:) , pointer :: gzitakh => null( )
  real(rkx) , dimension(:) , pointer :: xknu => null( )
  real(rkx) , dimension(:,:) , pointer :: p2d => null( )
  real(rkx) , dimension(:,:) , pointer :: xlat => null( )
  real(rkx) , dimension(:,:) , pointer :: xlon => null( )
  real(rkx) , dimension(:,:) , pointer :: coru => null( )
  real(rkx) , dimension(:,:) , pointer :: corv => null( )
  real(rkx) , dimension(:,:) , pointer :: mu => null( )
  real(rkx) , dimension(:,:) , pointer :: hx => null( )
  real(rkx) , dimension(:,:) , pointer :: mx => null( )
  real(rkx) , dimension(:,:) , pointer :: mv => null( )
  real(rkx) , dimension(:,:) , pointer :: hy => null( )
  real(rkx) , dimension(:,:) , pointer :: ps => null( )
  real(rkx) , dimension(:,:) , pointer :: ht => null( )
  real(rkx) , dimension(:,:,:) , pointer :: fmz => null( )
  real(rkx) , dimension(:,:,:) , pointer :: fmzf => null( )
  real(rkx) , dimension(:,:,:) , pointer :: pai => null( )
  real(rkx) , dimension(:,:,:) , pointer :: pf => null( )
  real(rkx) , dimension(:,:,:) , pointer :: tetav => null( )
  real(rkx) , dimension(:,:,:) , pointer :: tf => null( )
  real(rkx) , dimension(:,:,:) , pointer :: tvirt => null( )
  real(rkx) , dimension(:,:,:) , pointer :: zeta => null( )
  real(rkx) , dimension(:,:,:) , pointer :: u => null( )
  real(rkx) , dimension(:,:,:) , pointer :: v => null( )
  real(rkx) , dimension(:,:,:) , pointer :: w => null( )
  real(rkx) , dimension(:,:,:) , pointer :: ux => null( )
  real(rkx) , dimension(:,:,:) , pointer :: vx => null( )
  real(rkx) , dimension(:,:,:) , pointer :: ud => null( )
  real(rkx) , dimension(:,:,:) , pointer :: vd => null( )
  real(rkx) , dimension(:,:,:) , pointer :: p => null( )
  real(rkx) , dimension(:,:,:) , pointer :: t => null( )
  real(rkx) , dimension(:,:,:) , pointer :: rho => null( )
  real(rkx) , dimension(:,:,:) , pointer :: qv => null( )
  real(rkx) , dimension(:,:,:) , pointer :: qc => null( )
  real(rkx) , dimension(:,:,:) , pointer :: qi => null( )
  real(rkx) , dimension(:,:,:) , pointer :: qr => null( )
  real(rkx) , dimension(:,:,:) , pointer :: qs => null( )
  real(rkx) , dimension(:,:,:) , pointer :: qsat => null( )
  real(rkx) , dimension(:,:,:) , pointer :: qwltot => null( )
  real(rkx) , dimension(:,:,:) , pointer :: qwitot => null( )
  real(rkx) , dimension(:,:,:,:) , pointer :: qx => null( )

  public :: allocate_moloch , init_moloch , moloch

  real(rkx) , parameter :: minden = 1.0e-30_rkx
  real(rkx) , parameter :: xdamp = 0.0625_rkx

  logical , parameter :: do_fulleq       = .true.
  logical , parameter :: do_divdamp      = .true.
  logical , parameter :: do_vadvtwice    = .true.
  logical , parameter :: do_filterpai    = .false.
  logical , parameter :: do_filtertheta  = .false.

  logical :: lrotllr

  real(rkx) , parameter :: nupaitq = 0.05_rkx

  real(rkx) :: dzita
  integer(ik4) :: jmin , jmax , imin , imax

  contains

  ! Computes saturation pressure
  ! Reference:  Polynomial approximations from:
  !             Piotr J. Flatau, et al.,1992:  Polynomial fits to saturation
  !             vapor pressure.  Journal of Applied Meteorology, 31, 1507-1513.
  !
!DIR$ ATTRIBUTES FORCEINLINE :: pfesat
  pure elemental real(rkx) function pfesat(t,p) result(es)
!$acc routine seq
    implicit none
    real(rkx) , intent(in) :: t , p ! Temperature (K) , Pressure (Pa)

    real(rk8) :: td , t_limit , esat
    !
    ! For water vapor (temperature range 0C-100C)
    !
    real(rk8) , parameter :: a0 =  0.611213476e+03_rk8
    real(rk8) , parameter :: a1 =  0.444007856e+02_rk8
    real(rk8) , parameter :: a2 =  0.143064234e+01_rk8
    real(rk8) , parameter :: a3 =  0.264461437e-01_rk8
    real(rk8) , parameter :: a4 =  0.305903558e-03_rk8
    real(rk8) , parameter :: a5 =  0.196237241e-05_rk8
    real(rk8) , parameter :: a6 =  0.892344772e-08_rk8
    real(rk8) , parameter :: a7 = -0.373208410e-10_rk8
    real(rk8) , parameter :: a8 =  0.209339997e-13_rk8
    !
    ! For ice (temperature range -75C-0C)
    !
    real(rk8) , parameter :: c0 =  0.611123516e+03_rk8
    real(rk8) , parameter :: c1 =  0.503109514e+02_rk8
    real(rk8) , parameter :: c2 =  0.188369801e+01_rk8
    real(rk8) , parameter :: c3 =  0.420547422e-01_rk8
    real(rk8) , parameter :: c4 =  0.614396778e-03_rk8
    real(rk8) , parameter :: c5 =  0.602780717e-05_rk8
    real(rk8) , parameter :: c6 =  0.387940929e-07_rk8
    real(rk8) , parameter :: c7 =  0.149436277e-09_rk8
    real(rk8) , parameter :: c8 =  0.262655803e-12_rk8

    t_limit = t - tzero
    if ( t_limit > 100.0_rk8 ) t_limit = 100.0_rk8
    if ( t_limit < -75.0_rk8 ) t_limit = -75.0_rk8
    td = t_limit
    if ( td >= 0.0_rk8 ) then
      esat = a0 + td*(a1 + td*(a2 + td*(a3 + td*(a4 &
         + td*(a5 + td*(a6 + td*(a7 + td*a8)))))))
    else
      esat = c0 + td*(c1 + td*(c2 + td*(c3 + td*(c4 &
         + td*(c5 + td*(c6 + td*(c7 + td*c8)))))))
    end if
    es = real(min(esat,0.15_rk8*p),rkx) ! pa
  end function pfesat

!DIR$ ATTRIBUTES FORCEINLINE :: pfwsat
  pure elemental real(rkx) function pfwsat(t,p,e) result(ws) ! In kg/kg
!$acc routine seq
    implicit none
    real(rkx) , intent(in) :: t             ! Temperature (K)
    real(rkx) , intent(in) :: p             ! Pressure (Pa)
    real(rkx) , intent(in) , optional :: e  ! Saturated vapor pressure (Pa)
    real(rkx) :: es ! , qs
    if ( present(e) ) then
      es = e
    else
      es = pfesat(t,p)
    end if
    ws = ep2 * (es / (p - es))
    ! Bolton 1980
    ! qs = ep2 * es / (p - (0.378_rkx * es))
    ! ws = qs * ( d_one - qs)
  end function pfwsat

  subroutine allocate_moloch
    implicit none
    integer(ik4) :: k
    call getmem1d(gzitak,1,kzp1,'moloch:gzitak')
    call getmem1d(gzitakh,1,kz,'moloch:gzitakh')
    call getmem2d(p2d,jdi1,jdi2,idi1,idi2,'moloch:p2d')
    call getmem3d(deltaw,jce1ga,jce2ga,ice1ga,ice2ga,1,kzp1,'moloch:dw')
    call getmem3d(s,jci1,jci2,ici1,ici2,1,kzp1,'moloch:s')
    call getmem3d(wx,jce1,jce2,ice1,ice2,1,kz,'moloch:wx')
    call getmem3d(zdiv2,jce1ga,jce2ga,ice1ga,ice2ga,1,kz,'moloch:zdiv2')
    call getmem3d(wz,jci1,jci2,ice1gb,ice2gb,1,kz,'moloch:wz')
    call getmem3d(p0,jce1gb,jce2gb,ici1,ici2,1,kz,'moloch:p0')
    call getmem2d(mx2,jde1,jde2,ide1,ide2,'moloch:mx2')
    call getmem2d(rmu,jde1ga,jde2ga,ide1,ide2,'moloch:rmu')
    call getmem2d(rmv,jde1,jde2,ide1ga,ide2ga,'moloch:rmv')
    call getmem2d(coru,jde1,jde2,ice1,ice2,'moloch:coru')
    call getmem2d(corv,jce1,jce2,ide1,ide2,'moloch:corv')
    call getmem3d(ud,jde1ga,jde2ga,ice1ga,ice2ga,1,kz,'moloch:ud')
    call getmem3d(vd,jce1ga,jce2ga,ide1ga,ide2ga,1,kz,'moloch:vd')
    if ( do_fulleq ) then
      call getmem3d(qwltot,jci1,jci2,ici1,ici2,1,kz,'moloch:qwltot')
      call getmem3d(qwitot,jci1,jci2,ici1,ici2,1,kz,'moloch:qwitot')
    end if
    call getmem1d(xknu,1,kz,'moloch:xknu')
    do concurrent ( k = 1:kz )
      xknu(k) = xdamp + &
        (1.0_rkx-xdamp) * sin(d_half*mathpi*(1.0_rkx-real(k-1,rkx)/kzm1))
    end do
    if ( do_filterpai ) then
      call getmem3d(pf,jce1,jce2,ice1,ice2,1,kz,'moloch:pf')
    end if
    if ( do_filtertheta ) then
      call getmem3d(tf,jce1,jce2,ice1,ice2,1,kz,'moloch:tf')
    end if
  end subroutine allocate_moloch

  subroutine init_moloch
    implicit none
    integer(ik4) :: i , j
    call assignpnt(mddom%msfu,mu)
    call assignpnt(mddom%msfv,mv)
    call assignpnt(mddom%msfx,mx)
    call assignpnt(mddom%hx,hx)
    call assignpnt(mddom%hy,hy)
    call assignpnt(mddom%xlat,xlat)
    call assignpnt(mddom%xlon,xlon)
    call assignpnt(mddom%ht,ht)
    call assignpnt(sfs%ps,ps)
    call assignpnt(mo_atm%fmz,fmz)
    call assignpnt(mo_atm%fmzf,fmzf)
    call assignpnt(mo_atm%pai,pai)
    call assignpnt(mo_atm%tetav,tetav)
    call assignpnt(mo_atm%u,u)
    call assignpnt(mo_atm%ux,ux)
    call assignpnt(mo_atm%v,v)
    call assignpnt(mo_atm%vx,vx)
    call assignpnt(mo_atm%w,w)
    call assignpnt(mo_atm%tvirt,tvirt)
    call assignpnt(mo_atm%zeta,zeta)
    call assignpnt(mo_atm%p,p)
    call assignpnt(mo_atm%t,t)
    call assignpnt(mo_atm%rho,rho)
    call assignpnt(mo_atm%qx,qx)
    call assignpnt(mo_atm%qs,qsat)
    call assignpnt(mo_atm%qx,qv,iqv)
    call assignpnt(mo_atm%qx,qc,iqc)
    call assignpnt(mo_atm%qx,qi,iqi)
    call assignpnt(mo_atm%qx,qr,iqr)
    call assignpnt(mo_atm%qx,qs,iqs)
    coru = eomeg2*sin(mddom%ulat(jde1:jde2,ice1:ice2)*degrad)
    corv = eomeg2*sin(mddom%vlat(jce1:jce2,ide1:ide2)*degrad)
    mx2 = mx * mx
    rmu = d_one/mu
    rmv = d_one/mv
    gzitak = gzita(zita,mo_ztop,mo_a0)
    gzitakh = gzita(zitah,mo_ztop,mo_a0)
    dzita = mo_dzita
    do concurrent ( j = jce1:jce2, i = ice1:ice2 )
      w(j,i,1) = d_zero
    end do
    lrotllr = (iproj == 'ROTLLR')
    jmin = jcross1 - 2
    jmax = jcross2 + 2
    imin = icross1 - 2
    imax = icross2 + 2
  end subroutine init_moloch
  !
  ! Moloch dynamical integration engine
  !
  subroutine moloch(istep, lprint)
    implicit none
    integer(ik4) , intent(in) :: istep
    logical , intent(in) :: lprint
    real(rkx) :: dtsound , dtstepa
    real(rkx) :: maxp , minp , pmax , pmin
    real(rkx) :: maxu , minu , umax , umin
    real(rkx) :: maxv , minv , vmax , vmin
    real(rkx) :: zdgz , lrt , tv
    !real(rk8) :: jday
    integer(ik4) :: i , j , k , nadv
    integer(ik4) :: iconvec

    dtstepa = dt / real(mo_nadv,rkx)
    dtsound = dtstepa / real(mo_nsound,rkx)
    iconvec = 0

    call reset_tendencies

    do concurrent ( j = jce1:jce2, i = ice1:ice2, k = 1:kz )
      p(j,i,k) = (pai(j,i,k)**cpovr) * p00
      qsat(j,i,k) = pfwsat(t(j,i,k),p(j,i,k))
    end do

    do concurrent ( j = jce1:jce2, i = ice1:ice2, k = 1:kz )
      tvirt(j,i,k) = t(j,i,k) * (d_one + ep1*qv(j,i,k) - &
                                 qc(j,i,k) - qi(j,i,k) - &
                                 qr(j,i,k) - qs(j,i,k))
    end do
    if ( do_fulleq ) then
      do concurrent ( j = jci1:jci2, i = ici1:ici2, k = 1:kz )
        qwltot(j,i,k) = qc(j,i,k) + qr(j,i,k)
        qwitot(j,i,k) = qi(j,i,k) + qs(j,i,k)
      end do
    end if

    do concurrent ( j = jce1:jce2, i = ice1:ice2, k = 1:kz )
      tetav(j,i,k) = tvirt(j,i,k)/pai(j,i,k)
    end do

    if ( do_filterpai ) then
      do concurrent ( j = jce1:jce2, i = ice1:ice2, k = 1:kz )
        pf(j,i,k) = pai(j,i,k)
      end do
    end if
    if ( do_fulleq ) then
      if ( do_filtertheta ) then
        do concurrent ( j = jce1:jce2, i = ice1:ice2, k = 1:kz )
          tf(j,i,k) = tetav(j,i,k)
        end do
      end if
    end if

    do nadv = 1 , mo_nadv

      call sound(dtsound)

      call advection(dtstepa)

    end do ! Advection loop

    if ( do_filterpai ) then
      do concurrent ( j = jce1:jce2, i = ice1:ice2, k = 1:kz )
        pai(j,i,k) = pai(j,i,k) - pf(j,i,k)
      end do
      call filtpai
      do concurrent ( j = jce1:jce2, i = ice1:ice2, k = 1:kz )
        pai(j,i,k) = pai(j,i,k) + pf(j,i,k)
      end do
    end if

    if ( do_fulleq ) then
      if ( do_filtertheta ) then
        do concurrent ( j = jce1:jce2, i = ice1:ice2, k = 1:kz )
          tetav(j,i,k) = tetav(j,i,k) - tf(j,i,k)
        end do
        call filttheta
        do concurrent ( j = jce1:jce2, i = ice1:ice2, k = 1:kz )
          tetav(j,i,k) = tetav(j,i,k) + tf(j,i,k)
        end do
      end if
    end if

    do concurrent ( j = jci1:jci2, i = ici1:ici2, k = 1:kz )
      tvirt(j,i,k) = tetav(j,i,k)*pai(j,i,k)
    end do

    do concurrent ( j = jci1:jci2, i = ici1:ici2, k = 1:kz )
      t(j,i,k) = tvirt(j,i,k) / (d_one + ep1*qv(j,i,k) - &
                     qc(j,i,k) - qi(j,i,k) - qr(j,i,k) - qs(j,i,k))
    end do

    do concurrent ( j = jce1:jce2, i = ice1:ice2, k = 1:kz )
      p(j,i,k) = (pai(j,i,k)**cpovr) * p00
      rho(j,i,k) = p(j,i,k)/(rgas*t(j,i,k))
    end do

    !jday = yeardayfrac(rcmtimer%idate)

#ifdef STDPAR
    do concurrent ( j = jce1:jce2, i = ice1:ice2 ) local(zdgz,lrt,tv)
#else
    do i = ice1 , ice2
      do j = jce1 , jce2
#endif
        zdgz = zeta(j,i,kz)*egrav
        lrt = (tvirt(j,i,kz-1)-tvirt(j,i,kz))/(zeta(j,i,kz-1)-zeta(j,i,kz))
        lrt = 0.65_rkx*lrt - 0.35_rkx*lrate
        tv = tvirt(j,i,kz) - 0.5_rkx*zeta(j,i,kz)*lrt ! Mean temperature
        ps(j,i) = p(j,i,kz) * exp(zdgz/(rgas*tv))
#ifndef STDPAR
      end do
#endif
    end do
    !
    ! Recompute saturation
    !
    do concurrent ( j = jce1:jce2, i = ice1:ice2, k = 1:kz )
      qsat(j,i,k) = pfwsat(t(j,i,k),p(j,i,k))
    end do
    !
    ! Diagnostic and end timestep
    !
    if ( lprint ) then
      maxp = maxval(ps(jci1:jci2,ici1:ici2))
      minp = minval(ps(jci1:jci2,ici1:ici2))
      maxu = maxval(u(jci1:jci2,ici1:ici2,1:kz))
      minu = minval(u(jci1:jci2,ici1:ici2,1:kz))
      maxv = maxval(v(jci1:jci2,ici1:ici2,1:kz))
      minv = minval(v(jci1:jci2,ici1:ici2,1:kz))
      call maxall(maxp,pmax)
      call minall(minp,pmin)
      call maxall(maxu,umax)
      call minall(minu,umin)
      call maxall(maxv,vmax)
      call minall(minv,vmin)
      if ( is_nan(pmax) .or. is_nan(pmin) .or. &
           is_inf(pmax) .or. is_inf(pmin) ) then
        write (error_unit,*) 'WHUUUUBBBASAAAGASDDWD!!!!!!!!!!!!!!!!'
        write (error_unit,*) 'No more atmosphere here....'
        write (error_unit,*) 'CFL violation detected, so model STOP'
        write (error_unit,*) '#####################################'
        write (error_unit,*) '#            DECREASE DT !!!!       #'
        write (error_unit,*) '#####################################'
        call fatal(__FILE__,__LINE__,'CFL VIOLATION')
      end if
      if ( myid == 0 ) then
        write(output_unit,*) '$$$  STEP IS : ', istep
        write(output_unit,'(a,2f8.2)') &
            ' $$$ max, min of ps (mb) = ', pmax*d_r100 , pmin*d_r100
        write(output_unit,'(a,2f8.2)') &
            ' $$$ max, min of u (m/s) = ', umax , umin
        write(output_unit,'(a,2f8.2)') &
            ' $$$ max, min of v (m/s) = ', vmax , vmin
      end if
    end if
    !
  end subroutine moloch

    subroutine divergence_filter( )
      implicit none
      integer(ik4) :: j , i , k
      call exchange_lrbt(zdiv2,1,jce1,jce2,ice1,ice2,1,kz)
      do k = 1 , kz
        do concurrent ( j = jci1:jci2, i = ici1:ici2 )
          p2d(j,i) = 0.125_rkx * (zdiv2(j-1,i,k) + zdiv2(j+1,i,k) + &
                                  zdiv2(j,i-1,k) + zdiv2(j,i+1,k)) - &
                       d_half   * zdiv2(j,i,k)
        end do
        do concurrent ( j = jci1:jci2, i = ici1:ici2 )
          zdiv2(j,i,k) = zdiv2(j,i,k) + mo_anu2 * xknu(k) * p2d(j,i)
        end do
      end do
    end subroutine divergence_filter

    subroutine filtpai
      implicit none
      integer(ik4) :: j , i , k

      call exchange_lrbt(pai,1,jce1,jce2,ice1,ice2,1,kz)

      do k = 1 , kz
        do concurrent ( j = jci1:jci2, i = ici1:ici2 )
          p2d(j,i) = 0.125_rkx * (pai(j-1,i,k) + pai(j+1,i,k) + &
                                  pai(j,i-1,k) + pai(j,i+1,k)) - &
                       d_half   * pai(j,i,k)
        end do
        do concurrent ( j = jci1:jci2, i = ici1:ici2 )
          pai(j,i,k) = pai(j,i,k) + nupaitq * p2d(j,i)
        end do
      end do
    end subroutine filtpai

    subroutine filttheta
      implicit none
      integer(ik4) :: j , i , k

      call exchange_lrbt(tetav,1,jce1,jce2,ice1,ice2,1,kz)

      do k = 1 , kz
        do concurrent ( j = jci1:jci2, i = ici1:ici2 )
          p2d(j,i) = 0.125_rkx * (tetav(j-1,i,k) + tetav(j+1,i,k) + &
                                  tetav(j,i-1,k) + tetav(j,i+1,k)) - &
                       d_half   * tetav(j,i,k)
        end do
        do concurrent ( j = jci1:jci2, i = ici1:ici2 )
          tetav(j,i,k) = tetav(j,i,k) + nupaitq * p2d(j,i)
        end do
      end do
    end subroutine filttheta

    subroutine sound(dts)
      implicit none
      real(rkx) , intent(in) :: dts
      integer(ik4) :: i , j , k , nsound
      real(rkx) :: dtrdx , dtrdy , dtrdz , zcs2
      real(rkx) , dimension(jci1:jci2,ici1:ici2,2:kzp1) :: wwkw
      real(rkx) :: zrfmzum , zrfmzvm , zrfmzup , zrfmzvp
      real(rkx) :: zum , zup , zvm , zvp , zuh , zvh
      real(rkx) :: zrom1w , zwexpl , zqs , zdth , zu , zd , zrapp
      real(rkx) :: zcx , zcy , zfz
      real(rkx) :: zrom1u , zcor1u , zrom1v , zcor1v

      dtrdx = dts/dx
      dtrdy = dts/dx
      dtrdz = dts/dzita
      zcs2 = dtrdz**2*rdrcv

      do concurrent ( j = jci1:jci2, i = ici1:ici2 )
        wwkw(j,i,kzp1) = d_zero
      end do

      !  sound waves

      if ( .not. do_fulleq ) then
        call exchange_lrbt(tetav,1,jce1,jce2,ice1,ice2,1,kz)
      end if

      do nsound = 1 , mo_nsound

        call exchange(u,1,jde1,jde2,ice1,ice2,1,kz)
        call exchange(v,1,jce1,jce2,ide1,ide2,1,kz)

        ! partial definition of the generalized vertical velocity

#ifdef STDPAR
        do concurrent ( j = jci1:jci2, i = ici1:ici2 ) local(zuh,zvh)
#else
        do i = ici1 , ici2
          do j = jci1 , jci2
#endif
            zuh = u(j,i,kz) * hx(j,i) + u(j+1,i,kz) * hx(j+1,i)
            zvh = v(j,i,kz) * hy(j,i) + v(j,i+1,kz) * hy(j,i+1)
            w(j,i,kzp1) = d_half * (zuh+zvh)
#ifndef STDPAR
          end do
#endif
        end do

        do concurrent ( j = jci1:jci2, i = ici1:ici2 )
          s(j,i,kzp1) = -w(j,i,kzp1)
        end do

        ! Equation 10, generalized vertical velocity

#ifdef STDPAR
        do concurrent ( j = jci1:jci2, i = ici1:ici2, &
                        k = 2:kz ) local(zuh,zvh)
#else
        do k = 2 , kz
          do i = ici1 , ici2
            do j = jci1 , jci2
#endif
              zuh = (u(j,i,k)   + u(j,i,k-1))   * hx(j,i) +    &
                    (u(j+1,i,k) + u(j+1,i,k-1)) * hx(j+1,i)
              zvh = (v(j,i,k)   + v(j,i,k-1))   * hy(j,i) +    &
                    (v(j,i+1,k) + v(j,i+1,k-1)) * hy(j,i+1)
              s(j,i,k) = -0.25_rkx * (zuh+zvh) * gzitak(k)
#ifndef STDPAR
            end do
          end do
#endif
        end do

        ! Part of divergence (except w contribution) put in zdiv2
        ! Equation 16

        if ( lrotllr ) then
#ifdef STDPAR
          do concurrent ( j = jci1:jci2, i = ici1:ici2, k = 1:kz ) &
            local(zrfmzum,zrfmzvm,zrfmzup,zrfmzvp,zum,zup,zvm,zvp)
#else
          do k = 1 , kz
            do i = ici1 , ici2
              do j = jci1 , jci2
#endif
                zrfmzum = d_two / (fmz(j,i,k) + fmz(j-1,i,k))
                zrfmzvm = d_two / (fmz(j,i,k) + fmz(j,i-1,k))
                zrfmzup = d_two / (fmz(j,i,k) + fmz(j+1,i,k))
                zrfmzvp = d_two / (fmz(j,i,k) + fmz(j,i+1,k))
                zum = dtrdx * u(j,i,k) * zrfmzum
                zup = dtrdx * u(j+1,i,k) * zrfmzup
                zvm = dtrdy * v(j,i,k) * zrfmzvm * rmv(j,i)
                zvp = dtrdy * v(j,i+1,k) * zrfmzvp * rmv(j,i+1)
                zdiv2(j,i,k) = fmz(j,i,k) * mx(j,i) * ((zup-zum) + (zvp-zvm))
#ifndef STDPAR
              end do
            end do
#endif
          end do
        else
#ifdef STDPAR
          do concurrent ( j = jci1:jci2, i = ici1:ici2, k = 1:kz ) &
            local(zrfmzum,zrfmzvm,zrfmzup,zrfmzvp,zum,zup,zvm,zvp)
#else
          do k = 1 , kz
            do i = ici1 , ici2
              do j = jci1 , jci2
#endif
                zrfmzum = d_two / (fmz(j,i,k) + fmz(j-1,i,k))
                zrfmzvm = d_two / (fmz(j,i,k) + fmz(j,i-1,k))
                zrfmzup = d_two / (fmz(j,i,k) + fmz(j+1,i,k))
                zrfmzvp = d_two / (fmz(j,i,k) + fmz(j,i+1,k))
                zum = dtrdx * u(j,i,k)   * rmu(j,i)   * zrfmzum
                zup = dtrdx * u(j+1,i,k) * rmu(j+1,i) * zrfmzup
                zvm = dtrdy * v(j,i,k)   * rmv(j,i)   * zrfmzvm
                zvp = dtrdy * v(j,i+1,k) * rmv(j,i+1) * zrfmzvp
                zdiv2(j,i,k) = mx2(j,i) * fmz(j,i,k) * ((zup-zum)+(zvp-zvm))
#ifndef STDPAR
              end do
            end do
#endif
          end do
        end if

        if ( do_divdamp ) then
          call divdamp(dts)
        end if

        do concurrent ( j = jci1:jci2, i = ici1:ici2, k = 1:kz )
          zdiv2(j,i,k) = zdiv2(j,i,k) + fmz(j,i,k) * &
                     dtrdz * (s(j,i,k) - s(j,i,k+1))
        end do

        ! new w (implicit scheme) from Equation 19

        do k = kz , 2 , -1
#ifdef STDPAR
          do concurrent ( j = jci1:jci2, i = ici1:ici2 ) &
            local(zrom1w,zwexpl,zqs,zdth,zu,zd,zrapp)
#else
          do i = ici1 , ici2
            do j = jci1 , jci2
#endif
              deltaw(j,i,k) = -w(j,i,k)
              ! explicit w:
              !    it must be consistent with the initialization of pai
              zrom1w = d_half * cpd * fmzf(j,i,k) * &
                       (tetav(j,i,k-1)+tetav(j,i,k))
              zrom1w = zrom1w - cpd * w(j,i,k) * &
                       fmzf(j,i,k)*fmzf(j,i,k) * &
                       real(nsound,rkx) * dtrdz * &
                       (tetav(j,i,k-1)-tetav(j,i,k)) !! GW
              if ( qv(j,i,k) > 0.96_rkx*qsat(j,i,k) .and. &
                    w(j,i,k) > 0.1_rkx ) then
                zqs = d_half*(qsat(j,i,k)+qsat(j,i,k-1))
                zdth = egrav*w(j,i,k)*real(nsound-1,rkx)*dts*wlhv*wlhv* &
                  zqs/(cpd*pai(j,i,k-1)*rwat*t(j,i,k-1)*t(j,i,k-1))
                zrom1w = zrom1w + zdth*fmzf(j,i,k)
              end if
              ! explicit w:
              !    it must be consistent with the initialization of pai
              zwexpl = w(j,i,k) - zrom1w * dtrdz * &
                       (pai(j,i,k-1) - pai(j,i,k)) - egrav*dts
              zwexpl = zwexpl + rdrcv * zrom1w * dtrdz * &
                       (pai(j,i,k-1) * zdiv2(j,i,k-1) - &
                        pai(j,i,k)   * zdiv2(j,i,k))
              ! computation of the tridiagonal matrix coefficients
              ! -zu*w(k+1) + (1+zu+zd)*w(k) - zd*w(k-1) = zwexpl
              zu = zcs2 * fmz(j,i,k-1) * zrom1w * pai(j,i,k-1) + ffilt(k)
              zd = zcs2 * fmz(j,i,k)   * zrom1w * pai(j,i,k)   + ffilt(k)
              ! 1st loop for the tridiagonal inversion
              ! a = -zd ; b = (1+zu+zd) ; c = -zu
              zrapp = d_one / (d_one + zd + zu - zd*wwkw(j,i,k+1))
              w(j,i,k) = zrapp * (zwexpl + zd * w(j,i,k+1))
              wwkw(j,i,k) = zrapp * zu
#ifndef STDPAR
            end do
#endif
          end do
        end do

        ! 2nd loop for the tridiagonal inversion
        do k = 2 , kz
          do concurrent ( j = jci1:jci2, i = ici1:ici2 )
            w(j,i,k) = w(j,i,k) + wwkw(j,i,k)*w(j,i,k-1)
            deltaw(j,i,k) = deltaw(j,i,k) + w(j,i,k)
          end do
        end do

        ! new Exner function (Equation 19)

        do concurrent ( j = jci1:jci2, i = ici1:ici2, k = 1:kz )
          zdiv2(j,i,k) = zdiv2(j,i,k) + dtrdz * fmz(j,i,k) * &
                    (w(j,i,k) - w(j,i,k+1))
        end do

        if ( do_fulleq ) then
          do concurrent ( j = jci1:jci2, i = ici1:ici2, k = 1:kz )
            zdiv2(j,i,k) = zdiv2(j,i,k) * &
                   (d_one + 0.86_rkx * qv(j,i,k) + &
                            3.2_rkx * qc(j,i,k)) / &
                   (d_one + 0.96_rkx * qv(j,i,k) + &
                            4.8_rkx * qc(j,i,k))
            tetav(j,i,k) = tetav(j,i,k) * &
                   (d_one + rdrcv*zdiv2(j,i,k) * &
                    (0.25_rkx * qv(j,i,k) +      &
                     4.2_rkx * qwltot(j,i,k) +   &
                     2.1_rkx * qwitot(j,i,k)))
          end do
          call exchange_lrbt(tetav,1,jce1,jce2,ice1,ice2,1,kz)
        end if

        if ( mo_divfilter ) call divergence_filter( )

        ! horizontal momentum equations
        do concurrent ( j = jci1:jci2, i = ici1:ici2, k = 1:kz )
          pai(j,i,k) = pai(j,i,k) * (d_one - rdrcv*zdiv2(j,i,k))
        end do

        call exchange_lrbt(pai,1,jce1,jce2,ice1,ice2,1,kz)
        call exchange_lrbt(deltaw,1,jce1,jce2,ice1,ice2,1,kzp1)

        do concurrent ( j = jde1ga:jde2ga , i = ice1ga:ice2ga , k = 1:kz )
          ud(j,i,k) = u(j,i,k)
        end do
        do concurrent ( j = jce1ga:jce2ga , i = ide1ga:ide2ga , k = 1:kz )
          vd(j,i,k) = v(j,i,k)
        end do

        if ( lrotllr ) then
          ! Equation 17
#ifdef STDPAR
          do concurrent ( j = jdi1:jdi2, i = ici1:ici2, k = 1:kz ) &
            local(zcx,zfz,zrom1u,zcor1u)
#else
          do k = 1 , kz
            do i = ici1 , ici2
              do j = jci1 , jci2
#endif
                zcx = dtrdx * mu(j,i)
                zfz = egrav * dts + 0.25_rkx * &
                    (deltaw(j-1,i,k) + deltaw(j-1,i,k+1) + &
                     deltaw(j,i,k)   + deltaw(j,i,k+1))
                zrom1u = d_half * cpd * (tetav(j-1,i,k) + tetav(j,i,k))
                zcor1u = coru(j,i) * dts * 0.25_rkx * &
                     (vd(j,i,k) + vd(j-1,i,k) + vd(j-1,i+1,k) + vd(j,i+1,k))
                ! Equation 17
                u(j,i,k) = u(j,i,k) + zcor1u - &
                           zfz * hx(j,i) * gzitakh(k) - &
                           zcx * zrom1u * (pai(j,i,k) - pai(j-1,i,k))
#ifndef STDPAR
              end do
            end do
#endif
          end do
          ! Equation 18
#ifdef STDPAR
          do concurrent ( j = jci1:jci2, i = idi1:idi2, k = 1:kz ) &
            local(zcy,zfz,zrom1v,zcor1v)
#else
          do k = 1 , kz
            do i = idi1 , idi2
              do j = jci1 , jci2
#endif
                zcy = dtrdy
                zfz = egrav * dts + 0.25_rkx * &
                    (deltaw(j,i-1,k) + deltaw(j,i-1,k+1) + &
                     deltaw(j,i,k)   + deltaw(j,i,k+1))
                zrom1v = d_half * cpd * (tetav(j,i-1,k) + tetav(j,i,k))
                zcor1v = corv(j,i) * dts * 0.25_rkx * &
                       (ud(j,i,k) + ud(j,i-1,k) + ud(j+1,i,k) + ud(j+1,i-1,k))
                ! Equation 18
                v(j,i,k) = v(j,i,k) - zcor1v - &
                           zfz * hy(j,i) * gzitakh(k) -  &
                           zcy * zrom1v * (pai(j,i,k) - pai(j,i-1,k))
#ifndef STDPAR
              end do
            end do
#endif
          end do
        else
#ifdef STDPAR
          do concurrent ( j = jdi1:jdi2, i = ici1:ici2, k = 1:kz ) &
            local(zcx,zfz,zrom1u,zcor1u)
#else
          do k = 1 , kz
            do i = ici1 , ici2
              do j = jci1 , jci2
#endif
                zcx = dtrdx * mu(j,i)
                zfz = egrav * dts + 0.25_rkx * &
                    (deltaw(j-1,i,k) + deltaw(j-1,i,k+1) + &
                     deltaw(j,i,k)   + deltaw(j,i,k+1))
                zrom1u = d_half * cpd * (tetav(j-1,i,k) + tetav(j,i,k))
                zcor1u = coru(j,i) * dts * 0.25_rkx * &
                     (vd(j,i,k) + vd(j-1,i,k) + vd(j-1,i+1,k) + vd(j,i+1,k))
                ! Equation 17
                u(j,i,k) = u(j,i,k) + zcor1u - &
                           zfz * hx(j,i) * gzitakh(k) - &
                           zcx * zrom1u * (pai(j,i,k) - pai(j-1,i,k))
#ifndef STDPAR
              end do
            end do
#endif
          end do
#ifdef STDPAR
          do concurrent ( j = jci1:jci2, i = idi1:idi2, k = 1:kz ) &
            local(zcy,zfz,zrom1v,zcor1v)
#else
          do k = 1 , kz
            do i = idi1 , idi2
              do j = jci1 , jci2
#endif
                zcy = dtrdy * mv(j,i)
                zfz = egrav * dts + 0.25_rkx * &
                    (deltaw(j,i-1,k) + deltaw(j,i-1,k+1) + &
                     deltaw(j,i,k)   + deltaw(j,i,k+1))
                zrom1v = d_half * cpd * (tetav(j,i-1,k) + tetav(j,i,k))
                zcor1v = corv(j,i) * dts * 0.25_rkx * &
                       (ud(j,i,k) + ud(j,i-1,k) + ud(j+1,i,k) + ud(j+1,i-1,k))
                ! Equation 18
                v(j,i,k) = v(j,i,k) - zcor1v - &
                           zfz * hy(j,i) * gzitakh(k) -  &
                           zcy * zrom1v * (pai(j,i,k) - pai(j,i-1,k))
#ifndef STDPAR
              end do
            end do
#endif
          end do
        end if

      end do ! sound loop

      ! complete computation of generalized vertical velocity
      ! Complete Equation 10
      do concurrent ( j = jci1:jci2, i = ici1:ici2, k = 2:kz )
        s(j,i,k) = (w(j,i,k) + s(j,i,k)) * fmzf(j,i,k)
      end do
      do concurrent ( j = jci1:jci2, i = ici1:ici2 )
        s(j,i,1) = d_zero
        s(j,i,kzp1) = d_zero
      end do

    end subroutine sound

    subroutine divdamp(dts)
      implicit none
      real(rkx) , intent(in) :: dts
      integer(ik4) :: i , j , k
      real(rkx) :: ddamp

      call exchange_lrbt(zdiv2,1,jce1,jce2,ice1,ice2,1,kz)

      ddamp = 0.125_rkx * (dx/dts)
      if ( lrotllr ) then
        do concurrent ( j = jdi1:jdi2, i = ici1:ici2, k = 1:kz )
          u(j,i,k) = u(j,i,k) + &
                  xknu(k)*ddamp*mu(j,i)*(zdiv2(j,i,k)-zdiv2(j-1,i,k))
        end do
        do concurrent ( j = jci1:jci2, i = idi1:idi2, k = 1:kz )
          v(j,i,k) = v(j,i,k) + &
                  xknu(k)*ddamp*(zdiv2(j,i,k)-zdiv2(j,i-1,k))
        end do
      else
        do concurrent ( j = jdi1:jdi2, i = ici1:ici2, k = 1:kz )
          u(j,i,k) = u(j,i,k) + &
                  xknu(k)*ddamp*mu(j,i)*(zdiv2(j,i,k)-zdiv2(j-1,i,k))
        end do
        do concurrent ( j = jci1:jci2, i = idi1:idi2, k = 1:kz )
          v(j,i,k) = v(j,i,k) + &
                  xknu(k)*ddamp*mv(j,i)*(zdiv2(j,i,k)-zdiv2(j,i-1,k))
        end do
      end if
      ! Horizontal diffusion
      ddamp = xdamp * 0.015625_rkx/dts
      do k = 1 , kz
        do concurrent ( j = jci1:jci2, i = ici1:ici2 )
          p2d(j,i) = 0.125_rkx * (zdiv2(j-1,i,k) + zdiv2(j+1,i,k) + &
                                  zdiv2(j,i-1,k) + zdiv2(j,i+1,k)) - &
                       d_half   * zdiv2(j,i,k)
        end do
        do concurrent ( j = jci1:jci2, i = ici1:ici2 )
          zdiv2(j,i,k) = zdiv2(j,i,k) + ddamp * p2d(j,i)
        end do
      end do
    end subroutine divdamp

    subroutine advection(dta)
      implicit none
      real(rkx) , intent(in) :: dta
      integer(ik4) :: n
      real(rkx) , pointer , dimension(:,:,:) :: ptr => null( )

      ! Compute U,V on cross points

      call uvstagtox(u,v,ux,vx)

      ! Compute W (and TKE if required) on zita levels

      call wstagtox(w,wx)

      call wafone(tetav,dta)
      call wafone(pai,dta)
      call wafone(ux,dta)
      call wafone(vx,dta)
      call wafone(wx,dta)
      call wafone(qv,dta)
      do n = iqfrst , nqx
        call assignpnt(qx,ptr,n)
        call wafone(ptr,dta)
      end do

      ! Interpolate on staggered points
      call xtouvstag(ux,vx,u,v)

      ! Back to half-levels
      call xtowstag(wx,w)
    end subroutine advection

    subroutine wafone(pp,dta,pfac,pmin)
      implicit none
      real(rkx) , dimension(:,:,:) , pointer , intent(inout) :: pp
      real(rkx) , intent(in) :: dta
      real(rkx) , optional , intent(in) :: pfac , pmin
      integer(ik4) :: j , i , k
      real(rkx) :: dtrdx , dtrdy , dtrdz
      real(rkx) , parameter :: wlow  = 0.0_rkx
      real(rkx) , parameter :: whigh = 2.0_rkx
      real(rkx) , dimension(jci1:jci2,ici1:ici2,1:kzp1) :: wfw
      real(rkx) , dimension(jci1:jci2,ici1:ice2ga,1:kz) :: zpby
      real(rkx) , dimension(jci1:jce2ga,ici1:ici2,1:kz) :: zpbw
      real(rkx) :: zamu , is , r , b , zphi , zzden , zdv
      real(rkx) :: zhxvtn , zhxvts , zcostx
      real(rkx) :: zrfmu , zrfmd
      real(rkx) :: zrfmn , zrfms
      real(rkx) :: zrfme , zrfmw
      real(rkx) :: zfac , zmin
      integer(ik4) :: k1 , k1p1
      integer(ik4) :: ih , ihm1
      integer(ik4) :: jh , jhm1

      dtrdx = dta/dx
      dtrdy = dta/dx
      dtrdz = dta/dzita
      if ( do_vadvtwice ) then
        dtrdz = 0.5_rkx * dtrdz
      end if

      if ( present(pfac) ) then
        zfac = pfac
        do concurrent ( j = jce1:jce2, i = ice1:ice2, k = 1:kz )
          pp(j,i,k) = pp(j,i,k) * zfac
        end do
      end if

      ! Vertical advection
      do concurrent ( j = jci1:jci2 , i = ici1:ici2 )
        wfw(j,i,1) = d_zero
        wfw(j,i,kzp1) = d_zero
      end do

#ifdef STDPAR
      do concurrent ( j = jci1:jci2, i = ici1:ici2, k = 1:kzm1 ) &
        local(zamu,is,r,b,zphi,zzden,k1,k1p1)
#else
      do k = 1 , kzm1
        do i = ici1 , ici2
          do j = jci1 , jci2
#endif
            zamu = s(j,i,k+1) * dtrdz
            if ( zamu >= d_zero ) then
              is = d_one
              k1 = k + 1
              k1p1 = k1 + 1
              if ( k1p1 > kz ) k1p1 = kz
            else
              is = -d_one
              k1 = k - 1
              k1p1 = k
              if ( k1 < 1 ) k1 = 1
            end if
            zzden = pp(j,i,k)-pp(j,i,k+1)
            zzden = sign(max(abs(zzden),minden),zzden)
            r = (pp(j,i,k1)-pp(j,i,k1p1))/zzden
            b = max(wlow, min(whigh, max(r, min(d_two*r,d_one))))
            zphi = is + zamu * b - is * b
            wfw(j,i,k+1) = d_half * s(j,i,k+1) * ((d_one+zphi)*pp(j,i,k+1) + &
                                                  (d_one-zphi)*pp(j,i,k))
#ifndef STDPAR
          end do
        end do
#endif
      end do
#ifdef STDPAR
      do concurrent ( j = jci1:jci2, i = ici1:ici2, k = 1:kz ) &
        local(zrfmu,zrfmd,zdv)
#else
      do k = 1 , kz
        do i = ici1 , ici2
          do j = jci1 , jci2
#endif
            zrfmu = dtrdz * fmz(j,i,k)/fmzf(j,i,k)
            zrfmd = dtrdz * fmz(j,i,k)/fmzf(j,i,k+1)
            zdv = (s(j,i,k)*zrfmu - s(j,i,k+1)*zrfmd) * pp(j,i,k)
            wz(j,i,k) = pp(j,i,k) - &
              wfw(j,i,k)*zrfmu + wfw(j,i,k+1)*zrfmd + zdv
#ifndef STDPAR
          end do
        end do
#endif
      end do

      if ( do_vadvtwice ) then
#ifdef STDPAR
        do concurrent ( j = jci1:jci2, i = ici1:ici2, k = 1:kzm1 ) &
          local(zamu,is,r,b,zphi,zzden,k1,k1p1)
#else
        do k = 1 , kzm1
          do i = ici1 , ici2
            do j = jci1 , jci2
#endif
              zamu = s(j,i,k+1) * dtrdz
              if ( zamu >= d_zero ) then
                is = d_one
                k1 = k + 1
                k1p1 = k1 + 1
                if ( k1p1 > kz ) k1p1 = kz
              else
                is = -d_one
                k1 = k - 1
                k1p1 = k
                if ( k1 < 1 ) k1 = 1
              end if
              zzden = wz(j,i,k)-wz(j,i,k+1)
              zzden = sign(max(abs(zzden),minden),zzden)
              r = (wz(j,i,k1)-wz(j,i,k1p1))/zzden
              b = max(wlow, min(whigh, max(r, min(d_two*r,d_one))))
              zphi = is + zamu * b - is * b
              wfw(j,i,k+1) = d_half * s(j,i,k+1) * &
                ((d_one+zphi)*wz(j,i,k+1) + (d_one-zphi)*wz(j,i,k))
#ifndef STDPAR
            end do
          end do
#endif
        end do
#ifdef STDPAR
        do concurrent ( j = jci1:jci2, i = ici1:ici2, k = 1:kz ) &
          local(zrfmu,zrfmd,zdv)
#else
        do k = 1 , kz
          do i = ici1 , ici2
            do j = jci1 , jci2
#endif
              zrfmu = dtrdz * fmz(j,i,k)/fmzf(j,i,k)
              zrfmd = dtrdz * fmz(j,i,k)/fmzf(j,i,k+1)
              zdv = (s(j,i,k)*zrfmu - s(j,i,k+1)*zrfmd) * wz(j,i,k)
              wz(j,i,k) = wz(j,i,k) - wfw(j,i,k)*zrfmu + &
                                      wfw(j,i,k+1)*zrfmd + zdv
#ifndef STDPAR
            end do
          end do
#endif
        end do
      end if

      call exchange_bt(wz,2,jci1,jci2,ice1,ice2,1,kz)

      if ( lrotllr ) then

        ! Meridional advection
#ifdef STDPAR
        do concurrent ( j = jci1:jci2, i = ici1:ice2ga, k = 1:kz ) &
          local(zamu,is,r,b,zphi,zzden,ih,ihm1)
#else
        do k = 1 , kz
          do i = ici1 , ice2ga
            do j = jci1 , jci2
#endif
              zamu = v(j,i,k) * dtrdy
              if ( zamu > d_zero ) then
                is = d_one
                ih = i-1
              else
                is = -d_one
                ih = min(i+1,imax)
              end if
              ihm1 = max(ih-1,imin)
              zzden = wz(j,i,k)-wz(j,i-1,k)
              zzden = sign(max(abs(zzden),minden),zzden)
              r = (wz(j,ih,k)-wz(j,ihm1,k))/zzden
              b = max(wlow, min(whigh, max(r, min(d_two*r,d_one))))
              zphi = is + zamu*b - is*b
              zpby(j,i,k) = d_half * v(j,i,k) * &
                  ((d_one+zphi)*wz(j,i-1,k) + (d_one-zphi)*wz(j,i,k))
#ifndef STDPAR
            end do
          end do
#endif
        end do
#ifdef STDPAR
        do concurrent ( j = jci1:jci2, i = ici1:ici2, k = 1:kz ) &
          local(zhxvtn,zhxvts,zrfmn,zrfms,zdv)
#else
        do k = 1 , kz
          do i = ici1 , ici2
            do j = jci1 , jci2
#endif
              zhxvtn = dtrdy * rmv(j,i+1) * mx(j,i)
              zhxvts = dtrdy * rmv(j,i) * mx(j,i)
              zrfmn = zhxvtn * d_two * fmz(j,i,k)/(fmz(j,i,k)+fmz(j,i+1,k))
              zrfms = zhxvts * d_two * fmz(j,i,k)/(fmz(j,i,k)+fmz(j,i-1,k))
              zdv = (v(j,i+1,k) * zrfmn - v(j,i,k) * zrfms) * pp(j,i,k)
              p0(j,i,k) = wz(j,i,k) + &
                    zpby(j,i,k)*zrfms - zpby(j,i+1,k)*zrfmn + zdv
#ifndef STDPAR
            end do
          end do
#endif
        end do

        call exchange_lr(p0,2,jce1,jce2,ici1,ici2,1,kz)

        ! Zonal advection

#ifdef STDPAR
        do concurrent ( j = jci1:jce2ga, i = ici1:ici2, k = 1:kz ) &
          local(zamu,is,r,b,zphi,zzden,jh,jhm1)
#else
        do k = 1 , kz
          do i = ici1 , ici2
            do j = jci1 , jce2ga
#endif
              zamu = u(j,i,k) * mu(j,i) * dtrdx
              if ( zamu > d_zero ) then
                is = d_one
                jh = j-1
              else
                is = -d_one
                jh = min(j+1,jmax)
              end if
              jhm1 = max(jh-1,jmin)
              zzden = p0(j,i,k)-p0(j-1,i,k)
              zzden = sign(max(abs(zzden),minden),zzden)
              r = (p0(jh,i,k)-p0(jhm1,i,k))/zzden
              b = max(wlow, min(whigh, max(r, min(d_two*r,d_one))))
              zphi = is + zamu*b - is*b
              zpbw(j,i,k) = d_half * u(j,i,k) * &
                     ((d_one+zphi)*p0(j-1,i,k) + (d_one-zphi)*p0(j,i,k))
#ifndef STDPAR
            end do
          end do
#endif
        end do
#ifdef STDPAR
        do concurrent ( j = jci1:jci2, i = ici1:ici2, k = 1:kz ) &
          local(zcostx,zrfme,zrfmw,zdv)
#else
        do k = 1 , kz
          do i = ici1 , ici2
            do j = jci1 , jci2
#endif
              zcostx = dtrdx * mu(j,i)
              zrfme = zcostx * d_two * fmz(j,i,k)/(fmz(j,i,k)+fmz(j+1,i,k))
              zrfmw = zcostx * d_two * fmz(j,i,k)/(fmz(j,i,k)+fmz(j-1,i,k))
              zdv = (u(j+1,i,k) * zrfme - u(j,i,k) * zrfmw) * pp(j,i,k)
              pp(j,i,k) = p0(j,i,k) + &
                     zpbw(j,i,k)*zrfmw - zpbw(j+1,i,k)*zrfme + zdv
#ifndef STDPAR
            end do
          end do
#endif
        end do

      else

        ! Meridional advection

#ifdef STDPAR
        do concurrent ( j = jci1:jci2, i = ici1:ice2ga, k = 1:kz ) &
          local(zamu,is,r,b,zphi,zzden,ihm1,ih)
#else
        do k = 1 , kz
          do i = ici1 , ice2ga
            do j = jci1 , jci2
#endif
              zamu = v(j,i,k) * rmv(j,i) * dtrdy
              if ( zamu > d_zero ) then
                is = d_one
                ih = i-1
              else
                is = -d_one
                ih = min(i+1,imax)
              end if
              ihm1 = max(ih-1,imin)
              zzden = wz(j,i,k)-wz(j,i-1,k)
              zzden = sign(max(abs(zzden),minden),zzden)
              r = (wz(j,ih,k)-wz(j,ihm1,k))/zzden
              b = max(wlow, min(whigh, max(r, min(d_two*r,d_one))))
              zphi = is + zamu*b - is*b
              zpby(j,i,k) = d_half * v(j,i,k) * rmv(j,i) * &
                  ((d_one+zphi)*wz(j,i-1,k) + (d_one-zphi)*wz(j,i,k))
#ifndef STDPAR
            end do
          end do
#endif
        end do
#ifdef STDPAR
        do concurrent ( j = jci1:jci2, i = ici1:ici2, k = 1:kz ) &
          local(zrfmn,zrfms,zdv)
#else
        do k = 1 , kz
          do i = ici1 , ici2
            do j = jci1 , jci2
#endif
              zrfmn = dtrdy * d_two * fmz(j,i,k)/(fmz(j,i,k)+fmz(j,i+1,k))
              zrfms = dtrdy * d_two * fmz(j,i,k)/(fmz(j,i,k)+fmz(j,i-1,k))
              zdv = (v(j,i+1,k) * rmv(j,i+1) * zrfmn - &
                     v(j,i,k)   * rmv(j,i)   * zrfms) * pp(j,i,k)
              p0(j,i,k) = wz(j,i,k) + &
                mx2(j,i) * (zpby(j,i,k)*zrfms - zpby(j,i+1,k)*zrfmn + zdv)
#ifndef STDPAR
            end do
          end do
#endif
        end do

        call exchange_lr(p0,2,jce1,jce2,ici1,ici2,1,kz)

        ! Zonal advection

#ifdef STDPAR
        do concurrent ( j = jci1:jce2ga, i = ici1:ici2, k = 1:kz ) &
          local(zamu,is,r,b,zphi,zzden,jh,jhm1)
#else
        do k = 1 , kz
          do i = ici1 , ici2
            do j = jci1 , jce2ga
#endif
              zamu = u(j,i,k) * rmu(j,i) * dtrdx
              if ( zamu > d_zero ) then
                is = d_one
                jh = j-1
              else
                is = -d_one
                jh = min(j+1,jmax)
              end if
              jhm1 = max(jh-1,jmin)
              zzden = p0(j,i,k)-p0(j-1,i,k)
              zzden = sign(max(abs(zzden),minden),zzden)
              r = (p0(jh,i,k)-p0(jhm1,i,k))/zzden
              b = max(wlow, min(whigh, max(r, min(d_two*r,d_one))))
              zphi = is + zamu*b - is*b
              zpbw(j,i,k) = d_half * u(j,i,k) * rmu(j,i) * &
                     ((d_one+zphi)*p0(j-1,i,k) + (d_one-zphi)*p0(j,i,k))
#ifndef STDPAR
            end do
          end do
#endif
        end do

#ifdef STDPAR
        do concurrent ( j = jci1:jci2, i = ici1:ici2, k = 1:kz ) &
          local(zrfme,zrfmw,zdv)
#else
        do k = 1 , kz
          do i = ici1 , ici2
            do j = jci1 , jci2
#endif
              zrfme = dtrdx * d_two * fmz(j,i,k)/(fmz(j,i,k)+fmz(j+1,i,k))
              zrfmw = dtrdx * d_two * fmz(j,i,k)/(fmz(j,i,k)+fmz(j-1,i,k))
              zdv = (u(j+1,i,k) * rmu(j+1,i) * zrfme - &
                     u(j,i,k)   * rmu(j,i)   * zrfmw) * pp(j,i,k)
              pp(j,i,k) = p0(j,i,k) + &
                  mx2(j,i) * (zpbw(j,i,k)*zrfmw - zpbw(j+1,i,k)*zrfme + zdv)
#ifndef STDPAR
            end do
          end do
#endif
        end do
      end if

      if ( present(pfac) ) then
        do concurrent ( j = jce1:jce2, i = ice1:ice2, k = 1:kz )
          pp(j,i,k) = pp(j,i,k) / zfac
        end do
      end if
      if ( present(pmin) ) then
        zmin = pmin
        do concurrent ( j = jce1:jce2, i = ice1:ice2, k = 1:kz )
          pp(j,i,k) = max(pp(j,i,k),zmin)
        end do
      end if
    end subroutine wafone

    subroutine reset_tendencies
      implicit none
      integer(ik4) :: i , j , k
      do concurrent ( j = jci1:jci2, i = ici1:ici2, k = 1:kzp1 )
        s(j,i,k) = d_zero
      end do
      do concurrent ( j = jce1ga:jce2ga, i = ice1ga:ice2ga, k = 1:kz )
        zdiv2(j,i,k) = d_zero
      end do
      do concurrent ( j = jce1ga:jce2ga, i = ice1ga:ice2ga, k = 1:kzp1 )
        deltaw(j,i,k) = d_zero
      end do
    end subroutine reset_tendencies

  subroutine wstagtox(w,wx)
    implicit none
    real(rkx) , intent(in) , dimension(:,:,:) , pointer :: w
    real(rkx) , intent(inout) , dimension(:,:,:) , pointer :: wx
    integer(ik4) :: i , j , k

    do concurrent ( j = jce1:jce2, i = ice1:ice2, k = 2:kzm1 )
      wx(j,i,k) = 0.5625_rkx * (w(j,i,k+1)+w(j,i,k)) - &
                  0.0625_rkx * (w(j,i,k+2)+w(j,i,k-1))
    end do
    do concurrent ( j = jce1:jce2, i = ice1:ice2 )
      wx(j,i,1)  = d_half * (w(j,i,2)+w(j,i,1))
      wx(j,i,kz) = d_half * (w(j,i,kzp1)+w(j,i,kz))
    end do
  end subroutine wstagtox

  subroutine xtowstag(wx,w)
    implicit none
    real(rkx) , intent(in) , dimension(:,:,:) , pointer :: wx
    real(rkx) , intent(inout) , dimension(:,:,:) , pointer :: w
    integer(ik4) :: i , j , k

    do concurrent ( j = jce1:jce2, i = ice1:ice2, k = 3:kzm1 )
      w(j,i,k) = 0.5625_rkx * (wx(j,i,k)  +wx(j,i,k-1)) - &
                 0.0625_rkx * (wx(j,i,k+1)+wx(j,i,k-2))
    end do
    do concurrent ( j = jce1:jce2, i = ice1:ice2 )
      w(j,i,2) = d_half * (wx(j,i,2)  +wx(j,i,1))
      w(j,i,kz) = d_half * (wx(j,i,kz)+wx(j,i,kzm1))
    end do
  end subroutine xtowstag

  subroutine xtouvstag(ux,vx,u,v)
    implicit none
    real(rkx) , intent(inout) , dimension(:,:,:) , pointer :: ux , vx
    real(rkx) , intent(inout) , dimension(:,:,:) , pointer :: u , v
    integer(ik4) :: i , j , k

    call exchange_lr(ux,2,jce1,jce2,ice1,ice2,1,kz)
    call exchange_bt(vx,2,jce1,jce2,ice1,ice2,1,kz)

    ! Back to wind points: U (fourth order)

    do concurrent ( j = jdii1:jdii2, i = ici1:ici2, k = 1:kz )
      u(j,i,k) = 0.5625_rkx * (ux(j,i,k)  +ux(j-1,i,k)) - &
                 0.0625_rkx * (ux(j+1,i,k)+ux(j-2,i,k))
    end do

    ! Back to wind points: V (fourth order)

    do concurrent ( j = jci1:jci2, i = idii1:idii2, k = 1:kz )
      v(j,i,k) = 0.5625_rkx * (vx(j,i,k)  +vx(j,i-1,k)) - &
                 0.0625_rkx * (vx(j,i+1,k)+vx(j,i-2,k))
    end do
  end subroutine xtouvstag

  subroutine uvstagtox(u,v,ux,vx)
    implicit none
    real(rkx) , intent(inout) , dimension(:,:,:) , pointer :: u , v
    real(rkx) , intent(inout) , dimension(:,:,:) , pointer :: ux , vx
    integer(ik4) :: i , j , k

    call exchange_lr(u,2,jde1,jde2,ice1,ice2,1,kz)
    call exchange_bt(v,2,jce1,jce2,ide1,ide2,1,kz)

    ! Compute U-wind on T points

    do concurrent ( j = jci1:jci2, i = ice1:ice2, k = 1:kz )
      ux(j,i,k) = 0.5625_rkx * (u(j+1,i,k)+u(j,i,k)) - &
                  0.0625_rkx * (u(j+2,i,k)+u(j-1,i,k))
    end do

    ! Compute V-wind on T points

    do concurrent ( j = jce1:jce2, i = ici1:ici2, k = 1:kz )
      vx(j,i,k) = 0.5625_rkx * (v(j,i+1,k)+v(j,i,k)) - &
                  0.0625_rkx * (v(j,i+2,k)+v(j,i-1,k))
    end do
  end subroutine uvstagtox

end module mod_moloch

! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
