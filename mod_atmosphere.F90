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
!    along with ICTP RegCM.  If not, see <http://www.gnu.org/licenses/>.
!
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

module mod_atmosphere

  use mod_constants
  use mod_dynparam
  use mod_runparams
  use mod_regcm_types
  use mod_memutil

  implicit none

  private

  type(domain) , public :: mddom
  type(surfstate) , public :: sfs
  type(atmosphere) , public :: mo_atm

  public :: allocate_mod_atm_interface
  public :: setup_model_indexes
  public :: fill_atmosphere

  contains

    subroutine setup_model_indexes
      implicit none
      jde1  = global_dot_jstart
      jdi1  = global_dot_jstart
      jdii1 = global_dot_jstart
      jde2  = global_dot_jend
      jdi2  = global_dot_jend
      jdii2 = global_dot_jend
      ide1  = global_dot_istart
      idi1  = global_dot_istart
      idii1 = global_dot_istart
      ide2  = global_dot_iend
      idi2  = global_dot_iend
      idii2 = global_dot_iend
      jce1  = global_cross_jstart
      jci1  = global_cross_jstart
      jcii1 = global_cross_jstart
      jce2  = global_cross_jend
      jci2  = global_cross_jend
      jcii2 = global_cross_jend
      ice1  = global_cross_istart
      ici1  = global_cross_istart
      icii1 = global_cross_istart
      ice2  = global_cross_iend
      ici2  = global_cross_iend
      icii2 = global_cross_iend
      idi1ga = idi1 - 1
      idi2ga = idi2 + 1
      jdi1ga = jdi1 - 1
      jdi2ga = jdi2 + 1
      ide1ga = ide1 - 1
      ide2ga = ide2 + 1
      jde1ga = jde1 - 1
      jde2ga = jde2 + 1
      ici1ga = ici1 - 1
      ici2ga = ici2 + 1
      jci1ga = jci1 - 1
      jci2ga = jci2 + 1
      ice1ga = ice1 - 1
      ice2ga = ice2 + 1
      jce1ga = jce1 - 1
      jce2ga = jce2 + 1

      idi1gb = idi1 - 2
      idi2gb = idi2 + 2
      jdi1gb = jdi1 - 2
      jdi2gb = jdi2 + 2
      ide1gb = ide1 - 2
      ide2gb = ide2 + 2
      jde1gb = jde1 - 2
      jde2gb = jde2 + 2
      ici1gb = ici1 - 2
      ici2gb = ici2 + 2
      jci1gb = jci1 - 2
      jci2gb = jci2 + 2
      ice1gb = ice1 - 2
      ice2gb = ice2 + 2
      jce1gb = jce1 - 2
      jce2gb = jce2 + 2

    end subroutine setup_model_indexes

    subroutine allocate_atmosphere(atm)
      implicit none
      type(atmosphere) , intent(inout) :: atm
      call getmem3d(atm%u,jde1gb,jde2gb,ice1ga,ice2ga,1,kz,'atmstate:u')
      call getmem3d(atm%v,jce1ga,jce2ga,ide1gb,ide2gb,1,kz,'atmstate:v')
      call getmem3d(atm%ux,jce1gb,jce2gb,ice1ga,ice2ga,1,kz,'atmstate:ux')
      call getmem3d(atm%vx,jce1ga,jce2ga,ice1gb,ice2gb,1,kz,'atmstate:vx')
      call getmem3d(atm%t,jce1ga,jce2ga,ice1ga,ice2ga,1,kz,'atmstate:t')
      call getmem3d(atm%tetav,jce1ga,jce2ga,ice1ga,ice2ga,1,kz,'atmstate:tetav')
      call getmem3d(atm%w,jce1,jce2,ice1,ice2,1,kzp1,'atmstate:w')
      call getmem3d(atm%pai,jce1ga,jce2ga,ice1ga,ice2ga,1,kz,'atmstate:pai')
      call getmem4d(atm%qx,jce1ga,jce2ga,ice1ga,ice2ga,1,kz,1,nqx,'atmstate:qx')
      call getmem3d(atm%zeta,jce1gb,jce2gb,ice1gb,ice2gb,1,kz,'atmstate:zeta')

      call getmem3d(atm%rho,jce1,jce2,ice1,ice2,1,kz,'atmstate:rho')
      call getmem3d(atm%pf,jce1,jce2,ice1,ice2,1,kzp1,'atmstate:pf')
      call getmem3d(atm%p,jce1,jce2,ice1,ice2,1,kz,'atmstate:p')
      call getmem3d(atm%tvirt,jce1,jce2,ice1,ice2,1,kz,'atmstate:tvirt')
      call getmem3d(atm%zetaf,jce1,jce2,ice1,ice2,1,kzp1,'atmstate:zetaf')
      call getmem3d(atm%dz,jce1,jce2,ice1,ice2,1,kz,'atmstate:dz')
      call getmem3d(atm%qs,jce1,jce2,ice1,ice2,1,kz,'atmstate:qs')
      call getmem3d(atm%fmz,jce1ga,jce2ga,ice1ga,ice2ga,1,kz,'atmstate:fmz')
      call getmem3d(atm%fmzf,jce1,jce2,ice1,ice2,1,kzp1,'atmstate:fmzf')
    end subroutine allocate_atmosphere

    subroutine allocate_domain(dom)
      implicit none
      type(domain) , intent(inout) :: dom
      call getmem2d(dom%ht,jde1gb,jde2gb,ide1gb,ide2gb,'storage:ht')
      call getmem2d(dom%lndcat,jde1,jde2,ide1,ide2,'storage:lndcat')
      call getmem2d(dom%lndtex,jde1,jde2,ide1,ide2,'storage:lndtex')
      call getmem2d(dom%xlat,jde1ga,jde2ga,ide1ga,ide2ga,'storage:xlat')
      call getmem2d(dom%xlon,jde1ga,jde2ga,ide1ga,ide2ga,'storage:xlon')
      call getmem2d(dom%dlat,jde1,jde2,ide1,ide2,'storage:dlat')
      call getmem2d(dom%dlon,jde1,jde2,ide1,ide2,'storage:dlon')
      call getmem2d(dom%mask,jde1,jde2,ide1,ide2,'storage:mask')
      call getmem2d(dom%area,jde1,jde2,ide1,ide2,'storage:area')
      call getmem2d(dom%msfx,jde1,jde2,ide1,ide2,'storage:msfx')
      call getmem2d(dom%msfu,jde1ga,jde2ga,ide1,ide2,'storage:msfu')
      call getmem2d(dom%msfv,jde1,jde2,ide1ga,ide2ga,'storage:msfv')
      call getmem2d(dom%hx,jde1ga,jde2ga,ice1,ice2,'storage:hx')
      call getmem2d(dom%hy,jce1,jce2,ide1ga,ide2ga,'storage:hy')
      call getmem2d(dom%ulat,jde1,jde2,ide1,ide2,'storage:ulat')
      call getmem2d(dom%ulon,jde1,jde2,ide1,ide2,'storage:ulon')
      call getmem2d(dom%vlat,jde1,jde2,ide1,ide2,'storage:vlat')
      call getmem2d(dom%vlon,jde1,jde2,ide1,ide2,'storage:vlon')
      call getmem2d(dom%coriou,jde1,jde2,ice1,ice2,'storage:fu')
      call getmem2d(dom%coriov,jce1,jce2,ide1,ide2,'storage:fv')
      call getmem2d(dom%coriol,jde1,jde2,ide1,ide2,'storage:f')
      call getmem2d(dom%ldmsk,jci1,jci2,ici1,ici2,'storage:ldmsk')
      call getmem2d(dom%xmsf,jdi1,jdi2,idi1,idi2,'storage:xmsf')
      call getmem2d(dom%dmsf,jdi1,jdi2,idi1,idi2,'storage:dmsf')
    end subroutine allocate_domain

    subroutine allocate_surfstate(sfs)
      implicit none
      type(surfstate) , intent(inout) :: sfs
      call getmem2d(sfs%ps,jce1ga,jce2ga,ice1ga,ice2ga,'surf:ps')
    end subroutine allocate_surfstate

    subroutine allocate_mod_atm_interface
      implicit none
      call allocate_domain(mddom)
      call allocate_atmosphere(mo_atm)
      call allocate_surfstate(sfs)
    end subroutine allocate_mod_atm_interface

    subroutine fill_atmosphere
      implicit none
      real(rkx) , parameter :: ts = 300.0_rkx
      real(rkx) , parameter :: qs = 0.01865_rkx
      real(rkx) , parameter :: ps = 101480.0_rkx
      real(rkx) , parameter :: tropoh = 15000.0_rkx
      real(rkx) :: tvs , tvt , zi , qi
      real(rkx) :: tv1 , tv2 , lrt , tv , zz , zb , p , zdelta
      integer(ik4) :: i , j , k
      integer(ik4) :: mx , my
      real(rkx) :: mr , md
      tvs = ts * (1.0_rkx + ep1*qs)
      tvt = tvs - lrate * tropoh
      sfs%ps = ps
      do concurrent ( j = jce1:jce2 , i = ice1:ice2 )
        do k = 1 , kz
          if ( mo_atm%zeta(j,i,k) < tropoh ) then
            zi = mo_atm%zeta(j,i,k)
            qi = qs*exp(-zi/4000.0_rkx) * exp((-zi/7500.0_rkx)**2)
            mo_atm%qx(j,i,k,iqv) = qi
            mo_atm%t(j,i,k) = (tvs - lrate * zi)/(1.0_rkx + qi)
          else
            mo_atm%qx(j,i,k,iqv) = 1.0e-8_rkx
            mo_atm%t(j,i,k) = tvt
          end if
          mo_atm%qx(j,i,k,iqfrst:iqlst) = 0.0_rkx
        end do
      end do
      mx = jx / 2
      my = iy / 2
      mr = sqrt((real(jx)/5.0_rkx)**2+(real(iy)/5_rkx)**2)/2.0_rkx
      do concurrent ( j = jce1:jce2 , i = ice1:ice2 )
        md = sqrt(real(j-mx)**2+real(i-my)**2)
        if ( md < mr ) then
          mo_atm%t(j,i,kz) = mo_atm%t(j,i,kz) + (mr-md)/mr
        end if
      end do
      mo_atm%u = 1.0_rkx
      mo_atm%v = 1.0_rkx
      ! Hydrostatic initialization of pai
      do i = ice1 , ice2
        do j = jce1 , jce2
          zdelta = mo_atm%zeta(j,i,kz)*egrav
          tv1 = mo_atm%t(j,i,kz) * (d_one + ep1*mo_atm%qx(j,i,kz,iqv))
          tv2 = mo_atm%t(j,i,kz-1) * (d_one + ep1*mo_atm%qx(j,i,kz-1,iqv))
          lrt = (tv2-tv1)/(mo_atm%zeta(j,i,kz-1)-mo_atm%zeta(j,i,kz))
          lrt = 0.65_rkx*lrt - 0.35_rkx*lrate
          tv = tv1 - 0.5_rkx*mo_atm%zeta(j,i,kz)*lrt
          zz = d_one/(rgas*tv)
          p = sfs%ps(j,i) * exp(-zdelta*zz)
          mo_atm%pai(j,i,kz) = (p/p00)**rovcp
        end do
      end do
      do k = kzm1 , 1 , -1
        do i = ice1 , ice2
          do j = jce1 , jce2
            tv1 = mo_atm%t(j,i,k) * (d_one + ep1*mo_atm%qx(j,i,k,iqv))
            tv2 = mo_atm%t(j,i,k+1) * (d_one + ep1*mo_atm%qx(j,i,k+1,iqv))
            zb = d_two*egrav*mo_dzita/(mo_atm%fmzf(j,i,k+1)*cpd) + tv1 - tv2
            zdelta = sqrt(zb**2 + d_four * tv2 * tv1)
            mo_atm%pai(j,i,k) = -mo_atm%pai(j,i,k+1)/(d_two*tv2)*(zb-zdelta)
          end do
        end do
      end do
    end subroutine fill_atmosphere

end module mod_atmosphere

! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
