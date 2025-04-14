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
!
module mod_mppparam

  use mod_intkinds
  use mod_realkinds
  use mod_memutil
  use mod_dynparam
  use mod_runparams
  use mod_constants
  use netcdf
  use mod_regcm_types

  implicit none

  private

  public :: set_nproc

  type grid_nc_var2d
    character(len=64) :: varname
    integer(ik4) :: irec = -1
    integer(ik4) :: ncid = -1
    integer(ik4) :: varid = -1
    integer(ik4) :: nx = 0
    integer(ik4) :: ny = 0
    integer(ik4) :: mynx1 = 0
    integer(ik4) :: mynx2 = 0
    integer(ik4) :: myny1 = 0
    integer(ik4) :: myny2 = 0
    real(rk8) , pointer , dimension(:,:) :: val => null()
    real(rk8) , pointer , dimension(:,:) :: iobuf => null()
  end type grid_nc_var2d

  type grid_nc_var3d
    character(len=64) :: varname
    integer(ik4) :: irec = -1
    integer(ik4) :: ncid = -1
    integer(ik4) :: varid = -1
    integer(ik4) :: nx = 0
    integer(ik4) :: ny = 0
    integer(ik4) :: mynx1 = 0
    integer(ik4) :: mynx2 = 0
    integer(ik4) :: myny1 = 0
    integer(ik4) :: myny2 = 0
    integer(ik4) :: nz = 0
    real(rk8) , pointer , dimension(:,:,:) :: val => null()
    real(rk8) , pointer , dimension(:,:,:) :: iobuf => null()
  end type grid_nc_var3d

  type grid_nc_var4d
    character(len=64) :: varname
    integer(ik4) :: irec = -1
    integer(ik4) :: ncid = -1
    integer(ik4) :: varid = -1
    integer(ik4) :: nx = 0
    integer(ik4) :: ny = 0
    integer(ik4) :: mynx1 = 0
    integer(ik4) :: mynx2 = 0
    integer(ik4) :: myny1 = 0
    integer(ik4) :: myny2 = 0
    integer(ik4) :: nz = 0
    integer(ik4) :: nl = 0
    real(rk8) , pointer , dimension(:,:,:,:) :: val => null()
    real(rk8) , pointer , dimension(:,:,:,:) :: iobuf => null()
  end type grid_nc_var4d

  public :: grid_nc_var2d , grid_nc_var3d , grid_nc_var4d

  interface grid_nc_create
    module procedure grid_nc_create_var2d , &
                     grid_nc_create_var3d , &
                     grid_nc_create_var4d
  end interface grid_nc_create

  interface grid_nc_write
    module procedure grid_nc_write_var2d , &
                     grid_nc_write_var3d , &
                     grid_nc_write_var4d
  end interface grid_nc_write

  interface grid_nc_destroy
    module procedure grid_nc_destroy_var2d , &
                     grid_nc_destroy_var3d , &
                     grid_nc_destroy_var4d
  end interface grid_nc_destroy

  public :: grid_nc_create , grid_nc_write , grid_nc_destroy

  contains

  subroutine set_nproc
    implicit none
    global_dot_jstart = 1
    global_dot_istart = 1
    global_dot_jend = jx
    global_dot_iend = iy
    global_cross_jstart = 1
    global_cross_istart = 1
    global_cross_jend = jx
    global_cross_iend = iy
  end subroutine set_nproc

  subroutine grid_nc_create_var2d(varname,ldot,val,xvar)
    implicit none
    character(len=*) , intent(in) :: varname
    logical , intent(in) :: ldot
    real(rk8) , pointer , dimension(:,:) , intent(in) :: val
    type (grid_nc_var2d) , intent(inout) :: xvar
    integer(ik4) :: istat
    integer(ik4) , dimension(3) :: idims

    xvar%varname = varname
    call assignpnt(val,xvar%val)

    if ( ldot ) then
      xvar%nx = njdot
      xvar%ny = nidot
      xvar%mynx1 = jde1
      if ( lbound(val,1) > jde1 ) xvar%mynx1 = lbound(val,1)
      xvar%mynx2 = jde2
      if ( ubound(val,1) < jde2 ) xvar%mynx2 = ubound(val,1)
      xvar%myny1 = ide1
      if ( lbound(val,2) > ide1 ) xvar%myny1 = lbound(val,2)
      xvar%myny2 = ide2
      if ( ubound(val,2) < ide2 ) xvar%myny2 = ubound(val,2)
    else
      xvar%nx = njcross
      xvar%ny = nicross
      xvar%mynx1 = jce1
      if ( lbound(val,1) > jce1 ) xvar%mynx1 = lbound(val,1)
      xvar%mynx2 = jce2
      if ( ubound(val,1) < jce2 ) xvar%mynx2 = ubound(val,1)
      xvar%myny1 = ice1
      if ( lbound(val,2) > ice1 ) xvar%myny1 = lbound(val,2)
      xvar%myny2 = ice2
      if ( ubound(val,2) < ice2 ) xvar%myny2 = ubound(val,2)
    end if
    xvar%irec = 1
    istat = nf90_create(trim(xvar%varname)//'.nc',nf90_clobber,xvar%ncid)
    if ( istat /= nf90_noerr ) then
      write(error_unit, *) nf90_strerror(istat)
      return
    end if
    istat = nf90_def_dim(xvar%ncid, 'JX', xvar%nx, idims(1))
    if ( istat /= nf90_noerr ) then
      write(error_unit, *) nf90_strerror(istat)
      return
    end if
    istat = nf90_def_dim(xvar%ncid, 'IY', xvar%ny, idims(2))
    if ( istat /= nf90_noerr ) then
      write(error_unit, *) nf90_strerror(istat)
      return
    end if
    istat = nf90_def_dim(xvar%ncid, 'KTAU', nf90_unlimited, idims(3))
    if ( istat /= nf90_noerr ) then
      write(error_unit, *) nf90_strerror(istat)
      return
    end if
    istat = nf90_def_var(xvar%ncid,xvar%varname,nf90_double,idims,xvar%varid)
    if ( istat /= nf90_noerr ) then
      write(error_unit, *) nf90_strerror(istat)
      return
    end if
    istat = nf90_enddef(xvar%ncid)
    if ( istat /= nf90_noerr ) then
      write(error_unit, *) nf90_strerror(istat)
      return
    end if
    call getmem2d(xvar%iobuf,1,xvar%nx,1,xvar%ny,'var2d:iobuf')
  end subroutine grid_nc_create_var2d

  subroutine grid_nc_write_var2d(xvar)
    implicit none
    type (grid_nc_var2d) , intent(inout) :: xvar
    integer(ik4) :: istat
    integer(ik4) , dimension(3) :: istart , icount
    if ( .not. associated(xvar%val) .or. xvar%irec < 1 ) then
      return
    end if
    xvar%iobuf = xvar%val(1:xvar%nx,1:xvar%ny)
    istart(3) = xvar%irec
    istart(2) = 1
    istart(1) = 1
    icount(3) = 1
    icount(2) = xvar%ny
    icount(1) = xvar%nx
    istat = nf90_put_var(xvar%ncid,xvar%varid,xvar%iobuf,istart,icount)
    if ( istat /= nf90_noerr ) then
      write(error_unit, *) nf90_strerror(istat)
      return
    end if
    istat = nf90_sync(xvar%ncid)
    if ( istat /= nf90_noerr ) then
      write(error_unit, *) nf90_strerror(istat)
      return
    end if
    xvar%irec = xvar%irec + 1
  end subroutine grid_nc_write_var2d

  subroutine grid_nc_destroy_var2d(xvar)
    implicit none
    type (grid_nc_var2d) , intent(inout) :: xvar
    integer(ik4) :: istat
    istat = nf90_close(xvar%ncid)
    if ( istat /= nf90_noerr ) then
      write(error_unit, *) nf90_strerror(istat)
      return
    end if
    xvar%ncid = -1
    xvar%irec = -1
    nullify(xvar%val)
    call relmem2d(xvar%iobuf)
  end subroutine grid_nc_destroy_var2d

  subroutine grid_nc_create_var3d(varname,ldot,val,xvar)
    implicit none
    character(len=*) , intent(in) :: varname
    logical , intent(in) :: ldot
    real(rk8) , pointer , dimension(:,:,:) , intent(in) :: val
    type (grid_nc_var3d) , intent(inout) :: xvar
    integer(ik4) :: istat
    integer(ik4) , dimension(4) :: idims

    xvar%varname = varname
    call assignpnt(val,xvar%val)

    if ( ldot ) then
      xvar%nx = njdot
      xvar%ny = nidot
      xvar%mynx1 = jde1
      if ( lbound(val,1) > jde1 ) xvar%mynx1 = lbound(val,1)
      xvar%mynx2 = jde2
      if ( ubound(val,1) < jde2 ) xvar%mynx2 = ubound(val,1)
      xvar%myny1 = ide1
      if ( lbound(val,2) > ide1 ) xvar%myny1 = lbound(val,2)
      xvar%myny2 = ide2
      if ( ubound(val,2) < ide2 ) xvar%myny2 = ubound(val,2)
    else
      xvar%nx = njcross
      xvar%ny = nicross
      xvar%mynx1 = jce1
      if ( lbound(val,1) > jce1 ) xvar%mynx1 = lbound(val,1)
      xvar%mynx2 = jce2
      if ( ubound(val,1) < jce2 ) xvar%mynx2 = ubound(val,1)
      xvar%myny1 = ice1
      if ( lbound(val,2) > ice1 ) xvar%myny1 = lbound(val,2)
      xvar%myny2 = ice2
      if ( ubound(val,2) < ice2 ) xvar%myny2 = ubound(val,2)
    end if
    xvar%nz = size(xvar%val,3)
    xvar%irec = 1
    istat = nf90_create(trim(xvar%varname)//'.nc',nf90_clobber,xvar%ncid)
    if ( istat /= nf90_noerr ) then
      write(error_unit, *) nf90_strerror(istat)
      return
    end if
    istat = nf90_def_dim(xvar%ncid, 'JX', xvar%nx, idims(1))
    if ( istat /= nf90_noerr ) then
      write(error_unit, *) nf90_strerror(istat)
      return
    end if
    istat = nf90_def_dim(xvar%ncid, 'IY', xvar%ny, idims(2))
    if ( istat /= nf90_noerr ) then
      write(error_unit, *) nf90_strerror(istat)
      return
    end if
    istat = nf90_def_dim(xvar%ncid, 'KZ', xvar%nz, idims(3))
    if ( istat /= nf90_noerr ) then
      write(error_unit, *) nf90_strerror(istat)
      return
    end if
    istat = nf90_def_dim(xvar%ncid, 'KTAU', nf90_unlimited, idims(4))
    if ( istat /= nf90_noerr ) then
      write(error_unit, *) nf90_strerror(istat)
      return
    end if
    istat = nf90_def_var(xvar%ncid,xvar%varname,nf90_double,idims,xvar%varid)
    if ( istat /= nf90_noerr ) then
      write(error_unit, *) nf90_strerror(istat)
      return
    end if
    istat = nf90_enddef(xvar%ncid)
    if ( istat /= nf90_noerr ) then
      write(error_unit, *) nf90_strerror(istat)
      return
    end if
    call getmem3d(xvar%iobuf,1,xvar%nx,1,xvar%ny,1,xvar%nz,'var3d:iobuf')
  end subroutine grid_nc_create_var3d

  subroutine grid_nc_write_var3d(xvar)
    implicit none
    type (grid_nc_var3d) , intent(inout) :: xvar
    integer(ik4) :: istat
    integer(ik4) , dimension(4) :: istart , icount
    if ( .not. associated(xvar%val) .or. xvar%irec < 1 ) then
      return
    end if
    xvar%iobuf = xvar%val(1:xvar%nx,1:xvar%ny,1:xvar%nz)
    istart(4) = xvar%irec
    istart(3) = 1
    istart(2) = 1
    istart(1) = 1
    icount(4) = 1
    icount(3) = xvar%nz
    icount(2) = xvar%ny
    icount(1) = xvar%nx
    istat = nf90_put_var(xvar%ncid,xvar%varid,xvar%iobuf,istart,icount)
    if ( istat /= nf90_noerr ) then
      write(error_unit, *) nf90_strerror(istat)
      return
    end if
    istat = nf90_sync(xvar%ncid)
    if ( istat /= nf90_noerr ) then
      write(error_unit, *) nf90_strerror(istat)
      return
    end if
    xvar%irec = xvar%irec + 1
  end subroutine grid_nc_write_var3d

  subroutine grid_nc_destroy_var3d(xvar)
    implicit none
    type (grid_nc_var3d) , intent(inout) :: xvar
    integer(ik4) :: istat
    istat = nf90_close(xvar%ncid)
    if ( istat /= nf90_noerr ) then
      write(error_unit, *) nf90_strerror(istat)
      return
    end if
    xvar%ncid = -1
    xvar%irec = -1
    nullify(xvar%val)
    call relmem3d(xvar%iobuf)
  end subroutine grid_nc_destroy_var3d

  subroutine grid_nc_create_var4d(varname,ldot,val,xvar)
    implicit none
    character(len=*) , intent(in) :: varname
    logical , intent(in) :: ldot
    real(rk8) , pointer , dimension(:,:,:,:) , intent(in) :: val
    type (grid_nc_var4d) , intent(inout) :: xvar
    integer(ik4) :: istat
    integer(ik4) , dimension(5) :: idims

    xvar%varname = varname
    call assignpnt(val,xvar%val)

    if ( ldot ) then
      xvar%nx = njdot
      xvar%ny = nidot
      xvar%mynx1 = jde1
      if ( lbound(val,1) > jde1 ) xvar%mynx1 = lbound(val,1)
      xvar%mynx2 = jde2
      if ( ubound(val,1) < jde2 ) xvar%mynx2 = ubound(val,1)
      xvar%myny1 = ide1
      if ( lbound(val,2) > ide1 ) xvar%myny1 = lbound(val,2)
      xvar%myny2 = ide2
      if ( ubound(val,2) < ide2 ) xvar%myny2 = ubound(val,2)
    else
      xvar%nx = njcross
      xvar%ny = nicross
      xvar%mynx1 = jce1
      if ( lbound(val,1) > jce1 ) xvar%mynx1 = lbound(val,1)
      xvar%mynx2 = jce2
      if ( ubound(val,1) < jce2 ) xvar%mynx2 = ubound(val,1)
      xvar%myny1 = ice1
      if ( lbound(val,2) > ice1 ) xvar%myny1 = lbound(val,2)
      xvar%myny2 = ice2
      if ( ubound(val,2) < ice2 ) xvar%myny2 = ubound(val,2)
    end if
    xvar%nz = size(xvar%val,3)
    xvar%nl = size(xvar%val,4)
    xvar%irec = 1
    istat = nf90_create(trim(xvar%varname)//'.nc',nf90_clobber,xvar%ncid)
    if ( istat /= nf90_noerr ) then
      write(error_unit, *) nf90_strerror(istat)
      return
    end if
    istat = nf90_def_dim(xvar%ncid, 'JX', xvar%nx, idims(1))
    if ( istat /= nf90_noerr ) then
      write(error_unit, *) nf90_strerror(istat)
      return
    end if
    istat = nf90_def_dim(xvar%ncid, 'IY', xvar%ny, idims(2))
    if ( istat /= nf90_noerr ) then
      write(error_unit, *) nf90_strerror(istat)
      return
    end if
    istat = nf90_def_dim(xvar%ncid, 'KZ', xvar%nz, idims(3))
    if ( istat /= nf90_noerr ) then
      write(error_unit, *) nf90_strerror(istat)
      return
    end if
    istat = nf90_def_dim(xvar%ncid, 'LL', xvar%nl, idims(4))
    if ( istat /= nf90_noerr ) then
      write(error_unit, *) nf90_strerror(istat)
      return
    end if
    istat = nf90_def_dim(xvar%ncid, 'KTAU', nf90_unlimited, idims(5))
    if ( istat /= nf90_noerr ) then
      write(error_unit, *) nf90_strerror(istat)
      return
    end if
    istat = nf90_def_var(xvar%ncid,xvar%varname,nf90_double,idims,xvar%varid)
    if ( istat /= nf90_noerr ) then
      write(error_unit, *) nf90_strerror(istat)
      return
    end if
    istat = nf90_enddef(xvar%ncid)
    if ( istat /= nf90_noerr ) then
      write(error_unit, *) nf90_strerror(istat)
      return
    end if
    call getmem4d(xvar%iobuf,1,xvar%nx,1,xvar%ny,1,xvar%nz, &
                  1,xvar%nl,'var3d:iobuf')
  end subroutine grid_nc_create_var4d

  subroutine grid_nc_write_var4d(xvar)
    implicit none
    type (grid_nc_var4d) , intent(inout) :: xvar
    integer(ik4) :: istat
    integer(ik4) , dimension(5) :: istart , icount
    if ( .not. associated(xvar%val) .or. xvar%irec < 1 ) then
      return
    end if
    xvar%iobuf = xvar%val(1:xvar%nx,1:xvar%ny,1:xvar%nz,1:xvar%nl)
    istart(5) = xvar%irec
    istart(4) = 1
    istart(3) = 1
    istart(2) = 1
    istart(1) = 1
    icount(5) = 1
    icount(4) = xvar%nl
    icount(3) = xvar%nz
    icount(2) = xvar%ny
    icount(1) = xvar%nx
    istat = nf90_put_var(xvar%ncid,xvar%varid,xvar%iobuf,istart,icount)
    if ( istat /= nf90_noerr ) then
      write(error_unit, *) nf90_strerror(istat)
      return
    end if
    istat = nf90_sync(xvar%ncid)
    if ( istat /= nf90_noerr ) then
      write(error_unit, *) nf90_strerror(istat)
      return
    end if
    xvar%irec = xvar%irec + 1
  end subroutine grid_nc_write_var4d

  subroutine grid_nc_destroy_var4d(xvar)
    implicit none
    type (grid_nc_var4d) , intent(inout) :: xvar
    integer(ik4) :: istat
    istat = nf90_close(xvar%ncid)
    if ( istat /= nf90_noerr ) then
      write(error_unit, *) nf90_strerror(istat)
      return
    end if
    xvar%ncid = -1
    xvar%irec = -1
    nullify(xvar%val)
    call relmem4d(xvar%iobuf)
  end subroutine grid_nc_destroy_var4d

end module mod_mppparam
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
