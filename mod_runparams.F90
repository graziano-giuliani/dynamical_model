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

module mod_runparams

  use mod_intkinds
  use mod_realkinds
  use mod_constants
  use mod_dynparam
  use mod_memutil

  implicit none

  private

  integer(ik4) , public :: global_dot_istart
  integer(ik4) , public :: global_dot_iend
  integer(ik4) , public :: global_cross_istart
  integer(ik4) , public :: global_cross_iend
  integer(ik4) , public :: global_dot_jstart
  integer(ik4) , public :: global_dot_jend
  integer(ik4) , public :: global_cross_jstart
  integer(ik4) , public :: global_cross_jend

  integer(ik4) , public :: nqx , iqfrst , iqlst
  integer(ik4) , public , parameter :: iqv = 1
  integer(ik4) , public , parameter :: iqc = 2
  integer(ik4) , public , parameter :: iqi = 3
  integer(ik4) , public , parameter :: iqr = 4
  integer(ik4) , public , parameter :: iqs = 5
  integer(ik4) , public , parameter :: iqg = 6
  integer(ik4) , public , parameter :: iqh = 7
  integer(ik4) , public , parameter :: cqn = 8
  integer(ik4) , public , parameter :: cqc = 9
  integer(ik4) , public , parameter :: cqr = 10

  ! Moloch
  real(rkx) , public :: mo_dzita , mo_anu2
  logical , public :: mo_divfilter = .false.
  integer(ik4) , public :: mo_nzfilt
  integer(ik4) , public :: mo_nadv
  integer(ik4) , public :: mo_nsound

  real(rkx) , public :: dt , rdt
  real(rkx) , public :: dx , rdx

  real(rkx) , pointer , dimension(:) , public :: zita , zitah
  real(rkx) , pointer , dimension(:) , public :: sigma , hsigma
  real(rkx) , pointer , dimension(:) , public :: ffilt , ak , bk

  public :: allocate_mod_runparams

  contains

  subroutine allocate_mod_runparams
    implicit none
    call getmem1d(zita,1,kzp1,'mod_runparams:zita')
    call getmem1d(zitah,1,kz,'mod_runparams:zitah')
    call getmem1d(ffilt,1,kz,'mod_runparams:ffilt')
!$acc enter data create(ffilt)
    call getmem1d(ak,1,kz,'mod_runparams:ak')
    call getmem1d(bk,1,kz,'mod_runparams:bk')
    call getmem1d(hsigma,1,kz,'mod_runparams:hsigma')
    call getmem1d(sigma,1,kzp1,'mod_runparams:sigma')
  end subroutine allocate_mod_runparams

end module mod_runparams
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
