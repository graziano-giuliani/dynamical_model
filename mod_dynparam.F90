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
module mod_dynparam

  use, intrinsic :: iso_fortran_env
  use mod_intkinds
  use mod_realkinds
  use mpi

  implicit none

  private

  integer(ik4) , public  :: nstep = 2000
  integer(ik4) , public :: iprint = 20
  integer(ik4) , public :: iout = 100
  integer(ik4) , public :: iy = 100
  integer(ik4) , public :: jx = 100
  integer(ik4) , public :: kz = 30

  real(rkx) , public :: mo_ztop , mo_a0 , mo_h

  character(len=6) , public , parameter :: iproj = 'ROTLLR'

  real(rkx) , public :: ds

  character(len=512) :: aline
  character(len=8) :: cline

  integer(ik4) , public :: ide1 , ide2
  integer(ik4) , public :: jde1 , jde2
  integer(ik4) , public :: idi1 , idi2
  integer(ik4) , public :: jdi1 , jdi2
  integer(ik4) , public :: idii1 , idii2
  integer(ik4) , public :: jdii1 , jdii2

  integer(ik4) , public :: ice1 , ice2
  integer(ik4) , public :: jce1 , jce2
  integer(ik4) , public :: ici1 , ici2
  integer(ik4) , public :: jci1 , jci2
  integer(ik4) , public :: icii1 , icii2
  integer(ik4) , public :: jcii1 , jcii2

  integer(ik4) , public :: ici1ga , ici2ga , jci1ga , jci2ga
  integer(ik4) , public :: ice1ga , ice2ga , jce1ga , jce2ga
  integer(ik4) , public :: ide1ga , ide2ga , jde1ga , jde2ga
  integer(ik4) , public :: idi1ga , idi2ga , jdi1ga , jdi2ga
  integer(ik4) , public :: ici1gb , ici2gb , jci1gb , jci2gb
  integer(ik4) , public :: ice1gb , ice2gb , jce1gb , jce2gb
  integer(ik4) , public :: ide1gb , ide2gb , jde1gb , jde2gb
  integer(ik4) , public :: idi1gb , idi2gb , jdi1gb , jdi2gb

  integer(ik4) , public :: mycomm
  integer(ik4) , public :: nproc
  integer(ik4) , public :: myid
  integer(ik4) , public :: njxcpus , niycpus
  integer(ik4) , public :: iyp , jxp

  integer(ik4) , public :: idot1 , idot2
  integer(ik4) , public :: jdot1 , jdot2
  integer(ik4) , public :: icross1 , icross2
  integer(ik4) , public :: jcross1 , jcross2

  integer(ik4) , public :: kzp1 , kzp2 , kzm1 , kzm2
  integer(ik4) , public :: jxp1 , jxp2 , jxm1 , jxm2
  integer(ik4) , public :: iyp1 , iyp2 , iym1 , iym2

  integer(ik4) , public :: niout, njout
  integer(ik4) , public :: iout1, iout2
  integer(ik4) , public :: jout1, jout2

  integer(ik4) , public :: nicross , njcross
  integer(ik4) , public :: nidot , njdot

  public :: initparam , fatal

  contains

  subroutine initparam

    njxcpus = -1
    niycpus = -1

!   Setup all convenience dimensions

    iym1 = iy - 1
    iym2 = iy - 2
    jxm1 = jx - 1
    jxm2 = jx - 2
    kzm1 = kz - 1
    kzm2 = kz - 2
    kzp1 = kz + 1
    kzp2 = kz + 2
    jdot1 = 1
    jdot2 = jx
    jcross1 = 1
    jcross2 = jx
    jout1 = 1
    jout2 = jx
    idot1 = 1
    idot2 = iy
    icross1 = 1
    icross2 = iy
    iout1 = 1
    iout2 = iy
    njcross = jcross2-jcross1+1
    nicross = icross2-icross1+1
    njdot = jdot2-jdot1+1
    nidot = idot2-idot1+1
    njout = jout2-jout1+1
    niout = iout2-iout1+1
  end subroutine initparam

  subroutine fatal(filename,line,str)
    implicit none
    character(*) , intent(in) :: filename , str
    integer(ik4) , intent(in) :: line
    integer(ik4) :: ierr , myid
    character(len=8) :: date
    character(len=10) :: time
    character(len=5) :: zone
    write (cline,'(i8)') line
    write (error_unit,*) '-------------- FATAL CALLED ---------------'
    if ( line > 0 ) then
      write (aline,*) 'Fatal in file: '//filename//' at line: '//trim(cline)
      write (error_unit,*) trim(aline)
    end if
    write (error_unit,*) str
    write (error_unit,*) '-------------------------------------------'
    call date_and_time(date,time,zone)
    call mpi_comm_rank(mycomm, myid, ierr)
    write (error_unit,*) 'Abort called by computing node ', myid, 'at ', &
            date(1:4),'-',date(5:6),'-',date(7:8),' ', &
            time(1:2),':',time(3:4),':',time(5:10),' ',&
            zone
    write(error_unit,*) 'Execution terminated because of runtime error'
    call mpi_abort(mycomm,1,ierr)
  end subroutine fatal

end module mod_dynparam
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
