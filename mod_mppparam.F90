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
  use mod_dynparam
  use mod_runparams
  use mod_constants
  use mod_memutil
  use netcdf
  use mod_regcm_types
  use mpi

  implicit none

  private

  logical , parameter :: lreorder = .false.

  integer(ik4) , public , parameter :: iocpu = 0 ! The id of the cpu doing I/O
  integer(ik4) , public , parameter :: italk = 0 ! Who is doing the print ?

  public :: set_nproc

  integer(ik4) :: cartesian_communicator
  integer(ik4) :: ccio , ccid

  integer(ik4) , public :: ncout_mpi_info = mpi_info_null

  integer(ik4) :: location(2)
  integer(ik4) :: right , left
  integer(ik4) :: bottom , top
  integer(ik4) :: topright , topleft
  integer(ik4) :: bottomright , bottomleft

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

  interface grid_distribute
    module procedure real8_2d_distribute ,   &
                     real8_3d_distribute ,   &
                     real8_4d_distribute ,   &
                     real4_2d_distribute ,   &
                     real4_3d_distribute ,   &
                     real4_4d_distribute ,   &
                     integer_2d_distribute , &
                     integer_3d_distribute , &
                     integer_4d_distribute , &
                     logical_2d_distribute , &
                     logical_3d_distribute , &
                     logical_4d_distribute
  end interface grid_distribute

  interface grid_collect
    module procedure real8_2d_collect ,    &
                     real8_2d_3d_collect , &
                     real8_3d_collect ,    &
                     real8_3d_2d_collect , &
                     real8_4d_collect ,    &
                     real8_4d_2d_collect , &
                     real4_2d_collect ,    &
                     real4_2d_3d_collect , &
                     real4_3d_collect ,    &
                     real4_3d_2d_collect , &
                     real4_4d_collect ,    &
                     real4_4d_2d_collect , &
                     integer_2d_collect ,  &
                     integer_3d_collect ,  &
                     integer_4d_collect ,  &
                     logical_2d_collect
  end interface grid_collect

  interface exchange
    module procedure real8_2d_exchange , &
                     real8_3d_exchange , &
                     real8_4d_exchange , &
                     real4_2d_exchange , &
                     real4_3d_exchange , &
                     real4_4d_exchange
  end interface exchange

  interface exchange_lrbt
    module procedure real8_2d_exchange_left_right_bottom_top , &
                     real8_3d_exchange_left_right_bottom_top , &
                     real8_4d_exchange_left_right_bottom_top , &
                     real4_2d_exchange_left_right_bottom_top , &
                     real4_3d_exchange_left_right_bottom_top , &
                     real4_4d_exchange_left_right_bottom_top
  end interface exchange_lrbt

  interface exchange_lr
    module procedure real8_2d_exchange_left_right , &
                     real8_3d_exchange_left_right , &
                     real8_4d_exchange_left_right , &
                     real4_2d_exchange_left_right , &
                     real4_3d_exchange_left_right , &
                     real4_4d_exchange_left_right
  end interface exchange_lr

  interface exchange_bt
    module procedure real8_2d_exchange_bottom_top , &
                     real8_3d_exchange_bottom_top , &
                     real8_4d_exchange_bottom_top , &
                     real4_2d_exchange_bottom_top , &
                     real4_3d_exchange_bottom_top , &
                     real4_4d_exchange_bottom_top
  end interface exchange_bt

  interface exchange_lb
    module procedure real8_2d_exchange_left_bottom , &
                     real8_3d_exchange_left_bottom , &
                     real8_4d_exchange_left_bottom , &
                     real4_2d_exchange_left_bottom , &
                     real4_3d_exchange_left_bottom , &
                     real4_4d_exchange_left_bottom
  end interface exchange_lb

  interface exchange_rt
    module procedure real8_2d_exchange_right_top , &
                     real8_3d_exchange_right_top , &
                     real8_4d_exchange_right_top , &
                     real4_2d_exchange_right_top , &
                     real4_3d_exchange_right_top , &
                     real4_4d_exchange_right_top
  end interface exchange_rt

  interface bcast
    module procedure bcast_logical,           &
                     bcast_int4,              &
                     bcast_int8,              &
                     bcast_real4,             &
                     bcast_real8,             &
                     bcast_arr_logical,       &
                     bcast_arr_character,     &
                     bcast_arr_text_list,     &
                     bcast_arr_int4,          &
                     bcast_arr_int8,          &
                     bcast_arr_real4,         &
                     bcast_arr_real8,         &
                     bcast_matr_real8,        &
                     bcast_matr_real4
  end interface bcast

  interface sumall
    module procedure sumall_real8 , &
                     sumall_real4 , &
                     sumall_int4 ,  &
                     sumall_int4_array
  end interface sumall

  interface maxall
    module procedure maxall_real8
    module procedure maxall_real4
    module procedure maxall_integer4
  end interface maxall

  interface minall
    module procedure minall_real8
    module procedure minall_real4
    module procedure minall_integer4
  end interface minall

  interface meanall
    module procedure meanall_real8
    module procedure meanall_real4
  end interface meanall

  interface cross2dot
    module procedure cross2dot2d
    module procedure cross2dot3d
  end interface cross2dot

  real(rk8) , pointer , dimension(:) :: r8vector1
  real(rk8) , pointer , dimension(:) :: r8vector2
  real(rk4) , pointer , dimension(:) :: r4vector1
  real(rk4) , pointer , dimension(:) :: r4vector2
  integer(ik4) , pointer , dimension(:) :: i4vector1
  integer(ik4) , pointer , dimension(:) :: i4vector2
  logical , pointer , dimension(:) :: lvector1
  logical , pointer , dimension(:) :: lvector2

  integer(ik4) , dimension(4) :: window
  integer(ik4) , pointer , dimension(:) :: windows
  integer(ik4) , dimension(1) , target :: ifake
  integer(ik4) , pointer , dimension(:) :: wincount
  integer(ik4) , pointer , dimension(:) :: windispl

  integer(ik4) :: mpierr

  integer(ik4) , parameter :: tag_bt = 1 ! FROM bottom TO top
  integer(ik4) , parameter :: tag_tb = 2 ! FROM top TO bottom
  integer(ik4) , parameter :: tag_lr = 3 ! FROM left TO right
  integer(ik4) , parameter :: tag_rl = 4 ! FROM right TO left
  integer(ik4) , parameter :: tag_brtl = 5 ! FROM bottomrigth TO topleft
  integer(ik4) , parameter :: tag_tlbr = 6 ! FROM topleft TO bottomright
  integer(ik4) , parameter :: tag_bltr = 7 ! FROM bottomleft TO topright
  integer(ik4) , parameter :: tag_trbl = 8 ! FROM topright TO bottomleft

  public :: exchange , exchange_lb , exchange_rt
  public :: exchange_lr , exchange_bt , exchange_lrbt
  public :: grid_distribute , grid_collect
  public :: uvcross2dot , uvdot2cross , cross2dot
  public :: bcast , sumall , maxall , minall , meanall
  public :: gather_r , gather_i
  public :: allgather_r , allgather_i
  public :: trueforall
  public :: allsync

  logical, save, public :: on_device = .false.

  contains

  subroutine bcast_logical(lval)
    implicit none
    logical , intent(inout) :: lval
    call mpi_bcast(lval,1,mpi_logical,iocpu,mycomm,mpierr)
  end subroutine bcast_logical

  subroutine bcast_int4(ival)
    implicit none
    integer(ik4) , intent(inout) :: ival
    call mpi_bcast(ival,1,mpi_integer4,iocpu,mycomm,mpierr)
  end subroutine bcast_int4

  subroutine bcast_int8(ival)
    implicit none
    integer(rk8) , intent(inout) :: ival
    call mpi_bcast(ival,1,mpi_integer8,iocpu,mycomm,mpierr)
  end subroutine bcast_int8

  subroutine bcast_real4(rval)
    implicit none
    real(rk4) , intent(inout) :: rval
    call mpi_bcast(rval,1,mpi_real4,iocpu,mycomm,mpierr)
  end subroutine bcast_real4

  subroutine bcast_real8(rval)
    implicit none
    real(rk8) , intent(inout) :: rval
    call mpi_bcast(rval,1,mpi_real8,iocpu,mycomm,mpierr)
  end subroutine bcast_real8

  subroutine bcast_arr_logical(lval)
    implicit none
    logical , dimension(:) , intent(inout) :: lval
    call mpi_bcast(lval,size(lval),mpi_logical,iocpu,mycomm,mpierr)
  end subroutine bcast_arr_logical

  subroutine bcast_arr_character(cval,is)
    implicit none
    character(len=*) , intent(inout) :: cval
    integer(ik4) , intent(in) :: is
    call mpi_bcast(cval,is,mpi_character,iocpu,mycomm,mpierr)
  end subroutine bcast_arr_character

  subroutine bcast_arr_text_list(cval,is)
    implicit none
    character(len=*) , intent(inout) , dimension(:) :: cval
    integer(ik4) , intent(in) :: is
    call mpi_bcast(cval,is*size(cval),mpi_character,iocpu,mycomm,mpierr)
  end subroutine bcast_arr_text_list

  subroutine bcast_arr_int4(ival)
    implicit none
    integer(ik4) , dimension(:) , intent(inout) :: ival
    call mpi_bcast(ival,size(ival),mpi_integer4,iocpu,mycomm,mpierr)
  end subroutine bcast_arr_int4

  subroutine bcast_arr_int8(ival)
    implicit none
    integer(rk8) , dimension(:) , intent(inout) :: ival
    call mpi_bcast(ival,size(ival),mpi_integer8,iocpu,mycomm,mpierr)
  end subroutine bcast_arr_int8

  subroutine bcast_arr_real4(rval)
    implicit none
    real(rk4) , dimension(:) , intent(inout) :: rval
    call mpi_bcast(rval,size(rval),mpi_real4,iocpu,mycomm,mpierr)
  end subroutine bcast_arr_real4

  subroutine bcast_arr_real8(rval)
    implicit none
    real(rk8) , dimension(:) , intent(inout) :: rval
    call mpi_bcast(rval,size(rval),mpi_real8,iocpu,mycomm,mpierr)
  end subroutine bcast_arr_real8

  subroutine bcast_matr_real8(rval)
    implicit none
    real(rk8) , dimension(:,:) , intent(inout) :: rval
    call mpi_bcast(rval,size(rval,1)*size(rval,2), &
                   mpi_real8,iocpu,mycomm,mpierr)
  end subroutine bcast_matr_real8

  subroutine bcast_matr_real4(rval)
    implicit none
    real(rk4) , dimension(:,:) , intent(inout) :: rval
    call mpi_bcast(rval,size(rval,1)*size(rval,2), &
                   mpi_real4,iocpu,mycomm,mpierr)
  end subroutine bcast_matr_real4

  subroutine trueforall(rlval,rtval)
    implicit none
    logical , intent(in) :: rlval
    logical , intent(out) :: rtval
    call mpi_allreduce(rlval,rtval,1,mpi_logical,mpi_lor,mycomm,mpierr)
  end subroutine trueforall

  subroutine sumall_real8(rlval,rtval)
    implicit none
    real(rk8) , intent(in) :: rlval
    real(rk8) , intent(out) :: rtval
    call mpi_allreduce(rlval,rtval,1,mpi_real8,mpi_sum,mycomm,mpierr)
  end subroutine sumall_real8

  subroutine sumall_real4(rlval,rtval)
    implicit none
    real(rk4) , intent(in) :: rlval
    real(rk4) , intent(out) :: rtval
    call mpi_allreduce(rlval,rtval,1,mpi_real4,mpi_sum,mycomm,mpierr)
  end subroutine sumall_real4

  subroutine maxall_real8(rlval,rtval)
    implicit none
    real(rk8) , intent(in) :: rlval
    real(rk8) , intent(out) :: rtval
    call mpi_allreduce(rlval,rtval,1,mpi_real8,mpi_max,mycomm,mpierr)
  end subroutine maxall_real8

  subroutine maxall_real4(rlval,rtval)
    implicit none
    real(rk4) , intent(in) :: rlval
    real(rk4) , intent(out) :: rtval
    call mpi_allreduce(rlval,rtval,1,mpi_real4,mpi_max,mycomm,mpierr)
  end subroutine maxall_real4

  subroutine maxall_integer4(rlval,rtval)
    implicit none
    integer(ik4) , intent(in) :: rlval
    integer(ik4) , intent(out) :: rtval
    call mpi_allreduce(rlval,rtval,1,mpi_integer4,mpi_max,mycomm,mpierr)
  end subroutine maxall_integer4

  subroutine meanall_real8(rlval,rtval)
    implicit none
    real(rk8) , intent(in) :: rlval
    real(rk8) , intent(out) :: rtval
    call mpi_allreduce(rlval,rtval,1,mpi_real8,mpi_sum,mycomm,mpierr)
    rtval = rtval/real(nproc,rk8)
  end subroutine meanall_real8

  subroutine meanall_real4(rlval,rtval)
    implicit none
    real(rk4) , intent(in) :: rlval
    real(rk4) , intent(out) :: rtval
    call mpi_allreduce(rlval,rtval,1,mpi_real4,mpi_sum,mycomm,mpierr)
    rtval = rtval/real(nproc,rk4)
  end subroutine meanall_real4

  subroutine minall_real8(rlval,rtval)
    implicit none
    real(rk8) , intent(in) :: rlval
    real(rk8) , intent(out) :: rtval
    call mpi_allreduce(rlval,rtval,1,mpi_real8,mpi_min,mycomm,mpierr)
  end subroutine minall_real8

  subroutine minall_real4(rlval,rtval)
    implicit none
    real(rk4) , intent(in) :: rlval
    real(rk4) , intent(out) :: rtval
    call mpi_allreduce(rlval,rtval,1,mpi_real4,mpi_min,mycomm,mpierr)
  end subroutine minall_real4

  subroutine minall_integer4(rlval,rtval)
    implicit none
    integer(ik4) , intent(in) :: rlval
    integer(ik4) , intent(out) :: rtval
    call mpi_allreduce(rlval,rtval,1,mpi_integer4,mpi_min,mycomm,mpierr)
  end subroutine minall_integer4

  subroutine sumall_int4(ilval,itval)
    implicit none
    integer(ik4) , intent(in) :: ilval
    integer(ik4) , intent(out) :: itval
    call mpi_allreduce(ilval,itval,1,mpi_integer4,mpi_sum,mycomm,mpierr)
  end subroutine sumall_int4

  subroutine sumall_int4_array(ilval,itval)
    implicit none
    integer(ik4) , dimension(:) , intent(in) :: ilval
    integer(ik4) , dimension(:) , intent(out) :: itval
    call mpi_allreduce(ilval,itval,size(itval),mpi_integer4, &
                       mpi_sum,mycomm,mpierr)
  end subroutine sumall_int4_array

  subroutine set_nproc
    implicit none
    integer(ik4) , dimension(2) :: cpus_per_dim
    logical , dimension(2) :: dim_period
    integer(ik4) , dimension(2) :: isearch
    integer(ik4) :: imaxcpus , imax1 , imax2 , imiss
    integer(ik4) :: maximum_buffer_size
    data dim_period /.false.,.false./

    if ( nproc == 1 ) then
      cpus_per_dim(1) = 1
      cpus_per_dim(2) = 1
      jxp =  jx
      iyp =  iy
      global_dot_jstart = 1
      global_dot_istart = 1
      global_dot_jend = jx
      global_dot_iend = iy
      global_cross_jstart = 1
      global_cross_istart = 1
      global_cross_jend = jx
      global_cross_iend = iy
      dim_period(1) = .true.
      dim_period(2) = .true.
      call mpi_cart_create(mycomm,2,cpus_per_dim,dim_period,lreorder, &
                           cartesian_communicator,mpierr)
      call mpi_comm_rank(cartesian_communicator,ccid,mpierr)
      call mpi_cart_coords(cartesian_communicator,ccid,2,location,mpierr)
      ccio = iocpu
    else
      dim_period(1) = .true.
      dim_period(2) = .true.
      if ( nproc < 4 ) then
        cpus_per_dim(1) = nproc
        cpus_per_dim(2) = 1
      else if ( nproc >= 4 ) then
        cpus_per_dim(1) = (nint(sqrt(dble(nproc)))/2)*2
        if ( iy > int(1.5*dble(jx)) ) then
          cpus_per_dim(1) = cpus_per_dim(1) - 1
          do while ( mod(nproc,cpus_per_dim(1)) /= 0 )
            cpus_per_dim(1) = cpus_per_dim(1) - 1
          end do
        else if ( jx > int(1.5*dble(iy)) ) then
          cpus_per_dim(1) = cpus_per_dim(1) + 1
          do while ( mod(nproc,cpus_per_dim(1)) /= 0 )
            cpus_per_dim(1) = cpus_per_dim(1) + 1
          end do
        else
          do while ( mod(nproc,cpus_per_dim(1)) /= 0 )
            cpus_per_dim(1) = cpus_per_dim(1) + 1
          end do
        end if
        cpus_per_dim(2) = nproc/cpus_per_dim(1)
        imaxcpus = cpus_per_dim(1)*cpus_per_dim(2)
        if ( mod(nproc,imaxcpus) /= 0 ) then
          write(error_unit,*) 'Work does not evenly divide.'
          write(error_unit,*) 'I have calculated : '
          write(error_unit,*) 'CPUS DIM1 = ', cpus_per_dim(1)
          write(error_unit,*) 'CPUS DIM2 = ', cpus_per_dim(2)
          imax1 = ((jx/3)/2)*2
          imax2 = ((iy/3)/2)*2
          write(error_unit,*) 'Suggested maximum number of CPUS jx: ', imax1
          write(error_unit,*) 'Suggested maximum number of CPUS iy: ', imax2
          write(error_unit,*) 'Closest number : ' , imaxcpus
          call fatal(__FILE__,__LINE__,'CPU/WORK mismatch')
        end if
      end if

      call mpi_cart_create(mycomm,2,cpus_per_dim,dim_period,lreorder, &
                           cartesian_communicator,mpierr)
      call mpi_comm_rank(cartesian_communicator,ccid,mpierr)
      call mpi_cart_coords(cartesian_communicator,ccid,2,location,mpierr)
      if ( myid == iocpu ) ccio = ccid
      call bcast(ccio)
      call mpi_cart_shift(cartesian_communicator, 0, 1, left, right, mpierr)
      call mpi_cart_shift(cartesian_communicator, 1, 1, bottom, top, mpierr)

      isearch(1) = location(1)+1
      isearch(2) = location(2)+1
      call mpi_cart_rank(cartesian_communicator,isearch,topright,mpierr)
      isearch(1) = location(1)-1
      isearch(2) = location(2)+1
      call mpi_cart_rank(cartesian_communicator,isearch,topleft,mpierr)
      isearch(1) = location(1)+1
      isearch(2) = location(2)-1
      call mpi_cart_rank(cartesian_communicator,isearch,bottomright,mpierr)
      isearch(1) = location(1)-1
      isearch(2) = location(2)-1
      call mpi_cart_rank(cartesian_communicator,isearch,bottomleft,mpierr)

      jxp =  jx/cpus_per_dim(1)
      iyp =  iy/cpus_per_dim(2)

      global_dot_jstart = location(1)*jxp+1
      global_dot_istart = location(2)*iyp+1

      if ( jxp * cpus_per_dim(1) < jx ) then
        imiss = jx - jxp * cpus_per_dim(1)
        if ( location(1) < imiss ) then
          global_dot_jstart = global_dot_jstart + location(1)
          jxp = jxp + 1
        else
          global_dot_jstart = global_dot_jstart + imiss
        end if
      end if
      if ( iyp * cpus_per_dim(2) < iy ) then
        imiss = iy - iyp * cpus_per_dim(2)
        if ( location(2) < imiss ) then
          global_dot_istart = global_dot_istart + location(2)
          iyp = iyp + 1
        else
          global_dot_istart = global_dot_istart + imiss
        end if
      end if

      global_dot_jend = global_dot_jstart+jxp-1
      global_dot_iend = global_dot_istart+iyp-1
      if ( global_dot_iend > iy .or. global_dot_jend > jx ) then
        write(error_unit,*) 'Cannot evenly divide!!!!'
        write(error_unit,*) 'Processor ',myid,' has I : ', global_dot_istart, &
                                                       global_dot_iend
        write(error_unit,*) 'Processor ',myid,' has J : ', global_dot_jstart, &
                                                       global_dot_jend
        call fatal(__FILE__,__LINE__,'DECOMPOSITION ERROR')
      end if

      global_cross_istart = global_dot_istart
      global_cross_iend = global_dot_iend

      global_cross_jstart = global_dot_jstart
      global_cross_jend = global_dot_jend
    end if
    !
    ! Check the results to be fit (minum for the advection is to have 3 points
    !
    if ( jxp < 3 .or. iyp < 3 ) then
      write(error_unit,*) 'CPUS DIM1 = ', cpus_per_dim(1)
      write(error_unit,*) 'CPUS DIM2 = ', cpus_per_dim(2)
      write(error_unit,*) 'Cannot have one processor with less than 3x3 points.'
      write(error_unit,*) 'Processor ',myid,' has ',jxp*iyp,' (',jxp,'x',iyp,')'
      call fatal(__FILE__,__LINE__,'Too much processors')
    end if

    if ( myid == italk ) then
      write(output_unit,*) 'CPUS DIM1 = ', cpus_per_dim(1)
      write(output_unit,*) 'CPUS DIM2 = ', cpus_per_dim(2)
      write(output_unit,*)
    end if

    if ( nproc > 1 ) then
      if ( myid == ccio ) then
        call getmem1d(windows,1,nproc*4,'set_nproc:windows')
      else
        windows => ifake
      end if
      call getmem1d(wincount,1,nproc*4,'set_nproc:wincount')
      call getmem1d(windispl,1,nproc*4,'set_nproc:windispl')
      ! Allocate to something should fit all
      maximum_buffer_size = jx*iy
      maximum_buffer_size = max(maximum_buffer_size,jxp*iyp*kzp1)
      call getmem1d(r8vector1,1,maximum_buffer_size,'set_nproc:r8vector1')
      call getmem1d(r8vector2,1,maximum_buffer_size,'set_nproc:r8vector2')
      call getmem1d(r4vector1,1,maximum_buffer_size,'set_nproc:r4vector1')
      call getmem1d(r4vector2,1,maximum_buffer_size,'set_nproc:r4vector2')
      call getmem1d(i4vector1,1,maximum_buffer_size,'set_nproc:i4vector1')
      call getmem1d(i4vector2,1,maximum_buffer_size,'set_nproc:i4vector2')
      call getmem1d(lvector1,1,maximum_buffer_size,'set_nproc:lvector1')
      call getmem1d(lvector2,1,maximum_buffer_size,'set_nproc:lvector2')
    end if
  end subroutine set_nproc

  integer(ik4) function glosplitw(j1,j2,i1,i2,ls) result(tsize)
    implicit none
    integer(ik4) , intent(in) :: j1 , j2 , i1 , i2
    logical , intent(in) , optional :: ls
    integer(ik4) :: isize , jsize , lsize , icpu , isub
    tsize = 0
    if ( nproc == 1 ) return
    isub = 1
    isize = i2-i1+1
    jsize = j2-j1+1
    tsize = isize*jsize*isub
    window(1) = i1
    window(2) = window(1)+isize-1
    window(3) = j1
    window(4) = window(3)+jsize-1
    call mpi_gather(window,4,mpi_integer4, &
                    windows,4,mpi_integer4,ccio,mycomm,mpierr)
    if ( ccid == ccio ) then
      do icpu = 0 , nproc-1
        window = windows(icpu*4+1:icpu*4+4)
        isize = window(2)-window(1)+1
        jsize = window(4)-window(3)+1
        lsize = isize*jsize*isub
        wincount(icpu+1) = lsize
        windispl(icpu+1) = sum(wincount(1:icpu))
      end do
    end if
  end function glosplitw

  subroutine real8_2d_do_distribute(mg,ml,j1,j2,i1,i2,tsize)
    implicit none
    real(rk8) , pointer , dimension(:,:) , intent(in) :: mg  ! model global
    real(rk8) , pointer , dimension(:,:) , intent(inout) :: ml ! model local
    integer(ik4) , intent(in) :: j1 , j2 , i1 , i2 , tsize
    integer(ik4) :: ib , i , j , icpu
    if ( nproc == 1 ) then
      do concurrent ( j = j1:j2, i = i1:i2 )
        ml(j,i) = mg(j,i)
      end do
      return
    end if
    if ( ccid == ccio ) then
      ib = 1
      do icpu = 0 , nproc-1
        window = windows(icpu*4+1:icpu*4+4)
        do i = window(1) , window(2)
          do j = window(3) , window(4)
            r8vector1(ib) = mg(j,i)
            ib = ib + 1
          end do
        end do
      end do
    end if
    call mpi_scatterv(r8vector1,wincount,windispl,mpi_real8, &
                      r8vector2,tsize,mpi_real8,ccio,mycomm,mpierr)
    ib = 1
    do i = i1 , i2
      do j = j1 , j2
        ml(j,i) = r8vector2(ib)
        ib = ib + 1
      end do
    end do
  end subroutine real8_2d_do_distribute

  subroutine real4_2d_do_distribute(mg,ml,j1,j2,i1,i2,tsize)
    implicit none
    real(rk4) , pointer , dimension(:,:) , intent(in) :: mg  ! model global
    real(rk4) , pointer , dimension(:,:) , intent(inout) :: ml ! model local
    integer(ik4) , intent(in) :: j1 , j2 , i1 , i2 , tsize
    integer(ik4) :: ib , i , j , icpu
    if ( nproc == 1 ) then
      do concurrent ( j = j1:j2, i = i1:i2 )
        ml(j,i) = mg(j,i)
      end do
      return
    end if
    if ( ccid == ccio ) then
      ib = 1
      do icpu = 0 , nproc-1
        window = windows(icpu*4+1:icpu*4+4)
        do i = window(1) , window(2)
          do j = window(3) , window(4)
            r4vector1(ib) = mg(j,i)
            ib = ib + 1
          end do
        end do
      end do
    end if
    call mpi_scatterv(r4vector1,wincount,windispl,mpi_real4, &
                      r4vector2,tsize,mpi_real4,ccio,mycomm,mpierr)
    ib = 1
    do i = i1 , i2
      do j = j1 , j2
        ml(j,i) = r4vector2(ib)
        ib = ib + 1
      end do
    end do
  end subroutine real4_2d_do_distribute

  subroutine integer4_2d_do_distribute(mg,ml,j1,j2,i1,i2,tsize)
    implicit none
    integer(ik4) , pointer , dimension(:,:) , intent(in) :: mg  ! model global
    integer(ik4) , pointer , dimension(:,:) , intent(inout) :: ml ! model local
    integer(ik4) , intent(in) :: j1 , j2 , i1 , i2 , tsize
    integer(ik4) :: ib , i , j , icpu
    if ( nproc == 1 ) then
      do concurrent ( j = j1:j2, i = i1:i2 )
        ml(j,i) = mg(j,i)
      end do
      return
    end if
    if ( ccid == ccio ) then
      ib = 1
      do icpu = 0 , nproc-1
        window = windows(icpu*4+1:icpu*4+4)
        do i = window(1) , window(2)
          do j = window(3) , window(4)
            i4vector1(ib) = mg(j,i)
            ib = ib + 1
          end do
        end do
      end do
    end if
    call mpi_scatterv(i4vector1,wincount,windispl,mpi_integer4, &
                      i4vector2,tsize,mpi_integer4,ccio,mycomm,mpierr)
    ib = 1
    do i = i1 , i2
      do j = j1 , j2
        ml(j,i) = i4vector2(ib)
        ib = ib + 1
      end do
    end do
  end subroutine integer4_2d_do_distribute

  subroutine logical_2d_do_distribute(mg,ml,j1,j2,i1,i2,tsize)
    implicit none
    logical , pointer , dimension(:,:) , intent(in) :: mg  ! model global
    logical , pointer , dimension(:,:) , intent(inout) :: ml ! model local
    integer(ik4) , intent(in) :: j1 , j2 , i1 , i2 , tsize
    integer(ik4) :: ib , i , j , icpu
    if ( nproc == 1 ) then
      do concurrent ( j = j1:j2, i = i1:i2 )
        ml(j,i) = mg(j,i)
      end do
      return
    end if
    if ( ccid == ccio ) then
      ib = 1
      do icpu = 0 , nproc-1
        window = windows(icpu*4+1:icpu*4+4)
        do i = window(1) , window(2)
          do j = window(3) , window(4)
            lvector1(ib) = mg(j,i)
            ib = ib + 1
          end do
        end do
      end do
    end if
    call mpi_scatterv(lvector1,wincount,windispl,mpi_logical, &
                      lvector2,tsize,mpi_logical,ccio,mycomm,mpierr)
    ib = 1
    do i = i1 , i2
      do j = j1 , j2
        ml(j,i) = lvector2(ib)
        ib = ib + 1
      end do
    end do
  end subroutine logical_2d_do_distribute

  subroutine real8_2d_distribute(mg,ml,j1,j2,i1,i2)
    implicit none
    real(rk8) , pointer , dimension(:,:) , intent(in) :: mg  ! model global
    real(rk8) , pointer , dimension(:,:) , intent(inout) :: ml ! model local
    integer(ik4) , intent(in) :: j1 , j2 , i1 , i2
    integer(ik4) :: tsize
    tsize = glosplitw(j1,j2,i1,i2)
    call real8_2d_do_distribute(mg,ml,j1,j2,i1,i2,tsize)
  end subroutine real8_2d_distribute

  subroutine real8_3d_distribute(mg,ml,j1,j2,i1,i2,k1,k2)
    implicit none
    real(rk8) , pointer , dimension(:,:,:) , intent(in) :: mg  ! model global
    real(rk8) , pointer , dimension(:,:,:) , intent(inout) :: ml ! model local
    integer(ik4) , intent(in) :: j1 , j2 , i1 , i2 , k1 , k2
    real(rk8) , pointer , dimension(:,:) :: mg2 => null()
    real(rk8) , pointer , dimension(:,:) :: ml2 => null()
    integer(ik4) :: tsize , k
    tsize = glosplitw(j1,j2,i1,i2)
    do k = k1 , k2
      if ( ccid == ccio ) call assignpnt(mg,mg2,k)
      call assignpnt(ml,ml2,k)
      call real8_2d_do_distribute(mg2,ml2,j1,j2,i1,i2,tsize)
    end do
  end subroutine real8_3d_distribute

  subroutine real8_4d_distribute(mg,ml,j1,j2,i1,i2,k1,k2,n1,n2)
    implicit none
    real(rk8) , pointer , dimension(:,:,:,:) , intent(in) :: mg  ! model global
    real(rk8) , pointer , dimension(:,:,:,:) , intent(inout) :: ml ! model local
    integer(ik4) , intent(in) :: j1 , j2 , i1 , i2 , k1 , k2 , n1 , n2
    real(rk8) , pointer , dimension(:,:) :: mg2 => null()
    real(rk8) , pointer , dimension(:,:) :: ml2 => null()
    integer(ik4) :: tsize , k , n
    tsize = glosplitw(j1,j2,i1,i2)
    do n = n1 , n2
      do k = k1 , k2
        if ( ccid == ccio ) call assignpnt(mg,mg2,k,n)
        call assignpnt(ml,ml2,k,n)
        call real8_2d_do_distribute(mg2,ml2,j1,j2,i1,i2,tsize)
      end do
    end do
  end subroutine real8_4d_distribute

  subroutine real4_2d_distribute(mg,ml,j1,j2,i1,i2)
    implicit none
    real(rk4) , pointer , dimension(:,:) , intent(in) :: mg  ! model global
    real(rk4) , pointer , dimension(:,:) , intent(inout) :: ml ! model local
    integer(ik4) , intent(in) :: j1 , j2 , i1 , i2
    integer(ik4) :: tsize
    tsize = glosplitw(j1,j2,i1,i2)
    call real4_2d_do_distribute(mg,ml,j1,j2,i1,i2,tsize)
  end subroutine real4_2d_distribute

  subroutine real4_3d_distribute(mg,ml,j1,j2,i1,i2,k1,k2)
    implicit none
    real(rk4) , pointer , dimension(:,:,:) , intent(in) :: mg  ! model global
    real(rk4) , pointer , dimension(:,:,:) , intent(inout) :: ml ! model local
    integer(ik4) , intent(in) :: j1 , j2 , i1 , i2 , k1 , k2
    real(rk4) , pointer , dimension(:,:) :: mg2 => null( )
    real(rk4) , pointer , dimension(:,:) :: ml2 => null( )
    integer(ik4) :: tsize , k
    tsize = glosplitw(j1,j2,i1,i2)
    do k = k1 , k2
      if ( ccid == ccio ) call assignpnt(mg,mg2,k)
      call assignpnt(ml,ml2,k)
      call real4_2d_do_distribute(mg2,ml2,j1,j2,i1,i2,tsize)
    end do
  end subroutine real4_3d_distribute

  subroutine real4_4d_distribute(mg,ml,j1,j2,i1,i2,k1,k2,n1,n2)
    implicit none
    real(rk4) , pointer , dimension(:,:,:,:) , intent(in) :: mg  ! model global
    real(rk4) , pointer , dimension(:,:,:,:) , intent(inout) :: ml ! model local
    integer(ik4) , intent(in) :: j1 , j2 , i1 , i2 , k1 , k2 , n1 , n2
    real(rk4) , pointer , dimension(:,:) :: mg2 => null( )
    real(rk4) , pointer , dimension(:,:) :: ml2 => null( )
    integer(ik4) :: tsize , k , n
    tsize = glosplitw(j1,j2,i1,i2)
    do n = n1 , n2
      do k = k1 , k2
        if ( ccid == ccio ) call assignpnt(mg,mg2,k,n)
        call assignpnt(ml,ml2,k,n)
        call real4_2d_do_distribute(mg2,ml2,j1,j2,i1,i2,tsize)
      end do
    end do
  end subroutine real4_4d_distribute

  subroutine integer_2d_distribute(mg,ml,j1,j2,i1,i2)
    implicit none
    integer(ik4) , pointer , dimension(:,:) , intent(in) :: mg  ! model global
    integer(ik4) , pointer , dimension(:,:) , intent(inout) :: ml ! model local
    integer(ik4) , intent(in) :: j1 , j2 , i1 , i2
    integer(ik4) :: tsize
    tsize = glosplitw(j1,j2,i1,i2)
    call integer4_2d_do_distribute(mg,ml,j1,j2,i1,i2,tsize)
  end subroutine integer_2d_distribute

  subroutine integer_3d_distribute(mg,ml,j1,j2,i1,i2,k1,k2)
    implicit none
    integer(ik4) , pointer , dimension(:,:,:) , intent(in) :: mg  ! model global
    integer(ik4) , pointer , dimension(:,:,:) , intent(inout) :: ml !model local
    integer(ik4) , intent(in) :: j1 , j2 , i1 , i2 , k1 , k2
    integer(ik4) , pointer , dimension(:,:) :: mg2 => null()
    integer(ik4) , pointer , dimension(:,:) :: ml2 => null()
    integer(ik4) :: tsize , k
    tsize = glosplitw(j1,j2,i1,i2)
    do k = k1 , k2
      if ( ccid == ccio ) call assignpnt(mg,mg2,k)
      call assignpnt(ml,ml2,k)
      call integer4_2d_do_distribute(mg2,ml2,j1,j2,i1,i2,tsize)
    end do
  end subroutine integer_3d_distribute

  subroutine integer_4d_distribute(mg,ml,j1,j2,i1,i2,k1,k2,n1,n2)
    implicit none
    integer(ik4) , pointer , dimension(:,:,:,:) , intent(in) :: mg  ! model glob
    integer(ik4) , pointer , dimension(:,:,:,:) , intent(inout) :: ml !model loc
    integer(ik4) , intent(in) :: j1 , j2 , i1 , i2 , k1 , k2 , n1 , n2
    integer(ik4) , pointer , dimension(:,:) :: mg2 => null( )
    integer(ik4) , pointer , dimension(:,:) :: ml2 => null( )
    integer(ik4) :: tsize , k , n
    tsize = glosplitw(j1,j2,i1,i2)
    do n = n1 , n2
      do k = k1 , k2
        if ( ccid == ccio ) call assignpnt(mg,mg2,k,n)
        call assignpnt(ml,ml2,k,n)
        call integer4_2d_do_distribute(mg2,ml2,j1,j2,i1,i2,tsize)
      end do
    end do
  end subroutine integer_4d_distribute

  subroutine logical_2d_distribute(mg,ml,j1,j2,i1,i2)
    implicit none
    logical , pointer , dimension(:,:) , intent(in) :: mg  ! model global
    logical , pointer , dimension(:,:) , intent(inout) :: ml ! model local
    integer(ik4) , intent(in) :: j1 , j2 , i1 , i2
    integer(ik4) :: tsize
    tsize = glosplitw(j1,j2,i1,i2)
    call logical_2d_do_distribute(mg,ml,j1,j2,i1,i2,tsize)
  end subroutine logical_2d_distribute

  subroutine logical_3d_distribute(mg,ml,j1,j2,i1,i2,k1,k2)
    implicit none
    logical , pointer , dimension(:,:,:) , intent(in) :: mg  ! model global
    logical , pointer , dimension(:,:,:) , intent(inout) :: ml !model local
    integer(ik4) , intent(in) :: j1 , j2 , i1 , i2 , k1 , k2
    logical , pointer , dimension(:,:) :: mg2 => null()
    logical , pointer , dimension(:,:) :: ml2 => null()
    integer(ik4) :: tsize , k
    tsize = glosplitw(j1,j2,i1,i2)
    do k = k1 , k2
      if ( ccid == ccio ) call assignpnt(mg,mg2,k)
      call assignpnt(ml,ml2,k)
      call logical_2d_do_distribute(mg2,ml2,j1,j2,i1,i2,tsize)
    end do
  end subroutine logical_3d_distribute

  subroutine logical_4d_distribute(mg,ml,j1,j2,i1,i2,k1,k2,n1,n2)
    implicit none
    logical , pointer , dimension(:,:,:,:) , intent(in) :: mg  ! model glob
    logical , pointer , dimension(:,:,:,:) , intent(inout) :: ml !model loc
    integer(ik4) , intent(in) :: j1 , j2 , i1 , i2 , k1 , k2 , n1 , n2
    logical , pointer , dimension(:,:) :: mg2 => null( )
    logical , pointer , dimension(:,:) :: ml2 => null( )
    integer(ik4) :: tsize , k , n
    tsize = glosplitw(j1,j2,i1,i2)
    do n = n1 , n2
      do k = k1 , k2
        if ( ccid == ccio ) call assignpnt(mg,mg2,k,n)
        call assignpnt(ml,ml2,k,n)
        call logical_2d_do_distribute(mg2,ml2,j1,j2,i1,i2,tsize)
      end do
    end do
  end subroutine logical_4d_distribute

  subroutine real8_2d_do_collect(ml,mg,j1,j2,i1,i2,tsize)
    implicit none
    real(rk8) , pointer , dimension(:,:) , intent(in) :: ml  ! model local
    real(rk8) , pointer , dimension(:,:) , intent(inout) :: mg ! model global
    integer(ik4) , intent(in) :: j1 , j2 , i1 , i2 , tsize
    integer(ik4) :: ib , i , j , icpu
    if ( nproc == 1 ) then
      do concurrent ( j = j1:j2, i = i1:i2 )
        mg(j,i) = ml(j,i)
      end do
      return
    end if
    ib = 1
    do i = i1 , i2
      do j = j1 , j2
        r8vector2(ib) = ml(j,i)
        ib = ib + 1
      end do
    end do
    call mpi_gatherv(r8vector2,tsize,mpi_real8, &
                     r8vector1,wincount,windispl,mpi_real8, &
                     ccio,mycomm,mpierr)
    if ( ccid == ccio ) then
      ib = 1
      do icpu = 0 , nproc-1
        window = windows(icpu*4+1:icpu*4+4)
        do i = window(1) , window(2)
          do j = window(3) , window(4)
            mg(j,i) = r8vector1(ib)
            ib = ib + 1
          end do
        end do
      end do
    end if
  end subroutine real8_2d_do_collect

  subroutine real4_2d_do_collect(ml,mg,j1,j2,i1,i2,tsize)
    implicit none
    real(rk4) , pointer , dimension(:,:) , intent(in) :: ml  ! model local
    real(rk4) , pointer , dimension(:,:) , intent(inout) :: mg ! model global
    integer(ik4) , intent(in) :: j1 , j2 , i1 , i2 , tsize
    integer(ik4) :: ib , i , j , icpu
    if ( nproc == 1 ) then
      do concurrent ( j = j1:j2, i = i1:i2 )
        mg(j,i) = ml(j,i)
      end do
      return
    end if
    ib = 1
    do i = i1 , i2
      do j = j1 , j2
        r4vector2(ib) = ml(j,i)
        ib = ib + 1
      end do
    end do
    call mpi_gatherv(r4vector2,tsize,mpi_real4, &
                     r4vector1,wincount,windispl,mpi_real4, &
                     ccio,mycomm,mpierr)
    if ( ccid == ccio ) then
      ib = 1
      do icpu = 0 , nproc-1
        window = windows(icpu*4+1:icpu*4+4)
        do i = window(1) , window(2)
          do j = window(3) , window(4)
            mg(j,i) = r4vector1(ib)
            ib = ib + 1
          end do
        end do
      end do
    end if
  end subroutine real4_2d_do_collect

  subroutine integer4_2d_do_collect(ml,mg,j1,j2,i1,i2,tsize)
    implicit none
    integer(ik4) , pointer , dimension(:,:) , intent(in) :: ml  ! model local
    integer(ik4) , pointer , dimension(:,:) , intent(inout) :: mg ! model global
    integer(ik4) , intent(in) :: j1 , j2 , i1 , i2 , tsize
    integer(ik4) :: ib , i , j , icpu
    if ( nproc == 1 ) then
      do concurrent ( j = j1:j2, i = i1:i2 )
        mg(j,i) = ml(j,i)
      end do
      return
    end if
    ib = 1
    do i = i1 , i2
      do j = j1 , j2
        i4vector2(ib) = ml(j,i)
        ib = ib + 1
      end do
    end do
    call mpi_gatherv(i4vector2,tsize,mpi_integer4, &
                     i4vector1,wincount,windispl,mpi_integer4, &
                     ccio,mycomm,mpierr)
    if ( ccid == ccio ) then
      ib = 1
      do icpu = 0 , nproc-1
        window = windows(icpu*4+1:icpu*4+4)
        do i = window(1) , window(2)
          do j = window(3) , window(4)
            mg(j,i) = i4vector1(ib)
            ib = ib + 1
          end do
        end do
      end do
    end if
  end subroutine integer4_2d_do_collect

  subroutine logical_2d_do_collect(ml,mg,j1,j2,i1,i2,tsize)
    implicit none
    logical , pointer , dimension(:,:) , intent(in) :: ml  ! model local
    logical , pointer , dimension(:,:) , intent(inout) :: mg ! model global
    integer(ik4) , intent(in) :: j1 , j2 , i1 , i2 , tsize
    integer(ik4) :: ib , i , j , icpu
    if ( nproc == 1 ) then
      do concurrent ( j = j1:j2, i = i1:i2 )
        mg(j,i) = ml(j,i)
      end do
      return
    end if
    ib = 1
    do i = i1 , i2
      do j = j1 , j2
        lvector2(ib) = ml(j,i)
        ib = ib + 1
      end do
    end do
    call mpi_gatherv(lvector2,tsize,mpi_logical, &
                     lvector1,wincount,windispl,mpi_logical, &
                     ccio,mycomm,mpierr)
    if ( ccid == ccio ) then
      ib = 1
      do icpu = 0 , nproc-1
        window = windows(icpu*4+1:icpu*4+4)
        do i = window(1) , window(2)
          do j = window(3) , window(4)
            mg(j,i) = lvector1(ib)
            ib = ib + 1
          end do
        end do
      end do
    end if
  end subroutine logical_2d_do_collect

  subroutine real8_2d_collect(ml,mg,j1,j2,i1,i2)
    implicit none
    real(rk8) , pointer , dimension(:,:) , intent(in) :: ml  ! model local
    real(rk8) , pointer , dimension(:,:) , intent(inout) :: mg ! model global
    integer(ik4) , intent(in) :: j1 , j2 , i1 , i2
    integer(ik4) :: tsize
    tsize = glosplitw(j1,j2,i1,i2)
    call real8_2d_do_collect(ml,mg,j1,j2,i1,i2,tsize)
  end subroutine real8_2d_collect

  subroutine real8_2d_3d_collect(ml,mg,j1,j2,i1,i2,k)
    implicit none
    real(rk8) , pointer , dimension(:,:) , intent(in) :: ml    ! model local
    real(rk8) , pointer , dimension(:,:,:) , intent(inout) :: mg ! model global
    integer(ik4) , intent(in) :: j1 , j2 , i1 , i2
    integer(ik4) , intent(in) , optional :: k
    real(rk8) , pointer , dimension(:,:) :: mg2 => null( )
    integer(ik4) :: kk
    integer(ik4) :: tsize
    tsize = glosplitw(j1,j2,i1,i2)
    kk = 1
    if ( present(k) ) kk = k
    if ( ccid == ccio )  call assignpnt(mg,mg2,kk)
    call real8_2d_do_collect(ml,mg2,j1,j2,i1,i2,tsize)
  end subroutine real8_2d_3d_collect

  subroutine real8_3d_collect(ml,mg,j1,j2,i1,i2,k1,k2)
    implicit none
    real(rk8) , pointer , dimension(:,:,:) , intent(in) :: ml  ! model local
    real(rk8) , pointer , dimension(:,:,:) , intent(inout) :: mg ! model global
    integer(ik4) , intent(in) :: j1 , j2 , i1 , i2 , k1 , k2
    real(rk8) , pointer , dimension(:,:) :: ml2 => null( )
    real(rk8) , pointer , dimension(:,:) :: mg2 => null( )
    integer(ik4) :: tsize , k
    tsize = glosplitw(j1,j2,i1,i2)
    do k = k1 , k2
      if ( ccid == ccio )  call assignpnt(mg,mg2,k)
      call assignpnt(ml,ml2,k)
      call real8_2d_do_collect(ml2,mg2,j1,j2,i1,i2,tsize)
    end do
  end subroutine real8_3d_collect

  subroutine real8_3d_2d_collect(ml,mg,j1,j2,i1,i2,k)
    implicit none
    real(rk8) , pointer , dimension(:,:,:) , intent(in) :: ml  ! model local
    real(rk8) , pointer , dimension(:,:) , intent(inout) :: mg   ! model global
    integer(ik4) , intent(in) :: j1 , j2 , i1 , i2 , k
    real(rk8) , pointer , dimension(:,:) :: ml2 => null( )
    integer(ik4) :: tsize
    tsize = glosplitw(j1,j2,i1,i2)
    call assignpnt(ml,ml2,k)
    call real8_2d_do_collect(ml2,mg,j1,j2,i1,i2,tsize)
  end subroutine real8_3d_2d_collect

  subroutine real8_4d_collect(ml,mg,j1,j2,i1,i2,k1,k2,n1,n2)
    implicit none
    real(rk8) , pointer , dimension(:,:,:,:) , intent(in) :: ml  ! model local
    real(rk8) , pointer , dimension(:,:,:,:) , intent(inout) :: mg ! model glob
    integer(ik4) , intent(in) :: j1 , j2 , i1 , i2 , k1 , k2 , n1 , n2
    real(rk8) , pointer , dimension(:,:) :: ml2 => null( )
    real(rk8) , pointer , dimension(:,:) :: mg2 => null( )
    integer(ik4) :: tsize , k , n
    tsize = glosplitw(j1,j2,i1,i2)
    do n = n1 , n2
      do k = k1 , k2
        if ( ccid == ccio ) call assignpnt(mg,mg2,k,n)
        call assignpnt(ml,ml2,k,n)
        call real8_2d_do_collect(ml2,mg2,j1,j2,i1,i2,tsize)
      end do
    end do
  end subroutine real8_4d_collect

  subroutine real8_4d_2d_collect(ml,mg,j1,j2,i1,i2,k,n)
    implicit none
    real(rk8) , pointer , dimension(:,:,:,:) , intent(in) :: ml ! model local
    real(rk8) , pointer , dimension(:,:) , intent(inout) :: mg ! model global
    integer(ik4) , intent(in) :: j1 , j2 , i1 , i2 , k , n
    real(rk8) , pointer , dimension(:,:) :: ml2 => null( )
    integer(ik4) :: tsize
    tsize = glosplitw(j1,j2,i1,i2)
    call assignpnt(ml,ml2,k,n)
    call real8_2d_do_collect(ml2,mg,j1,j2,i1,i2,tsize)
  end subroutine real8_4d_2d_collect

  subroutine real4_2d_collect(ml,mg,j1,j2,i1,i2)
    implicit none
    real(rk4) , pointer , dimension(:,:) , intent(in) :: ml  ! model local
    real(rk4) , pointer , dimension(:,:) , intent(inout) :: mg ! model global
    integer(ik4) , intent(in) :: j1 , j2 , i1 , i2
    integer(ik4) :: tsize
    tsize = glosplitw(j1,j2,i1,i2)
    call real4_2d_do_collect(ml,mg,j1,j2,i1,i2,tsize)
  end subroutine real4_2d_collect

  subroutine real4_2d_3d_collect(ml,mg,j1,j2,i1,i2,k)
    implicit none
    real(rk4) , pointer , dimension(:,:) , intent(in) :: ml    ! model local
    real(rk4) , pointer , dimension(:,:,:) , intent(inout) :: mg ! model global
    integer(ik4) , intent(in) :: j1 , j2 , i1 , i2
    integer(ik4) , intent(in) , optional :: k
    real(rk4) , pointer , dimension(:,:) :: mg2 => null( )
    integer(ik4) :: kk
    integer(ik4) :: tsize
    tsize = glosplitw(j1,j2,i1,i2)
    kk = 1
    if ( present(k) ) kk = k
    if ( ccid == ccio )  call assignpnt(mg,mg2,kk)
    call real4_2d_do_collect(ml,mg2,j1,j2,i1,i2,tsize)
  end subroutine real4_2d_3d_collect

  subroutine real4_3d_collect(ml,mg,j1,j2,i1,i2,k1,k2)
    implicit none
    real(rk4) , pointer , dimension(:,:,:) , intent(in) :: ml  ! model local
    real(rk4) , pointer , dimension(:,:,:) , intent(inout) :: mg ! model global
    integer(ik4) , intent(in) :: j1 , j2 , i1 , i2 , k1 , k2
    real(rk4) , pointer , dimension(:,:) :: ml2 => null( )
    real(rk4) , pointer , dimension(:,:) :: mg2 => null( )
    integer(ik4) :: tsize , k
    tsize = glosplitw(j1,j2,i1,i2)
    do k = k1 , k2
      if ( ccid == ccio )  call assignpnt(mg,mg2,k)
      call assignpnt(ml,ml2,k)
      call real4_2d_do_collect(ml2,mg2,j1,j2,i1,i2,tsize)
    end do
  end subroutine real4_3d_collect

  subroutine real4_3d_2d_collect(ml,mg,j1,j2,i1,i2,k)
    implicit none
    real(rk4) , pointer , dimension(:,:,:) , intent(in) :: ml  ! model local
    real(rk4) , pointer , dimension(:,:) , intent(inout) :: mg   ! model global
    integer(ik4) , intent(in) :: j1 , j2 , i1 , i2 , k
    real(rk4) , pointer , dimension(:,:) :: ml2 => null( )
    integer(ik4) :: tsize
    tsize = glosplitw(j1,j2,i1,i2)
    call assignpnt(ml,ml2,k)
    call real4_2d_do_collect(ml2,mg,j1,j2,i1,i2,tsize)
  end subroutine real4_3d_2d_collect

  subroutine real4_4d_collect(ml,mg,j1,j2,i1,i2,k1,k2,n1,n2)
    implicit none
    real(rk4) , pointer , dimension(:,:,:,:) , intent(in) :: ml  ! model local
    real(rk4) , pointer , dimension(:,:,:,:) , intent(inout) :: mg ! model glob
    integer(ik4) , intent(in) :: j1 , j2 , i1 , i2 , k1 , k2 , n1 , n2
    real(rk4) , pointer , dimension(:,:) :: ml2 => null( )
    real(rk4) , pointer , dimension(:,:) :: mg2 => null( )
    integer(ik4) :: tsize , k , n
    tsize = glosplitw(j1,j2,i1,i2)
    do n = n1 , n2
      do k = k1 , k2
        if ( ccid == ccio ) call assignpnt(mg,mg2,k,n)
        call assignpnt(ml,ml2,k,n)
        call real4_2d_do_collect(ml2,mg2,j1,j2,i1,i2,tsize)
      end do
    end do
  end subroutine real4_4d_collect

  subroutine real4_4d_2d_collect(ml,mg,j1,j2,i1,i2,k,n)
    implicit none
    real(rk4) , pointer , dimension(:,:,:,:) , intent(in) :: ml ! model local
    real(rk4) , pointer , dimension(:,:) , intent(inout) :: mg ! model global
    integer(ik4) , intent(in) :: j1 , j2 , i1 , i2 , k , n
    real(rk4) , pointer , dimension(:,:) :: ml2 => null( )
    integer(ik4) :: tsize
    tsize = glosplitw(j1,j2,i1,i2)
    call assignpnt(ml,ml2,k,n)
    call real4_2d_do_collect(ml2,mg,j1,j2,i1,i2,tsize)
  end subroutine real4_4d_2d_collect

  subroutine logical_2d_collect(ml,mg,j1,j2,i1,i2)
    implicit none
    logical , pointer , dimension(:,:) , intent(in) :: ml  ! model local
    logical , pointer , dimension(:,:) , intent(inout) :: mg ! model global
    integer(ik4) , intent(in) :: j1 , j2 , i1 , i2
    integer(ik4) :: tsize
    tsize = glosplitw(j1,j2,i1,i2)
    call logical_2d_do_collect(ml,mg,j1,j2,i1,i2,tsize)
  end subroutine logical_2d_collect

  subroutine integer_2d_collect(ml,mg,j1,j2,i1,i2)
    implicit none
    integer(ik4) , pointer , dimension(:,:) , intent(in) :: ml  ! model local
    integer(ik4) , pointer , dimension(:,:) , intent(inout) :: mg ! model global
    integer(ik4) , intent(in) :: j1 , j2 , i1 , i2
    integer(ik4) :: tsize
    tsize = glosplitw(j1,j2,i1,i2)
    call integer4_2d_do_collect(ml,mg,j1,j2,i1,i2,tsize)
  end subroutine integer_2d_collect

  subroutine integer_3d_collect(ml,mg,j1,j2,i1,i2,k1,k2)
    implicit none
    integer(ik4) , pointer , dimension(:,:,:) , intent(in) :: ml  ! model local
    integer(ik4) , pointer , dimension(:,:,:) , intent(inout) :: mg ! model glbl
    integer(ik4) , intent(in) :: j1 , j2 , i1 , i2 , k1 , k2
    integer(ik4) , pointer , dimension(:,:) :: ml2 => null()
    integer(ik4) , pointer , dimension(:,:) :: mg2 => null()
    integer(ik4) :: tsize , k
    tsize = glosplitw(j1,j2,i1,i2)
    do k = k1 , k2
      if ( ccid == ccio )  call assignpnt(mg,mg2,k)
      call assignpnt(ml,ml2,k)
      call integer4_2d_do_collect(ml2,mg2,j1,j2,i1,i2,tsize)
    end do
  end subroutine integer_3d_collect

  subroutine integer_4d_collect(ml,mg,j1,j2,i1,i2,k1,k2,n1,n2)
    implicit none
    integer(ik4) , pointer , dimension(:,:,:,:) , intent(in) :: ml  ! mdl local
    integer(ik4) , pointer , dimension(:,:,:,:) , intent(inout) :: mg ! mdl glob
    integer(ik4) , intent(in) :: j1 , j2 , i1 , i2 , k1 , k2 , n1 , n2
    integer(ik4) , pointer , dimension(:,:) :: ml2 => null( )
    integer(ik4) , pointer , dimension(:,:) :: mg2 => null( )
    integer(ik4) :: tsize , k , n
    tsize = glosplitw(j1,j2,i1,i2)
    do n = n1 , n2
      do k = k1 , k2
        if ( ccid == ccio ) call assignpnt(mg,mg2,k,n)
        call assignpnt(ml,ml2,k,n)
        call integer4_2d_do_collect(ml2,mg2,j1,j2,i1,i2,tsize)
      end do
    end do
  end subroutine integer_4d_collect

  subroutine real8_2d_exchange(ml,nex,j1,j2,i1,i2)
    implicit none
    real(rk8) , pointer , dimension(:,:) , intent(inout) :: ml
    integer(ik4) , intent(in) :: nex , j1 , j2  , i1 , i2
    integer(ik4) :: nx , ny
    integer(ik4) :: lb , rb
    integer(ik4) :: ndx , ndy , tx , ty , sizex , sizey

    nx = j2-j1+1
    ny = i2-i1+1
    lb = nex
    rb = nex
    if ( left  == mpi_proc_null ) lb = 0
    if ( right == mpi_proc_null ) rb = 0
    tx = ny
    ty = nx+lb+rb
    sizex = nex*tx
    sizey = nex*ty
    ndx = 2*sizex
    ndy = 2*sizey

    exchange : block

    integer(ik4), dimension(4) :: counts , displs
    real(rk8), dimension(ndx) :: sdatax
    real(rk8), dimension(ndy) :: sdatay
    real(rk8), dimension(ndx), volatile :: rdatax
    real(rk8), dimension(ndy), volatile :: rdatay
    integer(ik4) :: ib1 , ib2 , iex

    ib2 = 0
    do iex = 1 , nex
      ib1 = ib2 + 1
      ib2 = ib1 + tx - 1
      sdatax(ib1:ib2) = ml(j1+(iex-1),i1:i2)
    end do
    do iex = 1 , nex
      ib1 = ib2 + 1
      ib2 = ib1 + tx - 1
      sdatax(ib1:ib2) = ml(j2-(iex-1),i1:i2)
    end do

    counts = [ sizex, sizex, 0, 0 ]
    displs = [ 0, sizex, 2*sizex, 2*sizex ]
    call mpi_neighbor_alltoallv(sdatax, counts, displs, mpi_real8, &
        rdatax, counts, displs, mpi_real8, cartesian_communicator, mpierr)

    ib2 = 0
    if ( left /= mpi_proc_null ) then
      do iex = 1 , nex
        ib1 = ib2 + 1
        ib2 = ib1 + tx - 1
        ml(j1-iex,i1:i2) = rdatax(ib1:ib2)
      end do
    else
      ib2 = ib2 + sizex
    end if
    if ( right /= mpi_proc_null ) then
      do iex = 1 , nex
        ib1 = ib2 + 1
        ib2 = ib1 + tx - 1
        ml(j2+iex,i1:i2) = rdatax(ib1:ib2)
      end do
    else
      ib2 = ib2 + sizex
    end if

    ib2 = 0
    do iex = 1 , nex
      ib1 = ib2 + 1
      ib2 = ib1 + ty - 1
      sdatay(ib1:ib2) = ml(j1-lb:j2+rb,i1+(iex-1))
    end do
    do iex = 1 , nex
      ib1 = ib2 + 1
      ib2 = ib1 + ty - 1
      sdatay(ib1:ib2) = ml(j1-lb:j2+rb,i2-(iex-1))
    end do
    counts = [ 0, 0, sizey, sizey ]
    displs = [ 0, 0, 0, sizey ]
    call mpi_neighbor_alltoallv(sdatay, counts, displs, mpi_real8, &
        rdatay, counts, displs, mpi_real8, cartesian_communicator, mpierr)

    ib2 = 0
    if ( bottom /= mpi_proc_null ) then
      do iex = 1 , nex
        ib1 = ib2 + 1
        ib2 = ib1 + ty - 1
        ml(j1-lb:j2+rb,i1-iex) = rdatay(ib1:ib2)
      end do
    else
      ib2 = ib2 + sizey
    end if
    if ( top /= mpi_proc_null ) then
      do iex = 1 , nex
        ib1 = ib2 + 1
        ib2 = ib1 + ty - 1
        ml(j1-lb:j2+rb,i2+iex) = rdatay(ib1:ib2)
      end do
    else
      ib2 = ib2 + sizey
    end if

    end block exchange
  end subroutine real8_2d_exchange

  subroutine real8_3d_exchange(ml,nex,j1,j2,i1,i2,k1,k2)
    implicit none
    real(rk8) , pointer , dimension(:,:,:) , intent(inout) :: ml
    integer(ik4) , intent(in) :: nex , j1 , j2  , i1 , i2 , k1 , k2
    integer(ik4) :: nx , ny , nk
    integer(ik4) :: lb , rb
    integer(ik4) :: ndx , ndy , tx , ty , sizex , sizey

    nx = j2-j1+1
    ny = i2-i1+1
    nk = k2-k1+1
    lb = nex
    rb = nex
    if ( left  == mpi_proc_null ) lb = 0
    if ( right == mpi_proc_null ) rb = 0
    tx = ny
    ty = nx+lb+rb
    sizex = nex*tx*nk
    sizey = nex*ty*nk
    ndx = 2*sizex
    ndy = 2*sizey

    exchange : block

    integer(ik4), dimension(4) :: counts , displs
    real(rk8), dimension(ndx) :: sdatax
    real(rk8), dimension(ndy) :: sdatay
    real(rk8), dimension(ndx), volatile :: rdatax
    real(rk8), dimension(ndy), volatile :: rdatay
    integer(ik4) :: ib1 , ib2 , iex , k

    ib2 = 0
    do k = k1 , k2
      do iex = 1 , nex
        ib1 = ib2 + 1
        ib2 = ib1 + tx - 1
        sdatax(ib1:ib2) = ml(j1+(iex-1),i1:i2,k)
      end do
    end do
    do k = k1 , k2
      do iex = 1 , nex
        ib1 = ib2 + 1
        ib2 = ib1 + tx - 1
        sdatax(ib1:ib2) = ml(j2-(iex-1),i1:i2,k)
      end do
    end do

    counts = [ sizex, sizex, 0, 0 ]
    displs = [ 0, sizex, 2*sizex, 2*sizex ]
    call mpi_neighbor_alltoallv(sdatax, counts, displs, mpi_real8, &
        rdatax, counts, displs, mpi_real8, cartesian_communicator, mpierr)

    ib2 = 0
    if ( left /= mpi_proc_null ) then
      do k = k1 , k2
        do iex = 1 , nex
          ib1 = ib2 + 1
          ib2 = ib1 + tx - 1
          ml(j1-iex,i1:i2,k) = rdatax(ib1:ib2)
        end do
      end do
    else
      ib2 = ib2 + sizex
    end if
    if ( right /= mpi_proc_null ) then
      do k = k1 , k2
        do iex = 1 , nex
          ib1 = ib2 + 1
          ib2 = ib1 + tx - 1
          ml(j2+iex,i1:i2,k) = rdatax(ib1:ib2)
        end do
      end do
    else
      ib2 = ib2 + sizex
    end if

    ib2 = 0
    do k = k1 , k2
      do iex = 1 , nex
        ib1 = ib2 + 1
        ib2 = ib1 + ty - 1
        sdatay(ib1:ib2) = ml(j1-lb:j2+rb,i1+(iex-1),k)
      end do
    end do
    do k = k1 , k2
      do iex = 1 , nex
        ib1 = ib2 + 1
        ib2 = ib1 + ty - 1
        sdatay(ib1:ib2) = ml(j1-lb:j2+rb,i2-(iex-1),k)
      end do
    end do
    counts = [ 0, 0, sizey, sizey ]
    displs = [ 0, 0, 0, sizey ]
    call mpi_neighbor_alltoallv(sdatay, counts, displs, mpi_real8, &
        rdatay, counts, displs, mpi_real8, cartesian_communicator, mpierr)

    ib2 = 0
    if ( bottom /= mpi_proc_null ) then
      do k = k1 , k2
        do iex = 1 , nex
          ib1 = ib2 + 1
          ib2 = ib1 + ty - 1
          ml(j1-lb:j2+rb,i1-iex,k) = rdatay(ib1:ib2)
        end do
      end do
    else
      ib2 = ib2 + sizey
    end if
    if ( top /= mpi_proc_null ) then
      do k = k1 , k2
        do iex = 1 , nex
          ib1 = ib2 + 1
          ib2 = ib1 + ty - 1
          ml(j1-lb:j2+rb,i2+iex,k) = rdatay(ib1:ib2)
        end do
      end do
    else
      ib2 = ib2 + sizey
    end if

    end block exchange

  end subroutine real8_3d_exchange

  subroutine real8_4d_exchange(ml,nex,j1,j2,i1,i2,k1,k2,n1,n2)
    implicit none
    real(rk8) , pointer , dimension(:,:,:,:) , intent(inout) :: ml
    integer(ik4) , intent(in) :: nex , j1 , j2  , i1 , i2 , k1 , k2 , n1 , n2
    integer(ik4) :: nx , ny , nk , nn
    integer(ik4) :: lb , rb
    integer(ik4) :: ndx , ndy , tx , ty , sizex , sizey

    nx = j2-j1+1
    ny = i2-i1+1
    nk = k2-k1+1
    nn = n2-n1+1
    lb = nex
    rb = nex
    if ( left == mpi_proc_null ) lb = 0
    if ( right == mpi_proc_null ) rb = 0
    tx = ny
    ty = nx+lb+rb
    sizex = nex*tx*nk*nn
    sizey = nex*ty*nk*nn
    ndx = 2*sizex
    ndy = 2*sizey

    exchange : block

    integer(ik4), dimension(4) :: counts , displs
    real(rk8), dimension(ndx) :: sdatax
    real(rk8), dimension(ndy) :: sdatay
    real(rk8), dimension(ndx), volatile :: rdatax
    real(rk8), dimension(ndy), volatile :: rdatay
    integer(ik4) :: ib1 , ib2 , iex , k , n

    ib2 = 0
    do n = n1 , n2
      do k = k1 , k2
        do iex = 1 , nex
          ib1 = ib2 + 1
          ib2 = ib1 + tx - 1
          sdatax(ib1:ib2) = ml(j1+(iex-1),i1:i2,k,n)
        end do
      end do
    end do
    do n = n1 , n2
      do k = k1 , k2
        do iex = 1 , nex
          ib1 = ib2 + 1
          ib2 = ib1 + tx - 1
          sdatax(ib1:ib2) = ml(j2-(iex-1),i1:i2,k,n)
        end do
      end do
    end do

    counts = [ sizex, sizex, 0, 0 ]
    displs = [ 0, sizex, 2*sizex, 2*sizex ]
    call mpi_neighbor_alltoallv(sdatax, counts, displs, mpi_real8, &
        rdatax, counts, displs, mpi_real8, cartesian_communicator, mpierr)

    ib2 = 0
    if ( left /= mpi_proc_null ) then
      do n = n1 , n2
        do k = k1 , k2
          do iex = 1 , nex
            ib1 = ib2 + 1
            ib2 = ib1 + tx - 1
            ml(j1-iex,i1:i2,k,n) = rdatax(ib1:ib2)
          end do
        end do
      end do
    else
      ib2 = ib2 + sizex
    end if
    if ( right /= mpi_proc_null ) then
      do n = n1 , n2
        do k = k1 , k2
          do iex = 1 , nex
            ib1 = ib2 + 1
            ib2 = ib1 + tx - 1
            ml(j2+iex,i1:i2,k,n) = rdatax(ib1:ib2)
          end do
        end do
      end do
    else
      ib2 = ib2 + sizex
    end if

    ib2 = 0
    do n = n1 , n2
      do k = k1 , k2
        do iex = 1 , nex
          ib1 = ib2 + 1
          ib2 = ib1 + ty - 1
          sdatay(ib1:ib2) = ml(j1-lb:j2+rb,i1+(iex-1),k,n)
        end do
      end do
    end do
    do n = n1 , n2
      do k = k1 , k2
        do iex = 1 , nex
          ib1 = ib2 + 1
          ib2 = ib1 + ty - 1
          sdatay(ib1:ib2) = ml(j1-lb:j2+rb,i2-(iex-1),k,n)
        end do
      end do
    end do
    counts = [ 0, 0, sizey, sizey ]
    displs = [ 0, 0, 0, sizey ]
    call mpi_neighbor_alltoallv(sdatay, counts, displs, mpi_real8, &
        rdatay, counts, displs, mpi_real8, cartesian_communicator, mpierr)

    ib2 = 0
    if ( bottom /= mpi_proc_null ) then
      do n = n1 , n2
        do k = k1 , k2
          do iex = 1 , nex
            ib1 = ib2 + 1
            ib2 = ib1 + ty - 1
            ml(j1-lb:j2+rb,i1-iex,k,n) = rdatay(ib1:ib2)
          end do
        end do
      end do
    else
      ib2 = ib2 + sizey
    end if
    if ( top /= mpi_proc_null ) then
      do n = n1 , n2
        do k = k1 , k2
          do iex = 1 , nex
            ib1 = ib2 + 1
            ib2 = ib1 + ty - 1
            ml(j1-lb:j2+rb,i2+iex,k,n) = rdatay(ib1:ib2)
          end do
        end do
      end do
    else
      ib2 = ib2 + sizey
    end if

    end block exchange

  end subroutine real8_4d_exchange

  subroutine real4_2d_exchange(ml,nex,j1,j2,i1,i2)
    implicit none
    real(rk4) , pointer , dimension(:,:) , intent(inout) :: ml
    integer(ik4) , intent(in) :: nex , j1 , j2  , i1 , i2
    integer(ik4) :: nx , ny
    integer(ik4) :: lb , rb
    integer(ik4) :: ndx , ndy , tx , ty , sizex , sizey

    nx = j2-j1+1
    ny = i2-i1+1
    lb = nex
    rb = nex
    if ( left  == mpi_proc_null ) lb = 0
    if ( right == mpi_proc_null ) rb = 0
    tx = ny
    ty = nx+lb+rb
    sizex = nex*tx
    sizey = nex*ty
    ndx = 2*sizex
    ndy = 2*sizey

    exchange : block

    integer(ik4), dimension(4) :: counts , displs
    real(rk4), dimension(ndx) :: sdatax
    real(rk4), dimension(ndy) :: sdatay
    real(rk4), dimension(ndx), volatile :: rdatax
    real(rk4), dimension(ndy), volatile :: rdatay
    integer(ik4) :: ib1 , ib2 , iex

    ib2 = 0
    do iex = 1 , nex
      ib1 = ib2 + 1
      ib2 = ib1 + tx - 1
      sdatax(ib1:ib2) = ml(j1+(iex-1),i1:i2)
    end do
    do iex = 1 , nex
      ib1 = ib2 + 1
      ib2 = ib1 + tx - 1
      sdatax(ib1:ib2) = ml(j2-(iex-1),i1:i2)
    end do

    counts = [ sizex, sizex, 0, 0 ]
    displs = [ 0, sizex, 2*sizex, 2*sizex ]
    call mpi_neighbor_alltoallv(sdatax, counts, displs, mpi_real4, &
        rdatax, counts, displs, mpi_real4, cartesian_communicator, mpierr)

    ib2 = 0
    if ( left /= mpi_proc_null ) then
      do iex = 1 , nex
        ib1 = ib2 + 1
        ib2 = ib1 + tx - 1
        ml(j1-iex,i1:i2) = rdatax(ib1:ib2)
      end do
    else
      ib2 = ib2 + sizex
    end if
    if ( right /= mpi_proc_null ) then
      do iex = 1 , nex
        ib1 = ib2 + 1
        ib2 = ib1 + tx - 1
        ml(j2+iex,i1:i2) = rdatax(ib1:ib2)
      end do
    else
      ib2 = ib2 + sizex
    end if

    ib2 = 0
    do iex = 1 , nex
      ib1 = ib2 + 1
      ib2 = ib1 + ty - 1
      sdatay(ib1:ib2) = ml(j1-lb:j2+rb,i1+(iex-1))
    end do
    do iex = 1 , nex
      ib1 = ib2 + 1
      ib2 = ib1 + ty - 1
      sdatay(ib1:ib2) = ml(j1-lb:j2+rb,i2-(iex-1))
    end do

    counts = [ 0, 0, sizey, sizey ]
    displs = [ 0, 0, 0, sizey ]
    call mpi_neighbor_alltoallv(sdatay, counts, displs, mpi_real4, &
        rdatay, counts, displs, mpi_real4, cartesian_communicator, mpierr)

    ib2 = 0
    if ( bottom /= mpi_proc_null ) then
      do iex = 1 , nex
        ib1 = ib2 + 1
        ib2 = ib1 + ty - 1
        ml(j1-lb:j2+rb,i1-iex) = rdatay(ib1:ib2)
      end do
    else
      ib2 = ib2 + sizey
    end if
    if ( top /= mpi_proc_null ) then
      do iex = 1 , nex
        ib1 = ib2 + 1
        ib2 = ib1 + ty - 1
        ml(j1-lb:j2+rb,i2+iex) = rdatay(ib1:ib2)
      end do
    else
      ib2 = ib2 + sizey
    end if

    end block exchange
  end subroutine real4_2d_exchange

  subroutine real4_3d_exchange(ml,nex,j1,j2,i1,i2,k1,k2)
    implicit none
    real(rk4) , pointer , dimension(:,:,:) , intent(inout) :: ml
    integer(ik4) , intent(in) :: nex , j1 , j2  , i1 , i2 , k1 , k2
    integer(ik4) :: nx , ny , nk
    integer(ik4) :: lb , rb
    integer(ik4) :: ndx , ndy , tx , ty , sizex , sizey

    nx = j2-j1+1
    ny = i2-i1+1
    nk = k2-k1+1
    lb = nex
    rb = nex
    if ( left  == mpi_proc_null ) lb = 0
    if ( right == mpi_proc_null ) rb = 0
    tx = ny
    ty = nx+lb+rb
    sizex = nex*tx*nk
    sizey = nex*ty*nk
    ndx = 2*sizex
    ndy = 2*sizey

    exchange : block

    integer(ik4), dimension(4) :: counts , displs
    real(rk4), dimension(ndx) :: sdatax
    real(rk4), dimension(ndy) :: sdatay
    real(rk4), dimension(ndx), volatile :: rdatax
    real(rk4), dimension(ndy), volatile :: rdatay
    integer(ik4) :: ib1 , ib2 , iex , k

    ib2 = 0
    do k = k1 , k2
      do iex = 1 , nex
        ib1 = ib2 + 1
        ib2 = ib1 + tx - 1
        sdatax(ib1:ib2) = ml(j1+(iex-1),i1:i2,k)
      end do
    end do
    do k = k1 , k2
      do iex = 1 , nex
        ib1 = ib2 + 1
        ib2 = ib1 + tx - 1
        sdatax(ib1:ib2) = ml(j2-(iex-1),i1:i2,k)
      end do
    end do

    counts = [ sizex, sizex, 0, 0 ]
    displs = [ 0, sizex, 2*sizex, 2*sizex ]
    call mpi_neighbor_alltoallv(sdatax, counts, displs, mpi_real4, &
        rdatax, counts, displs, mpi_real4, cartesian_communicator, mpierr)

    ib2 = 0
    if ( left /= mpi_proc_null ) then
      do k = k1 , k2
        do iex = 1 , nex
          ib1 = ib2 + 1
          ib2 = ib1 + tx - 1
          ml(j1-iex,i1:i2,k) = rdatax(ib1:ib2)
        end do
      end do
    else
      ib2 = ib2 + sizex
    end if
    if ( right /= mpi_proc_null ) then
      do k = k1 , k2
        do iex = 1 , nex
          ib1 = ib2 + 1
          ib2 = ib1 + tx - 1
          ml(j2+iex,i1:i2,k) = rdatax(ib1:ib2)
        end do
      end do
    else
      ib2 = ib2 + sizex
    end if

    ib2 = 0
    do k = k1 , k2
      do iex = 1 , nex
        ib1 = ib2 + 1
        ib2 = ib1 + ty - 1
        sdatay(ib1:ib2) = ml(j1-lb:j2+rb,i1+(iex-1),k)
      end do
    end do
    do k = k1 , k2
      do iex = 1 , nex
        ib1 = ib2 + 1
        ib2 = ib1 + ty - 1
        sdatay(ib1:ib2) = ml(j1-lb:j2+rb,i2-(iex-1),k)
      end do
    end do
    counts = [ 0, 0, sizey, sizey ]
    displs = [ 0, 0, 0, sizey ]
    call mpi_neighbor_alltoallv(sdatay, counts, displs, mpi_real8, &
        rdatay, counts, displs, mpi_real8, cartesian_communicator, mpierr)

    ib2 = 0
    if ( bottom /= mpi_proc_null ) then
      do k = k1 , k2
        do iex = 1 , nex
          ib1 = ib2 + 1
          ib2 = ib1 + ty - 1
          ml(j1-lb:j2+rb,i1-iex,k) = rdatay(ib1:ib2)
        end do
      end do
    else
      ib2 = ib2 + sizey
    end if
    if ( top /= mpi_proc_null ) then
      do k = k1 , k2
        do iex = 1 , nex
          ib1 = ib2 + 1
          ib2 = ib1 + ty - 1
          ml(j1-lb:j2+rb,i2+iex,k) = rdatay(ib1:ib2)
        end do
      end do
    else
      ib2 = ib2 + sizey
    end if

    end block exchange

  end subroutine real4_3d_exchange

  subroutine real4_4d_exchange(ml,nex,j1,j2,i1,i2,k1,k2,n1,n2)
    implicit none
    real(rk4) , pointer , dimension(:,:,:,:) , intent(inout) :: ml
    integer(ik4) , intent(in) :: nex , j1 , j2  , i1 , i2 , k1 , k2 , n1 , n2
    integer(ik4) :: nx , ny , nk , nn
    integer(ik4) :: lb , rb
    integer(ik4) :: ndx , ndy , tx , ty , sizex , sizey

    nx = j2-j1+1
    ny = i2-i1+1
    nk = k2-k1+1
    nn = n2-n1+1
    lb = nex
    rb = nex
    if ( left == mpi_proc_null ) lb = 0
    if ( right == mpi_proc_null ) rb = 0
    tx = ny
    ty = nx+lb+rb
    sizex = nex*tx*nk*nn
    sizey = nex*ty*nk*nn
    ndx = 2*sizex
    ndy = 2*sizey

    exchange : block

    integer(ik4), dimension(4) :: counts , displs
    real(rk4), dimension(ndx) :: sdatax
    real(rk4), dimension(ndy) :: sdatay
    real(rk4), dimension(ndx), volatile :: rdatax
    real(rk4), dimension(ndy), volatile :: rdatay
    integer(ik4) :: ib1 , ib2 , iex , k , n

    ib2 = 0
    do n = n1 , n2
      do k = k1 , k2
        do iex = 1 , nex
          ib1 = ib2 + 1
          ib2 = ib1 + tx - 1
          sdatax(ib1:ib2) = ml(j1+(iex-1),i1:i2,k,n)
        end do
      end do
    end do
    do n = n1 , n2
      do k = k1 , k2
        do iex = 1 , nex
          ib1 = ib2 + 1
          ib2 = ib1 + tx - 1
          sdatax(ib1:ib2) = ml(j2-(iex-1),i1:i2,k,n)
        end do
      end do
    end do

    counts = [ sizex, sizex, 0, 0 ]
    displs = [ 0, sizex, 2*sizex, 2*sizex ]
    call mpi_neighbor_alltoallv(sdatax, counts, displs, mpi_real4, &
        rdatax, counts, displs, mpi_real4, cartesian_communicator, mpierr)

    ib2 = 0
    if ( left /= mpi_proc_null ) then
      do n = n1 , n2
        do k = k1 , k2
          do iex = 1 , nex
            ib1 = ib2 + 1
            ib2 = ib1 + tx - 1
            ml(j1-iex,i1:i2,k,n) = rdatax(ib1:ib2)
          end do
        end do
      end do
    else
      ib2 = ib2 + sizex
    end if
    if ( right /= mpi_proc_null ) then
      do n = n1 , n2
        do k = k1 , k2
          do iex = 1 , nex
            ib1 = ib2 + 1
            ib2 = ib1 + tx - 1
            ml(j2+iex,i1:i2,k,n) = rdatax(ib1:ib2)
          end do
        end do
      end do
    else
      ib2 = ib2 + sizex
    end if

    ib2 = 0
    do n = n1 , n2
      do k = k1 , k2
        do iex = 1 , nex
          ib1 = ib2 + 1
          ib2 = ib1 + ty - 1
          sdatay(ib1:ib2) = ml(j1-lb:j2+rb,i1+(iex-1),k,n)
        end do
      end do
    end do
    do n = n1 , n2
      do k = k1 , k2
        do iex = 1 , nex
          ib1 = ib2 + 1
          ib2 = ib1 + ty - 1
          sdatay(ib1:ib2) = ml(j1-lb:j2+rb,i2-(iex-1),k,n)
        end do
      end do
    end do
    counts = [ 0, 0, sizey, sizey ]
    displs = [ 0, 0, 0, sizey ]
    call mpi_neighbor_alltoallv(sdatay, counts, displs, mpi_real4, &
        rdatay, counts, displs, mpi_real4, cartesian_communicator, mpierr)

    ib2 = 0
    if ( bottom /= mpi_proc_null ) then
      do n = n1 , n2
        do k = k1 , k2
          do iex = 1 , nex
            ib1 = ib2 + 1
            ib2 = ib1 + ty - 1
            ml(j1-lb:j2+rb,i1-iex,k,n) = rdatay(ib1:ib2)
          end do
        end do
      end do
    else
      ib2 = ib2 + sizey
    end if
    if ( top /= mpi_proc_null ) then
      do n = n1 , n2
        do k = k1 , k2
          do iex = 1 , nex
            ib1 = ib2 + 1
            ib2 = ib1 + ty - 1
            ml(j1-lb:j2+rb,i2+iex,k,n) = rdatay(ib1:ib2)
          end do
        end do
      end do
    else
      ib2 = ib2 + sizey
    end if

    end block exchange

  end subroutine real4_4d_exchange

  subroutine real8_2d_exchange_left_right_bottom_top(ml,nex,j1,j2,i1,i2)
    implicit none
    real(rk8) , pointer , dimension(:,:) , intent(inout) :: ml
    integer(ik4) , intent(in) :: nex , j1 , j2  , i1 , i2
    integer(ik4) :: ndx , ndy , nx , ny , tx , ty , sizex , sizey

    nx = j2-j1+1
    ny = i2-i1+1
    tx = ny
    ty = nx
    sizex = nex*tx
    sizey = nex*ty
    ndx = 2*sizex
    ndy = 2*sizey

    exchange : block

    integer(ik4), dimension(4) :: counts , displs
    real(rk8), dimension(ndx+ndy) :: sdata
    real(rk8), dimension(ndx+ndy), volatile :: rdata
    integer(ik4) :: ib1 , ib2 , iex

    ib2 = 0
    do iex = 1 , nex
      ib1 = ib2 + 1
      ib2 = ib1 + tx - 1
      sdata(ib1:ib2) = ml(j1+(iex-1),i1:i2)
    end do
    do iex = 1 , nex
      ib1 = ib2 + 1
      ib2 = ib1 + tx - 1
      sdata(ib1:ib2) = ml(j2-(iex-1),i1:i2)
    end do
    do iex = 1 , nex
      ib1 = ib2 + 1
      ib2 = ib1 + ty - 1
      sdata(ib1:ib2) = ml(j1:j2,i1+(iex-1))
    end do
    do iex = 1 , nex
      ib1 = ib2 + 1
      ib2 = ib1 + ty - 1
      sdata(ib1:ib2) = ml(j1:j2,i2-(iex-1))
    end do

    counts = [ sizex, sizex, sizey, sizey ]
    displs = [ 0, sizex, 2*sizex, 2*sizex+sizey ]
    call mpi_neighbor_alltoallv(sdata, counts, displs, mpi_real8, &
        rdata, counts, displs, mpi_real8, cartesian_communicator, mpierr)

    ib2 = 0
    if ( left /= mpi_proc_null ) then
      do iex = 1 , nex
        ib1 = ib2 + 1
        ib2 = ib1 + tx - 1
        ml(j1-iex,i1:i2) = rdata(ib1:ib2)
      end do
    else
      ib2 = ib2 + sizex
    end if
    if ( right /= mpi_proc_null ) then
      do iex = 1 , nex
        ib1 = ib2 + 1
        ib2 = ib1 + tx - 1
        ml(j2+iex,i1:i2) = rdata(ib1:ib2)
      end do
    else
      ib2 = ib2 + sizex
    end if
    if ( bottom /= mpi_proc_null ) then
      do iex = 1 , nex
        ib1 = ib2 + 1
        ib2 = ib1 + ty - 1
        ml(j1:j2,i1-iex) = rdata(ib1:ib2)
      end do
    else
      ib2 = ib2 + sizey
    end if
    if ( top /= mpi_proc_null ) then
      do iex = 1 , nex
        ib1 = ib2 + 1
        ib2 = ib1 + ty - 1
        ml(j1:j2,i2+iex) = rdata(ib1:ib2)
      end do
    else
      ib2 = ib2 + sizey
    end if

    end block exchange

  end subroutine real8_2d_exchange_left_right_bottom_top

  subroutine real8_3d_exchange_left_right_bottom_top(ml,nex,j1,j2,i1,i2,k1,k2)
    implicit none
    real(rk8) , pointer , dimension(:,:,:) , intent(inout) :: ml
    integer(ik4) , intent(in) :: nex , j1 , j2  , i1 , i2 , k1 , k2
    integer(ik4) :: ndx , ndy , nx , ny , nk , tx , ty , sizex , sizey

    nx = j2-j1+1
    ny = i2-i1+1
    nk = k2-k1+1
    tx = ny
    ty = nx
    sizex = nex*tx*nk
    sizey = nex*ty*nk
    ndx = 2*sizex
    ndy = 2*sizey

    exchange : block

    integer(ik4), dimension(4) :: counts , displs
    real(rk8), dimension(ndx+ndy) :: sdata
    real(rk8), dimension(ndx+ndy), volatile :: rdata
    integer(ik4) :: ib1 , ib2 , iex , k

    ib2 = 0
    do k = k1 , k2
      do iex = 1 , nex
        ib1 = ib2 + 1
        ib2 = ib1 + tx - 1
        sdata(ib1:ib2) = ml(j1+(iex-1),i1:i2,k)
      end do
    end do
    do k = k1 , k2
      do iex = 1 , nex
        ib1 = ib2 + 1
        ib2 = ib1 + tx - 1
        sdata(ib1:ib2) = ml(j2-(iex-1),i1:i2,k)
      end do
    end do
    do k = k1 , k2
      do iex = 1 , nex
        ib1 = ib2 + 1
        ib2 = ib1 + ty - 1
        sdata(ib1:ib2) = ml(j1:j2,i1+(iex-1),k)
      end do
    end do
    do k = k1 , k2
      do iex = 1 , nex
        ib1 = ib2 + 1
        ib2 = ib1 + ty - 1
        sdata(ib1:ib2) = ml(j1:j2,i2-(iex-1),k)
      end do
    end do

    counts = [ sizex, sizex, sizey, sizey ]
    displs = [ 0, sizex, 2*sizex, 2*sizex+sizey ]
    call mpi_neighbor_alltoallv(sdata, counts, displs, mpi_real8, &
        rdata, counts, displs, mpi_real8, cartesian_communicator, mpierr)

    ib2 = 0
    if ( left /= mpi_proc_null ) then
      do k = k1 , k2
        do iex = 1 , nex
          ib1 = ib2 + 1
          ib2 = ib1 + tx - 1
          ml(j1-iex,i1:i2,k) = rdata(ib1:ib2)
        end do
      end do
    else
      ib2 = ib2 + sizex
    end if
    if ( right /= mpi_proc_null ) then
      do k = k1 , k2
        do iex = 1 , nex
          ib1 = ib2 + 1
          ib2 = ib1 + tx - 1
          ml(j2+iex,i1:i2,k) = rdata(ib1:ib2)
        end do
      end do
    else
      ib2 = ib2 + sizex
    end if
    if ( bottom /= mpi_proc_null ) then
      do k = k1 , k2
        do iex = 1 , nex
          ib1 = ib2 + 1
          ib2 = ib1 + ty - 1
          ml(j1:j2,i1-iex,k) = rdata(ib1:ib2)
        end do
      end do
    else
      ib2 = ib2 + sizey
    end if
    if ( top /= mpi_proc_null ) then
      do k = k1 , k2
        do iex = 1 , nex
          ib1 = ib2 + 1
          ib2 = ib1 + ty - 1
          ml(j1:j2,i2+iex,k) = rdata(ib1:ib2)
        end do
      end do
    else
      ib2 = ib2 + sizey
    end if

    end block exchange

  end subroutine real8_3d_exchange_left_right_bottom_top

  subroutine real8_4d_exchange_left_right_bottom_top(ml,nex,j1,j2, &
                                                     i1,i2,k1,k2,n1,n2)
    implicit none
    real(rk8) , pointer , dimension(:,:,:,:) , intent(inout) :: ml
    integer(ik4) , intent(in) :: nex , j1 , j2  , i1 , i2 , k1 , k2 , n1 , n2
    integer(ik4) :: ndx , ndy , nx , ny , nk , nn , tx , ty , sizex , sizey

    nx = j2-j1+1
    ny = i2-i1+1
    nk = k2-k1+1
    nn = n2-n1+1
    tx = ny
    ty = nx
    sizex = nex*tx*nk*nn
    sizey = nex*ty*nk*nn
    ndx = 2*sizex
    ndy = 2*sizey

    exchange : block

    integer(ik4), dimension(4) :: counts , displs
    real(rk8), dimension(ndx+ndy) :: sdata
    real(rk8), dimension(ndx+ndy), volatile :: rdata
    integer(ik4) :: ib1 , ib2 , iex , k , n

    ib2 = 0
    do n = n1 , n2
      do k = k1 , k2
        do iex = 1 , nex
          ib1 = ib2 + 1
          ib2 = ib1 + tx - 1
          sdata(ib1:ib2) = ml(j1+(iex-1),i1:i2,k,n)
        end do
      end do
    end do
    do n = n1 , n2
      do k = k1 , k2
        do iex = 1 , nex
          ib1 = ib2 + 1
          ib2 = ib1 + tx - 1
          sdata(ib1:ib2) = ml(j2-(iex-1),i1:i2,k,n)
        end do
      end do
    end do
    do n = n1 , n2
      do k = k1 , k2
        do iex = 1 , nex
          ib1 = ib2 + 1
          ib2 = ib1 + ty - 1
          sdata(ib1:ib2) = ml(j1:j2,i1+(iex-1),k,n)
        end do
      end do
    end do
    do n = n1 , n2
      do k = k1 , k2
        do iex = 1 , nex
          ib1 = ib2 + 1
          ib2 = ib1 + ty - 1
          sdata(ib1:ib2) = ml(j1:j2,i2-(iex-1),k,n)
        end do
      end do
    end do

    counts = [ sizex, sizex, sizey, sizey ]
    displs = [ 0, sizex, 2*sizex, 2*sizex+sizey ]
    call mpi_neighbor_alltoallv(sdata, counts, displs, mpi_real8, &
        rdata, counts, displs, mpi_real8, cartesian_communicator, mpierr)

    ib2 = 0
    if ( left /= mpi_proc_null ) then
      do n = n1 , n2
        do k = k1 , k2
          do iex = 1 , nex
            ib1 = ib2 + 1
            ib2 = ib1 + tx - 1
            ml(j1-iex,i1:i2,k,n) = rdata(ib1:ib2)
          end do
        end do
      end do
    else
      ib2 = ib2 + sizex
    end if
    if ( right /= mpi_proc_null ) then
      do n = n1 , n2
        do k = k1 , k2
          do iex = 1 , nex
            ib1 = ib2 + 1
            ib2 = ib1 + tx - 1
            ml(j2+iex,i1:i2,k,n) = rdata(ib1:ib2)
          end do
        end do
      end do
    else
      ib2 = ib2 + sizex
    end if
    if ( bottom /= mpi_proc_null ) then
      do n = n1 , n2
        do k = k1 , k2
          do iex = 1 , nex
            ib1 = ib2 + 1
            ib2 = ib1 + ty - 1
            ml(j1:j2,i1-iex,k,n) = rdata(ib1:ib2)
          end do
        end do
      end do
    else
      ib2 = ib2 + sizey
    end if
    if ( top /= mpi_proc_null ) then
      do n = n1 , n2
        do k = k1 , k2
          do iex = 1 , nex
            ib1 = ib2 + 1
            ib2 = ib1 + ty - 1
            ml(j1:j2,i2+iex,k,n) = rdata(ib1:ib2)
          end do
        end do
      end do
    else
      ib2 = ib2 + sizey
    end if

    end block exchange

  end subroutine real8_4d_exchange_left_right_bottom_top

  subroutine real4_2d_exchange_left_right_bottom_top(ml,nex,j1,j2,i1,i2)
    implicit none
    real(rk4) , pointer , dimension(:,:) , intent(inout) :: ml
    integer(ik4) , intent(in) :: nex , j1 , j2  , i1 , i2
    integer(ik4) :: ndx , ndy , nx , ny , tx , ty , sizex , sizey

    nx = j2-j1+1
    ny = i2-i1+1
    tx = ny
    ty = nx
    sizex = nex*tx
    sizey = nex*ty
    ndx = 2*sizex
    ndy = 2*sizey

    exchange : block

    integer(ik4), dimension(4) :: counts , displs
    real(rk4), dimension(ndx+ndy) :: sdata
    real(rk4), dimension(ndx+ndy), volatile :: rdata
    integer(ik4) :: ib1 , ib2 , iex

    ib2 = 0
    do iex = 1 , nex
      ib1 = ib2 + 1
      ib2 = ib1 + tx - 1
      sdata(ib1:ib2) = ml(j1+(iex-1),i1:i2)
    end do
    do iex = 1 , nex
      ib1 = ib2 + 1
      ib2 = ib1 + tx - 1
      sdata(ib1:ib2) = ml(j2-(iex-1),i1:i2)
    end do
    do iex = 1 , nex
      ib1 = ib2 + 1
      ib2 = ib1 + ty - 1
      sdata(ib1:ib2) = ml(j1:j2,i1+(iex-1))
    end do
    do iex = 1 , nex
      ib1 = ib2 + 1
      ib2 = ib1 + ty - 1
      sdata(ib1:ib2) = ml(j1:j2,i2-(iex-1))
    end do

    counts = [ sizex, sizex, sizey, sizey ]
    displs = [ 0, sizex, 2*sizex, 2*sizex+sizey ]
    call mpi_neighbor_alltoallv(sdata, counts, displs, mpi_real4, &
        rdata, counts, displs, mpi_real4, cartesian_communicator, mpierr)

    ib2 = 0
    if ( left /= mpi_proc_null ) then
      do iex = 1 , nex
        ib1 = ib2 + 1
        ib2 = ib1 + tx - 1
        ml(j1-iex,i1:i2) = rdata(ib1:ib2)
      end do
    else
      ib2 = ib2 + sizex
    end if
    if ( right /= mpi_proc_null ) then
      do iex = 1 , nex
        ib1 = ib2 + 1
        ib2 = ib1 + tx - 1
        ml(j2+iex,i1:i2) = rdata(ib1:ib2)
      end do
    else
      ib2 = ib2 + sizex
    end if
    if ( bottom /= mpi_proc_null ) then
      do iex = 1 , nex
        ib1 = ib2 + 1
        ib2 = ib1 + ty - 1
        ml(j1:j2,i1-iex) = rdata(ib1:ib2)
      end do
    else
      ib2 = ib2 + sizey
    end if
    if ( top /= mpi_proc_null ) then
      do iex = 1 , nex
        ib1 = ib2 + 1
        ib2 = ib1 + ty - 1
        ml(j1:j2,i2+iex) = rdata(ib1:ib2)
      end do
    else
      ib2 = ib2 + sizey
    end if

    end block exchange

  end subroutine real4_2d_exchange_left_right_bottom_top

  subroutine real4_3d_exchange_left_right_bottom_top(ml,nex,j1,j2,i1,i2,k1,k2)
    implicit none
    real(rk4) , pointer , dimension(:,:,:) , intent(inout) :: ml
    integer(ik4) , intent(in) :: nex , j1 , j2  , i1 , i2 , k1 , k2
    integer(ik4) :: ndx , ndy , nx , ny , nk , tx , ty , sizex , sizey

    nx = j2-j1+1
    ny = i2-i1+1
    nk = k2-k1+1
    tx = ny
    ty = nx
    sizex = nex*tx*nk
    sizey = nex*ty*nk
    ndx = 2*sizex
    ndy = 2*sizey

    exchange : block

    integer(ik4), dimension(4) :: counts , displs
    real(rk4), dimension(ndx+ndy) :: sdata
    real(rk4), dimension(ndx+ndy), volatile :: rdata
    integer(ik4) :: ib1 , ib2 , iex , k

    ib2 = 0
    do k = k1 , k2
      do iex = 1 , nex
        ib1 = ib2 + 1
        ib2 = ib1 + tx - 1
        sdata(ib1:ib2) = ml(j1+(iex-1),i1:i2,k)
      end do
    end do
    do k = k1 , k2
      do iex = 1 , nex
        ib1 = ib2 + 1
        ib2 = ib1 + tx - 1
        sdata(ib1:ib2) = ml(j2-(iex-1),i1:i2,k)
      end do
    end do
    do k = k1 , k2
      do iex = 1 , nex
        ib1 = ib2 + 1
        ib2 = ib1 + ty - 1
        sdata(ib1:ib2) = ml(j1:j2,i1+(iex-1),k)
      end do
    end do
    do k = k1 , k2
      do iex = 1 , nex
        ib1 = ib2 + 1
        ib2 = ib1 + ty - 1
        sdata(ib1:ib2) = ml(j1:j2,i2-(iex-1),k)
      end do
    end do

    counts = [ sizex, sizex, sizey, sizey ]
    displs = [ 0, sizex, 2*sizex, 2*sizex+sizey ]
    call mpi_neighbor_alltoallv(sdata, counts, displs, mpi_real4, &
        rdata, counts, displs, mpi_real4, cartesian_communicator, mpierr)

    ib2 = 0
    if ( left /= mpi_proc_null ) then
      do k = k1 , k2
        do iex = 1 , nex
          ib1 = ib2 + 1
          ib2 = ib1 + tx - 1
          ml(j1-iex,i1:i2,k) = rdata(ib1:ib2)
        end do
      end do
    else
      ib2 = ib2 + sizex
    end if
    if ( right /= mpi_proc_null ) then
      do k = k1 , k2
        do iex = 1 , nex
          ib1 = ib2 + 1
          ib2 = ib1 + tx - 1
          ml(j2+iex,i1:i2,k) = rdata(ib1:ib2)
        end do
      end do
    else
      ib2 = ib2 + sizex
    end if
    if ( bottom /= mpi_proc_null ) then
      do k = k1 , k2
        do iex = 1 , nex
          ib1 = ib2 + 1
          ib2 = ib1 + ty - 1
          ml(j1:j2,i1-iex,k) = rdata(ib1:ib2)
        end do
      end do
    else
      ib2 = ib2 + sizey
    end if
    if ( top /= mpi_proc_null ) then
      do k = k1 , k2
        do iex = 1 , nex
          ib1 = ib2 + 1
          ib2 = ib1 + ty - 1
          ml(j1:j2,i2+iex,k) = rdata(ib1:ib2)
        end do
      end do
    else
      ib2 = ib2 + sizey
    end if

    end block exchange

  end subroutine real4_3d_exchange_left_right_bottom_top

  subroutine real4_4d_exchange_left_right_bottom_top(ml,nex,j1,j2, &
                                                     i1,i2,k1,k2,n1,n2)
    implicit none
    real(rk4) , pointer , dimension(:,:,:,:) , intent(inout) :: ml
    integer(ik4) , intent(in) :: nex , j1 , j2  , i1 , i2 , k1 , k2 , n1 , n2
    integer(ik4) :: ndx , ndy , nx , ny , nk , nn , tx , ty , sizex , sizey

    nx = j2-j1+1
    ny = i2-i1+1
    nk = k2-k1+1
    nn = n2-n1+1
    tx = ny
    ty = nx
    sizex = nex*tx*nk*nn
    sizey = nex*ty*nk*nn
    ndx = 2*sizex
    ndy = 2*sizey

    exchange : block

    integer(ik4), dimension(4) :: counts , displs
    real(rk4), dimension(ndx+ndy) :: sdata
    real(rk4), dimension(ndx+ndy), volatile :: rdata
    integer(ik4) :: ib1 , ib2 , iex , k , n

    ib2 = 0
    do n = n1 , n2
      do k = k1 , k2
        do iex = 1 , nex
          ib1 = ib2 + 1
          ib2 = ib1 + tx - 1
          sdata(ib1:ib2) = ml(j1+(iex-1),i1:i2,k,n)
        end do
      end do
    end do
    do n = n1 , n2
      do k = k1 , k2
        do iex = 1 , nex
          ib1 = ib2 + 1
          ib2 = ib1 + tx - 1
          sdata(ib1:ib2) = ml(j2-(iex-1),i1:i2,k,n)
        end do
      end do
    end do
    do n = n1 , n2
      do k = k1 , k2
        do iex = 1 , nex
          ib1 = ib2 + 1
          ib2 = ib1 + ty - 1
          sdata(ib1:ib2) = ml(j1:j2,i1+(iex-1),k,n)
        end do
      end do
    end do
    do n = n1 , n2
      do k = k1 , k2
        do iex = 1 , nex
          ib1 = ib2 + 1
          ib2 = ib1 + ty - 1
          sdata(ib1:ib2) = ml(j1:j2,i2-(iex-1),k,n)
        end do
      end do
    end do

    counts = [ sizex, sizex, sizey, sizey ]
    displs = [ 0, sizex, 2*sizex, 2*sizex+sizey ]
    call mpi_neighbor_alltoallv(sdata, counts, displs, mpi_real4, &
        rdata, counts, displs, mpi_real4, cartesian_communicator, mpierr)

    ib2 = 0
    if ( left /= mpi_proc_null ) then
      do n = n1 , n2
        do k = k1 , k2
          do iex = 1 , nex
            ib1 = ib2 + 1
            ib2 = ib1 + tx - 1
            ml(j1-iex,i1:i2,k,n) = rdata(ib1:ib2)
          end do
        end do
      end do
    else
      ib2 = ib2 + sizex
    end if
    if ( right /= mpi_proc_null ) then
      do n = n1 , n2
        do k = k1 , k2
          do iex = 1 , nex
            ib1 = ib2 + 1
            ib2 = ib1 + tx - 1
            ml(j2+iex,i1:i2,k,n) = rdata(ib1:ib2)
          end do
        end do
      end do
    else
      ib2 = ib2 + sizex
    end if
    if ( bottom /= mpi_proc_null ) then
      do n = n1 , n2
        do k = k1 , k2
          do iex = 1 , nex
            ib1 = ib2 + 1
            ib2 = ib1 + ty - 1
            ml(j1:j2,i1-iex,k,n) = rdata(ib1:ib2)
          end do
        end do
      end do
    else
      ib2 = ib2 + sizey
    end if
    if ( top /= mpi_proc_null ) then
      do n = n1 , n2
        do k = k1 , k2
          do iex = 1 , nex
            ib1 = ib2 + 1
            ib2 = ib1 + ty - 1
            ml(j1:j2,i2+iex,k,n) = rdata(ib1:ib2)
          end do
        end do
      end do
    else
      ib2 = ib2 + sizey
    end if

    end block exchange

  end subroutine real4_4d_exchange_left_right_bottom_top

  subroutine real8_2d_exchange_left_right(ml,nex,j1,j2,i1,i2)
    implicit none
    real(rk8) , pointer , dimension(:,:) , intent(inout) :: ml
    integer(ik4) , intent(in) :: nex , j1 , j2  , i1 , i2
    integer(ik4) :: ndx , ny , tx , sizex

    ny = i2-i1+1
    tx = ny
    sizex = nex*tx
    ndx = 2*sizex

    exchange : block

    integer(ik4), dimension(4) :: counts , displs
    real(rk8), dimension(ndx) :: sdata
    real(rk8), dimension(ndx), volatile :: rdata
    integer(ik4) :: ib1 , ib2 , iex

    ib2 = 0
    do iex = 1 , nex
      ib1 = ib2 + 1
      ib2 = ib1 + tx - 1
      sdata(ib1:ib2) = ml(j1+(iex-1),i1:i2)
    end do
    do iex = 1 , nex
      ib1 = ib2 + 1
      ib2 = ib1 + tx - 1
      sdata(ib1:ib2) = ml(j2-(iex-1),i1:i2)
    end do

    counts = [ sizex, sizex, 0, 0 ]
    displs = [ 0, sizex, 2*sizex, 2*sizex ]
    call mpi_neighbor_alltoallv(sdata, counts, displs, mpi_real8, &
        rdata, counts, displs, mpi_real8, cartesian_communicator, mpierr)

    ib2 = 0
    if ( left /= mpi_proc_null ) then
      do iex = 1 , nex
        ib1 = ib2 + 1
        ib2 = ib1 + tx - 1
        ml(j1-iex,i1:i2) = rdata(ib1:ib2)
      end do
    else
      ib2 = ib2 + sizex
    end if
    if ( right /= mpi_proc_null ) then
      do iex = 1 , nex
        ib1 = ib2 + 1
        ib2 = ib1 + tx - 1
        ml(j2+iex,i1:i2) = rdata(ib1:ib2)
      end do
    else
      ib2 = ib2 + sizex
    end if

    end block exchange

  end subroutine real8_2d_exchange_left_right

  subroutine real8_3d_exchange_left_right(ml,nex,j1,j2,i1,i2,k1,k2)
    implicit none
    real(rk8) , pointer , dimension(:,:,:) , intent(inout) :: ml
    integer(ik4) , intent(in) :: nex , j1 , j2  , i1 , i2 , k1 , k2
    integer(ik4) :: ndx , ny , nk , tx , sizex

    ny = i2-i1+1
    nk = k2-k1+1
    tx = ny
    sizex = nex*tx*nk
    ndx = 2*sizex

    exchange : block

    integer(ik4), dimension(4) :: counts , displs
    real(rk8), dimension(ndx) :: sdata
    real(rk8), dimension(ndx), volatile :: rdata
    integer(ik4) :: ib1 , ib2 , iex , k

    ib2 = 0
    do k = k1 , k2
      do iex = 1 , nex
        ib1 = ib2 + 1
        ib2 = ib1 + tx - 1
        sdata(ib1:ib2) = ml(j1+(iex-1),i1:i2,k)
      end do
    end do
    do k = k1 , k2
      do iex = 1 , nex
        ib1 = ib2 + 1
        ib2 = ib1 + tx - 1
        sdata(ib1:ib2) = ml(j2-(iex-1),i1:i2,k)
      end do
    end do

    counts = [ sizex, sizex, 0, 0 ]
    displs = [ 0, sizex, 2*sizex, 2*sizex ]
    call mpi_neighbor_alltoallv(sdata, counts, displs, mpi_real8, &
        rdata, counts, displs, mpi_real8, cartesian_communicator, mpierr)

    ib2 = 0
    if ( left /= mpi_proc_null ) then
      do k = k1 , k2
        do iex = 1 , nex
          ib1 = ib2 + 1
          ib2 = ib1 + tx - 1
          ml(j1-iex,i1:i2,k) = rdata(ib1:ib2)
        end do
      end do
    else
      ib2 = ib2 + sizex
    end if
    if ( right /= mpi_proc_null ) then
      do k = k1 , k2
        do iex = 1 , nex
          ib1 = ib2 + 1
          ib2 = ib1 + tx - 1
          ml(j2+iex,i1:i2,k) = rdata(ib1:ib2)
        end do
      end do
    else
      ib2 = ib2 + sizex
    end if

    end block exchange

  end subroutine real8_3d_exchange_left_right

  subroutine real8_4d_exchange_left_right(ml,nex,j1,j2,i1,i2,k1,k2,n1,n2)
    implicit none
    real(rk8) , pointer , dimension(:,:,:,:) , intent(inout) :: ml
    integer(ik4) , intent(in) :: nex , j1 , j2  , i1 , i2 , k1 , k2 , n1 , n2
    integer(ik4) :: ndx , ny , nk , nn , tx , sizex

    ny = i2-i1+1
    nk = k2-k1+1
    nn = n2-n1+1
    tx = ny
    sizex = nex*tx*nk*nn
    ndx = 2*sizex

    exchange : block

    integer(ik4), dimension(4) :: counts , displs
    real(rk8), dimension(ndx) :: sdata
    real(rk8), dimension(ndx), volatile :: rdata
    integer(ik4) :: ib1 , ib2 , iex , k , n

    ib2 = 0
    do n = n1 , n2
      do k = k1 , k2
        do iex = 1 , nex
          ib1 = ib2 + 1
          ib2 = ib1 + tx - 1
          sdata(ib1:ib2) = ml(j1+(iex-1),i1:i2,k,n)
        end do
      end do
    end do
    do n = n1 , n2
      do k = k1 , k2
        do iex = 1 , nex
          ib1 = ib2 + 1
          ib2 = ib1 + tx - 1
          sdata(ib1:ib2) = ml(j2-(iex-1),i1:i2,k,n)
        end do
      end do
    end do

    counts = [ sizex, sizex, 0, 0 ]
    displs = [ 0, sizex, 2*sizex, 2*sizex ]
    call mpi_neighbor_alltoallv(sdata, counts, displs, mpi_real8, &
        rdata, counts, displs, mpi_real8, cartesian_communicator, mpierr)

    ib2 = 0
    if ( left /= mpi_proc_null ) then
      do n = n1 , n2
        do k = k1 , k2
          do iex = 1 , nex
            ib1 = ib2 + 1
            ib2 = ib1 + tx - 1
            ml(j1-iex,i1:i2,k,n) = rdata(ib1:ib2)
          end do
        end do
      end do
    else
      ib2 = ib2 + sizex
    end if
    if ( right /= mpi_proc_null ) then
      do n = n1 , n2
        do k = k1 , k2
          do iex = 1 , nex
            ib1 = ib2 + 1
            ib2 = ib1 + tx - 1
            ml(j2+iex,i1:i2,k,n) = rdata(ib1:ib2)
          end do
        end do
      end do
    else
      ib2 = ib2 + sizex
    end if

    end block exchange

  end subroutine real8_4d_exchange_left_right

  subroutine real4_2d_exchange_left_right(ml,nex,j1,j2,i1,i2)
    implicit none
    real(rk4) , pointer , dimension(:,:) , intent(inout) :: ml
    integer(ik4) , intent(in) :: nex , j1 , j2  , i1 , i2
    integer(ik4) :: ndx , ny , tx , sizex

    ny = i2-i1+1
    tx = ny
    sizex = nex*tx
    ndx = 2*sizex

    exchange : block

    integer(ik4), dimension(4) :: counts , displs
    real(rk4), dimension(ndx) :: sdata
    real(rk4), dimension(ndx), volatile :: rdata
    integer(ik4) :: ib1 , ib2 , iex

    ib2 = 0
    do iex = 1 , nex
      ib1 = ib2 + 1
      ib2 = ib1 + tx - 1
      sdata(ib1:ib2) = ml(j1+(iex-1),i1:i2)
    end do
    do iex = 1 , nex
      ib1 = ib2 + 1
      ib2 = ib1 + tx - 1
      sdata(ib1:ib2) = ml(j2-(iex-1),i1:i2)
    end do

    counts = [ sizex, sizex, 0, 0 ]
    displs = [ 0, sizex, 2*sizex, 2*sizex ]
    call mpi_neighbor_alltoallv(sdata, counts, displs, mpi_real4, &
        rdata, counts, displs, mpi_real4, cartesian_communicator, mpierr)

    ib2 = 0
    if ( left /= mpi_proc_null ) then
      do iex = 1 , nex
        ib1 = ib2 + 1
        ib2 = ib1 + tx - 1
        ml(j1-iex,i1:i2) = rdata(ib1:ib2)
      end do
    else
      ib2 = ib2 + sizex
    end if
    if ( right /= mpi_proc_null ) then
      do iex = 1 , nex
        ib1 = ib2 + 1
        ib2 = ib1 + tx - 1
        ml(j2+iex,i1:i2) = rdata(ib1:ib2)
      end do
    else
      ib2 = ib2 + sizex
    end if

    end block exchange

  end subroutine real4_2d_exchange_left_right

  subroutine real4_3d_exchange_left_right(ml,nex,j1,j2,i1,i2,k1,k2)
    implicit none
    real(rk4) , pointer , dimension(:,:,:) , intent(inout) :: ml
    integer(ik4) , intent(in) :: nex , j1 , j2  , i1 , i2 , k1 , k2
    integer(ik4) :: ndx , ny , nk , tx , sizex

    ny = i2-i1+1
    nk = k2-k1+1
    tx = ny
    sizex = nex*tx*nk
    ndx = 2*sizex

    exchange : block

    integer(ik4), dimension(4) :: counts , displs
    real(rk4), dimension(ndx) :: sdata
    real(rk4), dimension(ndx), volatile :: rdata
    integer(ik4) :: ib1 , ib2 , iex , k

    ib2 = 0
    do k = k1 , k2
      do iex = 1 , nex
        ib1 = ib2 + 1
        ib2 = ib1 + tx - 1
        sdata(ib1:ib2) = ml(j1+(iex-1),i1:i2,k)
      end do
    end do
    do k = k1 , k2
      do iex = 1 , nex
        ib1 = ib2 + 1
        ib2 = ib1 + tx - 1
        sdata(ib1:ib2) = ml(j2-(iex-1),i1:i2,k)
      end do
    end do

    counts = [ sizex, sizex, 0, 0 ]
    displs = [ 0, sizex, 2*sizex, 2*sizex ]
    call mpi_neighbor_alltoallv(sdata, counts, displs, mpi_real4, &
        rdata, counts, displs, mpi_real4, cartesian_communicator, mpierr)

    ib2 = 0
    if ( left /= mpi_proc_null ) then
      do k = k1 , k2
        do iex = 1 , nex
          ib1 = ib2 + 1
          ib2 = ib1 + tx - 1
          ml(j1-iex,i1:i2,k) = rdata(ib1:ib2)
        end do
      end do
    else
      ib2 = ib2 + sizex
    end if
    if ( right /= mpi_proc_null ) then
      do k = k1 , k2
        do iex = 1 , nex
          ib1 = ib2 + 1
          ib2 = ib1 + tx - 1
          ml(j2+iex,i1:i2,k) = rdata(ib1:ib2)
        end do
      end do
    else
      ib2 = ib2 + sizex
    end if

    end block exchange

  end subroutine real4_3d_exchange_left_right

  subroutine real4_4d_exchange_left_right(ml,nex,j1,j2,i1,i2,k1,k2,n1,n2)
    implicit none
    real(rk4) , pointer , dimension(:,:,:,:) , intent(inout) :: ml
    integer(ik4) , intent(in) :: nex , j1 , j2  , i1 , i2 , k1 , k2 , n1 , n2
    integer(ik4) :: ndx , ny , nk , nn , tx , sizex

    ny = i2-i1+1
    nk = k2-k1+1
    nn = n2-n1+1
    tx = ny
    sizex = nex*tx*nk*nn
    ndx = 2*sizex

    exchange : block

    integer(ik4), dimension(4) :: counts , displs
    real(rk4), dimension(ndx) :: sdata
    real(rk4), dimension(ndx), volatile :: rdata
    integer(ik4) :: ib1 , ib2 , iex , k , n

    ib2 = 0
    do n = n1 , n2
      do k = k1 , k2
        do iex = 1 , nex
          ib1 = ib2 + 1
          ib2 = ib1 + tx - 1
          sdata(ib1:ib2) = ml(j1+(iex-1),i1:i2,k,n)
        end do
      end do
    end do
    do n = n1 , n2
      do k = k1 , k2
        do iex = 1 , nex
          ib1 = ib2 + 1
          ib2 = ib1 + tx - 1
          sdata(ib1:ib2) = ml(j2-(iex-1),i1:i2,k,n)
        end do
      end do
    end do

    counts = [ sizex, sizex, 0, 0 ]
    displs = [ 0, sizex, 2*sizex, 2*sizex ]
    call mpi_neighbor_alltoallv(sdata, counts, displs, mpi_real4, &
        rdata, counts, displs, mpi_real4, cartesian_communicator, mpierr)

    ib2 = 0
    if ( left /= mpi_proc_null ) then
      do n = n1 , n2
        do k = k1 , k2
          do iex = 1 , nex
            ib1 = ib2 + 1
            ib2 = ib1 + tx - 1
            ml(j1-iex,i1:i2,k,n) = rdata(ib1:ib2)
          end do
        end do
      end do
    else
      ib2 = ib2 + sizex
    end if
    if ( right /= mpi_proc_null ) then
      do n = n1 , n2
        do k = k1 , k2
          do iex = 1 , nex
            ib1 = ib2 + 1
            ib2 = ib1 + tx - 1
            ml(j2+iex,i1:i2,k,n) = rdata(ib1:ib2)
          end do
        end do
      end do
    else
      ib2 = ib2 + sizex
    end if

    end block exchange

  end subroutine real4_4d_exchange_left_right

  subroutine real8_2d_exchange_bottom_top(ml,nex,j1,j2,i1,i2)
    implicit none
    real(rk8) , pointer , dimension(:,:) , intent(inout) :: ml
    integer(ik4) , intent(in) :: nex , j1 , j2  , i1 , i2
    integer(ik4) :: ndy , nx , ty , sizey

    nx = j2-j1+1
    ty = nx
    sizey = nex*ty
    ndy = 2*sizey

    exchange : block

    integer(ik4), dimension(4) :: counts , displs
    real(rk8), dimension(ndy) :: sdata
    real(rk8), dimension(ndy), volatile :: rdata
    integer(ik4) :: ib1 , ib2 , iex

    ib2 = 0
    do iex = 1 , nex
      ib1 = ib2 + 1
      ib2 = ib1 + ty - 1
      sdata(ib1:ib2) = ml(j1:j2,i1+(iex-1))
    end do
    do iex = 1 , nex
      ib1 = ib2 + 1
      ib2 = ib1 + ty - 1
      sdata(ib1:ib2) = ml(j1:j2,i2-(iex-1))
    end do

    counts = [ 0, 0, sizey, sizey ]
    displs = [ 0, 0, 0, sizey ]
    call mpi_neighbor_alltoallv(sdata, counts, displs, mpi_real8, &
        rdata, counts, displs, mpi_real8, cartesian_communicator, mpierr)

    ib2 = 0
    if ( bottom /= mpi_proc_null ) then
      do iex = 1 , nex
        ib1 = ib2 + 1
        ib2 = ib1 + ty - 1
        ml(j1:j2,i1-iex) = rdata(ib1:ib2)
      end do
    else
      ib2 = ib2 + sizey
    end if
    if ( top /= mpi_proc_null ) then
      do iex = 1 , nex
        ib1 = ib2 + 1
        ib2 = ib1 + ty - 1
        ml(j1:j2,i2+iex) = rdata(ib1:ib2)
      end do
    else
      ib2 = ib2 + sizey
    end if

    end block exchange

  end subroutine real8_2d_exchange_bottom_top

  subroutine real8_3d_exchange_bottom_top(ml,nex,j1,j2,i1,i2,k1,k2)
    implicit none
    real(rk8) , pointer , dimension(:,:,:) , intent(inout) :: ml
    integer(ik4) , intent(in) :: nex , j1 , j2  , i1 , i2 , k1 , k2
    integer(ik4) :: ndy , nx , nk , ty , sizey

    nx = j2-j1+1
    nk = k2-k1+1
    ty = nx
    sizey = nex*ty*nk
    ndy = 2*sizey

    exchange : block

    integer(ik4), dimension(4) :: counts , displs
    real(rk8), dimension(ndy) :: sdata
    real(rk8), dimension(ndy), volatile :: rdata
    integer(ik4) :: ib1 , ib2 , iex , k

    ib2 = 0
    do k = k1 , k2
      do iex = 1 , nex
        ib1 = ib2 + 1
        ib2 = ib1 + ty - 1
        sdata(ib1:ib2) = ml(j1:j2,i1+(iex-1),k)
      end do
    end do
    do k = k1 , k2
      do iex = 1 , nex
        ib1 = ib2 + 1
        ib2 = ib1 + ty - 1
        sdata(ib1:ib2) = ml(j1:j2,i2-(iex-1),k)
      end do
    end do

    counts = [ 0, 0, sizey, sizey ]
    displs = [ 0, 0, 0, sizey ]
    call mpi_neighbor_alltoallv(sdata, counts, displs, mpi_real8, &
        rdata, counts, displs, mpi_real8, cartesian_communicator, mpierr)

    ib2 = 0
    if ( bottom /= mpi_proc_null ) then
      do k = k1 , k2
        do iex = 1 , nex
          ib1 = ib2 + 1
          ib2 = ib1 + ty - 1
          ml(j1:j2,i1-iex,k) = rdata(ib1:ib2)
        end do
      end do
    else
      ib2 = ib2 + sizey
    end if
    if ( top /= mpi_proc_null ) then
      do k = k1 , k2
        do iex = 1 , nex
          ib1 = ib2 + 1
          ib2 = ib1 + ty - 1
          ml(j1:j2,i2+iex,k) = rdata(ib1:ib2)
        end do
      end do
    else
      ib2 = ib2 + sizey
    end if

    end block exchange

  end subroutine real8_3d_exchange_bottom_top

  subroutine real8_4d_exchange_bottom_top(ml,nex,j1,j2,i1,i2,k1,k2,n1,n2)
    implicit none
    real(rk8) , pointer , dimension(:,:,:,:) , intent(inout) :: ml
    integer(ik4) , intent(in) :: nex , j1 , j2  , i1 , i2 , k1 , k2 , n1 , n2
    integer(ik4) :: ndy , nx , nk , nn , ty , sizey

    nx = j2-j1+1
    nk = k2-k1+1
    nn = n2-n1+1
    ty = nx
    sizey = nex*ty*nk*nn
    ndy = 2*sizey

    exchange : block

    integer(ik4), dimension(4) :: counts , displs
    real(rk8), dimension(ndy) :: sdata
    real(rk8), dimension(ndy), volatile :: rdata
    integer(ik4) :: ib1 , ib2 , iex , k , n

    ib2 = 0
    do n = n1 , n2
      do k = k1 , k2
        do iex = 1 , nex
          ib1 = ib2 + 1
          ib2 = ib1 + ty - 1
          sdata(ib1:ib2) = ml(j1:j2,i1+(iex-1),k,n)
        end do
      end do
    end do
    do n = n1 , n2
      do k = k1 , k2
        do iex = 1 , nex
          ib1 = ib2 + 1
          ib2 = ib1 + ty - 1
          sdata(ib1:ib2) = ml(j1:j2,i2-(iex-1),k,n)
        end do
      end do
    end do

    counts = [ 0, 0, sizey, sizey ]
    displs = [ 0, 0, 0, sizey ]
    call mpi_neighbor_alltoallv(sdata, counts, displs, mpi_real8, &
        rdata, counts, displs, mpi_real8, cartesian_communicator, mpierr)

    ib2 = 0
    if ( bottom /= mpi_proc_null ) then
      do n = n1 , n2
        do k = k1 , k2
          do iex = 1 , nex
            ib1 = ib2 + 1
            ib2 = ib1 + ty - 1
            ml(j1:j2,i1-iex,k,n) = rdata(ib1:ib2)
          end do
        end do
      end do
    else
      ib2 = ib2 + sizey
    end if
    if ( top /= mpi_proc_null ) then
      do n = n1 , n2
        do k = k1 , k2
          do iex = 1 , nex
            ib1 = ib2 + 1
            ib2 = ib1 + ty - 1
            ml(j1:j2,i2+iex,k,n) = rdata(ib1:ib2)
          end do
        end do
      end do
    else
      ib2 = ib2 + sizey
    end if

    end block exchange

  end subroutine real8_4d_exchange_bottom_top

  subroutine real4_2d_exchange_bottom_top(ml,nex,j1,j2,i1,i2)
    implicit none
    real(rk4) , pointer , dimension(:,:) , intent(inout) :: ml
    integer(ik4) , intent(in) :: nex , j1 , j2  , i1 , i2
    integer(ik4) :: ndy , nx , ty , sizey

    nx = j2-j1+1
    ty = nx
    sizey = nex*ty
    ndy = 2*sizey

    exchange : block

    integer(ik4), dimension(4) :: counts , displs
    real(rk4), dimension(ndy) :: sdata
    real(rk4), dimension(ndy), volatile :: rdata
    integer(ik4) :: ib1 , ib2 , iex

    ib2 = 0
    do iex = 1 , nex
      ib1 = ib2 + 1
      ib2 = ib1 + ty - 1
      sdata(ib1:ib2) = ml(j1:j2,i1+(iex-1))
    end do
    do iex = 1 , nex
      ib1 = ib2 + 1
      ib2 = ib1 + ty - 1
      sdata(ib1:ib2) = ml(j1:j2,i2-(iex-1))
    end do

    counts = [ 0, 0, sizey, sizey ]
    displs = [ 0, 0, 0, sizey ]
    call mpi_neighbor_alltoallv(sdata, counts, displs, mpi_real4, &
        rdata, counts, displs, mpi_real4, cartesian_communicator, mpierr)

    ib2 = 0
    if ( bottom /= mpi_proc_null ) then
      do iex = 1 , nex
        ib1 = ib2 + 1
        ib2 = ib1 + ty - 1
        ml(j1:j2,i1-iex) = rdata(ib1:ib2)
      end do
    else
      ib2 = ib2 + sizey
    end if
    if ( top /= mpi_proc_null ) then
      do iex = 1 , nex
        ib1 = ib2 + 1
        ib2 = ib1 + ty - 1
        ml(j1:j2,i2+iex) = rdata(ib1:ib2)
      end do
    else
      ib2 = ib2 + sizey
    end if

    end block exchange

  end subroutine real4_2d_exchange_bottom_top

  subroutine real4_3d_exchange_bottom_top(ml,nex,j1,j2,i1,i2,k1,k2)
    implicit none
    real(rk4) , pointer , dimension(:,:,:) , intent(inout) :: ml
    integer(ik4) , intent(in) :: nex , j1 , j2  , i1 , i2 , k1 , k2
    integer(ik4) :: ndy , nx , nk , ty , sizey

    nx = j2-j1+1
    nk = k2-k1+1
    ty = nx
    sizey = nex*ty*nk
    ndy = 2*sizey

    exchange : block

    integer(ik4), dimension(4) :: counts , displs
    real(rk4), dimension(ndy) :: sdata
    real(rk4), dimension(ndy), volatile :: rdata
    integer(ik4) :: ib1 , ib2 , iex , k

    ib2 = 0
    do k = k1 , k2
      do iex = 1 , nex
        ib1 = ib2 + 1
        ib2 = ib1 + ty - 1
        sdata(ib1:ib2) = ml(j1:j2,i1+(iex-1),k)
      end do
    end do
    do k = k1 , k2
      do iex = 1 , nex
        ib1 = ib2 + 1
        ib2 = ib1 + ty - 1
        sdata(ib1:ib2) = ml(j1:j2,i2-(iex-1),k)
      end do
    end do

    counts = [ 0, 0, sizey, sizey ]
    displs = [ 0, 0, 0, sizey ]
    call mpi_neighbor_alltoallv(sdata, counts, displs, mpi_real4, &
        rdata, counts, displs, mpi_real4, cartesian_communicator, mpierr)

    ib2 = 0
    if ( bottom /= mpi_proc_null ) then
      do k = k1 , k2
        do iex = 1 , nex
          ib1 = ib2 + 1
          ib2 = ib1 + ty - 1
          ml(j1:j2,i1-iex,k) = rdata(ib1:ib2)
        end do
      end do
    else
      ib2 = ib2 + sizey
    end if
    if ( top /= mpi_proc_null ) then
      do k = k1 , k2
        do iex = 1 , nex
          ib1 = ib2 + 1
          ib2 = ib1 + ty - 1
          ml(j1:j2,i2+iex,k) = rdata(ib1:ib2)
        end do
      end do
    else
      ib2 = ib2 + sizey
    end if

    end block exchange

  end subroutine real4_3d_exchange_bottom_top

  subroutine real4_4d_exchange_bottom_top(ml,nex,j1,j2,i1,i2,k1,k2,n1,n2)
    implicit none
    real(rk4) , pointer , dimension(:,:,:,:) , intent(inout) :: ml
    integer(ik4) , intent(in) :: nex , j1 , j2  , i1 , i2 , k1 , k2 , n1 , n2
    integer(ik4) :: ndy , nx , nk , nn , ty , sizey

    nx = j2-j1+1
    nk = k2-k1+1
    nn = n2-n1+1
    ty = nx
    sizey = nex*ty*nk*nn
    ndy = 2*sizey

    exchange : block

    integer(ik4), dimension(4) :: counts , displs
    real(rk4), dimension(ndy) :: sdata
    real(rk4), dimension(ndy), volatile :: rdata
    integer(ik4) :: ib1 , ib2 , iex , k , n

    ib2 = 0
    do n = n1 , n2
      do k = k1 , k2
        do iex = 1 , nex
          ib1 = ib2 + 1
          ib2 = ib1 + ty - 1
          sdata(ib1:ib2) = ml(j1:j2,i1+(iex-1),k,n)
        end do
      end do
    end do
    do n = n1 , n2
      do k = k1 , k2
        do iex = 1 , nex
          ib1 = ib2 + 1
          ib2 = ib1 + ty - 1
          sdata(ib1:ib2) = ml(j1:j2,i2-(iex-1),k,n)
        end do
      end do
    end do

    counts = [ 0, 0, sizey, sizey ]
    displs = [ 0, 0, 0, sizey ]
    call mpi_neighbor_alltoallv(sdata, counts, displs, mpi_real4, &
        rdata, counts, displs, mpi_real4, cartesian_communicator, mpierr)

    ib2 = 0
    if ( bottom /= mpi_proc_null ) then
      do n = n1 , n2
        do k = k1 , k2
          do iex = 1 , nex
            ib1 = ib2 + 1
            ib2 = ib1 + ty - 1
            ml(j1:j2,i1-iex,k,n) = rdata(ib1:ib2)
          end do
        end do
      end do
    else
      ib2 = ib2 + sizey
    end if
    if ( top /= mpi_proc_null ) then
      do n = n1 , n2
        do k = k1 , k2
          do iex = 1 , nex
            ib1 = ib2 + 1
            ib2 = ib1 + ty - 1
            ml(j1:j2,i2+iex,k,n) = rdata(ib1:ib2)
          end do
        end do
      end do
    else
      ib2 = ib2 + sizey
    end if

    end block exchange

  end subroutine real4_4d_exchange_bottom_top

  subroutine real8_2d_exchange_left_bottom(ml,nex,j1,j2,i1,i2)
    implicit none
    real(rk8) , pointer , dimension(:,:) , intent(inout) :: ml
    integer(ik4) , intent(in) :: nex , j1 , j2  , i1 , i2
    integer(ik4) :: nx , ny
    integer(ik4) :: lb
    integer(ik4) :: ndx , ndy , tx , ty , sizex , sizey

    nx = j2-j1+1
    ny = i2-i1+1
    lb = nex
    if ( left == mpi_proc_null ) lb = 0
    tx = ny
    ty = nx+lb
    sizex = nex*tx
    sizey = nex*ty
    ndx = 2*sizex
    ndy = 2*sizey

    exchange : block

    integer(ik4), dimension(4) :: counts , displs
    real(rk8), dimension(ndx) :: sdatax
    real(rk8), dimension(ndy) :: sdatay
    real(rk8), dimension(ndx), volatile :: rdatax
    real(rk8), dimension(ndy), volatile :: rdatay
    integer(ik4) :: ib1 , ib2 , iex

    ib2 = 0
    do iex = 1 , nex
      ib1 = ib2 + 1
      ib2 = ib1 + tx - 1
      sdatax(ib1:ib2) = ml(j1+(iex-1),i1:i2)
    end do
    do iex = 1 , nex
      ib1 = ib2 + 1
      ib2 = ib1 + tx - 1
      sdatax(ib1:ib2) = ml(j2-(iex-1),i1:i2)
    end do

    counts = [ sizex, sizex, 0, 0 ]
    displs = [ 0, sizex, 2*sizex, 2*sizex ]
    call mpi_neighbor_alltoallv(sdatax, counts, displs, mpi_real8, &
        rdatax, counts, displs, mpi_real8, cartesian_communicator, mpierr)

    ib2 = 0
    if ( left /= mpi_proc_null ) then
      do iex = 1 , nex
        ib1 = ib2 + 1
        ib2 = ib1 + tx - 1
        ml(j1-iex,i1:i2) = rdatax(ib1:ib2)
      end do
    else
      ib2 = ib2 + sizex
    end if

    ib2 = 0
    do iex = 1 , nex
      ib1 = ib2 + 1
      ib2 = ib1 + ty - 1
      sdatay(ib1:ib2) = ml(j1-lb:j2,i1+(iex-1))
    end do
    do iex = 1 , nex
      ib1 = ib2 + 1
      ib2 = ib1 + ty - 1
      sdatay(ib1:ib2) = ml(j1-lb:j2,i2-(iex-1))
    end do

    counts = [ 0, 0, sizey, sizey ]
    displs = [ 0, 0, 0, sizey ]
    call mpi_neighbor_alltoallv(sdatay, counts, displs, mpi_real8, &
        rdatay, counts, displs, mpi_real8, cartesian_communicator, mpierr)

    ib2 = 0
    if ( bottom /= mpi_proc_null ) then
      do iex = 1 , nex
        ib1 = ib2 + 1
        ib2 = ib1 + ty - 1
        ml(j1-lb:j2,i1-iex) = rdatay(ib1:ib2)
      end do
    else
      ib2 = ib2 + sizey
    end if

    end block exchange

  end subroutine real8_2d_exchange_left_bottom

  subroutine real8_3d_exchange_left_bottom(ml,nex,j1,j2,i1,i2,k1,k2)
    implicit none
    real(rk8) , pointer , dimension(:,:,:) , intent(inout) :: ml
    integer(ik4) , intent(in) :: nex , j1 , j2  , i1 , i2 , k1 , k2
    integer(ik4) :: nx , ny , nk
    integer(ik4) :: lb
    integer(ik4) :: ndx , ndy , tx , ty , sizex , sizey

    nx = j2-j1+1
    ny = i2-i1+1
    nk = k2-k1+1
    lb = nex
    if ( left == mpi_proc_null ) lb = 0
    tx = ny
    ty = nx+lb
    sizex = nex*tx*nk
    sizey = nex*ty*nk
    ndx = 2*sizex
    ndy = 2*sizey

    exchange : block

    integer(ik4), dimension(4) :: counts , displs
    real(rk8), dimension(ndx) :: sdatax
    real(rk8), dimension(ndy) :: sdatay
    real(rk8), dimension(ndx), volatile :: rdatax
    real(rk8), dimension(ndy), volatile :: rdatay
    integer(ik4) :: ib1 , ib2 , iex , k

    ib2 = 0
    do k = k1 , k2
      do iex = 1 , nex
        ib1 = ib2 + 1
        ib2 = ib1 + tx - 1
        sdatax(ib1:ib2) = ml(j1+(iex-1),i1:i2,k)
      end do
    end do
    do k = k1 , k2
      do iex = 1 , nex
        ib1 = ib2 + 1
        ib2 = ib1 + tx - 1
        sdatax(ib1:ib2) = ml(j2-(iex-1),i1:i2,k)
      end do
    end do

    counts = [ sizex, sizex, 0, 0 ]
    displs = [ 0, sizex, 2*sizex, 2*sizex ]
    call mpi_neighbor_alltoallv(sdatax, counts, displs, mpi_real8, &
        rdatax, counts, displs, mpi_real8, cartesian_communicator, mpierr)

    ib2 = 0
    if ( left /= mpi_proc_null ) then
      do k = k1 , k2
        do iex = 1 , nex
          ib1 = ib2 + 1
          ib2 = ib1 + tx - 1
          ml(j1-iex,i1:i2,k) = rdatax(ib1:ib2)
        end do
      end do
    else
      ib2 = ib2 + sizex
    end if

    ib2 = 0
    do k = k1 , k2
      do iex = 1 , nex
        ib1 = ib2 + 1
        ib2 = ib1 + ty - 1
        sdatay(ib1:ib2) = ml(j1-lb:j2,i1+(iex-1),k)
      end do
    end do
    do k = k1 , k2
      do iex = 1 , nex
        ib1 = ib2 + 1
        ib2 = ib1 + ty - 1
        sdatay(ib1:ib2) = ml(j1-lb:j2,i2-(iex-1),k)
      end do
    end do
    counts = [ 0, 0, sizey, sizey ]
    displs = [ 0, 0, 0, sizey ]
    call mpi_neighbor_alltoallv(sdatay, counts, displs, mpi_real8, &
        rdatay, counts, displs, mpi_real8, cartesian_communicator, mpierr)

    ib2 = 0
    if ( bottom /= mpi_proc_null ) then
      do k = k1 , k2
        do iex = 1 , nex
          ib1 = ib2 + 1
          ib2 = ib1 + ty - 1
          ml(j1-lb:j2,i1-iex,k) = rdatay(ib1:ib2)
        end do
      end do
    else
      ib2 = ib2 + sizey
    end if

    end block exchange

  end subroutine real8_3d_exchange_left_bottom

  subroutine real8_4d_exchange_left_bottom(ml,nex,j1,j2,i1,i2,k1,k2,n1,n2)
    implicit none
    real(rk8) , pointer , dimension(:,:,:,:) , intent(inout) :: ml
    integer(ik4) , intent(in) :: nex , j1 , j2  , i1 , i2 , k1 , k2 , n1 , n2
    integer(ik4) :: nx , ny , nk , nn
    integer(ik4) :: lb
    integer(ik4) :: ndx , ndy , tx , ty , sizex , sizey

    nx = j2-j1+1
    ny = i2-i1+1
    nk = k2-k1+1
    nn = n2-n1+1
    lb = nex
    if ( left == mpi_proc_null ) lb = 0
    tx = ny
    ty = nx+lb
    sizex = nex*tx*nk
    sizey = nex*ty*nk
    ndx = 2*sizex
    ndy = 2*sizey

    exchange : block

    integer(ik4), dimension(4) :: counts , displs
    real(rk8), dimension(ndx) :: sdatax
    real(rk8), dimension(ndy) :: sdatay
    real(rk8), dimension(ndx), volatile :: rdatax
    real(rk8), dimension(ndy), volatile :: rdatay
    integer(ik4) :: ib1 , ib2 , iex , k , n

    ib2 = 0
    do n = n1 , n2
      do k = k1 , k2
        do iex = 1 , nex
          ib1 = ib2 + 1
          ib2 = ib1 + tx - 1
          sdatax(ib1:ib2) = ml(j1+(iex-1),i1:i2,k,n)
        end do
      end do
    end do
    do n = n1 , n2
      do k = k1 , k2
        do iex = 1 , nex
          ib1 = ib2 + 1
          ib2 = ib1 + tx - 1
          sdatax(ib1:ib2) = ml(j2-(iex-1),i1:i2,k,n)
        end do
      end do
    end do

    counts = [ sizex, sizex, 0, 0 ]
    displs = [ 0, sizex, 2*sizex, 2*sizex ]
    call mpi_neighbor_alltoallv(sdatax, counts, displs, mpi_real8, &
        rdatax, counts, displs, mpi_real8, cartesian_communicator, mpierr)

    ib2 = 0
    if ( left /= mpi_proc_null ) then
      do n = n1 , n2
        do k = k1 , k2
          do iex = 1 , nex
            ib1 = ib2 + 1
            ib2 = ib1 + tx - 1
            ml(j1-iex,i1:i2,k,n) = rdatax(ib1:ib2)
          end do
        end do
      end do
    else
      ib2 = ib2 + sizex
    end if

    ib2 = 0
    do n = n1 , n2
      do k = k1 , k2
        do iex = 1 , nex
          ib1 = ib2 + 1
          ib2 = ib1 + ty - 1
          sdatay(ib1:ib2) = ml(j1-lb:j2,i1+(iex-1),k,n)
        end do
      end do
    end do
    do n = n1 , n2
      do k = k1 , k2
        do iex = 1 , nex
          ib1 = ib2 + 1
          ib2 = ib1 + ty - 1
          sdatay(ib1:ib2) = ml(j1-lb:j2,i2-(iex-1),k,n)
        end do
      end do
    end do

    counts = [ 0, 0, sizey, sizey ]
    displs = [ 0, 0, 0, sizey ]
    call mpi_neighbor_alltoallv(sdatay, counts, displs, mpi_real8, &
        rdatay, counts, displs, mpi_real8, cartesian_communicator, mpierr)

    ib2 = 0
    if ( bottom /= mpi_proc_null ) then
      do n = n1 , n2
        do k = k1 , k2
          do iex = 1 , nex
            ib1 = ib2 + 1
            ib2 = ib1 + ty - 1
            ml(j1-lb:j2,i1-iex,k,n) = rdatay(ib1:ib2)
          end do
        end do
      end do
    else
      ib2 = ib2 + sizey
    end if

    end block exchange

  end subroutine real8_4d_exchange_left_bottom

  subroutine real4_2d_exchange_left_bottom(ml,nex,j1,j2,i1,i2)
    implicit none
    real(rk4) , pointer , dimension(:,:) , intent(inout) :: ml
    integer(ik4) , intent(in) :: nex , j1 , j2  , i1 , i2
    integer(ik4) :: nx , ny
    integer(ik4) :: lb
    integer(ik4) :: ndx , ndy , tx , ty , sizex , sizey

    nx = j2-j1+1
    ny = i2-i1+1
    lb = nex
    if ( left == mpi_proc_null ) lb = 0
    tx = ny
    ty = nx+lb
    sizex = nex*tx
    sizey = nex*ty
    ndx = 2*sizex
    ndy = 2*sizey

    exchange : block

    integer(ik4), dimension(4) :: counts , displs
    real(rk4), dimension(ndx) :: sdatax
    real(rk4), dimension(ndy) :: sdatay
    real(rk4), dimension(ndx), volatile :: rdatax
    real(rk4), dimension(ndy), volatile :: rdatay
    integer(ik4) :: ib1 , ib2 , iex

    ib2 = 0
    do iex = 1 , nex
      ib1 = ib2 + 1
      ib2 = ib1 + tx - 1
      sdatax(ib1:ib2) = ml(j1+(iex-1),i1:i2)
    end do
    do iex = 1 , nex
      ib1 = ib2 + 1
      ib2 = ib1 + tx - 1
      sdatax(ib1:ib2) = ml(j2-(iex-1),i1:i2)
    end do

    counts = [ sizex, sizex, 0, 0 ]
    displs = [ 0, sizex, 2*sizex, 2*sizex ]
    call mpi_neighbor_alltoallv(sdatax, counts, displs, mpi_real4, &
        rdatax, counts, displs, mpi_real4, cartesian_communicator, mpierr)

    ib2 = 0
    if ( left /= mpi_proc_null ) then
      do iex = 1 , nex
        ib1 = ib2 + 1
        ib2 = ib1 + tx - 1
        ml(j1-iex,i1:i2) = rdatax(ib1:ib2)
      end do
    else
      ib2 = ib2 + sizex
    end if

    ib2 = 0
    do iex = 1 , nex
      ib1 = ib2 + 1
      ib2 = ib1 + ty - 1
      sdatay(ib1:ib2) = ml(j1-lb:j2,i1+(iex-1))
    end do
    do iex = 1 , nex
      ib1 = ib2 + 1
      ib2 = ib1 + ty - 1
      sdatay(ib1:ib2) = ml(j1-lb:j2,i2-(iex-1))
    end do

    counts = [ 0, 0, sizey, sizey ]
    displs = [ 0, 0, 0, sizey ]
    call mpi_neighbor_alltoallv(sdatay, counts, displs, mpi_real4, &
        rdatay, counts, displs, mpi_real4, cartesian_communicator, mpierr)

    ib2 = 0
    if ( bottom /= mpi_proc_null ) then
      do iex = 1 , nex
        ib1 = ib2 + 1
        ib2 = ib1 + ty - 1
        ml(j1-lb:j2,i1-iex) = rdatay(ib1:ib2)
      end do
    else
      ib2 = ib2 + sizey
    end if

    end block exchange

  end subroutine real4_2d_exchange_left_bottom

  subroutine real4_3d_exchange_left_bottom(ml,nex,j1,j2,i1,i2,k1,k2)
    implicit none
    real(rk4) , pointer , dimension(:,:,:) , intent(inout) :: ml
    integer(ik4) , intent(in) :: nex , j1 , j2  , i1 , i2 , k1 , k2
    integer(ik4) :: nx , ny , nk
    integer(ik4) :: lb
    integer(ik4) :: ndx , ndy , tx , ty , sizex , sizey

    nx = j2-j1+1
    ny = i2-i1+1
    nk = k2-k1+1
    lb = nex
    if ( left == mpi_proc_null ) lb = 0
    tx = ny
    ty = nx+lb
    sizex = nex*tx*nk
    sizey = nex*ty*nk
    ndx = 2*sizex
    ndy = 2*sizey

    exchange : block

    integer(ik4), dimension(4) :: counts , displs
    real(rk4), dimension(ndx) :: sdatax
    real(rk4), dimension(ndy) :: sdatay
    real(rk4), dimension(ndx), volatile :: rdatax
    real(rk4), dimension(ndy), volatile :: rdatay
    integer(ik4) :: ib1 , ib2 , iex , k

    ib2 = 0
    do k = k1 , k2
      do iex = 1 , nex
        ib1 = ib2 + 1
        ib2 = ib1 + tx - 1
        sdatax(ib1:ib2) = ml(j1+(iex-1),i1:i2,k)
      end do
    end do
    do k = k1 , k2
      do iex = 1 , nex
        ib1 = ib2 + 1
        ib2 = ib1 + tx - 1
        sdatax(ib1:ib2) = ml(j2-(iex-1),i1:i2,k)
      end do
    end do

    counts = [ sizex, sizex, 0, 0 ]
    displs = [ 0, sizex, 2*sizex, 2*sizex ]
    call mpi_neighbor_alltoallv(sdatax, counts, displs, mpi_real4, &
        rdatax, counts, displs, mpi_real4, cartesian_communicator, mpierr)

    ib2 = 0
    if ( left /= mpi_proc_null ) then
      do k = k1 , k2
        do iex = 1 , nex
          ib1 = ib2 + 1
          ib2 = ib1 + tx - 1
          ml(j1-iex,i1:i2,k) = rdatax(ib1:ib2)
        end do
      end do
    else
      ib2 = ib2 + sizex
    end if

    ib2 = 0
    do k = k1 , k2
      do iex = 1 , nex
        ib1 = ib2 + 1
        ib2 = ib1 + ty - 1
        sdatay(ib1:ib2) = ml(j1-lb:j2,i1+(iex-1),k)
      end do
    end do
    do k = k1 , k2
      do iex = 1 , nex
        ib1 = ib2 + 1
        ib2 = ib1 + ty - 1
        sdatay(ib1:ib2) = ml(j1-lb:j2,i2-(iex-1),k)
      end do
    end do

    counts = [ 0, 0, sizey, sizey ]
    displs = [ 0, 0, 0, sizey ]
    call mpi_neighbor_alltoallv(sdatay, counts, displs, mpi_real4, &
        rdatay, counts, displs, mpi_real4, cartesian_communicator, mpierr)

    ib2 = 0
    if ( bottom /= mpi_proc_null ) then
      do k = k1 , k2
        do iex = 1 , nex
          ib1 = ib2 + 1
          ib2 = ib1 + ty - 1
          ml(j1-lb:j2,i1-iex,k) = rdatay(ib1:ib2)
        end do
      end do
    else
      ib2 = ib2 + sizey
    end if

    end block exchange

  end subroutine real4_3d_exchange_left_bottom

  subroutine real4_4d_exchange_left_bottom(ml,nex,j1,j2,i1,i2,k1,k2,n1,n2)
    implicit none
    real(rk4) , pointer , dimension(:,:,:,:) , intent(inout) :: ml
    integer(ik4) , intent(in) :: nex , j1 , j2  , i1 , i2 , k1 , k2 , n1 , n2
    integer(ik4) :: nx , ny , nk , nn
    integer(ik4) :: lb
    integer(ik4) :: ndx , ndy , tx , ty , sizex , sizey

    nx = j2-j1+1
    ny = i2-i1+1
    nk = k2-k1+1
    nn = n2-n1+1
    lb = nex
    if ( left == mpi_proc_null ) lb = 0
    tx = ny
    ty = nx+lb
    sizex = nex*tx*nk
    sizey = nex*ty*nk
    ndx = 2*sizex
    ndy = 2*sizey

    exchange : block

    integer(ik4), dimension(4) :: counts , displs
    real(rk4), dimension(ndx) :: sdatax
    real(rk4), dimension(ndy) :: sdatay
    real(rk4), dimension(ndx), volatile :: rdatax
    real(rk4), dimension(ndy), volatile :: rdatay
    integer(ik4) :: ib1 , ib2 , iex , k , n

    ib2 = 0
    do n = n1 , n2
      do k = k1 , k2
        do iex = 1 , nex
          ib1 = ib2 + 1
          ib2 = ib1 + tx - 1
          sdatax(ib1:ib2) = ml(j1+(iex-1),i1:i2,k,n)
        end do
      end do
    end do
    do n = n1 , n2
      do k = k1 , k2
        do iex = 1 , nex
          ib1 = ib2 + 1
          ib2 = ib1 + tx - 1
          sdatax(ib1:ib2) = ml(j2-(iex-1),i1:i2,k,n)
        end do
      end do
    end do

    counts = [ sizex, sizex, 0, 0 ]
    displs = [ 0, sizex, 2*sizex, 2*sizex ]
    call mpi_neighbor_alltoallv(sdatax, counts, displs, mpi_real4, &
        rdatax, counts, displs, mpi_real4, cartesian_communicator, mpierr)

    ib2 = 0
    if ( left /= mpi_proc_null ) then
      do n = n1 , n2
        do k = k1 , k2
          do iex = 1 , nex
            ib1 = ib2 + 1
            ib2 = ib1 + tx - 1
            ml(j1-iex,i1:i2,k,n) = rdatax(ib1:ib2)
          end do
        end do
      end do
    else
      ib2 = ib2 + sizex
    end if

    ib2 = 0
    do n = n1 , n2
      do k = k1 , k2
        do iex = 1 , nex
          ib1 = ib2 + 1
          ib2 = ib1 + ty - 1
          sdatay(ib1:ib2) = ml(j1-lb:j2,i1+(iex-1),k,n)
        end do
      end do
    end do
    do n = n1 , n2
      do k = k1 , k2
        do iex = 1 , nex
          ib1 = ib2 + 1
          ib2 = ib1 + ty - 1
          sdatay(ib1:ib2) = ml(j1-lb:j2,i2-(iex-1),k,n)
        end do
      end do
    end do

    counts = [ 0, 0, sizey, sizey ]
    displs = [ 0, 0, 0, sizey ]
    call mpi_neighbor_alltoallv(sdatay, counts, displs, mpi_real4, &
        rdatay, counts, displs, mpi_real4, cartesian_communicator, mpierr)

    ib2 = 0
    if ( bottom /= mpi_proc_null ) then
      do n = n1 , n2
        do k = k1 , k2
          do iex = 1 , nex
            ib1 = ib2 + 1
            ib2 = ib1 + ty - 1
            ml(j1-lb:j2,i1-iex,k,n) = rdatay(ib1:ib2)
          end do
        end do
      end do
    else
      ib2 = ib2 + sizey
    end if

    end block exchange

  end subroutine real4_4d_exchange_left_bottom

  subroutine real8_2d_exchange_right_top(ml,nex,j1,j2,i1,i2)
    implicit none
    real(rk8) , pointer , dimension(:,:) , intent(inout) :: ml
    integer(ik4) , intent(in) :: nex , j1 , j2  , i1 , i2
    integer(ik4) :: nx , ny
    integer(ik4) :: rb
    integer(ik4) :: ndx , ndy , tx , ty , sizex , sizey

    nx = j2-j1+1
    ny = i2-i1+1
    rb = nex
    if ( right == mpi_proc_null ) rb = 0
    tx = ny
    ty = nx+rb
    sizex = nex*tx
    sizey = nex*ty
    ndx = 2*sizex
    ndy = 2*sizey

    exchange : block

    integer(ik4), dimension(4) :: counts , displs
    real(rk8), dimension(ndx) :: sdatax
    real(rk8), dimension(ndy) :: sdatay
    real(rk8), dimension(ndx), volatile :: rdatax
    real(rk8), dimension(ndy), volatile :: rdatay
    integer(ik4) :: ib1 , ib2 , iex

    ib2 = 0
    do iex = 1 , nex
      ib1 = ib2 + 1
      ib2 = ib1 + tx - 1
      sdatax(ib1:ib2) = ml(j1+(iex-1),i1:i2)
    end do
    do iex = 1 , nex
      ib1 = ib2 + 1
      ib2 = ib1 + tx - 1
      sdatax(ib1:ib2) = ml(j2-(iex-1),i1:i2)
    end do

    counts = [ sizex, sizex, 0, 0 ]
    displs = [ 0, sizex, 2*sizex, 2*sizex ]
    call mpi_neighbor_alltoallv(sdatax, counts, displs, mpi_real8, &
        rdatax, counts, displs, mpi_real8, cartesian_communicator, mpierr)

    ib2 = sizex
    if ( right /= mpi_proc_null ) then
      do iex = 1 , nex
        ib1 = ib2 + 1
        ib2 = ib1 + tx - 1
        ml(j2+iex,i1:i2) = rdatax(ib1:ib2)
      end do
    else
      ib2 = ib2 + sizex
    end if

    ib2 = 0
    do iex = 1 , nex
      ib1 = ib2 + 1
      ib2 = ib1 + ty - 1
      sdatay(ib1:ib2) = ml(j1:j2+rb,i1+(iex-1))
    end do
    do iex = 1 , nex
      ib1 = ib2 + 1
      ib2 = ib1 + ty - 1
      sdatay(ib1:ib2) = ml(j1:j2+rb,i2-(iex-1))
    end do

    counts = [ 0, 0, sizey, sizey ]
    displs = [ 0, 0, 0, sizey ]
    call mpi_neighbor_alltoallv(sdatay, counts, displs, mpi_real8, &
        rdatay, counts, displs, mpi_real8, cartesian_communicator, mpierr)

    ib2 = sizey
    if ( top /= mpi_proc_null ) then
      do iex = 1 , nex
        ib1 = ib2 + 1
        ib2 = ib1 + ty - 1
        ml(j1:j2+rb,i2+iex) = rdatay(ib1:ib2)
      end do
    else
      ib2 = ib2 + sizey
    end if

    end block exchange

  end subroutine real8_2d_exchange_right_top

  subroutine real8_3d_exchange_right_top(ml,nex,j1,j2,i1,i2,k1,k2)
    implicit none
    real(rk8) , pointer , dimension(:,:,:) , intent(inout) :: ml
    integer(ik4) , intent(in) :: nex , j1 , j2  , i1 , i2 , k1 , k2
    integer(ik4) :: nx , ny , nk
    integer(ik4) :: rb
    integer(ik4) :: ndx , ndy , tx , ty , sizex , sizey

    nx = j2-j1+1
    ny = i2-i1+1
    nk = k2-k1+1
    rb = nex
    if ( right == mpi_proc_null ) rb = 0
    tx = ny
    ty = nx+rb
    sizex = nex*tx*nk
    sizey = nex*ty*nk
    ndx = 2*sizex
    ndy = 2*sizey

    exchange : block

    integer(ik4), dimension(4) :: counts , displs
    real(rk8), dimension(ndx) :: sdatax
    real(rk8), dimension(ndy) :: sdatay
    real(rk8), dimension(ndx), volatile :: rdatax
    real(rk8), dimension(ndy), volatile :: rdatay
    integer(ik4) :: ib1 , ib2 , iex , k

    ib2 = 0
    do k = k1 , k2
      do iex = 1 , nex
        ib1 = ib2 + 1
        ib2 = ib1 + tx - 1
        sdatax(ib1:ib2) = ml(j1+(iex-1),i1:i2,k)
      end do
    end do
    do k = k1 , k2
      do iex = 1 , nex
        ib1 = ib2 + 1
        ib2 = ib1 + tx - 1
        sdatax(ib1:ib2) = ml(j2-(iex-1),i1:i2,k)
      end do
    end do

    counts = [ sizex, sizex, 0, 0 ]
    displs = [ 0, sizex, 2*sizex, 2*sizex ]
    call mpi_neighbor_alltoallv(sdatax, counts, displs, mpi_real8, &
        rdatax, counts, displs, mpi_real8, cartesian_communicator, mpierr)

    ib2 = sizex
    if ( right /= mpi_proc_null ) then
      do k = k1 , k2
        do iex = 1 , nex
          ib1 = ib2 + 1
          ib2 = ib1 + tx - 1
          ml(j2+iex,i1:i2,k) = rdatax(ib1:ib2)
        end do
      end do
    else
      ib2 = ib2 + sizex
    end if

    ib2 = 0
    do k = k1 , k2
      do iex = 1 , nex
        ib1 = ib2 + 1
        ib2 = ib1 + ty - 1
        sdatay(ib1:ib2) = ml(j1:j2+rb,i1+(iex-1),k)
      end do
    end do
    do k = k1 , k2
      do iex = 1 , nex
        ib1 = ib2 + 1
        ib2 = ib1 + ty - 1
        sdatay(ib1:ib2) = ml(j1:j2+rb,i2-(iex-1),k)
      end do
    end do

    counts = [ 0, 0, sizey, sizey ]
    displs = [ 0, 0, 0, sizey ]
    call mpi_neighbor_alltoallv(sdatay, counts, displs, mpi_real8, &
        rdatay, counts, displs, mpi_real8, cartesian_communicator, mpierr)

    ib2 = sizey
    if ( top /= mpi_proc_null ) then
      do k = k1 , k2
        do iex = 1 , nex
          ib1 = ib2 + 1
          ib2 = ib1 + ty - 1
          ml(j1:j2+rb,i2+iex,k) = rdatay(ib1:ib2)
        end do
      end do
    else
      ib2 = ib2 + sizey
    end if

    end block exchange

  end subroutine real8_3d_exchange_right_top

  subroutine real8_4d_exchange_right_top(ml,nex,j1,j2,i1,i2,k1,k2,n1,n2)
    implicit none
    real(rk8) , pointer , dimension(:,:,:,:) , intent(inout) :: ml
    integer(ik4) , intent(in) :: nex , j1 , j2  , i1 , i2 , k1 , k2 , n1 , n2
    integer(ik4) :: nx , ny , nk , nn
    integer(ik4) :: rb
    integer(ik4) :: ndx , ndy , tx , ty , sizex , sizey

    nx = j2-j1+1
    ny = i2-i1+1
    nk = k2-k1+1
    nn = n2-n1+1
    rb = nex
    if ( right == mpi_proc_null ) rb = 0
    tx = ny
    ty = nx+rb
    sizex = nex*tx*nk*nn
    sizey = nex*ty*nk*nn
    ndx = 2*sizex
    ndy = 2*sizey

    exchange : block

    integer(ik4), dimension(4) :: counts , displs
    real(rk8), dimension(ndx) :: sdatax
    real(rk8), dimension(ndy) :: sdatay
    real(rk8), dimension(ndx), volatile :: rdatax
    real(rk8), dimension(ndy), volatile :: rdatay
    integer(ik4) :: ib1 , ib2 , iex , k , n

    ib2 = 0
    do n = n1 , n2
      do k = k1 , k2
        do iex = 1 , nex
          ib1 = ib2 + 1
          ib2 = ib1 + tx - 1
          sdatax(ib1:ib2) = ml(j1+(iex-1),i1:i2,k,n)
        end do
      end do
    end do
    do n = n1 , n2
      do k = k1 , k2
        do iex = 1 , nex
          ib1 = ib2 + 1
          ib2 = ib1 + tx - 1
          sdatax(ib1:ib2) = ml(j2-(iex-1),i1:i2,k,n)
        end do
      end do
    end do

    counts = [ sizex, sizex, 0, 0 ]
    displs = [ 0, sizex, 2*sizex, 2*sizex ]
    call mpi_neighbor_alltoallv(sdatax, counts, displs, mpi_real8, &
        rdatax, counts, displs, mpi_real8, cartesian_communicator, mpierr)

    ib2 = sizex
    if ( right /= mpi_proc_null ) then
      do n = n1 , n2
        do k = k1 , k2
          do iex = 1 , nex
            ib1 = ib2 + 1
            ib2 = ib1 + tx - 1
            ml(j2+iex,i1:i2,k,n) = rdatax(ib1:ib2)
          end do
        end do
      end do
    else
      ib2 = ib2 + sizex
    end if

    ib2 = 0
    do n = n1 , n2
      do k = k1 , k2
        do iex = 1 , nex
          ib1 = ib2 + 1
          ib2 = ib1 + ty - 1
          sdatay(ib1:ib2) = ml(j1:j2+rb,i1+(iex-1),k,n)
        end do
      end do
    end do
    do n = n1 , n2
      do k = k1 , k2
        do iex = 1 , nex
          ib1 = ib2 + 1
          ib2 = ib1 + ty - 1
          sdatay(ib1:ib2) = ml(j1:j2+rb,i2-(iex-1),k,n)
        end do
      end do
    end do

    counts = [ 0, 0, sizey, sizey ]
    displs = [ 0, 0, 0, sizey ]
    call mpi_neighbor_alltoallv(sdatay, counts, displs, mpi_real8, &
        rdatay, counts, displs, mpi_real8, cartesian_communicator, mpierr)

    ib2 = sizey
    if ( top /= mpi_proc_null ) then
      do n = n1 , n2
        do k = k1 , k2
          do iex = 1 , nex
            ib1 = ib2 + 1
            ib2 = ib1 + ty - 1
            ml(j1:j2+rb,i2+iex,k,n) = rdatay(ib1:ib2)
          end do
        end do
      end do
    else
      ib2 = ib2 + sizey
    end if

    end block exchange

  end subroutine real8_4d_exchange_right_top

  subroutine real4_2d_exchange_right_top(ml,nex,j1,j2,i1,i2)
    implicit none
    real(rk4) , pointer , dimension(:,:) , intent(inout) :: ml
    integer(ik4) , intent(in) :: nex , j1 , j2  , i1 , i2
    integer(ik4) :: nx , ny
    integer(ik4) :: rb
    integer(ik4) :: ndx , ndy , tx , ty , sizex , sizey

    nx = j2-j1+1
    ny = i2-i1+1
    rb = nex
    if ( right == mpi_proc_null ) rb = 0
    tx = ny
    ty = nx+rb
    sizex = nex*tx
    sizey = nex*ty
    ndx = 2*sizex
    ndy = 2*sizey

    exchange : block

    integer(ik4), dimension(4) :: counts , displs
    real(rk4), dimension(ndx) :: sdatax
    real(rk4), dimension(ndy) :: sdatay
    real(rk4), dimension(ndx), volatile :: rdatax
    real(rk4), dimension(ndy), volatile :: rdatay
    integer(ik4) :: ib1 , ib2 , iex

    ib2 = 0
    do iex = 1 , nex
      ib1 = ib2 + 1
      ib2 = ib1 + tx - 1
      sdatax(ib1:ib2) = ml(j1+(iex-1),i1:i2)
    end do
    do iex = 1 , nex
      ib1 = ib2 + 1
      ib2 = ib1 + tx - 1
      sdatax(ib1:ib2) = ml(j2-(iex-1),i1:i2)
    end do

    counts = [ sizex, sizex, 0, 0 ]
    displs = [ 0, sizex, 2*sizex, 2*sizex ]
    call mpi_neighbor_alltoallv(sdatax, counts, displs, mpi_real4, &
        rdatax, counts, displs, mpi_real4, cartesian_communicator, mpierr)

    ib2 = sizex
    if ( right /= mpi_proc_null ) then
      do iex = 1 , nex
        ib1 = ib2 + 1
        ib2 = ib1 + tx - 1
        ml(j2+iex,i1:i2) = rdatax(ib1:ib2)
      end do
    else
      ib2 = ib2 + sizex
    end if

    ib2 = 0
    do iex = 1 , nex
      ib1 = ib2 + 1
      ib2 = ib1 + ty - 1
      sdatay(ib1:ib2) = ml(j1:j2+rb,i1+(iex-1))
    end do
    do iex = 1 , nex
      ib1 = ib2 + 1
      ib2 = ib1 + ty - 1
      sdatay(ib1:ib2) = ml(j1:j2+rb,i2-(iex-1))
    end do

    counts = [ 0, 0, sizey, sizey ]
    displs = [ 0, 0, 0, sizey ]
    call mpi_neighbor_alltoallv(sdatay, counts, displs, mpi_real4, &
        rdatay, counts, displs, mpi_real4, cartesian_communicator, mpierr)

    ib2 = sizey
    if ( top /= mpi_proc_null ) then
      do iex = 1 , nex
        ib1 = ib2 + 1
        ib2 = ib1 + ty - 1
        ml(j1:j2+rb,i2+iex) = rdatay(ib1:ib2)
      end do
    else
      ib2 = ib2 + sizey
    end if

    end block exchange

  end subroutine real4_2d_exchange_right_top

  subroutine real4_3d_exchange_right_top(ml,nex,j1,j2,i1,i2,k1,k2)
    implicit none
    real(rk4) , pointer , dimension(:,:,:) , intent(inout) :: ml
    integer(ik4) , intent(in) :: nex , j1 , j2  , i1 , i2 , k1 , k2
    integer(ik4) :: nx , ny , nk
    integer(ik4) :: rb
    integer(ik4) :: ndx , ndy , tx , ty , sizex , sizey

    nx = j2-j1+1
    ny = i2-i1+1
    nk = k2-k1+1
    rb = nex
    if ( right == mpi_proc_null ) rb = 0
    tx = ny
    ty = nx+rb
    sizex = nex*tx*nk
    sizey = nex*ty*nk
    ndx = 2*sizex
    ndy = 2*sizey

    exchange : block

    integer(ik4), dimension(4) :: counts , displs
    real(rk4), dimension(ndx) :: sdatax
    real(rk4), dimension(ndy) :: sdatay
    real(rk4), dimension(ndx), volatile :: rdatax
    real(rk4), dimension(ndy), volatile :: rdatay
    integer(ik4) :: ib1 , ib2 , iex , k

    ib2 = 0
    do k = k1 , k2
      do iex = 1 , nex
        ib1 = ib2 + 1
        ib2 = ib1 + tx - 1
        sdatax(ib1:ib2) = ml(j1+(iex-1),i1:i2,k)
      end do
    end do
    do k = k1 , k2
      do iex = 1 , nex
        ib1 = ib2 + 1
        ib2 = ib1 + tx - 1
        sdatax(ib1:ib2) = ml(j2-(iex-1),i1:i2,k)
      end do
    end do

    counts = [ sizex, sizex, 0, 0 ]
    displs = [ 0, sizex, 2*sizex, 2*sizex ]
    call mpi_neighbor_alltoallv(sdatax, counts, displs, mpi_real4, &
        rdatax, counts, displs, mpi_real4, cartesian_communicator, mpierr)

    ib2 = sizex
    if ( right /= mpi_proc_null ) then
      do k = k1 , k2
        do iex = 1 , nex
          ib1 = ib2 + 1
          ib2 = ib1 + tx - 1
          ml(j2+iex,i1:i2,k) = rdatax(ib1:ib2)
        end do
      end do
    else
      ib2 = ib2 + sizex
    end if

    ib2 = 0
    do k = k1 , k2
      do iex = 1 , nex
        ib1 = ib2 + 1
        ib2 = ib1 + ty - 1
        sdatay(ib1:ib2) = ml(j1:j2+rb,i1+(iex-1),k)
      end do
    end do
    do k = k1 , k2
      do iex = 1 , nex
        ib1 = ib2 + 1
        ib2 = ib1 + ty - 1
        sdatay(ib1:ib2) = ml(j1:j2+rb,i2-(iex-1),k)
      end do
    end do

    counts = [ 0, 0, sizey, sizey ]
    displs = [ 0, 0, 0, sizey ]
    call mpi_neighbor_alltoallv(sdatay, counts, displs, mpi_real4, &
        rdatay, counts, displs, mpi_real4, cartesian_communicator, mpierr)

    ib2 = sizey
    if ( top /= mpi_proc_null ) then
      do k = k1 , k2
        do iex = 1 , nex
          ib1 = ib2 + 1
          ib2 = ib1 + ty - 1
          ml(j1:j2+rb,i2+iex,k) = rdatay(ib1:ib2)
        end do
      end do
    else
      ib2 = ib2 + sizey
    end if

    end block exchange

  end subroutine real4_3d_exchange_right_top

  subroutine real4_4d_exchange_right_top(ml,nex,j1,j2,i1,i2,k1,k2,n1,n2)
    implicit none
    real(rk4) , pointer , dimension(:,:,:,:) , intent(inout) :: ml
    integer(ik4) , intent(in) :: nex , j1 , j2  , i1 , i2 , k1 , k2 , n1 , n2
    integer(ik4) :: nx , ny , nk , nn
    integer(ik4) :: rb
    integer(ik4) :: ndx , ndy , tx , ty , sizex , sizey

    nx = j2-j1+1
    ny = i2-i1+1
    nk = k2-k1+1
    nn = n2-n1+1
    rb = nex
    if ( right == mpi_proc_null ) rb = 0
    tx = ny
    ty = nx+rb
    sizex = nex*tx*nk*nn
    sizey = nex*ty*nk*nn
    ndx = 2*sizex
    ndy = 2*sizey

    exchange : block

    integer(ik4), dimension(4) :: counts , displs
    real(rk4), dimension(ndx) :: sdatax
    real(rk4), dimension(ndy) :: sdatay
    real(rk4), dimension(ndx), volatile :: rdatax
    real(rk4), dimension(ndy), volatile :: rdatay
    integer(ik4) :: ib1 , ib2 , iex , k , n

    ib2 = 0
    do n = n1 , n2
      do k = k1 , k2
        do iex = 1 , nex
          ib1 = ib2 + 1
          ib2 = ib1 + tx - 1
          sdatax(ib1:ib2) = ml(j1+(iex-1),i1:i2,k,n)
        end do
      end do
    end do
    do n = n1 , n2
      do k = k1 , k2
        do iex = 1 , nex
          ib1 = ib2 + 1
          ib2 = ib1 + tx - 1
          sdatax(ib1:ib2) = ml(j2-(iex-1),i1:i2,k,n)
        end do
      end do
    end do

    counts = [ sizex, sizex, 0, 0 ]
    displs = [ 0, sizex, 2*sizex, 2*sizex ]
    call mpi_neighbor_alltoallv(sdatax, counts, displs, mpi_real4, &
        rdatax, counts, displs, mpi_real4, cartesian_communicator, mpierr)

    ib2 = sizex
    if ( right /= mpi_proc_null ) then
      do n = n1 , n2
        do k = k1 , k2
          do iex = 1 , nex
            ib1 = ib2 + 1
            ib2 = ib1 + tx - 1
            ml(j2+iex,i1:i2,k,n) = rdatax(ib1:ib2)
          end do
        end do
      end do
    else
      ib2 = ib2 + sizex
    end if

    ib2 = 0
    do n = n1 , n2
      do k = k1 , k2
        do iex = 1 , nex
          ib1 = ib2 + 1
          ib2 = ib1 + ty - 1
          sdatay(ib1:ib2) = ml(j1:j2+rb,i1+(iex-1),k,n)
        end do
      end do
    end do
    do n = n1 , n2
      do k = k1 , k2
        do iex = 1 , nex
          ib1 = ib2 + 1
          ib2 = ib1 + ty - 1
          sdatay(ib1:ib2) = ml(j1:j2+rb,i2-(iex-1),k,n)
        end do
      end do
    end do

    counts = [ 0, 0, sizey, sizey ]
    displs = [ 0, 0, 0, sizey ]
    call mpi_neighbor_alltoallv(sdatay, counts, displs, mpi_real4, &
        rdatay, counts, displs, mpi_real4, cartesian_communicator, mpierr)

    ib2 = sizey
    if ( top /= mpi_proc_null ) then
      do n = n1 , n2
        do k = k1 , k2
          do iex = 1 , nex
            ib1 = ib2 + 1
            ib2 = ib1 + ty - 1
            ml(j1:j2+rb,i2+iex,k,n) = rdatay(ib1:ib2)
          end do
        end do
      end do
    else
      ib2 = ib2 + sizey
    end if

    end block exchange

  end subroutine real4_4d_exchange_right_top

  ! Takes u and v tendencies on the cross grid (as t, qv, qc, etc.)
  ! and interpolates the u and v to the dot grid.
  ! This routine sheilds the user of the function from the need to worry
  ! about the details of the domain decomposition.
  !
  ! Written by Travis A. O'Brien 01/04/11.
  !
  subroutine uvcross2dot(ux,vx,ud,vd)
    implicit none
    real(rkx) , pointer , dimension(:,:,:) , intent(inout) :: ux , vx
    real(rkx) , pointer , dimension(:,:,:) , intent(inout) :: ud , vd
    integer(ik4) :: i , j , k

    call exchange(ux,1,jce1,jce2,ice1,ice2,1,kz)
    call exchange(vx,1,jce1,jce2,ice1,ice2,1,kz)

    !
    !  o     o     o     o     o     o     o
    !
    !     x     x     x     x     x     x
    !
    !  o     o     o     o     o     o     o
    !         (i-1,j-1)     (i,j-1)
    !     x     x     x-----x     x     x
    !                 |(i,j)|
    !  o     o     o  |  o  |  o     o     o
    !                 |     |
    !     x     x     x-----x     x     x
    !           (i-1,j)     (i,j)
    !  o     o     o     o     o     o     o
    !
    !     x     x     x     x     x     x
    !
    !  o     o     o     o     o     o     o

    ! Perform the bilinear interpolation necessary
    ! to put the u and v variables on the dot grid.

    do concurrent ( j = jdi1:jdi2, i = idi1:idi2, k = 1:kz )
      ud(j,i,k) =  d_rfour*(ux(j,i,k) + ux(j-1,i,k) +   &
                            ux(j,i-1,k) + ux(j-1,i-1,k))
      vd(j,i,k) =  d_rfour*(vx(j,i,k) + vx(j-1,i,k) +   &
                            vx(j,i-1,k) + vx(j-1,i-1,k))
    end do
  end subroutine uvcross2dot

  subroutine uvdot2cross(ud,vd,ux,vx)
    implicit none
    real(rkx) , pointer , dimension(:,:,:) , intent(inout) :: ud , vd
    real(rkx) , pointer , dimension(:,:,:) , intent(inout) :: ux , vx
    integer(ik4) :: i , j , k

    call exchange(ud,1,jdi1,jdi2,idi1,idi2,1,kz)
    call exchange(vd,1,jdi1,jdi2,idi1,idi2,1,kz)

    !
    !     o     o     o     o     o     o
    !
    !        x     x     x     x     x
    !           (i+1,j)   (i+1,j+1)
    !     o     o     o-----o     o     o
    !                 |(i,j)|
    !        x     x  |  x  |  x     x
    !                 |     |
    !     o     o     o-----o     o     o
    !             (i,j)     (i,j+1)
    !        x     x     x     x     x
    !
    !     o     o     o     o     o     o
    !

    ! Perform the bilinear interpolation necessary
    ! to put the u and v variables on the cross grid.

    do concurrent ( j = jci1:jci2, i = ici1:ici2, k = 1:kz )
      ux(j,i,k) =  d_rfour*(ud(j,i,  k) + ud(j+1,i,  k) + &
                            ud(j,i+1,k) + ud(j+1,i+1,k))
      vx(j,i,k) =  d_rfour*(vd(j,i  ,k) + vd(j+1,i  ,k) + &
                            vd(j,i+1,k) + vd(j+1,i+1,k))
    end do
  end subroutine uvdot2cross

  subroutine cross2dot2d(x,d)
    implicit none
    real(rkx) , pointer , dimension(:,:) , intent(inout) :: x
    real(rkx) , pointer , dimension(:,:) , intent(inout) :: d
    integer(ik4) :: i , j

    call exchange(x,1,jci1,jci2,ici1,ici2)

    do concurrent ( j = jdii1:jdii2, i = idii1:idii2 )
      d(j,i) =  d_rfour*(x(j,i)   + x(j-1,i)   + &
                           x(j,i-1) + x(j-1,i-1))
    end do
  end subroutine cross2dot2d

  subroutine cross2dot3d(x,d)
    implicit none
    real(rkx) , pointer , dimension(:,:,:) , intent(inout) :: x
    real(rkx) , pointer , dimension(:,:,:) , intent(inout) :: d
    integer(ik4) :: i , j , k

    call exchange(x,1,jci1,jci2,ici1,ici2,1,kz)

    do concurrent ( j = jdii1:jdii2, i = idii1:idii2, k = 1:kz )
      d(j,i,k) =  d_rfour*(x(j,i,k)   + x(j-1,i,k)   + &
                           x(j,i-1,k) + x(j-1,i-1,k))
    end do
  end subroutine cross2dot3d

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
    if ( myid /= iocpu ) return
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
    call grid_collect(xvar%val,xvar%iobuf, &
                      xvar%mynx1,xvar%mynx2,xvar%myny1,xvar%myny2)
    if ( myid == iocpu ) then
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
    end if
    xvar%irec = xvar%irec + 1
  end subroutine grid_nc_write_var2d

  subroutine grid_nc_destroy_var2d(xvar)
    implicit none
    type (grid_nc_var2d) , intent(inout) :: xvar
    integer(ik4) :: istat
    if ( myid == iocpu ) then
      istat = nf90_close(xvar%ncid)
      if ( istat /= nf90_noerr ) then
        write(error_unit, *) nf90_strerror(istat)
        return
      end if
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
    if ( myid /= iocpu ) return
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
    call grid_collect(xvar%val,xvar%iobuf, &
                      xvar%mynx1,xvar%mynx2,xvar%myny1,xvar%myny2,1,xvar%nz)
    if ( myid == iocpu ) then
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
    end if
    xvar%irec = xvar%irec + 1
  end subroutine grid_nc_write_var3d

  subroutine grid_nc_destroy_var3d(xvar)
    implicit none
    type (grid_nc_var3d) , intent(inout) :: xvar
    integer(ik4) :: istat
    if ( myid == iocpu ) then
      istat = nf90_close(xvar%ncid)
      if ( istat /= nf90_noerr ) then
        write(error_unit, *) nf90_strerror(istat)
        return
      end if
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
    if ( myid /= iocpu ) return
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
    call grid_collect(xvar%val,xvar%iobuf, &
                      xvar%mynx1,xvar%mynx2,xvar%myny1,xvar%myny2, &
                      1,xvar%nz,1,xvar%nl)
    if ( myid == iocpu ) then
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
    end if
    xvar%irec = xvar%irec + 1
  end subroutine grid_nc_write_var4d

  subroutine grid_nc_destroy_var4d(xvar)
    implicit none
    type (grid_nc_var4d) , intent(inout) :: xvar
    integer(ik4) :: istat
    if ( myid == iocpu ) then
      istat = nf90_close(xvar%ncid)
      if ( istat /= nf90_noerr ) then
        write(error_unit, *) nf90_strerror(istat)
        return
      end if
    end if
    xvar%ncid = -1
    xvar%irec = -1
    nullify(xvar%val)
    call relmem4d(xvar%iobuf)
  end subroutine grid_nc_destroy_var4d

  subroutine gather_r(f_collect,f_sub)
    implicit none
    real(rk8) , dimension(:) , intent(out) :: f_collect
    real(rk8) , intent(in) :: f_sub
    real(rk8) , dimension(1) :: tmp
    tmp(1) = f_sub
    call mpi_gather(tmp,      1,mpi_real8, &
                    f_collect,1,mpi_real8,iocpu,mycomm,mpierr)
  end subroutine gather_r

  subroutine gather_i(i_collect,i_sub)
    implicit none
    integer(ik4) , dimension(:) , intent(out) :: i_collect
    integer(ik4) , intent(in) :: i_sub
    integer(ik4) , dimension(1) :: tmp
    tmp(1) = i_sub
    call mpi_gather(tmp,      1,mpi_integer4, &
                    i_collect,1,mpi_integer4,iocpu,mycomm,mpierr)
  end subroutine gather_i

  subroutine allgather_r(f_collect,f_sub)
    implicit none
    real(rk8) , dimension(:) , intent(out) :: f_collect
    real(rk8) , intent(in) :: f_sub
    real(rk8) , dimension(1) :: tmp
    tmp(1) = f_sub
    call mpi_allgather(tmp,      1,mpi_real8, &
                       f_collect,1,mpi_real8,mycomm,mpierr)
  end subroutine allgather_r

  subroutine allgather_i(i_collect,i_sub)
    implicit none
    integer(ik4) , dimension(:) , intent(out) :: i_collect
    integer(ik4) , intent(in) :: i_sub
    integer(ik4) , dimension(1) :: tmp
    tmp(1) = i_sub
    call mpi_allgather(tmp,      1,mpi_integer4, &
                       i_collect,1,mpi_integer4,mycomm,mpierr)
  end subroutine allgather_i

  subroutine allsync
    implicit none
    call mpi_barrier(mycomm,mpierr)
  end subroutine allsync

end module mod_mppparam
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
