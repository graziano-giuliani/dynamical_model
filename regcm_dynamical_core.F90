program regcm_dynamical_core
  use, intrinsic :: iso_fortran_env
  use mpi
  use mod_atmosphere
  use mod_intkinds
  use mod_dynparam
  use mod_runparams
  use mod_regcm_types
  use mod_realkinds
  use mod_constants
  use mod_mppparam
  use mod_moloch
  use mod_memutil
  use mod_zita
  implicit none

  logical :: lprint
  integer(ik4) :: i
  integer(ik4) :: ierr
  integer(ik4) :: iprov
  integer(ik4) , dimension(8) :: tval
  real(rkx) :: start_time , finish_time , last_time
  character (len=32) :: cdata='?'
  character (len=5) :: czone='?'
  type(grid_nc_var2d) :: psout
  type(grid_nc_var3d) :: tout
  character(len=*) , parameter :: f99001 = &
        '(2x," GIT Revision: ",a," compiled at: data : ",a,"  time: ",a,/)'

  call mpi_init_thread(mpi_thread_single,iprov,ierr)
  if ( ierr /= mpi_success ) then
    write(error_unit,*) 'Cannot initilize MPI'
    stop
  end if

  call mpi_comm_dup(MPI_COMM_WORLD,mycomm,ierr)
  if ( ierr /= 0 ) then
    call fatal(__FILE__,__LINE__,'Cannot get communicator!')
  end if
  call mpi_comm_rank(mycomm, myid, ierr)
  if ( ierr /= 0 ) then
    call fatal(__FILE__,__LINE__,'mpi_comm_rank Failure!')
  end if
  call mpi_comm_size(mycomm, nproc, ierr)
  if ( ierr /= 0 ) then
    call fatal(__FILE__,__LINE__,'mpi_comm_size Failure!')
  end if

  if ( myid == iocpu ) then
    call cpu_time(start_time)
    last_time = start_time
    write (output_unit,"(/,2x,'This is RegCM development version')")
    write (output_unit,f99001)  '1.0', __DATE__ , __TIME__
  end if

  call memory_init
  call initparam
  call set_nproc

  call setup_model_indexes

  call param
  call init

  if ( myid == iocpu ) print *, 'Starting Integration!'

  do i = 1 , nstep
    call output(i)
    lprint = mod(i,iprint) == 0
    call moloch(i,lprint)
  end do

  call memory_destroy
  call finalize
  call mpi_finalize(ierr)

  contains

    subroutine param
      implicit none
      integer(ik4) :: i , j , k
      real(rkx) :: dl
      real(rkx) :: clat = 0.0
      real(rkx) :: clon = 0.0
      ds = 3.0
      dx = ds * 1000.0_rkx
      rdx = 1.0_rkx/dx
      dt = 30.0_rkx
      rdt = 1.0_rkx/dt
      mo_ztop = 36000_rkx
      mo_h = 12000_rkx
      mo_a0 = 0.0_rkx
      mo_nadv = 3
      mo_nsound = 5
      mo_nzfilt = kz/5
      mo_anu2 = 0.6_rkx
      nqx = 5
      iqfrst = iqc
      iqlst = iqs
      call allocate_mod_runparams
      call allocate_mod_atm_interface
      call allocate_moloch

      call model_zitaf(zita,mo_ztop)
      call model_zitah(zitah,mo_ztop)
      mo_dzita = zita(kz)
      sigma = sigmazita(zita,mo_ztop)
      hsigma = sigmazita(zitah,mo_ztop)
      ak = md_ak(zitah,mo_ztop,mo_h)
      bk = md_bk(zitah,mo_ztop,mo_a0)

      mddom%ht = 0.0_rkx
      mddom%lndcat = 15.0_rkx
      mddom%mask = 0.0_rkx
      mddom%msfx = 1.0_rkx
      mddom%msfu = 1.0_rkx
      mddom%msfv = 1.0_rkx

      dl = raddeg * dx/earthrad

      do concurrent ( j = jde1:jde2, i = ide1:ide2 )
        mddom%xlat(j,i) = clat + dl * (i-real(iy,rkx)*d_half-0.5_rkx)
        mddom%xlon(j,i) = clon + dl * (j-real(jx,rkx)*d_half-0.5_rkx)
        mddom%coriol(j,i) = eomeg2*sin(mddom%xlat(j,i)*degrad)
        mddom%ulat(j,i) = mddom%xlat(j,i)
        mddom%ulon(j,i) = clon + dl * (j-real(jx,rkx)*d_half)
        mddom%vlat(j,i) = clat + dl * (i-real(iy,rkx)*d_half)
        mddom%vlon(j,i) = mddom%xlon(j,i)
      end do
      call exchange_lr(mddom%msfu,1,jde1,jde2,ide1,ide2)
      call exchange_bt(mddom%msfv,1,jde1,jde2,ide1,ide2)
      do concurrent ( j = jdi1:jdi2, i = ice1:ice2 )
        mddom%hx(j,i) = (mddom%ht(j,i) - mddom%ht(j-1,i)) * &
                         mddom%msfu(j,i) * rdx * regrav
      end do
      if ( iproj == 'ROTLLR' ) then
        do concurrent ( j = jce1:jce2, i = idi1:idi2 )
          mddom%hy(j,i) = (mddom%ht(j,i) - mddom%ht(j,i-1)) * &
                           rdx * regrav
        end do
      else
        do concurrent ( j = jce1:jce2, i = idi1:idi2 )
          mddom%hy(j,i) = (mddom%ht(j,i) - mddom%ht(j,i-1)) * &
                           mddom%msfv(j,i) * rdx * regrav
        end do
      end if
      call exchange_lr(mddom%hx,1,jde1,jde2,ice1,ice2)
      call exchange_bt(mddom%hy,1,jce1,jce2,ide1,ide2)
      do concurrent ( j = jce1:jce2, i = ice1:ice2, k = 1:kz )
        mo_atm%zeta(j,i,k) = md_zeta(zitah(k), &
          mddom%ht(j,i),mo_ztop,mo_h,mo_a0)
        mo_atm%fmz(j,i,k) = md_fmz(zitah(k), &
          mddom%ht(j,i),mo_ztop,mo_h,mo_a0)
      end do
      call exchange_lrbt(mo_atm%fmz,1,jce1,jce2,ice1,ice2,1,kz)
      call exchange_lrbt(mo_atm%zeta,2,jce1,jce2,ice1,ice2,1,kz)
      do concurrent ( j = jce1:jce2, i = ice1:ice2, k = 1:kzp1 )
        mo_atm%fmzf(j,i,k) = md_fmz(zita(k), &
                 mddom%ht(j,i),mo_ztop,mo_h,mo_a0)
        mo_atm%zetaf(j,i,k) = md_zeta(zita(k), &
                 mddom%ht(j,i),mo_ztop,mo_h,mo_a0)
      end do
      do concurrent ( j = jce1:jce2 , i = ice1:ice2 , k = 1:kz )
        mo_atm%dz(j,i,k) = mo_atm%zetaf(j,i,k) - mo_atm%zetaf(j,i,k+1)
      end do
    end subroutine param

    subroutine init
      implicit none
      call fill_atmosphere
      call init_moloch
      return
    end subroutine init

    subroutine output(istep)
      implicit none
      integer(ik4) , intent(in) :: istep
      if ( istep == 1 ) then
        call grid_nc_create('ps',.false.,sfs%ps,psout)
        call grid_nc_create('t',.false.,mo_atm%t,tout)
        call grid_nc_write(psout)
        call grid_nc_write(tout)
      end if
      if ( mod(istep, iout) == 0 ) then
        call grid_nc_write(psout)
        call grid_nc_write(tout)
      end if
    end subroutine output

    subroutine finalize
      implicit none
      call grid_nc_destroy(psout)
      call grid_nc_destroy(tout)
      if ( myid == iocpu ) then
        call cpu_time(finish_time)
        call date_and_time(zone=czone,values=tval)
        write(cdata,'(i0.4,"-",i0.2,"-",i0.2," ",i0.2,":",i0.2,":",i0.2,a)') &
              tval(1), tval(2), tval(3), tval(5), tval(6), tval(7), czone
        write (error_unit,*) ': this run stops at  : ', trim(cdata)
        write (error_unit,*) ': Run has been completed using ', &
                nproc, ' processors.'
        write (error_unit,*) ': Total elapsed seconds of run : ', &
                (finish_time - start_time)
      end if
    end subroutine

end program regcm_dynamical_core
