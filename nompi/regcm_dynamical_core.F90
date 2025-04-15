program regcm_dynamical_core
  use, intrinsic :: iso_fortran_env
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
  integer(ik4) :: ipunit
  integer(ik4) :: i
  integer(ik4) :: ierr
  integer(ik4) , dimension(8) :: tval
  real(rkx) :: start_time , finish_time , last_time
  character (len=32) :: cdata='?'
  character (len=5) :: czone='?'
  character (len=256) :: namelistfile
  type(grid_nc_var2d) :: psout
  type(grid_nc_var3d) :: tout
  character(len=*) , parameter :: f99001 = &
        '(2x," GIT Revision: ",a," compiled at: data : ",a,"  time: ",a,/)'

  call cpu_time(start_time)
  last_time = start_time
  write (output_unit,"(/,2x,'This is RegCM development version')")
  write (output_unit,f99001)  '1.0', __DATE__ , __TIME__

  call get_command_argument(1,value=namelistfile)
  call readinput(namelistfile)

  call memory_init
  call initparam
  call set_nproc

  call setup_model_indexes

  call param
  call init

  print *, 'Starting Integration!'

  do i = 1 , nstep
    call output(i)
    lprint = mod(i,iprint) == 0
    call moloch(i,lprint)
  end do

  call memory_destroy
  call finalize

  contains

    subroutine readinput(namelistfile)
      implicit none
      character(len=*) , intent(in) :: namelistfile
      namelist /dimensions/ jx , iy , kz
      namelist /timeparams/ nstep , iprint , iout

      open(newunit=ipunit, file=namelistfile, status='old', &
           action='read', iostat=ierr)
      if ( ierr /= 0 ) then
        write(output_unit, *) 'Namelist file error, using defaults.'
      else
        rewind(ipunit)
        read (ipunit, nml=dimensions, iostat=ierr, err=100)
        if ( ierr /= 0 ) then
           write(error_unit, *) 'Error reading namelist dimensions'
           stop
        end if
        rewind(ipunit)
        read (ipunit, nml=timeparams, iostat=ierr, err=200)
        if ( ierr /= 0 ) then
          write(error_unit, *) 'Error reading namelist timeparams'
          stop
        end if
      end if

      return

 100  write(error_unit, *) 'Error reading namelist dimensions'
      stop
 200  write(error_unit, *) 'Error reading namelist timeparams'
      stop
    end subroutine readinput

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
      mddom%msfu(jde1ga,ide1:ide2) = mddom%msfu(jde2,ide1:ide2)
      mddom%msfu(jde2ga,ide1:ide2) = mddom%msfu(jde1,ide1:ide2)
      mddom%msfv(jde1:jde2,ide1ga) = mddom%msfv(jde1:jde2,ide2)
      mddom%msfv(jde1:jde2,ide2ga) = mddom%msfv(jde1:jde2,ide1)
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
      mddom%hx(jde1ga,ice1:ice2) = mddom%hx(jde2,ice1:ice2)
      mddom%hx(jde2ga,ice1:ice2) = mddom%hx(jde1,ice1:ice2)
      mddom%hy(jce1:jce2,ide1ga) = mddom%hy(jce1:jce2,ide2)
      mddom%hy(jce1:jce2,ide2ga) = mddom%hy(jce1:jce2,ide1)
      do concurrent ( j = jce1:jce2, i = ice1:ice2, k = 1:kz )
        mo_atm%zeta(j,i,k) = md_zeta(zitah(k), &
          mddom%ht(j,i),mo_ztop,mo_h,mo_a0)
        mo_atm%fmz(j,i,k) = md_fmz(zitah(k), &
          mddom%ht(j,i),mo_ztop,mo_h,mo_a0)
      end do
      mo_atm%fmz(jce1ga,ice1:ice2,:) = mo_atm%fmz(jce2,ice1:ice2,:)
      mo_atm%fmz(jce2ga,ice1:ice2,:) = mo_atm%fmz(jce1,ice1:ice2,:)
      mo_atm%fmz(jce1:jce2,ice1ga,:) = mo_atm%fmz(jce1:jce2,ice2,:)
      mo_atm%fmz(jce1:jce2,ice2ga,:) = mo_atm%fmz(jce1:jce2,ice1,:)
      mo_atm%zeta(jce1ga,ice1:ice2,:) = mo_atm%zeta(jce2,ice1:ice2,:)
      mo_atm%zeta(jce2ga,ice1:ice2,:) = mo_atm%zeta(jce1,ice1:ice2,:)
      mo_atm%zeta(jce1gb,ice1:ice2,:) = mo_atm%zeta(jce2-1,ice1:ice2,:)
      mo_atm%zeta(jce2gb,ice1:ice2,:) = mo_atm%zeta(jce1+1,ice1:ice2,:)
      mo_atm%zeta(jce1:jce2,ice1ga,:) = mo_atm%zeta(jce1:jce2,ice2,:)
      mo_atm%zeta(jce1:jce2,ice2ga,:) = mo_atm%zeta(jce1:jce2,ice1,:)
      mo_atm%zeta(jce1:jce2,ice1gb,:) = mo_atm%zeta(jce1:jce2,ice2-1,:)
      mo_atm%zeta(jce1:jce2,ice2gb,:) = mo_atm%zeta(jce1:jce2,ice1+1,:)
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
      call cpu_time(finish_time)
      call date_and_time(zone=czone,values=tval)
      write(cdata,'(i0.4,"-",i0.2,"-",i0.2," ",i0.2,":",i0.2,":",i0.2,a)') &
            tval(1), tval(2), tval(3), tval(5), tval(6), tval(7), czone
      write (error_unit,*) ': this run stops at  : ', trim(cdata)
      write (error_unit,*) ': Total elapsed seconds of run : ', &
              (finish_time - start_time)
    end subroutine

end program regcm_dynamical_core
