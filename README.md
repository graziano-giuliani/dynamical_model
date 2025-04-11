# RegCM-MOLOCH dynamical model

Requires MPI and NetCDF libraries. Edit the Makefile to select one of the
Fortran compilers (Intel ifx, GNU gfortran, Nvidia nvfortran).

To compile, just type make.

Change model settings in mod\_dynparams.F90:

    integer(ik4) , public , parameter :: nstep = 2000
    integer(ik4) , public , parameter :: iprint = 20
    integer(ik4) , public , parameter :: iout = 100
    integer(ik4) , public , parameter :: iy = 100
    integer(ik4) , public , parameter :: jx = 100
    integer(ik4) , public , parameter :: kz = 30

 * nstep is the total number of timesteps (longer simulation)
 * iprint is the printout frequency on screen (in timesteps units)
 * iout is the output netCDF frequency (in timesteps units)
 * iy , jx are the horizontal dimensions lenght
 * kz is the vertical dimensions lenght

