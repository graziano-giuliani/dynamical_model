
FC = gfortran
#FC = ifx
#FC = nvfortran

#FCFLAGS = `nf-config --fflags` -O0 -g -Wall -pedantic -fcheck=all -fbacktrace
FCFLAGS = `nf-config --fflags` -O3
LDFLAGS = `nf-config --flibs`

.SUFFIXES: .F90 .f90 .o

%.o: %.F90
	$(FC) $(FCFLAGS) -c $<

%.o: %.f90
	$(FC) $(FCFLAGS) -c $<

OBJS = mod_intkinds.o mod_realkinds.o mod_space.o mod_memutil.o \
       mod_constants.o mod_regcm_types.o mod_dynparam.o mod_runparams.o \
       mod_atmosphere.o mod_zita.o mod_moloch.o mod_mppparam.o

all :: regcm_dynamical_core

regcm_dynamical_core : regcm_dynamical_core.F90 $(OBJS)
	$(FC) $(FCFLAGS) -o $@ $< $(OBJS) $(LDFLAGS)

clean :
	rm -f *.o *.mod regcm_dynamical_core

mod_intkinds.o : mod_intkinds.F90
mod_realkinds.o : mod_realkinds.F90
mod_space.o : mod_space.F90 mod_realkinds.o mod_intkinds.o
mod_memutil.o : mod_memutil.F90 mod_space.o
mod_constants.o : mod_constants.F90 mod_realkinds.o
mod_regcm_types.o : mod_regcm_types.F90 mod_realkinds.o
mod_dynparam.o : mod_dynparam.F90 mod_realkinds.o mod_intkinds.o
mod_runparams.o : mod_runparams.F90 mod_dynparam.o mod_memutil.o mod_constants.o
mod_atmosphere.o : mod_atmosphere.F90 mod_runparams.o mod_regcm_types.o
mod_zita.o : mod_zita.F90 mod_dynparam.o mod_constants.o
mod_moloch.o : mod_moloch.F90 mod_atmosphere.o mod_zita.o mod_mppparam.o
mod_mppparam.o : mod_mppparam.F90 mod_runparams.o mod_regcm_types.o
