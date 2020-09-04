UNAME_S := $(shell uname -s)
UNAME_N := $(shell uname -n | sed 's/[0-9]*//g')
ARCH := $(shell uname -p)


ifeq ($(UNAME_S),Linux)
ifeq ($(ARCH),ppc64le)
LIB =
F2PY = f2py3
OPTFLAG = -O
FFLAGS= -fopenmp -cpp
CMD = ${F2PY} --fcompiler=gfortran -c  --link-lapack_opt --opt="-cpp -heap-arrays ${OPTFLAG}" ${LIB} --f90flags="${FFLAGS}"
F90=gfortran
F77=${F90}
F77FLAGS = ${FFLAGS}
else
#LIB = -lmkl_rt 
LIB = -lmkl_rt -L/opt/intel-14.0/mkl/lib/intel64
# if --link-lapack_opt fails, try to link to mkl manually with mkl_rt library
F2PY = f2py3
OPTFLAG = -fast -xHost
FFLAGS= -DIntelMKL -qopenmp -fpp
CMD = ${F2PY} --fcompiler=intelem   -c  --link-lapack_opt --opt="-cpp -heap-arrays ${OPTFLAG}" ${LIB} --f90flags="${FFLAGS}"
F90=ifort
F77=${F90}
F77FLAGS = -qopenmp -O3 -warn nousage -xHost -no-prec-div -mtune=pentium4 -noautomatic -heap-arrays 0 -fpp
endif
endif

ifeq ($(UNAME_S),Darwin)
LIB = -L${HOME}/lib
F2PY = f2py
OPTFLAG = -O
FFLAGS = -ffree-line-length-none -std=f2008 -fopenmp -cpp
CMD = ${F2PY} --fcompiler=gfortran -c  --link-lapack_opt --opt="-cpp ${OPTFLAG}" ${LIB} --f90flags="${FFLAGS}"
F90 = gfortran
# WARNING: DO NOT use -O* options otherwise code behavior is crazy, seemingly due to the associate statement of fortran 2003
F77=${F90}
F77FLAGS =-fopenmp -O -cpp
endif
#FFLAGS+= -g -fbacktrace -fbounds-check -ffree-line-length-none -ftrapv

.SUFFIXES:

#BINARIES = bregman.so f_phonon.so f_util.so bcs_driver.so
BINARIES = bregman.so f_phonon.so f_util.so

BCSSOLVER= cssolve/num_types.f90  cssolve/matrix_sets.f90  cssolve/laplace.f90  cssolve/bcs.f90
BCSSOLVEROBJ=$(BCSSOLVER:.f90=.o)


TET = csld/phonon/dostet.f90 csld/phonon/pdstet.f90 csld/phonon/setk06.f90
TET_OBJ=$(TET:.f90=.o)
DMdipole = csld/phonon/DM_dipole_dipole.f csld/phonon/gbox.f csld/phonon/tripl.f csld/phonon/derfc.f csld/phonon/find_Ewald_eta_screened.f csld/phonon/DMstandard.f csld/phonon/lattc.f csld/phonon/cross.f csld/phonon/Fewald.f
DMdipole_OBJ=$(DMdipole:.f=.o)


.PHONY: all
all : $(BINARIES)

bregman.so: cssolve/bregman.f90

f_phonon.so: csld/phonon/f_phonon.f90 ${TET_OBJ} ${DMdipole_OBJ}

f_util.so: compile/f_util/f_util.f90

clean:
	rm -f cssolve/*.o cssolve/*.mod cssolve/*.a *.o *.mod $(BINARIES) csld/phonon/*.o


%.o: %.f90
	$(F90) -fPIC -c $(FFLAGS) $^ -o $@

%.o: %.f
	$(F77) -fPIC -c $(F77FLAGS) $^ -o $@

cssolve/bcs.a: $(BCSSOLVEROBJ)
	ar ru $@ $?
	ranlib $@

cssolve/bcs.so: $(BCSSOLVEROBJ)
	$(F90) -shared -fPIC $(FFLAGS) -o $@ $^


$(BINARIES):
	rm -f $@
	${CMD} -m $(basename $@) $^
	if [ ! -f $@ ]; then mv -f $(basename $@).*.so $@; fi


bcs_driver.so: cssolve/bcs_driver.f90 $(BCSSOLVEROBJ)
	rm -f $@
	echo "dict(real=dict(dp='double'))" > .f2py_f2cmap
	${CMD} -m $(basename $@) $^
	if [ ! -f $@ ]; then mv -f $(basename $@).*.so $@; fi
	rm .f2py_f2cmap
