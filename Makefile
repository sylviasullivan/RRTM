# Pulling the template partially from
# https://www.lmd.jussieu.fr/~lguez/templates-of-makefiles-for-fortran-programs-and-libraries.html
# and partially from
# https://github.com/ZedThree/fort_depend.py
#
# 1. Source file
# Assume that the source files are all the f90 files in this directory.
VPATH = .

# 2. Objects and executable files.
srtm := $(shell cat ${VPATH}/srtm_kgb_list)

# Must be the name of the file containing the main program, without
# suffix (else linking will not work).
execut = rrtm_driver

# 3. Compiler-dependent part
FC = gfortran
# Adding in low optimization for debugging.
FFLAGS = -ggdb -Wall -Wmaybe-uninitialized -finit-local-zero -ffpe-trap=invalid,zero -fcheck=all -I/sw/rhel6-x64/netcdf/netcdf_fortran-4.4.3-gcc71/include -fbacktrace
LDLIBS = -L/sw/rhel6-x64/netcdf/netcdf_fortran-4.4.3-gcc71/lib -lnetcdff
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/sw/rhel6-x64/netcdf/netcdf_fortran-4.4.3-gcc71/lib

# 4. Rules
# Use the Fortran compiler for linking.
# First rule
all: ${execut}

${execut}:
	$(FC) $(FFLAGS) $^ -o $@ $(LDLIBS) $(srtm)

%.o: %.f90 
	$(FC) $(FFLAGS) -c $< $(LDLIBS) 

include ${VPATH}/rrtm.dep

clean:
	rm -f ${execut} *.o *.mod

