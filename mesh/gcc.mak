#==============================================================================
#  Makefile for mesh utilites (GCC)
#
#  Author:  Scott Collis
#
#  Revised: 2-9-2020 
#==============================================================================
DEBUG  = -O2
FFLAGS   = -cpp -freal-4-real-8 -fdefault-real-8 -std=legacy -ffixed-line-length-120 \
					 -Wno-align-commons $(DEFINES) $(DEBUG)
F90FLAGS = -cpp -freal-4-real-8 -fdefault-real-8 -Wno-align-commons $(DEFINES) $(DEBUG)
OFLAGS = $(DEBUG)
COMP   = gfortran 
F77    = gfortran
CC     = gcc-11
LIB    = 
#
#  Define Fortran 90 suffix
#
.SUFFIXES: 
.SUFFIXES: .f90 .f .c .o
#
#==============================================================================
#
BSLIB = bslib1.o bslib2.o
#
# Numerical Recipese routines cannot be shared since they are commercial
# licenced you will need a valid license for these to use.
#
ifdef USE_NR
  NRLIB = nr_odeint.o
  NRSPLINE = nr_spline.o
endif

MODS = stencil.o

ALL = ogrid cyl circ reverse cyl_org 

ifdef USE_NR
  ALL += mse msecurve confpc interpc cinterpc pc testpc pcurve level 
endif

all: $(ALL)

ogrid: ogrid.o
	$(COMP) $(OFLAGS) -o ogrid ogrid.o

mse: mse.o $(NRLIB) $(NRSPLINE)
	$(COMP) $(OFLAGS) $(NRLIB) $(NRSPLINE) -o mse mse.o

msecurve: msecurve.o $(NRLIB)
	$(COMP) $(OFLAGS) $(NRLIB) -o msecurve msecurve.o

confpc: wgrid.o  wdata.o  confpc.o $(NRLIB)
	$(COMP) $(OFLAGS) $(NRLIB) -o confpc wgrid.o wdata.o confpc.o 

interpc: interpc.o calcd.o grad.o $(BSLIB) $(NRLIB) $(NRSPLINE)
	$(COMP) $(OFLAGS) -o interpc interpc.o calcd.o grad.o $(NRLIB) $(NRSPLINE) $(BSLIB) $(LIB)

cinterpc: cinterpc.o calcd.o grad.o $(BSLIB) $(NRLIB) $(NRSPLINE) 
	$(COMP) $(OFLAGS) -o cinterpc cinterpc.o calcd.o grad.o $(NRLIB) $(NRSPLINE) $(BSLIB) $(LIB) 

pc: pc.o wdata.o wgrid.o $(NRLIB)
	$(COMP) $(OFLAGS) $(NRLIB) -o pc pc.o wdata.o wgrid.o

testpc: testpc.o $(NRLIB)
	$(COMP) $(OFLAGS) $(NRLIB) -o testpc testpc.o

pcurve: pcurve.o $(NRLIB)
	$(COMP) $(OFLAGS) $(NRLIB) -o pcurve pcurve.o

level: level.o splcrv.o crvdist.o wgrid.o wdata.o elliptic.o conformal.o $(BSLIB) $(NRLIB)
	$(COMP) $(OFLAGS) -o level level.o splcrv.o crvdist.o wgrid.o wdata.o \
		elliptic.o conformal.o $(LIB) $(BSLIB) $(NRLIB) 

grad.o: stencil.o

cyl: cyl.o
	$(COMP) $(OFLAGS) -o cyl cyl.o

circ: circ.o
	$(COMP) $(OFLAGS) -o circ circ.o

$(MODS):

clean:
	$(RM) *.o *.mod $(ALL)

.c.o:
	$(CC) $(DEBUG) -c $*.c

.f.o:
	$(F77) $(FFLAGS) -c $*.f

.f:
	$(F77) $(FFLAGS) -o $* $*.f

.f90.o:
	$(COMP) $(F90FLAGS) -c $*.f90 
