#==============================================================================
#  Makefile for mesh utilites (GCC)
#
#  Author:  Scott Collis
#
#  Revised: 2-9-2020 
#==============================================================================
DEBUG  = -O2 
#-ffpe-trap=invalid,zero,overflow 
DEFINES += -DUSE_BSLIB
FFLAGS = -cpp -fdefault-real-8 -fdefault-double-8 -std=legacy \
-ffixed-line-length-120 -Wno-align-commons $(DEFINES) $(DEBUG)
F90FLAGS = -cpp -fdefault-real-8 -fdefault-double-8 -Wno-align-commons \
$(DEFINES) $(DEBUG)
OFLAGS = $(DEBUG)
FC     = gfortran 
F77    = gfortran
CC     = gcc-12
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
  ifeq ($(LIBNR_DIR),)
    LIBNR_DIR = $(HOME)/git/NR-utilities
  endif
  LIB += -L$(LIBNR_DIR) -lnr 
endif

MODS = stencil.o

ALL = ogrid cyl circ reverse cyl_org 

ifdef USE_NR
  ALL += mse msecurve confpc interpc cinterpc pc testpc pcurve level 
endif

all: $(ALL)

ogrid: ogrid.o
	$(FC) $(OFLAGS) -o ogrid ogrid.o

mse: mse.o
	$(FC) mse.o $(OFLAGS) $(LIB) -o mse

msecurve: msecurve.o
	$(FC) msecurve.o $(OFLAGS) $(LIB) -o msecurve

confpc: wgrid.o  wdata.o  confpc.o
	$(FC) $(OFLAGS) -o confpc wgrid.o wdata.o confpc.o $(LIB) 

interpc: interpc.o calcd.o grad.o $(BSLIB)
	$(FC) $(OFLAGS) -o interpc interpc.o calcd.o grad.o $(BSLIB) $(LIB)

cinterpc: cinterpc.o calcd.o grad.o $(BSLIB)
	$(FC) $(OFLAGS) -o cinterpc cinterpc.o calcd.o grad.o $(BSLIB) $(LIB) 

pc: pc.o wdata.o wgrid.o
	$(FC) $(OFLAGS) -o pc pc.o wdata.o wgrid.o $(LIB)

testpc: testpc.o
	$(FC) $(OFLAGS) -o testpc testpc.o $(LIB)

pcurve: pcurve.o
	$(FC) $(OFLAGS) -o pcurve pcurve.o $(LIB)

level: level.o splcrv.o crvdist.o wgrid.o wdata.o elliptic.o conformal.o $(BSLIB)
	$(FC) $(OFLAGS) -o level level.o splcrv.o crvdist.o wgrid.o wdata.o \
		elliptic.o conformal.o $(LIB) $(BSLIB)

grad.o: stencil.o

cyl: cyl.o
	$(FC) $(OFLAGS) -o cyl cyl.o

circ: circ.o
	$(FC) $(OFLAGS) -o circ circ.o

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
	$(FC) $(F90FLAGS) -c $*.f90 
