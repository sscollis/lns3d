#==============================================================================
#  Makefile for mesh utilites (SGI)
#
#  Author:  Scott Collis
#
#  Revised: 6-15-97
#==============================================================================
DEBUG  = -g
FFLAGS = -r8 $(DEBUG) -c -extend_source
OFLAGS = -r8 $(DEBUG)
COMP   = f90
LIB    = /usr/people/collis/lib/bslib.a
#
#  Define Fortran 90 suffix
#
.SUFFIXES: .f90 
#
#==============================================================================
#
ogrid: ogrid.o
	f90 $(OFLAGS) -o ogrid ogrid.o

confpc: wgrid.o  wdata.o  confpc.o
	$(COMP) $(OFLAGS) -o confpc wgrid.o  wdata.o  confpc.o 

interpc: interpc.o calcd.o grad.o
	f90 $(OFLAGS) -o interpc interpc.o calcd.o grad.o $(LIB)

interpc.o: interpc.f90
	f90 $(FFLAGS) interpc.f90 -c

cinterpc: cinterpc.o calcd.o grad.o
	f90 $(OFLAGS) -o cinterpc cinterpc.o calcd.o grad.o $(LIB) 

cinterpc.o: cinterpc.f90
	f90 $(FFLAGS) cinterpc.f90 -c

grad.o: grad.f90 stencil.o
	f90 $(FFLAGS) grad.f90 stencil.o -c

stencil.o: stencil.f90
	f90 $(FFLAGS) stencil.f90 -c

pc: pc.o wdata.o wgrid.o
	$(COMP) $(OFLAGS) -o pc pc.o wdata.o wgrid.o

level: level.o splcrv.o crvdist.o wgrid.o wdata.o elliptic.o conformal.o
	$(COMP) $(OFLAGS) -o level level.o splcrv.o crvdist.o wgrid.o wdata.o \
        elliptic.o conformal.o $(LIB) 

cyl: cyl.o
	$(COMP) $(OFLAGS) -o cyl cyl.o

circ: circ.o
	$(COMP) $(OFLAGS) -o circ circ.o

clean:
	/bin/rm *.o *.mod

all: interpc cinterpc grad.o stencil.o pc

.c.o:
	cc $(DEBUG) -c $*.c

.f.o:
	f90 $(FFLAGS) $*.f

.f90.o:
	f90 $(FFLAGS) $*.f90 
