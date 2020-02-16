#==============================================================================
#  Makefile for mesh utilites
#
#  Author:   Scott Collis
#
#  Revised:  6-15-97
#==============================================================================
NPROC  = 2
DEBUG  =
FFLAGS = -N80 $(DEBUG)
OFLAGS = $(DEBUG) -o $(NAME)
COMP   = cf77 
#==============================================================================

confpc: wgrid.o  wdata.o  confpc.o
	$(COMP) $(DEBUG) wgrid.o  wdata.o  confpc.o -limsl -o confpc

interpc: interpc.o calcd.o grad.o
	f90 $(DEBUG) interpc.o calcd.o grad.o -limsl -o interpc 

interpc.o: interpc.f90
	f90 $(DEBUG) interpc.f90 -c

cinterpc: cinterpc.o calcd.o grad.o
	f90 $(DEBUG) cinterpc.o calcd.o grad.o -limsl -o cinterpc 

cinterpc.o: cinterpc.f90
	f90 $(DEBUG) cinterpc.f90 -c

grad.o: grad.f90 stencil.o
	f90 grad.f90 -p stencil.o -c

stencil.o: stencil.f90
	f90 stencil.f90 -c

pc: pc.o wdata.o wgrid.o
	$(COMP) $(DEBUG) pc.o wdata.o wgrid.o -o pc

level: level.o splcrv.o crvdist.o wgrid.o wdata.o elliptic.o conformal.o
	$(COMP) $(DEBUG) level.o splcrv.o crvdist.o wgrid.o wdata.o elliptic.o conformal.o -limsl -o level 

clean:
	/bin/rm *.o

all: interpc cinterpc grad.o stencil.o pc level
