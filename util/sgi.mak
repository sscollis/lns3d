#==============================================================================
#  Makefile for LNS utility programs (SGI)
#
#  Author:  Scott Collis
#
#  Revised: 10-17-95
#==============================================================================
DEBUG  = -g
FFLAGS = -r8 -c $(DEBUG)
OFLAGS = -r8 $(DEBUG) -o $(NAME)
LIB    = /home/collis/lib/bslib/bslib.a
ARPACK = /home/collis/disk4/ARPACK/libarpack_SGI.a
COMP   = f90
#
# Define the Fortran 90 suffix
#
.SUFFIXES: .f90 
#==============================================================================
GRAD = grad.o  grad2.o

GRAD_JI = grad_ji.o  grad2_ji.o

conv.sgi: const.o conv.sgi.o wgrid.o wdata.o $(GRAD)
	$(COMP) $(FFLAGS) conv.sgi.o const.o wgrid.o wdata.o -o conv.sgi

lpost: const.o lpost.o  $(GRAD_JI)
	$(COMP) $(FFLAGS) lpost.o const.o $(GRAD_JI) -o lpost

spost: const.o spost.o $(GRAD_JI) filter.o 
	$(COMP) spost.o $(GRAD_JI) filter.o const.o $(LIB) -o spost

npost: const.o npost.o $(GRAD) filter.o error.o
	$(COMP) $(FFLAGS) npost.o $(GRAD) \
        filter.o const.o error.o $(LIB) -o npost

npost_ji: const.o npost_ji.o $(GRAD_JI) filter.o
	$(COMP) npost_ji.o $(GRAD_JI) filter.o const.o $(LIB) -o npost_ji

lpost3d: const.o fmax.o cfilter.o rtsafe.o lpost3d.o 
	$(COMP) lpost3d.o cfilter.o fmax.o rtsafe.o const.o $(LIB) -o lpost3d

subwave: const.o subwave.o 
	$(COMP) subwave.o const.o -o subwave

csubwave: const.o csubwave.o  
	$(COMP) $(FFLAGS) csubwave.o const.o -o csubwave

mkamp: const.o mkamp.o  
	$(COMP) mkamp.o const.o -o mkamp 

mkini: const.o mkini.o  
	$(COMP) mkini.o const.o -o mkini

mkdist: const.o mkdist.o 
	$(COMP) $(FFLAGS) mkdist.o const.o -o mkdist

mkdist3d: const.o mkdist3d.o 
	$(COMP) $(FFLAGS) mkdist3d.o const.o -o mkdist3d

mkmean: const.o mkmean.o spline.o  
	$(COMP) $(FFLAGS) mkmean.o spline.o const.o -o mkmean

inter: const.o inter.o
	$(COMP) inter.o const.o $(LIB) -o inter

genmesh: const.o genmesh.o 
	$(COMP) $(FFLAGS) genmesh.o const.o -o genmesh

stat: const.o stat.o rtsafe.o 
	$(COMP) stat.o const.o rtsafe.o $(LIB) -o stat 

initial: const.o initial.o 
	$(COMP) initial.o const.o -o initial

nconvert: const.o nconvert.o
	$(COMP) nconvert.o const.o -o nconvert

mkvortex: const.o mkvortex.o 
	$(COMP) $(FFLAGS) mkvortex.o const.o /usr/people/collis/lib/slatec.a \
        -lcomplib.sgimath -o mkvortex

mkvortex_v2: const.o mkvortex_v2.o 
	$(COMP) $(FFLAGS) mkvortex_v2.o const.o \
        /usr/people/collis/lib/slatec.a \
        -lcomplib.sgimath -o mkvortex_v2

getevec: getevec.o
	$(COMP) $(FFLAGS) getevec.o $(ARPACK) -o getevec

$(GRAD): stencil.o
	$(COMP) $(FFLAGS) $*.f90

$(GRAD_JI): stencil.o
	$(COMP) $(FFLAGS) $*.f90

all: conv lpost spost npost npost2 lpost3d subwave csubwave mkamp \
     mkini mkdist mkdist3d mkmean genmesh initial nconvert 

clean:
	/bin/rm *.o *.mod

.f90.o:
	$(COMP) $(FFLAGS) $*.f90 

.f.o:
	f77 -col120 $(FFLAGS) $*.f

.c.o:
	cc -n32 $(DEBUG) -c $*.c
