#==============================================================================
#  Makefile for LNS utility programs (GCC)
#
#  Author:  Scott Collis
#
#  Revised: 10-17-95
#==============================================================================
DEBUG  = -O2 -fopenmp
FFLAGS = -cpp -fdefault-real-8 -freal-4-real-8 $(DEBUG)
OFLAGS = $(DEBUG) -o $(NAME)
LIB    = -L$(HOME)/local/OpenBLAS/lib -lopenblas
ARPACK = -L/usr/local/lib -larpack
COMP   = gfortran
FC     = gfortran
F77    = gfortran
CC     = gcc-10
#
# Define the Fortran 90 suffix
#
.SUFFIXES: 
.SUFFIXES: .f90 .f .c .o

#==============================================================================

GRAD = grad.o  grad2.o

GRAD_JI = grad_ji.o  grad2_ji.o

BSLIB = bslib1.o bslib2.o

NRLIB = rtsafe.o

ALL = conv.sgi lpost spost npost lpost3d subwave csubwave mkamp \
mkini mkdist mkdist3d mkmean genmesh initial nconvert getevec \
mkvortex ij2ji ji2ij mkmean_ji lpost3d_ji mkdist3d_ji mkdist_ji \
r4tor8 dirp3d stat p3dlns3d unipot

all: $(ALL) 

conv.sgi: const.o conv.sgi.o wgrid.o wdata.o $(GRAD)
	$(COMP) $(FFLAGS) conv.sgi.o const.o wgrid.o wdata.o -o conv.sgi

lpost: const.o lpost.o  $(GRAD_JI)
	$(COMP) $(FFLAGS) lpost.o const.o $(GRAD_JI) -o lpost

spost: const.o spost.o $(GRAD_JI) filter.o $(BSLIB) $(NRLIB)
	$(COMP) spost.o $(GRAD_JI) filter.o const.o $(LIB) $(BSLIB) $(NRLIB) -o spost

npost: const.o npost.o $(GRAD) filter.o error.o $(BSLIB) $(NRLIB)
	$(COMP) $(FFLAGS) npost.o $(GRAD) $(BSLIB) $(NRLIB) \
        filter.o const.o error.o $(LIB) -o npost

npost_ji: const.o npost_ji.o $(GRAD_JI) filter.o
	$(COMP) npost_ji.o $(GRAD_JI) filter.o const.o $(LIB) -o npost_ji

lpost3d: const.o fmax.o cfilter.o rtsafe.o lpost3d.o $(BSLIB)
	$(COMP) lpost3d.o cfilter.o fmax.o rtsafe.o const.o $(LIB) $(BSLIB) -o lpost3d

lpost3d_ji: const.o fmax.o cfilter.o rtsafe.o lpost3d_ji.o $(BSLIB)
	$(COMP) lpost3d_ji.o cfilter.o fmax.o rtsafe.o const.o $(LIB) $(BSLIB) -o lpost3d_ji

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

mkdist_ji: const.o mkdist_ji.o 
	$(COMP) $(FFLAGS) mkdist_ji.o const.o -o mkdist_ji

mkdist3d: const.o mkdist3d.o 
	$(COMP) $(FFLAGS) mkdist3d.o const.o -o mkdist3d

mkdist3d_ji: const.o mkdist3d_ji.o 
	$(COMP) $(FFLAGS) mkdist3d_ji.o const.o -o mkdist3d_ji

mkmean: const.o mkmean.o spline.o  
	$(COMP) $(FFLAGS) mkmean.o spline.o const.o -o mkmean

mkmean_ji: const.o mkmean_ji.o spline.o  
	$(COMP) $(FFLAGS) mkmean_ji.o spline.o const.o -o mkmean_ji

inter: const.o inter.o $(BSLIB)
	$(COMP) inter.o const.o $(LIB) $(BSLIB) -o inter

genmesh: const.o genmesh.o 
	$(COMP) $(FFLAGS) genmesh.o const.o -o genmesh

stat: const.o stat.o rtsafe.o $(BSLIB) 
	$(COMP) stat.o const.o rtsafe.o $(LIB) $(BSLIB) -o stat 

initial: const.o initial.o 
	$(COMP) initial.o const.o -o initial

nconvert: const.o nconvert.o
	$(COMP) nconvert.o const.o -o nconvert

mkvortex: const.o mkvortex.o 
	$(COMP) $(FFLAGS) mkvortex.o const.o \
		$(HOME)/lib/libslatec.a $(LIB) -o mkvortex

mkvortex_v2: const.o mkvortex_v2.o 
	$(COMP) $(FFLAGS) mkvortex_v2.o const.o \
		$(HOME)/lib/libslatec.a $(LIB) -o mkvortex_v2

getevec: getevec.o
	$(COMP) $(FFLAGS) getevec.o $(ARPACK) -o getevec

$(GRAD): stencil.o

$(GRAD_JI): stencil.o

clean:
	$(RM) *.o *.mod $(ALL)

.f90.o:
	$(COMP) $(FFLAGS) -c $*.f90 

.f90:
	$(COMP) $(FFLAGS) -o $* $*.f90 

.o:
	$(COMP) -o $* $*.o

.f:
	$(F77) $(FFLAGS) -std=legacy -ffixed-line-length-120 -o $* $*.f

.f.o:
	$(F77) $(FFLAGS) -std=legacy -ffixed-line-length-120 -c $*.f

.c.o:
	$(CC) $(DEBUG) -c $*.c
