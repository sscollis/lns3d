#==============================================================================
#  Makefile for LNS utility programs (GCC)
#
#  Author:  Scott Collis
#
#  Revised: 4/30/2022 
#==============================================================================
DEBUG  = -O2 -fopenmp
FFLAGS = -cpp -fdefault-real-8 -freal-4-real-8 $(DEBUG)
OFLAGS = $(DEBUG) -o $(NAME)
LIB    = -L$(HOME)/local/OpenBLAS/lib -lopenblas
ARPACK = -L/usr/local/lib -larpack
SLATEC = -L$(HOME)/local/slatec/lib -lslatec
FC     = gfortran
F77    = gfortran
CC     = gcc-11
#
# Define the Fortran 90 suffix
#
.SUFFIXES: 
.SUFFIXES: .f90 .f .c .o

#==============================================================================

GRAD = grad.o  grad2.o

GRAD_JI = grad_ji.o  grad2_ji.o

BSLIB = bslib1.o bslib2.o

ifdef USE_NR
  NRLIB = rtsafe.o
endif

ALL = conv.sgi lpost subwave csubwave mkamp mkini mkdist mkdist3d mkmean \
genmesh initial nconvert getevec mkvortex ij2ji ji2ij mkmean_ji mkdist3d_ji \
mkdist_ji r4tor8 dirp3d p3dlns3d unipot

ifdef USE_NR
  ALL += spost npost lpost3d lpost3d_ji stat
endif

all: $(ALL) 

conv.sgi: const.o conv.sgi.o wgrid.o wdata.o $(GRAD)
	$(FC) $(FFLAGS) conv.sgi.o const.o wgrid.o wdata.o -o conv.sgi

lpost: const.o lpost.o  $(GRAD_JI)
	$(FC) $(FFLAGS) lpost.o const.o $(GRAD_JI) -o lpost

spost: const.o spost.o $(GRAD_JI) filter.o $(BSLIB) $(NRLIB)
	$(FC) spost.o $(GRAD_JI) filter.o const.o $(LIB) $(BSLIB) $(NRLIB) \
	-o spost

npost: const.o npost.o $(GRAD) filter.o error.o $(BSLIB) $(NRLIB)
	$(FC) $(FFLAGS) npost.o $(GRAD) $(BSLIB) $(NRLIB) \
	filter.o const.o error.o $(LIB) -o npost

npost_ji: const.o npost_ji.o $(GRAD_JI) filter.o
	$(FC) npost_ji.o $(GRAD_JI) filter.o const.o $(LIB) -o npost_ji

lpost3d: const.o fmax.o cfilter.o $(NRLIB) lpost3d.o $(BSLIB)
	$(FC) lpost3d.o cfilter.o fmax.o $(NRLIB) const.o $(LIB) $(BSLIB) \
	-o lpost3d

lpost3d_ji: const.o fmax.o cfilter.o $(NRLIB) lpost3d_ji.o $(BSLIB)
	$(FC) lpost3d_ji.o cfilter.o fmax.o $(NRLIB) const.o $(LIB) $(BSLIB) \
	-o lpost3d_ji

subwave: const.o subwave.o 
	$(FC) subwave.o const.o -o subwave

csubwave: const.o csubwave.o  
	$(FC) $(FFLAGS) csubwave.o const.o -o csubwave

mkamp: const.o mkamp.o  
	$(FC) mkamp.o const.o -o mkamp 

mkini: const.o mkini.o  
	$(FC) mkini.o const.o -o mkini

mkdist: const.o mkdist.o 
	$(FC) $(FFLAGS) mkdist.o const.o -o mkdist

mkdist_ji: const.o mkdist_ji.o 
	$(FC) $(FFLAGS) mkdist_ji.o const.o -o mkdist_ji

mkdist3d: const.o mkdist3d.o 
	$(FC) $(FFLAGS) mkdist3d.o const.o -o mkdist3d

mkdist3d_ji: const.o mkdist3d_ji.o 
	$(FC) $(FFLAGS) mkdist3d_ji.o const.o -o mkdist3d_ji

mkmean: const.o mkmean.o spline.o  
	$(FC) $(FFLAGS) mkmean.o spline.o const.o -o mkmean

mkmean_ji: const.o mkmean_ji.o spline.o  
	$(FC) $(FFLAGS) mkmean_ji.o spline.o const.o -o mkmean_ji

inter: const.o inter.o $(BSLIB)
	$(FC) inter.o const.o $(LIB) $(BSLIB) -o inter

genmesh: const.o genmesh.o 
	$(FC) $(FFLAGS) genmesh.o const.o -o genmesh

stat: const.o stat.o $(NRLIB) $(BSLIB) 
	$(FC) stat.o const.o $(NRLIB) $(LIB) $(BSLIB) -o stat 

initial: const.o initial.o 
	$(FC) initial.o const.o -o initial

nconvert: const.o nconvert.o
	$(FC) nconvert.o const.o -o nconvert

mkvortex: const.o mkvortex.o 
	$(FC) $(FFLAGS) mkvortex.o const.o $(SLATEC) $(LIB) -o mkvortex

mkvortex_v2: const.o mkvortex_v2.o 
	$(FC) $(FFLAGS) mkvortex_v2.o const.o $(SLATEC) $(LIB) -o mkvortex_v2

getevec: getevec.o
	$(FC) $(FFLAGS) getevec.o $(ARPACK) -o getevec

$(GRAD): stencil.o

$(GRAD_JI): stencil.o

clean:
	$(RM) *.o *.mod $(ALL)

.f90.o:
	$(FC) $(FFLAGS) -c $*.f90 

.f90:
	$(FC) $(FFLAGS) -o $* $*.f90 

.f.o:
	$(F77) $(FFLAGS) -std=legacy -ffixed-line-length-120 -c $*.f

.f:
	$(F77) $(FFLAGS) -std=legacy -ffixed-line-length-120 -o $* $*.f

.c.o:
	$(CC) $(DEBUG) -c $*.c

.o:
	$(FC) -o $* $*.o
