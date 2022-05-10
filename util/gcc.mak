#==============================================================================
#  Makefile for LNS utility programs (GCC)
#
#  Author:  Scott Collis
#
#  Revised: 4/30/2022 
#==============================================================================
DEBUG  = -O2 -fopenmp
FFLAGS = -cpp -fdefault-real-8 -fdefault-double-8 -std=legacy \
         -ffixed-line-length-120 $(DEBUG)
OFLAGS = $(DEBUG) -o $(NAME)
LIB    = -L$(HOME)/local/OpenBLAS/lib -lopenblas
ARPACK = -L/usr/local/lib -larpack
SLATEC = -L../slatec/lib -lslatec
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
#
# Optionally use Numerical-Recipes (commercial licensed)
#
ifdef USE_NR
  ifeq ($(LIBNR_DIR),)
    LIBNR_DIR = $(HOME)/git/NR-utilities
  endif
  LIB += -L$(LIBNR_DIR) -lnr 
endif

MATHLIB = zeroin.o d1mach.o

ALL = conv-sgi lpost subwave csubwave mkamp mkini mkdist mkdist3d mkmean \
genmesh initial nconvert getevec mkvortex ij2ji ji2ij mkmean_ji mkdist3d_ji \
mkdist_ji r4tor8 dirp3d p3dlns3d unipot npost spost lpost3d lpost3d_ji stat

all: $(ALL) 

conv-sgi: const.o conv-sgi.o wgrid.o wdata.o $(GRAD)
	$(FC) $(FFLAGS) conv-sgi.o const.o wgrid.o wdata.o -o conv-sgi

lpost: const.o lpost.o $(GRAD_JI)
	$(FC) $(FFLAGS) lpost.o const.o $(GRAD_JI) -o lpost

spost: const.o spost.o $(GRAD_JI) filter.o $(BSLIB)
	$(FC) spost.o $(GRAD_JI) filter.o const.o $(LIB) $(BSLIB) \
	-o spost

npost: const.o npost.o $(GRAD) filter.o error.o $(BSLIB) $(MATHLIB)
	$(FC) $(FFLAGS) npost.o $(GRAD) $(BSLIB) $(MATHLIB) \
	filter.o const.o error.o $(LIB) -o npost

npost_ji: const.o npost_ji.o $(GRAD_JI) filter.o
	$(FC) npost_ji.o $(GRAD_JI) filter.o const.o $(LIB) -o npost_ji

lpost3d: const.o fmax.o cfilter.o $(MATHLIB) lpost3d.o $(BSLIB)
	$(FC) lpost3d.o cfilter.o fmax.o $(MATHLIB) const.o $(LIB) $(BSLIB) \
	-o lpost3d

lpost3d_ji: const.o fmax.o cfilter.o $(MATHLIB) lpost3d_ji.o $(BSLIB)
	$(FC) lpost3d_ji.o cfilter.o fmax.o $(MATHLIB) const.o $(LIB) $(BSLIB) \
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

stat: const.o stat.o $(BSLIB) 
	$(FC) stat.o const.o $(LIB) $(BSLIB) -o stat 

initial: const.o initial.o 
	$(FC) initial.o const.o -o initial

nconvert: const.o nconvert.o
	$(FC) nconvert.o const.o -o nconvert

mkvortex: const.o mkvortex.o 
	$(FC) $(FFLAGS) mkvortex.o const.o $(SLATEC) -Xlinker -rpath -Xlinker ../slatec/lib \
	$(LIB) -o mkvortex

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
	$(F77) $(FFLAGS) -c $*.f

.f:
	$(F77) $(FFLAGS) -o $* $*.f

.c.o:
	$(CC) $(DEBUG) -c $*.c

.o:
	$(FC) -o $* $*.o
