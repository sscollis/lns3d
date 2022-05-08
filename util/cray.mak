#==============================================================================
#  Makefile for LNS utility programs (Cray)
#
#  Author:  Scott Collis
#
#  Revised: 10-17-95
#==============================================================================
NPROC  = 1
DEBUG  =
FFLAGS = -c $(DEBUG)
OFLAGS = $(DEBUG) -o $(NAME)
LIB    = /u/va/collis/lib/bslib.a
COMP   = f90
#
# Define the Fortran 90 suffix
#
.SUFFIXES: .f90 
#==============================================================================
GRAD = grad.o  grad2.o

lpost: const.o lpost.o $(GRAD)
	$(COMP) lpost.o $(GRAD) const.o  -o lpost

spost: const.o spost.o $(GRAD) filter.o
#	$(COMP) spost.o $(GRAD) filter.o const.o -limsl -o spost
	$(COMP) spost.o $(GRAD) filter.o const.o $(LIB) -o spost

npost: const.o npost.o $(GRAD) filter.o
	$(COMP) npost.o $(GRAD) filter.o const.o $(LIB) -o npost
#	$(COMP) npost.o $(GRAD) filter.o const.o -limsl -o npost

npost2: const.o npost2.o $(GRAD) filter.o
	$(COMP) npost2.o $(GRAD) filter.o const.o $(LIB) -o npost2
#	$(COMP) npost2.o $(GRAD) filter.o const.o -limsl -o npost2

lpost3d: const.o fmax.o nr_rtsafe.o lpost3d.o cfilter.o 
#	$(COMP) lpost3d.o cfilter.o fmax.o nr_rtsafe.o const.o -limsl -o lpost3d
	$(COMP) lpost3d.o cfilter.o fmax.o nr_rtsafe.o const.o $(LIB) -o lpost3d

subwave: const.o subwave.o
	$(COMP) subwave.o const.o -o subwave

csubwave: const.o csubwave.o
	$(COMP) csubwave.o const.o -o csubwave

mkamp: const.o mkamp.o
	$(COMP) mkamp.o const.o -o mkamp 

mkini: const.o mkini.o
	$(COMP) mkini.o const.o -o mkini

mkdist: const.o mkdist.o
	$(COMP) mkdist.o const.o -o mkdist

mkmean: const.o mkmean.o spline.o
#	$(COMP) mkmean.o spline.o const.o -limsl -o mkmean
	$(COMP) mkmean.o spline.o const.o $(LIB) -o mkmean

inter: const.o inter.o
#	$(COMP) inter.o const.o -limsl -o inter
	$(COMP) inter.o const.o $(LIB) -o inter

genmesh: const.o genmesh.o
	$(COMP) genmesh.o const.o -o genmesh

stat: const.o stat.o nr_rtsafe.o
#	$(COMP) stat.o const.o nr_rtsafe.o -limsl -o stat 
	$(COMP) stat.o const.o nr_rtsafe.o $(LIB) -o stat 

initial: const.o initial.o
	$(COMP) initial.o const.o -o initial

mkdist3d: mkdist3d.o
	$(COMP) mkdist3d.o -o mkdist3d

nconvert: const.o nconvert.o
	$(COMP) nconvert.o const.o -o nconvert

$(GRAD): stencil.o
	$(COMP) $(FFLAGS) $*.f90

all: lpost spost npost npost2 lpost3d subwave csubwave mkamp \
mkini mkdist mkdist3d mkmean inter genmesh stat initial nconvert 

clean:
	/bin/rm *.o

.f90.o: const.o
	$(COMP) $(FFLAGS) $*.f90 

.f.o:
	$(COMP) -N80 $(FFLAGS) $*.f
