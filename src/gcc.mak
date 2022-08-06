#=============================================================================
#
#  Makefile for lns3d (GCC)
#
#  Author:  Scott Collis
#
#  Revised: 12-19-00
#
#=============================================================================
NAME     = lns3d
DEBUG    = -O2 -fopenmp -g
#
# Activates unstead inflow BC's and sponge for linearied 3d runs
#
DEFINES = 
DEFINES += -DUSE_TRANSIENT_EIGENFUNCTION
DEFINES += $(ADDONS)
F77FLAGS = -cpp -fdefault-real-8 -fdefault-double-8 -ffixed-line-length-120 \
-std=legacy $(DEFINES) -c $(DEBUG)
F90FLAGS = -cpp -fdefault-real-8 -fdefault-double-8 $(DEFINES) -c $(DEBUG)
OFLAGS   = $(DEBUG) -o $(NAME)
LIB      = -L$(HOME)/local/OpenBLAS/lib -lopenblas
FC  = gfortran
F77 = gfortran
CC  = gcc
#
#  Optionally use ARpack
#
ifdef USE_ARPACK
  DEFINES += -DUSE_ARPACK
  ifeq ($(ARPACK_DIR),)
    ARPACK_DIR = /usr/local/lib
  endif
  LIB += -L$(ARPACK_DIR) -larpack
endif
#
#  Define Fortran 90 suffix
#
.SUFFIXES: .f90 
#
#  All modules are listed here
#
MODS = global.o stencil.o local2d.o buffer.o sponge.o local3d.o pot.o	\
bump.o stencil.o

MODULES = $(MODS:.o=.mod)
#
#  All objects are listed here
#
OBJS = lhs.o penta1p.o lhsbc1.o penta2bc.o checkbc.o gendist.o		\
lhsbc1f.o penta2p.o damper.o gengrid.o lhsbc2.o resid.o genini.o	\
lhsbc2f.o resstat.o lhsbt1.o temporal.o genmean.o lhsbs1f.o error.o	\
genmtrx.o lhsbt1f.o expdrv.o lhsbt2f.o bcfix.o calcp.o filter.o lhsf.o	\
rhsbc.o fstat2d.o getmat.o lns3d.o g1.o getver.o lhsbt2.o g2.o grad.o	\
lrhs.o smoother.o grad2.o rk.o rkbc.o gradbc.o nrk.o sstat.o		\
nsolve2p.o impdrv.o nsolve1bc.o lhsbs1.o nsolve1p.o traces.o input.o	\
nsolve2bc.o itrbc.o penta1bc.o wdata.o reimann.o genini3d.o	\
cgrad.o impdrv3d.o resstat3d.o cgrad2.o cgradbc.o itrbc3d.o lrhs3d.o	\
rhsbc3d.o cpenta1p.o cpenta2p.o cpenta1bc.o cpenta2bc.o init.o		\
lhs1f3d.o lhsbt1f3d.o lhsbt2f3d.o lhsbc1f3d.o lhsbc2f3d.o lhsbs1f3d.o	\
expdrv3d.o rk3d.o rkbc3d.o lhs2f3d.o fstat3d.o potential.o genbump.o	\
smoother3d.o misc.o rhs_l.o rpenta1p.o rpenta2p.o rpenta1bc.o rpenta2bc.o \
lwallbc.o wallbc.o gena_n.o genb_n.o gena_l.o genb_l.o

ifdef USE_ARPACK
  OBJS += ZnaupdClass.o ZneupdClass.o si_eigdrv3d_new.o eigbc3d.o
endif

#adi.o
#la_eigdrv3d.o

$(NAME): $(MODS) $(OBJS)
	$(FC) $(OFLAGS) $(MODS) $(OBJS) $(LIB)
	\cp $(NAME) $(NAME).exe

help:
	@echo $(MODULES)

getmat.o gengrid.o genmean.o resid.o genini.o gendist.o gradbc.o	\
sstat.o checkbc.o damper.o traces.o genini3d.o cgradbc.o rkbc3d.o	\
spg_it.o buffer.o: global.o

grad.o  grad2.o  g1.o  g2.o  cgrad.o  cgrad2.o: stencil.o

wallbc.o itrbc.o rhsbc.o reimann.o potential.o lwallbc.o: global.o	\
stencil.o pot.o

lhsbc1f.o lhsbc1.o lhsbc2.o lhsbt1.o lhsbt2.o lhsbs1.o lhsbc2f.o	\
lhsbt1f.o lhsbt2f.o lhsbs1f.o smoother.o lhs.o lhsf.o lhs1f3d.o		\
lhsbt1f3d.o lhsbt2f3d.o lhsbc1f3d.o lhsbc2f3d.o lhsbs1f3d.o lhs2f3d.o	\
smoother3d.o: global.o stencil.o buffer.o pot.o

local2d.o local3d.o: global.o

nrk.o expdrv.o impdrv.o lrhs.o rk.o rkbc.o: global.o local2d.o

expdrv3d.o impdrv3d.o lrhs3d.o rk3d.o: global.o local2d.o

lns3d.o: global.o local2d.o local3d.o

sponge.o: global.o

input.o: global.o sponge.o buffer.o

genmtrx.o: global.o local2d.o local3d.o pot.o

genbump.o: global.o bump.o

itrbc3d.o  rhsbc3d.o: global.o stencil.o bump.o

temporal.o: global.o local2d.o

init.o: global.o buffer.o sponge.o

resstat.o resstat3d.o: global.o

fstat2d.o fstat3d.o: global.o

calcp.o: global.o pot.o

rhs_l.o: global.o local2d.o

ZneupdClass.o: ZnaupdClass.o

eigdrv3d.o: ZneupdClass.o ZnaupdClass.o global.o

la_eigdrv3d.o: global.o local2d.o local3d.o

si_eigdrv3d_new.o: global.o ZneupdClass.o ZnaupdClass.o local2d.o local3d.o

eigbc3d.o: global.o stencil.o bump.o

adi.o: global.o local2d.o

.f90.mod:
	$(FC) $(F90FLAGS) $*.f90

clean:
	$(RM) *.o *.mod

.f90.o:
	$(FC) $(F90FLAGS) $*.f90 

.f.o:
	$(F77) $(F77FLAGS) $*.f

.c.o:
	$(CC) -O -c $*.c
