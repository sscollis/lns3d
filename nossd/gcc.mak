#=============================================================================
#
#  Makefile for lns3d (SGI)
#
#  Author:  Scott Collis
#
#  Revised: 1-29-97
#
#=============================================================================
NAME     = lns3d 
DEBUG    = -O2 -fopenmp 
FFLAGS   = -cpp -fdefault-real-8 -std=legacy -ffixed-line-length-120 -c $(DEBUG)
F90FLAGS = -cpp -fdefault-real-8 -c $(DEBUG)
OFLAGS   = -fdefault-real-8 $(DEBUG) -o $(NAME)
LIB      = -L$(HOME)/local/OpenBLAS/lib -lopenblas
COMP     = gfortran
F77      = gfortran
CC       = gcc-9
#
#  Define Fortran 90 suffix
#
.SUFFIXES: .f90 
#
#  All modules are listed here
#
MODS = \
global.o   stencil.o  local.o   buff_mod.o  spg_mod.o  \
local3d.o  pot.o      bump.o
#
#  All objects are listed here
#
#  csolve1p.o  csolve1bc.o    csolve2bc.o    ! these aren't used?
#
OBJS = \
gend.o         lhs.o          penta1p.o      gend_l.o \
buffer.o       spg_it.o       lhsbc1.o       penta2bc.o \
checkbc.o      gendist.o      lhsbc1f.o      penta2p.o \
damper.o       gengrid.o      lhsbc2.o       resid.o \
genini.o       lhsbc2f.o      resstat.o      lhsbt1.o \
dtcfl.o        genmean.o      lhsbs1f.o      rhs.o \
error.o        genmtrx.o      lhsbt1f.o      rhs_l.o \
expdrv.o       lhsbt2f.o      bcfix.o \
filter.o       genvij.o       lhsf.o         rhsbc.o \
fstat2d.o      getmat.o       lns3d.o        rk.o \
g1.o           getver.o       rkbc.o         lhsbt2.o \
g2.o           grad.o         lrhs.o         smoother.o \
gena.o         grad2.o        mlocal.o       sponge.o \
gradbc.o       nrk.o          sstat.o \
genb.o         impdrv.o       nsolve1bc.o    lhsbs1.o \
imprk.o        nsolve1p.o     traces.o \
genc.o         input.o        nsolve2bc.o    wallbc.o \
itrbc.o        penta1bc.o     wdata.o \
gend_n.o       gena_l.o       gena_n.o       reimann.o \
genini3d.o     cgrad.o        impdrv3d.o     resstat3d.o \
cgrad2.o       cgradbc.o      itrbc3d.o      lrhs3d.o \
rhsbc3d.o      cpenta1p.o     cpenta2p.o \
cpenta1bc.o    cpenta2bc.o    init.o \
lhs1f3d.o      lhsbt1f3d.o    lhsbt2f3d.o    lhsbc1f3d.o \
lhsbc2f3d.o    lhsbs1f3d.o    expdrv3d.o     rk3d.o \
rkbc3d.o       mlocal3d.o     lhs2f3d.o      fstat3d.o \
potential.o    genb_l.o       genb_n.o	     genbump.o \
lwallbc.o      smoother3d.o   misc.o
#
#  These objects depend on only on global.f90
#
OBJS1 = \
getmat.o     gengrid.o   genmean.o    resid.o    \
genini.o     gendist.o   gradbc.o     sstat.o    \
checkbc.o    damper.o    traces.o     genini3d.o \
cgradbc.o    rkbc3d.o
#
#  These objects depend on stencil.f90
#
OBJS2 = grad.o  grad2.o  g1.o  g2.o  cgrad.o  cgrad2.o
#
#  These objects depend on global.f90 and local.f90
#
OBJS3 = \
gena.o     genb.o      genc.o     gend.o     genvij.o    \
rhs.o      mlocal.o    expdrv.o   rk.o       lrhs.o      \
impdrv.o   rkbc.o      imprk.o    rhs_l.o     \
nrk.o      gend_l.o    gend_n.o   gena_l.o   gena_n.o    \
impdrv3d.o expdrv3d.o  genb_l.o   genb_n.o   spg_it.o
#
#  These objects depend on global.f90 and stencil.f90 and pot.f90
#
OBJS4 = \
wallbc.o  itrbc.o  rhsbc.o  reimann.o  potential.o  lwallbc.o
#
#  These objects depend on global.f90, stencil.f90, and buff_mod and pot
#
OBJS5 = \
lhsbc1f.o    lhsbc1.o     lhsbc2.o    lhsbt1.o      lhsbt2.o     lhsbs1.o    \
lhsbc2f.o    lhsbt1f.o    lhsbt2f.o   lhsbs1f.o     buffer.o     smoother.o  \
lhs.o        lhsf.o       lhs1f3d.o   lhsbt1f3d.o   lhsbt2f3d.o  lhsbc1f3d.o \
lhsbc2f3d.o  lhsbs1f3d.o  lhs2f3d.o   smoother3d.o
#
#  These objects depend on global.f90, local.f90, and local3d.f90
#
OBJS6 = \
rk3d.o   mlocal3d.o   lrhs3d.o

$(NAME): $(MODS) $(OBJS)
	$(COMP) $(OFLAGS) $(MODS) $(OBJS) $(LIB)

$(OBJS1): global.o

$(OBJS2): stencil.o

$(OBJS3): global.o local.o

$(OBJS4): global.o stencil.o pot.o

$(OBJS5): global.o stencil.o buff_mod.o pot.o

$(OBJS6): global.o local.o local3d.o

lns3d.o: global.o local.o

sponge.o: global.o spg_mod.o

input.o: global.o spg_mod.o buff_mod.o

genmtrx.o: global.o local.o pot.o

genbump.o: global.o bump.o

itrbc3d.o  rhsbc3d.o: global.o stencil.o pot.o bump.o

dtcfl.o: global.o local.o

init.o: global.o buff_mod.o spg_it.o 

resstat.o resstat3d.o: global.o

fstat2d.o fstat3d.o: global.o

$(MODS):

clean:
	$(RM) *.o *.mod

.f90.o:
	$(COMP) $(F90FLAGS) $*.f90 

.f.o:
	$(F77) $(FFLAGS) $*.f

.c.o:
	$(CC) -O2 -c $*.c
