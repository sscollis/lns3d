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
DEBUG    = -O2 -OPT:fold_arith_limit=2000:const_copy_limit=10126:fprop_limit=1500
FFLAGS   = -64 -r8 -col120 -c $(DEBUG)
F90FLAGS = -64 -r8 -c $(DEBUG)
OFLAGS   = -64 -r8 $(DEBUG) -o $(NAME)
LIB      = 
COMP     = f90
#
#  Define Fortran 90 suffix
#
.SUFFIXES: .f90 
#
#  All modules are listed here
#
MODS = \
stuff.o    diff.o  local.o   buff_stuff.o  spg_mod.o  \
local3d.o  pot.o   bump.o
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
fstat.o        getmat.o       lns3d.o        rk.o \
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
cpenta1bc.o    cpenta2bc.o    \
lhs1f3d.o      lhsbt1f3d.o    lhsbt2f3d.o    lhsbc1f3d.o \
lhsbc2f3d.o    lhsbs1f3d.o    expdrv3d.o     rk3d.o \
rkbc3d.o       mlocal3d.o     lhs2f3d.o      fstat3d.o \
potential.o    genb_l.o       genb_n.o	     genbump.o \
lwallbc.o      smoother3d.o   misc.o
#
#  These objects depend on only on stuff.f90
#
OBJS1 = \
getmat.o     gengrid.o   genmean.o    resid.o    \
resstat.o    genini.o    fstat.o                 \
gendist.o    gradbc.o    sstat.o      checkbc.o  \
damper.o     traces.o    genini3d.o   cgradbc.o  \
resstat3d.o  rkbc3d.o    fstat3d.o
#
#  These objects depend on diff.f90
#
OBJS2 = grad.o  grad2.o  g1.o  g2.o  cgrad.o  cgrad2.o
#
#  These objects depend on stuff.f90 and local.f90
#
OBJS3 = \
gena.o     genb.o      genc.o     gend.o     genvij.o    \
rhs.o      mlocal.o    expdrv.o   rk.o       lrhs.o      \
impdrv.o   rkbc.o     imprk.o    rhs_l.o     \
nrk.o      gend_l.o    gend_n.o   gena_l.o   gena_n.o    \
impdrv3d.o expdrv3d.o  genb_l.o   genb_n.o   spg_it.o
#
#  These objects depend on stuff.f90 and diff.f90 and pot.f90
#
OBJS4 = \
wallbc.o  itrbc.o  rhsbc.o  reimann.o  potential.o  lwallbc.o
#
#  These objects depend on stuff.f90, diff.f90, and buff_stuff and pot
#
OBJS5 = \
lhsbc1f.o    lhsbc1.o     lhsbc2.o    lhsbt1.o      lhsbt2.o     lhsbs1.o    \
lhsbc2f.o    lhsbt1f.o    lhsbt2f.o   lhsbs1f.o     buffer.o     smoother.o  \
lhs.o        lhsf.o       lhs1f3d.o   lhsbt1f3d.o   lhsbt2f3d.o  lhsbc1f3d.o \
lhsbc2f3d.o  lhsbs1f3d.o  lhs2f3d.o   smoother3d.o
#
#  These objects depend on stuff.f90, local.f90, and local3d.f90
#
OBJS6 = \
rk3d.o   mlocal3d.o   lrhs3d.o

$(NAME): $(MODS) $(OBJS)
	$(COMP) $(OFLAGS) $(MODS) $(OBJS) $(LIB)

$(OBJS1): stuff.o
	$(COMP) $(F90FLAGS) $*.f90

$(OBJS2): diff.o
	$(COMP) $(F90FLAGS) $*.f90

$(OBJS3): stuff.o local.o
	$(COMP) $(F90FLAGS) $*.f90

$(OBJS4): stuff.o diff.o pot.o
	$(COMP) $(F90FLAGS) $*.f90

$(OBJS5): stuff.o diff.o buff_stuff.o pot.o
	$(COMP) $(F90FLAGS) $*.f90

$(OBJS6): stuff.o local.o local3d.o
	$(COMP) $(F90FLAGS) $*.f90

lns3d.o: stuff.o local.o
	$(COMP) $(F90FLAGS) -DSSD $*.f90

sponge.o: stuff.o spg_mod.o
	$(COMP) $(F90FLAGS) $*.f90

input.o: stuff.o spg_mod.o buff_stuff.o
	$(COMP) $(F90FLAGS) $*.f90

genmtrx.o: stuff.o local.o pot.o
	$(COMP) $(F90FLAGS) $*.f90

genbump.o: stuff.o bump.o
	$(COMP) $(F90FLAGS) $*.f90

itrbc3d.o  rhsbc3d.o: stuff.o diff.o pot.o bump.o
	$(COMP) $(F90FLAGS) $*.f90

dtcfl.o: stuff.o local.o
	$(COMP) $(F90FLAGS) $*.f90

$(MODS):
	$(COMP) $(F90FLAGS) $*.f90

clean:
	/bin/rm *.o *.kmo

.f90.o:
	$(COMP) $(F90FLAGS) $*.f90 

.f.o:
	f90 $(FFLAGS) $*.f

.c.o:
	cc -n32 -O -c $*.c
