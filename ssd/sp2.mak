#==============================================================
#  Makefile for lns3d on the SP2
#
#  Author:  Scott Collis
#
#  Revised: 7-10-96
#==============================================================
NAME   = lns3d 
DEBUG  = 
FFLAGS = -O2 -Q -ghot -qalias=noaryovrlp -qarch=pwr2 \
         -qtune=pwr2 -qrealsize=8 -c $(DEBUG)
OFLAGS = $(DEBUG) -o $(NAME)
LIB    = /a/raid1fr2sw/u14/collis/lib/misc.o
COMP   = xlf90
#
#  Define Fortran 90 suffix
#
.SUFFIXES: .mod .f90 
#
#  All modules are listed here
#
MODS = stuff.o  diff.o  local.o   buff_stuff.o  spg_mod.o  \
local3d.o  pot.o  bump.o
#
#  All objects are listed here
#
OBJS = \
gend.o         lhs.o          penta1p.o      gend_l.o \
buffer.o       gend_p.o       lhsbc1.o       penta2bc.o \
checkbc.o      gendist.o      lhsbc1f.o      penta2p.o \
damper.o       gengrid.o      lhsbc2.o       resid.o \
genini.o       lhsbc2f.o      resstat.o      lhsbt1.o \
dtcfl.o        genmean.o      lhsbs1f.o      rhs.o \
error.o        genmtrx.o      lhsbt1f.o      rhs_l.o \
expdrv.o       genmtrx_p.o    lhsbt2f.o      rhs_p.o \
filter.o       genvij.o       lhsf.o         rhsbc.o \
fstat.o        getmat.o       lns3d.o        rk.o \
g1.o           getver.o       rkbc.o         lhsbt2.o \
g2.o           grad.o         lrhs.o         smoother.o \
gena.o         grad2.o        mlocal.o       sponge.o \
gena_p.o       gradbc.o       nrk.o          sstat.o \
genb.o         impdrv.o       nsolve1bc.o    lhsbs1.o \
genb_p.o       imprk.o        nsolve1p.o     traces.o \
genc.o         input.o        nsolve2bc.o    wallbc.o \
genc_p.o       itrbc.o        penta1bc.o     wdata.o \
gend_n.o       gena_l.o       gena_n.o       reimann.o \
genini3d.o     cgrad.o        impdrv3d.o     resstat3d.o \
cgrad2.o       cgradbc.o      itrbc3d.o      lrhs3d.o \
rhsbc3d.o      cpenta1p.o     cpenta2p.o     csolve1p.o \
cpenta1bc.o    cpenta2bc.o    csolve1bc.o    csolve2bc.o \
lhs1f3d.o      lhsbt1f3d.o    lhsbt2f3d.o    lhsbc1f3d.o \
lhsbc2f3d.o    lhsbs1f3d.o    expdrv3d.o     rk3d.o \
rkbc3d.o       mlocal3d.o     lhs2f3d.o      fstat3d.o \
potential.o    genb_l.o       genb_n.o	     genbump.o \
lwallbc.o      smoother3d.o
#
# These objects depend on only on stuff.f90
#
OBJS1 = \
getmat.o    gengrid.o   genmean.o    resid.o     \
resstat.o   genini.o    fstat.o      gendist.o   \
gradbc.o    sstat.o     checkbc.o    damper.o    \
traces.o    genini3d.o  cgradbc.o    resstat3d.o \
rkbc3d.o    fstat3d.o
#
# These objects depend on diff.f90
#
OBJS2 = grad.o  grad2.o  g1.o  g2.o  cgrad.o  cgrad2.o
#
# These objects depend on stuff.f90 and local.f90
#
OBJS3 = \
gena.o     genb.o     genc.o      gend.o   \
genvij.o   rhs.o      mlocal.o    expdrv.o    rk.o        \
lrhs.o     impdrv.o   lns3d.o     rkbc.o      imprk.o     \
gena_p.o   genb_p.o   genc_p.o    gend_p.o    genmtrx_p.o \
rhs_p.o    rhs_l.o    nrk.o       gend_l.o    gend_n.o    \
gena_l.o   gena_n.o   impdrv3d.o  expdrv3d.o  genb_l.o    \
genb_n.o
#
# These objects depend on stuff.f90 and diff.f90 and pot.f90
#
OBJS4 = wallbc.o  itrbc.o  rhsbc.o  reimann.o  potential.o  lwallbc.o
#
# These objects depend on stuff.f90, diff.f90, 
# buff_stuff, and pot
#
OBJS5 = \
lhsbc1f.o    lhsbc1.o     lhsbc2.o     lhsbt1.o     \
lhsbt2.o     lhsbs1.o     lhsbc2f.o    lhsbt1f.o    \
lhsbt2f.o    lhsbs1f.o    buffer.o     smoother.o   \
lhs.o        lhsf.o       lhs1f3d.o    lhsbt1f3d.o  \
lhsbt2f3d.o  lhsbc1f3d.o  lhsbc2f3d.o  lhsbs1f3d.o  \
lhs2f3d.o    smoother3d.o
#
# These objects depend on stuff.f90, local.f90, and local3d.f90
#
OBJS6 = rk3d.o   mlocal3d.o   lrhs3d.o

$(NAME): $(MODS) $(OBJS)
	$(COMP) $(OFLAGS) $(OBJS) $(LIB)

$(OBJS1): stuff.mod
	$(COMP) $(FFLAGS) $*.f

$(OBJS2): diff.mod
	$(COMP) $(FFLAGS) $*.f

$(OBJS3): stuff.mod local.mod
	$(COMP) $(FFLAGS) $*.f

$(OBJS4): stuff.mod diff.mod pot.mod
	$(COMP) $(FFLAGS) $*.f

$(OBJS5): stuff.mod diff.mod buff_stuff.mod pot.mod
	$(COMP) $(FFLAGS) $*.f

$(OBJS6): stuff.mod local.mod local3d.mod
	$(COMP) $(FFLAGS) $*.f

sponge.o: stuff.mod spg_mod.mod
	$(COMP) $(FFLAGS) $*.f

input.o: stuff.mod spg_mod.mod buff_stuff.mod
	$(COMP) $(FFLAGS) $*.f

genmtrx.o: stuff.mod local.mod pot.mod
	$(COMP) $(FFLAGS) $*.f

genbump.o: stuff.mod bump.mod
	$(COMP) $(FFLAGS) $*.f

itrbc3d.o  rhsbc3d.o: stuff.mod diff.mod pot.mod bump.mod
	$(COMP) $(FFLAGS) $*.f

dtcfl.o: stuff.mod local.mod
	$(COMP) $(FFLAGS) $*.f

$(MODS):
	$(COMP) $(FFLAGS) $*.f

clean:
	/bin/rm *.o

.f90.f:
	mv $*.f90 $*.f

.f.o:
	$(COMP) $(FFLAGS) $*.f

.f.mod:
	$(COMP) $(FFLAGS) $*.f
