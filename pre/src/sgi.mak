#=============================================================================
#
#  Makefile for pre (SGI)
#
#  Author:  Scott Collis
#
#  Revised: 4-5-95
#
#=============================================================================
NAME   = pre.sgi
DEBUG  =  
FFLAGS = -O -r8 -c $(DEBUG)
OFLAGS = -O -r8 $(DEBUG) -o $(NAME)
LIB    = 
COMP   = f90
#
#  Define Fortran90 suffix
#
.SUFFIXES: .f90 
#
#  These are standalone objects
#
OBJS = grad6.o  wdata.o  error.o 
#
#  Objects that depend on modules
#
MODS = const.o
OBJECTS = pre.o

$(NAME): $(OBJS) $(MODS) $(OBJECTS)
	 $(COMP) $(OFLAGS) $(MODS) $(OBJS) $(OBJECTS) $(LIB)

$(OBJECTS): $(MODS)
	$(COMP) $(FFLAGS) $*.f90

.f90.o:
	$(COMP) $(FFLAGS) $*.f90 

.f.o:
	f77 -col120 $(FFLAGS) $*.f

clean:
	/bin/rm *.o *.mod
