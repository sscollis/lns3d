#=============================================================================
#
#  Makefile for pre (SGI)
#
#  Author:  Scott Collis
#
#  Revised: 4-5-95
#
#=============================================================================
NAME   = pre.exe
DEBUG  =  
FFLAGS = -g -r8 -c $(DEBUG)
OFLAGS = -g -r8 $(DEBUG) -o $(NAME)
LIB    = 
COMP   = f90
#
#  Define Fortran90 suffix
#
.SUFFIXES: .f90 
#
#  These are standalone objects
#
OBJS = wdata.o  error.o 
#
#  Objects that depend on modules
#
MODS = const.o stencil.o
OBJECTS = pre.o grad.o grad2.o 

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
