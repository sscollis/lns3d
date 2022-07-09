#=============================================================================
#
#  Makefile for pre (GCC)
#
#  Author:  Scott Collis
#
#  Revised: 2-12-2020
#
#=============================================================================
NAME     = pre.exe
DEBUG    = -g 
FFLAGS   = -cpp -fdefault-real-8 -O2 -c $(DEBUG)
F77FLAGS = -cpp -fdefault-real-8 -ffixed-line-length-120 -std=legacy \
-O2 -c $(DEBUG)
OFLAGS   = -fdefault-real-8 $(DEBUG) -o $(NAME)
LIB      = 
COMP     = gfortran
F77      = gfortran
CC       = gcc-11
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

all: $(NAME)

$(NAME): $(OBJS) $(MODS) $(OBJECTS)
	 $(COMP) $(OFLAGS) $(MODS) $(OBJS) $(OBJECTS) $(LIB)

$(OBJECTS): $(MODS)
	$(COMP) $(FFLAGS) $*.f90

.f90.o:
	$(COMP) $(FFLAGS) $*.f90 

.f.o:
	$(F77) $(F77FLAGS) $*.f

clean:
	$(RM) *.o *.mod $(NAME)
