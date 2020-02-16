#
#  Makefile for pre
#
#  Author:  Scott Collis
#
#  Revised: 4-5-95
#
NPROC  =  4
NAME   =  pre 
DEBUG  = 
FFLAGS =  -c $(DEBUG)
OFLAGS =  $(DEBUG) -o $(NAME)
LIB    = -limsl
COMP   =  f90
#
#  Define Fortran90 suffix
#
.SUFFIXES: .f90 
#
#  These are standalone objects
#
OBJS = \
grad6.o  wdata.o  pre.o  error.o

$(NAME): $(OBJS) 
	 $(COMP) $(OFLAGS) $(OBJS) $(LIB)

stuff:
	$(COMP) $(FFLAGS) stuff.f90

clean:
	/bin/rm *.o

pre.o:const.o
	$(COMP) -p const.o $(FFLAGS) $*.f90

.f90.o:
	 $(COMP) $(FFLAGS) $*.f90 

.f.o:
	 $(COMP) $(FFLAGS) $*.f

	
