#This is the fluxgrid and sol makefile.

NIMHOME = ..

include $(NIMHOME)/make.inc

LIBS = $(NIMHOME)/nimlib/libnimlib.a $(ODE_LIB)
MODDIR=-I$(NIMHOME)/nimlib

OBJS = \
	nim_stop.o \
	physdat.o \
	time_stat.o \
	global.o \
	input.o \
	grid.o \
	direct.o \
	inverse.o \
	analyze.o \
	read_efit.o \
	read_galkin.o \
	read_miller.o \
	read_miller8.o \
	read_millasc.o \
	read_chum.o \
	read_chease.o \
	read_rsteq.o \
	read_soloviev.o \
	output.o \
	codes.o \
	fluxgrid.o

OBJS_SOL = \
	input_sol.o \
	nim_stop.o \
	sol.o

all: fluxgrid sol

fluxgrid: $(OBJS)
	$(LDR) -o fluxgrid $(LDR_FLAGS) $(OBJS) $(LIBS)

sol: $(OBJS_SOL)
	$(LDR) -o sol $(LDR_FLAGS) $(OBJS_SOL) $(LIBS)

.f.o: 
	$(F90) -c $(FFLAGS) $(MODDIR) $*.f

read_miller.o:read_miller.f
	$(F90) -c $(FFLAGS_NODBL) $(MODDIR) $*.f

clean:
	-rm *.o *.mod *.d

realclean: clean cleanlinks
	-rm sol fluxgrid 

# module dependencies

grid.o: input.o global.o
direct.o: input.o global.o grid.o
fluxgrid.o: input.o global.o output.o analyze.o
inverse.o: input.o global.o 
output.o: global.o input.o physdat.o analyze.o
read_efit.o: direct.o input.o global.o physdat.o
read_galkin.o: inverse.o input.o global.o physdat.o
read_miller.o: inverse.o input.o global.o physdat.o
read_miller8.o: inverse.o input.o global.o physdat.o
read_millasc.o: inverse.o input.o global.o physdat.o
read_chum.o: inverse.o input.o global.o physdat.o
read_soloviev.o: direct.o input.o global.o physdat.o
read_rsteq.o: direct.o input.o global.o physdat.o
read_chease.o: direct.o input.o global.o physdat.o
sol.o: input_sol.o 
analyze.o : global.o input.o physdat.o

# links

physdat.f:
	ln -f -s $(NIMHOME)/nimset/physdat.f physdat.f

cleanlinks:
	rm physdat.f

