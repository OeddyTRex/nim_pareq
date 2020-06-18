#This is the nimfl makefile.

#The dependencies listed here assume that nimlib
#and nimcore modules (in the nimlib and nimcore directories)
#are already current.

NIMHOME = ..
COREDIR = $(NIMHOME)/nimcore
NIMLIBDIR = $(NIMHOME)/nimlib

include $(NIMHOME)/make.inc

BINFILE = nimfl

LIBS =  $(COREDIR)/libnimpar.a \
	$(ODE_LIB) \
	$(NIMLIBDIR)/libnimlib.a \
	$(NIMLIBDIR)/libnimmpi.a \
	$(COREDIR)/libnimbltype.a \
	$(BLAS_HOME) $(BLAS_LIB)

MODDIR = -I$(NIMLIBDIR) -I$(COREDIR)

OBJS = \
        input.o \
        parbcast_input.o \
        input0.o \
        global.o \
        dump.o \
        dump_reset.o \
        onefl.o \
        int_segment.o \
        plane_cross.o \
        flder.o \
        get_bfield.o \
	magnetic_axis.o\
	rfft.o\
	start_positions.o\
	psigrid.o\
	poincare.o\
	threshold.o\
	polfft.o\
	nimfl.o

CORE_OBJS = \
	$(COREDIR)/pardata.o \
	$(COREDIR)/time.o \
	$(COREDIR)/fields.o \
	$(COREDIR)/seam_storage.o

$(BINFILE): $(OBJS)
	$(LDR) -o $(BINFILE) $(LDR_FLAGS) $(CORE_OBJS) $(OBJS) $(LIBS)

clean: cleanlinks
	-rm *.o *.mod *.d

realclean: clean cleanlinks
	-rm $(BINFILE)

.f.o: 
	$(F90) -c $(FFLAGS) $(MODDIR) $*.f

# Module dependencies

dump.o: input.o
magnetic_axis.o: input.o global.o dump.o dump_reset.o
start_positions.o: input0.o global.o dump.o \
	dump_reset.o magnetic_axis.o
poincare.o: input.o global.o dump.o dump_reset.o \
	start_positions.o input0.o
psigrid.o: input.o global.o dump.o dump_reset.o \
	input0.o magnetic_axis.o start_positions.o
polfft.o: input.o global.o dump.o dump_reset.o \
	input0.o magnetic_axis.o start_positions.o rfft.o
threshold.o: input.o global.o dump.o dump_reset.o \
	magnetic_axis.o start_positions.o input0.o
nimfl.o: input.o global.o dump.o \
	dump_reset.o magnetic_axis.o start_positions.o input0.o
onefl.o: input.o dump_reset.o
plane_cross.o: input.o input0.o
int_segment.o: plane_cross.o dump_reset.o input0.o
dump_reset.o: input.o dump.o
get_bfield.o: input.o global.o dump_reset.o

# links

input.f:
	ln -f -s $(NIMHOME)/nimset/input.f input.f
parbcast_input.f:
	ln -f -s $(NIMHOME)/nimrod/parbcast_input.f parbcast_input.f
global.f:
	ln -f -s $(NIMHOME)/nimrod/global.f global.f
dump.f:
	ln -f -s $(NIMHOME)/nimset/dump.f dump.f
dump_reset.f:
	ln -f -s $(NIMHOME)/nimset/dump_reset.f dump_reset.f

cleanlinks:
	-rm global.f dump.f dump_reset.f input.f parbcast_input.f
