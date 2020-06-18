#This is the xlog and nimplot makefile.

#The dependencies listed here assume that nimlib, nimmpi,
#and nimcore modules (in the nimlib and nimcore directories)
#are already current.

NIMHOME = ..
COREDIR = $(NIMHOME)/nimcore
NIMLIBDIR = $(NIMHOME)/nimlib

include $(NIMHOME)/make.inc

LIBS = $(COREDIR)/libnimpar.a $(NIMLIBDIR)/libnimmpi.a \
	$(COREDIR)/libnimblock.a $(COREDIR)/libnimbltype.a \
	$(COREDIR)/libnimiter.a $(COREDIR)/libnimmatrix.a \
	$(NIMLIBDIR)/libnimlib.a \
	$(NIMHOME)/externals/liblapack_dummy.a \
	$(NIMHOME)/externals/libsslu_dummy.a \
	$(NIMHOME)/externals/libslud_dummy.a
	
MODDIR = -I$(NIMLIBDIR) -I$(COREDIR)

OBJSXL = \
	nim_stop.o \
	xlog.o

OBJS = \
	physdat.o \
	global.o \
	nim_stop.o \
	input.o \
	matrix_storage.o \
	dump.o \
	nimplot_ints.o \
	nimplot_init.o \
	plot_data.o \
	data_dealloc.o \
	nimplot_mgt.o \
	diagnose.o \
	contour.o \
	plot_io.o \
	nimplot.o

# The following "core" objects are not in a library.

CORE_OBJS = \
        $(COREDIR)/seam_storage.o \
        $(COREDIR)/fields.o \
        $(COREDIR)/computation_pointers.o \
        $(COREDIR)/pardata.o \
        $(COREDIR)/time.o \
        $(COREDIR)/edge.o \
        $(COREDIR)/boundary.o \
        $(COREDIR)/regularity.o \
        $(COREDIR)/block_create_tang.o \
        $(COREDIR)/generic_evals.o \
        $(COREDIR)/fft_mod.o \
        $(COREDIR)/iter_externals.o

all: nimplot xlog

xlog:	$(OBJSXL)
	$(LDR) -o xlog $(LDR_FLAGS) $(OBJSXL) $(LIBS)

nimplot: $(OBJS)
	$(LDR) -o nimplot $(LDR_FLAGS) $(CORE_OBJS) $(OBJS) $(LIBS)

clean:
	-rm *.o *.mod *.d

realclean: clean cleanlinks
	-rm nimplot xlog

.f.o: 
	$(F90) -c $(FFLAGS) $(MODDIR) $*.f

xlog.o:xlog.f
	$(F90) -c $(FFLAGS_NODBL) $(MODDIR) $*.f

# object dependencies

dump.o: input.o
nimplot_ints.o: global.o input.o physdat.o
nimplot_init.o: input.o matrix_storage.o nimplot_ints.o
plot_data.o: input.o dump.o global.o physdat.o
data_dealloc.o: input.o global.o
nimplot_mgt.o: input.o physdat.o nimplot_ints.o \
	global.o contour.o plot_data.o matrix_storage.o
diagnose.o: input.o
plot_io.o: plot_data.o
nimplot.o: dump.o diagnose.o contour.o plot_data.o \
	global.o nimplot_mgt.o nimplot_init.o physdat.o

# links

physdat.f:
	ln -f -s $(NIMHOME)/nimset/physdat.f .
global.f:
	ln -f -s $(NIMHOME)/nimrod/global.f .
input.f:
	ln -f -s $(NIMHOME)/nimset/input.f .
dump.f:
	ln -f -s $(NIMHOME)/nimset/dump.f .

cleanlinks:
	rm input.f dump.f physdat.f global.f
