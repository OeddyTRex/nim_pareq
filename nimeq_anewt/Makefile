#This is the makefile for the nimeq_anewt version of nimeq.
#Source files that are not modified are linked from the nimeq
#directory.

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
	$(NIMHOME)/nimcore/libnimlocate.a \
	$(NIMHOME)/externals/liblapack_dummy.a \
	$(SLU_NEQ_LIB) $(ODE_LIB) $(BLAS_LIB)

MODDIR = -I$(NIMLIBDIR) -I$(COREDIR)

OBJS = \
	physdat.o \
	global.o \
	input.o \
	parbcast_input.o \
	nim_stop.o \
	nimeq_input.o \
	matrix_storage.o \
	dump.o \
	dump_reset.o \
	nimeq_mod.o \
	nimeq_ints.o \
	nimeq_srf_ints.o \
	nimeq_init.o \
	nimeq_all.o \
	nimeq_data.o \
	data_dealloc.o \
	eq_swap.o \
	nimeq_btr.o \
	nimeq_free.o \
	nimeq_mgt.o \
	diagnose.o \
	contour.o \
	plot_io.o \
	nimeq_dot_mgt.o \
	nimeq.o

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
	$(COREDIR)/surface.o \
	$(COREDIR)/aph_eval.o \
	$(COREDIR)/iter_externals.o

all: nimeq 

nimeq: $(OBJS)
	$(LDR) -o nimeq $(LDR_FLAGS) $(CORE_OBJS) $(OBJS) $(LIBS)

clean:
	-rm *.o *.mod *.d

realclean: clean cleanlinks
	-rm nimeq

.f.o: 
	$(F90) -c $(FFLAGS) $(MODDIR) $*.f

# object dependencies

nim_stop.o: input.o nimeq_input.o
dump.o: input.o
dump_reset.o: input.o dump.o
nimeq_ints.o: global.o physdat.o nimeq_input.o
nimeq_srf_ints.o: global.o physdat.o nimeq_input.o
nimeq_init.o: input.o matrix_storage.o nimeq_ints.o nimeq_input.o
nimeq_data.o: input.o dump.o global.o physdat.o nimeq_all.o
data_dealloc.o: input.o global.o nimeq_all.o nimeq_mod.o
eq_swap.o: physdat.o
nimeq_btr.o: nimeq_input.o nimeq_all.o
nimeq_free.o: input.o nimeq_input.o physdat.o nimeq_mod.o 
nimeq_mgt.o: input.o physdat.o nimeq_ints.o global.o \
	contour.o nimeq_data.o matrix_storage.o nimeq_srf_ints.o \
        nimeq_input.o nimeq_mod.o nimeq_btr.o nimeq_free.o 
nimeq_dot_mgt.o: input.o physdat.o nimeq_ints.o global.o \
	nimeq_data.o nimeq_input.o nimeq_mod.o nimeq_free.o
diagnose.o: input.o
plot_io.o: nimeq_data.o
nimeq.o: dump.o dump_reset.o diagnose.o contour.o nimeq_data.o \
	global.o nimeq_mgt.o nimeq_init.o \
	physdat.o nimeq_input.o nimeq_mod.o \
	nimeq_all.o nimeq_dot_mgt.o

# links

physdat.f:
	ln -f -s $(NIMHOME)/nimset/physdat.f .
global.f:
	ln -f -s $(NIMHOME)/nimrod/global.f .
input.f:
	ln -f -s $(NIMHOME)/nimset/input.f .
parbcast_input.f:
	ln -f -s $(NIMHOME)/nimrod/parbcast_input.f .
dump.f:
	ln -f -s $(NIMHOME)/nimrod/dump.f .
dump_reset.f:
	ln -f -s $(NIMHOME)/nimset/dump_reset.f .
eq_swap.f:
	ln -f -s $(NIMHOME)/nimset/eq_swap.f .
contour.f:
	ln -f -s $(NIMHOME)/nimplot/contour.f
diagnose.f:
	ln -f -s $(NIMHOME)/nimplot/diagnose.f
plot_io.f:
	ln -f -s $(NIMHOME)/nimplot/plot_io.f
nim_stop.f:
	ln -f -s $(NIMHOME)/nimeq/nim_stop.f
nimeq_all.f:
	ln -f -s $(NIMHOME)/nimeq/nimeq_all.f
nimeq_btr.f:
	ln -f -s $(NIMHOME)/nimeq/nimeq_btr.f
nimeq_data.f:
	ln -f -s $(NIMHOME)/nimeq/nimeq_data.f
nimeq_srf_ints.f:
	ln -f -s $(NIMHOME)/nimeq/nimeq_srf_ints.f
nimeq_free.f:
	ln -f -s $(NIMHOME)/nimeq/nimeq_free.f
nimeq.f:
	ln -f -s $(NIMHOME)/nimeq/nimeq.f
matrix_storage.f:
	ln -f -s $(NIMHOME)/nimeq/matrix_storage.f


cleanlinks:
	rm input.f dump.f physdat.f global.f contour.f diagnose.f \
	plot_io.f eq_swap.f dump_reset.f parbcast_input.f nim_stop.f \
	nimeq_all.f nimeq_btr.f nimeq_srf_ints.f nimeq_data.f nimeq_free.f \
	nimeq.f matrix_storage.f
