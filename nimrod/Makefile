#This is the nimrod makefile.

#The dependencies listed here assume that nimlib, nimmpi,
#and nimcore modules (in the nimlib and nimcore directories)
#are already current.

NIMHOME = ..
COREDIR = $(NIMHOME)/nimcore
NIMLIBDIR = $(NIMHOME)/nimlib

include $(NIMHOME)/make.inc

MODDIR = -I$(NIMLIBDIR) -I$(COREDIR)

LIBS = $(COREDIR)/libnimpar.a $(NIMLIBDIR)/libnimmpi.a \
	$(COREDIR)/libnimblock.a $(COREDIR)/libnimbltype.a \
	$(COREDIR)/libnimiter.a $(COREDIR)/libnimmatrix.a \
	$(NIMLIBDIR)/libnimlib.a $(LAP_LIB) $(SLU_NIM_LIB) \
	$(MPI_LIB) $(BLAS_LIB)

OBJS = \
	physdat.o \
	input.o \
	parbcast_input.o \
	global.o \
        matrix_storage.o \
	dump.o \
	integrands_rhs.o \
	integrands_dot.o \
	integrands_mat.o \
	integrands.o \
	surface_ints.o \
	diagnostic_ints.o \
	extrap_mod.o \
	finite_element.o \
	threed_dot_mgt.o \
	nimrod_init.o \
	history.o \
	diagnose.o \
	utilities.o \
	closure_model.o \
	field_comps.o \
	threed_precon.o \
	nimrod.o

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
	$(COREDIR)/surface.o \
	$(COREDIR)/iter_externals.o

nimrod: $(OBJS)
	$(MPILDR) -o nimrod $(LDR_FLAGS) $(CORE_OBJS) $(OBJS) $(LIBS)

clean:
	-rm *.o *.mod *.d

realclean: clean cleanlinks
	-rm nimrod

.f.o: 
	$(MPIF90) -c $(FFLAGS) $(MODDIR) $*.f

# Module dependencies

parbcast_input.o: input.o
dump.o: input.o
integrands_rhs.o: physdat.o input.o global.o
integrands_dot.o: physdat.o input.o global.o
integrands_mat.o: physdat.o input.o global.o
integrands.o: integrands_rhs.o integrands_dot.o integrands_mat.o
surface_ints.o: physdat.o input.o global.o
diagnostic_ints.o: physdat.o input.o global.o 
extrap_mod.o: global.o input.o
finite_element.o: input.o global.o
threed_dot_mgt.o: input.o global.o integrands.o finite_element.o \
	surface_ints.o matrix_storage.o
nimrod_init.o: input.o global.o integrands.o matrix_storage.o \
	physdat.o
history.o: input.o global.o physdat.o diagnostic_ints.o
diagnose.o: input.o
utilities.o: physdat.o input.o global.o history.o dump.o \
	diagnose.o diagnostic_ints.o extrap_mod.o matrix_storage.o
closure_model.o: physdat.o input.o global.o
field_comps.o: physdat.o input.o global.o closure_model.o
threed_precon.o: matrix_storage.o input.o global.o
nimrod.o: physdat.o nimrod_init.o input.o global.o integrands.o \
	history.o dump.o diagnose.o extrap_mod.o utilities.o \
	finite_element.o matrix_storage.o \
	surface_ints.o field_comps.o 

# links (Note that mpi_link is a machine-dependent link.)

physdat.f:
	ln -f -s $(NIMHOME)/nimset/physdat.f physdat.f
input.f:
	ln -f -s $(NIMHOME)/nimset/input.f input.f

cleanlinks:
	rm physdat.f input.f
