#This is the nimlib makefile.

NIMHOME = ..

include $(NIMHOME)/make.inc

NIMLIB_OBJS = \
	io.o \
	local.o \
	bessel.o \
	spline.o \
	polynomials.o \
        compl_ell.o \
	bicube.o \
	lagrange_quad.o \
	lagr_disc_quad.o \
	modal_disc_quad.o \
	tri_linear.o \
	math_tran.o \
        vector_type_mod.o \
	edge_type_mod.o \
	cell_type_mod.o \
	eldata_type_mod.o \
	tecplot_mod.o \
	vtk_mod.o

NIMMPI_OBJS = \
	mpi_link.o

COBJS = local_system.o

all: nimlib nimmpi

nimlib: $(NIMLIB_OBJS) $(COBJS)
	$(AR) lib$(@).a $(NIMLIB_OBJS) $(COBJS)
	$(RANLIB) lib$(@).a

nimmpi: $(NIMMPI_OBJS) $(COBJS)
	$(AR) lib$(@).a $(NIMMPI_OBJS) $(COBJS)
	$(RANLIB) lib$(@).a

clean:
	-rm *.o *.mod *.d

realclean: clean cleanlinks
	-rm libnimlib.a
	-rm libnimmpi.a

.c.o: 
	$(CC) -c $(CFLAGS) $*.c

.f.o: 
	$(F90) -c $(FFLAGS) $*.f

tecplot_mod.o:tecplot_mod.f
	$(F90) -c $(FFLAGS_NODBL) $*.f

vtk_mod.o:vtk_mod.f
	$(F90) -c $(FFLAGS_NODBL) $*.f

# Module dependencies

local.o: io.o
bessel.o: local.o
spline.o: local.o
polynomials.o: local.o
compl_ell.o: local.o
bicube.o: local.o spline.o
lagrange_quad.o: local.o
lagr_disc_quad.o: lagrange_quad.o
modal_disc_quad.o: local.o
tri_linear.o: local.o
math_tran.o: local.o
vector_type_mod.o: local.o bicube.o lagrange_quad.o modal_disc_quad.o
edge_type_mod.o: local.o
cell_type_mod.o: local.o
eldata_type_mod.o: local.o
tecplot_mod.o: local.o eldata_type_mod.o
vtk_mod.o: local.o eldata_type_mod.o

# links (Note that local is a machine-dependent link.)

local.f:
	ln -f -s $(LOCAL_FILE) local.f

mpi_link.f:
	ln -f -s $(MPILINK) mpi_link.f

cleanlinks:
	rm local.f mpi_link.f
