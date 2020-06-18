#This is the makefile for nimuw's nimcore libraries.

#The dependencies listed here assume that nimlib and nimmpi
#modules (both in the nimlib directory) are already current.

NIMHOME = ..

include $(NIMHOME)/make.inc

LIB_BLTYPE=nimbltype
LIB_BLOCK=nimblock
LIB_MAT=nimmatrix
LIB_ITER=nimiter
LIB_LOCATE=nimlocate
LIB_PAR=nimpar

BLTYPE_OBJS = \
	rblock_type_mod.o \
	tblock_type_mod.o

BLOCK_OBJS = \
	rblock.o \
	tblock.o

MAT_OBJS = \
	matrix_type_mod.o \
	factor_type_mod.o \
	matrix_mod.o

ITER_OBJS = \
	iter_utils.o \
	iter_precon_real.o \
	iter_precon_comp.o \
	iter_cg_real.o \
	iter_cg_comp.o \
	iter_cg.o \
	iter_gmres_c2d.o \
	iter_gmres_r2d.o \
	iter_precon_c3dto2d.o \
	iter_krylov_c3d.o

LOCATE_OBJS = \
	nim_locate.o

PAR_OBJS = \
	parallel.o \
	parallel_io.o

MODDIR = -I$(NIMHOME)/nimlib

all: $(LIB_BLTYPE) $(LIB_BLOCK) fields.o $(LIB_MAT) $(LIB_ITER)  \
	$(LIB_PAR) $(LIB_LOCATE) \
	iter_externals.o seam_storage.o pardata.o time.o edge.o \
	computation_pointers.o regularity.o boundary.o fft_mod.o \
	generic_evals.o surface.o aph_eval.o block_create_tang.o

$(LIB_BLTYPE): $(BLTYPE_OBJS)
	   $(AR) lib$(@).a $(BLTYPE_OBJS)
	   $(RANLIB) lib$(@).a

$(LIB_BLOCK): $(LIB_BLTYPE) $(LIB_MAT) $(BLOCK_OBJS)
	   $(AR) lib$(@).a $(BLOCK_OBJS)
	   $(RANLIB) lib$(@).a

$(LIB_MAT): $(MAT_OBJS)
	   $(AR) lib$(@).a $(MAT_OBJS)
	   $(RANLIB) lib$(@).a

$(LIB_ITER): $(ITER_OBJS)
	   $(AR) lib$(@).a $(ITER_OBJS)
	   $(RANLIB) lib$(@).a

$(LIB_LOCATE): $(LOCATE_OBJS)
	   $(AR) lib$(@).a $(LOCATE_OBJS)
	   $(RANLIB) lib$(@).a

$(LIB_PAR): $(PAR_OBJS)
	   $(AR) lib$(@).a $(PAR_OBJS)
	   $(RANLIB) lib$(@).a

clean:
	-rm *.o *.mod *.d

realclean: clean
	-rm libnimbltype.a
	-rm libnimblock.a
	-rm libnimmatrix.a
	-rm libnimiter.a
	-rm libnimpar.a
	-rm libnimlocate.a

.c.o: 
	$(CC) -c $(CFLAGS_SLU) $(CINCL_SLU) $*.c

.f.o: 
	$(F90) -c $(FFLAGS) $(MODDIR) $*.f

# module dependencies

nim_locate.o: rblock_type_mod.o tblock_type_mod.o

fields.o: rblock_type_mod.o tblock_type_mod.o
 
matrix_mod.o: matrix_type_mod.o edge.o seam_storage.o time.o
factor_type_mod.o: matrix_type_mod.o

iter_utils.o: matrix_mod.o pardata.o seam_storage.o factor_type_mod.o
iter_precon_real.o: matrix_type_mod.o matrix_mod.o edge.o seam_storage.o \
        pardata.o time.o iter_utils.o factor_type_mod.o
iter_cg_real.o: matrix_type_mod.o matrix_mod.o edge.o seam_storage.o \
        pardata.o time.o iter_utils.o iter_precon_real.o
iter_precon_comp.o: matrix_type_mod.o matrix_mod.o edge.o seam_storage.o \
        pardata.o time.o iter_utils.o factor_type_mod.o
iter_cg_comp.o: matrix_type_mod.o matrix_mod.o edge.o seam_storage.o \
        pardata.o time.o iter_utils.o iter_precon_comp.o
iter_cg.o: iter_cg_real.o iter_cg_comp.o
iter_gmres_c2d.o: iter_precon_comp.o pardata.o time.o edge.o \
        iter_utils.o seam_storage.o matrix_mod.o factor_type_mod.o
iter_gmres_r2d.o: iter_precon_real.o pardata.o time.o edge.o \
        iter_utils.o seam_storage.o matrix_mod.o factor_type_mod.o
iter_precon_c3dto2d.o: time.o matrix_mod.o factor_type_mod.o \
        seam_storage.o edge.o iter_precon_comp.o
iter_krylov_c3d.o: iter_cg.o pardata.o time.o matrix_mod.o \
        seam_storage.o edge.o matrix_type_mod.o time.o iter_precon_comp.o \
	factor_type_mod.o

iter_externals.o: seam_storage.o pardata.o fields.o edge.o

parallel_io.o: rblock_type_mod.o tblock_type_mod.o pardata.o

parallel.o: rblock_type_mod.o fields.o pardata.o seam_storage.o edge.o

edge.o: rblock_type_mod.o tblock_type_mod.o fields.o pardata.o time.o \
	computation_pointers.o seam_storage.o

time.o: pardata.o

fft_mod.o: time.o

generic_evals.o: rblock_type_mod.o tblock_type_mod.o

boundary.o: seam_storage.o
regularity.o: seam_storage.o matrix_type_mod.o
block_create_tang.o: seam_storage.o edge.o fields.o

surface.o: rblock_type_mod.o tblock_type_mod.o time.o
