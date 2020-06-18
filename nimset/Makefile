#This is the nimset makefile.

NIMHOME = ..
COREDIR = $(NIMHOME)/nimcore

include $(NIMHOME)/make.inc

LIBS = $(NIMHOME)/nimlib/libnimlib.a \
       $(COREDIR)/libnimbltype.a

MODDIR = -I$(NIMHOME)/nimlib -I$(COREDIR)

OBJS = \
	physdat.o \
	input.o \
	dump.o \
	dump_reset.o \
	polar_init.o \
	nimset_init.o \
	physics_init.o \
	nim_stop.o \
	diagnose.o \
        eq_swap.o \
	nimset.o

OBJS_MS = \
        nim_stop.o \
        input.o \
        physdat.o \
        dump.o \
        mode_scramble.o

OBJS_ST = \
	physdat.o \
	input.o \
	dump.o \
	polar_init.o \
	nimset_init.o \
	physics_init.o \
	nim_stop.o \
	diagnose.o \
	region_type.o \
        eq_swap.o \
	stitch.o

# The following "core" objects are not in a library.

CORE_OBJS = \
	$(COREDIR)/seam_storage.o \
	$(COREDIR)/aph_eval.o \
	$(COREDIR)/fields.o

nimset: $(OBJS)
	$(LDR) -o nimset $(LDR_FLAGS) $(CORE_OBJS) $(OBJS) $(LIBS)

mode_scramble: $(OBJS_MS)
	$(LDR) -o mode_scramble $(LDR_FLAGS) $(CORE_OBJS) $(OBJS_MS) $(LIBS)

stitch: $(OBJS_ST)
	$(LDR) -o stitch $(LDR_FLAGS) $(CORE_OBJS) $(OBJS_ST) $(LIBS)

clean:
	-rm *.o *.mod *.d

realclean: clean
	-rm nimset mode_scramble stitch

.f.o: 
	$(F90) -c $(FFLAGS) $(MODDIR) $*.f

# Module dependencies within this directory

dump.o: input.o
dump_reset.o: input.o
polar_init.o: input.o physdat.o
nimset_init.o: input.o physdat.o polar_init.o
poly_init.o: input.o physdat.o
physics_init.o: input.o physdat.o polar_init.o 
diagnose.o: input.o
poly.o: input.o
eq_swap.o: physdat.o
nimset.o: input.o dump.o nimset_init.o  \
	physics_init.o polar_init.o diagnose.o \
	dump_reset.o eq_swap.o
stitch.o: input.o dump.o \
	nimset_init.o diagnose.o region_type.o eq_swap.o
