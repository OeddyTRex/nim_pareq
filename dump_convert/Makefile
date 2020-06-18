#This is the dump_convert makefile.

NIMHOME = ..
COREDIR = $(NIMHOME)/nimcore

include $(NIMHOME)/make.inc

LIBS = $(NIMHOME)/nimlib/libnimlib.a
MODDIR = -I$(NIMHOME)/nimlib -I$(COREDIR)

OBJS = \
        nim_stop.o \
	rblock31_type_mod.o \
	tblock_type_mod.o \
	fields.o \
        seam_storage.o

OBJS2 = \
	dump_convert.o

OBJS3 = \
	input.o \
	physdat.o \
	dump_conv30to31.o

OBJS4 = \
	input.o \
	physdat.o \
	dump_conv31to32.o

all: dump_convert dump_conv30to31 dump_conv31to32

dump_convert: $(OBJS) $(OBJS2)
	$(LDR) -o dump_convert $(LDR_FLAGS) $(OBJS) $(OBJS2) $(LIBS)

dump_conv30to31: $(OBJS) $(OBJS3)
	$(LDR) -o dump_conv30to31 $(LDR_FLAGS) $(OBJS) $(OBJS3) $(LIBS)

dump_conv31to32: $(OBJS) $(OBJS4)
	$(LDR) -o dump_conv31to32 $(LDR_FLAGS) $(OBJS) $(OBJS4) $(LIBS)

clean:
	-rm *.o *.mod *.d

realclean: clean cleanlinks
	-rm dump_convert
	-rm dump_conv30to31
	-rm dump_conv31to32

.f.o: 
	$(F90) -c $(FFLAGS) $(MODDIR) $*.f

# Module dependencies

fields.o: rblock31_type_mod.o tblock_type_mod.o
dump_convert.o: rblock31_type_mod.o tblock_type_mod.o \
	fields.o seam_storage.o
dump_conv30to31.o: rblock31_type_mod.o tblock_type_mod.o \
	fields.o seam_storage.o input.o physdat.o
dump_conv31to32.o: rblock31_type_mod.o tblock_type_mod.o \
	fields.o seam_storage.o input.o physdat.o

input.f:
	ln -f -s $(NIMHOME)/nimset/input.f input.f
physdat.f:
	ln -f -s $(NIMHOME)/nimset/physdat.f physdat.f
nim_stop.f:
	ln -f -s $(NIMHOME)/nimset/nim_stop.f nim_stop.f
tblock_type_mod.f:
	ln -f -s $(COREDIR)/tblock_type_mod.f tblock_type_mod.f
fields.f:
	ln -f -s $(COREDIR)/fields.f fields.f
seam_storage.f:
	ln -f -s $(COREDIR)/seam_storage.f seam_storage.f

cleanlinks:
	-rm nim_stop.f tblock_type_mod.f \
	   fields.f seam_storage.f input.f physdat.f

