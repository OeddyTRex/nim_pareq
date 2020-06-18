#This is the nimbase makefile.

#The dependencies listed here assume that nimlib and
#nimcore modules (in the nimlib and nimcore directories)
#are already current.

NIMHOME = ..
COREDIR = $(NIMHOME)/nimcore
NIMLIBDIR = $(NIMHOME)/nimlib

include $(NIMHOME)/make.inc

LIBS = $(NIMLIBDIR)/libnimlib.a $(COREDIR)/libnimbltype.a
MODDIR = -I$(NIMLIBDIR) -I$(COREDIR)

OBJS = \
	input.o \
	dump.o \
	dumpb.o \
	nim_stop.o \
	nimbase.o

OBJS_COMB = \
	input.o \
	dump.o \
	dumpb.o \
	nim_stop.o \
	nimcomb.o

CORE_OBJS = \
	$(COREDIR)/fields.o \
	$(COREDIR)/seam_storage.o

nimbase: $(OBJS)
	$(LDR) -o nimbase $(LDR_FLAGS) $(CORE_OBJS) $(OBJS) $(LIBS)

nimcomb: $(OBJS_COMB)
	$(LDR) -o nimcomb $(LDR_FLAGS) $(CORE_OBJS) $(OBJS_COMB) $(LIBS)

clean:
	-rm *.o *.mod *.d

realclean: clean cleanlinks
	-rm nimbase

.f.o: 
	$(F90) -c $(FFLAGS) $(MODDIR) $*.f

# Module dependencies within this directory

dump.o: input.o
dumpb.o: dump.o input.o
nimbase.o: dump.o input.o
nimcomb.o: dump.o input.o


# links
 
input.f:
	ln -f -s $(NIMHOME)/nimset/input.f input.f
dump.f:
	ln -f -s $(NIMHOME)/nimset/dump.f dump.f
nim_stop.f:
	ln -f -s $(NIMHOME)/nimset/nim_stop.f nim_stop.f
 
cleanlinks:
	rm input.f nim_stop.f dump.f
