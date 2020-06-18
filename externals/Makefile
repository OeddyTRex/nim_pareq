#This is the externals makefile.

NIMHOME = ..

include $(NIMHOME)/make.inc

SSLU_OBJS = \
	c_fortran_zgssv.o \
	c_fortran_dgssv.o

SLUD_OBJS = \
	c_fortran_slugrid.o \
	c_fortran_pzglob.o \
	c_fortran_pdglob.o \
	c_fortran_pzloc.o \
	c_fortran_pdloc.o

DUM_SSLU_OBJS = \
	c_dum_sslu.o

DUM_SLUD_OBJS = \
	c_dum_slud.o

DUM_LAP_OBJS = \
	f90_dum_lapack.o

dummies: sslu_dummy slud_dummy lapack_dummy

sslu_link: slu_config  $(SSLU_OBJS)
	   $(AR) lib$(@).a $(SSLU_OBJS)
	   $(RANLIB) lib$(@).a

slud_link: slu_config  $(SLUD_OBJS)
	   $(AR) lib$(@).a $(SLUD_OBJS)
	   $(RANLIB) lib$(@).a

sslu_dummy: $(DUM_SSLU_OBJS)
	   $(AR) lib$(@).a $(DUM_SSLU_OBJS)
	   $(RANLIB) lib$(@).a

slud_dummy: $(DUM_SLUD_OBJS)
	   $(AR) lib$(@).a $(DUM_SLUD_OBJS)
	   $(RANLIB) lib$(@).a

lapack_dummy: $(DUM_LAP_OBJS)
	   $(AR) lib$(@).a $(DUM_LAP_OBJS)
	   $(RANLIB) lib$(@).a

clean:
	-rm *.o *.mod *.d

realclean: clean
	-rm libsslu_link.a
	-rm libslud_link.a
	-rm libsslu_dummy.a
	-rm libslud_dummy.a
	-rm liblapack_dummy.a
	-rm config.h

slu_config:
	sed s/defndefn/$(SLUD6_DEFN)/g<config.tmpl>config.h

.c.o: 
	$(CC) -c $(CFLAGS_SLU) $(CINCL_SLU) $*.c

.f.o: 
	$(F90) -c $(FFLAGS) $*.f

