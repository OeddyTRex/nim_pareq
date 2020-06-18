#Global makefile for the nimrod suite of codes.

NIMHOME='..'

include make.inc

.IGNORE:

RUNDIR=rundir

all: solex fluxgridex nimsetex nimrodex nimplotex xlogex \
     nimflex nimeqex stitchex runlinks

# executables:

solex: nimliblib
	cd fluxgrid; make sol

fluxgridex: nimliblib
	cd fluxgrid; make fluxgrid

nimsetex: nimliblib coreset
	cd nimset; make 

stitchex: nimliblib coreset
	cd nimset; make stitch 

nimrodex: nimliblib externallibs coreobjs
	cd nimrod; make 

nimplotex: nimliblib externalstubs coreobjs
	cd nimplot; make nimplot

xlogex: nimliblib
	cd nimplot; make xlog

nimeqex: nimliblib externalforeq coreobjs
	cd nimeq; make

nimflex: nimliblib coreflex
	cd nimfl; make 

# included libraries and other common object files:

nimliblib:
	cd nimlib; make 

nimlocate:
	cd nimcore; make nimlocate 

externallibs:
	cd externals; make dummies $(SLU_LINK)

externalforeq:
	cd externals; make dummies sslu_link $(SLU_LINK)

externalstubs:
	cd externals; make dummies

coreobjs:
	cd nimcore; make

coreset:
	cd nimcore; make nimbltype seam_storage.o fields.o aph_eval.o

coreflex:
	cd nimcore; make nimbltype nimpar seam_storage.o fields.o time.o

# run links

$(RUNDIR):
	mkdir $(RUNDIR)

runlinks: executablel

executablel: $(RUNDIR)
	cd $(RUNDIR); ln -f -s $(NIMHOME)/fluxgrid/fluxgrid .
	cd $(RUNDIR); ln -f -s $(NIMHOME)/fluxgrid/sol .
	cd $(RUNDIR); ln -f -s $(NIMHOME)/nimset/nimset .
	cd $(RUNDIR); ln -f -s $(NIMHOME)/nimrod/nimrod .
	cd $(RUNDIR); ln -f -s $(NIMHOME)/nimplot/nimplot .
	cd $(RUNDIR); ln -f -s $(NIMHOME)/nimplot/xlog .
	cd $(RUNDIR); ln -f -s $(NIMHOME)/nimfl/nimfl .
	cd $(RUNDIR); ln -f -s $(NIMHOME)/nimeq/nimeq .

# cleaners

clean:
	cd nimlib; make clean
	cd fluxgrid; make clean
	cd nimset; make clean
	cd nimrod; make clean
	cd nimplot; make clean
	cd nimeq; make clean
	cd nimfl; make clean
	cd externals; make clean
	cd nimcore; make clean
	cd nimbase; make clean
	cd dump_convert; make clean

realcleanex:
	cd nimlib; make realclean
	cd fluxgrid; make realclean
	cd nimset; make realclean
	cd nimrod; make realclean
	cd nimplot; make realclean
	cd nimeq; make realclean
	cd nimfl; make realclean
	cd externals; make realclean
	cd nimcore; make realclean
	cd nimbase; make realclean
	cd dump_convert; make realclean

cleanlinks:
	-rm $(RUNDIR)/sol
	-rm $(RUNDIR)/fluxgrid
	-rm $(RUNDIR)/nimset
	-rm $(RUNDIR)/nimrod
	-rm $(RUNDIR)/xlog
	-rm $(RUNDIR)/nimplot
	-rm $(RUNDIR)/nimfl
	-rm $(RUNDIR)/nimeq

realclean: realcleanex cleanlinks

