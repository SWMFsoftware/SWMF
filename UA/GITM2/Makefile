default : GITM

include Makefile.def

ABDIR   = ${UADIR}/srcSphereAB
EDDIR   = ${UADIR}/srcRudy
IEDIR   = ${UADIR}/srcIE
IODIR   = ${UADIR}/srcIO
MAINDIR = ${UADIR}/src

install: Makefile.def.orig MAKEFILE_DEF
	@make install_cont;

Makefile.def.orig:
	mv Makefile.def Makefile.def.orig
	cp Makefile.def.orig Makefile.def

MAKEFILE_DEF:
	@(if [ "$(STANDALONE)" != "NO" ]; then \
		echo UADIR=`pwd`                        >  Makefile.def; \
		echo OS=`uname`                         >> Makefile.def; \
		echo STANDALONE=${STANDALONE}           >> Makefile.def; \
		cat srcMake/Makefile.def                >> Makefile.def; \
	fi);

PLANET=earth

install_cont: 
	@(if [ "$(STANDALONE)" != "NO" ]; then \
		cp -f share/build/Makefile.${OS}${COMPILER} Makefile.conf; \
		cd share; make install;\
	else \
		echo include $(DIR)/Makefile.conf > Makefile.conf; \
	fi);
	touch src/Makefile.DEPEND srcInterface/Makefile.DEPEND
	cd src; make STATIC
	./config.pl -${PLANET}

#
#       General Housekeeping
#

GITM:
	@cd ${SHAREDIR}; make LIB
	@cd $(ABDIR); make LIB
	@cd $(IEDIR); make LIB
	@cd $(EDDIR); make LIB
	@cd $(IODIR); make LIB
	@cd $(MAINDIR); make GITM.exe

LIB:
	cd $(ABDIR)     ; make                                LIB
	cd $(EDDIR)     ; make LIBPREV=${ABDIR}/libSphere.a   LIBADD
	cd $(IODIR)     ; make LIBPREV=${EDDIR}/libABIEED.a   LIBADD
	cd $(MAINDIR)   ; make LIBPREV=${IODIR}/libABIEEDIO.a LIB
	cd srcInterface ; make LIBPREV=${MAINDIR}/libUA.a     LIB

nompirun:
	make GITM
	cd run; ./GITM.exe

clean:
	@cd $(ABDIR); make clean
	@cd $(IEDIR); make clean
	@cd $(EDDIR); make clean
	@cd $(IODIR); make clean
	@cd $(MAINDIR); make clean
	@(if [ -d share ]; then cd share; make clean; fi);

distclean: clean
	@(if [ -d share ]; then cd share; make distclean; fi);
	rm -f Makefile.conf Makefile.def *~
	mv Makefile.def.orig Makefile.def

#
#       Create run directories
#
rundir:
	mkdir -p ${RUNDIR}/UA
	cd ${RUNDIR}/UA; mkdir data restartOUT
	ln -s ${UADIR}/srcData ${RUNDIR}/UA/DataIn
	@(if [ "$(STANDALONE)" != "NO" ]; then \
		cd ${RUNDIR} ; \
		ln -s ${BINDIR}/GITM.exe . ; \
		cp UA/DataIn/UAM.in . ; \
		touch core ; chmod 444 core ; \
		ln -s UA/* .; \
	fi);

