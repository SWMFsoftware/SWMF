include Makefile.def

install: Makefile.def.orig MAKEFILE_DEF
	@make install_cont;

Makefile.def.orig:
	mv Makefile.def Makefile.def.orig
	cp Makefile.def.orig Makefile.def

MAKEFILE_DEF:
	@(if [ "$(STANDALONE)" != "NO" ]; then \
		echo IEDIR=`pwd`                        >  Makefile.def; \
		echo OS=`uname`                         >> Makefile.def; \
		echo STANDALONE=YES                     >> Makefile.def; \
		echo include `pwd`/src/Makefile.def     >> Makefile.def; \
	fi);

install_cont:
	@(if [ "$(STANDALONE)" != "NO" ]; then \
		cp -f share/build/Makefile.${OS}${COMPILER} Makefile.conf; \
	else \
		echo include $(DIR)/Makefile.conf > Makefile.conf; \
	fi);
	touch src/Makefile.DEPEND

LIB:
	cd src; make LIB

PIONO:
	cd src; make PostIONO.exe

rundir:
	mkdir -p ${RUNDIR}/IE/ionosphere
	cd ${RUNDIR}/IE; ln -s ${IEDIR}/Scripts/pION .
	cd ${RUNDIR}/IE; ln -s ${BINDIR}/PostIONO.exe .

clean:
	touch src/Makefile.DEPEND
	cd src; make clean

distclean:
	touch src/Makefile.DEPEND
	cd src; make distclean
	rm -f Makefile.conf Makefile.def *~
	mv Makefile.def.orig Makefile.def
