
include Makefile.def 

install: Makefile.def.orig
	touch src/Makefile.DEPEND srcInterface/Makefile.DEPEND

Makefile.def.orig:
	mv Makefile.def Makefile.def.orig
	cp Makefile.def.orig Makefile.def

LIB:
	cd src; make LIB
	cd srcInterface; make LIB

test:
	echo "SP/Kota test is incomplete" > notest.diff

rundir: 
	mkdir -p ${RUNDIR}/SP ${RUNDIR}/SP/restartOUT ${RUNDIR}/SP/restartIN
	@(if [ "${STANDALONE}" != 'NO' ]; then 	\
		cd ${RUNDIR}; 			\
		ln -s ../bin/SEP.exe .; 	\
		ln -s ../SP/Kota/src/a.out .;	\
	fi)

clean: install
	@(if [ -r "Makefile.conf" ]; then  \
		cd src; make clean; \
		cd ../srcInterface; make clean; \
	fi)

distclean: install
	@(if [ -r "Makefile.conf" ]; then \
		cd src; make distclean; \
		cd ../srcInterface; make distclean; \
		rm -f Makefile.conf; \
	fi)
	rm -f Makefile.def *~
	mv Makefile.def.orig Makefile.def
