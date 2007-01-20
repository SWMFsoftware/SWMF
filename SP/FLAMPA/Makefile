
include Makefile.def 

install: 
	touch src/Makefile.DEPEND srcInterface/Makefile.DEPEND

LIB:    install
	cd src;          make LIB
	cd srcInterface; make LIB

rundir: 
	mkdir -p ${RUNDIR}/SP ${RUNDIR}/SP/restartOUT ${RUNDIR}/SP/restartIN ${RUNDIR}/SP/MHDATA ${RUNDIR}/SP/IO2
	cd ${RUNDIR}; ln -s ${DIR}/bin/SEP.exe .

clean:  install
	@(if [ -r "Makefile.conf" ]; then \
		cd src; make clean; \
		cd ../srcInterface; make clean; \
	fi)

distclean: install
	@(if [ -r "Makefile.conf" ]; then \
		cd src; make distclean; \
		cd ../srcInterface; make distclean; \
		rm -f ../Makefile.conf; \
	fi)
	rm -f Makefile.def *~
