include Makefile.def

install:
	touch src/Makefile.DEPEND

LIB:
	cd src; make LIB

rundir:
	mkdir -p ${RUNDIR}/IE/Output
	cd ${RUNDIR}/IE; ln -s ${DIR}/UA/GITM2/srcData Input; cd ${DIR}
	cd ${RUNDIR}; ln -s ${DIR}/UA/GITM2/src/PostProcess.exe PostGITM.exe; cd ${DIR}
	cd ${DIR}/UA/GITM2 ; make POST ; cd ${DIR}
	cd ${RUNDIR}/IE; ln -s ${IEDIR}/src/pIE .; cd ${DIR}

clean:
	@touch src/Makefile.DEPEND src/Makefile.RULES
	cd src; make clean

distclean: 
	./Config.pl -uninstall

allclean:
	@touch src/Makefile.DEPEND src/Makefile.RULES
	cd src; make distclean
	rm -f *~

test:
	echo "There is no test for IE/RIM" > notest.diff

test_tmp:
	cd ${RUNDIR}; cp ${DIR}/Param/LAYOUT.in.test.RIM.Weimer LAYOUT.in; cd ${DIR}
	cd ${RUNDIR}; cp ${DIR}/Param/PARAM.in.test.RIM.Weimer PARAM.in; cd ${DIR}
	cd ${RUNDIR}; ./SWMF.exe
