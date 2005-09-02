include ../../Makefile.def

install:
	touch src/Makefile.DEPEND

rundir: ${RUNDIR}/ESMF_SWMF.input ${RUNDIR}/ESMF_SWMF.exe
	cd ${RUNDIR}; 

${RUNDIR}/ESMF_SWMF.input:
	cp src/ESMF_SWMF.input ${RUNDIR}/

${RUNDIR}/ESMF_SWMF.exe:
	@(if perl -e 'exit(-l "${RUNDIR}/ESMF_SWMF.exe")' ; then \
	   echo "Creating link to ${BINDIR}/ESMF_SWMF.exe"; \
	   ln -s ${BINDIR}/ESMF_SWMF.exe ${RUNDIR}/ESMF_SWMF.exe; \
	fi)

clean:
	touch src/Makefile.DEPEND
	cd src; make clean

distclean:
	touch src/Makefile.DEPEND
	cd src; make distclean
