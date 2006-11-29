include ../../Makefile.def

install:
	touch src/Makefile.DEPEND

rundir: ${RUNDIR}/ESMF_SWMF.input ${RUNDIR}/ESMF_SWMF.exe

${RUNDIR}/ESMF_SWMF.input:
	cp src/ESMF_SWMF.input ${RUNDIR}/

${RUNDIR}/ESMF_SWMF.exe:
	@(if perl -e 'exit(-l "${RUNDIR}/ESMF_SWMF.exe")' ; then \
	   echo "Creating link to ${BINDIR}/ESMF_SWMF.exe"; \
	   ln -s ${BINDIR}/ESMF_SWMF.exe ${RUNDIR}/ESMF_SWMF.exe; \
	fi)

clean:
	@(if([ "${ESMF_DIR}" != "" ]); then \
		cd src; touch Makefile.DEPEND; make clean; \
	fi)

distclean:
	@(if([ "${ESMF_DIR}" != "" ]); then \
		cd src; touch Makefile.DEPEND; make distclean; \
	fi)
	@rm -f src/Makefile.DEPEND
