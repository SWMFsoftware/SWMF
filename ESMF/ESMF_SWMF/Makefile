include ../../Makefile.def

ESMF3 = `cd ${ESMF_DIR}; make info | grep ESMF_VERSION_STRING | grep 3`

# Change source code for ESMF3 compatibility
install:
	touch src/Makefile.DEPEND
	@(if [ "$(ESMF3)" != "" ]; then \
	   perl -ni -e 'print unless /\.f90\.o|FC_MOD|FREENO/' src/Makefile; \
	   perl -pi -e \
	     's/type=ESMF_DATA_REAL, kind=ESMF_R8/typekind=ESMF_TYPEKIND_R8/' \
		src/ESMF_SWMF_Mod.f90; \
	   perl -pi -e \
	     's/(ESMF_(Grid|Cpl)CompCreate)\((\w+),/$$1\(parentvm=$$3,/' \
		src/ESMF_SWMF_Driver.f90 src/ESMF_SWMF_GridComp.f90; \
	fi)

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
