

include ../../Makefile.def


check_src:
	@(if([ ! -f src/Makefile ]); then \
		echo "The src/Makefile is missing!"; \
		echo "make SCBATSRUS"; \
		exit 1; \
	fi)

install: check_src
	touch src/Makefile.DEPEND srcInterface/Makefile.DEPEND
	cp -f src/Makefile.RULES.${OS}${COMPILER} src/Makefile.RULES
	cp -f src/Makefile.OPTIONS.ORIG src/Makefile.OPTIONS
	cp -f src/ModSize_${OS}.f90 src/ModSize.f90
	cd src; make user_routines.f90


LIB: check_src
	cd src; make LIB
	cd srcInterface; make LIB

PLOT = ???

rundir: check_src
	mkdir ${RUNDIR}/SC
	cd ${RUNDIR}/SC; mkdir restartIN restartOUT IO2
	cd ${RUNDIR}/SC; ln -s ${SCDIR}/Scripts/Run/pp_${PLOT}/[Pp]* .

clean:  
	@(if([ -f src/Makefile ]); then \
	   touch src/Makefile.DEPEND src/Makefile.RULES src/Makefile.OPTIONS; \
	   cd src; make clean; cd ..;\
	   touch srcInterface/Makefile.DEPEND; cd srcInterface; make clean; \
	fi)

distclean: clean
	@(if([ -f src/Makefile ]); then \
	  cd src; make distclean; cd ..; \
	  cd srcInterface; make distclean; \
	cd ..; rm -f Makefile.*;\
	fi)

veryclean: distclean
	@(if([ -f src/Makefile ]); then \
	rm -rf  src \
	srcInterface/SC_wrapper.f90 \
	srcInterface/SC_get_for_sp.f90 \
	srcInterface/ModGridDescriptor.f90 \
	srcInterface/update_lagrangian_grid.f90 \
	srcInterface/Makefile.DEPEND \
	PARAM.XML \
	*.pl Test* \
	Idl Param \
	Scripts srcPostProc *~; \
	fi)


