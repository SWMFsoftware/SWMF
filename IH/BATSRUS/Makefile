include ../../Makefile.def

check_src:
	@(if([ ! -f src/Makefile ]); then \
		echo "The src/Makefile is missing!"; \
		echo "make IHBATSRUS"; \
		exit 1; \
	fi)

install: check_src
	touch src/Makefile.DEPEND srcInterface/Makefile.DEPEND
	cp -f src/Makefile.RULES.${OS}${COMPILER} src/Makefile.RULES
	cp -f src/Makefile.OPTIONS.ORIG src/Makefile.OPTIONS
	cp -f src/ModSize_${OS}.f90 src/ModSize.f90

LIB: check_src
	cd src; make LIB
	cd srcInterface; make LIB

PLOT = ???

rundir: check_src
	mkdir ${RUNDIR}/IH
	cd ${RUNDIR}/IH; mkdir restartIN restartOUT IO2
	cd ${RUNDIR}/IH; \
		ln -s ${DIR}/GM/BATSRUS/Scripts/Run/pp_${PLOT}/P*.exe .; \
		cp    ${DIR}/GM/BATSRUS/Scripts/Run/pp_${PLOT}/p??? .

clean:  
	@(if([ -f src/Makefile ]); then \
	   touch src/Makefile.DEPEND src/Makefile.RULES src/Makefile.OPTIONS; \
	   cd src; make clean; cd ..; \
	   touch srcInterface/Makefile.DEPEND; cd srcInterface; make clean; \
	fi)

distclean: clean
	@(if([ -f Makefile.def ]); then \
		touch srcInterface/Makefile.DEPEND; \
		cd srcInterface; make distclean; \
	fi)
	rm -rf  src Makefile.* PARAM.XML PARAM.pl GridSize.pl *~

veryclean: distclean
