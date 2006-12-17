

include ../../Makefile.def


check_src:
	@(if([ ! -f src/Makefile ]); then \
		echo "The src/Makefile is missing!"; \
		echo "make SCBATSRUS"; \
		exit 1; \
	fi)

install:
	@(if([ -f src/Makefile ]); then \
		touch src/Makefile.DEPEND srcInterface/Makefile.DEPEND; \
		cp -f src/Makefile.RULES.${OS}${COMPILER} src/Makefile.RULES; \
		cp -f src/ModSize_orig.f90 src/ModSize.f90; \
	fi)

LIB: check_src
	cd src; make LIB
	cd srcInterface; make LIB

test:
	touch see_gm_batsrus.diff

rundir: check_src
	mkdir ${RUNDIR}/SC
	cd ${RUNDIR}/SC; \
		mkdir restartIN restartOUT IO2; \
		ln -s ${BINDIR}/PostIDL.exe .; \
		ln -s ${BINDIR}/PostSPH.exe .; \
		cp    ${DIR}/GM/BATSRUS/Scripts/IDL/pIDL .; \
		cp    ${DIR}/GM/BATSRUS/Scripts/TEC/pTEC .
	cp ../../Param/CR1935_WSO.dat ${RUNDIR}/SC/

clean:  
	@(if([ -f src/Makefile ]); then \
	   touch src/Makefile.DEPEND src/Makefile.RULES; \
	   cd src; make clean; cd ..;\
	   touch srcInterface/Makefile.DEPEND; cd srcInterface; make clean; \
	fi)

distclean: clean
	@(if([ -f Makefile.def ]); then \
		touch srcInterface/Makefile.DEPEND; \
		cd srcInterface; make distclean; \
	fi)
	rm -rf src srcUser
	rm -f Makefile.* PARAM.XML RenameList.pl GridSize.pl Options.pl *~

veryclean: distclean
