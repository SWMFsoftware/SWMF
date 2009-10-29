include ../../Makefile.def

check_src:
	@(if([ ! -f src/Makefile ]); then \
		echo "The src/Makefile is missing!"; \
		echo "make OHBATSRUS"; \
		exit 1; \
	fi)

INSTALLFILES = src/Makefile.DEPEND src/Makefile.RULES srcInterface/Makefile.DEPEND

install:
	@(if([ -f src/Makefile ]); then \
		touch ${INSTALLFILES}; \
		cp -f src/ModSize_orig.f90 src/ModSize.f90; \
		cd src; make STATIC; \
	fi)

LIB: check_src
	cd src; make LIB
	cd srcInterface; make LIB

test:
	@echo OH/BATSRUS is tested in GM/BATSRUS

rundir: check_src
	mkdir ${RUNDIR}/OH
	cd ${RUNDIR}/OH; \
		mkdir restartIN restartOUT IO2; \
		ln -s ${BINDIR}/PostIDL.exe .; \
		ln -s ${BINDIR}/PostSPH.exe .; \
		cp    ${DIR}/GM/BATSRUS/Scripts/IDL/pIDL .; \
		cp    ${DIR}/GM/BATSRUS/Scripts/TEC/pTEC .

clean:  
	@(if([ -f src/Makefile ]); then \
	   touch ${INSTALLFILES}; \
	   cd src; make clean; cd ..; \
	   cd srcInterface; make clean; \
	fi)

distclean: 
	./Config.pl -uninstall

allclean:
	@(if([ -f Makefile.def ]); then \
		touch srcInterface/Makefile.DEPEND; \
		cd srcInterface; make distclean; \
	fi)
	rm -rf src srcBATL srcUser srcEquation \
		PARAM.XML RenameList.pl Config.pl *~
