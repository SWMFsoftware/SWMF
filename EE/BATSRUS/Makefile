include ../../Makefile.def

check_src:
	@(if([ ! -f src/Makefile ]); then \
		echo "The src/Makefile is missing!"; \
		echo "make EEBATSRUS"; \
		exit 1; \
	fi)

install:
	@(if([ -f src/Makefile ]); then \
		cp -f src/ModOptimizeParam_orig.f90 src/ModOptimizeParam.f90; \
		cp -f srcBATL/BATL_size_orig.f90 srcBATL/BATL_size.f90; \
	fi)

LIB: check_src
	cd srcBATL; $(MAKE) LIB
	cd src; $(MAKE) LIB
	cd srcInterface; $(MAKE) LIB

test:
	@echo EE/BATSRUS is tested in GM/BATSRUS

rundir: check_src
	mkdir ${RUNDIR}/EE
	cd ${RUNDIR}/EE; \
		mkdir restartIN restartOUT IO2; \
		ln -s IO2 plots; \
		ln -s ${BINDIR}/PostIDL.exe .; \
		cp    ${DIR}/share/Scripts/pIDL .; \
		cp    ${DIR}/GM/BATSRUS/Scripts/TEC/pTEC .; \
		ln -s ${DIR}/GM/BATSRUS/Param/FLUXEMERGENCE Param

clean:  
	@(if([ -f src/Makefile ]); then \
	   cd src; make clean; cd ..; \
	   cd srcBATL; make clean; cd ..; \
	   cd srcInterface; make clean; \
	fi)

distclean: 
	./Config.pl -uninstall

allclean:
	@(if([ -f Makefile.def ]); then \
		cd srcInterface; make distclean; \
	fi)
	rm -rf src srcBATL srcUser srcEquation \
		PARAM.XML RenameList.pl Config.pl *~
