include ../../Makefile.def

check_src:
	@(if([ ! -f src/Makefile ]); then \
		echo "The src/Makefile is missing!"; \
		echo "make SCBATSRUS"; \
		exit 1; \
	fi)

install:
	@(if([ -f src/Makefile ]); then \
		cp -f src/ModOptimizeParam_orig.f90 src/ModOptimizeParam.f90; \
		cp -f srcBATL/BATL_size_orig.f90 srcBATL/BATL_size.f90; \
	fi)

newtr:
	@(if([ -f src/ModThreadedLC_New.f90 ]); then \
		cd src; mv ModThreadedLC.f90 ModThreadedLC_Old.f90; \
		mv ModFieldLineThread.f90 ModFieldLineThread_Old.f90; \
		mv ModThreadedLC_New.f90 ModThreadedLC.f90; \
		mv ModFieldLineThread_New.f90 ModFieldLineThread.f90; \
		touch ModFieldLineThread.f90; \
	fi)

oldtr:
	@(if([ -f src/ModThreadedLC_Old.f90 ]); then \
		cd src; mv ModThreadedLC.f90 ModThreadedLC_New.f90; \
		mv ModFieldLineThread.f90 ModFieldLineThread_New.f90; \
		mv ModThreadedLC_Old.f90 ModThreadedLC.f90; \
		mv ModFieldLineThread_Old.f90 ModFieldLineThread.f90; \
		touch ModFieldLineThread.f90; \
	fi)

LIB: check_src
	cd srcBATL; $(MAKE) LIB
	cd src; $(MAKE) LIB
	cd srcInterface; $(MAKE) LIB

test:
	@echo SC/BATSRUS is tested in GM/BATSRUS

rundir: check_src
	mkdir ${RUNDIR}/SC
	cd ${RUNDIR}/SC; \
		mkdir restartIN restartOUT IO2; \
		ln -s IO2 plots; \
		ln -s ${BINDIR}/PostIDL.exe .; \
		ln -s ${BINDIR}/FDIPS.exe .; \
		ln -s ${BINDIR}/HARMONICS.exe .; \
		ln -s ${BINDIR}/CONVERTHARMONICS.exe .; \
		ln -s ${BINDIR}/FRMAGNETOGRAM.exe .; \
		ln -s ${BINDIR}/ELIMINATECMEBR.exe .; \
		cp ${MAGNETOGRAMDIR}/HARMONICS.in.orig \
		HARMONICS.in; \
		cp ${MAGNETOGRAMDIR}/HARMONICSGRID.in.orig \
		HARMONICSGRID.in; \
		ln -s ${MAGNETOGRAMDIR}/remap_magnetogram.py .; \
		ln -s ${MAGNETOGRAMDIR}/read_magnetogram.pro .; \
		ln -s ${MAGNETOGRAMDIR}/readfits.pro .; \
		ln -s ${MAGNETOGRAMDIR}/*SETUP*py .; \
		ln -s ${MAGNETOGRAMDIR}/*SETUP*pro .; \
		cp    ${DIR}/share/Scripts/pIDL .; \
		cp    ${DIR}/GM/BATSRUS/Scripts/TEC/pTEC .; \
		ln -s ${DIR}/GM/BATSRUS/Param/CORONA Param; \
		ln -s ${DIR}/GM/BATSRUS/data/TRAJECTORY .

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
