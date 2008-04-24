include ../../Makefile.def

check_src:
	@(if([ ! -f src/Makefile ]); then \
		echo "The src/Makefile is missing!"; \
		echo "make IHBATSRUS"; \
		exit 1; \
	fi)

INSTALLFILES = src/Makefile.DEPEND src/Makefile.RULES srcInterface/Makefile.DEPEND

install:
	@(if([ -f src/Makefile ]); then \
		touch ${INSTALLFILES}; \
		cp -f src/ModSize_orig.f90 src/ModSize.f90; \
	fi)

LIB: check_src
	cd src; make LIB
	cd srcInterface; make LIB

test:
	@echo IH/BATSRUS is tested in GM/BATSRUS

rundir: check_src
	mkdir ${RUNDIR}/IH
	cd ${RUNDIR}/IH; \
		mkdir restartIN restartOUT IO2; \
		ln -s ${BINDIR}/PostIDL.exe .; \
		ln -s ${BINDIR}/PostSPH.exe .; \
		ln -s ${BINDIR}/EARTH_TRAJ.exe .; \
		cp    ${DIR}/GM/BATSRUS/Scripts/IDL/pIDL .; \
		cp    ${DIR}/GM/BATSRUS/Scripts/TEC/pTEC .; \
		cp    ${DIR}/GM/BATSRUS/Param/HELIOSPHERE/EARTH_TRAJ.in .

clean:  
	@(if([ -f src/Makefile ]); then \
	   touch ${INSTALLFILES}; \
	   cd src; make clean; cd ..; \
	   cd srcInterface; make clean; \
	fi)

distclean: 
	@(if([ -f Makefile.def ]); then \
		touch srcInterface/Makefile.DEPEND; \
		cd srcInterface; make distclean; \
	fi)
	rm -rf src srcUser PARAM.XML RenameList.pl Config.pl *~
