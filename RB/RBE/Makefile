default : RBE

include Makefile.def

install: Makefile.def.orig MAKEFILE_DEF
	@make install_cont;

Makefile.def.orig:
	mv Makefile.def Makefile.def.orig
	cp Makefile.def.orig Makefile.def

MAKEFILE_DEF:
	@(if [ "$(STANDALONE)" != "NO" ]; then \
		echo RBDIR=`pwd`                        >  Makefile.def; \
		echo OS=`uname`                         >> Makefile.def; \
		echo STANDALONE=${STANDALONE}           >> Makefile.def; \
		cat src/Makefile.def                    >> Makefile.def; \
	fi);

install_cont: 
	@(if [ "$(STANDALONE)" != "NO" ]; then \
		cp -f share/build/Makefile.${OS}${COMPILER} Makefile.conf; \
		cd share; make install;\
		cd util;  make install;\
	else \
		echo include $(DIR)/Makefile.conf > Makefile.conf; \
	fi);
	@(if [ -f src/Makefile.RULES.${OS}${COMPILER} ]; then                \
		cp -f src/Makefile.RULES.${OS}${COMPILER} src/Makefile.RULES;\
	else \
		rm -f src/Makefile.RULES; touch src/Makefile.RULES; \
	fi);
	touch src/Makefile.DEPEND

#
#       General Housekeeping
#

RBE:
	@cd ${SHAREDIR};  make LIB
	@cd ${TIMINGDIR}; make LIB
	@cd src;          make RBE

LIB:
	cd src; make LIB

test:	RBE
	cd run; ./rbe.exe
	make test_compare

test_compare:
	gunzip -c test/2000f223_e.fls.gz > test/2000f223_e.fls
	${SCRIPTDIR}/DiffNum.pl -r=0.001 -a=1e-10 \
		test/2000f223_e.fls run/2000f223_e.fls

clean:
	@touch src/Makefile.DEPEND src/Makefile.RULES
	@cd src; make clean
	@(if [ -d util ];  then cd util;  make clean; fi);
	@(if [ -d share ]; then cd share; make clean; fi);

distclean: clean
	@(if [ -d util  ]; then cd util;  make distclean; fi);
	@(if [ -d share ]; then cd share; make distclean; fi);
	@cd src; make distclean
	rm -f Makefile.conf Makefile.def *~
	mv Makefile.def.orig Makefile.def

#
#       Create run directories
#
rundir:
	mkdir -p ${RUNDIR}/RB
	cd ${RUNDIR}/RB; \
		mkdir restartOUT output; \
		ln -s restartOUT restartIN;
	@(if [ "$(STANDALONE)" != "NO" ]; then \
		cd ${RUNDIR} ; \
		ln -s ${BINDIR}/rbe.exe .   ; \
		ln -s ../input/2000_223.* . ; \
		ln -s ../input/w2k.dat .    ;\
		ln -s ../input/rbe_e.fin .  ; \
		cp ../input/rbe_swmf.dat .  ; \
		touch core ; chmod 444 core ; \
		ln -s RB/* .; \
	fi);

