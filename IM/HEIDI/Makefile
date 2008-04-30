#default : HEIDI

include Makefile.def

INSTALLFILES =  src/Makefile.DEPEND \
		src/Makefile.RULES 

install: 
	touch ${INSTALLFILES}

#
#       General Housekeeping
#

HEIDI:
	@cd ${NOMPIDIR};	make LIB
	@cd ${SHAREDIR};  	make LIB
	@cd ${TIMINGDIR}; 	make LIB 
	@cd ${EMPIRICALIEDIR};	make LIB
	@cd ${EMPIRICALGMDIR};	make LIB
	@cd src;	        make HEIDI

TESTDIR = run_test

test:
	@echo "There is no test for HEIDI" > notest.diff

clean: install
	@(if [ -r "Makefile.conf" ]; then  \
	cd src;                      make clean;\
	cd srcInterface;             make clean;\
	cd srcInterface/srcIONO;     make clean;\
	cd srcInterface/srcIE;       make clean;\
	fi)

distclean: install
	@(if [ -r "Makefile.conf" ]; then \
		cd src;                      make distclean;\
		cd srcInterface;             make distclean;\
		cd srcInterface/srcIONO;     make distclean;\
		cd srcInterface/srcIE;       make distclean;\
	fi)
	rm -f *~
#
#       Create run directories
#
rundir:
	mkdir -p ${RUNDIR}/IM
	cd ${RUNDIR}/IM; \
		mkdir restartOUT restartIN plots
	@(if [ "$(STANDALONE)" != "NO" ]; then \
		cd ${RUNDIR} ; \
		ln -s ${BINDIR}/HEIDI.exe .
		touch core ; chmod 444 core;\
	fi);
