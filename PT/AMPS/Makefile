SHELL=/bin/sh

DEFAULT_TARGET : amps

# These definitions may be overwritten by Makefile.def
SOURCES=src
WSD=srcTemp
InputFileAMPS=moon.input
SPICE=nospice

#Compiling with the CCMC's Kameleon
KAMELEON=nokameleon
BOOST=noboost
BATL=nobatl
TESTMODE=off
INTERFACE=off

include Makefile.def
include Makefile.conf


#include the local makefile (defined the AMPS' compiling variables)  
include Makefile.local

#the default value of the c++ compiler flags
SEARCH_C=-DMPI_ON -LANG:std -I${CWD}/${WSD}/pic -I${CWD}/${WSD}/main  -I${CWD}/${WSD}/meshAMR -I${CWD}/${WSD}/interface -I${CWD}/${WSD}/general -I${CWD}/${WSD}/models/electron_impact -I${CWD}/${WSD}/models/sputtering -I${CWD}/${WSD}/models/dust -I${CWD}/${WSD}/models/charge_exchange -I${CWD}/${WSD}/models/photolytic_reactions -I${CWD}/${WSD}/species -I${CWD}/${WSD}/models/exosphere -I${SPICE}/include -I${BOOST}/include -I${KAMELEON}/src -I${CWD}

#define the "compile kameleon' flag only when KAMELEON is used (to exclude including of the KAMELEON headers on machimes where KAMELEON is not installed) 
ifneq ($(KAMELEON),nokameleon)   
	SEARCH_C+=-D_PIC_COMPILE__KAMELEON_ 
endif
	

#the additional argument string for the fortran compiler
SEARCH_F=
#-fdefault-real-8 

# These definitions are inherited from Makefile.def and Makefile.conf
CC=${COMPILE.mpicxx}
CWD=${MYDIR}
AMPSLINKER=${CC}


AMPSLINKLIB= 

#include BATL-related libraries for linking
ifneq ($(BATL),nobatl)	
	AMPSLINKLIB+=${BATL}/lib/libREADAMR.a
	AMPSLINKLIB+=${BATL}/lib/libSHARE.a
	AMPSLINKER=${LINK.f90}

ifeq ($(COMPILE.f90),gfortran)  
	SEARCH_F+= -J${BATL}/share/include 
else ifeq ($(COMPILE.f90),mpif90)
	SEARCH_F+= -I${BATL}/share/include
else ifeq ($(COMPILE.f90),pgf90)
	SEARCH_F+= -module ${BATL}/share/include
else ifeq  ($(COMPILE.f90),ifort)
	SEARCH_F+= -module ${BATL}/share/include
else ifeq  ($(COMPILE.f90),nagfor)
	SEARCH_F+= -I${BATL}/share/include
endif

endif

# include interface with external FORTRAN subroutines
ifeq ($(INTERFACE),on)
	AMPSLINKER=${LINK.f90}
	AMPSLINKLIB+=${WSD}/interface/interface.a	
endif 


# when linking mixed C/C++ and FORTRAN code mpif90 is used for linking
# certain compilers (Intel, PGI) require an additional flag
# if main() subroutine is written in C/C++

ifeq (${AMPSLINKER},${LINK.f90})

ifeq ($(COMPILE.f90),pgf90)
	AMPSLINKER+= -Mnomain
else ifeq  ($(COMPILE.f90),ifort)
	AMPSLINKER+= -nofor-main
endif

endif

# This looks wrong
ifeq ($(STANDALONE), NO)
	AMPSLINKLIB+=${LIBDIR}/*
	AMPSLINKER=${LINK.f90}
endif

ifneq ($(KAMELEON),nokameleon)
	AMPSLINKLIB+=${KAMELEON}/lib/ccmc/libccmc.dylib 
endif

ifneq ($(SPICE),nospice)
	AMPSLINKLIB+=${SPICE}/lib/cspice.a
endif

ifeq ($(TESTMODE),on)  
	SEARCH_C+=-D _PIC_NIGHTLY_TEST_MODE_=_PIC_MODE_ON_
endif



install:
	@echo " " > Makefile.local
	@rm -f .ampsConfig.Settings
	@./Config.pl -application=Moon
	@echo "AMPS installed"

distclean:
	./Config.pl -uninstall

allclean: clean
	rm -rf main srcTemp *.input* amps Makefile.local .ampsConfig.Settings

rundir:
	mkdir -p ${RUNDIR}/PT
	cd ${RUNDIR}/PT; mkdir restartIN restartOUT plots



EXE=amps

#Lib=  -lm

#MPIRUN=mpirun -np 4
#RUNDIR=run

LIB_AMPS = ${WSD}/libAMPS.a

clean:
	rm -rf ${LIB_AMPS} ${WSD}
	@(if [ -d srcInterface ]; then cd srcInterface; $(MAKE) clean; fi);

tar:
	cd ../pic-tower/sources/general; rm -f *.o *.a
	cd ../pic-tower/sources/dsmc; rm -f *.o *.a
	tar -cvf sources.tar sources

${WSD}:
	./ampsConfig.pl -input ${InputFileAMPS} -no-compile

LIB: 
	@(if [ -d ${WSD} ]; then rm -rf ${WSD}; fi)
	$(MAKE) ${WSD}
	make LIB_after_build
	(if [ "$(STANDALONE)" == "NO" ]; then \
		cd srcInterface; make LIB SEARCH_C="${SEARCH_C}"; fi)

LIB_after_build: 
ifeq ($(INTERFACE),on)
	cd ${WSD}/interface; make SEARCH_C="${SEARCH_C}" SEARCH="${SEARCH_F}" 
endif
	cd ${WSD}/general;                     make SEARCH_C=
	cd ${WSD}/meshAMR;                     make SEARCH_C="${SEARCH_C}" 
	cd ${WSD}/pic;                         make SEARCH_C="${SEARCH_C}" SEARCH="${SEARCH_F}" 
	cd ${WSD}/species;                     make SEARCH_C="${SEARCH_C}"
	cd ${WSD}/models/electron_impact;      make SEARCH_C="${SEARCH_C}"
	cd ${WSD}/models/sputtering;           make SEARCH_C="${SEARCH_C}"
	cd ${WSD}/models/dust;                 make SEARCH_C="${SEARCH_C}"
	cd ${WSD}/models/charge_exchange;      make SEARCH_C="${SEARCH_C}"
	cd ${WSD}/models/photolytic_reactions; make SEARCH_C="${SEARCH_C}" 
#compile external modules
	$(foreach src, $(ExternalModules), (cd ${WSD}/$(src); make SEARCH_C="${SEARCH_C}")) 
	cd ${WSD}/main; make SEARCH_C="${SEARCH_C}"
	cp -f ${WSD}/main/mainlib.a ${WSD}/libAMPS.a
ifeq ($(SPICE),nospice)
	cd ${WSD}; ${AR} libAMPS.a general/*.o meshAMR/*.o pic/*.o species/*.o models/electron_impact/*.o models/sputtering/*.o models/dust/*.o models/charge_exchange/*.o models/photolytic_reactions/*.o
	$(foreach src, $(ExternalModules), (cd ${WSD}; ${AR} libAMPS.a $(src)/*.o))
else
	rm -rf ${WSD}/tmpSPICE
	mkdir ${WSD}/tmpSPICE
	cp ${SPICE}/lib/cspice.a ${WSD}/tmpSPICE
	cd ${WSD}/tmpSPICE; ar -x cspice.a
	cd ${WSD}; ${AR} libAMPS.a general/*.o meshAMR/*.o pic/*.o species/*.o models/electron_impact/*.o models/sputtering/*.o models/dust/*.o models/charge_exchange/*.o models/photolytic_reactions/*.o tmpSPICE/*.o 
	$(foreach src, $(ExternalModules), (cd ${WSD}; ${AR} libAMPS.a $(src)/*.o))
endif

.PHONY: amps
amps:
	@(if [ -d ${WSD} ]; then rm -rf ${WSD}; fi)
	$(MAKE) $(WSD)
	$(MAKE) amps_after_build

amps_after_build: LIB_after_build 
	@rm -f amps
	cd ${WSD}/main; make amps SEARCH_C="${SEARCH_C}"
	make amps_link

amps_link:
	${AMPSLINKER} -o amps srcTemp/main/main.a srcTemp/libAMPS.a \
		${CPPLIB} ${AMPSLINKLIB}

.PHONY: test
test:
	echo "Use make test_all to run the AMPS tests"
