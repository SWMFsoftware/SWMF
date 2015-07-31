#  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
#  For more information, see http://csem.engin.umich.edu/tools/swmf


SHELL =/bin/sh

default : HEIDI.exe

include ../Makefile.def

SEARCHDIR = -I${SHAREDIR} -I${DATAREADINDICESDIR} -I${EMPIRICALIEDIR}

include ../Makefile.conf
include Makefile.DEPEND

include Makefile.RULES

OBJECTS = \
	ModHeidiSize.o\
	ModHeidiMain.o\
	ModHeidiCurrents.o\
	ModHeidiDGCPM.o\
	ModHeidiDrifts.o\
	ModHeidiIO.o\
	ModHeidiWaves.o\
	ModInit.o\
	ModIonoHeidi.o\
	ModHeidiSatellites.o\
	ModProcHEIDI.o\
	ModHeidiInput.o\
	ModHeidiBField.o\
	ModHeidiNeutralH.o\
	ModHeidiBACoefficients.o\
	ModHeidiAsymStretchedB.o\
	heidi_emudt.o\
	heidi_integrals.o\
	heidi_setup.o\
	heidi_icbc.o\
	heidi_coefficients.o\
	heidi_operators.o\
	heidi_currents.o\
	heidi_output.o\
	heidi_waves.o\
	recipes.o\
	pbo.o\
	wgrdar.o\
	get_potential.o \
	heidi_check.o \
	heidi_init.o \
	heidi_run.o \
	read_ring_current.o\
	write_ring_current.o\
	epencalc.o\
	iono_init.o\
	HModel.o

EXEOBJECTS = heidi_main.o ${OBJECTS}

DEPEND:
	@perl ${SCRIPTDIR}/depend.pl ${SEARCHDIR} ${EXEOBJECTS}

MY_LIB = libHEIDI.a



LIB : DEPEND
	make ${MY_LIB}
	@echo
	@echo ${MY_LIB} has been brought up to date.
	@echo

${MY_LIB}: ${OBJECTS}
	rm -f ${MY_LIB}
	${AR} ${MY_LIB} ${OBJECTS} 

# main.o: ${LIBSHARE}

EXE = ${BINDIR}/HEIDI.exe

HEIDI:  DEPEND
	@make ${EXE}
	@echo "${EXE} has been created"

${EXE}: ${EXEOBJECTS}
	rm -rf Tmp_; mkdir Tmp_
	cd Tmp_; \
		ar -x ${LIBDIR}/libEMPIRICALIE.a; \
		ar -x ${LIBDIR}/libINDICES.a; \
		ar -x ${LIBDIR}/libSHARE.a;\
		ar -x ${LIBDIR}/libTIMING.a
	${LINK.f90}  -o ${EXE} ${EXEOBJECTS} Tmp_/*.o ${LBLAS} ${Lflag1}
	rm -rf Tmp_

clean: cleanfiles

distclean: clean
	rm -f Makefile.DEPEND Makefile.RULES



