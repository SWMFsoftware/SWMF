
include ../Makefile.def

SEARCHDIR = -I${SHAREDIR} -I${EMPIRICALGMDIR} -I${EMPIRICALIEDIR} -I${DATAREADINDICESDIR} -I../srcSAMI3

include ../Makefile.conf
include Makefile.DEPEND
include Makefile.RULES

OBJECTS = \
        ModPlanet.o\
	ModCimi.o\
	ModCimiPlot.o\
	ModFieldTrace.o\
	ModGmCIMI.o\
	ModGrid.o\
	ModCoupleSAMI.o\
	ModIeCimi.o\
	ModImSat.o\
	ModImTime.o\
	ModInitialize.o\
	ModRestart.o\
	ModCimiBoundary.o\
	ModTsyInput.o\
	ModPrerunField.o\
	cimi.o\
	set_parameters.o\
	plot_fieldline.o\
	trace_dipole.o\
	geopack.o\
        ModWaveDiff.o\
        ModLstar.o

DEPEND:
	@${SCRIPTDIR}/depend.pl ${SEARCHDIR} ${OBJECTS}


MY_LIB = libCIMI.a

LIB: DEPEND
	make ${MY_LIB}
	@echo
	@echo ${MY_LIB} has been brought up to date.
	@echo

${MY_LIB}: ${OBJECTS}
	rm -f ${MY_LIB}
	${AR} ${MY_LIB} ${OBJECTS} 

# Object files not included into the library
OBJECTS_EXE = cimi_main.o

# Other requireed libraries
LIBSHARE  = ${LIBDIR}/libSHARE.a
LIBTIMING = ${LIBDIR}/libTIMING.a
LIBINDICES  = ${LIBDIR}/libINDICES.a
LIBEMPIRICALIE   = ${LIBDIR}/libEMPIRICALIE.a
LIBEMPIRICALGM = ${LIBDIR}/libEMPIRICALGM.a
LIBSAMI3 = ${LIBDIR}/libSAMI3.a

# Libraries should be compiled first, because modules are used in main
${OBJECTS_EXE}: ${LIBSHARE} ${LIBTIMING} ${MY_LIB}   

EXE = ${BINDIR}/cimi.exe

CIMI:
	@echo ' '
	@make ${EXE}
	@echo "Program ${EXE} has been brought up to date."
	@echo ' '

${EXE}: ${OBJECTS_EXE}
	rm -rf Tmp_; mkdir Tmp_
	cd Tmp_; \
		ar -x ../${MY_LIB}; \
		ar -x ${LIBTIMING}; \
		ar -x ${LIBSHARE}; \
		ar -x ${LIBINDICES}; \
		ar -x ${LIBEMPIRICALGM};\
		ar -x ${LIBEMPIRICALIE};
	${LINK.f90} -o ${EXE} ${OBJECTS_EXE} Tmp_/*.o ${LBLAS} ${Lflag1}
	rm -rf Tmp_

CIMI_orig:
	${COMPILE.f90} ${Lflag1} -o ${EXE} cimi_main.f90 ${OBJECTS} \
	-L${LIBDIR} -lTIMING -lSHARE -lNOMPI -lINDICES

#SAMI-CIMI targets
# Object files not included into the library
CIMI_SAMI_OBJECTS_EXE = main_cimi_sami.o

${CIMI_SAMI_OBJECTS_EXE}: ${LIBSHARE} ${LIBTIMING} ${LIBSAMI3} ${MY_LIB}   

CIMI_SAMI_EXE = ${BINDIR}/cimi_sami.exe

CIMI_SAMI:
	@echo ' '
	@make ${CIMI_SAMI_EXE}
	@echo "Program ${CIMI_SAMI_EXE} has been brought up to date."
	@echo ' '

${CIMI_SAMI_EXE}: ${CIMI_SAMI_OBJECTS_EXE}
	rm -rf Tmp_; mkdir Tmp_
	cd Tmp_; \
		ar -x ../${MY_LIB}; \
		ar -x ${LIBSAMI3}; \
		ar -x ${LIBTIMING}; \
		ar -x ${LIBSHARE}; \
		ar -x ${LIBINDICES}; \
		ar -x ${LIBEMPIRICALGM};\
		ar -x ${LIBEMPIRICALIE};
	${LINK.f90} -o ${CIMI_SAMI_EXE} ${CIMI_SAMI_OBJECTS_EXE} Tmp_/*.o ${LBLAS} ${Lflag1}
	rm -rf Tmp_

CIMI_SAMI_orig:
	${COMPILE.f90} ${Lflag1} -o ${EXE} main_cimi_sami.f90 ${OBJECTS} \
	-L${LIBDIR} -lTIMING -lSHARE -lNOMPI -lINDICES -lSAMI3

clean: cleanfiles

distclean: clean
	rm -f ModPlanet.f90 ModGrid.f90
