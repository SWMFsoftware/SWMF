include ../Makefile.def
include ../Makefile.conf
include Makefile.DEPEND
include Makefile.RULES

MODULES = \
	ModFSAM.o \
	ModPar.o \
	ModSolar.o \
	ModPhysunits.o \
	ModDomain.o \
	ModOutfile.o \
	ModCoeff.o \
	ModSundry.o \
	ModControl.o \
	ModGrid.o \
	ModBack.o \
	ModBlk.o \
	ModBoundary.o \
	ModField.o \
	ModDel.o \
	ModInterp.o \
	ModIoFSAM.o \
	ModGetQtys.o \
	ModBval.o \
	ModGetP.o \
	ModRHS.o \
	ModIteration.o \
	ModUserSetup.o \
	ModInitialization.o \
	ModFish.o

OBJ = main.o

DEPEND:
	@${SCRIPTDIR}/depend.pl ${SEARCHDIR} ${OBJ} ${MODULES}

${MODULES}: ${LIBSHARE}

$(OBJ): ${MODULES} ${LIBSHARE}

EXE = xfsam.exe

xfsam: DEPEND
	make EXE
	@echo "xfsam.exe has been created"

EXE: $(OBJ)
	${LINK.f90} $(FFLAGS) -o ${EXE} $(OBJ) $(MODULES) ${LIBFISH} \
		-L${LIBDIR} -L${SHAREDIR} -lSHARE ${Lflag1}

MY_LIB = libFSAM.a

LIB: DEPEND ${MODULES} 
	rm -f ${MY_LIB}
	${AR} ${MY_LIB} ${MODULES}
	@echo
	@echo ${MY_LIB} has been brought up to date.
	@echo

distclean: clean
	rm -f ModUserSetup.f90 ModPar.f90
