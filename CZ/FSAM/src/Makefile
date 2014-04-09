include ../Makefile.def
include ../Makefile.conf
include Makefile.DEPEND
include Makefile.RULES

LIBFISH = -L../input -lfishpak 

EXTRA = -fp-model strict
CFLAGS = ${Cflag3} ${EXTRA}
FFLAGS = ${Lflag1} ${OPT3} ${EXTRA}

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
	ModIO.o \
	ModGetQtys.o \
	ModBval.o \
	ModGetP.o \
	ModRHS.o \
	ModIteration.o \
	ModUserSetup.o \
	ModInitialization.o 

OBJ = main.o

DEPEND:
	@${SCRIPTDIR}/depend.pl ${SEARCH} ${OBJ} ${MODULES}

${MODULES}: ${LIBSHARE}

$(OBJ): ${MODULES} ${LIBSHARE}

EXE = xfsam.exe

.f90.o:
	${COMPILE.f90} ${CFLAGS} ${LIBFISH} $<

xfsam: $(OBJ)
	${LINK.f90} $(FFLAGS) -o ${EXE} $(OBJ) $(MODULES) ${LIBFISH} \
		-L${LIBDIR} -L${SHAREDIR} -lSHARE
	@echo "xfsam.exe has been created"

MY_LIB = libFSAM.a

LIB: DEPEND ${MODULES} 
	rm -f ${MY_LIB}
	${AR} ${MY_LIB} ${MODULES}
	@echo
	@echo ${MY_LIB} has been brought up to date.
	@echo


distclean: clean
	rm -f Makefile.DEPEND Makefile.RULES