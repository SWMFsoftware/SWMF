#^CFG COPYRIGHT UM
#
#	Specific Rules for Darwin fortran NAG compiler
#

calc_chemistry.o: calc_chemistry.f90\
	ModSize.o \
	ModGITM.o \
	ModPlanet.o \
	ModRates.o \
	ModEUV.o \
	ModSources.o \
	ModInputs.o \
	ModConstants.o
	${COMPILE.f90} ${Cflag1} calc_chemistry.f90

calc_electrodynamics.o: calc_electrodynamics.f90 \
	ModElectrodynamics.o \
	ModGITM.o \
	ModInputs.o \
	${SHAREDIR}/ModMpi.o
	${COMPILE.f90} ${Cflag0} calc_electrodynamics.f90

iri90.o: iri90.f
	${COMPILE.f77} ${Cflag1} -132 iri90.f 

apex.o: apex.f
	${COMPILE.f77} ${Cflag1} -132 apex.f 
