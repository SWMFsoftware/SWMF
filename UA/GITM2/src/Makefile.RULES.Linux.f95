#^CFG COPYRIGHT UM
#
#	Specific Rules for NAG compiler on grendel (Linux)
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

