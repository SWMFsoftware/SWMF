#^CFG COPYRIGHT UM
#
#	Specific Rules for Darwin fortran NAG compiler
#

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
