include ../Makefile.def

SEARCHDIR = -I${SHAREDIR} -I${EMPIRICALGMDIR} -I${EMPIRICALIEDIR} -I${DATAREADINDICESDIR}

include ../Makefile.conf
include Makefile.DEPEND
include Makefile.RULES


OBJECTS= \
	ModSAMI.o\
	ModCoupleCimi.o\
	apexcord.o\
	grid3_mpi-1.98.o\
	hwm93.o\
	madala_sevp_dp.o\
	nrlmsise00.o\
	sami3_mpi-1.98.o\
	sevp13_dp.o\
	hwm14.o

###f90 = mpif90 -shared-intel -fp-model precise -mcmodel=large  -save -O3 -vec_report0

#f90 = mpif90_intel -shared-intel -fp-model precise -mcmodel=large  -save -O3 -vec_report0
#f90 = ompif90_intel -shared-intel -fp-model precise -mcmodel=large  -save -O3 -vec_report0

#f90 = mpif90 -mcmodel=large -save -C -traceback -shared-intel -vec_report0

#f90 = /opt/openmpi-1.4.3_intel/bin/mpif90 -mcmodel=medium -s -Rc -Rb -g -w

DEPEND:
	@${SCRIPTDIR}/depend.pl ${SEARCHDIR} ${OBJECTS}


MY_LIB = ${LIBDIR}/libSAMI3.a

LIB: DEPEND
	make ${MY_LIB}
	@echo
	@echo ${MY_LIB} has been brought up to date.
	@echo

${MY_LIB}: ${OBJECTS}
	rm -f ${MY_LIB}
	${AR} ${MY_LIB} ${OBJECTS} 

#Object files not included into the library
OBJECTS_EXE = sami3_main.o

LIBSHARE  = ${LIBDIR}/libSHARE.a
LIBTIMING = ${LIBDIR}/libTIMING.a

# Libraries should be compiled first, because modules are used in main
${OBJECTS_EXE}: ${LIBSHARE} ${MY_LIB}   



EXE = SAMI3.exe

SAMI3:
	@echo ' '
	@make ${EXE}
	@echo "Program ${EXE} has been brought up to date."
	@echo ' '

${EXE}: ${OBJECTS_EXE}
	rm -rf Tmp_; mkdir Tmp_
	cd Tmp_; \
		ar -x ${MY_LIB}; \
		ar -x ${LIBSHARE};  \
		ar -x ${LIBTIMING};  
	${LINK.f90} -o ${EXE} ${OBJECTS_EXE} Tmp_/*.o  ${Lflag1}
	rm -rf Tmp_

SAMI3_orig: ${OBJECTS_EXE}
	${COMPILE.f90} ${Lflag1} -o ${EXE} sami3_main.f90 ${OBJECTS} \
	-L${LIBDIR} -lTIMING

clean: cleanfiles

distclean: clean


#SAMI3 requires save flag  for origional .f files
# right now this only works for gfortran and would have to get 
# modified for other codes. Probably this hase to get fixed at Config.pl
include Makefile.suff


###.SUFFIXES : .o .f90 .f 
###
###
###%.o:%.mod
###
###.f.o:
###	$(f90) -c $*.f 
###
###.f90.o:
###	$(f90) -c $*.f90
###
###clean:
###	rm *.o *.exe 
###
###sami3.x:	$(OBJ)
###	$(f90) -o sami3.x $(OBJ)
###
###$(OBJ): com3_mpi-1.98.inc
###$(OBJ): param3_mpi-1.98.inc

