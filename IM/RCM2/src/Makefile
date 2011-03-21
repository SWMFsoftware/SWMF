#^CFG COPYRIGHT UM
SHELL=/bin/sh

include ../Makefile.def
include ../Makefile.conf

include Makefile.DEPEND

CLAWDIR = claw

#
#  Makefile for RCM code:
#
#  Type "make" to compile and build executable
#
#

CLAWOBJECTS = \
  $(CLAWDIR)/qinit.o \
  $(CLAWDIR)/rpn2ad1.o \
  $(CLAWDIR)/rpt2ad1.o \
  $(CLAWDIR)/b4step2.o \
  $(CLAWDIR)/setaux.o \
  $(CLAWDIR)/psi.o \
  $(CLAWDIR)/claw2ez.o \
  $(CLAWDIR)/bc2.o \
  $(CLAWDIR)/claw2.o \
  $(CLAWDIR)/step2.o \
  $(CLAWDIR)/step2ds.o \
  $(CLAWDIR)/dimsp2.o \
  $(CLAWDIR)/flux2.o \
  $(CLAWDIR)/copyq2.o \
  $(CLAWDIR)/limiter.o \
  $(CLAWDIR)/philim.o \
  $(CLAWDIR)/src2.o 

RCMLIBRARY = \
  IM_wrapper.o

RCMOBJECTS = \
  ModProcIM.o\
  rcm_variables.o \
  rcm_io.o \
  rcm_routines.o \
  rcm_comput.o \
  gmresm.o \
  rcm_efield.o\
  impl_bicgstab.o \
  gamma.o\
  RCM_advec.o

OBJECTS = $(CLAWOBJECTS) $(RCMOBJECTS) $(RCMLIBRARY)

SEARCH_EXTRA = -I${LIBRARYDIR} -I${COUPLERDIR}

DEPEND:
	@perl ${SCRIPTDIR}/depend.pl ${SEARCH} ${OBJECTS}

MY_LIB = ${LIBDIR}/libIM.a

LIB:	DEPEND
	cd $(CLAWDIR); make lib
	make ${MY_LIB}
	@echo
	@echo ${MY_LIB} has been brought up to date.
	@echo

${MY_LIB}:	$(OBJECTS) ${MAKEFILE_COMP_SELECT}
		rm -f ${MY_LIB}
		${AR} ${MY_LIB} $(OBJECTS)

# RCM needs its own rules

.SUFFIXES:
.SUFFIXES: .f90 .F90 .f .o

.f90.o:
	${COMPILE.f90} ${CFLAGS} ${PRECISION} $<

.f.o:
	${COMPILE.f77} ${CFLAGS} ${PRECISION} $<

# IM_wrapper needs to be compiled with the same precision as others
IM_wrapper.o: IM_wrapper.f90 \
	ModProcIM.o \
	rcm_variables.o \
	rcm_io.o \
	${SHAREDIR}/ModReadParam.o \
	${SHAREDIR}/ModUtilities.o \
	${SHAREDIR}/ModNumConst.o \
	${LIBRARYDIR}/CON_comp_info.o \
	${COUPLERDIR}/CON_coupler.o \
	${SHAREDIR}/ModIoUnit.o
		${COMPILE.f90} $(CFLAGS) $(PRECISION) IM_wrapper.f90

distclean: clean
	-@(rm -f ${LIBDIR}/libRCM.a) 
	-@(rm -f ${BINDIR}/RCM.exe)

##############################################################################
# SWMF Makefile ends here, the rest is for standalone executables
##############################################################################

RCMEXTRACT = \
  rcm_variables.o \
  rcm_io.o \
  extract_arrays.o

RCMCONTROL = \
  control.o

RCM:
	cd $(CLAWDIR); make lib
	make RCMEXE

Extract: Extract.exe

all:	ALL

ALL:	RCM Extract


MPI_INC = ../Common mpif90.h

#
# stand-alone RCM executable
#
RCMEXE:		$(CLAWOBJECTS) $(RCMOBJECTS) $(RCMCONTROL)
		${LINK.f90} -o ${BINDIR}/RCM.exe $(Lmpi) $(CLAWOBJECTS) $(RCMOBJECTS) $(RCMCONTROL)

Extract.exe:	$(RCMEXTRACT)
		${COMPILE.f90} -o Extract.exe $(RCMEXTRACT)

setup_grid.out:		setup_grid.f90 $(RCMOBJECTS) $(CLAWOBJECTS)
			$(COMPILE.f90) $(Lmpi) setup_grid.f90 $(RCMOBJECTS) $(CLAWOBJECTS) -o setup_grid.out

setup_bf.out:		setup_bf.f90 $(RCMOBJECTS) $(CLAWOBJECTS)
			$(COMPILE.f90) $(Lmpi) setup_bf.f90 $(RCMOBJECTS) $(CLAWOBJECTS) -o setup_bf.out

setup_plasma.out:	setup_plasma.f90 $(RCMOBJECTS) $(CLAWOBJECTS)
			$(COMPILE.f90) $(Lmpi) setup_plasma.f90 $(RCMOBJECTS) $(CLAWOBJECTS) \
                        -o setup_plasma.out -lpgplot -lX11

read_mhd.out:		read_mhd.f90 $(RCMOBJECTS) $(CLAWOBJECTS)
			$(COMPILE.f90) $(Lmpi) read_mhd.f90 $(RCMOBJECTS) $(CLAWOBJECTS)\
                         -o read_mhd.out $(PGLIB) -lX11

plot_v.out:		plot_v.f90 $(RCMOBJECTS) $(CLAWOBJECTS)
			$(COMPILE.f90) $(Lmpi) plot_v.f90 $(RCMOBJECTS) $(CLAWOBJECTS)\
                         -o plot_v.out $(PGLIB) -lX11

plot_veff.out:		plot_veff.f90 $(RCMOBJECTS) $(CLAWOBJECTS)
			$(COMPILE.f90) plot_veff.f90 $(RCMOBJECTS) $(CLAWOBJECTS)\
                         -o plot_veff.out $(PGLIB) -lX11

plot_eeta.out:		plot_eeta.f90 $(RCMOBJECTS) $(CLAWOBJECTS)
			$(COMPILE.f90) $(Lmpi) plot_eeta.f90 $(RCMOBJECTS) $(CLAWOBJECTS)\
                         -o plot_eeta.out $(PGLIB) -lX11

plot_pvgamma.out:	plot_pvgamma.f90 $(RCMOBJECTS) $(CLAWOBJECTS)
			$(COMPILE.f90) $(Lmpi) plot_pvgamma.f90 $(RCMOBJECTS) $(CLAWOBJECTS)\
                         -o plot_pvgamma.out $(PGLIB) -lX11

plot_birk.out:		plot_birk.f90 $(CLAWOBJECTS) $(RCMOBJECTS)
			$(COMPILE.f90) $(Lmpi) plot_birk.f90 $(RCMOBJECTS) $(CLAWOBJECTS)\
                         -o plot_birk.out $(PGLIB) -lX11

plot_density.out:	plot_density.f90 $(CLAWOBJECTS) $(RCMOBJECTS)
			$(COMPILE.f90) $(Lmpi) plot_density.f90 $(RCMOBJECTS) $(CLAWOBJECTS)\
                         -o plot_density.out $(PGLIB) -lX11

plot_pressure.out:	plot_pressure.f90 $(CLAWOBJECTS) $(RCMOBJECTS)
			$(COMPILE.f90) $(Lmpi) plot_pressure.f90 $(RCMOBJECTS) $(CLAWOBJECTS)\
                         -o plot_pressure.out $(PGLIB) -lX11

plot_temperature.out:	plot_temperature.f90 $(CLAWOBJECTS) $(RCMOBJECTS)
			$(COMPILE.f90) $(Lmpi) plot_temperature.f90 $(RCMOBJECTS) $(CLAWOBJECTS)\
                         -o plot_temperature.out $(PGLIB) -lX11

plot_vm.out:		plot_vm.f90 $(CLAWOBJECTS) $(RCMOBJECTS)
			$(COMPILE.f90) $(Lmpi) plot_vm.f90 $(RCMOBJECTS) $(CLAWOBJECTS)\
                         -o plot_vm.out $(PGLIB) -lX11

### DO NOT remove this line - make depends on it ###

