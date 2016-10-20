include ../../Makefile.def

SEARCH_EXTRA = -I../ -I${DATAREADINDICESDIR} -I${EMPIRICALIEDIR} -I${EMPIRICALUADIR} \
	-I${LIBDIR}

include ../../Makefile.conf

help:
	@echo "Tests for STET modules: unit tests and combinations"
	@echo
	@echo "   unit_test_grid_potential"
	@echo "   unit_test_update_plasmasphere"
	@echo "   unit_test_update_plasmasphere_pot"
	@echo "   unit_test_update_iono"
	@echo "   unit_test_update_iono_pot"
	@echo "   unit_test_update_iono_nosources"
	@echo "   unit_test_update_iono_nosources_pot"
	@echo "   unit_test_background"
	@echo "   unit_test_res_scatter"
	@echo ""

test:   unit_test_update_plasmasphere 
	@echo ---------------------------------------------------------
	@ls -l *.diff

distclean: clean
	rm -rf *.tmp *.diff test*.out test_initial_state.ref \
		test_process_var_name.log

#######################################################################
OBJ_GRID_POT = \
	../ModMath.o\
	../ModGrid.o \
	../ModSeMpi.o \
	../ModSePlot.o 



unit_test_grid_potential_exe:
	@(cd .. ;make ModMath.o ModGrid.o ModSeMpi.o ModSePlot.o)
	@make unit_test_grid_potential.o 
	${LINK.f90} -o unit_test_grid_potential.exe ${OBJ_GRID_POT} \
		unit_test_grid_potential.o -L${LIBDIR}  -lEMPIRICALIE \
		-lEMPIRICALUA -lINDICES -lSHARE -lTIMING ${Lflag1}

unit_test_grid_potential: 
	@echo unit_test_grid_potential starting
	@make unit_test_grid_potential_exe 
	@echo ./unit_test_grid_potential_exe
	-@(./unit_test_grid_potential_exe )

#######################################################################

#######################################################################
OBJ_UPDATE_PLAS = \
	../ModMath.o\
	../ModGrid.o \
	../ModSeMpi.o \
	../ModSePlot.o \
	../ModSeState.o 


unit_test_update_plasmasphere_exe:
	@(cd .. ;make ModMath.o ModGrid.o ModSeMpi.o ModSePlot.o ModSeState.o)
	@make unit_test_update_plasmasphere.o 
	${LINK.f90} -o unit_test_update_plasmasphere.exe ${OBJ_UPDATE_PLAS} \
		unit_test_update_plasmasphere.o -L${LIBDIR}  -lEMPIRICALIE \
		-lEMPIRICALUA -lINDICES -lSHARE -lTIMING ${Lflag1}

unit_test_update_plasmasphere: 
	@echo unit_test_update_plasmasphere starting
	@make unit_test_update_plasmasphere_exe 
	@echo ./unit_test_update_plasmasphere_exe
	-@(./unit_test_update_plasmasphere_exe )

#######################################################################
#######################################################################
OBJ_UPDATE_PLAS_POT = \
	../ModMath.o\
	../ModGrid.o \
	../ModSeMpi.o \
	../ModSePlot.o \
	../ModSeState.o 


unit_test_update_plasmasphere_pot_exe:
	@(cd .. ;make ModMath.o ModGrid.o ModSeMpi.o ModSePlot.o ModSeState.o)
	@make unit_test_update_plasmasphere_pot.o 
	${LINK.f90} -o unit_test_update_plasmasphere_pot.exe \
		${OBJ_UPDATE_PLAS_POT} unit_test_update_plasmasphere_pot.o \
		-L${LIBDIR}  -lEMPIRICALIE -lEMPIRICALUA -lINDICES -lSHARE \
		-lTIMING ${Lflag1}

unit_test_update_plasmasphere_pot: 
	@echo unit_test_update_plasmasphere_pot starting
	@make unit_test_update_plasmasphere_pot_exe 
	@echo ./unit_test_update_plasmasphere_pot_exe
	-@(./unit_test_update_plasmasphere_pot_exe )

#######################################################################

#######################################################################
OBJ_UPDATE_IONO = \
	../ModMath.o\
	../ModGrid.o \
	../ModSeMpi.o \
	../ModSePlot.o \
	../ModSeCross.o \
	../ModSeProduction.o\
	../geomag.o\
	../iri90.o\
	../ModBackground.o \
	../ModSeState.o 


unit_test_update_iono_exe:
	@(cd .. ;make ModMath.o ModGrid.o ModSeMpi.o ModSePlot.o ModSeCross.o ModSeProduction.o geomag.o iri90.o ModBackground.o ModSeState.o)
	@make unit_test_update_iono.o 
	${LINK.f90} -o unit_test_update_iono.exe ${OBJ_UPDATE_IONO} \
		unit_test_update_iono.o -L${LIBDIR}  -lEMPIRICALIE \
		-lEMPIRICALUA -lINDICES -lSHARE -lTIMING ${Lflag1}

unit_test_update_iono: 
	@echo unit_test_update_iono starting
	@make unit_test_update_iono_exe 
	@echo ./unit_test_update_iono_exe
	-@(./unit_test_update_iono_exe )

#######################################################################

#######################################################################
OBJ_UPDATE_IONO_POT = \
	../ModMath.o\
	../ModGrid.o \
	../ModSeMpi.o \
	../ModSePlot.o \
	../ModSeCross.o \
	../ModSeProduction.o\
	../geomag.o\
	../iri90.o\
	../ModBackground.o \
	../ModSeState.o 


unit_test_update_iono_pot_exe:
	@(cd .. ;make ModMath.o ModGrid.o ModSeMpi.o ModSePlot.o ModSeCross.o ModSeProduction.o geomag.o iri90.o ModBackground.o ModSeState.o)
	@make unit_test_update_iono_pot.o 
	${LINK.f90} -o unit_test_update_iono_pot.exe ${OBJ_UPDATE_IONO_POT} \
		unit_test_update_iono_pot.o -L${LIBDIR}  -lEMPIRICALIE \
		-lEMPIRICALUA -lINDICES -lSHARE -lTIMING ${Lflag1}

unit_test_update_iono_pot: 
	@echo unit_test_update_iono_pot starting
	@make unit_test_update_iono_pot_exe 
	@echo ./unit_test_update_iono_pot_exe
	-@(./unit_test_update_iono_pot_exe )

#######################################################################
#######################################################################
OBJ_UPDATE_IONO_PLAS_POT = \
	../ModMath.o\
	../ModGrid.o \
	../ModSeMpi.o \
	../ModSePlot.o \
	../ModSeCross.o \
	../ModSeProduction.o\
	../geomag.o\
	../iri90.o\
	../ModBackground.o \
	../ModSeState.o 


unit_test_update_iono_plas_pot_exe:
	@(cd .. ;make ModMath.o ModGrid.o ModSeMpi.o ModSePlot.o ModSeCross.o ModSeProduction.o geomag.o iri90.o ModBackground.o ModSeState.o)
	@make unit_test_update_iono_plas_pot.o 
	${LINK.f90} -o unit_test_update_iono_plas_pot.exe ${OBJ_UPDATE_IONO_PLAS_POT} \
		unit_test_update_iono_plas_pot.o -L${LIBDIR}  -lEMPIRICALIE \
		-lEMPIRICALUA -lINDICES -lSHARE -lTIMING ${Lflag1}

unit_test_update_iono_plas_pot: 
	@echo unit_test_update_iono_plas_pot starting
	@make unit_test_update_iono_plas_pot_exe 
	@echo ./unit_test_update_iono_plas_pot_exe
	-@(./unit_test_update_iono_plas_pot_exe )

#######################################################################

#######################################################################
OBJ_UPDATE_IONO_NOSOURCES = \
	../ModMath.o\
	../ModGrid.o \
	../ModSeMpi.o \
	../ModSePlot.o \
	../ModSeState.o 


unit_test_update_iono_nosources_exe:
	@(cd .. ;make ModMath.o ModGrid.o ModSeMpi.o ModSePlot.o ModSeState.o)
	@make unit_test_update_iono_nosources.o 
	${LINK.f90} -o unit_test_update_iono_nosources.exe ${OBJ_UPDATE_IONO_NOSOURCES} \
		unit_test_update_iono_nosources.o -L${LIBDIR}  -lEMPIRICALIE \
		-lEMPIRICALUA -lINDICES -lSHARE -lTIMING ${Lflag1}

unit_test_update_iono_nosources: 
	@echo unit_test_update_iono_nosources starting
	@make unit_test_update_iono_nosources_exe 
	@echo ./unit_test_update_iono_nosources_exe
	-@(./unit_test_update_iono_nosources_exe )

#######################################################################
#######################################################################
OBJ_UPDATE_IONO_NOSOURCES_POT = \
	../ModMath.o\
	../ModGrid.o \
	../ModSeMpi.o \
	../ModSePlot.o \
	../ModSeState.o 


unit_test_update_iono_nosources_pot_exe:
	@(cd .. ;make ModMath.o ModGrid.o ModSeMpi.o ModSePlot.o ModSeState.o)
	@make unit_test_update_iono_nosources_pot.o 
	${LINK.f90} -o unit_test_update_iono_nosources_pot.exe ${OBJ_UPDATE_IONO_NOSOURCES_POT} \
		unit_test_update_iono_nosources_pot.o -L${LIBDIR}  -lEMPIRICALIE \
		-lEMPIRICALUA -lINDICES -lSHARE -lTIMING ${Lflag1}

unit_test_update_iono_nosources_pot: 
	@echo unit_test_update_iono_nosources_pot starting
	@make unit_test_update_iono_nosources_pot_exe 
	@echo ./unit_test_update_iono_nosources_pot_exe
	-@(./unit_test_update_iono_nosources_pot_exe )

#######################################################################
#######################################################################
OBJ_BACKGROUND = \
	../ModMath.o\
	../ModGrid.o \
	../ModSeMpi.o \
	../ModSeCross.o \
	../ModSeProduction.o\
	../geomag.o\
	../iri90.o\
	../ModBackground.o 


unit_test_background_exe:
	@(cd .. ;make ModMath.o ModGrid.o ModSeMpi.o ModSeCross.o ModSeProduction.o geomag.o iri90.o ModBackground.o)
	@make unit_test_background.o 
	${LINK.f90} -o unit_test_background.exe ${OBJ_BACKGROUND} \
		unit_test_background.o -L${LIBDIR}  -lEMPIRICALIE \
		-lEMPIRICALUA -lINDICES -lSHARE -lTIMING ${Lflag1}

unit_test_background: 
	@echo unit_test_background starting
	@make unit_test_background_exe 
	@echo ./unit_test_background_exe
	-@(./unit_test_background_exe )

#######################################################################
OBJ_RESSCATTER = \
	../ModMath.o\
	../ModGrid.o \
	../ModSeMpi.o \
	../ModSeProduction.o



unit_test_res_scatter_exe:
	@(cd .. ;make ModMath.o ModGrid.o ModSeMpi.o ModSeProduction.o)
	@make unit_test_res_scatter.o 
	${LINK.f90} -o unit_test_res_scatter.exe ${OBJ_RESSCATTER} \
		unit_test_res_scatter.o -L${LIBDIR}  -lEMPIRICALIE \
		-lEMPIRICALUA -lINDICES -lSHARE -lTIMING ${Lflag1}

unit_test_res_scatter: 
	@echo unit_test_res_scatter starting
	@make unit_test_res_scatter_exe 
	@echo ./unit_test_res_scatter_exe
	-@(./unit_test_background_exe )

#######################################################################
