default : CIMI

include Makefile.def

INSTALLFILES =  src/Makefile.DEPEND \
		src/Makefile.RULES \
		srcSAMI3/Makefile.DEPEND \
		srcSAMI3/Makefile.RULES \
		srcInterface/Makefile.DEPEND


install: 
	touch ${INSTALLFILES}
	./Config.pl -EarthHO -GridDefault
	@(if [ ! -d input ];  then ln -s data/input  input;  fi)
	@(if [ ! -d output ]; then ln -s data/output output; fi)

#
#       General Housekeeping
#

CIMI:
	@cd ${SHAREDIR};  	make LIB
	@cd ${NOMPIDIR};	make LIB
	@cd ${TIMINGDIR}; 	make LIB 
	@cd ${EMPIRICALIEDIR};	make LIB
	@cd ${EMPIRICALGMDIR};	make LIB
	@cd ${DATAREADINDICESDIR};make LIB
	@cd src;	make LIB
	@cd src;	make CIMI

CIMI_SAMI:
	@cd ${SHAREDIR};  	make LIB
	@cd ${NOMPIDIR};	make LIB
	@cd ${TIMINGDIR}; 	make LIB 
	@cd ${EMPIRICALIEDIR};	make LIB
	@cd ${EMPIRICALGMDIR};	make LIB
	@cd ${DATAREADINDICESDIR};make LIB
	@cd srcSAMI3;	make LIB
	@cd src;	make LIB
	@cd src;	make CIMI_SAMI


SAMI3:
	@cd ${SHAREDIR};  	make LIB
	@cd ${TIMINGDIR}; 	make LIB 
	@cd srcSAMI3;	make LIB
	@cd srcSAMI3;	make SAMI3

NOMPI:
	cd util/NOMPI/src; $(MAKE) LIB

LIB:
	cd src; make LIB
	cd srcInterface; make LIB

TESTDIR = run_test

test:
	@echo "test_compile..." > test_cimi.diff
	make   test_compile
	@echo "test_rundir..." >> test_cimi.diff
	make   test_rundir
	@echo "test_run..."    >> test_cimi.diff
	make   test_run
	@echo "test_check..."  >> test_cimi.diff
	make   test_check

test_all:
	@echo "test_compile..." > test_cimi.diff
	make   test_compile
	@echo "test_rundir..." >> test_cimi.diff
	make   test_rundir
	@echo "test_run..."    >> test_cimi.diff
	make   test_run
	@echo "test_check..."  >> test_cimi.diff
	make   test_check_all
	@echo "test_rundir_WAVES..." >> test_cimi.diff
	make   test_rundir_WAVES
	@echo "test_run..."    >> test_cimi.diff
	make   test_run
	@echo "test_check_WAVES..."  >> test_cimi.diff
	make   test_check_WAVES
	@echo "test_rundir_dipole..." >> test_cimi.diff
	make   test_rundir_dipole
	@echo "test_run..."    >> test_cimi.diff
	make   test_run
	@echo "test_check_dipole..."  >> test_cimi.diff
	make   test_check_dipole
	@echo "test_rundir_Prerun..." >> test_cimi.diff
	make   test_rundir_Prerun
	@echo "test_run..."    >> test_cimi.diff
	make   test_run
	@echo "test_check_Prerun..."  >> test_cimi.diff
	make   test_check_Prerun
	ls -l test_*.diff

test_WAVES:
	@echo "test_compile..." > test_cimi.diff
	make   test_compile
	@echo "test_rundir_WAVES..." >> test_cimi.diff
	make   test_rundir_WAVES
	@echo "test_run..."    >> test_cimi.diff
	make   test_run
	@echo "test_check_WAVES..."  >> test_cimi.diff
	make   test_check_WAVES
	ls -l test_*.diff

test_dipole:
	@echo "test_compile..." > test_cimi.diff
	make   test_compile
	@echo "test_rundir_dipole..." >> test_cimi.diff
	make   test_rundir_dipole
	@echo "test_run..."    >> test_cimi.diff
	make   test_run
	@echo "test_check_dipole..."  >> test_cimi.diff
	make   test_check_dipole
	ls -l test_*.diff

test_flux:
	@echo "test_compile..." > test_cimi.diff
	make   test_compile
	@echo "test_rundir_flux..." >> test_cimi.diff
	make   test_rundir_flux
	@echo "test_run..."    >> test_cimi.diff
	make   test_run
	@echo "test_check_flux..."  >> test_cimi.diff
	make   test_check_flux
	@echo "test_check_eq..."  >> test_cimi.diff
	make   test_check_eq
	ls -l test_*.diff

test_drift:
	@echo "test_compile..." > test_cimi.diff
	make   test_compile
	@echo "test_rundir_drift..." >> test_cimi.diff
	make   test_rundir_drift
	@echo "test_run..."    >> test_cimi.diff
	make   test_run
	@echo "test_check_eq..."  >> test_cimi.diff
	make   test_check_eq
	@echo "test_check_drift..."  >> test_cimi.diff
	make   test_check_drift
	ls -l test_*.diff

test_Prerun:
	@echo "test_compile..." > test_cimi.diff
	make   test_compile
	@echo "test_rundir_Prerun..." >> test_cimi.diff
	make   test_rundir_Prerun
	@echo "test_run..."    >> test_cimi.diff
	make   test_run
	@echo "test_check_Prerun..."  >> test_cimi.diff
	make   test_check_Prerun
	ls -l test_*.diff

test_compile:
	./Config.pl -EarthHO -GridDefault -show
	make CIMI

test_rundir:
	rm -rf ${TESTDIR}
	make rundir RUNDIR=${TESTDIR} STANDALONE=YES IMDIR=`pwd`
	cp input/testfiles/*.dat ${TESTDIR}/
	cp input/testfiles/PARAM.in.test.NOWAVES ${TESTDIR}/PARAM.in

test_rundir_WAVES:
	rm -rf ${TESTDIR}
	make rundir RUNDIR=${TESTDIR} STANDALONE=YES IMDIR=`pwd`
	cp input/testfiles/*.dat ${TESTDIR}/
	cp input/testfiles/PARAM.in.test.WAVES ${TESTDIR}/PARAM.in

test_rundir_flux:
	rm -rf ${TESTDIR}
	make rundir RUNDIR=${TESTDIR} STANDALONE=YES IMDIR=`pwd`
	cp input/testfiles/*.dat ${TESTDIR}/
	cp input/testfiles/PARAM.in.test.flux ${TESTDIR}/PARAM.in

test_rundir_drift:
	rm -rf ${TESTDIR}
	make rundir RUNDIR=${TESTDIR} STANDALONE=YES IMDIR=`pwd`
	cp input/testfiles/*.dat ${TESTDIR}/
	cp input/testfiles/PARAM.in.test.drift ${TESTDIR}/PARAM.in

test_rundir_dipole:
	rm -rf ${TESTDIR}
	make rundir RUNDIR=${TESTDIR} STANDALONE=YES IMDIR=`pwd`
	cp input/testfiles/*.dat ${TESTDIR}/
	cp input/testfiles/PARAM.in.test.dipole ${TESTDIR}/PARAM.in

test_rundir_Prerun:
	rm -rf ${TESTDIR}
	make rundir RUNDIR=${TESTDIR} STANDALONE=YES IMDIR=`pwd`
	cp input/testfiles/imf.dat.Prerun ${TESTDIR}/imf.dat
	cp input/testfiles/Indices.dat.Prerun ${TESTDIR}/Indices.dat
	cp input/testfiles/Prerun/* ${TESTDIR}/IM/.
	cp input/testfiles/PARAM.in.test.Prerun ${TESTDIR}/PARAM.in

test_run:
	cd ${TESTDIR}; ${MPIRUN} ./cimi.exe > runlog 

# reduced set of checks for the SWMF nightly tests
test_check:
	-make test_check_eq
	ls -l test_cimi*.diff

# complete set of checks
test_check_all:
	-make test_check_flux
	-make test_check_psd
	-make test_check_eq
	-make test_check_drift
	-make test_check_log
	ls -l test_cimi*.diff

test_check_WAVES:
	-${SCRIPTDIR}/DiffNum.pl -r=0.001 -a=1e-10 \
		${TESTDIR}/IM/plots/CimiFlux_e.fls \
		output/CimiFlux_e.fls.WAVES \
		> test_cimi_waves.diff
	-${SCRIPTDIR}/DiffNum.pl -r=0.001 -a=1e-10 \
		${TESTDIR}/IM/plots/CimiPSD_e.psd \
		output/CimiPSD_e.psd.WAVES \
		>> test_cimi_waves.diff
	-${SCRIPTDIR}/DiffNum.pl -r=0.001 -a=1e-10 \
		${TESTDIR}/IM/plots/CIMIeq.outs \
		output/CIMIeq.outs.WAVES \
		>> test_cimi_waves.diff
	-${SCRIPTDIR}/DiffNum.pl -r=0.001 -a=1e-10 \
		${TESTDIR}/IM/plots/CIMI.log \
		output/CIMI.log.WAVES \
		>> test_cimi_waves.diff

test_check_flux:
	-${SCRIPTDIR}/DiffNum.pl -r=0.001 -a=1e-10 \
		${TESTDIR}/IM/plots/CimiFlux_h.fls \
		output/CimiFlux_h.fls \
		> test_cimi_flux.diff
	-${SCRIPTDIR}/DiffNum.pl -r=0.001 -a=1e-10 \
		${TESTDIR}/IM/plots/CimiFlux_o.fls \
		output/CimiFlux_o.fls \
		>> test_cimi_flux.diff
	-${SCRIPTDIR}/DiffNum.pl -r=0.001 -a=1e-10 \
		${TESTDIR}/IM/plots/CimiFlux_e.fls \
		output/CimiFlux_e.fls \
		>> test_cimi_flux.diff

test_check_psd:
	-${SCRIPTDIR}/DiffNum.pl -t -r=0.001 -a=1e-10 \
		${TESTDIR}/IM/plots/CimiPSD_h.psd \
		output/CimiPSD_h.psd \
		> test_cimi_psd.diff
	-${SCRIPTDIR}/DiffNum.pl -t -r=0.001 -a=1e-10 \
		${TESTDIR}/IM/plots/CimiPSD_o.psd \
		output/CimiPSD_o.psd \
		>> test_cimi_psd.diff
	-${SCRIPTDIR}/DiffNum.pl -t -r=0.001 -a=1e-10 \
		${TESTDIR}/IM/plots/CimiPSD_e.psd \
		output/CimiPSD_e.psd \
		>> test_cimi_psd.diff

test_check_drift:
	-${SCRIPTDIR}/DiffNum.pl -r=0.001 -a=1e-10 \
		${TESTDIR}/IM/plots/CimiDrift_h.vp \
		output/CimiDrift_h.vp \
		> test_cimi_drift.diff
	-${SCRIPTDIR}/DiffNum.pl -r=0.001 -a=1e-10 \
		${TESTDIR}/IM/plots/CimiDrift_o.vp \
		output/CimiDrift_o.vp \
		>> test_cimi_drift.diff
	-${SCRIPTDIR}/DiffNum.pl -r=0.001 -a=1e-10 \
		${TESTDIR}/IM/plots/CimiDrift_e.vp \
		output/CimiDrift_e.vp \
		>> test_cimi_drift.diff

test_check_dipole:
	-${SCRIPTDIR}/DiffNum.pl -r=0.001 -a=1e-10 \
		${TESTDIR}/IM/plots/CimiDrift_h.vp \
		output/CimiDrift_h.vp.dipole \
		> test_cimi_dipole.diff
	-${SCRIPTDIR}/DiffNum.pl -r=0.001 -a=1e-10 \
		${TESTDIR}/IM/plots/CimiDrift_o.vp \
		output/CimiDrift_o.vp.dipole \
		>> test_cimi_dipole.diff
	-${SCRIPTDIR}/DiffNum.pl -r=0.001 -a=1e-10 \
		${TESTDIR}/IM/plots/CimiDrift_e.vp \
		output/CimiDrift_e.vp.dipole \
		>> test_cimi_dipole.diff
	-${SCRIPTDIR}/DiffNum.pl -r=0.001 -a=1e-10 \
		${TESTDIR}/IM/plots/CIMI.log \
		output/CIMI.log.dipole \
		>> test_cimi_dipole.diff

test_check_eq:
	-${SCRIPTDIR}/DiffNum.pl -r=0.001 -a=1e-10 \
		${TESTDIR}/IM/plots/CIMIeq.outs \
		output/CIMIeq.outs \
		> test_cimi.diff

test_check_log:
	-${SCRIPTDIR}/DiffNum.pl -r=0.001 -a=1e-10 \
		${TESTDIR}/IM/plots/CIMI.log \
		output/CIMI.log \
		> test_cimi_log.diff

test_check_Prerun:
	-${SCRIPTDIR}/DiffNum.pl -r=0.001 -a=1e-10 \
		${TESTDIR}/IM/plots/CimiFlux_h.fls \
		output/CimiFlux_h.fls.Prerun \
		> test_cimi_Prerun.diff
	-${SCRIPTDIR}/DiffNum.pl -r=0.001 -a=1e-10 \
		${TESTDIR}/IM/plots/CimiFlux_o.fls \
		output/CimiFlux_o.fls.Prerun \
		>> test_cimi_Prerun.diff
	-${SCRIPTDIR}/DiffNum.pl -r=0.001 -a=1e-10 \
		${TESTDIR}/IM/plots/CimiFlux_e.fls \
		output/CimiFlux_e.fls.Prerun \
		>> test_cimi_Prerun.diff
	-${SCRIPTDIR}/DiffNum.pl -r=0.001 -a=1e-10 \
		${TESTDIR}/IM/plots/CimiPSD_h.psd \
		output/CimiPSD_h.psd.Prerun \
		>> test_cimi_Prerun.diff
	-${SCRIPTDIR}/DiffNum.pl -r=0.001 -a=1e-10 \
		${TESTDIR}/IM/plots/CimiPSD_o.psd \
		output/CimiPSD_o.psd.Prerun \
		>> test_cimi_Prerun.diff
	-${SCRIPTDIR}/DiffNum.pl -r=0.001 -a=1e-10 \
		${TESTDIR}/IM/plots/CimiPSD_e.psd \
		output/CimiPSD_e.psd.Prerun \
		>> test_cimi_Prerun.diff
	-${SCRIPTDIR}/DiffNum.pl -r=0.001 -a=1e-10 \
		${TESTDIR}/IM/plots/CIMIeq.outs \
		output/CIMIeq.outs.Prerun \
		>> test_cimi_Prerun.diff
	-${SCRIPTDIR}/DiffNum.pl -r=0.001 -a=1e-10 \
		${TESTDIR}/IM/plots/CIMI.log \
		output/CIMI.log.Prerun \
		>> test_cimi_Prerun.diff

PDF:
	@cd doc/Tex; make PDF

clean:
	touch ${INSTALLFILES}
	cd src; make clean
	cd srcSAMI3; make clean
	cd srcInterface; make clean
	cd doc/Tex; make clean
	(if [ -d util ];  then cd util;  make clean; fi);
	(if [ -d share ]; then cd share; make clean; fi);

distclean:
	./Config.pl -uninstall

allclean:
	@touch ${INSTALLFILES}
	cd src; make distclean
	cd srcInterface; make distclean
	rm -f config.log *~
	if [ -h input ]; then rm -f input; fi
	if [ -h output ]; then rm -f output; fi

#
#       Create run directories
#
rundir:
	mkdir -p ${RUNDIR}/IM
	@(cd ${RUNDIR}; \
		if [ ! -e "EIE/README" ]; then \
			ln -s ${EMPIRICALIEDIR}/data EIE;\
		fi;)
	cd ${RUNDIR}/IM; \
		cp ${IMDIR}/input/quiet*fin . ;\
		cp ${IMDIR}/input/WaveData/*dat . ;\
		mkdir plots restartIN restartOUT
	@(if [ "$(STANDALONE)" != "NO" ]; then \
		cd ${RUNDIR} ; \
		ln -s ${BINDIR}/cimi.exe .   ; \
		touch core ; chmod 444 core;\
	fi);

rundir_cimi_sami:
	mkdir -p ${RUNDIR}/IM
	@(cd ${RUNDIR}; \
		if [ ! -e "EIE/README" ]; then \
			ln -s ${EMPIRICALIEDIR}/data EIE;\
		fi;)
	cd ${RUNDIR}/IM; \
		cp ${IMDIR}/input/quiet*fin . ;\
		cp ${IMDIR}/input/WaveData/*dat . ;\
		mkdir plots restartIN restartOUT plotsSAMI3
	@(if [ "$(STANDALONE)" != "NO" ]; then \
		cp input/testfiles/*.dat ${RUNDIR}/ ;\
		cp input/testfiles/PARAM.in.test.WAVES ${RUNDIR}/;\
		cp srcSAMI3/sami3_mpi-1.98.namelist ${RUNDIR}/ ;\
		cp srcSAMI3/*.inp ${RUNDIR}/ ;\
		cd ${RUNDIR} ; \
		ln -s ${BINDIR}/cimi_sami.exe .   ; \
		touch core ; chmod 444 core;\
	fi);


rundir_sami:
	mkdir -p ${RUNDIR}/IM
	cd ${RUNDIR}/IM; \
		mkdir restartIN restartOUT plotsSAMI3
	@(if [ "$(STANDALONE)" != "NO" ]; then \
		cp srcSAMI3/sami3_mpi-1.98.namelist ${RUNDIR}/ ;\
		cp srcSAMI3/*.inp ${RUNDIR}/ ;\
		ln -s srcSAMI3/sami3.exe ${RUNDIR}/SAMI3.exe ;\
		cd ${RUNDIR} ; \
		touch core ; chmod 444 core;\
	fi);
