# Makefile for AMPS stand-alone nightly test for the Moon application

TESTDIRMOON = run_test_moon
WSDMOON     = srcTestMoon

### Moon application general test ###
		
test_moon:
	$(MAKE) test_moon_compile
	$(MAKE) test_moon_rundir
	$(MAKE) test_moon_run
	$(MAKE) test_moon_check

test_moon_compile:
	@echo "test_moon_compile..." > test_moon.diff
	./Config.pl -application=Moon
	rm -rf srcTestMoon
	./ampsConfig.pl -no-compile 
	$(MAKE) amps

test_moon_rundir:
	@echo "test_moon_rundir..." >> test_moon.diff
	rm -rf   ${TESTDIRMOON}
	mkdir -p ${TESTDIRMOON}
	mv amps  ${TESTDIRMOON}

test_moon_run:
	@echo "test_moon_run..." >> test_moon.diff
	cd ${TESTDIRMOON}; ${MPIRUN} ./amps

test_moon_check:
	@echo "test_moon_check..." >> test_moon.diff
	-(${SCRIPTDIR}/DiffNum.pl ${TESTDIRMOON}/PT/plots/amps.dat \
	output/test_amps.ref_np`ls ${TESTDIRMOON}/PT/thread* | wc -l | tr -d ' '` \
	> test_moon.diff)
	@ls -l test_moon.diff


