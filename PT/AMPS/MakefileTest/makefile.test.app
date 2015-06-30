test_<APP>:
	@($(MAKE) test_<APP>_compile)
	@($(MAKE) test_<APP>_rundir)
	@($(MAKE) test_<APP>_run)
	@($(MAKE) test_<APP>_check)

test_<APP>_compile:
	@echo "test_<APP>_compile..." > test_<APP>.diff
	./Config.pl -application=<APPPATH> -spice-path=nospice -spice-kernels=nospice -model-data-path=$(MYDIR)/data/input/<APP>  -amps-test=on $(TEST<APP>KEYS)
	rm -rf srcTemp
	./ampsConfig.pl -no-compile
	-@($(MAKE) amps)

test_<APP>_rundir:
	@echo "test_<APP>_rundir..." >> test_<APP>.diff
	rm -rf   $(TEST<APP>DIR)
	mkdir -p $(TEST<APP>DIR)
	mv amps  $(TEST<APP>DIR)

test_<APP>_run:
	@echo "test_<APP>_run..." >> test_<APP>.diff
	cd $(TEST<APP>DIR); ${MPIRUN} ./amps

test_<APP>_check:
	@echo "test_<APP>_check..." >> test_<APP>.diff
	-@$(foreach OUT,$(TEST<APP>OUTFILES),                                 \
	-(${SCRIPTDIR}/DiffNum.pl $(TEST<APP>DIR)/PT/plots/$(OUT).dat         \
	output/$(OUT).ref_np`ls $(TEST<APP>DIR)/PT/thread* |wc -l |tr -d ' '` \
	> test_<APP>.diff);)
	@ls -l test_<APP>.diff


