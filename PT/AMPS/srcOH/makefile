SHELL =/bin/sh

DEFAULT_TARGET : lib

include ../../Makefile.conf

MAINLIBOBJ = OH.o  main_lib.o OH_Coupler.o OH_sample_distribution_function.o
MAINOBJ = main.o

lib : $(MAINLIBOBJ)
	ar -rc mainlib.a $(MAINLIBOBJ)
	@echo "Target 'main' is done."

amps : $(MAINOBJ)
	ar -rc main.a $(MAINOBJ)
