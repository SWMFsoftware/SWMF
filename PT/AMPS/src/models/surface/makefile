SHELL =/bin/sh

DEFAULT_TARGET : surface

include ../../../Makefile.conf

SURFACEGOBJ = Surface_particle_interaction.o Surface_temperature.o Surface.o

surface : $(SURFACEGOBJ)
	ar -src surface.a $(SURFACEGOBJ)
	@echo "Target 'surface' is done." 
