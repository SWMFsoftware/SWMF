# defines OS
-include ../Makefile.def

install:
	cd Control/src; cp show_git_info.h.orig show_git_info.h

clean:
	cd Control/src;   make clean
	cd Library/src;   make clean
	cd Coupler/src;   make clean
	cd Interface/src; make clean
	cd Stubs/src;     make clean

distclean: clean
