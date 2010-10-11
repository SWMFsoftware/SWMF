
SHELL =/bin/sh

include ../Makefile.def

SEARCH_EXTRA = -I../src

include ../Makefile.conf

########################################################################

OBJ_BLOCK_DATA = \
	test_block_data.o \
	../src/ModBlockData.o

test_block_data.exe: ${OBJ_BLOCK_DATA}
	${LINK.f90} -o test_block_data.exe ${OBJ_BLOCK_DATA} ${Lflag2}

test_block_data: 
	@echo -----------------------------------------------------------------
	@echo test_block_data starting
	@echo
	@(cd ../src; make ModBlockData.o)
	@make test_block_data.exe
	@echo
	@echo './test_block_data.exe'
	@./test_block_data.exe
	@echo
	@echo test_block_data done
	@echo -----------------------------------------------------------------

distclean: clean
