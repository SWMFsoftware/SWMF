SHELL =/bin/sh

DEFAULT_TARGET : mesh  

clean:
	rm -f *.o *.a

CC=gcc
Flags=-g   
IncludeList=

.SUFFIXES: .cpp .o 
.cpp.o: 
	${CC} -c ${Flags} ${IncludeList} $< 


MESHOBJ = meshAMRgeneric.o  meshAMRcutcell.o
mesh : $(MESHOBJ)
	ar -src mesh.a $(MESHOBJ)
	@echo "Target 'mesh' is done." 
