SHELL =/bin/sh

DEFAULT_TARGET : species  

clean:
	rm -f *.o *.a

CC=gcc
Flags=-g   
IncludeList=

.SUFFIXES: .cpp .o 
.cpp.o: 
	${CC} -c ${Flags} ${IncludeList} $< 


SPECIESOBJ = Na.o He.o 
species : $(SPECIESOBJ)
	ar -src species.a $(SPECIESOBJ)
	@echo "Target 'species' is done." 
