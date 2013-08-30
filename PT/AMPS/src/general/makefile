SHELL =/bin/sh

DEFAULT_TARGET : general  

clean:
	rm -f *.o *.a

CC=gcc
Flags=-g   
IncludeList=

.SUFFIXES: .cpp .o 
.cpp.o: 
	${CC} -c ${Flags} ${IncludeList} $< 
 
GENERALOBJ = specfunc.o rnd.o global.o 

general : $(GENERALOBJ)
	ar -src general.a $(GENERALOBJ)
	@echo "Target 'general' is done." 
