# Compiling lapack.f and/or blas.f with -O2 causes run time error. 
# Both -O4 and -O1 seem to work.
lapack.o: lapack.f
	${COMPILE.f77} ${Cflag4} -132 lapack.f

blas.o: blas.f
	${COMPILE.f77} ${Cflag4} -132 blas.f
