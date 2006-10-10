
PW3:
	cd src; make

UNSPLIT:
	cd src; make UNSPLIT

install: bin src/PLANET

bin:
	mkdir bin

include Makefile.planet

src/PLANET:
	mv Makefile.planet Makefile.planet.orig
	./config.pl -p=${PLANET}

test:	install run PW3
	cd run; pw3.exe
	diff -b run/restart/restart_iline0001.out  output/${PLANET}/restart_iline0001.out
	diff -b run/plots/plots_iline0001.out  output/${PLANET}/plots_iline0001.out
	diff -b run/log.out output/${PLANET}/log.out

testrun: PW3
	cd run; pw3.exe

rundir: run

run:
	mkdir run
	cd run; \
		ln -s ../bin/pw3.exe    .; \
		cp ../input/${PLANET}/pw.input pw.input; \
		mkdir restart;\
		cp ../input/${PLANET}/restart_iline0001.in restart/;\
		cp ../input/${PLANET}/North.dat .;\
		mkdir plots

clean:
	cd src; make clean

uninstall: 
	cd src; make distclean
	mv Makefile.planet.orig Makefile.planet
