#^CFG COPYRIGHT UM

include ../../../Makefile.def
include ../../../Makefile.conf

../src/IONIZ.exe:
	cd ../src/;make IONIZ

Table1.tex: ../src/IONIZ.exe
	cd ../src;IONIZ.exe
ionization_equilibrium.tex: Table1.tex Table2.tex

ionization_equilibrium.pdf: ionization_equilibrium.tex
	latex ionization_equilibrium.tex
	latex ionization_equilibrium.tex
	dvipdf ionization_equilibrium

PDF: ionization_equilibrium.pdf

distclean: clean
	rm -f *.pdf *.aux *.log *.dvi Table* ../src/IONIZ.exe


