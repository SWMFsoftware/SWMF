#^CFG COPYRIGHT UM

include ../../../../Makefile.def
include ${DIR}/share/build/Makefile.doc

# Note that the ionization_table*.tex files were produced
# by the util/CRASH/src/IONIZ.exe code

PDF: ionization_equilibrium.pdf GreyDiffusion.pdf EvidenceFile.pdf NTemperature.pdf FermiStatistics.pdf Coulomb.pdf

clean:
	rm -f *.log *.aux *~

distclean: clean
	rm -f ../*.pdf *.dvi *.ps
