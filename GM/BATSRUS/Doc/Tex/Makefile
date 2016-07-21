#  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
#  For more information, see http://csem.engin.umich.edu/tools/swmf
include ../../Makefile.def
include ${DIR}/share/build/Makefile.doc

help:
	@echo "       Make manuals for BATSRUS"
	@echo ""
	@echo "Targets:"
	@echo ""
	@echo "PDF          - produces USERMANUAL.pdf DESIGN.pdf STANDARDS.pdf"
	@echo ""
	@echo "clean        - rm -rf *.aux *.log *.toc *.dvi *~"
	@echo ""
	@echo "cleanpdf     - clean + remove PDF manuals"

USERMANUAL.dvi: \
	USERMANUAL.tex \
	HEADER.tex \
	introduction.tex \
	quickstart.tex \
	hardware_software.tex \
	installation.tex \
	runs.tex \
	param.tex \
	../../PARAM.xmltex \
	output.tex \
	idl.tex \
	tecplot.tex

PDF:  USERMANUAL.pdf DESIGN.pdf STANDARDS.pdf

clean:
	rm -rf *.aux *.log *.toc *.dvi *~ ../../PARAM.xmltex

cleanpdf: clean
	cd ..;     rm -rf USERMANUAL.pdf DESIGN.pdf STANDARDS.pdf


