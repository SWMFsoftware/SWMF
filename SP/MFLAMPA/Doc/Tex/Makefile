#  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
#  For more information, see http://csem.engin.umich.edu/tools/swmf
include ../../Makefile.def
include ${DIR}/share/build/Makefile.doc

help:
	@echo "       Make manuals for MFLAMPA"
	@echo ""
	@echo "Targets:"
	@echo ""
	@echo "PDF          - produces USERMANUAL.pdf"
	@echo ""
	@echo "clean        - rm -rf *.aux *.log *.toc *.dvi *~"
	@echo ""
	@echo "cleanpdf     - clean + remove PDF manuals"

USERMANUAL.dvi: \
	USERMANUAL.tex \
	../../PARAM.xmltex 

PDF:  USERMANUAL.pdf

clean:
	rm -rf *.aux *.log *.toc *.dvi *~ ../../PARAM.xmltex

cleanpdf: clean
	cd ..;     rm -rf USERMANUAL.pdf 


