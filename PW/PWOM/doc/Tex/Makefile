#^CMP COPYRIGHT UM
#BOP
#!ROUTINE: doc/Tex/Makefile - produce documentation from LaTex, XML, Protex
#!DESCRIPTION:
# This Makefile produces documentation for PWOM from Latex, PARAM.XML
# and Protex documented source code and scripts. 
#EOP
#BOC
include ../../Makefile.def
include ${DIR}/share/build/Makefile.doc

help:
	@echo "       Make manuals for PWOM"
	@echo ""
	@echo "Targets:"
	@echo ""
	@echo "PDF          - produces PDF manuals in the ../ directory"
	@echo ""
	@echo "clean        - remove temporary files"
	@echo ""
	@echo "cleanpdf     - clean + remove PDF manuals"
#EOC


PWOMXML:
	@rm -f PWOMXML.tex
	@touch PWOMXML.tex
	@perl -e 'for $$i (glob("../../PARAM.XML ../../[A-Z][A-Z]/*/PARAM.XML")){`${SCRIPTDIR}/XmlToTex.pl $$i >> PWOMXML.tex`}'
	@echo "${SCRIPTDIR}/XmlToTex.pl */PARAM.XML > PWOMXML.tex"

PWOM.dvi: \
	PWOM.tex \
	HEADER.tex \
	PWOM_introduction.tex \
	PWOM_configuration.tex \
	PWOM_param.tex \
	PWOMXML

PDF_FILES = PWOM.pdf
PDF: ${PDF_FILES}

clean:
	rm -f *.out *.aux *.log *.toc *.ind *.idx *.ilg *.dvi *.ps *~ *.protex
	rm -f PWOMXML.tex 
	@echo 'rm -f PWOMREF PWOMXML'
	@rm -f ${PWOMXML}

distclean: clean
	cd ..; rm -f ${PDF_FILES}
