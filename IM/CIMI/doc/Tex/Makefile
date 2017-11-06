# This Makefile produces documentation for CIMI from Latex, PARAM.XML
# and Protex documented source code and scripts. 
#EOP
#BOC
include ../../Makefile.def
include ${DIR}/share/build/Makefile.doc

help:
	@echo "       Make manuals for CIMI"
	@echo ""
	@echo "Targets:"
	@echo ""
	@echo "PDF          - produces PDF manuals in the ../ directory"
	@echo ""
	@echo "clean        - remove temporary files"
	@echo ""
	@echo "cleanpdf     - clean + remove PDF manuals"
#EOC


CIMIXML:
	@rm -f CIMIXML.tex
	@touch CIMIXML.tex
	@perl -e 'for $$i (glob("../../PARAM.XML ../../[A-Z][A-Z]/*/PARAM.XML")){`${SCRIPTDIR}/XmlToTex.pl $$i >> CIMIXML.tex`}'
	@echo "${SCRIPTDIR}/XmlToTex.pl */PARAM.XML > CIMIXML.tex"

CIMI.dvi: \
	CIMI.tex \
	HEADER.tex \
	CIMI_introduction.tex \
	CIMI_configuration.tex \
	CIMI_param.tex \
	CIMIXML

PDF_FILES = CIMI.pdf
PDF: ${PDF_FILES}

clean:
	rm -f *.out *.aux *.log *.toc *.ind *.idx *.ilg *.dvi *.ps *~ *.protex
	rm -f CIMIXML.tex 
	@echo 'rm -f CIMIREF CIMIXML'
	@rm -f ${CIMIXML}

distclean: clean
	cd ..; rm -f ${PDF_FILES}
