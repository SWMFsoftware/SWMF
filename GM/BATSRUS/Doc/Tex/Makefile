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
	@#^CFG IF DOCHTML BEGIN
	@echo ""
	@echo "HTML         - produces online manuals in ../HTML directory"
	@#^CFG END DOCHTML
	@echo ""
	@echo "clean        - rm -rf *.aux *.log *.toc *.dvi *~"
	@echo ""
	@echo "cleanpdf     - clean + remove PDF manuals"
	@#^CFG IF DOCHTML BEGIN
	@echo ""
	@echo "cleanhtml    - clean + remove HTML manuals"
	@echo ""
	@echo "cleanall     - clean + remove PDF and HTML manuals"
	@#^CFG END DOCHTML

#^CFG IF DOCHTML BEGIN
.tex.html:
	mv $*.tex $*orig.tex
	perl -pe 's/\\tableofcontents/\%\\tableofcontents/' $*orig.tex > $*.tex
	latex $*
	latex $*
	latex2html -image_type gif -local_icons \
		  -address "(c) 2002 University of Michigan" \
		  -prev_url   "../../index.html"  \
		  -up_url     "../../index.html" \
		  -up_title   "BATSRUS Documents" \
		  -prev_title "BATSRUS Documents" \
		  $*
	cd $*; \
	perl -pi -e 's/<BODY >/<BODY bgcolor="#CECFCE">/' $*.html node*.html index.html; \
	perl -pi -e 's/^\[\]//' $*.html index.html; \
	perl -pi -e 's/SRC="prev.png"\>\<\/A\>/SRC="prev.png"\>\<\/A\>\n \
		\<A HREF="http:\/\/csem.engin.umich.edu\/"\>\n\<IMG WIDTH="100" \
		HEIGHT="24" ALIGN="BOTTOM" BORDER="0" ALT="CSEM Home" \
		SRC="csem_home.png"\>\<\/A\>/' $*.html node*.html index.html; \
	cp ../../HTML/csem_home.png ./ 
	rm -f $*.tex
	mv $*orig.tex $*.tex
	mv $* ../HTML
#^CFG END DOCHTML

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

USERMANUAL.html: ../../PARAM.xmltex

PDF:  USERMANUAL.pdf DESIGN.pdf STANDARDS.pdf

HTML: cleanhtml USERMANUAL.html DESIGN.html STANDARDS.html #^CFG IF DOCHTML

clean:
	rm -rf *.aux *.log *.toc *.dvi *~ ../../PARAM.xmltex

cleanpdf: clean
	cd ..;     rm -rf USERMANUAL.pdf DESIGN.pdf STANDARDS.pdf

#^CFG IF DOCHTML BEGIN
cleanhtml: clean
	cd ../HTML;rm -rf USERMANUAL     DESIGN     STANDARDS

cleanall: cleanpdf cleanhtml
#^CFG END DOCHTML

