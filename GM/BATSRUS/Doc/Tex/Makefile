#^CFG COPYRIGHT UM
DIR	= `pwd`

.SUFFIXES:
.SUFFIXES: .tex .dvi .ps .pdf .html

.tex.dvi:
	latex $<
	latex $<
	latex $<

.dvi.ps: 
	dvips -o $*.ps $<

.ps.pdf:
	ps2pdf $*.ps ../$*.pdf

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

PDF:  USERMANUAL.pdf DESIGN.pdf STANDARDS.pdf

HTML: cleanhtml USERMANUAL.html DESIGN.html STANDARDS.html #^CFG IF DOCHTML

clean:
	rm -rf *.aux *.log *.toc *.dvi *~

cleanpdf: clean
	cd ..;     rm -rf USERMANUAL.pdf DESIGN.pdf STANDARDS.pdf

#^CFG IF DOCHTML BEGIN
cleanhtml: clean
	cd ../HTML;rm -rf USERMANUAL     DESIGN     STANDARDS

cleanall: cleanpdf cleanhtml
#^CFG END DOCHTML

