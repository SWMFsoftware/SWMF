#^CFG COPYRIGHT UM
include ../../Makefile.def

SWMFREF = \
	${CONTROLDIR}/CON_main.protex \
	${CONTROLDIR}/CON_variables.protex \
	${LIBRARYDIR}/CON_world.protex \
	${LIBRARYDIR}/CON_comp_param.protex \
	${LIBRARYDIR}/CON_comp_info.protex \
	${CONTROLDIR}/CON_io.protex \
	${SHAREDIR}/ModReadParam.protex \
	${SHAREDIR}/ModIoUnit.protex \
	${CONTROLDIR}/CON_session.protex \
	${INTERFACEDIR}/CON_wrapper.protex \
	${INTERFACEDIR}/CON_couple_all.protex \
	${INTERFACEDIR}/CON_couple_ih_gm.protex \
	${INTERFACEDIR}/CON_couple_gm_ie_swmf.protex \
	${INTERFACEDIR}/CON_couple_gm_im.protex \
	${INTERFACEDIR}/CON_couple_ie_im.protex \
	${INTERFACEDIR}/CON_couple_ie_ua.protex \
	${COUPLERDIR}/CON_domain_decomposition.protex \
	${COUPLERDIR}/CON_grid_descriptor.protex \
	${COUPLERDIR}/CON_router.protex \
	${COUPLERDIR}/CON_global_message_pass.protex \
	${SHAREDIR}/CON_axes.protex \
	${SHAREDIR}/ModCoordTransform.protex \
	${SHAREDIR}/CON_time.protex \
	${SHAREDIR}/ModTimeType.protex \
	${SHAREDIR}/ModKind.protex \
	${SHAREDIR}/CON_planet.protex \
	${SHAREDIR}/CON_planet_field.protex \
	${SHAREDIR}/CON_physics.protex \
	${CONTROLDIR}/CON_methods.protex \
	${SHAREDIR}/ModMpi.protex \
	${SCRIPTDIR}/depend.protex \
	${SCRIPTDIR}/Rename.protex \
	${SCRIPTDIR}/Methods.protex \
	${GMDIR}/TestBatsrus.protex \
	${GMDIR}/TestSuite.protex \
	${GMDIR}/TestCompare.protex

.SUFFIXES:
.SUFFIXES: .f90 .pl .protex .tex .dvi .ps .pdf .html .XML .xmltex

.XML.xmltex:
	${SCRIPTDIR}/XmlToTex.pl $*.XML > $*.xmltex

.f90.protex:
	${SCRIPTDIR}/protex -b -n -f $*.f90 > $*.protex

.pl.protex:
	${SCRIPTDIR}/protex -b -n -f -S $*.pl > $*.protex

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
		  -address "(c) 2003 University of Michigan" \
		  -prev_url   "../../index.html"  \
		  -up_url     "../../index.html" \
		  -up_title   "SWMF Documents" \
		  -prev_title "SWMF Documents" \
		  $*
	cd $*; \
	perl -pi -e 's/<BODY >/<BODY bgcolor="#CECFCE">/' \
		$*.html node*.html index.html; \
	perl -pi -e 's/^\[\]//' $*.html index.html; \
	perl -pi -e 's/SRC="prev.png"\>\<\/A\>/SRC="prev.png"\>\<\/A\>\n \
		\<A HREF="http:\/\/csem.engin.umich.edu\/"\>\n \
		\<IMG WIDTH="100" HEIGHT="24" ALIGN="BOTTOM" BORDER="0" \
		ALT="CSEM Home" SRC="csem_home.png"\>\<\/A\>/' \
		$*.html node*.html index.html; \
	cp ../../HTML/csem_home.png ./ 
	rm -f $*.tex
	mv $*orig.tex $*.tex
	mv $* ../HTML
#^CFG END DOCHTML

help:
	@echo "       Make manuals for SWMF"
	@echo ""
	@echo "Targets:"
	@echo ""
	@echo "PDF          - produces SWMF.pdf"
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

SWMFREF.tex: ${SWMFREF}
	cat ${SWMFREF} > SWMFREF.tex

SWMFXML = \
	${DIR}/Param/PARAM.xmltex \
	${DIR}/GM/BATSRUS/PARAM.xmltex

SWMF_xml.tex: ${SWMFXML}
	cat ${SWMFXML} > SWMF_xml.tex

SWMF.dvi: \
	SWMF.tex \
	HEADER.tex \
	SWMF_introduction.tex \
	SWMF_configuration.tex \
	SWMF_xml.tex \
	SWMF_couple.tex \
	SWMF_new_comp.tex \
	SWMFREF.tex

PDF: SWMF.pdf

#^CFG IF DOCHTML BEGIN
SWMF.html: SWMFREF.tex
HTML: cleanhtml SWMF.html
#^CFG END DOCHTML

clean:
	rm -rf *.aux *.log *.toc *.dvi *~ *.f90 SWMFREF.tex
	@echo 'rm -f SWMFREF SWMFXML'
	@rm -f ${SWMFREF} ${SWMFXML}

cleanpdf: clean
	cd ..;     rm -f  SWMF.pdf

#^CFG IF DOCHTML BEGIN
cleanhtml: clean
	cd ../HTML;rm -rf SWMF

cleanall: cleanpdf cleanhtml
#^CFG END DOCHTML

