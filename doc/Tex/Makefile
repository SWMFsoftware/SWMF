#^CFG COPYRIGHT UM
include ../../Makefile.def

#^CMP IF GM BEGIN
BATSRUS_SCRIPTS = \
	${DIR}/GM/BATSRUS/TestBatsrus.protex \
	${DIR}/GM/BATSRUS/TestSuite.protex \
	${DIR}/GM/BATSRUS/TestCompare.protex \
	${DIR}/GM/BATSRUS/Scripts/ConvertRestart.protex
#^CMP END GM

SCRIPT_CONF = ${SCRIPTDIR}/Configure.protex #^CMP IF CONFIGURE

SWMFREF = \
	${DIR}/SetSWMF.protex \
	${CONTROLDIR}/CON_main.protex \
	${CONTROLDIR}/CON_variables.protex \
	${CONTROLDIR}/CON_session.protex \
	${CONTROLDIR}/CON_io.protex \
	${CONTROLDIR}/CON_methods.protex \
	${LIBRARYDIR}/CON_time.protex \
	${LIBRARYDIR}/CON_physics.protex \
	${LIBRARYDIR}/CON_world.protex \
	${LIBRARYDIR}/CON_comp_param.protex \
	${LIBRARYDIR}/CON_comp_info.protex \
	${CONDIR}/Interface/src/CON_wrapper.protex \
	${CONDIR}/Interface/src/CON_couple_all.protex \
	${CONDIR}/Interface/src/CON_couple_ih_gm.protex \
	${CONDIR}/Interface/src/CON_couple_gm_ie_swmf.protex \
	${CONDIR}/Interface/src/CON_couple_gm_im.protex \
	${CONDIR}/Interface/src/CON_couple_ie_im.protex \
	${CONDIR}/Interface/src/CON_couple_ie_ua.protex \
	${CONDIR}/Stubs/src/CON_wrapper.protex \
	${CONDIR}/Stubs/src/CON_couple_all.protex \
	${COUPLERDIR}/CON_domain_decomposition.protex \
	${COUPLERDIR}/CON_grid_descriptor.protex \
	${COUPLERDIR}/CON_router.protex \
	${COUPLERDIR}/CON_global_message_pass.protex \
	${SHAREDIR}/ModReadParam.protex \
	${SHAREDIR}/ModIoUnit.protex \
	${SHAREDIR}/ModKind.protex \
	${SHAREDIR}/ModMpi.protex \
	${SHAREDIR}/ModTimeConvert.protex \
	${SHAREDIR}/ModFreq.protex \
	${SHAREDIR}/ModUtilities.protex \
	${SHAREDIR}/ModCoordTransform.protex \
	${SHAREDIR}/CON_axes.protex \
	${SHAREDIR}/CON_planet.protex \
	${SHAREDIR}/CON_planet_field.protex \
	${SHAREDIR}/CON_ray_trace.protex \
	${SCRIPTDIR}/CheckParam.protex \
	${SCRIPTDIR}/ExpandParam.protex ${SCRIPT_CONF} \
	${SCRIPTDIR}/depend.protex \
	${SCRIPTDIR}/FixMakefileDef.protex \
	${SCRIPTDIR}/FixMpi.protex \
	${SCRIPTDIR}/Rename.protex \
	${SCRIPTDIR}/Methods.protex \
	${SCRIPTDIR}/StaticToDynamic.protex \
	${SCRIPTDIR}/ChangeIndexOrder.protex \
	${SCRIPTDIR}/CompareDirs.protex \
	${SCRIPTDIR}/Fix4Endian.protex \
	${SCRIPTDIR}/Fix8Endian.protex \
	${SCRIPTDIR}/FixI4toI8.protex \
	${SCRIPTDIR}/FixI8toI4.protex \
	${SCRIPTDIR}/XmlToTex.protex \
	${SCRIPTDIR}/protex.protex \
	${UTILDIR}/NOMPI/src/NOMPI.protex \
	${UTILDIR}/TIMING/src/ModTiming.protex \
	${UTILDIR}/TIMING/src/timing.protex \
	${BATSRUS_SCRIPTS}

.SUFFIXES:
.SUFFIXES: .f90 .F90 .pl .protex .tex .dvi .ps .pdf .html .XML .xmltex

.XML.xmltex:
	${SCRIPTDIR}/XmlToTex.pl $*.XML > $*.xmltex

.f90.protex:
	${SCRIPTDIR}/protex.pl -b -n -f $*.f90 > $*.protex

.F90.protex:
	${SCRIPTDIR}/protex.pl -b -n -f $*.F90 > $*.protex


${SCRIPTDIR}/CompareDirs.protex: ${SCRIPTDIR}/CompareDirs
	${SCRIPTDIR}/protex.pl -b -n -f -S \
		${SCRIPTDIR}/CompareDirs > ${SCRIPTDIR}/CompareDirs.protex

.pl.protex:
	${SCRIPTDIR}/protex.pl -b -n -f -S $*.pl > $*.protex

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

GMXML = ${DIR}/GM/BATSRUS/PARAM.xmltex #	^CMP IF GM
IMXML = ${DIR}/IM/RCM/PARAM.xmltex #		^CMP IF IM
IEXML = ${DIR}/IE/Ridley_serial/PARAM.xmltex #	^CMP IF IE

SWMFXML = ${DIR}/Param/PARAM.xmltex ${GMXML} ${IMXML} ${IEXML}

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
SWMF.html: SWMFREF.tex SWMF_xml.tex
HTML: cleanhtml
	make SWMF.html
#^CFG END DOCHTML

clean:
	rm -f *.aux *.log *.toc *.dvi *~ *.f90 SWMFREF.tex SWMF_xml.tex
	@echo 'rm -f SWMFREF SWMFXML'
	@rm -f ${SWMFREF} ${SWMFXML}

cleanpdf: clean
	cd ..;     rm -f  SWMF.pdf

#^CFG IF DOCHTML BEGIN
cleanhtml: clean
	cd ../HTML;rm -rf SWMF

cleanall: cleanpdf cleanhtml
#^CFG END DOCHTML
