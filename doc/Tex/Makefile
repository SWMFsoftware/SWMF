#^CMP COPYRIGHT UM
#BOP
#!ROUTINE: doc/Tex/Makefile - produce documentation from LaTex, XML, Protex
#!DESCRIPTION:
# This Makefile produces documentation for SWMF from Latex, PARAM.XML
# and Protex documented source code and scripts. 
#EOP
#BOC
include ../../Makefile.def
include ${DIR}/share/build/Makefile.doc

help:
	@echo "       Make manuals for SWMF"
	@echo ""
	@echo "Targets:"
	@echo ""
	@echo "PDF          - produces SWMF.pdf and TESTING.pdf"
	@#^CMP IF DOCHTML BEGIN
	@echo ""
	@echo "HTML         - produces online manuals in ../HTML directory"
	@#^CMP END DOCHTML
	@echo ""
	@echo "clean        - rm -rf *.aux *.log *.toc *.dvi *~"
	@echo ""
	@echo "cleanpdf     - clean + remove PDF manuals"
	@#^CMP IF DOCHTML BEGIN
	@echo ""
	@echo "cleanhtml    - clean + remove HTML manuals"
	@echo ""
	@echo "cleanall     - clean + remove PDF and HTML manuals"
	@#^CMP END DOCHTML
#EOC

#^CMP IF GM BEGIN
BATSRUS_SCRIPTS = \
	${DIR}/GM/BATSRUS/TestBatsrus.protex \
	${DIR}/GM/BATSRUS/TestSuite.protex \
	${DIR}/GM/BATSRUS/TestCompare.protex \
	${DIR}/GM/BATSRUS/Scripts/ConvertRestart.protex

COUPLE_IHGM = ${CONDIR}/Interface/src/CON_couple_ih_gm.protex #     ^CMP IF IH
COUPLE_GMIE = ${CONDIR}/Interface/src/CON_couple_gm_ie_swmf.protex #^CMP IF IE
COUPLE_GMIM = ${CONDIR}/Interface/src/CON_couple_gm_im.protex #     ^CMP IF IM
#^CMP END GM

#^CMP IF IE BEGIN
COUPLE_IEIM = ${CONDIR}/Interface/src/CON_couple_ie_im.protex #     ^CMP IF IM
COUPLE_IEUA = ${CONDIR}/Interface/src/CON_couple_ie_ua.protex #     ^CMP IF UA
#^CMP END IE

SCRIPT_CONF = ${SCRIPTDIR}/Configure.protex #^CMP IF CONFIGURE

SWMFREF = \
	${DIR}/SetSWMF.protex ${SCRIPT_CONF} \
	${DIR}/share/Makefile.protex \
	${DIR}/Makefile.protex \
	${DIR}/CON/Makefile.protex \
	${DIR}/share/build/Makefile.protex \
	${DIR}/Scripts/TestParam.protex \
	${SCRIPTDIR}/CheckParam.protex \
	${SCRIPTDIR}/ExpandParam.protex \
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
        ${COUPLE_IHGM} ${COUPLE_GMIE} ${COUPLE_GMIM} \
	${COUPLE_IEIM} ${COUPLE_IEUA} \
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
	${SCRIPTDIR}/depend.protex \
	${SCRIPTDIR}/FixMakefileDef.protex \
	${SCRIPTDIR}/FixMpi.protex \
	${SCRIPTDIR}/Rename.protex \
	${SCRIPTDIR}/Methods.protex \
	${SCRIPTDIR}/StopToConStop.protex \
	${SCRIPTDIR}/StaticToDynamic.protex \
	${SCRIPTDIR}/ChangeIndexOrder.protex \
	${SCRIPTDIR}/CompareDirs.protex \
	${SCRIPTDIR}/Fix4Endian.protex \
	${SCRIPTDIR}/Fix8Endian.protex \
	${SCRIPTDIR}/FixI4toI8.protex \
	${SCRIPTDIR}/FixI8toI4.protex \
	${SCRIPTDIR}/XmlToTex.protex \
	${SCRIPTDIR}/protex.protex \
	${DIR}/share/build/Makefile.doc.protex \
	${DIR}/doc/Tex/Makefile.protex \
	${UTILDIR}/NOMPI/src/NOMPI.protex \
	${UTILDIR}/TIMING/src/ModTiming.protex \
	${UTILDIR}/TIMING/src/timing.protex \
	${BATSRUS_SCRIPTS}

#^CMP IF DOCHTML BEGIN
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
#^CMP END DOCHTML

#
# Rules for files without a standard extension
#

${DIR}/Makefile.protex: ${DIR}/Makefile
	${SCRIPTDIR}/protex.pl -b -n -f -S \
		${DIR}/Makefile > ${DIR}/Makefile.protex

${DIR}/share/Makefile.protex: ${DIR}/share/Makefile
	${SCRIPTDIR}/protex.pl -b -n -f -S \
		${DIR}/share/Makefile > ${DIR}/share/Makefile.protex

${DIR}/CON/Makefile.protex: ${DIR}/CON/Makefile.def
	${SCRIPTDIR}/protex.pl -b -n -f -S \
		${DIR}/CON/Makefile.def > ${DIR}/CON/Makefile.protex

${DIR}/share/build/Makefile.protex: ${DIR}/share/build/Makefile.Linuxpgf90
	${SCRIPTDIR}/protex.pl -b -n -f -S \
		${DIR}/share/build/Makefile.Linuxpgf90 > \
		${DIR}/share/build/Makefile.protex

${DIR}/share/build/Makefile.doc.protex: ${DIR}/share/build/Makefile.doc
	${SCRIPTDIR}/protex.pl -b -n -f -S \
		${DIR}/share/build/Makefile.doc > \
		${DIR}/share/build/Makefile.doc.protex

${DIR}/doc/Tex/Makefile.protex: ${DIR}/doc/Tex/Makefile
	${SCRIPTDIR}/protex.pl -b -n -f -S \
		${DIR}/doc/Tex/Makefile > ${DIR}/doc/Tex/Makefile.protex

${SCRIPTDIR}/CompareDirs.protex: ${SCRIPTDIR}/CompareDirs
	${SCRIPTDIR}/protex.pl -b -n -f -S \
		${SCRIPTDIR}/CompareDirs > ${SCRIPTDIR}/CompareDirs.protex

${SWMFREF} : ${SCRIPTDIR}/protex.pl

SWMFREF.tex: ${SWMFREF}
	cat ${SWMFREF} > SWMFREF.tex

GMXML = ${DIR}/GM/BATSRUS/PARAM.xmltex #			  ^CMP IF GM
IMXML = ${DIR}/IM/RCM/PARAM.xmltex #				  ^CMP IF IM
RBXML = ${DIR}/RB/RiceV5/PARAM.xmltex #				  ^CMP IF RB
IEXML = ${DIR}/IE/Ridley_serial/PARAM.xmltex #			  ^CMP IF IE
SPXML = ${DIR}/SP/Kota/PARAM.xmltex #           		  ^CMP IF SP
UAXML = ${DIR}/UA/GITM/PARAM.xmltex ${DIR}/UA/GITM2/PARAM.xmltex #^CMP IF UA

SWMFXML = ${DIR}/Param/PARAM.xmltex \
	${GMXML} ${IEXML} ${UAXML} ${IMXML} ${RBXML} ${SPXML}

${SWMFXML}: ${SCRIPTDIR}/XmlToTex.pl

SWMFXML.tex: ${SWMFXML}
	cat ${SWMFXML} > SWMFXML.tex

SWMF.dvi: \
	SWMF.tex \
	HEADER.tex \
	SWMF_introduction.tex \
	SWMF_configuration.tex \
	SWMF_param.tex \
	SWMF_control.tex \
	SWMFXML.tex \
	SWMF_couple.tex \
	SWMF_new_comp.tex \
	SWMFREF.tex

PDF: SWMF.pdf TESTING.pdf

#^CMP IF DOCHTML BEGIN
SWMF.html: SWMFREF.tex SWMFXML.tex
HTML: cleanhtml
	make SWMF.html
#^CMP END DOCHTML

clean:
	rm -f *.aux *.log *.toc *.dvi *~ *.f90 SWMFREF.tex SWMFXML.tex
	@echo 'rm -f SWMFREF SWMFXML'
	@rm -f ${SWMFREF} ${SWMFXML}

cleanpdf: clean
	cd ..;     rm -f  SWMF.pdf TESTING.pdf

#^CMP IF DOCHTML BEGIN
cleanhtml: clean
	cd ../HTML;rm -rf SWMF

cleanall: cleanpdf cleanhtml
#^CMP END DOCHTML
