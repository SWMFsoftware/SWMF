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
	@echo "PDF          - produces PDF manuals in the ../ directory"
	@#^CMP IF DOCHTML BEGIN
	@echo ""
	@echo "HTML         - produces online manuals in ../HTML directory"
	@#^CMP END DOCHTML
	@echo ""
	@echo "clean        - remove temporary files"
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
COUPLE_IHGM = ${CONDIR}/Interface/src/CON_couple_ih_gm.protex #     ^CMP IF IH
COUPLE_GMIE = ${CONDIR}/Interface/src/CON_couple_gm_ie.protex #     ^CMP IF IE
COUPLE_GMIM = ${CONDIR}/Interface/src/CON_couple_gm_im.protex #     ^CMP IF IM
#^CMP END GM

#^CMP IF IE BEGIN
COUPLE_IEIM = ${CONDIR}/Interface/src/CON_couple_ie_im.protex #     ^CMP IF IM
COUPLE_IEUA = ${CONDIR}/Interface/src/CON_couple_ie_ua.protex #     ^CMP IF UA
#^CMP END IE

SCRIPT_CONF = ${SCRIPTDIR}/Configure.protex #^CMP IF CONFIGURE

BATSRUS_SCRIPTS = \
	${DIR}/GM/BATSRUS/TestBatsrus.protex \
	${DIR}/GM/BATSRUS/TestSuite.protex \
	${DIR}/GM/BATSRUS/TestCompare.protex \
	${DIR}/GM/BATSRUS/Scripts/ConvertRestart.protex

BATSRUS_SCRIPTS:
	@rm -f BATSRUS_SCRIPTS.tex
	@touch BATSRUS_SCRIPTS.tex
	@(if [ -d ${DIR}/GM/BATSRUS ]; then \
		echo "make GM/BATSRUS/*.protex"; \
		make ${BATSRUS_SCRIPTS}; \
		echo "cat GM/BATSRUS/*.protex > BATSRUS_SCRIPTS.tex"; \
		cat ${BATSRUS_SCRIPTS} > BATSRUS_SCRIPTS.tex; \
	fi)

SWMFREF = \
	${DIR}/Config.protex ${SCRIPT_CONF} \
	${DIR}/Scripts/TestParam.protex \
	${DIR}/Scripts/Performance.protex \
	BATSRUS_SCRIPTS.tex \
	${SCRIPTDIR}/CheckParam.protex \
	${SCRIPTDIR}/ParamConvert.protex \
	${SCRIPTDIR}/ParamEditor.protex \
	${SCRIPTDIR}/PostProc.protex \
	${SCRIPTDIR}/Restart.protex \
	${SCRIPTDIR}/depend.protex \
	${SCRIPTDIR}/FixMakefileDef.protex \
	${SCRIPTDIR}/CheckDataName.protex \
	${SCRIPTDIR}/Rename.protex \
	${SCRIPTDIR}/Methods.protex \
	${SCRIPTDIR}/FixMpi.protex \
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
	${DIR}/Makefile.protex \
	${DIR}/share/Makefile.protex \
	${DIR}/share/build/Makefile.protex \
	${DIR}/share/build/Makefile.doc.protex \
	${DIR}/doc/Tex/Makefile.protex \
	${DIR}/CON/Makefile.protex \
	${CONTROLDIR}/swmf.protex \
	${CONTROLDIR}/swmf_interface.protex \
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
	${UTILDIR}/NOMPI/src/NOMPI.protex \
	${UTILDIR}/TIMING/src/ModTiming.protex \
	${UTILDIR}/TIMING/src/timing.protex

#^CMP IF DOCHTML BEGIN
.tex.html:
	touch $*.ind $*.idx
	latex $*
	latex $*
	makeindex $*.idx
	latex2html -image_type gif -local_icons \
		  -address "(c) 2003 University of Michigan" \
		  -prev_url   "../../index.html"  \
		  -up_url     "../../index.html" \
		  -up_title   "SWMF Documents" \
		  -prev_title "SWMF Documents" \
		  -noauto_navigation \
		  -top_navigation \
		  $*
	cd $*; \
	perl -pi -e 's/^\[\]//' $*.html index.html; \
	cp ../../HTML/*.png ./ ; \
	rm -f $*.ind $*.idx $*.ilg
	rm -rf ../HTML/$*
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

${DIR}/share/build/Makefile.protex: ${DIR}/share/build/Makefile.Linux.pgf90
	${SCRIPTDIR}/protex.pl -b -n -f -S \
		${DIR}/share/build/Makefile.Linux.pgf90 > \
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

SWMFREF.tex: ${SWMFREF} BATSRUS_SCRIPTS
	cat ${SWMFREF} > SWMFREF.tex

SWMFXML:
	@rm -f SWMFXML.tex
	@touch SWMFXML.tex
	@perl -e 'for $$i (glob("../../Param/PARAM.XML ../../[A-Z][A-Z]/*/PARAM.XML")){`${SCRIPTDIR}/XmlToTex.pl $$i >> SWMFXML.tex`}'
	@echo "${SCRIPTDIR}/XmlToTex.pl */PARAM.XML > SWMFXML.tex"

SWMF.dvi: \
	SWMF.tex \
	HEADER.tex \
	SWMF_introduction.tex \
	SWMF_configuration.tex \
	SWMF_param.tex \
	SWMF_control.tex \
	SWMF_example_runs.tex \
	SWMFXML

REFERENCE.dvi:\
	SWMFREF.tex

MAINTENANCE.dvi:\
	PM_wrapper.protex \
	PM_XY_coupler.protex

PDF_FILES = SWMF.pdf TESTING.pdf REFERENCE.pdf MAINTENANCE.pdf PHYSICS.pdf \
		SOFTWARE_STANDARD.pdf

PDF: ${PDF_FILES}

#^CMP IF DOCHTML BEGIN
HTML_FILES = SWMF.html REFERENCE.html
HTML_DIRS  = SWMF REFERENCE

SWMF.html: SWMFXML
REFERENCE.html: SWMFREF.tex
HTML: cleanhtml
	make ${HTML_FILES}
#^CMP END DOCHTML

clean:
	rm -f *.out *.aux *.log *.toc *.ind *.idx *.ilg *.dvi *.ps *~ *.protex
	rm -f SWMFREF.tex SWMFXML.tex BATSRUS_SCRIPTS.tex
	@echo 'rm -f SWMFREF SWMFXML BATSRUS_SCRIPTS'
	@rm -f ${SWMFREF} ${SWMFXML} ${BATSRUS_SCRIPTS}

cleanpdf: clean
	cd ..; rm -f ${PDF_FILES}

#^CMP IF DOCHTML BEGIN
cleanhtml: clean
	cd ../HTML;rm -rf ${HTML_DIRS}

cleanall: cleanpdf cleanhtml
#^CMP END DOCHTML
