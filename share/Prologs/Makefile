include ../../Makefile.def
include ../../Makefile.conf

# This Makefile is a template for creating documentation from 
# Protex documented F90 source code and scripts

# SCRIPTDIR = ../Scripts

PROTEXFILES = \
	program.protex \
	module.protex \
	subroutine.protex \
	internal_subroutine.protex

OBJECTS = \
	module.o \
	subroutine.o \
	internal_subroutine.o \
	program.o

README.dvi: README.tex PROTEX.tex

PROTEX.tex: ${PROTEXFILES}
	cat ${PROTEXFILES} > PROTEX.tex

# TEST TARGETS

help:
	@echo 'The following targets are available:'
	@echo
	@echo 'test      - make COMPILE DVI PS PDF'
	@echo 'DVI       - create README.dvi'
	@echo 'PS        - create README.ps'
	@echo 'DVI       - create README.pdf'
	@echo 'COMPILE   - compile the Fortran files into program_example.exe'
	@echo 'RUN       - run program_example.exe'
	@echo 'clean     - remove *~ and *.protex files'
	@echo 'distclean - remove all temporary files'
	@echo

test:
	make COMPILE
	make DVI
	make PS
	make PDF

DVI: README.dvi

PS:  README.ps

PDF: README.pdf

RUN: COMPILE
	./program_example.exe

COMPILE: program_example.exe

program_example.exe: ${OBJECTS}
	${LINK.f90} -o program_example.exe ${OBJECTS}


# CLEANING

distclean: clean
	rm -f *.dvi *.ps *.pdf PROTEX.tex program_example.exe

# RULES

.SUFFIXES:
.SUFFIXES: .f90 .o .pl .protex .tex .dvi .ps .pdf

.f90.protex:
	${SCRIPTDIR}/protex -b -n -f $*.f90 > $*.protex

.pl.protex:
	${SCRIPTDIR}/protex -b -n -f -S $*.pl > $*.protex

# NOTE: Repeat latex $*.tex twice to get references right
# NOTE: Repeat latex $*.tex three times to get table of contents right
.tex.dvi:
	latex $*.tex

.dvi.ps: 
	dvips -o $*.ps $*.dvi

.ps.pdf:
	ps2pdf $*.ps $*.pdf

