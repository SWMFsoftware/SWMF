#^CMP COPYRIGHT UM
SHELL = /bin/sh

help:
	@echo 'Targets to make for the TIMING module:'
	@echo ' '
	@echo 'help              - this help (default target)'
	@echo ' '
	@echo 'MAN               - MAN_TIMING.ps.gz manual'
	@echo ' '
	@echo 'clean             - /bin/rm -f *.o *.exe'
	@echo 'distclean         - /bin/rm -f *.o *.exe *.a *~'
	@echo ' '
	@echo 'dist              - /bin/tar -cf ../TIMING.tar .'

MAN:	MAN_TIMING.ps.gz

MAN_TIMING.dvi : MAN_TIMING.tex
	latex MAN_TIMING.tex
	latex MAN_TIMING.tex

MAN_TIMING.ps : MAN_TIMING.dvi
	dvips -o MAN_TIMING.ps MAN_TIMING.dvi

MAN_TIMING.ps.gz : MAN_TIMING.ps
	gzip -f MAN_TIMING.ps

clean:
	/bin/rm -f *.log *.aux *.dvi *~

distclean: clean

