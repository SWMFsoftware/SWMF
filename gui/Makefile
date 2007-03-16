# Default port number unless set otherwise
PORT=8080

#
# The default target is 'install' so it is listed first
#
default : install

#
# Install
#
install: 
	./install.sh ${PORT}
	@echo 'GUI install finished.'

#
# Code documentation
#
doc: 
	./makemanual.sh
	@echo 'SWMF manual created.'

#
# Cleanup
#
clean: 
	@echo ' ... removing *~ files ...'
	@find . -name "*~" -print -exec rm -f {} \;

cleanplots: 
	rm -rf htdocs/plots

cleanruns: 
	rm -rf htdocs/runs

cleancodes: 
	rm -rf htdocs/codes

distclean: clean cleanplots cleanruns cleancodes
	rm -rf apache PORT* htdocs/Manuals htdocs/config.php security/passwords

