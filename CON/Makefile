# defines OS
include ../Makefile.def

CONTROLRULES = Control/src/Makefile.RULES

INSTALL_FILES = \
	Control/src/Makefile.DEPEND \
	Library/src/Makefile.DEPEND \
	Coupler/src/Makefile.DEPEND \
	Interface/src/Makefile.DEPEND \
	Stubs/src/Makefile.DEPEND \
	${CONTROLRULES}


install:
	touch ${INSTALL_FILES}
	@(if [ -f ${CONTROLRULES}.${OS}.${COMPILER} ]; then   \
		cp -f ${CONTROLRULES}.${OS}.${COMPILER} ${CONTROLRULES};\
	else \
                rm -f ${CONTROLRULES}; touch ${CONTROLRULES}; \
	fi);

clean:
	cd Control/src;   make clean
	cd Library/src;   make clean
	cd Coupler/src;   make clean
	cd Interface/src; make clean
	cd Stubs/src;     make clean

distclean: clean
	rm -f ${INSTALL_FILES} *~
