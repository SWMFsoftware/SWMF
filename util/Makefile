INSTALLFILE = \

touch_install_files:
	-@touch DATAREAD/srcIndices/Makefile.DEPEND
	-@touch EMPIRICAL/srcEE/Makefile.DEPEND
	-@touch EMPIRICAL/srcIE/Makefile.DEPEND
	-@touch EMPIRICAL/srcUA/Makefile.RULES
	-@touch NOMPI/src/Makefile.RULES
	-@touch CRASH/src/Makefile.DEPEND
	-@touch CRASH/src/Makefile.RULES

install: touch_install_files
	@(if [ -d HYPRE ]; then cd HYPRE;  make install; fi);

clean:: touch_install_files
	@(if [ -d NOMPI ];     then cd NOMPI/src;  make clean; fi)
	@(if [ -d TIMING ];    then cd TIMING;     make clean; fi)
	@(if [ -d DATAREAD ];  then cd DATAREAD;   make clean; fi)
	@(if [ -d EMPIRICAL ]; then cd EMPIRICAL;  make clean; fi)
	@(if [ -d CRASH ];     then cd CRASH;      make clean; fi)
	@(if [ -d HYPRE ];     then cd HYPRE;      make clean; fi)

distclean:: touch_install_files
	@(if [ -d NOMPI ];     then cd NOMPI/src;  make distclean; fi)
	@(if [ -d TIMING ];    then cd TIMING;     make distclean; fi)
	@(if [ -d DATAREAD ];  then cd DATAREAD;   make distclean; fi)
	@(if [ -d EMPIRICAL ]; then cd EMPIRICAL;  make distclean; fi)
	@(if [ -d CRASH ];     then cd CRASH;      make distclean; fi)
	@(if [ -d HYPRE ];     then cd HYPRE;      make distclean; fi)
	rm -f *~ */src*/Makefile.DEPEND */src*/Makefile.RULES
