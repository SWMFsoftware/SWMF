lib/libFISHPAK.a:
	tar -xzf fishpak.tar.gz
	mkdir lib
	cd src; make all

install: lib/libFISHPAK.a

test:	install
	cd test; make test

clean:	
	@(if [ -d src  ];   then cd src;    make clean; fi)
	@(if [ -d test ];   then cd test;   make clean; fi)

distclean:
	@(if [ -d src  ];   then rm -rf src;  fi)
	@(if [ -d test ];   then rm -rf test; fi)
	@(if [ -d lib  ];   then rm -rf lib;  fi)

uninstall: distclean
