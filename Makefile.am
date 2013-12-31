SUBDIRS = src
ACLOCAL_AMFLAGS = -I m4

DIST_SUBDIRS = src src/optional src/cgi src/test
EXTRA_DIST = scripts default.fusionseqrc

optional:
	cd src/optional && $(MAKE) $(AM_MAKEFLAGS)

cgi:
	cd src/cgi && $(MAKE) $(AM_MAKEFLAGS) cgi

test:
	cd src/test && $(MAKE) $(AM_MAKEFLAGS)


clean-local:
	cd src/optional && $(MAKE) clean
	cd src/test && $(MAKE) clean
	cd src/cgi && $(MAKE) clean

install-exec-local:
	if test -e src/optional/bp2wig; then cd src/optional && $(MAKE) install; fi
	if test -e src/test/bp2fasta; then cd src/test && $(MAKE) install; fi

install-data-hook:
	scripts/install_rc.sh

deploy:
	scripts/install_deploy.pl


.PHONY: clean-local install-exec-local install-data-hook deploy
