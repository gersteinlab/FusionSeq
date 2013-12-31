EXTRA_PROGRAMS = geneFusions_cgi \
		 showDetails_cgi \
		 seqViz_cgi \
		 findFusionPartner.cgi

cgi: $(EXTRA_PROGRAMS)

CLEANFILES = $(EXTRA_PROGRAMS)

.PHONY: cgi

AM_CPPFLAGS = -I ../

# -----------------------------------------------------------------------------
# CGI
# -----------------------------------------------------------------------------

geneFusions_cgi_SOURCES = geneFusions_cgi.c
geneFusions_cgi_LDADD = ../libfusionseq.la -lbios

showDetails_cgi_SOURCES = showDetails_cgi.c
showDetails_cgi_LDADD = ../libfusionseq.la -lbios

seqViz_cgi_SOURCES = seqViz_cgi.c
seqViz_cgi_LDADD = ../libfusionseq.la -lbios -lmrf

findFusionPartner_cgi_SOURCES = findFusionPartner_cgi.c
findFusionPartner_cgi_LDADD = ../libfusionseq.la -lbios
