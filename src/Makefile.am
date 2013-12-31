# -----------------------------------------------------------------------------
# Common library
# -----------------------------------------------------------------------------

lib_LTLIBRARIES = libfusionseq.la
libfusionseq_la_SOURCES = bp.c \
			  gfr.c \
			  sqvCircos.c \
			  sqvUtil.c \
			  sqvWeb.c \
			  util.c

libfusionseq_la_LIBADD = -lm -lbios

include_HEADERS = bp.h \
		  gfr.h \
		  sqvCircos.h \
		  sqvUtil.h \
		  sqvWeb.h \
		  util.h

# -----------------------------------------------------------------------------
# FusionSeq programs
# -----------------------------------------------------------------------------

bin_PROGRAMS = geneFusions \
	       gfrAnnotationConsistencyFilter \
	       gfrLargeScaleHomologyFilter \
	       gfrSmallScaleHomologyFilter \
	       gfrAbnormalInsertSizeFilter \
	       gfrAddInfo \
	       gfrProximityFilter \
	       bpFilter \
	       bp2alignment \
	       gfrPCRFilter \
	       gfrBlackListFilter \
	       gfr2bpJunctions \
	       bowtie2bp \
	       gfrConfidenceValues \
	       validateBpJunctions \
	       gfrExpressionConsistencyFilter \
	       gfrMitochondrialFilter \
	       gfrRepeatMaskerFilter \
	       gfrCountPairTypes \
	       gfrRibosomalFilter \
	       gfrSpliceJunctionFilter \
	       gfrClassify \
	       gfrWhiteListFilter

# -----------------------------------------------------------------------------
# CORE: Finding candidate fusions
# -----------------------------------------------------------------------------

geneFusions_SOURCES = geneFusions.c
geneFusions_LDADD = libfusionseq.la -lbios -lmrf -lm

gfrExpressionConsistencyFilter_SOURCES = gfrExpressionConsistencyFilter.c
gfrExpressionConsistencyFilter_LDADD = libfusionseq.la -lbios -lm

gfrAnnotationConsistencyFilter_SOURCES = gfrAnnotationConsistencyFilter.c
gfrAnnotationConsistencyFilter_LDADD = libfusionseq.la -lbios -lm

gfrRibosomalFilter_SOURCES = gfrRibosomalFilter.c
gfrRibosomalFilter_LDADD = libfusionseq.la -lbios

gfrMitochondrialFilter_SOURCES = gfrMitochondrialFilter.c
gfrMitochondrialFilter_LDADD = libfusionseq.la -lbios

gfrLargeScaleHomologyFilter_SOURCES = gfrLargeScaleHomologyFilter.c
gfrLargeScaleHomologyFilter_LDADD = libfusionseq.la -lbios

gfrSmallScaleHomologyFilter_SOURCES = gfrSmallScaleHomologyFilter.c
gfrSmallScaleHomologyFilter_LDADD = libfusionseq.la -lbios

gfrAbnormalInsertSizeFilter_SOURCES = gfrAbnormalInsertSizeFilter.c
gfrAbnormalInsertSizeFilter_LDADD = libfusionseq.la -lbios

gfrPCRFilter_SOURCES = gfrPCRFilter.c
gfrPCRFilter_LDADD = libfusionseq.la -lbios

gfrBlackListFilter_SOURCES = gfrBlackListFilter.c
gfrBlackListFilter_LDADD = libfusionseq.la -lbios

gfrProximityFilter_SOURCES = gfrProximityFilter.c
gfrProximityFilter_LDADD = libfusionseq.la -lbios

gfrAddInfo_SOURCES = gfrAddInfo.c
gfrAddInfo_LDADD = libfusionseq.la -lbios

gfrConfidenceValues_SOURCES = gfrConfidenceValues.c
gfrConfidenceValues_LDADD = libfusionseq.la -lbios -lmrf

gfrRepeatMaskerFilter_SOURCES = gfrRepeatMaskerFilter.c
gfrRepeatMaskerFilter_LDADD = libfusionseq.la -lbios -lm

gfrCountPairTypes_SOURCES = gfrCountPairTypes.c
gfrCountPairTypes_LDADD = libfusionseq.la -lbios -lm

gfrSpliceJunctionFilter_SOURCES = gfrSpliceJunctionFilter.c
gfrSpliceJunctionFilter_LDADD = libfusionseq.la -lbios

gfrClassify_SOURCES = gfrClassify.c
gfrClassify_LDADD = libfusionseq.la -lbios -lm

gfrWhiteListFilter_SOURCES = gfrWhiteListFilter.c
gfrWhiteListFilter_LDADD = libfusionseq.la -lbios

# -----------------------------------------------------------------------------
# CORE: Identifying sequences of the junction
# -----------------------------------------------------------------------------

gfr2bpJunctions_SOURCES = gfr2bpJunctions.c
gfr2bpJunctions_LDADD = libfusionseq.la -lbios -lm

bowtie2bp_SOURCES = bowtie2bp.c
bowtie2bp_LDADD = libfusionseq.la -lbios

validateBpJunctions_SOURCES = validateBpJunctions.c
validateBpJunctions_LDADD = libfusionseq.la -lbios

bpFilter_SOURCES = bpFilter.cpp
bpFilter_CXXFLAGS = -D_REENTRANT -pthread -m64
bpFilter_LDADD = libfusionseq.la -lstdc++ -lbios -L@ROOTLIBDIR@ @ROOTGLIBS@  @ROOTLIBS@ -ldl

bp2alignment_SOURCES = bp2alignment.c
bp2alignment_LDADD = libfusionseq.la -lbios

