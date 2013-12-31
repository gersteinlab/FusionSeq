#ifndef DEF_UTIL_H
#define DEF_UTIL_H

#include <bios/blastParser.h>
#include <bios/blatParser.h>

typedef struct {
  char* transcriptName;
  char* swissProt;
  char* uniprotId;
  char* geneSymbol;
  char* refseqId;
  char* refseqDescription;
} KgXref;



typedef struct {
  char* transcriptName;
  char* treeFamId;
} KgTreeFam;


extern int getNucleotideOverlap ( BlatQuery* blQ );
extern Array util_readKnownGeneXrefs (char* fileName);
extern Array util_readKnownGeneTreeFams (char* fileName);
extern int sortKgXrefsByTranscriptName (KgXref *a, KgXref *b);
extern void transcript2geneSymbolAndGeneDescription (Array kgXrefs, char *transcriptName, char** geneSymbol, char **description);
 

#endif
