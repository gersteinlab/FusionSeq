#include <bios/log.h>
#include <bios/format.h>
#include <bios/linestream.h>
#include <bios/bits.h>

#include "util.h"

int getNucleotideOverlap ( BlatQuery* blQ ) {
  int l;
  PslEntry* blE=NULL;     
  int qSize = arrp( blQ->entries, arrayMax(blQ->entries)-1, PslEntry)->qSize;
  int maxOverlap=0;
  //Bits* bitOverlap = bitAlloc( qSize  );	
  //bitClear( bitOverlap, qSize );
  //warn( "%d", arrayMax( blQ->entries ) );
  for( l=0; l < arrayMax ( blQ->entries ); l++ ) {
    blE = arrp( blQ->entries, l, PslEntry );
    if( blE->qSize != qSize ) die("Query size different from PslEntry query size: qSize(%d) - blE->qSize(%d)", qSize, blE->qSize);    
    if( ((blE->qEnd - blE->qStart)+1) > maxOverlap ) maxOverlap = (blE->qEnd - blE->qStart)+1;
    /*    if( (blE->qEnd > blE->qStart)  ) { // to ensure that at least two nucleotides matches and only significant hits are considered // 
      Bits* currBitEntry = bitAlloc( blE->qSize );
      bitClear ( currBitEntry, blE->qSize );
      if( (blE->qEnd - blE->qStart)> blE->qSize ) die("The query match is bigger than the read size: name(%s) - qSize(%d) - actualSize(%d)", blQ->qName, blE->qSize, (blE->qEnd - blE->qStart));
      bitSetRange( currBitEntry, (blE->qStart - 1), (blE->qEnd - blE->qStart));
      bitOr( bitOverlap, currBitEntry, qSize );
      bitFree( &currBitEntry );
    }//*/ 
  }
  //int overlap = 50; //bitCountRange( bitOverlap, 0, readSize);
  //bitFree( &bitOverlap);
  return maxOverlap;
}

Array util_readKnownGeneXrefs (char* fileName)
{
  WordIter w;
  LineStream ls;
  char *line,*pos;
  Array kgXrefs;
  KgXref *currKgXref;

  kgXrefs = arrayCreate (50000,KgXref);
  ls = ls_createFromFile (fileName);
  while (line = ls_nextLine (ls)) {
    if (line[0] == '\0') {
      continue;
    }
    currKgXref = arrayp (kgXrefs,arrayMax (kgXrefs),KgXref);
    w = wordIterCreate (line,"\t",0);
    currKgXref->transcriptName = hlr_strdup (wordNext (w));
    wordNext (w);
    currKgXref->swissProt = hlr_strdup (wordNext (w));
    currKgXref->uniprotId = hlr_strdup (wordNext (w));
    currKgXref->geneSymbol = hlr_strdup (wordNext (w));
    currKgXref->refseqId = hlr_strdup (wordNext (w));
    wordNext (w);
    currKgXref->refseqDescription = hlr_strdup (wordNext (w));
    wordIterDestroy (w);
    if (pos = strchr (currKgXref->uniprotId,'-')) {
      *pos = '\0';
    }
    if (pos = strchr (currKgXref->swissProt,'-')) {
      *pos = '\0';
    }
  }
  ls_destroy (ls);
  return kgXrefs;
}



Array util_readKnownGeneTreeFams (char* fileName) 
{
  LineStream ls;
  char* line;
  char* pos;
  Array kgTreeFams;
  KgTreeFam *currKgTreeFam;

  kgTreeFams = arrayCreate (30000,KgTreeFam);
  ls = ls_createFromFile (fileName);
  while (line = ls_nextLine (ls)) {
    if (line[0] == '\0') {
      continue;
    }
    if (pos = strchr (line,'\t')) {
      *pos = '\0';
      currKgTreeFam = arrayp (kgTreeFams,arrayMax (kgTreeFams),KgTreeFam);
      currKgTreeFam->transcriptName = hlr_strdup (line);
      currKgTreeFam->treeFamId = hlr_strdup (pos + 1);
    }
  }
  ls_destroy (ls);
  return kgTreeFams;
}



int sortKgXrefsByTranscriptName (KgXref *a, KgXref *b) 
{
  return strcmp (a->transcriptName,b->transcriptName);
}



static char* convert2string (Texta t)
{
  static Stringa buffer = NULL;
  int i;

  stringCreateClear (buffer,100);
  for (i = 0; i < arrayMax (t); i++) {
    stringAppendf (buffer,"%s%s",textItem (t,i),i < arrayMax (t) - 1 ? "|" : "");
  }
  return string (buffer);
}



void transcript2geneSymbolAndGeneDescription (Array kgXrefs, char *transcriptName, char** geneSymbol, char **description)
{
  Texta tokens;
  int i;
  KgXref testKX,*currKX;
  int index; 
  static Texta descriptions = NULL;
  static Texta geneSymbols = NULL;
  
  textCreateClear (descriptions,100);
  textCreateClear (geneSymbols,100);
  tokens = textFieldtokP (transcriptName,"|");
  for (i = 0; i < arrayMax (tokens); i++) {
    testKX.transcriptName = hlr_strdup (textItem (tokens,i));
    if (!arrayFind (kgXrefs,&testKX,&index,(ARRAYORDERF)sortKgXrefsByTranscriptName)) {
      die ("Expected to find KgXref: %s",testKX.transcriptName);
    }
    currKX = arrp (kgXrefs,index,KgXref);
    if (currKX->refseqDescription[0] != '\0') {
      textAdd (descriptions,currKX->refseqDescription);
    }
    if (currKX->geneSymbol[0] != '\0') {
      textAdd (geneSymbols,currKX->geneSymbol);
    }
    hlr_free (testKX.transcriptName);
  }
  textDestroy (tokens);
  textUniqKeepOrder (descriptions);
  textUniqKeepOrder (geneSymbols);
  *geneSymbol = hlr_strdup (convert2string (geneSymbols));
  *description = hlr_strdup (convert2string (descriptions));
}

