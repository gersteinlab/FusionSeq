#include <stdlib.h>
#include <unistd.h>

#include <bios/confp.h>
#include <bios/log.h>
#include <bios/format.h>
#include <bios/bowtieParser.h>
#include <bios/common.h>
#include <bios/fasta.h>

#include "bp.h"

typedef struct {
  char* chromosome1;  
  int start1;
  int end1;
  char* chromosome2;
  int start2;  
  int end2; 
} BreakPointJunction;

static config *Conf = NULL;

static char* getBreakPointSequence (char *tileCoordinate1, char *tileCoordinate2)
{
  Stringa buffer;
  Stringa targetsFile;
  FILE *fp;
  Array targetSeqs;
  int i;
  Seq *currSeq;
  static Stringa sequence = NULL;

  buffer = stringCreate (100);
  targetsFile = stringCreate (100);
  stringPrintf (targetsFile,"targets_%d.txt",getpid ());
  if (!(fp = fopen (string (targetsFile),"w")) ){
    die ("Unable to open target file: %s",string (targetsFile));
  }
  fprintf (fp,"%s\n%s",tileCoordinate1,tileCoordinate2);
  fclose (fp);
  stringPrintf (buffer,"%s %s/%s stdout -noMask -seqList=%s",
		confp_get(Conf, "BLAT_TWO_BIT_TO_FA"),
		confp_get(Conf, "BLAT_DATA_DIR"),
		confp_get(Conf, "BLAT_TWO_BIT_DATA_FILENAME"),
		string (targetsFile));
  fasta_initFromPipe (string (buffer));
  targetSeqs = fasta_readAllSequences (0);
  fasta_deInit ();
  if (arrayMax (targetSeqs) != 2) {
    die ("Expected only two target sequences");
  } 
  stringCreateClear (sequence,100);
  for (i = 0; i < arrayMax (targetSeqs); i++) {
    currSeq = arrp (targetSeqs,i,Seq);
    stringAppendf (sequence,"%s",currSeq->sequence);
    hlr_free (currSeq->name);
    hlr_free (currSeq->sequence);
  }
  arrayDestroy (targetSeqs);
  stringPrintf (buffer,"rm -rf %s",string (targetsFile));
  hlr_system (string (buffer),0);
  stringDestroy (targetsFile);
  stringDestroy (buffer);
  return string (sequence);
}

static void getCoordinates( char* str, BreakPointJunction* bpJunction) 
{
  Texta txt = textFieldtokP(str, ":|-"); 
  bpJunction->chromosome1 = hlr_strdup( textItem( txt, 0 ));
  bpJunction->start1 = atoi(textItem( txt, 1 ) );
  bpJunction->end1 = atoi( textItem( txt, 2 ) );
  bpJunction->chromosome2 = hlr_strdup( textItem( txt, 3 ) );
  bpJunction->start2 = atoi(textItem( txt, 4 ) );
  bpJunction->end2 = atoi( textItem( txt, 5 ) );  
  textDestroy( txt );
}
static int sortHits( BowtieEntry *a, BowtieEntry *b) 
{
  BreakPointJunction* BPJa, *BPJb;  
  AllocVar( BPJa );
  AllocVar( BPJb );
  getCoordinates( a->chromosome, BPJa);
  getCoordinates( b->chromosome, BPJb);
  int diff =  ( (BPJa->start1 - BPJb->start1) + (BPJa->start2 - BPJb->start2) );
  hlr_free( BPJa->chromosome1 );
  hlr_free( BPJa->chromosome2 );
  hlr_free( BPJb->chromosome1 );
  hlr_free( BPJb->chromosome2 );
  freeMem( BPJa );
  freeMem( BPJb );
  return diff;
}

static int sortBreakPointJunctions( BreakPointJunction *a, BreakPointJunction *b) 
{
  int diff =  ( (a->start1 - b->start1) + (a->start2 - b->start2) );
  return diff;
}

int main (int argc, char *argv[])
{
  BowtieQuery *currBQ;
  BowtieEntry *currBE;
  Array bowtieQueries;
  Array breakPointJunctions;
  BreakPointJunction *currBPJ;
  int i,j;

  if ((Conf = confp_open(getenv("FUSIONSEQ_CONFPATH"))) == NULL)
    return EXIT_FAILURE;

  bowtieParser_initFromFile ("-");
  bowtieQueries = bowtieParser_getAllQueries ();
  bowtieParser_deInit (); 
  breakPointJunctions = arrayCreate (10000,BreakPointJunction);
  for (i = 0; i < arrayMax (bowtieQueries); i++) {
    currBQ = arrp (bowtieQueries,i,BowtieQuery);
    currBPJ = arrayp (breakPointJunctions,arrayMax(breakPointJunctions), BreakPointJunction);
    if( arrayMax( currBQ->entries ) == 1 ) {
      // print junction, unique
      getCoordinates( currBQ->sequenceName, currBPJ);      
    } else { // junction not unique 
      arraySort( currBQ->entries, (ARRAYORDERF) sortHits);
      for( j=0; j<arrayMax( currBQ->entries ); j++ ) {
	currBE = arrp (currBQ->entries,j,BowtieEntry);
	BreakPointJunction *bpJ;
	if( j==0 )
	  getCoordinates( currBE->chromosome, currBPJ);
	else { // to check if close enough
	  AllocVar( bpJ );
	  getCoordinates( currBE->chromosome, bpJ);
	  currBPJ->end2 = bpJ->end2;
	  hlr_free(bpJ->chromosome1);
	  hlr_free(bpJ->chromosome2);
	  freeMem( bpJ);
	}
      }
    }
  }
  arraySort( breakPointJunctions, (ARRAYORDERF) sortBreakPointJunctions);
  arrayUniq( breakPointJunctions, NULL, (ARRAYORDERF) sortBreakPointJunctions );

  for( i=0; i<arrayMax( breakPointJunctions); i++ ) {
    currBPJ = arrp (breakPointJunctions, i, BreakPointJunction);

    Stringa tileCoordinate1 = stringCreate( 100 );
    stringPrintf( tileCoordinate1, "%s:%d-%d", currBPJ->chromosome1, currBPJ->start1, currBPJ->end1 );
    Stringa tileCoordinate2 = stringCreate( 100 );
    stringPrintf( tileCoordinate2, "%s:%d-%d", currBPJ->chromosome2, currBPJ->start2, currBPJ->end2 );
    char* breakPointSequence = getBreakPointSequence ( string(tileCoordinate1), string(tileCoordinate2) );

    printf(">%s:%d-%d|%s:%d-%d\n%s\n", 
	   currBPJ->chromosome1, 
	   currBPJ->start1, 
	   currBPJ->end1, 
	   currBPJ->chromosome2, 
	   currBPJ->start2, 
	   currBPJ->end2, 
	   breakPointSequence);
    stringDestroy( tileCoordinate1 );
    stringDestroy( tileCoordinate2 );
  }
  warn("bpClustering: Number of breakpoint junctions:\t%d",--i);
  confp_close(Conf);

  return EXIT_SUCCESS;
}
