#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>

#include <bios/confp.h>
#include <bios/log.h>
#include <bios/format.h>
#include <bios/fasta.h>
#include <bios/stringUtil.h>

#include "bp.h"

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

static int sortBreakPointsByTargetAndOffset (BreakPoint *a, BreakPoint *b) 
{
  int diff;
  diff = strcmp (a->tileCoordinate1, b->tileCoordinate1) + strcmp(a->tileCoordinate2, b->tileCoordinate2);
  return diff;
}

int main (int argc, char *argv[])
{
  Array breakPoints;
  BreakPoint *currBP;
  int i;
  char *breakPointSequence;

  if ((Conf = confp_open(getenv("FUSIONSEQ_CONFPATH"))) == NULL)
    return EXIT_FAILURE;

  bp_init ("-");
  breakPoints = bp_getBreakPoints ();
  arraySort (breakPoints,(ARRAYORDERF)sortBreakPointsByTargetAndOffset);
  
  for (i = 0; i < arrayMax (breakPoints); i++) {
    currBP = arrp (breakPoints,i,BreakPoint);
    breakPointSequence = getBreakPointSequence (currBP->tileCoordinate1,currBP->tileCoordinate2);
    printf( ">%s|%s\n%s\n", currBP->tileCoordinate1, currBP->tileCoordinate2, breakPointSequence);
    warn(">%s|%s\n%s", 
	 currBP->tileCoordinate1, 
	 currBP->tileCoordinate2, 
	 subString(breakPointSequence, 10, strlen(breakPointSequence)-10));
  }
  bp_deInit();
  confp_close(Conf);

  return EXIT_SUCCESS;
}
