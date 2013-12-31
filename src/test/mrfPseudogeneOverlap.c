#include <stdio.h>

#include <bios/log.h>
#include <bios/format.h>
#include <bios/intervalFind.h>
#include <mrf/mrf.h>



static int readIsContainedInPseudogene (MrfRead *currMrfRead)
{
  Array intervals;
  MrfBlock *currMrfBlock;
  int i;
  Interval *currInterval;

  if (arrayMax (currMrfRead->blocks) != 1) {
    die ("Expected only one alignment block for each end of paired-end read!");
  }
  currMrfBlock = arrp (currMrfRead->blocks,0,MrfBlock);
  intervals = intervalFind_getOverlappingIntervals (currMrfBlock->targetName,currMrfBlock->targetStart,currMrfBlock->targetEnd);
  for (i = 0; i < arrayMax (intervals); i++) {
    currInterval = arru (intervals,i,Interval*);
    if (currInterval->start <= currMrfBlock->targetStart && currMrfBlock->targetEnd <= currInterval->end) {
      return 1;
    }
  }
  return 0;
}



int main (int argc, char *argv[])
{
  MrfEntry *currMrfEntry;
  Stringa buffer;
  FILE *fp1,*fp2;
  int i;
 
  if (argc < 3) {
    usage ("%s <prefix> <pseudogene.annotation> [-eland]",argv[0]);
  }
  buffer = stringCreate (100);
  stringPrintf (buffer,"%s_1.fa",argv[1]);
  fp1 = fopen (string (buffer),"w");
  stringPrintf (buffer,"%s_2.fa",argv[1]);
  fp2 = fopen (string (buffer),"w");
  if (fp1 == NULL || fp2 == NULL) {
    die ("Unable to open output files");
  }
  i = 0;
  intervalFind_addIntervalsToSearchSpace (argv[2],0);
  mrf_init ("-");
  while (currMrfEntry = mrf_nextEntry ()) {
    if (readIsContainedInPseudogene (&currMrfEntry->read1) || readIsContainedInPseudogene (&currMrfEntry->read2)) {
      i++;
      if( argc==4 ) {
	if( strEqual( argv[3], "-eland") ) {
	  stringPrintf( buffer, "%s:1:1:10:10#%d", argv[1],i);
	}
      } else {
	stringPrintf( buffer, "%s_%d", argv[1], i);
      } 
      fprintf (fp1,">%s/1\n%s\n",string(buffer),currMrfEntry->read1.sequence);
      fprintf (fp2,">%s/2\n%s\n",string(buffer),currMrfEntry->read2.sequence);
    }
  }
  mrf_deInit ();
  stringDestroy (buffer);
  fclose (fp1);
  fclose (fp2);
  return 0;
}
