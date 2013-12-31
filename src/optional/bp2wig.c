#include <stdio.h>

#include <bios/log.h>
#include <bios/format.h>

#include "bp.h"



static int sortBreakPointsByTileCoordinate1 (BreakPoint *a, BreakPoint *b) 
{
  return strcmp (a->tileCoordinate1,b->tileCoordinate1);
}



static int sortBreakPointsByTileCoordinate2 (BreakPoint *a, BreakPoint *b) 
{
  return strcmp (a->tileCoordinate2,b->tileCoordinate2);
}



int main (int argc, char *argv[])
{
  Array breakPoints;
  BreakPoint *currBP,*nextBP;
  int i,j;
  char *pos;
  Stringa buffer;
  FILE *fp;
  char *tileCopy;
  int count;

  if (argc != 2) {
    usage ("%s <fileName.bp>",argv[0]);
  }

  bp_init (argv[1]);
  breakPoints = bp_getBreakPoints ();
  bp_deInit ();
  
  buffer = stringCreate (100);
  pos = strchr (argv[1],'.');
  *pos = '\0';
  stringPrintf (buffer,"%s_breakPointsTranscript%d.wig",argv[1],strstr(argv[1], "_AB") ? 1 : 2); // NB: TRANSCRIPT 1 always corresponds to the leftmost gene as in the BP results; hence, if the orientation is BA, the TileCoordinate1 corresponds to transcript2
  fp = fopen (string (buffer),"w");
  if (fp == NULL) {
    die ("Unable to open file: %s",string (buffer));
  }
  arraySort (breakPoints,(ARRAYORDERF)sortBreakPointsByTileCoordinate1);
  tileCopy = NULL;
  i = 0; 
  while (i < arrayMax (breakPoints)) {
    currBP = arrp (breakPoints,i,BreakPoint);
    count = arrayMax (currBP->breakPointReads);
    j = i + 1;
    while (j < arrayMax (breakPoints)) {
      nextBP = arrp (breakPoints,j,BreakPoint);
      if (strEqual (currBP->tileCoordinate1,nextBP->tileCoordinate1)) {
        count += arrayMax (nextBP->breakPointReads);
      }
      else {
        break;
      }
      j++;
    }
    i = j;
    if (tileCopy == NULL) {
      strReplace (&tileCopy,currBP->tileCoordinate1);
      pos = strchr (tileCopy,':');
      *pos = '\0';
      fprintf (fp,"track type=wiggle_0 name=\"Breakpoints %s_%d\"\n",argv[1], strstr(argv[1], "_AB") ? 1 : 2); // see comment above
      fprintf (fp,"variableStep chrom=%s span=1\n",tileCopy);
    }
    pos = strchr (currBP->tileCoordinate1,'-');
    fprintf (fp,"%d\t%d\n",atoi (pos + 1),count);
  }
  fprintf (fp, "%d\t%d\n",atoi (pos + 1)+1,0); //add 0 for scaling issues in the Genome Browser
  fclose (fp);
  
  stringPrintf (buffer,"%s_breakPointsTranscript%d.wig",argv[1],strstr(argv[1], "_AB") ? 2 : 1);  // see comment above
  fp = fopen (string (buffer),"w");
  if (fp == NULL) {
    die ("Unable to open file: %s",string (buffer));
  }
  arraySort (breakPoints,(ARRAYORDERF)sortBreakPointsByTileCoordinate2);
  tileCopy = NULL;
  i = 0; 
  while (i < arrayMax (breakPoints)) {
    currBP = arrp (breakPoints,i,BreakPoint);
    count = arrayMax (currBP->breakPointReads);
    j = i + 1;
    while (j < arrayMax (breakPoints)) {
      nextBP = arrp (breakPoints,j,BreakPoint);
      if (strEqual (currBP->tileCoordinate2,nextBP->tileCoordinate2)) {
        count += arrayMax (nextBP->breakPointReads);
      }
      else {
        break;
      }
      j++;
    }
    i = j;
    if (tileCopy == NULL) {
      strReplace (&tileCopy,currBP->tileCoordinate2);
      pos = strchr (tileCopy,':');
      *pos = '\0';
      fprintf (fp,"track type=wiggle_0 name=\"Breakpoints %s_%d\"\n",argv[1],strstr(argv[1], "_AB") ? 2 : 1);  // see comment above
      fprintf (fp,"variableStep chrom=%s span=1\n",tileCopy);
    }
    pos = strchr (currBP->tileCoordinate2,':');
    fprintf (fp,"%d\t%d\n",atoi (pos + 1),count);
  }    
  fprintf (fp, "%d\t%d\n",atoi (pos + 1)+1,0); //add 0 for scaling issues in the Genome Browser
  fclose (fp);
  return 0;
}

 
