#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/types.h>

#include <bios/confp.h>
#include <bios/log.h>
#include <bios/format.h>
#include <bios/bowtieParser.h>

#include "bp.h"


int main (int argc, char *argv[])
{
  Array breakPoints;
  BreakPoint *currBP;
  BreakPointRead *currBPR;
  int i,j;
  Stringa buffer,cmd;
  FILE *fp;
  BowtieQuery *currBQ;
  Texta invalidJunctions;
  int index;
  char *str = NULL;

  config *conf;

  if ((conf = confp_open(getenv("FUSIONSEQ_CONFPATH"))) == NULL)
    return EXIT_FAILURE;

  buffer = stringCreate (100);
  stringPrintf (buffer,"bpJunctionReads_%d.fa",getpid ());
  fp = fopen (string (buffer),"w");
  if (fp == NULL) {
    die ("Unable to open file: %s",string (buffer));
  }
  bp_init ("-");
  breakPoints = bp_getBreakPoints ();
  for (i = 0; i < arrayMax (breakPoints); i++) {
    currBP = arrp (breakPoints,i,BreakPoint);
    for (j = 0; j < arrayMax (currBP->breakPointReads); j++) {
      currBPR = arrp (currBP->breakPointReads,j,BreakPointRead);
      fprintf (fp,">%s|%s\n%s\n",currBP->tileCoordinate1,currBP->tileCoordinate2,currBPR->read);
    }
  }
  bp_deInit ();
  fclose (fp);

  invalidJunctions = textCreate (1000);
  cmd = stringCreate (100);
  stringPrintf (cmd,"bowtie --quiet -p 4 -n 0 -f %s/%s/%s %s",
		confp_get(conf, "BOWTIE_INDEXES"),
		confp_get(conf, "BOWTIE_GENOME"), 
		confp_get(conf, "BOWTIE_GENOME"), 
		string (buffer));
  bowtieParser_initFromPipe (string (cmd));
  while (currBQ = bowtieParser_nextQuery ()) {
    textAdd (invalidJunctions,currBQ->sequenceName);
  }
  bowtieParser_deInit ();
  stringPrintf (cmd,"rm -f %s",string (buffer));
  hlr_system (string (cmd),0);

  arraySort (invalidJunctions,(ARRAYORDERF)arrayStrcmp);
  arrayUniq (invalidJunctions,NULL,(ARRAYORDERF)arrayStrcmp);
  for (i = 0; i < arrayMax (breakPoints); i++) {
    currBP = arrp (breakPoints,i,BreakPoint);
    stringPrintf (buffer,"%s|%s",currBP->tileCoordinate1,currBP->tileCoordinate2);
    strReplace (&str,string (buffer));
    if (!arrayFind (invalidJunctions,&str,&index,(ARRAYORDERF)arrayStrcmp)) {
      puts (bp_writeBreakPoint (currBP));
    }
  }
  stringDestroy (buffer);
  stringDestroy (cmd);
  confp_close(conf);

  return EXIT_SUCCESS;
}
