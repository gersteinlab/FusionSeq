#include <bios/log.h>
#include <bios/format.h>
#include <bios/linestream.h>

#include "bp.h"



static LineStream lsBp = NULL;



void bp_init (const char* fileName)
{
  lsBp = ls_createFromFile (fileName);
}



void bp_deInit (void)
{
  ls_destroy (lsBp);
}



Array bp_getBreakPoints (void)
{
  Array breakPoints;
  BreakPoint *currBP;
  char *line,*pos,*token;
  Texta tokens;
  BreakPointRead *currBPR;
  WordIter w;

  breakPoints = arrayCreate (100,BreakPoint);
  while (line = ls_nextLine (lsBp)) {
    if (line[0] == '\0') {
      continue;
    }
    tokens = textStrtok (line,",");
    currBP = arrayp (breakPoints,arrayMax (breakPoints),BreakPoint);
    currBP->tileCoordinate1 = hlr_strdup (textItem (tokens,0));
    currBP->tileCoordinate2 = hlr_strdup (textItem (tokens,1)); 
    currBP->breakPointReads = arrayCreate (100,BreakPointRead);
    w = wordIterCreate (textItem (tokens,2),"|",0);
    while (token = wordNext (w)) {
      pos = strchr (token,':');
      *pos = '\0';
      currBPR = arrayp (currBP->breakPointReads,arrayMax (currBP->breakPointReads),BreakPointRead);
      currBPR->offset = atoi (token);
      currBPR->read = hlr_strdup (pos + 1);
    }
    wordIterDestroy (w);
    textDestroy (tokens);
  }
  return breakPoints;
}



char* bp_writeBreakPoint (BreakPoint *currBP)
{
  static Stringa buffer = NULL;
  BreakPointRead *currBPR;
  int i;

  stringCreateClear (buffer,1000);
  stringPrintf (buffer,"%s,%s,",currBP->tileCoordinate1,currBP->tileCoordinate2);
  for (i = 0; i < arrayMax (currBP->breakPointReads); i++) {
    currBPR = arrp (currBP->breakPointReads,i,BreakPointRead);
    stringAppendf (buffer,"%d:%s%s",currBPR->offset,currBPR->read, i < arrayMax (currBP->breakPointReads) - 1 ? "|" : "");
  }
  return string (buffer);
}
