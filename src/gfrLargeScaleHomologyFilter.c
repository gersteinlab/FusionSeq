#include <stdlib.h>

#include <bios/confp.h>
#include <bios/log.h>
#include <bios/format.h>

#include "util.h"
#include "gfr.h"



static int sortKgTreeFamsByTranscriptName (KgTreeFam *a, KgTreeFam *b) 
{
  return strcmp (a->transcriptName,b->transcriptName);
}



static char* lookUpTreeFam (Array kgTreeFams, char *transcript) 
{
  KgTreeFam testKGTF;
  int index;
  int foundIt;
   
  foundIt = 0;
  testKGTF.transcriptName = hlr_strdup (transcript);
  foundIt = arrayFind (kgTreeFams,&testKGTF,&index,(ARRAYORDERF)sortKgTreeFamsByTranscriptName);
  hlr_free (testKGTF.transcriptName);
  if (foundIt) {
    return  arrp (kgTreeFams,index,KgTreeFam)->treeFamId;
  }
  return NULL;
}



static int isHomologous (Array kgTreeFams, char *transcript1, char *transcript2)
{
  Texta tokens;
  int i,j;
  char *treeFamId;
  static Texta treeFamIdsTranscript1 = NULL;
  static Texta treeFamIdsTranscript2 = NULL;

  textCreateClear (treeFamIdsTranscript1,100);
  textCreateClear (treeFamIdsTranscript2,100);
  tokens = textFieldtokP (transcript1,"|");
  for (i = 0; i < arrayMax (tokens); i++) {
    if (treeFamId = lookUpTreeFam (kgTreeFams,textItem (tokens,i))) {
      textAdd (treeFamIdsTranscript1,treeFamId);
    } 
  }
  textDestroy (tokens);
  tokens = textFieldtokP (transcript2,"|");
  for (i = 0; i < arrayMax (tokens); i++) {
    if (treeFamId = lookUpTreeFam (kgTreeFams,textItem (tokens,i))) {
      textAdd (treeFamIdsTranscript2,treeFamId);
    } 
  }
  textDestroy (tokens);
  for (i = 0; i < arrayMax (treeFamIdsTranscript1); i++) {
     for (j = 0; j < arrayMax (treeFamIdsTranscript2); j++) {
       if (strEqual (textItem (treeFamIdsTranscript1,i),textItem (treeFamIdsTranscript2,j))) {
         return 1;
       }
     }
  }
  return 0;
}



int main (int argc, char *argv[])
{
  GfrEntry *currGE;
  Array kgTreeFams;
  Stringa buffer;
  int count;
  int countRemoved;

  config *conf;

  if ((conf = confp_open(getenv("FUSIONSEQ_CONFPATH"))) == NULL)
    return EXIT_FAILURE;

  buffer = stringCreate (100);
  stringPrintf (buffer,"%s/%s",
                confp_get(conf, "ANNOTATION_DIR"), 
		confp_get(conf, "KNOWN_GENE_TREE_FAM_FILENAME"));
  kgTreeFams = util_readKnownGeneTreeFams (string (buffer));
  arraySort (kgTreeFams,(ARRAYORDERF)sortKgTreeFamsByTranscriptName);
  stringDestroy (buffer);

  count = 0;
  countRemoved = 0;
  gfr_init ("-");
  puts (gfr_writeHeader ());
  while (currGE = gfr_nextEntry ()){
    if (isHomologous (kgTreeFams,currGE->nameTranscript1,currGE->nameTranscript2)) {
      countRemoved++;
      continue;
    }
    puts (gfr_writeGfrEntry (currGE));
    count++;
  }
  gfr_deInit ();
  warn ("%s_numRemoved: %d",argv[0],countRemoved);
  warn ("%s_numGfrEntries: %d",argv[0],count);

  confp_close(conf);

  return EXIT_SUCCESS;
}

