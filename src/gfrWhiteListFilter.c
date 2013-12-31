#include <stdio.h>

#include <bios/log.h>
#include <bios/format.h>
#include <bios/linestream.h>

#include "gfr.h"

typedef struct {
  char* gene1;
  char* gene2;
} WLEntry;

static int sortWhiteListByName1 (WLEntry *a, WLEntry *b) 
{
  int res = strcmp ( a->gene1, b->gene1);
  if( res==0 ) res = strcmp ( a->gene2, b->gene2 );
  return (res); //(strcmp ( a->gene1, b->gene1));
}

int main (int argc, char *argv[])
{
  GfrEntry *currGE;
  WLEntry *currWLE;
  WLEntry currQuery;
  FILE *fp;
  char *line;
  int count;

  int index;
  WordIter w;
  Array whiteList = arrayCreate(20, WLEntry);

  if (argc != 2) {
    usage ("%s <whiteList.txt>",argv[0]);
  }  
  fp = fopen( argv[1], "r" );
  
  if( !fp )  die("Unable to open file: %s", argv[1]);
  // reading whitelist file
  LineStream ls = ls_createFromFile( argv[1] );
  while( line = ls_nextLine(ls) ) {
    w = wordIterCreate( line, "\t", 1);
    currWLE = arrayp( whiteList, arrayMax(whiteList), WLEntry);
    currWLE->gene1 = hlr_strdup ( wordNext(w) );
    currWLE->gene2 = hlr_strdup ( wordNext(w) );    
    wordIterDestroy(w);
  }
  fclose(fp);
  arraySort( whiteList, (ARRAYORDERF) sortWhiteListByName1);

  // beginFiltering
  count = 0;
  gfr_init ("-");
  puts (gfr_writeHeader ());
  while (currGE = gfr_nextEntry ()) { // reading the gfr
    // creating a new query to the white list
    currQuery.gene1 = currGE->geneSymbolTranscript1;
    currQuery.gene2 = currGE->geneSymbolTranscript2;
    // searching against read_1/read_2
    int res = arrayFind( whiteList, &currQuery, &index, (ARRAYORDERF) sortWhiteListByName1);  
    if( res ) { // found, write the instance to stdout, update the counts 
      puts (gfr_writeGfrEntry (currGE));
      count++;
    } else { //not found: then searching against read_2/read_1
      currQuery.gene1 = currGE->geneSymbolTranscript2;
      currQuery.gene2 = currGE->geneSymbolTranscript1;
      
      res =  arrayFind( whiteList, &currQuery, 
			&index, (ARRAYORDERF) sortWhiteListByName1 );
      
      if( res ) { // found, write the instance to stdout, update the counts
	puts (gfr_writeGfrEntry (currGE));
	count++;	
      } 
    }
  }	           
  gfr_deInit ();
  arrayDestroy( whiteList );
  warn ("%s_WhiteListFilter: %s",argv[0], argv[1]);
  warn ("%s_numGfrEntries: %d",argv[0],count);
  return 0;
}

