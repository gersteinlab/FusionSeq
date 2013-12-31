#include <stdio.h>

#include <bios/log.h>
#include <bios/format.h>
#include <bios/linestream.h>

#include "gfr.h"

typedef struct {
  char* gene1;
  char* gene2;
} BLEntry;

static int sortBlackListByName1 (BLEntry *a, BLEntry *b) 
{
  int res = strcmp ( a->gene1, b->gene1);
  if( res==0 ) res = strcmp ( a->gene2, b->gene2 );
  return (res); //(strcmp ( a->gene1, b->gene1));
}

int main (int argc, char *argv[])
{
  GfrEntry *currGE;
  BLEntry *currBLE;
  BLEntry currQuery;
  FILE *fp;
  char *line;
  int count;
  int countRemoved;
  
  int index;
  WordIter w;
  Array blackList = arrayCreate(20, BLEntry);

  if (argc != 2) {
    usage ("%s <blackList.txt>",argv[0]);
  }  
  fp = fopen( argv[1], "r" );
  
  if( !fp )  die("Unable to open file: %s", argv[1]);
  // reading blacklist file
  LineStream ls = ls_createFromFile( argv[1] );
  while( line = ls_nextLine(ls) ) {
    w = wordIterCreate( line, "\t", 1);
    currBLE = arrayp( blackList, arrayMax(blackList), BLEntry);
    currBLE->gene1 = hlr_strdup ( wordNext(w) );
    currBLE->gene2 = hlr_strdup ( wordNext(w) );    
    wordIterDestroy(w);
  }
  fclose(fp);
  arraySort( blackList, (ARRAYORDERF) sortBlackListByName1);

  // beginFiltering
  count = 0;
  countRemoved = 0;
  gfr_init ("-");
  puts (gfr_writeHeader ());
  while (currGE = gfr_nextEntry ()) { // reading the gfr
    // creating a new query to the black list
    currQuery.gene1 = currGE->geneSymbolTranscript1;
    currQuery.gene2 = currGE->geneSymbolTranscript2;
    // searching against read_1/read_2
    int res = arrayFind( blackList, &currQuery, 
			 &index,  (ARRAYORDERF) sortBlackListByName1);  
    
    if( !res ) { // not found, then searching against read_2/read_1
      currQuery.gene1 = currGE->geneSymbolTranscript2;
      currQuery.gene2 = currGE->geneSymbolTranscript1;
      
      res =  arrayFind( blackList, &currQuery, 
			&index, (ARRAYORDERF) sortBlackListByName1 );
      
      if( !res ) { // not found, write the instance to stdout, update the counts
	puts (gfr_writeGfrEntry (currGE));
	count++;	
      } else { // found: read2/read1
	countRemoved++;
      }	
    } else { //found: read1/read2
      countRemoved++;
    }
  }	           
  gfr_deInit ();
  arrayDestroy( blackList );
  warn ("%s_BlackListFilter: %s",argv[0], argv[1]);
  warn ("%s_numRemoved: %d",argv[0],countRemoved);
  warn ("%s_numGfrEntries: %d",argv[0],count);
  return 0;
}

