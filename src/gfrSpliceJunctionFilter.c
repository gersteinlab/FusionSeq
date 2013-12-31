#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/types.h>
#include <time.h>

#include <bios/confp.h>
#include <bios/log.h>
#include <bios/format.h>
#include <bios/bowtieParser.h>
#include <bios/fasta.h>
#include <bios/intervalFind.h>
#include <bios/blatParser.h>

#include "gfr.h"
#include "util.h"

#define MAX_FRACTION_SPLICES 0.05

int sortGfrById( GfrEntry* a, GfrEntry* b) {
  return strcmp( a->id, b->id);
}

int sortGfrByDASPER( GfrEntry* a, GfrEntry* b) {
  return (int) ( b->DASPER*1e06 - a->DASPER*1e06);
}

int main (int argc, char *argv[])
{
  GfrEntry *currGE;
  int i,j,l;
  Stringa cmd;
  FILE *freads;
  Array gfrEntries;
  int count;
  int countRemoved;
  int readSize1,readSize2;
  BlatQuery *blQ = NULL;
  
  if( argc != 2 ) {
    usage("%s <splice_junction_library>\nNB: the splice junction library should be in 2bit format.\nEx: %s ucsc_nh_sj75.2bit", argv[0], argv[0]);
  }
  char* spliceJunctionLibrary = argv[1];

  gfr_init ("-");
  gfrEntries = arrayCreate( 100, GfrEntry );
  gfrEntries =  gfr_parse ();
  if (arrayMax (gfrEntries) == 0){
    puts (gfr_writeHeader ());
    gfr_deInit ();
    return 0;
  }

  // creating the fasta files with the reads 
  Stringa readsFA = stringCreate( 100 ); 
 
  stringPrintf( readsFA, "%d_reads.fa", (int) getpid());
  freads = fopen ( string(readsFA) ,"w");
  if (freads == NULL) {
    die ("Unable to open file: %s",string (readsFA));
  }     
  cmd = stringCreate (100);
  puts (gfr_writeHeader ());
  count = 0;
  countRemoved = 0;
  for (i = 0; i < arrayMax (gfrEntries); i++) {
    currGE = arrp (gfrEntries,i,GfrEntry);  
    if (strEqual( currGE->fusionType, "read-through")) {
      continue;
    }
    // creating one fasta files with the reads
    if (arrayMax(currGE->readsTranscript1) != arrayMax(currGE->readsTranscript2))
      die("Error: different number of inter-transcript reads %d vs. %d", 
          arrayMax( currGE->readsTranscript1), 
          arrayMax(currGE->readsTranscript2));
    // writing the reads into file
    for (l = 0; l < arrayMax (currGE->readsTranscript1); l++) {      
      char* currRead1 = hlr_strdup( textItem (currGE->readsTranscript1,l)); // read1
      char* currRead2 = hlr_strdup( textItem (currGE->readsTranscript2,l)); // read2
      fprintf( freads, ">%s_1_%d\n%s\n>%s_2_%d\n%s\n", currGE->id, l+1, currRead1, currGE->id, l+1, currRead2 );
      readSize1 = strlen( currRead1 );
      readSize2 = strlen( currRead2 );
      if(readSize1 != readSize2 ) die("The two reads have different lengths: 1:%d vs 2:%d", readSize1, readSize2);
      hlr_free( currRead1 );
      hlr_free( currRead2 );
    }
  }
  fclose( freads ); 
  freads=NULL;
  
  
  //blat of reads against the splice junction library
  stringPrintf( cmd, "blat -t=dna %s %s stdout", spliceJunctionLibrary, string (readsFA) );
  
  arraySort(gfrEntries, (ARRAYORDERF) sortGfrById);
  blatParser_initFromPipe(string(cmd));
  Texta toRemove = textCreate( 10 );
  while (blQ = blatParser_nextQuery()) {
    Texta tok  = textFieldtokP( blQ->qName, "_");
    Stringa readID = stringCreate( 20 );
    for (i = 0; i < (arrayMax(tok)-2); i++) {
      stringAppendf( readID, "%s%s", textItem( tok, i), i < (arrayMax(tok)-3) ? "_" : "" );
    }
    while(  l<arrayMax( toRemove) ) {
      if( strEqual( string(readID), textItem( toRemove, l ) ) ) 
	break;
      l++;
    }
    if( l<arrayMax( toRemove ) )
      continue;

    GfrEntry* gfrTest;
    AllocVar( gfrTest );    
    gfrTest->id = hlr_strdup( string(readID) );
    int index=-1;
    arrayFind( gfrEntries, gfrTest, &index, (ARRAYORDERF) sortGfrById);          
    GfrEntry* currGE = arrp( gfrEntries, index, GfrEntry );
    int numHits=0;
    for( j=0; j<arrayMax( blQ->entries ); j++ ) {
      PslEntry *currE = arrp( blQ->entries, j, PslEntry );           
      Texta readPos = textFieldtokP( currE->tName, "|");
      int readNum = atoi( textItem( tok, arrayMax(tok)-2 ) );
      if( readNum == 1 ) { // found hit for read1 
	if( strEqual( currGE->chromosomeTranscript2, textItem( readPos, 0 )) && 
	    currGE->startTranscript2 <= atoi(  textItem( readPos, 1 )) &&
	    currGE->endTranscript2 >= atoi( textItem( readPos, 2 )) ) { //found a proper hit in the splice junction
	  numHits++;
	}
      }
      if( readNum == 2 ) {// found hit for read2
	if( strEqual( currGE->chromosomeTranscript1, textItem( readPos, 0 )) && 
	    currGE->startTranscript1 <= atoi( textItem(readPos, 1 )) &&
	    currGE->endTranscript1 >= atoi(textItem( readPos, 2 )) ) { //found a proper hit in the splice junction
	  numHits++;
	}
      }      
      textDestroy( readPos);
    }    
    if( numHits > arrayMax( currGE->interReads ) * MAX_FRACTION_SPLICES )  {
      textAdd( toRemove, currGE->id);
    }
    
    hlr_free( gfrTest->id );
    freeMem( gfrTest );
    textDestroy( tok ); 
    stringDestroy( readID );
  }
  blatParser_deInit();
  arraySort( toRemove, (ARRAYORDERF) arrayStrcmp );
  arraySort( gfrEntries, (ARRAYORDERF) sortGfrByDASPER );
  for (i = 0; i < arrayMax(gfrEntries); i++) {
    GfrEntry *currGE = arrp( gfrEntries, i, GfrEntry);
    if( arrayFind( toRemove, &currGE->id, NULL, (ARRAYORDERF) arrayStrcmp ) ) {
      countRemoved++;
    } else {
      puts( gfr_writeGfrEntry( currGE ) );
      count++;
    }
  }

  // removing temporary files
  stringPrintf (cmd,"rm -rf %s", string(readsFA) );
  hlr_system( string(cmd) , 0);      

  gfr_deInit ();
  arrayDestroy ( gfrEntries );
  arrayDestroy ( toRemove );
  stringDestroy( cmd );
  stringDestroy( readsFA );
  warn ("%s_spliceJunctionLibrary: %s",argv[0],spliceJunctionLibrary);
  warn ("%s_numRemoved: %d",argv[0],countRemoved);  
  warn ("%s_numGfrEntries: %d",argv[0],count);
  return 0;
}

