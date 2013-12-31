#include <stdio.h>
#include <sys/types.h>
#include <unistd.h>

#include <bios/confp.h>
#include <bios/log.h>
#include <bios/format.h>
#include <bios/bowtieParser.h>
#include <bios/fasta.h>
#include <bios/intervalFind.h>
#include <bios/blatParser.h>

#include "gfr.h"
#include "util.h"

#define MAX_FRACTION_HOMOLOGOUS 0.05
#define MAX_OVERLAP_ALLOWED 0.75

int main (int argc, char *argv[])
{
  GfrEntry *currGE;
  int i,j,l;
  Stringa cmd;
  FILE *freads;
  Array gfrEntries;
  int ribosomalCount;
  int count;
  int countRemoved;
  int readSize1,readSize2;
  BlatQuery *blQ=NULL;

  config *conf;

  if ((conf = confp_open(getenv("FUSIONSEQ_CONFPATH"))) == NULL)
    return EXIT_FAILURE;
  
  gfr_init ("-");
  gfrEntries = arrayCreate( 100, GfrEntry );
  gfrEntries =  gfr_parse ();
  if (arrayMax (gfrEntries) == 0){
    puts (gfr_writeHeader ());
    gfr_deInit ();
    return 0;
  }

  cmd = stringCreate (100);
  count = 0;
  countRemoved = 0;
  puts (gfr_writeHeader ());
  j = 0;
  for (i = 0; i < arrayMax (gfrEntries); i++) {
    currGE = arrp (gfrEntries,i,GfrEntry);
    ribosomalCount = 0;
    
    // creating one fasta files with the reads
    Stringa readsFA = stringCreate( 100 ); 
 
    // creating the fasta files with the reads 
    stringPrintf( readsFA, "%s_reads.fa", currGE->id);
    freads = fopen ( string(readsFA) ,"w");
    if (freads == NULL) {
      die ("Unable to open file: %s",string (readsFA));
    }     
    if (arrayMax(currGE->readsTranscript1) != arrayMax(currGE->readsTranscript2))
      die("Error: different number of inter-transcript reads %d vs. %d", 
          arrayMax(currGE->readsTranscript1),
	  arrayMax( currGE->readsTranscript2));

    // writing the reads into file
    for (l = 0; l < arrayMax (currGE->readsTranscript1); l++) {      
      char* currRead1 = hlr_strdup( textItem (currGE->readsTranscript1,l)); // read1
      char* currRead2 = hlr_strdup( textItem (currGE->readsTranscript2,l)); // read2
      fprintf( freads, ">%d/1\n%s\n>%d/2\n%s\n", l+1, currRead1, l+1, currRead2 );
      readSize1 = strlen( currRead1 );
      readSize2 = strlen( currRead2 );
      if(readSize1 != readSize2 )
        die("The two reads have different lengths: 1:%d vs 2:%d", readSize1, readSize2);
      hlr_free( currRead1 );
      hlr_free( currRead2 );
    }
    fclose( freads ); 
    freads=NULL;
    stringDestroy( readsFA );

    //blat of reads against the ribosomal genes
    // 2*stepSize + tileSize - 1 = min_num_nt_to_trigger_alignment
    // tileSize = 15
    // int stepSize = (int)((readSize1*MAX_OVERLAP_ALLOWED + 1 + 11) * 0.5);
    stringPrintf(cmd, "blat -t=dna -q=dna -out=psl -fine -repMatch=1000000 -tileSize=15 %s/%s %s_reads.fa stdout", 
                 confp_get(conf, "RIBOSOMAL_DIR"), 
		 confp_get(conf, "RIBOSOMAL_FILENAME"), 
		 currGE->id);
    // reading the results of blast from Pipe
    blatParser_initFromPipe( string(cmd) );
    while( blQ = blatParser_nextQuery() ) {
      int nucleotideOverlap = getNucleotideOverlap ( blQ );
      if (nucleotideOverlap > (((double)readSize1)*MAX_OVERLAP_ALLOWED)) {
	ribosomalCount++;
      } 
    }
    blatParser_deInit();
    if (( (double)ribosomalCount / ( (double)currGE->numInter * 2.0 ) ) <= MAX_FRACTION_HOMOLOGOUS) {       
      // writing the gfrEntry
      puts (gfr_writeGfrEntry (currGE));
      count++;
    } else {
      countRemoved++;
    }
    // removing temporary files
    stringPrintf (cmd,"rm -rf %s_reads.fa", currGE->id );
    hlr_system( string(cmd) , 0);      
  }

  gfr_deInit ();
  arrayDestroy ( gfrEntries );
  stringDestroy( cmd );
  warn ("%s_numRemoved: %d",argv[0],countRemoved);  
  warn ("%s_numGfrEntries: %d",argv[0],count);

  confp_close(conf);
  return 0;
}

