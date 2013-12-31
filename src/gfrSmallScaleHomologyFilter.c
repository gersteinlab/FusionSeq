#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/types.h>

#include <bios/confp.h>
#include <bios/log.h>
#include <bios/format.h>
#include <bios/bowtieParser.h>
#include <bios/fasta.h>
#include <bios/intervalFind.h>
#include <bios/blastParser.h>

#include "gfr.h"
#include "util.h"

#define MAX_FRACTION_HOMOLOGOUS 0.05
#define MAX_OVERLAP_ALLOWED 0.75

static int sortBowtieQueriesBySequenceName (BowtieQuery *a, BowtieQuery *b)
{
  return strcmp (a->sequenceName,b->sequenceName);
}



static int sortFastaSequencesByName (Seq *a, Seq *b) 
{
  return strcmp (a->name,b->name);
}



static char getTranscriptNumber (char *seqName)
{
  char *pos;

  pos = strchr (seqName,'|');
  return *(pos + 1);
}

int main (int argc, char *argv[])
{
  GfrEntry *currGE;
  int i,j,k,l,index;
  Stringa buffer,cmd,fnSequencesToAlign;
  FILE *fp;
  FILE *fp1;
  FILE *fp2;
  FILE *freads1;
  FILE *freads2;
  Array gfrEntries;
  BowtieQuery *currBQ,testBQ;
  BowtieEntry *currBE;
  Texta seqNames;
  int readSize1, readSize2;
  Array bowtieQueries;
  char transcriptNumber;
  int isHomologous,homologousCount;
  int count;
  int countRemoved;
  BlatQuery *blQ;

  config *conf;

  if ((conf = confp_open(getenv("FUSIONSEQ_CONFPATH"))) == NULL)
    return EXIT_FAILURE;

  // read in the transcriptome sequences 
  Array transcriptome = arrayCreate( 20000, Seq );
  buffer = stringCreate(100);
  stringPrintf(buffer, "%s/%s", 
               confp_get(conf, "TRANSCRIPT_COMPOSITE_MODEL_DIR"), 
	       confp_get(conf, "TRANSCRIPT_COMPOSITE_MODEL_FA_FILENAME"));

  fasta_initFromFile( string(buffer) );
  transcriptome = fasta_readAllSequences( 0 );
  fasta_deInit();
  arraySort( transcriptome,(ARRAYORDERF)sortFastaSequencesByName ); 
  
  gfr_init ("-");
  gfrEntries =  gfr_parse ();
  if (arrayMax (gfrEntries) == 0){
    puts (gfr_writeHeader ());
    gfr_deInit ();
    return 0;
  }
  seqNames = textCreate (10000); 
  buffer = stringCreate (100);
  cmd = stringCreate (100);
  fnSequencesToAlign = stringCreate (100);
  stringPrintf (fnSequencesToAlign,"sequences_%d.fasta",getpid ());
  if (!(fp = fopen (string (fnSequencesToAlign),"w"))) {
    die ("Unable to open file: %s",string (fnSequencesToAlign));
  }
  for (i = 0; i < arrayMax (gfrEntries); i++) {
    currGE = arrp (gfrEntries,i,GfrEntry);
    for (j = 0; j < arrayMax (currGE->readsTranscript1); j++) {
      stringPrintf (buffer,"%s|1|%05d",currGE->id,j + 1);
      textAdd (seqNames,string (buffer));
      fprintf (fp,">%s\n%s\n",string (buffer),textItem (currGE->readsTranscript1,j));
    }
    for (j = 0; j < arrayMax (currGE->readsTranscript2); j++) {
      stringPrintf (buffer,"%s|2|%05d",currGE->id,j + 1);
      textAdd (seqNames,string (buffer));
      fprintf (fp,">%s\n%s\n",string (buffer),textItem (currGE->readsTranscript2,j));
    }
  } 
  fclose (fp);
  stringPrintf (cmd,"bowtie -a -v 3 --quiet -f %s/%s/%s %s", 
                confp_get(conf, "BOWTIE_INDEXES"), 
                confp_get(conf, "BOWTIE_COMPOSITE"), 
                confp_get(conf, "BOWTIE_COMPOSITE"), 
		string (fnSequencesToAlign) );
  bowtieParser_initFromPipe (string (cmd));
  bowtieQueries = bowtieParser_getAllQueries ();
  bowtieParser_deInit ();
  arraySort (bowtieQueries,(ARRAYORDERF)sortBowtieQueriesBySequenceName);

  count = 0;
  countRemoved = 0;
  puts (gfr_writeHeader ());
  j = 0;
  for (i = 0; i < arrayMax (gfrEntries); i++) {
    currGE = arrp (gfrEntries,i,GfrEntry);
    homologousCount = 0;
    while (j < arrayMax (seqNames)) {
      if (strStartsWith (textItem (seqNames,j),currGE->id)) {

        testBQ.sequenceName =  hlr_strdup (textItem (seqNames,j));
        if (arrayFind (bowtieQueries,&testBQ,&index,(ARRAYORDERF)sortBowtieQueriesBySequenceName)) {
          isHomologous = 0;
          transcriptNumber = getTranscriptNumber (textItem (seqNames,j));
          currBQ = arrp (bowtieQueries,index,BowtieQuery);
          k = 0; 
          while (k < arrayMax (currBQ->entries)) {
            currBE = arrp (currBQ->entries,k,BowtieEntry);
            if (transcriptNumber == '1') {
              if (strEqual (currBE->chromosome,currGE->nameTranscript2)) {
                isHomologous = 1;
                break;
              }
            }
            if (transcriptNumber == '2') {
              if (strEqual (currBE->chromosome,currGE->nameTranscript1)) {
                isHomologous = 1;
                break;
              }
            }
            k++;
          }
	 
          if (isHomologous == 1) {
            homologousCount++;
          }
        }
        hlr_free (testBQ.sequenceName);  
      }
      else {
        break;
      }
      j++;
    }
    //if (((double)homologousCount / currGE->numInter) <= MAX_FRACTION_HOMOLOGOUS) { 
    if (((double)homologousCount / arrayMax(currGE->readsTranscript1)) <= MAX_FRACTION_HOMOLOGOUS) { 
       // no homology for bowtie, then a blast analysis
      // creating two fasta files with the two genes
      //warn("%d %s", i, currGE->id);
      Stringa fa1 = stringCreate( 100 ); 
      Stringa fa2 = stringCreate( 100 );
      stringPrintf( fa1, "%s_transcript1.fa", currGE->id);
      if (!(fp1 = fopen ( string(fa1) ,"w"))) {
	die ("Unable to open file: %s",string (fa1));
      } 
      arrayFind(transcriptome, &currGE->nameTranscript1, &index,(ARRAYORDERF) sortFastaSequencesByName);
      fprintf( fp1, ">%s\n%s\n", arrp( transcriptome, index, Seq )->name,arrp( transcriptome, index, Seq )->sequence );
      fclose( fp1 );

      stringPrintf( fa2, "%s_transcript2.fa", currGE->id);
      if (!(fp2 = fopen ( string(fa2) ,"w"))) {
	die ("Unable to open file: %s",string (fa2));
      }
      arrayFind(transcriptome, &currGE->nameTranscript2, &index,(ARRAYORDERF) sortFastaSequencesByName);
      fprintf( fp2, ">%s\n%s\n", arrp( transcriptome, index, Seq )->name,arrp( transcriptome, index, Seq )->sequence );
      fclose( fp2 );
      
      // creating the two fasta files with the reads
      stringPrintf( fa1, "%s_reads1.fa", currGE->id);
      if (!(freads1 = fopen ( string(fa1) ,"w"))) {
	die ("Unable to open file: %s",string (fa1));
      }   
      // writing the reads of the first end into file
      for (l = 0; l < arrayMax (currGE->readsTranscript1); l++) {
	char* currRead1 = hlr_strdup( textItem (currGE->readsTranscript1,l)); // read1
	readSize1 = strlen( currRead1 );
	fprintf( freads1, ">%d\n%s\n", l, currRead1 );
	hlr_free( currRead1 );
      }
      fclose( freads1 );
      
      stringPrintf( fa2, "%s_reads2.fa", currGE->id);
      if (!(freads2 = fopen ( string(fa2) ,"w"))) {
	die ("Unable to open file: %s",string (fa2));
      } 
      // writing the reads of the second end into file
      for (l = 0; l < arrayMax (currGE->readsTranscript2); l++) {
	char* currRead2 = hlr_strdup( textItem (currGE->readsTranscript2,l)); // read2
	readSize2 = strlen( currRead2 );
	fprintf( freads2, ">%d\n%s\n", l, currRead2 );
	hlr_free( currRead2 );
      }
      fclose( freads2 );      
      
      //blat of reads2 against the first transcript
      stringPrintf( cmd, "blat -t=dna -out=psl -fine -tileSize=15 %s_transcript1.fa %s_reads2.fa stdout", currGE->id, currGE->id );
      
      // reading the results of blast from Pipe
      blatParser_initFromPipe( string(cmd) );
      while( blQ = blatParser_nextQuery() ) {
	int nucleotideOverlap = getNucleotideOverlap ( blQ );
	if ( nucleotideOverlap > ( ((double)readSize2)*MAX_OVERLAP_ALLOWED) )
	  homologousCount++;
      }
      blatParser_deInit();

      //blat of reads1 against the second transcript
      stringPrintf( cmd, "blat -t=dna -out=psl -fine -tileSize=15 %s_transcript2.fa %s_reads1.fa stdout", currGE->id, currGE->id  );
      blatParser_initFromPipe( string(cmd) );
      while( blQ = blatParser_nextQuery() ) {		
	int nucleotideOverlap = getNucleotideOverlap ( blQ );
	if ( nucleotideOverlap > ( ((double)readSize1)*MAX_OVERLAP_ALLOWED) ) {
	  homologousCount++;
	}
      }
      blatParser_deInit();

    //if (((double)homologousCount / (double)currGE->numInter) <= MAX_FRACTION_HOMOLOGOUS) {       
      if (((double)homologousCount / (double)arrayMax(currGE->readsTranscript1)) <= MAX_FRACTION_HOMOLOGOUS) { 
	// writing the gfrEntry
	puts (gfr_writeGfrEntry (currGE));
	count++;
      } else {
	countRemoved++;
      }
      // removing temporary files
      stringPrintf (cmd,"rm -rf %s_reads?.fa %s_transcript?.fa", currGE->id,currGE->id);
      hlr_system( string(cmd) , 0);      
    }
    else {
      countRemoved++;
    }
  }
  gfr_deInit ();
  stringPrintf (cmd,"rm -rf %s",string (fnSequencesToAlign));
  hlr_system (string (cmd),0);

  stringDestroy (fnSequencesToAlign);
  stringDestroy (cmd);
  stringDestroy (buffer);
  warn ("%s_numRemoved: %d",argv[0],countRemoved);  
  warn ("%s_numGfrEntries: %d",argv[0],count);

  confp_close(conf);

  return EXIT_SUCCESS;
}

