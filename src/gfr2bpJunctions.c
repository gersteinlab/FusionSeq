#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/types.h>

#include <bios/confp.h>
#include <bios/log.h>
#include <bios/format.h>
#include <bios/numUtil.h>
#include <bios/fasta.h>

#include "gfr.h"



#define MAX_NUMBER_OF_JUNCTIONS_PER_FILE 2000000

typedef struct {
  int start;
  int end;
} Region;


static config *Conf = NULL;

static int sortTilesByName (Seq *a, Seq *b)
{
  return strcmp (a->name,b->name);
}

static void updateCounts (Array positions, int start, int end) 
{
  int i;
  for (i = start; i <= end; i++) {
    array (positions,i,int)++;
  }
}

void createSignal( Array regions, Array positions ) {
  Region *currRegion;  
  int j = 0; 
  arrayClear (positions);
  while (j < arrayMax (regions)) {
    currRegion = arrp (regions,j,Region);    
    updateCounts (positions,currRegion->start,currRegion->end); 
    j++;
  }

}

static void getInterEndPoints (Array interReads, int *min1, int *max1, int *min2, int *max2)
{
  GfrInterRead *currGIR;
  int i;
  Array regions1, regions2;
  Region *currRegion1, *currRegion2;
  
  regions1 = arrayCreate( arrayMax( interReads), Region);
  regions2 = arrayCreate( arrayMax( interReads), Region);
  
  for (i = 0; i < arrayMax (interReads); i++) {
    currGIR = arrp (interReads,i,GfrInterRead);
    currRegion1 = arrayp( regions1, arrayMax( regions1 ), Region);
    currRegion2 = arrayp( regions2, arrayMax( regions2 ), Region);    
    currRegion1->start = currGIR->readStart1;
    currRegion1->end = currGIR->readEnd1;
    currRegion2->start = currGIR->readStart2;
    currRegion2->end = currGIR->readEnd2;
  }
  Array positions = arrayCreate (10000000,int);
  createSignal( regions1, positions );
  *min1=-1;*max1=-1;
  for( i=0; i<arrayMax( positions ); i++) {
    if( arru( positions, i, int) > 1 & (*min1)==-1 ) {
      *min1 = i+1;
      break;
    } 
  }
  for( i=arrayMax(positions)-1; i>0; i--) {
    if( arru( positions, i, int) > 1 & (*max1)==-1 ) {
      *max1 = i+1;
      break;
    } 
  }  
  createSignal( regions2, positions );
  *min2=-1;*max2=-1;
  for( i=0; i<arrayMax( positions ); i++) {
    if( arru( positions, i, int) > 1 & (*min2)==-1 ) {
      *min2 = i+1;
      break;
    } 
  }
  for( i=arrayMax(positions)-1; i>0; i--) {
    if( arru( positions, i, int) > 1 & (*max2)==-1 ) {
      *max2 = i+1;
      break;
    } 
  }
  arrayDestroy( regions1 );
  arrayDestroy( regions2 );
  arrayDestroy( positions);
}



static Array getTilesForRegion (int tileSize, 
				char *chromosome, 
				int start, 
				int end, 
				int junctionIsToTheRightOfTile) 
{
  Stringa buffer;
  Stringa targetsFile;
  FILE *fp;
  Array targetSeqs;
  int i;

  buffer = stringCreate (100);
  targetsFile = stringCreate (100);
  stringPrintf (targetsFile,"targets_%d.txt",getpid ());
  if (!(fp = fopen (string (targetsFile),"w")) ){
    die ("Unable to open target file: %s",string (targetsFile));
  }
  if (junctionIsToTheRightOfTile == 1) {
    for (i = start - tileSize; i <= end - tileSize; i++) {
      fprintf (fp,"%s:%d-%d\n",chromosome,i,i + tileSize);
    }
  }
  else {
    for (i = start; i <= end; i++) {
      fprintf (fp,"%s:%d-%d\n",chromosome,i,i + tileSize);
    }
  }
  fclose (fp);
  stringPrintf (buffer,"%s %s/%s stdout -noMask -seqList=%s",
		confp_get(Conf, "BLAT_TWO_BIT_TO_FA"),
		confp_get(Conf, "BLAT_DATA_DIR"),
		confp_get(Conf, "BLAT_TWO_BIT_DATA_FILENAME"),
		string (targetsFile));
  fasta_initFromPipe (string (buffer));
  targetSeqs = fasta_readAllSequences (0);
  fasta_deInit ();
  stringPrintf (buffer,"rm -rf %s",string (targetsFile));
  hlr_system (string (buffer),0);
  stringDestroy (targetsFile);
  stringDestroy (buffer);
  return targetSeqs;
}



static void collectTiles (Array totalTiles, Array regionalTiles)
{
  int i;
  Seq *currSeq1,*currSeq2;

  for (i = 0; i < arrayMax (regionalTiles); i++) {
    currSeq1 = arrp (regionalTiles,i,Seq);
    currSeq2 = arrayp (totalTiles,arrayMax (totalTiles),Seq);
    currSeq2->name = hlr_strdup (currSeq1->name);
    currSeq2->sequence = hlr_strdup (currSeq1->sequence);
    currSeq2->size = currSeq1->size;
    hlr_free (currSeq1->name);
    hlr_free (currSeq1->sequence);
  }
  arrayDestroy (regionalTiles);
}



static Array processRegion (Array exonCoordinates, 
			    char *chromosome, 
			    int regionStart, 
			    int regionEnd, 
			    int junctionIsToTheRightOfTile, 
			    int tileSize, 
			    int sizeFlankingRegion)
{
  GfrExonCoordinate *currEC;
  int i;
  Array totalTiles;
  int overlap;
  int exonOverlap=0;
  totalTiles = arrayCreate (10000,Seq);
  for (i = 0; i < arrayMax (exonCoordinates); i++) {
    currEC = arrp (exonCoordinates,i,GfrExonCoordinate);
    overlap = rangeIntersection (currEC->start,currEC->end,regionStart,regionEnd);
    if (overlap > 0) { // overlapping exons
      exonOverlap++;
      if (regionStart < currEC->start && currEC->end < regionEnd) {
        // contained exons, has a flanking region on each end
        collectTiles (totalTiles,getTilesForRegion (tileSize,chromosome,currEC->start - sizeFlankingRegion,currEC->end + sizeFlankingRegion,junctionIsToTheRightOfTile));
      }
      else if (regionStart >= currEC->start && currEC->end <= regionEnd) {
        // left endpoint, has one flanking region on the right
        collectTiles (totalTiles,getTilesForRegion (tileSize,chromosome,regionStart,currEC->end + sizeFlankingRegion,junctionIsToTheRightOfTile));
      }
      else if (regionStart <= currEC->start && currEC->end >= regionEnd) {
        // right endpoint, has one flanking region on the left
        collectTiles (totalTiles,getTilesForRegion (tileSize,chromosome,currEC->start - sizeFlankingRegion,regionEnd,junctionIsToTheRightOfTile));
      }
    }
  }
  if( exonOverlap == 0 ) { // no exon overlap, only intronic region
    warn( "Only intronic");
  }
  return totalTiles;
}



static void freeTiles (Array tiles)
{
  Seq *currSeq;
  int i;

  for (i = 0; i < arrayMax (tiles); i++) {
    currSeq = arrp (tiles,i,Seq);
    hlr_free (currSeq->name);
    hlr_free (currSeq->sequence);
  }
  arrayDestroy (tiles);
}



static void printCommands (char *id, 
			   char *orientation, 
			   int fileCount, 
			   FILE *fpJobList1, 
			   FILE *fpJobList2, 
			   char *gfrPrefix, 
			   int isFirst, 
			   int isLast)
{
  long size;
  char *buf;
  static char *currentDirectory = NULL;
  static int first = 1;

  if (first == 1) {
    size = pathconf(".", _PC_PATH_MAX);
    if ((buf = (char *)malloc((size_t)size)) != NULL) {
      currentDirectory = getcwd(buf, (size_t)size);
    }
    first = 0;
  }
  fprintf (fpJobList1,"cd %s; bowtie-build -q -f %s_%s_%d.fa %s/%s_%s_%d; zcat %s_allReads.txt.gz | bowtie --quiet -p 4 -r %s/%s_%s_%d - %s_%s_%d.bowtie; rm -rf %s_%s_%d.fa %s/%s_%s_%d.*.ebwt\n",
           currentDirectory,
	   id,
	   orientation,
	   fileCount,
	   confp_get(Conf, "BOWTIE_INDEXES"),
	   id,
	   orientation,
	   fileCount,
	   gfrPrefix,
	   confp_get(Conf, "BOWTIE_INDEXES"),
	   id,
	   orientation,
	   fileCount,
	   id,
	   orientation,
	   fileCount,
	   id,
	   orientation,
	   fileCount,
	   confp_get(Conf, "BOWTIE_INDEXES"),
	   id,
	   orientation,
	   fileCount);
  if (isFirst == 1) {
    fprintf (fpJobList2,"cd %s; cat ",currentDirectory);
  }
  fprintf (fpJobList2,"%s_%s_%d.bowtie ",id,orientation,fileCount);
  if (isLast == 1) {
    fprintf (fpJobList2,"| sort -n | bowtie2bp > %s_%s.bp\n",id,orientation);
  }
}



static void writeTiles (Array tiles1, 
			Array tiles2, 
			char *id, 
			char *orientation, 
			FILE *fpJobList1, 
			FILE *fpJobList2, 
			char *gfrPrefix)
{
  static Stringa buffer = NULL;
  int fileCount;
  int numJunctions;
  FILE *fp;
  int i,j;
  Seq *firstSeq,*secondSeq;

  arraySort (tiles1,(ARRAYORDERF)sortTilesByName);
  arrayUniq (tiles1,NULL,(ARRAYORDERF)sortTilesByName);
  arraySort (tiles2,(ARRAYORDERF)sortTilesByName);
  arrayUniq (tiles2,NULL,(ARRAYORDERF)sortTilesByName);
  warn ("%s, %s, num_tiles_transcript_1: %d, num_tiles_transcript_2: %d",
	id,
	orientation,
	arrayMax (tiles1),
	arrayMax (tiles2));
  fileCount = 1;
  stringCreateClear (buffer,100);
  stringPrintf (buffer,"%s_%s_%d.fa",id,orientation,fileCount);
  fp = fopen (string (buffer),"w");
  if (fp == NULL) {
    die ("Unable to open file: %s",string (buffer));
  }
  numJunctions = 0;
  for (i = 0; i < arrayMax (tiles1); i++) {
    firstSeq = arrp (tiles1,i,Seq);
    for (j = 0; j < arrayMax (tiles2); j++) {
      secondSeq = arrp (tiles2,j,Seq);
      numJunctions++; 
       if (numJunctions == MAX_NUMBER_OF_JUNCTIONS_PER_FILE) {
        printCommands (id,orientation,fileCount,fpJobList1,fpJobList2,gfrPrefix,fileCount == 1 ? 1 : 0,0);
        fclose (fp);
        numJunctions = 0;
        fileCount++;
        stringPrintf (buffer,"%s_%s_%d.fa",id,orientation,fileCount);
        fp = fopen (string (buffer),"w");
        if (fp == NULL) {
          die ("Unable to open file: %s",string (buffer));
        }
      }
      fprintf (fp,">%s|%s\n",firstSeq->name,secondSeq->name);
      fprintf (fp,"%s%s\n",firstSeq->sequence,secondSeq->sequence);
    }
  }
  printCommands (id,orientation,fileCount,fpJobList1,fpJobList2,gfrPrefix,fileCount == 1 ? 1 : 0,1);
  fclose (fp);
  freeTiles (tiles1);
  freeTiles (tiles2);
}



int main (int argc, char *argv[])
{
  GfrEntry *currGE;
  int min1,max1,min2,max2;
  int tileSize;
  int sizeFlankingRegion;
  Array tiles1;
  Array tiles2;
  char *gfrPrefix;
  char *pos;
  FILE *fpJobList1,*fpJobList2;
  Stringa buffer;
  double minDASPER;

  if (argc <= 4) {
    usage ("%s <file.gfr> <tileSize> <sizeFlankingRegion> [minDASPER]",argv[0]);
  }

  if ((Conf = confp_open(getenv("FUSIONSEQ_CONFPATH"))) == NULL)
    return EXIT_FAILURE;

  gfr_init (argv[1]);
  tileSize = atoi (argv[2]);
  sizeFlankingRegion = atoi (argv[3]);
  if(argc==5) {
    minDASPER = atof( argv[4] );
  } else {
    minDASPER = -100000000;
  }
  gfrPrefix = hlr_strdup (argv[1]);
  pos = strchr (gfrPrefix,'.');
  if (pos == NULL) {
    die ("Expected a '.gfr' extension: %s",argv[1]);
  }
  *pos = '\0';
  buffer = stringCreate (100);
  stringPrintf (buffer,"%s_jobList1.txt",gfrPrefix);
  fpJobList1 = fopen (string (buffer),"w");
  stringPrintf (buffer,"%s_jobList2.txt",gfrPrefix);
  fpJobList2 = fopen (string (buffer),"w");
  if (fpJobList1 == NULL || fpJobList2 == NULL) {
    die ("Unable to open job list files!");
  }
  while (currGE = gfr_nextEntry ()){
    if( currGE->DASPER < minDASPER ) continue;
    getInterEndPoints (currGE->interReads,&min1,&max1,&min2,&max2);
    tiles1 = processRegion (currGE->exonCoordinatesTranscript1,
		    	    currGE->chromosomeTranscript1,
			    min1,
			    currGE->endTranscript1,
			    1,
			    tileSize,
			    sizeFlankingRegion);
    tiles2 = processRegion (currGE->exonCoordinatesTranscript2,
		    	    currGE->chromosomeTranscript2,
			    currGE->startTranscript2,
			    max2,
			    0,
			    tileSize,
			    sizeFlankingRegion);
    writeTiles (tiles1,tiles2,currGE->id,"AB",fpJobList1,fpJobList2,gfrPrefix);
    tiles1 = processRegion (currGE->exonCoordinatesTranscript2,
		    	    currGE->chromosomeTranscript2,
			    min2,
			    currGE->endTranscript2,
			    1,
			    tileSize,
			    sizeFlankingRegion);
    tiles2 = processRegion (currGE->exonCoordinatesTranscript1,
		    	    currGE->chromosomeTranscript1,
			    currGE->startTranscript1,
			    max1,
			    0,
			    tileSize,
			    sizeFlankingRegion);    
    writeTiles (tiles1,tiles2,currGE->id,"BA",fpJobList1,fpJobList2,gfrPrefix);//*/
  }
  gfr_deInit ();
  fclose (fpJobList1);
  fclose (fpJobList2);
  hlr_free (gfrPrefix);
  stringDestroy (buffer);
  confp_close(Conf);

  return EXIT_SUCCESS;
}
