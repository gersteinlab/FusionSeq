#include <stdlib.h>
#include <time.h>
#include <math.h>

#include <bios/confp.h>
#include <bios/log.h>
#include <bios/format.h>
#include <bios/intervalFind.h>
#include <mrf/mrf.h>

#include "gfr.h"


#define SAMPLING_ITERATIONS 100000

static config *Conf = NULL;

typedef struct {
  Interval *transcript;
  char* read1;
  char* read2;
  int readStart1;
  int readStart2;
  int readEnd1;
  int readEnd2;
} Intra;



typedef struct {
  Interval *transcript1;
  Interval *transcript2;
  char* read1;
  char* read2;
  int readStart1;
  int readStart2;
  int readEnd1;
  int readEnd2;
  int pairType;
  int number1;
  int number2;
} Inter;



typedef struct {
  Interval *transcript;
  Array intras;  // of type Intra*
} SuperIntra;



typedef struct {
  Interval *transcript1;
  Interval *transcript2;
  Array inters;  // of type Inter*
} SuperInter;



typedef struct {
  char* chromosome;
  int genomic;
  int transcript;
} Coordinate; 


int getNumberOfIntras( SuperIntra* a ) {
  int i;
  double numberOfIntras=0.0;
  Intra* currIntra;
  for( i=0; i<arrayMax( a->intras); i++) {
    currIntra = arru( a->intras, i, Intra* );
    if( currIntra->read1==NULL | currIntra->read2==NULL | 
	(currIntra->readEnd1 - currIntra->readStart1 + 1) == 0 | 
	(currIntra->readEnd2 - currIntra->readStart2 + 1) == 0 ) 
      die("Something is wrong with the intra pairs: read1[%s](e%d-s%d +1 )=%d read2[%s](e%d-s%d + 1 )=%d", 
	  currIntra->read1, currIntra->readEnd1, currIntra->readStart1, (currIntra->readEnd1 - currIntra->readStart1 + 1),
	  currIntra->read2, currIntra->readEnd2, currIntra->readStart2, (currIntra->readEnd2 - currIntra->readStart2 + 1) );
    if( (currIntra->readEnd1 - currIntra->readStart1 + 1) != strlen(currIntra->read1) &
	(currIntra->readEnd2 - currIntra->readStart2 + 1) != strlen(currIntra->read2) ) {
      numberOfIntras += 0.25;
    } else if ( (currIntra->readEnd1 - currIntra->readStart1 + 1) != strlen(currIntra->read1) |
		(currIntra->readEnd2 - currIntra->readStart2 + 1) != strlen(currIntra->read2) ) {
      numberOfIntras += 0.5;
    } else {
      numberOfIntras += 1.0;
    }
  }
  return (int) rint( numberOfIntras );
}

int getNumberOfInters( SuperInter* a ) {
  int i;
  double numberOfInters=0.0;
  Inter* currInter;
  for( i=0; i<arrayMax( a->inters); i++) {
    currInter = arru( a->inters, i, Inter* );
    if( currInter->read1==NULL | currInter->read2==NULL | 
	(currInter->readEnd1 - currInter->readStart1 + 1) == 0 | 
	(currInter->readEnd2 - currInter->readStart2 + 1) == 0 ) 
      die("Something is wrong with the inter pairs: read1[%s](e%d-s%d +1 )=%d read2[%s](e%d-s%d + 1 )=%d", 
	  currInter->read1, currInter->readEnd1, currInter->readStart1, (currInter->readEnd1 - currInter->readStart1 + 1),
	  currInter->read2, currInter->readEnd2, currInter->readStart2, (currInter->readEnd2 - currInter->readStart2 + 1) );
    if( (currInter->readEnd1 - currInter->readStart1 + 1) != strlen(currInter->read1) &
	(currInter->readEnd2 - currInter->readStart2 + 1) != strlen(currInter->read2) ) {
      numberOfInters += 0.25;
    } else if ( (currInter->readEnd1 - currInter->readStart1 + 1) != strlen(currInter->read1) |
		(currInter->readEnd2 - currInter->readStart2 + 1) != strlen(currInter->read2) ) {
      numberOfInters += 0.5;
    } else {
      numberOfInters += 1.0;
    }
  }
  return (int) rint( numberOfInters );
}


static int sortIntersByTranscript (Inter *a, Inter *b)
{
  int diff;
  
  diff = a->transcript1 - b->transcript1;
  if (diff != 0) {
    return diff;
  }
  return a->transcript2 - b->transcript2;
}



static int sortIntrasByTranscript (Intra *a, Intra *b)
{
  return a->transcript - b->transcript;
}



static int sortSuperInters (SuperInter *a, SuperInter *b)
{
  //return arrayMax (b->inters) - arrayMax (a->inters);
  return getNumberOfInters(b) - getNumberOfInters(a);
}



static int sortSuperIntras (SuperIntra *a, SuperIntra *b)
{
  return a->transcript - b->transcript;
}



static int getExonNumber (Interval* transcript, int start, int end)
{
  SubInterval *currExon;
  int i;

  for (i = 0; i < arrayMax (transcript->subIntervals); i++) {
    currExon = arrp (transcript->subIntervals,i,SubInterval);
    if (start >= currExon->start && end <= currExon->end) {
      return i + 1;
    }
  }
  return 0;
}



static int getIntronNumber (Interval* transcript, int start, int end)
{
  SubInterval *currExon,*prevExon;
  int i;

  for (i = 1; i < arrayMax (transcript->subIntervals); i++) {
    prevExon = arrp (transcript->subIntervals,i - 1,SubInterval);
    currExon = arrp (transcript->subIntervals,i,SubInterval);
    if (start > prevExon->end && end < currExon->start) {
      return i;
    }
  }
  return 0;
}



static int getJunctionNumber (Interval* transcript, int start, int end, int exonNumber, int intronNumber)
{
  SubInterval *currExon;
  int i;
  
  if (exonNumber > 0 || intronNumber > 0) {
    return 0;
  }
  for (i = 0; i < arrayMax (transcript->subIntervals); i++) {
    currExon = arrp (transcript->subIntervals,i,SubInterval);
    if (start <= currExon->start && end >= currExon->start) {
      return i * 2 + 1;
    }
    if (start <= currExon->end && end >= currExon->end) {
      return i * 2 + 2;
    }
  }
  return 0;
}



static void assignPairType (Inter *currInter, int exon1, int intron1, int junction1, int exon2, int intron2, int junction2)
{
  if (exon1 > 0 && exon2 > 0) {
    currInter->pairType = GFR_PAIR_TYPE_EXONIC_EXONIC;
    currInter->number1 = exon1;
    currInter->number2 = exon2;
  }
  else if (exon1 > 0 && intron2 > 0) {
    currInter->pairType = GFR_PAIR_TYPE_EXONIC_INTRONIC;
    currInter->number1 = exon1;
    currInter->number2 = intron2;
  }
  else if (exon1 > 0 && junction2 > 0) {
    currInter->pairType = GFR_PAIR_TYPE_EXONIC_JUNCTION;
    currInter->number1 = exon1;
    currInter->number2 = junction2;
  }
  else if (intron1 > 0 && exon2 > 0) {
    currInter->pairType = GFR_PAIR_TYPE_INTRONIC_EXONIC;
    currInter->number1 = intron1;
    currInter->number2 = exon2;
  }
  else if (intron1 > 0 && intron2 > 0) {
    currInter->pairType = GFR_PAIR_TYPE_INTRONIC_INTRONIC;
    currInter->number1 = intron1;
    currInter->number2 = intron2;
  }
  else if (intron1 > 0 && junction2 > 0) {
    currInter->pairType = GFR_PAIR_TYPE_INTRONIC_JUNCTION;
    currInter->number1 = intron1;
    currInter->number2 = junction2;
  }
  else if (junction1 > 0 && junction2 > 0) {
    currInter->pairType = GFR_PAIR_TYPE_JUNCTION_JUNCTION;
    currInter->number1 = junction1;
    currInter->number2 = junction2;
  }
  else if (junction1 > 0 && exon2 > 0) {
    currInter->pairType = GFR_PAIR_TYPE_JUNCTION_EXONIC;
    currInter->number1 = junction1;
    currInter->number2 = exon2;
  }
  else if (junction1 > 0 && intron2 > 0) {
    currInter->pairType = GFR_PAIR_TYPE_JUNCTION_INTRONIC;
    currInter->number1 = junction1;
    currInter->number2 = intron2;
  }
  else {
    die ("Unexpected pair type: %d %d %d %d %d %d",exon1,intron1,junction1,exon2,intron2,junction2);
  }
}



static int sortCoordinates (Coordinate *a, Coordinate *b)
{
  int diff;
  
  diff = strcmp (a->chromosome,b->chromosome);
  if (diff != 0) {
    return diff;
  }
  return a->genomic - b->genomic;
}



static Array convertIntraCoordinates (Interval *transcript)
{
  int i,j,k;
  SubInterval *currExon;
  Coordinate *currCoordinate;
  Array coordinates;

  coordinates = arrayCreate (10000,Coordinate);
  k = 1;
  for (i = 0; i < arrayMax (transcript->subIntervals); i++) {
    currExon = arrp (transcript->subIntervals,i,SubInterval);
    for (j = currExon->start; j <= currExon->end; j++) {
      currCoordinate = arrayp (coordinates,arrayMax (coordinates),Coordinate);
      currCoordinate->genomic = j;
      currCoordinate->transcript = k;
      currCoordinate->chromosome = hlr_strdup (transcript->chromosome);
      k++;
    }
  }
  return coordinates;
  // Note: this array is already sorted by genomic coodinate
}



static void addInterCoordinates (Interval *transcript, Array coordinates, int *transcriptIndex, 
                                 int startFusionTranscript, int endFusionTranscript)
{
  SubInterval *currExon;
  Coordinate *currCoordinate;
  int i,j;

  for (i = 0; i < arrayMax (transcript->subIntervals); i++) {
    currExon = arrp (transcript->subIntervals,i,SubInterval);
    for (j = currExon->start; j <= currExon->end; j++) {
      if (j >= startFusionTranscript && j <= endFusionTranscript) {
        currCoordinate = arrayp (coordinates,arrayMax (coordinates),Coordinate);
        currCoordinate->genomic = j;
        currCoordinate->transcript = *transcriptIndex;
        currCoordinate->chromosome = hlr_strdup (transcript->chromosome);
        (*transcriptIndex)++;
      }
    }
  }
}



static Array convertInterCoordinates (SuperInter *currSuperInter, int isAB)
{
  int i,start;
  int startFusionTranscript1,endFusionTranscript1;
  int startFusionTranscript2,endFusionTranscript2;
  Inter *currInter;
  Array coordinates;
  int transcriptIndex;

  start = 0;
  while (start < arrayMax (currSuperInter->inters)) {
    currInter = arru (currSuperInter->inters,start,Inter*);
    if (currInter->pairType == GFR_PAIR_TYPE_EXONIC_EXONIC) {
      break;
    }
    start++;
  }
  startFusionTranscript1 = currInter->readStart1;
  endFusionTranscript1 = currInter->readEnd1;
  startFusionTranscript2 = currInter->readStart2;
  endFusionTranscript2 = currInter->readEnd2;
  for (i = start + 1; i < arrayMax (currSuperInter->inters); i++) {
    currInter = arru (currSuperInter->inters,i,Inter*);
    if (currInter->pairType != GFR_PAIR_TYPE_EXONIC_EXONIC) {
      continue;
    }
    if (currInter->readStart1 < startFusionTranscript1) {
      startFusionTranscript1 = currInter->readStart1;
    } 
    if (currInter->readEnd1 > endFusionTranscript1) {
      endFusionTranscript1 = currInter->readEnd1;
    } 
    if (currInter->readStart2 < startFusionTranscript2) {
      startFusionTranscript2 = currInter->readStart2;
    } 
    if (currInter->readEnd2 > endFusionTranscript2) {
      endFusionTranscript2 = currInter->readEnd2;
    } 
  }
  coordinates = arrayCreate (1000,Coordinate);
  transcriptIndex = 1;
  if (isAB == 1) {
    addInterCoordinates (currSuperInter->transcript1,coordinates,&transcriptIndex,startFusionTranscript1,endFusionTranscript1);
    addInterCoordinates (currSuperInter->transcript2,coordinates,&transcriptIndex,startFusionTranscript2,endFusionTranscript2);
  }
  else if (isAB == 0) {
    addInterCoordinates (currSuperInter->transcript2,coordinates,&transcriptIndex,startFusionTranscript2,endFusionTranscript2);
    addInterCoordinates (currSuperInter->transcript1,coordinates,&transcriptIndex,startFusionTranscript1,endFusionTranscript1);
  }
  else {
    die ("Unknown mode: %d",isAB);
  }
  arraySort (coordinates,(ARRAYORDERF)sortCoordinates);
  return coordinates;
}



static int getTranscriptCoordinate (Array coordinates, int genomicCoordinate, char* chromosome) 
{
  Coordinate testCoordinate;
  int index;
 
  testCoordinate.genomic = genomicCoordinate;
  testCoordinate.chromosome = hlr_strdup (chromosome);
  if (!arrayFind (coordinates,&testCoordinate,&index,(ARRAYORDERF)sortCoordinates)) {
    die ("Expected to find coordinate: %s %d",chromosome,genomicCoordinate); 
  }
  hlr_free (testCoordinate.chromosome);
  return arrp (coordinates,index,Coordinate)->transcript;
}



static void calculateIntraOffsets (Array coordinatesTanscript, SuperIntra *currSuperIntra, Array intraOffsets)
{
  int i;
  Intra *currIntra;
  int index1,index2;

  for (i = 0; i < arrayMax (currSuperIntra->intras); i++) {
    currIntra = arru (currSuperIntra->intras,i,Intra*);
    index2 = getTranscriptCoordinate (coordinatesTanscript,currIntra->readEnd2,currIntra->transcript->chromosome);
    index1 = getTranscriptCoordinate (coordinatesTanscript,currIntra->readStart1,currIntra->transcript->chromosome);
    array (intraOffsets,arrayMax (intraOffsets),int) = index2 - index1 + 1;
  }
}



static void calculateInterOffsets (Array interCoordinates, SuperInter *currSuperInter, Array interOffsets, int isAB)
{
  int i;
  Inter *currInter;
  int index1,index2;

  for (i = 0; i < arrayMax (currSuperInter->inters); i++) {
    currInter = arru (currSuperInter->inters,i,Inter*);
    if (currInter->pairType != GFR_PAIR_TYPE_EXONIC_EXONIC) {
      continue;
    }
    if (isAB == 1) {
      index2 = getTranscriptCoordinate (interCoordinates,currInter->readEnd2,currInter->transcript2->chromosome);
      index1 = getTranscriptCoordinate (interCoordinates,currInter->readStart1,currInter->transcript1->chromosome);
      array (interOffsets,arrayMax (interOffsets),int) = index2 - index1 + 1;
    }
    else if (isAB == 0) {
      index1 = getTranscriptCoordinate (interCoordinates,currInter->readEnd1,currInter->transcript1->chromosome);
      index2 = getTranscriptCoordinate (interCoordinates,currInter->readStart2,currInter->transcript2->chromosome);
      array (interOffsets,arrayMax (interOffsets),int) = index1 - index2 + 1;
    }
    else {
      die ("Unknown mode: %d",isAB);
    }
  }
}



static void destroyCoordinates (Array coordinates)
{
  Coordinate *currCoordinate;
  int i;

  for (i = 0; i < arrayMax (coordinates); i++) {
    currCoordinate = arrp (coordinates,i,Coordinate);
    hlr_free (currCoordinate->chromosome);
  }
  arrayDestroy (coordinates);
}



static double calculateMean (Array integers) 
{
  int i;
  int count;

  count = 0;
  for (i = 0; i < arrayMax (integers); i++) {
    count += arru (integers,i,int);
  }
  return 1.0 * count / arrayMax (integers);
}
int sortIntegers(int *a, int *b) {
  return *b - *a;
}
static double calculateMedian (Array integers) 
{
  arraySort( integers, (ARRAYORDERF) sortIntegers);
  return arru (integers, arrayMax( integers )/2,int);
}


static double compareDistributions (Array intraOffsets, Array interOffsets)
{
  double meanIntra, medianInter;
  int i,j;
  static Array sampledIntraOffsets = NULL;
  int count;

  if (sampledIntraOffsets == NULL) {
    sampledIntraOffsets = arrayCreate (1000,int);
  }
  count = 0;
  //meanInter = calculateMean (interOffsets);
  medianInter = calculateMedian( interOffsets );
  for (i = 0; i < SAMPLING_ITERATIONS; i++) {
    arrayClear (sampledIntraOffsets);
    for (j = 0; j < arrayMax (interOffsets); j++) {
      array (sampledIntraOffsets,arrayMax (sampledIntraOffsets),int) = arru (intraOffsets,rand () % arrayMax (intraOffsets),int);
    }
    meanIntra = calculateMean (sampledIntraOffsets);
    if (medianInter > meanIntra) { //replaced mean with median
      count++;
    }
  }
  return 1 - (1.0 * count / SAMPLING_ITERATIONS);
}



static SuperIntra* getSuperIntra (Array superIntras, Interval *testInterval)
{
   SuperIntra testSuperIntra;
   int index;

   testSuperIntra.transcript = testInterval;
   if (arrayFind (superIntras,&testSuperIntra,&index,(ARRAYORDERF)sortSuperIntras)) {
     return arrp (superIntras,index,SuperIntra);
   } 
   return NULL;
}



static char* exonCoordinates2string (Array exons)
{
  int i;
  static Stringa buffer = NULL;
  SubInterval *currExon;

  stringCreateClear (buffer,100);
  for (i = 0; i < arrayMax (exons); i++) {
    currExon = arrp (exons,i,SubInterval);
    stringAppendf (buffer,"%d,%d%s",currExon->start,currExon->end,i < arrayMax (exons) - 1 ? "|" : "");
  }
  return string (buffer);
}



static int isValidForInsertSizeFilter (SuperInter *currSuperInter)
{
  int i;
  Inter *currInter;

  for (i = 0; i < arrayMax (currSuperInter->inters); i++) {
    currInter = arru (currSuperInter->inters,i,Inter*);
    if (currInter->pairType == GFR_PAIR_TYPE_EXONIC_EXONIC) {
      return 1;
    } 
  }
  return 0;
}



int main (int argc, char *argv[])
{
  MrfEntry *currMrfEntry;
  MrfRead *currMrfRead1,*currMrfRead2;
  MrfBlock *currMrfBlock1,*currMrfBlock2;
  Stringa buffer;
  Array intervals1;
  Array intervals2;
  Interval *transcript1,*transcript2;
  Array inters;
  Inter *currInter,*nextInter;
  Array intras;
  Intra *currIntra,*nextIntra;
  int i,j;
  int exon1,exon2,intron1,intron2,junction1,junction2;
  SuperInter *currSuperInter;
  Array superInters;
  SuperIntra *currSuperIntra,*superIntra1,*superIntra2;
  Array superIntras;
  Array coordinatesTanscript;
  Array interCoordinatesAB,interCoordinatesBA;
  Array intraOffsets;
  Array interOffsetsAB,interOffsetsBA;
  double pvalueAB,pvalueBA;
  double meanInterAB,meanInterBA;
  FILE *fp;
  int mrfLines;
  char *exonCoordinates1,*exonCoordinates2;
  int minNumberOfPairedEndReads;

  if ((Conf = confp_open(getenv("FUSIONSEQ_CONFPATH"))) == NULL)
    return EXIT_FAILURE;

  if (argc != 3) {
    usage ("%s <prefix> <minNumberOfPairedEndReads>",argv[0]);
  }
  minNumberOfPairedEndReads = atoi (argv[2]);
  srand (time (0));
  buffer = stringCreate (100);
  stringPrintf (buffer,"%s/%s", 
                confp_get(Conf, "TRANSCRIPT_COMPOSITE_MODEL_DIR"), 
                confp_get(Conf, "TRANSCRIPT_COMPOSITE_MODEL_FILENAME"));
  intervalFind_addIntervalsToSearchSpace (string (buffer),0);
  inters = arrayCreate (1000000,Inter);
  intras = arrayCreate (1000000,Intra);
  mrfLines = 0;
  mrf_init ("-");
  while (currMrfEntry = mrf_nextEntry ()) {
    mrfLines++;
    currMrfRead1 = &currMrfEntry->read1;
    currMrfRead2 = &currMrfEntry->read2;
    for (i = 0; i < arrayMax (currMrfRead1->blocks); i++) {
      currMrfBlock1 = arrp (currMrfRead1->blocks,i,MrfBlock);
      for (j = 0; j < arrayMax (currMrfRead2->blocks); j++) {
        currMrfBlock2 = arrp (currMrfRead2->blocks,j,MrfBlock);  
        intervals1 = arrayCopy (intervalFind_getOverlappingIntervals (currMrfBlock1->targetName,currMrfBlock1->targetStart,currMrfBlock1->targetEnd));    
        intervals2 = arrayCopy (intervalFind_getOverlappingIntervals (currMrfBlock2->targetName,currMrfBlock2->targetStart,currMrfBlock2->targetEnd));    
        if (arrayMax (intervals1) == 1 && arrayMax (intervals2) == 1) {
          transcript1 = arru (intervals1,0,Interval*);
          transcript2 = arru (intervals2,0,Interval*);
          exon1 = getExonNumber (transcript1,currMrfBlock1->targetStart,currMrfBlock1->targetEnd);
          exon2 = getExonNumber (transcript2,currMrfBlock2->targetStart,currMrfBlock2->targetEnd);
          intron1 = getIntronNumber (transcript1,currMrfBlock1->targetStart,currMrfBlock1->targetEnd);
          intron2 = getIntronNumber (transcript2,currMrfBlock2->targetStart,currMrfBlock2->targetEnd);
          junction1 = getJunctionNumber (transcript1,currMrfBlock1->targetStart,currMrfBlock1->targetEnd,exon1,intron1);
          junction2 = getJunctionNumber (transcript2,currMrfBlock2->targetStart,currMrfBlock2->targetEnd,exon2,intron2);
          if (transcript1 != transcript2) {
            currInter = arrayp (inters,arrayMax (inters),Inter);
            currInter->transcript1 = transcript1;
            currInter->transcript2 = transcript2;
            currInter->readStart1 = currMrfBlock1->targetStart;
            currInter->readStart2 = currMrfBlock2->targetStart;
            currInter->readEnd1 = currMrfBlock1->targetEnd;
            currInter->readEnd2 = currMrfBlock2->targetEnd;
            currInter->read1 = hlr_strdup (currMrfRead1->sequence);
            currInter->read2 = hlr_strdup (currMrfRead2->sequence);
            assignPairType (currInter,exon1,intron1,junction1,exon2,intron2,junction2);
          }
          if (transcript1 == transcript2) {
            if (exon1 > 0 && exon2 > 0) {
              currIntra = arrayp (intras,arrayMax (intras),Intra);
              currIntra->transcript = transcript1;
              currIntra->readStart1 = currMrfBlock1->targetStart;
              currIntra->readStart2 = currMrfBlock2->targetStart;
              currIntra->readEnd1 = currMrfBlock1->targetEnd;
              currIntra->readEnd2 = currMrfBlock2->targetEnd;
              currIntra->read1 = hlr_strdup (currMrfRead1->sequence);
              currIntra->read2 = hlr_strdup (currMrfRead2->sequence);
            }
          }
        }
        arrayDestroy (intervals1);
        arrayDestroy (intervals2);
      }
    }
  }
  mrf_deInit ();
  arraySort (inters,(ARRAYORDERF)sortIntersByTranscript);
  arraySort (intras,(ARRAYORDERF)sortIntrasByTranscript);

  intraOffsets = arrayCreate (1000000,int);
  superIntras = arrayCreate (100000,SuperIntra);
  i = 0;
  while (i < arrayMax (intras)) {
    currIntra = arrp (intras,i,Intra);
    currSuperIntra = arrayp (superIntras,arrayMax (superIntras),SuperIntra);
    currSuperIntra->transcript = currIntra->transcript;
    currSuperIntra->intras = arrayCreate (100,Intra*);
    array (currSuperIntra->intras,arrayMax (currSuperIntra->intras),Intra*) = currIntra;
    j = i + 1;
    while (j < arrayMax (intras)) {
      nextIntra = arrp (intras,j,Intra);
      if (currIntra->transcript == nextIntra->transcript) {
        array (currSuperIntra->intras,arrayMax (currSuperIntra->intras),Intra*) = nextIntra;
      }
      else {
        break;
      }
      j++;
    }
    i = j;
    coordinatesTanscript = convertIntraCoordinates (currSuperIntra->transcript);
    calculateIntraOffsets (coordinatesTanscript,currSuperIntra,intraOffsets);
    destroyCoordinates (coordinatesTanscript);
  }
  arraySort (superIntras,(ARRAYORDERF)sortSuperIntras);
  
  superInters = arrayCreate (100000,SuperInter);
  i = 0;
  while (i < arrayMax (inters)) {
    currInter = arrp (inters,i,Inter);
    currSuperInter = arrayp (superInters,arrayMax (superInters),SuperInter);
    currSuperInter->transcript1 = currInter->transcript1;
    currSuperInter->transcript2 = currInter->transcript2;
    currSuperInter->inters = arrayCreate (100,Inter*);
    array (currSuperInter->inters,arrayMax (currSuperInter->inters),Inter*) = currInter;
    j = i + 1;
    while (j < arrayMax (inters)) {
      nextInter = arrp (inters,j,Inter);
      if (currInter->transcript1 == nextInter->transcript1 &&
          currInter->transcript2 == nextInter->transcript2) {
        array (currSuperInter->inters,arrayMax (currSuperInter->inters),Inter*) = nextInter;
      }
      else {
        break;
      }
      j++;
    }    
    i = j;
  }
  arraySort (superInters,(ARRAYORDERF)sortSuperInters);

  warn ("%s_numMrfLines: %d",argv[0],mrfLines);
  warn ("%s_numIntra: %d",argv[0],arrayMax (intras));
  warn ("%s_numInter: %d",argv[0],arrayMax (inters));
  warn ("%s_numSuperIntra: %d",argv[0],arrayMax (superIntras));
  warn ("%s_numSuperInter: %d",argv[0],arrayMax (superInters));
  printf ("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",
          GFR_COLUMN_NAME_NUM_INTER,
          GFR_COLUMN_NAME_INTER_MEAN_AB,
          GFR_COLUMN_NAME_INTER_MEAN_BA,
          GFR_COLUMN_NAME_PVALUE_AB,
          GFR_COLUMN_NAME_PVALUE_BA,
          GFR_COLUMN_NAME_NUM_INTRA1,
          GFR_COLUMN_NAME_NUM_INTRA2,
          GFR_COLUMN_NAME_FUSION_TYPE,
          GFR_COLUMN_NAME_NAME_TRANSCRIPT1,
          GFR_COLUMN_NAME_NUM_EXONS_TRANSCRIPT1, 
          GFR_COLUMN_NAME_EXON_COORDINATES_TRANSCRIPT1,
          GFR_COLUMN_NAME_CHROMOSOME_TRANSCRIPT1,
          GFR_COLUMN_NAME_STRAND_TRANSCRIPT1,
          GFR_COLUMN_NAME_START_TRANSCRIPT1,
          GFR_COLUMN_NAME_END_TRANSCRIPT1,
          GFR_COLUMN_NAME_NAME_TRANSCRIPT2,
          GFR_COLUMN_NAME_NUM_EXONS_TRANSCRIPT2,
          GFR_COLUMN_NAME_EXON_COORDINATES_TRANSCRIPT2,
          GFR_COLUMN_NAME_CHROMOSOME_TRANSCRIPT2,
          GFR_COLUMN_NAME_STRAND_TRANSCRIPT2,
          GFR_COLUMN_NAME_START_TRANSCRIPT2,
          GFR_COLUMN_NAME_END_TRANSCRIPT2,
          GFR_COLUMN_NAME_INTER_READS,
          GFR_COLUMN_NAME_ID,
          GFR_COLUMN_NAME_READS_TRANSCRIPT1,
          GFR_COLUMN_NAME_READS_TRANSCRIPT2);

  interOffsetsAB = arrayCreate (1000,int);
  interOffsetsBA = arrayCreate (1000,int);
  i = 0; 
  while (i < arrayMax (superInters)) {
    currSuperInter = arrp (superInters,i,SuperInter);
    if ( getNumberOfInters (currSuperInter) < minNumberOfPairedEndReads) {
      break;
    }
    if (isValidForInsertSizeFilter (currSuperInter)) {
      interCoordinatesAB = convertInterCoordinates (currSuperInter,1);
      arrayClear (interOffsetsAB);
      calculateInterOffsets (interCoordinatesAB,currSuperInter,interOffsetsAB,1);
      pvalueAB = compareDistributions (intraOffsets,interOffsetsAB); 
      meanInterAB = calculateMedian (interOffsetsAB);
      destroyCoordinates (interCoordinatesAB);
      interCoordinatesBA = convertInterCoordinates (currSuperInter,0);
      arrayClear (interOffsetsBA);
      calculateInterOffsets (interCoordinatesBA,currSuperInter,interOffsetsBA,0);
      pvalueBA = compareDistributions (intraOffsets,interOffsetsBA);
      meanInterBA = calculateMedian (interOffsetsBA);
      destroyCoordinates (interCoordinatesBA); 
    }
    else {
      meanInterAB = -1;
      meanInterBA = -1;
      pvalueAB = -1;
      pvalueBA = -1;
    }
    superIntra1 = getSuperIntra (superIntras,currSuperInter->transcript1);
    superIntra2 = getSuperIntra (superIntras,currSuperInter->transcript2);
    exonCoordinates1 = hlr_strdup (exonCoordinates2string (currSuperInter->transcript1->subIntervals));
    exonCoordinates2 = hlr_strdup (exonCoordinates2string (currSuperInter->transcript2->subIntervals));
    printf ("%d\t%.2f\t%.2f\t%.5f\t%.5f\t%d\t%d\t%s\t",
            getNumberOfInters ( currSuperInter ),
	    meanInterAB,meanInterBA,pvalueAB,pvalueBA,
            superIntra1 ? getNumberOfIntras(superIntra1) : 0,
	    superIntra2 ? getNumberOfIntras(superIntra2) : 0,
            strEqual (currSuperInter->transcript1->chromosome,currSuperInter->transcript2->chromosome) ? "cis" : "trans");
    printf ("%s\t%d\t%s\t%s\t%c\t%d\t%d\t",
            currSuperInter->transcript1->name,
	    arrayMax (currSuperInter->transcript1->subIntervals),
	    exonCoordinates1,
	    currSuperInter->transcript1->chromosome,
	    currSuperInter->transcript1->strand,
	    currSuperInter->transcript1->start,
	    currSuperInter->transcript1->end);
    printf ("%s\t%d\t%s\t%s\t%c\t%d\t%d\t",
            currSuperInter->transcript2->name,
	    arrayMax (currSuperInter->transcript2->subIntervals),
	    exonCoordinates2,
	    currSuperInter->transcript2->chromosome,
	    currSuperInter->transcript2->strand,
	    currSuperInter->transcript2->start,
	    currSuperInter->transcript2->end);
    for (j = 0; j < arrayMax (currSuperInter->inters); j++) {
      currInter = arru (currSuperInter->inters,j,Inter*);
      printf ("%d,%d,%d,%d,%d,%d,%d%s",
              currInter->pairType,currInter->number1,currInter->number2,
              currInter->readStart1,currInter->readEnd1,
              currInter->readStart2,currInter->readEnd2, 
              j < arrayMax (currSuperInter->inters) - 1 ? "|" : "\t");
    }
    printf ("%s_%05d\t",argv[1],i + 1);
    for (j = 0; j < arrayMax (currSuperInter->inters); j++) {
      currInter = arru (currSuperInter->inters,j,Inter*);
      printf ("%s%s",currInter->read1,j < arrayMax (currSuperInter->inters) - 1 ? "|" : "\t");
    }
    for (j = 0; j < arrayMax (currSuperInter->inters); j++) {
      currInter = arru (currSuperInter->inters,j,Inter*);
      printf ("%s%s",currInter->read2,j < arrayMax (currSuperInter->inters) - 1 ? "|" : "\n");
    }
    hlr_free (exonCoordinates1);
    hlr_free (exonCoordinates2);
    i++;
  }    
  warn ("%s_numGfrEntries: %d",argv[0],i);
  stringPrintf (buffer,"%s.intraOffsets",argv[1]);
  fp = fopen (string (buffer),"w");
  for (i = 0; i < arrayMax (intraOffsets); i++) {
    fprintf (fp,"%d\n",arru (intraOffsets,i,int));
  }
  fclose (fp);
  stringDestroy (buffer);
  return EXIT_SUCCESS;
}

