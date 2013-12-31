#include <math.h>

#include <bios/confp.h>
#include <bios/log.h>
#include <bios/format.h>
#include <bios/array.h>

#include "gfr.h"

double getNonExonicSPER( GfrEntry* currGE ) {
  int i;
  //int minArea1, minArea2, maxArea1, maxArea2;
  Array starts1 = arrayCreate( 50, int);
  Array starts2 = arrayCreate( 50, int);
  for( i = 0; i<arrayMax(currGE->interReads); i++) {
    array(starts1, arrayMax(starts1), int) = arrp(currGE->interReads, i, GfrInterRead)->readStart1; 
    array(starts2, arrayMax(starts2), int) = arrp(currGE->interReads, i, GfrInterRead)->readStart2; 
    /*if( i == 0 ) {
      minArea1 = currIR->readStart1; maxArea1 = currIR->readStart1;
      minArea2 = currIR->readStart2; maxArea2 = currIR->readStart2;
    }
    if(minArea1 > currIR->readStart1) minArea1 = currIR->readStart1;
    if(minArea2 > currIR->readStart2) minArea2 = currIR->readStart2;
    if(maxArea1 < currIR->readStart1) maxArea1 = currIR->readStart1;
    if(maxArea2 < currIR->readStart2) maxArea2 = currIR->readStart2;*/
  } 
  arraySort( starts1, (ARRAYORDERF) arrayIntcmp);
  arraySort( starts2, (ARRAYORDERF) arrayIntcmp);
  int tenthQuantile  = (int) round( arrayMax( starts1 ) * 0.10 );
  int ninetiethQuantile  = (int) round( arrayMax( starts1 ) * 0.90 );
  double area1 = (double) (arru( starts1, ninetiethQuantile , int ) - arru( starts1, tenthQuantile, int ) );
  double area2 = (double) (arru( starts2, ninetiethQuantile , int ) - arru( starts2, tenthQuantile, int ) );
  int nonExonicPEReads = (int) round( arrayMax( starts1 ) * 0.80 );
  arrayDestroy( starts1 ); 
  arrayDestroy( starts2 );
  //return MIN( (double) nonExonicPEReads / (double)(maxArea1 - minArea1), (double) nonExonicPEReads /  (double)(maxArea2 - minArea2 )); 
  return MIN( (double) nonExonicPEReads / (double)(area1), (double) nonExonicPEReads /  (double)( area2 )); 
}

int getTranscriptLength1( GfrEntry* currGE) {
  int i;
  int lengthTranscript=0;
  GfrExonCoordinate* currEC;
  for( i=0; i < arrayMax(currGE->exonCoordinatesTranscript1); i++) {
    currEC = arrp( currGE->exonCoordinatesTranscript1, i, GfrExonCoordinate);
    lengthTranscript += currEC->end - currEC->start;
  }
  return lengthTranscript;
}
int getTranscriptLength2( GfrEntry* currGE) {
  int i;
  int lengthTranscript = 0;
  GfrExonCoordinate* currEC;
  for( i=0; i < arrayMax(currGE->exonCoordinatesTranscript2); i++) {
    currEC = arrp( currGE->exonCoordinatesTranscript2, i, GfrExonCoordinate);
    lengthTranscript += currEC->end - currEC->start;
  }
  return lengthTranscript;
}

int main (int argc, char *argv[])
{
  GfrEntry *currGE;
  int count;
  int countRemoved;
  int exonicPEReads, nonExonicPEReads,lengthTranscript1,lengthTranscript2;
  int i;
  GfrPairCount* currGEPC;
  double nonExonicSPER;

  count = 0;
  countRemoved = 0;
  gfr_init ("-");
  puts (gfr_writeHeader ());
  while (currGE = gfr_nextEntry () ) {
    exonicPEReads=0; 
    nonExonicPEReads = 0;
    int allIntronic = 1;
  
    for( i=0; i<arrayMax( currGE->pairCounts ); i++ ) {
      currGEPC = arrp( currGE->pairCounts, i, GfrPairCount);
      if( currGEPC->pairType == GFR_PAIR_TYPE_EXONIC_EXONIC ) {
	exonicPEReads+= currGEPC->count; 
	allIntronic = 0;
      } /*else {
	nonExonicPEReads+= currGEPC->count;
	}//*/
    }   
    if( allIntronic == 1 ) {
      nonExonicSPER = getNonExonicSPER( currGE );
      lengthTranscript1 = getTranscriptLength1( currGE );
      lengthTranscript2 = getTranscriptLength2( currGE );
      if(   ( nonExonicSPER > ( (double) currGE->numIntra1 / (double) lengthTranscript1 ) ) ||
	    ( nonExonicSPER > ( (double) currGE->numIntra2 / (double) lengthTranscript2 ) ) ) {
	countRemoved++;
	continue;
      } else {
	puts (gfr_writeGfrEntry (currGE));
	count++;
      }
    } else {
      puts (gfr_writeGfrEntry (currGE));
      count++;
    }//*/
  }
  gfr_deInit ();
  warn ("%s_numRemoved: %d",argv[0],countRemoved);
  warn ("%s_numGfrEntries: %d",argv[0],count);
  return 0;
}
