#include <bios/log.h>
#include <bios/format.h>
#include <bios/bowtieParser.h>
#include <bios/common.h>

#include "bp.h"


static int sortBowtieQueries ( BowtieQuery *a, BowtieQuery *b) 
{
  return strcmp ( a->sequenceName, b->sequenceName );
}

int main (int argc, char *argv[])
{
  Array breakPoints;
  Array bowtieQueries;
  BreakPoint *currBP;
  Stringa buffer = stringCreate( 20 );
  int i;
  if( argc<2 ) {
    usage("Usage: %s <prefix>", argv[0]);
  }
  bp_init("-");
  breakPoints = bp_getBreakPoints();
  stringPrintf( buffer, "%s", argv[1]);
  bowtieParser_initFromFile (string(buffer));
  bowtieQueries = bowtieParser_getAllQueries ();
  bowtieParser_deInit (); 
  arraySort( bowtieQueries, (ARRAYORDERF) sortBowtieQueries);

  for( i=0; i<arrayMax( breakPoints ); i++ ) {
    currBP = arrp( breakPoints, i, BreakPoint );
    BowtieQuery *test;
    int pos;
    AllocVar( test );
    stringPrintf( buffer, "%s|%s", currBP->tileCoordinate1, currBP->tileCoordinate2 );
    test->sequenceName = hlr_strdup( string( buffer) );
    int found=arrayFind( bowtieQueries, test, &pos, (ARRAYORDERF) sortBowtieQueries );
    if( !found ) {
      printf( "%s\n", bp_writeBreakPoint( currBP ) );
    }
    hlr_free( test->sequenceName);
    freeMem( test );
  }
  stringDestroy( buffer );
  bp_deInit();  
  return 0;
}
