#include <stdlib.h>
#include <time.h>
#include <math.h>

#include <bios/log.h>
#include <bios/format.h>
#include <bios/intervalFind.h>
#include <bios/numUtil.h>
#include <bios/linestream.h>
#include <mrf/mrf.h>

#include "gfr.h"

#define MILLION 1000000.0


static int sortArrayByDASPER( GfrEntry* a, GfrEntry* b ) {
  return (int) (b->DASPER*MILLION - a->DASPER*MILLION);
}

int main (int argc, char *argv[])
{
  int i,j ;
  double inter, intra1, intra2, intras;
  double Npe = -1.0;
  char *line;

  Array gfrA = arrayCreate(20, GfrEntry); 
 
  if (argc != 2) {
    usage ("Usage: %s <prefix>\nNote:prefix.meta is supposed to exist.",argv[0]);
  }

  // reading GFR file 
  gfr_init ( "-" );
  gfr_addNewColumnType("SPER");
  gfr_addNewColumnType("DASPER");
  gfr_addNewColumnType("RESPER");
  gfrA = gfr_parse();

  Stringa buffer=stringCreate( 50 );
  stringPrintf( buffer, "%s.meta", argv[1]);
  LineStream ls = ls_createFromFile( string(buffer)  );
  if( ls == NULL ) {
    die( "File not found: (%s.meta)", argv[1]);
  } else {
    while( line = ls_nextLine( ls ) ) {
      WordIter w = wordIterCreate( line, "\t", 1);
      if( strEqual( wordNext ( w ), "Mapped_reads") ) {
	Npe = atof( wordNext( w ) );
      }
      wordIterDestroy( w );
    }
  }
  if( Npe <= 0 ) die("Number of mapped reads cannot be less or equal than zero: %1.3f", Npe);
  ls_destroy(ls);
  double expSPER;
  double resperAvg = 0.0;
  for( j = 0; j < arrayMax(gfrA); j++) resperAvg += (double) arrp( gfrA, j, GfrEntry)->numInter;
  resperAvg /= ( arrayMax(gfrA) * Npe );
  GfrEntry* gfrE;

  for( i=0; i<arrayMax(gfrA); i++) {
    gfrE = arrp( gfrA, i, GfrEntry );
    gfrE->SPER    = ( ((double) gfrE->numInter ) / Npe ) * MILLION;
    
    inter  = (double) gfrE->numInter;
    intra1 = (double) gfrE->numIntra1; 
    intra2 = (double) gfrE->numIntra2; 
    intras = intra1*intra2;
    
    expSPER = ( 2.0 * intra1 + inter) * (2.0 * intra2 + inter) / (4.0 * Npe * Npe ) * MILLION;
    gfrE->DASPER = (gfrE->SPER  - expSPER);
    
    gfrE->RESPER = gfrE->SPER / resperAvg / MILLION;

  }
  arraySort( gfrA, (ARRAYORDERF)sortArrayByDASPER );
  printf("%s\n", gfr_writeHeader() );
  for( i=0; i<arrayMax(gfrA); i++) {
    GfrEntry* gfrE =  arrp( gfrA, i, GfrEntry );
    printf("%s\n", gfr_writeGfrEntry( gfrE )) ;
  }
  gfr_deInit();
  
  warn("%s_numberMappedReads: %1.0f", argv[0], Npe);
  warn("%s_numGfrEntries: %d", argv[0], arrayMax( gfrA ));
  arrayDestroy(gfrA);
  return 0;
}

