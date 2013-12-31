#include <stdlib.h>
#include <math.h>


#include <bios/log.h>
#include <bios/format.h>
#include <bios/intervalFind.h>
#include <bios/linestream.h>
#include <bios/common.h>

#include "gfr.h"


#define MILLION 1000000.0

typedef struct {
	char* geneID;
	double value;
} GeneExpression;


static int sortGeneByID( GeneExpression* a, GeneExpression* b )
{
	return strcmp (b->geneID, a->geneID);
}

int main (int argc, char *argv[])
{
	int i ;
	char *line;
	Array geneExpression = arrayCreate(25000, GeneExpression);  

	if (argc != 2) {
		usage ("Usage: %s <prefix>\nNote:prefix.composite.expression is supposed to exist.",argv[0]);
	}
	Stringa buffer=stringCreate( 50 );
	stringPrintf( buffer, "%s.composite.expression", argv[1]);
	LineStream ls = ls_createFromFile( string(buffer)  );
	while (line = ls_nextLine(ls)) {
		WordIter w = wordIterCreate( line, "\t", 1);
		GeneExpression* currGE = arrayp( geneExpression, arrayMax( geneExpression ), GeneExpression);
		currGE->geneID = hlr_strdup( wordNext(w) );
		currGE->value = atof( wordNext(w) );
		wordIterDestroy( w );
	}
	arraySort( geneExpression, (ARRAYORDERF)sortGeneByID );
	gfr_init ( "-" );
	GfrEntry* gfrE;
	while (gfrE = gfr_nextEntry()) { 
		int idx1, idx2;
		GeneExpression* toSearchGE ;
		AllocVar(toSearchGE);
		toSearchGE->geneID = hlr_strdup( gfrE->nameTranscript1);
		arrayFind( geneExpression, toSearchGE, &idx1, (ARRAYORDERF)sortGeneByID );
		hlr_free( toSearchGE->geneID );
		toSearchGE->geneID = hlr_strdup( gfrE->nameTranscript2);
		arrayFind( geneExpression, toSearchGE, &idx2, (ARRAYORDERF)sortGeneByID );
		printf ("%s\t%f\n", 
			gfrE->id, 
			(double) gfrE->numInter / (arrp(geneExpression, idx1, GeneExpression)->value + 
				 		   arrp( geneExpression, idx2, GeneExpression)->value) / 2.0 );    
		hlr_free( toSearchGE->geneID );
		freeMem( toSearchGE );
	}
	gfr_deInit();
	for (i = 0; i < arrayMax( geneExpression ); i++)
		hlr_free( arrp( geneExpression, i, GeneExpression)->geneID);
	arrayDestroy(geneExpression);
	return 0;
}

