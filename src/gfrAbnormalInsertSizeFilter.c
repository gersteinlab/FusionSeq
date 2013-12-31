#include <bios/log.h>
#include <bios/format.h>
#include "gfr.h"



int main (int argc, char *argv[])
{
	GfrEntry *currGE;
	double pvalueCutOff;
	int count;
	int countRemoved;

	if (argc != 2) {
		usage ("%s <pvalueCutOff>",argv[0]);
	}
	count = 0;
	countRemoved = 0;
	pvalueCutOff = atof (argv[1]);
	gfr_init ("-");
	puts (gfr_writeHeader ());
	while (currGE = gfr_nextEntry ()){
		if ((MAX(currGE->pValueAB,currGE->pValueBA) < pvalueCutOff) && 
		    (currGE->pValueAB != -1) ) {
			countRemoved++;
			continue;
		} 
		puts (gfr_writeGfrEntry (currGE));
		count++;
	}
	gfr_deInit ();
	warn ("%s_pvalueCutOff: %f",argv[0],pvalueCutOff);
	warn ("%s_numRemoved: %d",argv[0],countRemoved);
	warn ("%s_numGfrEntries: %d",argv[0],count);
	return 0;
}

