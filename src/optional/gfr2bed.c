#include <stdio.h>

#include <bios/log.h>
#include <bios/format.h>

#include "gfr.h"



int main (int argc, char *argv[])
{
	GfrEntry *currGE;
	GfrInterRead *currGIR;
	int i;
	Stringa buffer;
	FILE *fp1,*fp2;
	int count;

	count = 0;
	buffer = stringCreate (100);
	gfr_init ("-");
	puts (gfr_writeHeader ());
	while (currGE = gfr_nextEntry ()) {
		stringPrintf (buffer,"%s_1.bed",currGE->id);
		fp1 = fopen (string (buffer),"w");
		stringPrintf (buffer,"%s_2.bed",currGE->id);
		fp2 = fopen (string (buffer),"w");
		if (fp1 == NULL || fp2 == NULL) {
			die ("Unable to open BED files");
		}
		fprintf (fp1,"browser full knownGene\n");
		fprintf (fp1,"track name=\"Inter paird-ends: %s_1\" visibility=2\n",currGE->id);
		fprintf (fp2,"browser full knownGene\n");
		fprintf (fp2,"track name=\"Inter paird-ends: %s_2\" visibility=2\n",currGE->id);
		for (i = 0; i < arrayMax (currGE->interReads); i++) {
			currGIR = arrp (currGE->interReads,i,GfrInterRead);
			fprintf (fp1,"%s\t%d\t%d\n",currGE->chromosomeTranscript1,currGIR->readStart1,currGIR->readEnd1);
			fprintf (fp2,"%s\t%d\t%d\n",currGE->chromosomeTranscript2,currGIR->readStart2,currGIR->readEnd2);
		}
		fclose (fp1);
		fclose (fp2);
		puts (gfr_writeGfrEntry (currGE));
		count++;
	}
	gfr_deInit ();
	stringDestroy (buffer);
	warn ("%s_numGfrEntries: %d",argv[0],count);
	return 0;
}

