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
	FILE *fp;
	int count;

	count = 0;
	buffer = stringCreate (100);
	gfr_init ("-");
	puts (gfr_writeHeader ());
  
	while (currGE = gfr_nextEntry ()){
		count++;
		if (strEqual (currGE->fusionType,"inter") || 
		    strEqual (currGE->fusionType,"trans") ) {
			puts (gfr_writeGfrEntry (currGE));
			continue;
		}
		stringPrintf (buffer,"%s.gff",currGE->id);
		fp = fopen (string (buffer),"w");
		if (fp == NULL) {
			die ("Unable to open GFF file");
		}
		fprintf (fp,"browser pack knownGene\n");
		fprintf (fp,"track name=\"%s\" visibility=2\n",currGE->id);
		for (i = 0; i < arrayMax (currGE->interReads); i++) {
			currGIR = arrp (currGE->interReads,i,GfrInterRead);     
			fprintf (fp, "%s\tRNA_Seq\texon\t%d\t%d\t.\t.\t.\tTG%d\n",
				currGE->chromosomeTranscript1,
				currGIR->readStart1,currGIR->readEnd1,i);
			fprintf (fp, "%s\tRNA_Seq\texon\t%d\t%d\t.\t.\t.\tTG%d\n",
				currGE->chromosomeTranscript2,
				currGIR->readStart2,currGIR->readEnd2,i);
		}
		fclose (fp);
		puts (gfr_writeGfrEntry (currGE));
	}
	gfr_deInit ();
	stringDestroy (buffer);
	warn ("%s_numGfrEntries: %d",argv[0],count);
	return 0;
}

