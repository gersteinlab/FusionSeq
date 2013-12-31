#include <stdio.h>

#include <bios/log.h>
#include <bios/format.h>

#include "gfr.h"



int main (int argc, char *argv[])
{
	GfrEntry *currGE;
	int i;
	Stringa buffer;
	FILE *fp1,*fp2;
	int count;

	buffer = stringCreate (100);
	count = 0;
	gfr_init ("-");
	puts (gfr_writeHeader ());
	while (currGE = gfr_nextEntry ()) {
		stringPrintf (buffer,"%s_1.fasta",currGE->id);
		fp1 = fopen (string (buffer),"w");
		stringPrintf (buffer,"%s_2.fasta",currGE->id);
		fp2 = fopen (string (buffer),"w");
		if (fp1 == NULL || fp2 == NULL) {
			die ("Unable to open FASTA files");
		}
		for (i = 0; i < arrayMax (currGE->readsTranscript1); i++) {
			fprintf (fp1,">%s_1_%d\n%s\n",currGE->id,i + 1,textItem (currGE->readsTranscript1,i));
		}
		for (i = 0; i < arrayMax (currGE->readsTranscript2); i++) {
			fprintf (fp2,">%s_2_%d\n%s\n",currGE->id,i + 1,textItem (currGE->readsTranscript2,i));
		}
		fclose (fp1);
		fclose (fp2);
		puts (gfr_writeGfrEntry (currGE));
		count++;
	}
	gfr_deInit ();
	warn ("%s_numGfrEntries: %d",argv[0],count);
	return 0;
}

