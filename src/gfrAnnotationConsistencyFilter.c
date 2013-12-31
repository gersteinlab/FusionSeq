#include <bios/log.h>
#include <bios/format.h>

#include "gfr.h"



int main (int argc, char *argv[])
{
	GfrEntry *currGE;
	int count;
	int countRemoved;
 
	if (argc != 2) {
		usage ("%s <string>",argv[0]);
	}
	count = 0;
	countRemoved = 0;
	gfr_init ("-");
	puts (gfr_writeHeader ());
	while (currGE = gfr_nextEntry ()){
		if (currGE->descriptionTranscript1 == NULL ||
				currGE->descriptionTranscript2 == NULL) {
			die ("Transcript description is missing");
		}
		if (strCaseStr (currGE->descriptionTranscript1,argv[1]) ||
		    strCaseStr (currGE->descriptionTranscript2,argv[1])) {
			countRemoved++;
			continue;
		}
		puts (gfr_writeGfrEntry (currGE));
		count++;
	}
	gfr_deInit ();
	warn ("%s_string: %s",argv[0],argv[1]);
	warn ("%s_numRemoved: %d",argv[0],countRemoved);
	warn ("%s_numGfrEntries: %d",argv[0],count);
	return 0;
}

