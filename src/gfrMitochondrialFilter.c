#include <bios/log.h>
#include <bios/format.h>
#include "gfr.h"



int main (int argc, char *argv[])
{
	GfrEntry *currGE;
	int count;
	int countRemoved;
 
	count = 0;
	countRemoved = 0;
	gfr_init ("-");
	puts (gfr_writeHeader ());
	while (currGE = gfr_nextEntry ()) {
		if (strEqual(currGE->chromosomeTranscript1, "chrM") || 
		    strEqual(currGE->chromosomeTranscript2, "chrM")) {
			countRemoved++;
			continue;
		} else {
			puts (gfr_writeGfrEntry (currGE));
			count++;
		}
	}
	gfr_deInit ();
	warn ("%s_numRemoved: %d",argv[0],countRemoved);
	warn ("%s_numGfrEntries: %d",argv[0],count);
	return 0;
}

