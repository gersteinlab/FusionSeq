#include <bios/log.h>
#include <bios/format.h>

#include "gfr.h"

int main (int argc, char *argv[])
{
	GfrEntry *currGE;
	int offset;
	int count;
	int countRemoved;

	if (argc != 2) {
		usage ("%s <offset>",argv[0]);
	}
	count = 0;
	countRemoved = 0;
	offset = atoi (argv[1]);
	gfr_init ("-");
	puts (gfr_writeHeader ());
	while (currGE = gfr_nextEntry ()) {
		if (strEqual (currGE->fusionType,"cis") && 
		    currGE->strandTranscript1 != currGE->strandTranscript2 &&
		    (currGE->startTranscript2 - currGE->endTranscript1) < offset) {
			countRemoved++;
			continue;
		}
		puts (gfr_writeGfrEntry (currGE));
		count++;
	}
	gfr_deInit ();
	warn ("%s_offset: %d",argv[0],offset);
	warn ("%s_numRemoved: %d",argv[0],countRemoved);
	warn ("%s_numGfrEntries: %d",argv[0],count);
	return 0;
}

