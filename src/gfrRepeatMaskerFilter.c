#include <bios/log.h>
#include <bios/format.h>
#include <bios/intervalFind.h>

#include "gfr.h"

int main (int argc, char *argv[])
{
	GfrEntry *currGE;
	int count,countRemoved;
	int i;
	Array intervals;
	GfrInterRead *currGIR;
	int totalOverlaps;
	int minNumInterReads;

	if (argc != 3) {
		usage ("%s <repeatMasker.interval> <minNumInterReads>",argv[0]);
	}
	intervalFind_addIntervalsToSearchSpace (argv[1],0);
	minNumInterReads = atoi (argv[2]);
	count = 0;
	countRemoved = 0;
	gfr_init ("-");
	puts (gfr_writeHeader ());
	while (currGE = gfr_nextEntry ()){
		for (i = 0; i < arrayMax (currGE->interReads); i++) {
			currGIR = arrp (currGE->interReads,i,GfrInterRead);
			if (currGIR->pairType == GFR_PAIR_TYPE_EXONIC_EXONIC) {
				continue;
			}
			totalOverlaps = 0;
			intervals = intervalFind_getOverlappingIntervals (currGE->chromosomeTranscript1,currGIR->readStart1,currGIR->readEnd1);
			totalOverlaps += arrayMax (intervals);
			intervals = intervalFind_getOverlappingIntervals (currGE->chromosomeTranscript2,currGIR->readStart2,currGIR->readEnd2);
			totalOverlaps += arrayMax (intervals);
			if (totalOverlaps > 0) {
				currGE->numInter--;
				currGIR->flag = 1;
			}
		}
		if (currGE->numInter < minNumInterReads) {
			countRemoved++;
			continue;
		}
		puts (gfr_writeGfrEntry (currGE));
		count++;
	}
	gfr_deInit ();
	warn ("%s_numRemoved: %d",argv[0],countRemoved);
	warn ("%s_numGfrEntries: %d",argv[0],count);
	return 0;
}

