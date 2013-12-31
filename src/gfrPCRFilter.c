#include <bios/log.h>
#include <bios/format.h>
#include "gfr.h"


int main (int argc, char *argv[])
{
	GfrEntry *currGE;
	int count;
	int countRemoved;
	int i;

	if (argc != 3) {
		usage ("%s <offsetCutoff> <minNumUniqueReads>",argv[0]);
	}
	count = 0;
	countRemoved = 0;

	int offsetCutOff = atoi (argv[1]);
	int minNumUniqueReads = atoi (argv[2]);

	gfr_init ("-");
	puts (gfr_writeHeader ());
	while (currGE = gfr_nextEntry ()) {
		Array starts = arrayCreate( 100, int);
		for (i = 0; i < arrayMax( currGE->interReads ); i++) {
			int currStart = arrp(currGE->interReads, i, GfrInterRead)->readStart1 + arrp(currGE->interReads, i, GfrInterRead)->readStart2;
			array(starts, arrayMax(starts), int) = currStart; 
		}
		arraySort( starts, (ARRAYORDERF) arrayIntcmp );
		arrayUniq( starts, NULL, (ARRAYORDERF) arrayIntcmp ) ;
		int numUniqeOffsets = arrayMax( starts );
		arrayDestroy( starts );

	if (arrayMax( currGE->readsTranscript1 ) != arrayMax( currGE->readsTranscript2 ) )
		die( "The two ends have a different number of reads");
	Texta reads = textCreate(arrayMax(currGE->readsTranscript1));
	for (i = 0; i < arrayMax(currGE->readsTranscript1); i++) {
		Stringa strA = stringCreate( strlen(textItem( currGE->readsTranscript1, i) ) * 2 + 1);
		stringAppendf( strA, textItem( currGE->readsTranscript1,i));
		stringAppendf( strA, textItem( currGE->readsTranscript2,i)); 
		textAdd( reads, string(strA));
		stringDestroy( strA );
	}
	textUniqKeepOrder( reads );
	int numRemaining = arrayMax( reads );
	textDestroy ( reads );

	if (numRemaining <= minNumUniqueReads || numUniqeOffsets <= offsetCutOff) {
		countRemoved++;
		continue;
	} 
	puts (gfr_writeGfrEntry (currGE));
	count++;
	}
	gfr_deInit ();
	warn("%s_PCRFilter: offset=%d minNumUniqueReads=%d",
	     argv[0],offsetCutOff, minNumUniqueReads);
	warn("%s_numRemoved: %d",argv[0],countRemoved);
	warn("%s_numGfrEntries: %d",argv[0],count);
	return 0;
}
