#include <stdlib.h>

#include <bios/confp.h>
#include <bios/log.h>
#include <bios/format.h>
#include <bios/intervalFind.h>

#include "util.h"
#include "gfr.h"

int main (int argc, char *argv[])
{
	int i;
	GfrEntry* gfrE;
	config *conf;
  
	if ((conf = confp_open(getenv("FUSIONSEQ_CONFPATH"))) == NULL)
		return EXIT_FAILURE;

	Array intervals = arrayCreate( 50, Interval);
	Stringa buffer = stringCreate (100);

	stringPrintf(buffer,"%s/%s",
		     confp_get(conf, "ANNOTATION_DIR"),
		     confp_get(conf, "TRANSCRIPT_COMPOSITE_MODEL_FILENAME"));
	intervalFind_addIntervalsToSearchSpace (string (buffer),0);
  
	gfr_init("-");
	printf( "%s\n", gfr_writeHeader());
	while( gfrE = gfr_nextEntry()) {
		if (strEqual(gfrE->fusionType, "cis" ) && 
		    (gfrE->strandTranscript1 == gfrE->strandTranscript2)) {
			intervals = arrayCopy (intervalFind_getOverlappingIntervals (gfrE->chromosomeTranscript1, gfrE->endTranscript1+1, gfrE->startTranscript2-1 ));
			stringPrintf(buffer, "read-through");
			if (arrayMax(intervals)>0) {
				for (i = 0; i < arrayMax(intervals); i++) {
					Interval *currInterval = arru( intervals, i, Interval*);
					if ((currInterval->strand == gfrE->strandTranscript1)) { 
						stringPrintf(buffer,"intra");
						i = arrayMax(intervals);
					}
				}
			}
		} else {
			if (strEqual(gfrE->fusionType, "trans")) {
				stringPrintf( buffer, "inter"); 
			} else { 
				stringPrintf( buffer, "cis");
			} 
		}
		gfrE->fusionType = hlr_strdup(string(buffer));
		printf("%s\n", gfr_writeGfrEntry( gfrE ));
    
	}
	gfr_deInit();
	arrayDestroy(intervals);
	stringDestroy (buffer);
	confp_close(conf);

  return 0;
}
