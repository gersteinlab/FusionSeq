#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include <bios/log.h>
#include <bios/format.h>
#include <bios/fasta.h>
#include <bios/confp.h>

#include "bp.h"


#define TILE_SEPARATOR "  "

static config *Conf = NULL;

static char* getBreakPointSequence (char *tileCoordinate1, char *tileCoordinate2)
{
	Stringa buffer;
	Stringa targetsFile;
	FILE *fp;
	Array targetSeqs;
	int i;
	Seq *currSeq;
	static Stringa sequence = NULL;

	buffer = stringCreate (100);
	targetsFile = stringCreate (100);
	stringPrintf (targetsFile,"targets_%d.txt",getpid ());
	if (!(fp = fopen (string (targetsFile),"w")) ){
		die ("Unable to open target file: %s",string (targetsFile));
	}
	fprintf (fp,"%s\n%s",tileCoordinate1,tileCoordinate2);
	fclose (fp);

	stringPrintf (buffer,"%s %s/%s stdout -noMask -seqList=%s",
		      confp_get(Conf, "BLAT_TWO_BIT_TO_FA"),
		      confp_get(Conf, "BLAT_DATA_DIR"),
		      confp_get(Conf, "BLAT_TWO_BIT_DATA_FILENAME"),
		      string (targetsFile));
	fasta_initFromPipe (string (buffer));
	targetSeqs = fasta_readAllSequences (0);
	fasta_deInit ();
	if (arrayMax (targetSeqs) != 2) {
		die ("Expected only two target sequences");
	} 
	stringCreateClear (sequence,100);
	for (i = 0; i < arrayMax (targetSeqs); i++) {
		currSeq = arrp (targetSeqs,i,Seq);
		stringAppendf (sequence,"%s",currSeq->sequence);
		hlr_free (currSeq->name);
		hlr_free (currSeq->sequence);
	}
	arrayDestroy (targetSeqs);
	stringPrintf (buffer,"rm -rf %s",string (targetsFile));
	hlr_system (string (buffer),0);
	stringDestroy (targetsFile);
	stringDestroy (buffer);
	return string (sequence);
}



static void extractCoordinates (char *string, int *start, int *end)
{
	static char* stringCopy = NULL; 
	char *pos1,*pos2;

	strReplace (&stringCopy,string);
	pos1 = strchr (stringCopy,':');
	pos2 = strchr (stringCopy,'-');
	*pos2 = '\0';
	*start = atoi (pos1 + 1);
	*end = atoi (pos2 + 1);
}



static int getTileSize (char *tileCoordinate1, char *tileCoordinate2)
{
	int startTile1,endTile1,tileSize1;
	int startTile2,endTile2,tileSize2;

	extractCoordinates (tileCoordinate1,&startTile1,&endTile1);
	tileSize1 = endTile1 - startTile1;
	extractCoordinates (tileCoordinate2,&startTile2,&endTile2);
	tileSize2 = endTile2 - startTile2;
	if (tileSize1 != tileSize2) {
		//die ("Expected the same tileSize for both tile coordinates: %s %s",tileCoordinate1,tileCoordinate2);
	}

		int junctionSize= tileSize1 + tileSize2;
	return (junctionSize/2 + junctionSize%2);
}



int main (int argc, char *argv[])
{
	Array breakPoints;
	BreakPoint *currBP;
	BreakPointRead *currBPR;
	int i,j,k;
	int readLength;
	int tileSize;
	char *breakPointSequence;

	if ((Conf = confp_open(getenv("FUSIONSEQ_CONFPATH"))) == NULL)
		return EXIT_FAILURE;

	bp_init ("-");
	breakPoints = bp_getBreakPoints ();
	for (i = 0; i < arrayMax (breakPoints); i++) {
		currBP = arrp (breakPoints,i,BreakPoint);
		tileSize = getTileSize (currBP->tileCoordinate1,currBP->tileCoordinate2);
		breakPointSequence = getBreakPointSequence (currBP->tileCoordinate1,currBP->tileCoordinate2);
		printf ("Tile 1: %s\n",currBP->tileCoordinate1);
		printf ("Tile 2: %s\n",currBP->tileCoordinate2);
		printf ("Number of reads spanning breakpoint: %d\n\n\n",arrayMax (currBP->breakPointReads));
		for (j = 0; j < arrayMax (currBP->breakPointReads); j++) {
			currBPR = arrp (currBP->breakPointReads,j,BreakPointRead);
			readLength = strlen (currBPR->read);
			for (k = 0; k < currBPR->offset; k++) {
				printf (" ");
			}
			for (k = 0; k < readLength; k++) {
				if (((currBPR->offset + k) % tileSize) == 0 && (currBPR->offset + k) != 0) {
					printf ("%s",TILE_SEPARATOR);
				}
				printf ("%c",currBPR->read[k]);
			}
			printf ("\n");
		}
		for (k = 0; k < (2 * tileSize); k++) {
			if ((k % tileSize) == 0 && k != 0) {
				printf ("%s",TILE_SEPARATOR);
			}
			printf ("%c",breakPointSequence[k]);
		}
		printf ("\n\n\n\n\n");
	}
	bp_deInit ();

	confp_close(Conf);

	return EXIT_SUCCESS;
}

 
