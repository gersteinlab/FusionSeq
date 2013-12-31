#include <cstdio>
#include <cstdlib>
#include <ctime>

extern "C" {
#include <bios/log.h>
#include <bios/array.h>
#include <bios/format.h>
#include "bp.h"
}

#ifndef __CINT__

//--- ROOT includes ---

#include <root/TROOT.h>
#include <root/TMath.h>
#include <root/TRandom.h>

TROOT root("Rint","The ROOT Interactive Interface");

#endif

int main(int argc, char *argv[])
{
	Array breakPoints;
	BreakPoint *currBP;
	BreakPointRead *currBPR;
	int minNumReads, minNumUniqueOffsets,
	    minNumReadsForKS,numPossibleOffsets;
	double pValueCutoffForKS;
	Array offsets;
	Array randomNumbers;
	double *observedOffsets;
	double *randomOffsets;

 	if (argc != 6) 
		usage("%s <minNumReads> <minNumUniqueOffsets> <minNumReadsForKS> <pValueCutoffForKS> <numPossibleOffsets>", argv[0]);
	
	minNumReads         = std::atoi(argv[1]);
	minNumUniqueOffsets = std::atoi(argv[2]);
	minNumReadsForKS    = std::atoi(argv[3]);
	pValueCutoffForKS   = std::atof(argv[4]);
	numPossibleOffsets  = std::atoi(argv[5]);
	bp_init("-");
	offsets = arrayCreate(100, int);
	randomNumbers = arrayCreate(100, int);
	breakPoints = bp_getBreakPoints();

	for (int i = 0; i < arrayMax(breakPoints); i++) {
		currBP = arrp(breakPoints, i, BreakPoint);
		arrayClear(offsets);
		for (int j = 0; j < arrayMax(currBP->breakPointReads); j++) {
			currBPR = arrp(currBP->breakPointReads, j, BreakPointRead);
			array(offsets, arrayMax(offsets), int) = currBPR->offset;
		}
		arraySort(offsets, (ARRAYORDERF) arrayIntcmp);
		arrayUniq(offsets, NULL, (ARRAYORDERF) arrayIntcmp);
		if (arrayMax(currBP->breakPointReads) >= minNumReads && 
		    arrayMax(currBP->breakPointReads) < minNumReadsForKS) {        
			if (arrayMax(offsets) >= minNumUniqueOffsets)
				std::puts(bp_writeBreakPoint(currBP));
		}
		else if (arrayMax(currBP->breakPointReads) >= minNumReads && 
			 arrayMax(currBP->breakPointReads) >= minNumReadsForKS) {
			arrayClear(randomNumbers);
			for (int j = 0; j < arrayMax(offsets); j++)
				array (randomNumbers, arrayMax(randomNumbers), int) = std::rand() % numPossibleOffsets;
			
			arraySort(randomNumbers, (ARRAYORDERF) arrayIntcmp);
			observedOffsets = (double *) hlr_malloc(arrayMax(offsets) * sizeof(double)); 
			randomOffsets = (double *) hlr_malloc(arrayMax(offsets) * sizeof(double)); 
			for (int j = 0; j < arrayMax(offsets); j++) {
				observedOffsets[j] = arru(offsets,j,int);
				randomOffsets[j] = arru(randomNumbers,j,int);
			}
			if (pValueCutoffForKS < TMath::KolmogorovTest(arrayMax(offsets), 
								      observedOffsets, 
								      arrayMax(offsets), 
								      randomOffsets, 
								      ""))
				std::puts(bp_writeBreakPoint(currBP));
			
			hlr_free(observedOffsets);
			hlr_free(randomOffsets);
		}
	}
	bp_deInit();
	return 0;
}

