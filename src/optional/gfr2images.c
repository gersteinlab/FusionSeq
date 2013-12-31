
#include <gd.h>
#include <gdfonts.h>
#include <gdfontt.h>
#include <gdfontmb.h>

#include <bios/log.h>
#include <bios/format.h>

#include "gfr.h"


#define IMAGE_HEIGHT 200
#define MIN_IMAGE_WIDTH 600
#define SIDE_MARGIN 10
#define Y_COORDINATE_TRANSCRIPT_1 50
#define Y_COORDINATE_TRANSCRIPT_2 150
#define SPACER 12 // number of pixels between exons
#define EXON_SIZE 6 // number of pixels along the side of a square 



static void drawTranscript (int numberExons, int yCoordinateTranscript, gdImagePtr im, int color)
{
	int i;
	int xEndCoordinateTranscript;

	xEndCoordinateTranscript = SIDE_MARGIN + 
				   numberExons * EXON_SIZE + 
				   (numberExons + 1) * SPACER;
	gdImageLine(im,
		    SIDE_MARGIN,
		    yCoordinateTranscript,
		    xEndCoordinateTranscript,
		    yCoordinateTranscript,
		    color);
	for (i = 1; i <= numberExons; i++) {
		gdImageFilledRectangle(im,
				       SIDE_MARGIN + i * SPACER + (i - 1) * EXON_SIZE,
				       yCoordinateTranscript - EXON_SIZE / 2,
				       SIDE_MARGIN + i * SPACER + (i - 1) * EXON_SIZE + EXON_SIZE,
				       yCoordinateTranscript + EXON_SIZE / 2,
				       color);
	}
}

void drawNumbers(gdImagePtr im, int numExons, int color)
{
	int i;
	Stringa buffer = stringCreate(3);
	for(i = 1; i <= numExons; i++) {
		stringPrintf( buffer, "%d", i);
		gdImageString(im, 
			      gdFontGetTiny(), 
			      SIDE_MARGIN + i * SPACER  + EXON_SIZE * i - EXON_SIZE / 2 - strlen(string(buffer)) * gdFontGetTiny()->w / 2,  
			      im->sy / 2,  
			      (unsigned char*) string(buffer), color);
	}
	stringDestroy(buffer);
}

static void makeImage (GfrEntry *currGE)
{
	gdImagePtr im;
	FILE *jpegout;
	int black, white, gray, blue, orange, green, color;
	int maxNumberExons;
	int imageWidth;
	int i;
	int xCoord1, xCoord2;
	//GfrInterRead *currGIR;
	static Stringa buffer = NULL;

	maxNumberExons = MAX (currGE->numExonsTranscript1,currGE->numExonsTranscript2);
	imageWidth = MAX (maxNumberExons * EXON_SIZE + (maxNumberExons + 1) * SPACER + 2 * SIDE_MARGIN,
			  MIN_IMAGE_WIDTH);
	im = gdImageCreate(imageWidth,IMAGE_HEIGHT);
	white = gdImageColorAllocate(im,255,255,255);
	gray  = gdImageColorAllocate(im,192,192,192);
	blue  = gdImageColorAllocate(im,0,0,255);
	black = gdImageColorAllocate(im,0,0,0);
	orange = gdImageColorAllocate(im,255,69,0);
	green = gdImageColorAllocate(im, 46,139,87);

	gdImageString (im,gdFontGetMediumBold (),SIDE_MARGIN,Y_COORDINATE_TRANSCRIPT_1 - 30,(unsigned char*)currGE->geneSymbolTranscript1,black);
	gdImageString (im,gdFontGetMediumBold (),SIDE_MARGIN,Y_COORDINATE_TRANSCRIPT_2 + 18,(unsigned char*)currGE->geneSymbolTranscript2,black);
	drawTranscript (currGE->numExonsTranscript1,Y_COORDINATE_TRANSCRIPT_1,im,black);
	drawTranscript (currGE->numExonsTranscript2,Y_COORDINATE_TRANSCRIPT_2,im,black);
	drawNumbers ( im, MAX( currGE->numExonsTranscript1, currGE->numExonsTranscript2), black);
	for (i = 0; i< arrayMax( currGE->pairCounts); i++) {
		GfrPairCount* currGEPC = arrp( currGE->pairCounts, i, GfrPairCount);
		switch (currGEPC->pairType) {
		case GFR_PAIR_TYPE_EXONIC_EXONIC:
			xCoord1 = SIDE_MARGIN + currGEPC->number1 * SPACER + currGEPC->number1 * EXON_SIZE - EXON_SIZE / 2;
			xCoord2 = SIDE_MARGIN + currGEPC->number2 * SPACER + currGEPC->number2 * EXON_SIZE - EXON_SIZE / 2;
			color = orange;
			break;
		case GFR_PAIR_TYPE_EXONIC_INTRONIC:
			xCoord1 = SIDE_MARGIN + currGEPC->number1 * SPACER + currGEPC->number1 * EXON_SIZE - EXON_SIZE / 2;
			xCoord2 = SIDE_MARGIN + currGEPC->number2 * SPACER + currGEPC->number2 * EXON_SIZE + SPACER / 2;
			color = orange;
			break;
		case GFR_PAIR_TYPE_EXONIC_JUNCTION:
			xCoord1 = SIDE_MARGIN + currGEPC->number1 * SPACER + currGEPC->number1 * EXON_SIZE - EXON_SIZE / 2;
			xCoord2 = SIDE_MARGIN + (currGEPC->number2/2) * SPACER + (currGEPC->number2)/2 * EXON_SIZE + (currGEPC->number2 % 2) * SPACER;
			color = orange;
			break;    
		case GFR_PAIR_TYPE_INTRONIC_EXONIC:
			xCoord1 = SIDE_MARGIN + currGEPC->number1 * SPACER + currGEPC->number1 * EXON_SIZE + SPACER / 2;
			xCoord2 = SIDE_MARGIN + currGEPC->number2 * SPACER + currGEPC->number2 * EXON_SIZE - EXON_SIZE / 2;
			color = green;
			break;
		case GFR_PAIR_TYPE_INTRONIC_INTRONIC:
			xCoord1 = SIDE_MARGIN + currGEPC->number1 * SPACER + currGEPC->number1 * EXON_SIZE + SPACER / 2;
			xCoord2 = SIDE_MARGIN + currGEPC->number2 * SPACER + currGEPC->number2 * EXON_SIZE + SPACER / 2;
			color = green;
			break;
		case GFR_PAIR_TYPE_INTRONIC_JUNCTION:
			xCoord1 = SIDE_MARGIN + currGEPC->number1 * SPACER + currGEPC->number1 * EXON_SIZE + SPACER / 2;
			xCoord2 = SIDE_MARGIN + (currGEPC->number2/2) * SPACER + (currGEPC->number2)/2 * EXON_SIZE + (currGEPC->number2 % 2) * SPACER;
			color = green;
			break;
		case GFR_PAIR_TYPE_JUNCTION_JUNCTION:
			xCoord1 = SIDE_MARGIN + (currGEPC->number1/2) * SPACER + (currGEPC->number1) / 2 * EXON_SIZE + (currGEPC->number1 % 2) * SPACER;
			xCoord2 = SIDE_MARGIN + (currGEPC->number2/2) * SPACER + (currGEPC->number2) / 2 * EXON_SIZE + (currGEPC->number2 % 2) * SPACER;
			color = gray;
			break;
		case GFR_PAIR_TYPE_JUNCTION_EXONIC:
			xCoord1 = SIDE_MARGIN + (currGEPC->number1/2) * SPACER + (currGEPC->number1) / 2 * EXON_SIZE + (currGEPC->number1 % 2) * SPACER;
			xCoord2 = SIDE_MARGIN + currGEPC->number2 * SPACER + currGEPC->number2 * EXON_SIZE - EXON_SIZE / 2;
			color = gray;
			break;     
		case GFR_PAIR_TYPE_JUNCTION_INTRONIC:
			xCoord1 = SIDE_MARGIN + (currGEPC->number1/2) * SPACER + (currGEPC->number1) / 2 * EXON_SIZE + (currGEPC->number1 % 2) * SPACER;
			xCoord2 = SIDE_MARGIN + currGEPC->number2 * SPACER + currGEPC->number2 * EXON_SIZE + SPACER / 2;
			color = gray;
			break;     
		default:
			xCoord1 = SIDE_MARGIN;
			xCoord2 = SIDE_MARGIN;
			color = white;
			break;
		}
		gdImageLine(im,
			    xCoord1,
			    Y_COORDINATE_TRANSCRIPT_1,
			    xCoord2,
			    Y_COORDINATE_TRANSCRIPT_2,
			    color);
	}
	/*
	for (i = 0; i < arrayMax (currGE->interReads); i++) {
	 	currGIR = arrp (currGE->interReads,i,GfrInterRead);
		gdImageLine(im,
			    xCoord1,
			    Y_COORDINATE_TRANSCRIPT_1,
			    xCoord2,
			    Y_COORDINATE_TRANSCRIPT_2,
			    color);
		gdImageLine(im,
			    SIDE_MARGIN + currGIR->exonNumber1 * SPACER + (currGIR->exonNumber1 - 1) * EXON_SIZE + EXON_SIZE / 2,
			    Y_COORDINATE_TRANSCRIPT_1,
			    SIDE_MARGIN + currGIR->exonNumber2 * SPACER + (currGIR->exonNumber2 - 1) * EXON_SIZE + EXON_SIZE / 2,
			    Y_COORDINATE_TRANSCRIPT_2,
			    black);
	}
	*/
	stringCreateClear (buffer,100);
	stringPrintf (buffer,"%s.jpg",currGE->id);
	jpegout = fopen(string (buffer), "wb");
	gdImageJpeg(im, jpegout, -1);
	fclose(jpegout);
	gdImageDestroy(im); 
}



int main (int argc, char *argv[])
{
	GfrEntry *currGE;
	int count;

	count = 0;
	gfr_init ("-");
	puts (gfr_writeHeader ());
	while (currGE = gfr_nextEntry ()) {
		makeImage (currGE);
		puts (gfr_writeGfrEntry (currGE));
		count++;
	}
	gfr_deInit ();
	warn ("%s_numGfrEntries: %d",argv[0],count);
	return 0;
}

