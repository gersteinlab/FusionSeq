#include <time.h>
#include <stdlib.h>


#include <bios/array.h>
#include <bios/log.h>
#include <bios/format.h>
#include <bios/numUtil.h>
#include <bios/html.h>
#include <bios/htmlLinker.h>
#include <bios/confp.h>
#include <mrf/mrf.h>
#include <mrf/segmentationUtil.h>

#include "sqvUtil.h"
#include "sqvWeb.h"
#include "sqvCircos.h"

#define RFILTER_THOLD 2
#define RFILTER_MAXGAP 5
#define RFILTER_MINRUN 20

static config *Conf = NULL;

/**
 * check if the start and end coordinate are valid
 */
void validateMappedRead( Locus *mr ) {
}

static int containedLocus ( MrfEntry* query, Locus* target )
{
  int overlap=0;
  
  if( query->isPairedEnd ) { // do this only if paired end
    int i;
    for( i=0; i<arrayMax( (query->read1).blocks ); i++) {
      MrfBlock* qBlock = arrp( (query->read1).blocks, i, MrfBlock);
      if( !strcmp( qBlock->targetName, target->chromosome) ) {
        if( positiveRangeIntersection( qBlock->targetStart, qBlock->targetEnd, target->start, target->end ) > 0 ) {
          overlap = 1; // found
          i = arrayMax( (query->read1).blocks ); // found, then stop
        }
      }
    }
    if( overlap == 0 ) { // not found in read 1, let's check read 2
      for( i=0; i<arrayMax( (query->read2).blocks ); i++) {
        MrfBlock* qBlock = arrp( (query->read2).blocks, i, MrfBlock);
        if( !strcmp( qBlock->targetName, target->chromosome) ) {
          if( positiveRangeIntersection( qBlock->targetStart, qBlock->targetEnd, target->start, target->end ) > 0 ) {
            overlap = 1; // found
            i = arrayMax( (query->read2).blocks ); // found, then stop
          }
        }
      }
    }
  }
  
  return overlap;
}

void printLocus( Locus l) {
  printf("%s:%d-%d\n", l.chromosome, l.start, l.end);
}
void sprintfLocus( char* str, Locus l) {
  sprintf(str, "%s\t%d\t%d", l.chromosome, l.start, l.end);
}

/** 
 * writes circos configuration file
 */
int write_circosConf (char* prefix,
		              Locus locus,
		              Array regions,
		              Chrdata_t *chromosomes,
		              SVCfg_t *settings)
{
  float rpos = 0.99;
  FILE *fp;
  Stringa buffer = stringCreate(50);
  int scale = getScale (regions);

  stringPrintf (buffer, "%s/test/circos_%s_%s_%d_%d.conf",
		confp_get(Conf, "WEB_CIRCOS_DIR"), 
		prefix, 
		locus.chromosome, 
		locus.start, 
		locus.end);

  if (!(fp = fopen (string (buffer), "w"))) {
    die ("Unable to open target file");
    return -1;
  };

  printf ("<h2>%i</h2>", scale);

  // write conf file
  conf_printHeader (fp, 
		    confp_get(Conf, "WEB_CIRCOS_DIR"), 
		    confp_get(Conf, "WEB_DATA_DIR"), 
		    confp_get(Conf, "WEB_SDATA_DIR"), 
		    chromosomes, 
		    prefix, 
		    locus, 
		    scale);
  conf_printUnits (fp, regions, chromosomes, scale);

  if (scale <= 10000) {
    conf_printDataTracks (fp, 
		    	  prefix, 
			  locus, 
			  confp_get(Conf, "WEB_SDATA_DIR"), 
			  confp_get(Conf, "WEB_DATA_DIR"), 
			  &rpos, 
			  regions, 
			  chromosomes, 
			  settings);
  }

  conf_printLinks (fp, confp_get(Conf, "WEB_DATA_DIR"), &rpos, prefix, locus, settings->readlim);
  conf_printFooter (fp);

  fclose (fp);
  stringDestroy (buffer);
  return 0;
}

static void randSelect (Array *arrPER, int readlim)
{
  Array new  = arrayCreate (readlim, PEreads);
  Array hash = arrayCreate (arrayMax (*arrPER), int);
  int rd;
  int max = arrayMax (*arrPER);
  int i, n;
  PEreads *newPER;

  // Initialize "hash table"
  for (i = 0; i < arrayMax(hash); i++)
    *(arrayp (hash, i, int)) = 0;

  n = 0;
  while (n < readlim) {
    rd = rand() % max;
    if (arru (hash, rd, int) == 0) {
      PEreads *oldPER = arrp (*arrPER, rd, PEreads);
      
      newPER = arrayp (new, n, PEreads);
      newPER->read1.chromosome = hlr_strdup (oldPER->read1.chromosome);
      newPER->read2.chromosome = hlr_strdup (oldPER->read2.chromosome);
      newPER->read1.start = oldPER->read1.start;
      newPER->read2.start = oldPER->read2.start;
      newPER->read1.end = oldPER->read1.end;
      newPER->read2.end = oldPER->read2.end;
      newPER->id = oldPER->id;

      *(arrayp (hash, rd, int)) = 1;
      n++;
    }
  }

  arrayDestroy (*arrPER);
  arrayDestroy (hash);
  *arrPER = new;
}

static void regionFilter (Array arrFiltered, Array arrPER, Array regions, int minspan)
{
  int i, j;
  int perspan;
  int chrnum1, chrnum2;
  int diffchr;

  PEreads *new;

  for (i = 0, j = 0; i < arrayMax (arrPER); i++) {
    PEreads* currPER = arrp (arrPER, i, PEreads);

    chrnum1 = getChrnum (currPER->read1.chromosome+3);
    chrnum2 = getChrnum (currPER->read2.chromosome+3);
    
    if (chrnum1 == chrnum2)
      diffchr = 0;
    else
      diffchr = 1;
    
    if (inRegions (regions, chrnum1, currPER->read1.start, currPER->read1.end) &&
        inRegions (regions, chrnum2, currPER->read2.start, currPER->read2.end)) {
      
      perspan = currPER->read2.start - currPER->read1.end;
      if (perspan < 0)
        perspan = 0 - perspan;
      
      if (diffchr || perspan >= minspan) {
        new = arrayp (arrFiltered, j, PEreads);
        new->read1.chromosome = hlr_strdup (currPER->read1.chromosome);
        new->read2.chromosome = hlr_strdup (currPER->read2.chromosome);
        new->read1.start = currPER->read1.start;
        new->read2.start = currPER->read2.start;
        new->read1.end = currPER->read1.end;
        new->read2.end = currPER->read2.end;
        new->id = currPER->id;
        j++;
      }
    }
  }
}

static void segmentFilter (char *prefix,
	                   Locus target,
	                   Array *arrFiltered,
	                   Array regions,
	                   int thold,
	                   int maxgap,
	                   int minrun,
	                   int showtars)
{
  Array newFiltered;

  Array arrWigs;
  Array arrTars;

  Array *tmpWigs;
  Array *tmpTars;

  SRegion_t *tmpReg;
  Wig *tmpW;
  Tar *tmpT;
  FILE *outTar;
  int i, j, k;
  int nregions = arrayMax (regions);
  int regionSpan;
  Stringa buffer = stringCreate (50);

  if (showtars) {
    stringPrintf (buffer, "%s/tmp/%s_%s_%d_%d.hlight.txt", 
		  confp_get(Conf, "WEB_DATA_DIR"), 
		  prefix, 
		  target.chromosome, 
		  target.start, 
		  target.end);
    if (!(outTar = fopen (string (buffer), "w"))) {
      die ("Unable to open target highlight file");
    }
  }

  tmpReg = arrayp (regions, 0, SRegion_t);
  if (tmpReg->chromosome == 0) {
    return;
  }

  //--------------------------------------------------------------------
  // Initialize arrays
  //--------------------------------------------------------------------
  newFiltered = arrayCreate (200, PEreads);
  arrWigs = arrayCreate (nregions, Array);
  arrTars = arrayCreate (nregions, Array);
  for (i = 0; i < nregions; i++) {
    tmpReg = arrayp (regions, i, SRegion_t);
    regionSpan = tmpReg->end - tmpReg->start + 1;

    tmpWigs = arrayp (arrWigs, i, Array);
    tmpTars = arrayp (arrTars, i, Array);

    *tmpWigs = arrayCreate (regionSpan, Wig);
    *tmpTars = arrayCreate (200, Tar);

    for (j = 0; j < regionSpan; j++) {
      tmpW = arrayp (*tmpWigs, j, Wig);
      tmpW->position = j + 1;
      tmpW->value    = 0;
    }
  }

  //--------------------------------------------------------------------
  // Make wigs
  //--------------------------------------------------------------------
  for (i = 0; i < arrayMax (*arrFiltered); i++) {
    PEreads* currPER = arrp (*arrFiltered, i, PEreads);

    // First read
    for (j = 0; j < nregions; j++) {
      tmpReg = arrayp (regions, j, SRegion_t);
      tmpWigs = arrayp (arrWigs, j, Array);

      if (tmpReg->chromosome == getChrnum (currPER->read1.chromosome+3) &&
          inRegion (currPER->read1.start, currPER->read1.end, tmpReg->start, tmpReg->end)) {
        break;
      }
    }

    int gstart = currPER->read1.start - tmpReg->start;
    int gend   = currPER->read1.end - tmpReg->start;
    for (k = gstart; k <= gend; k++) {
      tmpW = arrayp (*tmpWigs, k, Wig);
      tmpW->value += 1;
    }

    // Second read
    j = 0;
    for (j = 0; j < nregions; j++) {
      tmpReg = arrayp (regions, j, SRegion_t);
      tmpWigs = arrayp (arrWigs, j, Array);
      if (tmpReg->chromosome == getChrnum (currPER->read2.chromosome+3) &&
          inRegion (currPER->read2.start, currPER->read2.end, tmpReg->start, tmpReg->end)) {
        break;
      }
    }

	gstart = currPER->read2.start - tmpReg->start;
	gend   = currPER->read2.end - tmpReg->start;
	for (k = gstart; k <= gend; k++) {
	  tmpW = arrayp (*tmpWigs, k, Wig);
	  tmpW->value += 1;
	}
  }

  //--------------------------------------------------------------------
  // Perform segmentation, export tars to wigs
  //--------------------------------------------------------------------
  for (i = 0; i < nregions; i++) {
    tmpWigs = arrayp (arrWigs, i, Array);
    tmpTars = arrayp (arrTars, i, Array);
    tmpReg  = arrayp (regions, i, SRegion_t);

    Stringa chrname = stringCreate (3);
    if (tmpReg->chromosome == 23) {
      stringPrintf (chrname, "X");
    } else if (tmpReg->chromosome == 24) {
      stringPrintf (chrname, "Y");
    } else {
      stringPrintf (chrname, "%i", tmpReg->chromosome);
    }

    performSegmentation (*tmpTars, *tmpWigs, string (chrname), thold, maxgap, minrun);

    for (j = 0; j < arrayMax (*tmpWigs); j++) {
      tmpW = arrayp (*tmpWigs, j, Wig);
      tmpW->value = 0;
    }

    for (j = 0; j < arrayMax (*tmpTars); j++) {
      tmpT = arrayp (*tmpTars, j, Tar);
      int gstart = tmpT->start;
      int gend = tmpT->end;

      if (showtars)
        fprintf (outTar, "hs%s %i %i\n", string (chrname), gstart-1 + tmpReg->start, gend-1 + tmpReg->start);

      for (k = gstart-1; k <= gend-1; k++) {
        tmpW = arrayp (*tmpWigs, k, Wig);
        tmpW->value = 1;
      }
    }
    stringDestroy (chrname);
  }

  //--------------------------------------------------------------------
  // Filter reads
  //--------------------------------------------------------------------
  for (i = 0; i < arrayMax (*arrFiltered); i++) {
    PEreads* currPER = arrp (*arrFiltered, i, PEreads);
    int gstart, gend;

    // First read
    for (j = 0; j < nregions; j++) {
      tmpReg = arrayp (regions, j, SRegion_t);
      tmpWigs = arrayp (arrWigs, j, Array);
      if (tmpReg->chromosome == getChrnum (currPER->read1.chromosome+3) &&
          inRegion (currPER->read1.start, currPER->read1.end, tmpReg->start, tmpReg->end)) {
        break;
      }
    }

    gstart = currPER->read1.start - tmpReg->start;
    gend   = currPER->read1.end - tmpReg->start;
    if ((arrayp (*tmpWigs, gstart, Wig))->value > 0 ||
        (arrayp (*tmpWigs, gend, Wig))->value > 0) {
      for (j = 0; j < nregions; j++) {
        tmpReg = arrayp (regions, j, SRegion_t);
        tmpWigs = arrayp (arrWigs, j, Array);
        if (tmpReg->chromosome == getChrnum (currPER->read2.chromosome+3) &&
            inRegion (currPER->read2.start, currPER->read2.end, tmpReg->start, tmpReg->end)) {
          break;
        }
      }

      gstart = currPER->read2.start - tmpReg->start;
      gend   = currPER->read2.end - tmpReg->start;

      if ((arrayp (*tmpWigs, gstart, Wig))->value > 0 ||
    	  (arrayp (*tmpWigs, gend, Wig))->value > 0) {

        PEreads *newPER = arrayp (newFiltered, arrayMax (newFiltered), PEreads);
        newPER->id = currPER->id;
        newPER->read1.chromosome = hlr_strdup (currPER->read1.chromosome);
        newPER->read2.chromosome = hlr_strdup (currPER->read2.chromosome);
        newPER->read1.start = currPER->read1.start;
        newPER->read2.start = currPER->read2.start;
        newPER->read1.end = currPER->read1.end;
        newPER->read2.end = currPER->read2.end;
      }
    }
  }

  for (i = 0; i < nregions; i++) {
    tmpWigs = arrayp (arrWigs, i, Array);
    tmpTars = arrayp (arrTars, i, Array);
    arrayDestroy (*tmpWigs);
    arrayDestroy (*tmpTars);
  }

  arrayDestroy (arrWigs);
  arrayDestroy (arrTars);
  arrayDestroy (*arrFiltered);
  *arrFiltered = newFiltered;
  stringDestroy (buffer);
  if (showtars)
    fclose (outTar);
}

static void printLinks (FILE *fp, int *targetCounts, Array arrPER)
{
  int i, diffchr;
  for (i = 0; i < arrayMax (arrPER); i++) {
    PEreads* currPER = arrp (arrPER, i, PEreads);
    if (targetCounts[getChrnum (currPER->read1.chromosome+3) - 1] > 9 &&
        targetCounts[getChrnum (currPER->read2.chromosome+3) - 1] > 9 ) {
      Stringa chr1 = stringCreate (5);
      Stringa chr2 = stringCreate (5);

      Stringa entry1 = stringCreate (30);
      Stringa entry2 = stringCreate (30);

      stringPrintf (chr1, "hs%s", currPER->read1.chromosome+3);
      stringPrintf (chr2, "hs%s", currPER->read2.chromosome+3);

      diffchr = strDiffer (string (chr1), string (chr2));

      stringPrintf (entry1, "data%d %s %d %d",
	       currPER->id,
	       string(chr1),
	       currPER->read1.start,
	       currPER->read1.end);
      if (diffchr) {
        stringAppendf (entry1, " color=red");
      }

      stringPrintf (entry2, "data%d %s %d %d",
      	   currPER->id,
      	   string(chr2),
           currPER->read2.start,
           currPER->read2.end);
      if (diffchr) {
        stringAppendf (entry2, " color=red");
      }

      fprintf (fp, "%s\n", string (entry1));
      fprintf (fp, "%s\n", string (entry2));

      stringDestroy (entry1);
      stringDestroy (entry2);
      stringDestroy (chr1);
      stringDestroy (chr2);
    }
  }
}

static void generateOutput (char *prefix,
                            char *location,
                            Array regions,
                            Chrdata_t *chromosomes,
                            SVCfg_t *settings)
{ 
  int i;
  Stringa buffer = stringCreate (50);
  Stringa sbuffer = stringCreate (50);
  Stringa imgUrl = stringCreate (50);
  Stringa svgUrl = stringCreate (50);
  char *loc = strdup (location);

  WordIter w = wordIterCreate (location, ":-", 0);
  Locus target;
  target.chromosome = hlr_strdup (wordNext (w));

  char* tmp = wordNext (w);
  strTranslate (tmp, ",", "");
  target.start = atoi (tmp);

  tmp = wordNext (w);
  strTranslate (tmp, ",", "");
  target.end = atoi (tmp);

  if (target.end < target.start) {
    swapInt (&target.end, &target.start);
  }
  wordIterDestroy (w);
  stringPrintf (imgUrl, "%s/tmp/%s_%s_%d_%d.png", 
		confp_get(Conf, "WEB_DATA_LINK"), 
		prefix, target.chromosome, target.start, target.end );
  stringPrintf (svgUrl, "%s/tmp/%s_%s_%d_%d.svg", 
		confp_get(Conf, "WEB_DATA_LINK"), 
		prefix, target.chromosome, target.start, target.end );

  //--------------------------------------------------------
  // Display HTML header
  //--------------------------------------------------------
  web_printHTMLHead(regions, chromosomes, confp_get(Conf, "WEB_PUB_DIR"));
  web_printPageHeader(prefix, loc);

  web_printSidebar (confp_get(Conf, "WEB_URL_CGI"),
                    prefix,
                    loc,
                    string (imgUrl),
                    string (svgUrl),
                    regions,
                    chromosomes,
                    settings);
  fflush(stdout);


  //--------------------------------------------------------
  // Export links data from MRF
  //--------------------------------------------------------
  Stringa mrfFile = stringCreate (50);
  stringPrintf (mrfFile, "zcat %s/mrf/%s_%s.mrf.gz", confp_get(Conf, "WEB_DATA_DIR"), prefix, target.chromosome);
    //warn("here: %s", string(mrfFile));

  mrf_initFromPipe (string (mrfFile));

  MrfEntry* currMRFE;
  int count = 0;
  FILE *fp;
  FILE *sfp;
  Array arrPER = arrayCreate (1000, PEreads);
  Array arrFiltered = arrayCreate (1000, PEreads);

  int targetCounts[24];

  for (i = 0; i < 24; i++)
    targetCounts[i] = 0;

  while (currMRFE = mrf_nextEntry ()) {
	MrfBlock* currBlock1 = arrp ((currMRFE->read1).blocks, 0, MrfBlock);
	MrfBlock* currBlock2 = arrp ((currMRFE->read2).blocks, 0, MrfBlock);
    if (containedLocus (currMRFE, &target) &&
      strlen (currBlock1->targetName) < 6 &&
      strlen (currBlock2->targetName) < 6) {
      PEreads* currPER = arrayp (arrPER, arrayMax (arrPER), PEreads);
      if( !strcmp( currBlock1->targetName, currBlock2->targetName) ) { // same chromosome
        if( currBlock1->targetStart < currBlock2->targetStart ) {
          currPER->read1.chromosome = hlr_strdup (currBlock1->targetName);
          currPER->read2.chromosome = hlr_strdup (currBlock2->targetName);
          currPER->read1.start = currBlock1->targetStart;
          currPER->read2.start = currBlock2->targetStart;
          currPER->read1.end = currBlock1->targetEnd;
          currPER->read2.end = currBlock2->targetEnd;
        } else {
          currPER->read2.chromosome = hlr_strdup (currBlock1->targetName);
          currPER->read1.chromosome = hlr_strdup (currBlock2->targetName);
          currPER->read2.start = currBlock1->targetStart;
          currPER->read1.start = currBlock2->targetStart;
          currPER->read2.end = currBlock1->targetEnd;
          currPER->read1.end = currBlock2->targetEnd;
        }
      } else {
        int chr1, chr2;

        if( !strcmp(currBlock1->targetName,"X") ) {
          chr1=23;
        } else if(!strcmp(currBlock1->targetName,"Y") ) {
          chr1=24;
        } else {
          chr1=atoi( currBlock1->targetName+3 );
        }

        if( !strcmp(currBlock2->targetName,"X") ) {
          chr2=23;
        } else if(!strcmp(currBlock2->targetName,"Y") ) {
          chr2=24;
        } else {
          chr2=atoi( currBlock2->targetName+3 );
        }

        if( chr1 < chr2 ) {
          currPER->read1.chromosome = hlr_strdup (currBlock1->targetName);
          currPER->read2.chromosome = hlr_strdup (currBlock2->targetName);
          currPER->read1.start = currBlock1->targetStart;
          currPER->read2.start = currBlock2->targetStart;
          currPER->read1.end = currBlock1->targetEnd;
          currPER->read2.end = currBlock2->targetEnd;
        } else {
          currPER->read2.chromosome = hlr_strdup (currBlock1->targetName);
          currPER->read1.chromosome = hlr_strdup (currBlock2->targetName);
          currPER->read2.start = currBlock1->targetStart;
          currPER->read1.start = currBlock2->targetStart;
          currPER->read2.end = currBlock1->targetEnd;
          currPER->read1.end = currBlock2->targetEnd;
        }
      }
      currPER->id = count;
      targetCounts[getChrnum (currBlock1->targetName+3) - 1]++;
      targetCounts[getChrnum (currBlock2->targetName+3) - 1]++;
      count++;
    }
  }

  stringPrintf (buffer, "%s/tmp/%s_%s_%d_%d.data", 
		confp_get(Conf, "WEB_DATA_DIR"), 
		prefix, target.chromosome, target.start, target.end);
  if (!(fp = fopen (string (buffer), "w"))) {
    die ("Unable to open target file");
  };
  stringPrintf (sbuffer, "%s/tmp/%s_%s_%d_%d_s.data", 
		confp_get(Conf, "WEB_DATA_DIR"), 
		prefix, target.chromosome, target.start, target.end);
  if (!(sfp = fopen (string (sbuffer), "w"))) {
    die ("Unable to open target file");
  };

  arraySort (arrPER, (ARRAYORDERF) &sortPosition);
  arrayUniq (arrPER, NULL, (ARRAYORDERF) &sortPosition);
  
  regionFilter (arrFiltered, arrPER, regions, settings->minspan);
  if (settings->readlim != 0)
    randSelect(&arrFiltered, settings->readlim);

  if (settings->rfilter.on)
    segmentFilter (prefix,
    		       target,
    		       &arrFiltered,
    		       regions,
    		       settings->rfilter.thold,
    		       settings->rfilter.maxgap,
    		       settings->rfilter.minrun,
    		       settings->rfilter.showtars);

  printLinks (fp, targetCounts, arrPER);
  printLinks (sfp, targetCounts, arrFiltered);

  fclose (sfp);
  fclose (fp);

  //--------------------------------------------------------
  // Create Circos conf, run, make image
  //--------------------------------------------------------
  if (write_circosConf (prefix,
		                target,
		                regions,
		                chromosomes,
		                settings)) {
    die ("Error creating circos configuration file");
  }
  stringPrintf (buffer, "rm %s/tmp/%s_%s_%d_%d.png", 
		confp_get(Conf, "WEB_DATA_DIR"), 
		prefix, target.chromosome, target.start, target.end);
  //  hlr_system( string(buffer),0);


  //--------------------------------------------------------
  // Executing Circos
  //--------------------------------------------------------
  stringPrintf (buffer,"%s/bin/circos -conf %s/test/circos_%s_%s_%d_%d.conf -silent",
		confp_get(Conf, "WEB_CIRCOS_DIR"),
		confp_get(Conf, "WEB_CIRCOS_DIR"), 
		prefix, target.chromosome, target.start, target.end);

  hlr_system (string (buffer),0);
  hlr_free(target.chromosome);
  mrf_deInit ();

  printf ("<div class=\"img_circ\"><a href=\"%s\" target=\"_blank\">\n", string (imgUrl));
  printf ("<img width=\"650px\" src=\"%s\" id=\"circosImg\" border=\"1\" onload=\"setImageLoaded();\"/>", string (imgUrl));    
  puts ("</a></div>\n");

  //--------------------------------------------------------
  // Print rest of page
  //--------------------------------------------------------

  web_printBody (confp_get(Conf, "WEB_URL_CGI"),
                 prefix,
                 loc,
                 regions,
                 chromosomes,
                 settings);


  fflush (stdout);
  free (loc);
  arrayDestroy (arrPER);
  arrayDestroy (arrFiltered);
  stringDestroy (imgUrl);
  stringDestroy (svgUrl);
  stringDestroy (buffer);
  stringDestroy (sbuffer);
}


int main (int argc, char *argv[])
{
  if ((Conf = confp_open(getenv("FUSIONSEQ_CONFPATH"))) == NULL)
    return EXIT_FAILURE;

  cgiInit();
  cgiHeader("text/html");
  setenv ("PERL5LIB","/home/asboner/PerlModules/lib/perl5/site_perl/5.8.8/:/home/asboner/PerlModules/lib64/perl5/site_perl/5.8.8/:/home/asboner/PerlModules/lib/perl5/5.8.8/",1);  

  srand(time(NULL));
  char *querystr = cgiGet2Post(); // Look for params sent via GET

  Stringa item;
  Stringa value;
  char *iPtr, *vPtr;

  char *prefix = NULL;
  char *location = NULL;

  int chrnum, instnum, i;
  Array regions = arrayCreate(1, SRegion_t);
  Chrdata_t *chromosomes = (Chrdata_t *) calloc(25, sizeof(*chromosomes));
  SRegion_t *tmp;

  SVCfg_t settings;
  settings.minspan = 0;
  settings.readlim = 0;
  settings.rfilter.on = 0;
  settings.rfilter.thold = RFILTER_THOLD;
  settings.rfilter.maxgap = RFILTER_MAXGAP;
  settings.rfilter.minrun = RFILTER_MINRUN;
  settings.rfilter.showtars = 0;
  settings.tracks.expr = 0;
  settings.tracks.genes = 0;
  settings.tracks.exons = 0;

  int first = 1;

  tmp = arrayp(regions, 0, SRegion_t);
  tmp->chromosome = -1;
  tmp->instance   = -1;
  tmp->show       = 0;
  tmp->start      = 0;
  tmp->end        = -1;
  tmp->mstart     = 0;
  tmp->mend       = -1;

  for (i = 0; i < 25; i++) {
    chromosomes[i].show = 0;
    chromosomes[i].instances = 0;
  }

  if (strcmp(querystr, "") != 0) {
    item  = stringCreate (20);
	value = stringCreate (20);
	while (cgiGetNextPair (&first,item,value)) {
      iPtr = string (item);
      vPtr = string (value);
      if (strEqual (iPtr,"prefix")) {
        prefix = hlr_strdup ( vPtr );
      }
      if (strEqual (iPtr,"location")) {
        location = hlr_strdup ( vPtr );
      }

      if (strEqual (iPtr, "rfilter"))
        settings.rfilter.on = atoi (vPtr);
      if (strEqual (iPtr, "rthold"))
        settings.rfilter.thold = atoi (vPtr);
      if (strEqual (iPtr, "rminrun"))
        settings.rfilter.minrun = atoi (vPtr);
      if (strEqual (iPtr, "rmaxgap"))
        settings.rfilter.maxgap = atoi (vPtr);
      if (strEqual (iPtr, "rshowtars"))
        settings.rfilter.showtars = atoi (vPtr);
      if (strEqual (iPtr, "expr"))
        settings.tracks.expr = 1;
      if (strEqual (iPtr, "genes"))
        settings.tracks.genes = 1;
      if (strEqual (iPtr, "exons"))
        settings.tracks.exons = atoi (vPtr);
      if (strEqual (iPtr, "readlim"))
        settings.readlim = atoi (vPtr);
      if (strEqual (iPtr, "minspan"))
        settings.minspan = atoi (vPtr);

      if (!strncmp (iPtr, "hs", 2)) {
        if (!strcmp (&iPtr[2], "X")) {
          chrnum = 23;
         }
        else if (!strcmp (&iPtr[2], "Y")) {
          chrnum = 24;
         }
        else {
          chrnum = atoi (&iPtr[2]);
        }

        if (chromosomes[0].show != chrnum) {
          chromosomes[0].show = chrnum;
          chromosomes[0].instances++;
        }
        chromosomes[chrnum].show = 1;
      }

      //----
      // Note: this way of processing arguments is not "efficient" O(n^2).
      // However, most of the time, we are dealing with 1 region displayed.

      if (!strncmp(iPtr, "shs", 3) || !strncmp(iPtr, "mhs", 3)) {
        instnum = antoi(&iPtr[6], 6);
    	if (!strncmp(&iPtr[4], "X", 1))
    	  chrnum = 23;
    	else if (!strncmp(&iPtr[4], "Y", 1))
    	  chrnum = 24;
    	else
    	  chrnum  = antoi(&iPtr[3], 2);

    	for (i = 0; i < arrayMax(regions); i++) {
          tmp = arrayp(regions, i, SRegion_t);
          if ((chrnum == tmp->chromosome && instnum == tmp->instance) ||
        	  (chrnum == -1 && instnum && -1)) {
        	tmp->instance   = instnum;
        	tmp->chromosome = chrnum;
            break;
          }
    	}

    	if (i >= arrayMax(regions)) {
    	  tmp = arrayp(regions, i, SRegion_t);
    	  tmp->instance   = instnum;
    	  tmp->chromosome = chrnum;
    	  tmp->show       = 0;
    	  tmp->start      = 0;
    	  tmp->end        = -1;
    	  tmp->mend       = -1;
    	}

    	if (!strncmp(iPtr, "shs", 3)) {
    	  if (!strcmp(&iPtr[8], "show")) {
    	    tmp->show = 1;
    	   }
    	  else if (!strcmp(&iPtr[8], "start")) {
    	    tmp->start = atoi(vPtr);
    	   }
    	  else if (!strcmp(&iPtr[8], "end")) {
    	    tmp->end = atoi(vPtr);
    	  }
    	}
    	else {
    	  if (!strcmp(&iPtr[8], "start")) {
            tmp->mstart = atoi(vPtr);
    	   }
    	  else if (!strcmp(&iPtr[8], "end")) {
    	    tmp->mend = atoi(vPtr);
    	  }
    	}
      }
    }
	filterChrRegions (&regions, &chromosomes);
	if (!settings.rfilter.on) {
	  settings.rfilter.showtars = 0;
	}
  }

  if (!strcmp (querystr, "") || prefix == NULL || location == NULL) {
	web_printSingleLocusForm (confp_get(Conf, "WEB_URL_CGI"), confp_get(Conf, "WEB_PUB_DIR"));
    fflush (stdout);
   }
  else {
    generateOutput (prefix, location, regions, chromosomes, &settings);
  }

  cgiGet2PostReset();
  arrayDestroy(regions);
  free(chromosomes);

  confp_close(Conf);
  return 0;
}

