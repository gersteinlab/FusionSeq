#include <stdio.h>

#include <bios/log.h>
#include <bios/format.h>
#include <bios/linestream.h>
#include <bios/exportPEParser.h>
#include <mrf/mrf.h>

/**
 *  \file export2mrf.c Converts the export files of a paired-end experiment from Illumina to MRF.
 *  \author Andrea Sboner (andrea.sboner@yale.edu)
 */

static char convertStrand (char strand) 
{
  if (strand == 'F') {
    return '+';
  }
  if (strand != 'R') {
    die ("Unexpected strand: %c",strand);
  }
  return '-';
}


int main (int argc, char *argv[])
{
  int readLength;
  int chrEnd1, startEnd1;
  int chrEnd2, startEnd2;
  char *pos;
  Stringa id;
  FILE *fp;
  Stringa buffer, export1, export2;
  int countExport = 0;
  int countMRF = 0;
  int junction1, junction2;

  if (argc != 4) {
    usage ("%s <prefix> <exportFile1> <exportFile2>",argv[0]);
  }
  id = stringCreate (100);
  buffer = stringCreate (100);
  stringPrintf (buffer,"%s_allReads.txt",argv[1]);
  fp = fopen (string (buffer),"w");
  if (fp == NULL) {
    die ("Unable to open file: %s",string (buffer));
  }
  printf ("%s\t%s\t%s\n",MRF_COLUMN_NAME_BLOCKS,MRF_COLUMN_NAME_SEQUENCE,MRF_COLUMN_NAME_QUALITY_SCORES);

  export1 = stringCreate( 20 );
  export2 = stringCreate( 20 );
  if( strEndsWith( argv[2], ".gz" ) ) { // if export files are zipped, it uses zcat and createFromPipe
    stringPrintf( export1, "zcat %s", argv[2] );
    stringPrintf( export2, "zcat %s", argv[3] );
    exportPEParser_initFromPipe( string( export1 ), string( export2 ) );
  } else {
    exportPEParser_initFromFile( argv[2], argv[3] );
  }

  ExportPE* currExport;  
  while ( currExport = exportPEParser_nextEntry() ) {
    countExport++;
    junction1=0; 
    junction2=0;
    singleEnd* currEnd1 = currExport->end1;
    singleEnd* currEnd2 = currExport->end2;
    if(currEnd1->filter != 'Y' || currEnd2->filter != 'Y' ) {
      continue;
    }
    // skip reads with QC
    if (strStartsWith (currEnd1->chromosome,"QC") || strStartsWith (currEnd2->chromosome,"QC")  ) {
      continue;
    }
    // from here all reads should be kept
    fprintf (fp,"%s\n",currEnd1->sequence);
    fprintf (fp,"%s\n",currEnd2->sequence);

    // check for NM reads from exportFile1 or exportFile2
    if (strEqual (currEnd1->chromosome,"NM") || strEqual (currEnd2->chromosome,"NM" ) ) {
      continue;
    }
    // skip reads mapping to the mitochondrial chromosome
    if (strchr (currEnd1->chromosome,'M') || strchr (currEnd2->chromosome,'M')) {
      continue;
    }    
   if( strStartsWith( currEnd1->chromosome, "ribosomal") || strStartsWith( currEnd2->chromosome, "ribosomal")  ) {
      continue;
    }
    if( strchr( currEnd1->chromosome, ':' ) || strchr( currEnd2->chromosome, ':' ) ) { // multi-hits
      continue;
    }
    static Stringa end1 = NULL;
    static Stringa end2 = NULL;
    stringCreateClear( end1, 30);
    stringCreateClear( end2, 30);
   

    readLength = strlen( currEnd1->sequence );
    if( strStartsWith( currEnd1->chromosome, "chr") ) { // genomic
      pos = strchr ( currEnd1->chromosome,'.');	  
      *pos = '\0';
      chrEnd1 = atoi( currEnd1->chromosome+3 ); 
      startEnd1 = currEnd1->position;
      if( currEnd1->position<0 ) 
	continue;
      stringPrintf( end1, "%s:%c:%d:%d:1:%d", currEnd1->chromosome, convertStrand (currEnd1->strand), currEnd1->position, currEnd1->position + readLength - 1, readLength );
    } else if( strstr( currEnd1->chromosome, "junction" ) ) { // junction
      Texta location = textFieldtokP ( currEnd1->contig, "|" );
      int tilestart1 = atoi( textItem( location, 1) );
      int tilestart2 = atoi( textItem( location, 2) );
      int tileSize = atoi( textItem( location, 3 ) );
      if( (tileSize - currEnd1->position + 1 ) >= readLength  || ( (tileSize - currEnd1->position) < 0) ) 
	continue;
      int targetStart1 = tilestart1 + currEnd1->position; // tileStart1 is zero-based, thus no need to subtract 1
      int targetEnd1   = tilestart1 + tileSize; // similarly, no need to subtract 1
      int queryStart1  = 1;
      int queryEnd1    = (targetEnd1 - targetStart1 + 1); // targetEnd1 and targetStart1 are in the same space, thus need to add 1
      int basesOnSecondBlock = readLength - (queryEnd1 - queryStart1 + 1);

      int targetStart2 = tilestart2 + 1 ; // tileStart2 is zero-based
      int targetEnd2   = targetStart2 + basesOnSecondBlock - 1; 
      int queryStart2  = queryEnd1 + 1;
      int queryEnd2    = readLength;
      stringPrintf( end1, "%s:%c:%d:%d:%d:%d,%s:%c:%d:%d:%d:%d", 
		    textItem( location, 0),convertStrand (currEnd1->strand),targetStart1,targetEnd1, queryStart1, queryEnd1, // first block
		    textItem( location, 0),convertStrand (currEnd1->strand),targetStart2,targetEnd2, queryStart2, queryEnd2  // second block
		    );
      chrEnd1 = atoi( textItem( location, 0)+3  ); 
      startEnd1 = targetStart1;
      textDestroy( location );
      junction1 = 1;
    } else {
      warn( "End 1 not in a proper format:%s", exportPEParser_writeEntry( currEnd1 ) );
      continue;
    }
    readLength = strlen( currEnd2->sequence );
    if( strStartsWith( currEnd2->chromosome, "chr") ) { // genomic - inter
      pos = strchr ( currEnd2->chromosome,'.');	  
      *pos = '\0';
      chrEnd2 = atoi( currEnd2->chromosome+3 );
      startEnd2 = currEnd2->position;
      stringPrintf( end2, "%s:%c:%d:%d:1:%d", currEnd2->chromosome, convertStrand ( currEnd2->strand), currEnd2->position, currEnd2->position + readLength - 1, readLength );      
    } else if( strstr( currEnd2->chromosome, "junction" ) ) { // junction
      Texta location = textFieldtokP ( currEnd2->contig, "|" );
      //strReplace( &partnerChromosome,  textItem( location, 0));
      int tilestart1 = atoi( textItem( location, 1) );
      int tilestart2 = atoi( textItem( location, 2) );
      int tileSize = atoi( textItem( location, 3 ) );
      if( (tileSize - currEnd2->position + 1) >= readLength || ( (tileSize - currEnd2->position ) < 0) )
	continue;
      int targetStart1 = tilestart1 + currEnd2->position;
      int targetEnd1   = tilestart1 + tileSize;
      int queryStart1  = 1;
      int queryEnd1    = (targetEnd1 - targetStart1 + 1);
      int basesOnSecondBlock = readLength - (queryEnd1 - queryStart1 + 1);

      int targetStart2 = tilestart2 + 1 ; // tileStart2 is zero-based
      int targetEnd2   = targetStart2 + basesOnSecondBlock - 1;
      int queryStart2  = queryEnd1 + 1;
      int queryEnd2    = readLength;
      stringPrintf( end2, "%s:%c:%d:%d:%d:%d,%s:%c:%d:%d:%d:%d", 
		    textItem( location, 0),convertStrand (currEnd2->strand),targetStart1, targetEnd1, queryStart1, queryEnd1, // first block
		    textItem( location, 0),convertStrand (currEnd2->strand),targetStart2, targetEnd2, queryStart2, queryEnd2  // second block
		    );
      chrEnd2=atoi(textItem( location, 0)+3 );
      startEnd2 = targetStart2;
      textDestroy( location );
      junction2=1;
    } else {
      warn( "End 2 not in a proper format:%s", exportPEParser_writeEntry( currEnd2 ) );
      continue;
    }
    
    if( chrEnd1 < chrEnd2   ) {
      printf ("%s|%s\t%s|%s\t%s|%s\n", string(end1), string(end2), currEnd1->sequence, currEnd2->sequence, currEnd1->quality, currEnd2->quality);
    } else {
      if( chrEnd1 > chrEnd2 ) {
	printf ("%s|%s\t%s|%s\t%s|%s\n", string(end2), string(end1), currEnd2->sequence, currEnd1->sequence, currEnd2->quality, currEnd1->quality);
      } else { // same chromosome
	/*if( abs( startEnd2 - startEnd1 ) < readLength ) {
	  warn("Overlapping reads (1):\t%s", exportPEParser_writeEntry( currEnd1 ) );
	  warn("Overlapping reads (2):\t%s", exportPEParser_writeEntry( currEnd2 ) );
	  continue;
	  }*/
	if ((currEnd1->strand == 'R' && currEnd2->strand == 'F') || (currEnd1->strand == 'F' && currEnd2->strand == 'R') && 
	    ( junction1 == 0 && junction2 == 0) ) {
	  if( (startEnd2-startEnd1) > 0 ) {
	    printf ("%s|%s\t%s|%s\t%s|%s\n", string(end1), string(end2), currEnd1->sequence, currEnd2->sequence, currEnd1->quality, currEnd2->quality);
	  } else {
	    printf ("%s|%s\t%s|%s\t%s|%s\n", string(end2), string(end1), currEnd2->sequence, currEnd1->sequence, currEnd2->quality, currEnd1->quality);
	  }
	} else {
	  if( junction1 == 1 || junction2 == 1 ) {
	    if( (startEnd2-startEnd1) > 0 ) {
	      printf ("%s|%s\t%s|%s\t%s|%s\n", string(end1), string(end2), currEnd1->sequence, currEnd2->sequence, currEnd1->quality, currEnd2->quality);
	    } else {
	      printf ("%s|%s\t%s|%s\t%s|%s\n", string(end2), string(end1), currEnd2->sequence, currEnd1->sequence, currEnd2->quality, currEnd1->quality);
	    }
	  } else {
	    warn( "Same strand (1):\t%s", exportPEParser_writeEntry ( currEnd1 ) );
	    warn( "Same strand (2):\t%s", exportPEParser_writeEntry ( currEnd2 ) );
	    continue;
	  }
	}    	      	  
      }	
    }
    stringDestroy( end1 );
    stringDestroy( end2 );
    countMRF++;
  }
  fclose (fp);
  
  stringPrintf (buffer,"%s.meta",argv[1]);
  fp = fopen (string (buffer),"w");
  if (fp == NULL) {
    die ("Unable to open file: %s",string (buffer));
  } 
  fprintf(fp, "Total_number_reads\t%d\nMapped_reads\t%d", countExport,countMRF);
  fclose (fp);
  stringDestroy (buffer);
  stringDestroy (id);
  return 0;
}



