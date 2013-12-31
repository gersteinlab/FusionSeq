#include <bios/log.h>
#include <bios/format.h>
#include <bios/linestream.h>
#include <bios/common.h>
#include <bios/bits.h>

#include "gfr.h"



static LineStream lsGfr = NULL;
static Bits* presentColumnTypes = NULL;
static Array columnTypes = NULL;
static Texta columnHeaders = NULL;
static char* headerLine = NULL;



static void gfr_addColumnType (char *type)
{
  if (strEqual (type,GFR_COLUMN_NAME_NUM_INTER)) {
    bitSetOne (presentColumnTypes,GFR_COLUMN_TYPE_NUM_INTER);
    array (columnTypes,arrayMax (columnTypes),int) = GFR_COLUMN_TYPE_NUM_INTER;
    textAdd (columnHeaders,GFR_COLUMN_NAME_NUM_INTER);
  }
  else if (strEqual (type,GFR_COLUMN_NAME_INTER_MEAN_AB)) {
    bitSetOne (presentColumnTypes,GFR_COLUMN_TYPE_INTER_MEAN_AB);
    array (columnTypes,arrayMax (columnTypes),int) = GFR_COLUMN_TYPE_INTER_MEAN_AB;
    textAdd (columnHeaders,GFR_COLUMN_NAME_INTER_MEAN_AB);
  }
  else if (strEqual (type,GFR_COLUMN_NAME_INTER_MEAN_BA)) {
    bitSetOne (presentColumnTypes,GFR_COLUMN_TYPE_INTER_MEAN_BA);
    array (columnTypes,arrayMax (columnTypes),int) = GFR_COLUMN_TYPE_INTER_MEAN_BA;
    textAdd (columnHeaders,GFR_COLUMN_NAME_INTER_MEAN_BA);
  }
  else if (strEqual (type,GFR_COLUMN_NAME_PVALUE_AB)) {
    bitSetOne (presentColumnTypes,GFR_COLUMN_TYPE_PVALUE_AB);
    array (columnTypes,arrayMax (columnTypes),int) = GFR_COLUMN_TYPE_PVALUE_AB;
    textAdd (columnHeaders,GFR_COLUMN_NAME_PVALUE_AB);
  }
  else if (strEqual (type,GFR_COLUMN_NAME_PVALUE_BA)) {
    bitSetOne (presentColumnTypes,GFR_COLUMN_TYPE_PVALUE_BA);
    array (columnTypes,arrayMax (columnTypes),int) = GFR_COLUMN_TYPE_PVALUE_BA;
    textAdd (columnHeaders,GFR_COLUMN_NAME_PVALUE_BA);
  }
  else if (strEqual (type,GFR_COLUMN_NAME_NUM_INTRA1)) {
    bitSetOne (presentColumnTypes,GFR_COLUMN_TYPE_NUM_INTRA1);
    array (columnTypes,arrayMax (columnTypes),int) = GFR_COLUMN_TYPE_NUM_INTRA1;
    textAdd (columnHeaders,GFR_COLUMN_NAME_NUM_INTRA1);
  }
  else if (strEqual (type,GFR_COLUMN_NAME_NUM_INTRA2)) {
    bitSetOne (presentColumnTypes,GFR_COLUMN_TYPE_NUM_INTRA2);
    array (columnTypes,arrayMax (columnTypes),int) = GFR_COLUMN_TYPE_NUM_INTRA2;
    textAdd (columnHeaders,GFR_COLUMN_NAME_NUM_INTRA2);
  }
  else if (strEqual (type,GFR_COLUMN_NAME_FUSION_TYPE)) {
    bitSetOne (presentColumnTypes,GFR_COLUMN_TYPE_FUSION_TYPE);
    array (columnTypes,arrayMax (columnTypes),int) = GFR_COLUMN_TYPE_FUSION_TYPE;
    textAdd (columnHeaders,GFR_COLUMN_NAME_FUSION_TYPE);
  }
  else if (strEqual (type,GFR_COLUMN_NAME_NAME_TRANSCRIPT1)) {
    bitSetOne (presentColumnTypes,GFR_COLUMN_TYPE_NAME_TRANSCRIPT1);
    array (columnTypes,arrayMax (columnTypes),int) = GFR_COLUMN_TYPE_NAME_TRANSCRIPT1;
    textAdd (columnHeaders,GFR_COLUMN_NAME_NAME_TRANSCRIPT1);
  }
  else if (strEqual (type,GFR_COLUMN_NAME_CHROMOSOME_TRANSCRIPT1)) {
    bitSetOne (presentColumnTypes,GFR_COLUMN_TYPE_CHROMOSOME_TRANSCRIPT1);
    array (columnTypes,arrayMax (columnTypes),int) = GFR_COLUMN_TYPE_CHROMOSOME_TRANSCRIPT1;
    textAdd (columnHeaders,GFR_COLUMN_NAME_CHROMOSOME_TRANSCRIPT1);
  }
  else if (strEqual (type,GFR_COLUMN_NAME_STRAND_TRANSCRIPT1)) {
    bitSetOne (presentColumnTypes,GFR_COLUMN_TYPE_STRAND_TRANSCRIPT1);
    array (columnTypes,arrayMax (columnTypes),int) = GFR_COLUMN_TYPE_STRAND_TRANSCRIPT1;
    textAdd (columnHeaders,GFR_COLUMN_NAME_STRAND_TRANSCRIPT1);
  }
  else if (strEqual (type,GFR_COLUMN_NAME_NUM_EXONS_TRANSCRIPT1)) {
    bitSetOne (presentColumnTypes,GFR_COLUMN_TYPE_NUM_EXONS_TRANSCRIPT1);
    array (columnTypes,arrayMax (columnTypes),int) = GFR_COLUMN_TYPE_NUM_EXONS_TRANSCRIPT1;
    textAdd (columnHeaders,GFR_COLUMN_NAME_NUM_EXONS_TRANSCRIPT1);
  }
  else if (strEqual (type,GFR_COLUMN_NAME_EXON_COORDINATES_TRANSCRIPT1)) {
    bitSetOne (presentColumnTypes,GFR_COLUMN_TYPE_EXON_COORDINATES_TRANSCRIPT1);
    array (columnTypes,arrayMax (columnTypes),int) = GFR_COLUMN_TYPE_EXON_COORDINATES_TRANSCRIPT1;
    textAdd (columnHeaders,GFR_COLUMN_NAME_EXON_COORDINATES_TRANSCRIPT1);
  }
  else if (strEqual (type,GFR_COLUMN_NAME_START_TRANSCRIPT1)) {
    bitSetOne (presentColumnTypes,GFR_COLUMN_TYPE_START_TRANSCRIPT1);
    array (columnTypes,arrayMax (columnTypes),int) = GFR_COLUMN_TYPE_START_TRANSCRIPT1;
    textAdd (columnHeaders,GFR_COLUMN_NAME_START_TRANSCRIPT1);
  }
  else if (strEqual (type,GFR_COLUMN_NAME_END_TRANSCRIPT1)) {
    bitSetOne (presentColumnTypes,GFR_COLUMN_TYPE_END_TRANSCRIPT1);
    array (columnTypes,arrayMax (columnTypes),int) = GFR_COLUMN_TYPE_END_TRANSCRIPT1;
    textAdd (columnHeaders,GFR_COLUMN_NAME_END_TRANSCRIPT1);
  }
  else if (strEqual (type,GFR_COLUMN_NAME_GENE_SYMBOL_TRANSCRIPT1)) {
    bitSetOne (presentColumnTypes,GFR_COLUMN_TYPE_GENE_SYMBOL_TRANSCRIPT1);
    array (columnTypes,arrayMax (columnTypes),int) = GFR_COLUMN_TYPE_GENE_SYMBOL_TRANSCRIPT1;
    textAdd (columnHeaders,GFR_COLUMN_NAME_GENE_SYMBOL_TRANSCRIPT1);
  }
  else if (strEqual (type,GFR_COLUMN_NAME_DESCRIPTION_TRANSCRIPT1)) {
    bitSetOne (presentColumnTypes,GFR_COLUMN_TYPE_DESCRIPTION_TRANSCRIPT1);
    array (columnTypes,arrayMax (columnTypes),int) = GFR_COLUMN_TYPE_DESCRIPTION_TRANSCRIPT1;
    textAdd (columnHeaders,GFR_COLUMN_NAME_DESCRIPTION_TRANSCRIPT1);
  }
  else if (strEqual (type,GFR_COLUMN_NAME_NAME_TRANSCRIPT2)) {
    bitSetOne (presentColumnTypes,GFR_COLUMN_TYPE_NAME_TRANSCRIPT2);
    array (columnTypes,arrayMax (columnTypes),int) = GFR_COLUMN_TYPE_NAME_TRANSCRIPT2;
    textAdd (columnHeaders,GFR_COLUMN_NAME_NAME_TRANSCRIPT2);
  }
  else if (strEqual (type,GFR_COLUMN_NAME_CHROMOSOME_TRANSCRIPT2)) {
    bitSetOne (presentColumnTypes,GFR_COLUMN_TYPE_CHROMOSOME_TRANSCRIPT2);
    array (columnTypes,arrayMax (columnTypes),int) = GFR_COLUMN_TYPE_CHROMOSOME_TRANSCRIPT2;
    textAdd (columnHeaders,GFR_COLUMN_NAME_CHROMOSOME_TRANSCRIPT2);
  }
  else if (strEqual (type,GFR_COLUMN_NAME_STRAND_TRANSCRIPT2)) {
    bitSetOne (presentColumnTypes,GFR_COLUMN_TYPE_STRAND_TRANSCRIPT2);
    array (columnTypes,arrayMax (columnTypes),int) = GFR_COLUMN_TYPE_STRAND_TRANSCRIPT2;
    textAdd (columnHeaders,GFR_COLUMN_NAME_STRAND_TRANSCRIPT2);
  }
  else if (strEqual (type,GFR_COLUMN_NAME_NUM_EXONS_TRANSCRIPT2)) {
    bitSetOne (presentColumnTypes,GFR_COLUMN_TYPE_NUM_EXONS_TRANSCRIPT2);
    array (columnTypes,arrayMax (columnTypes),int) = GFR_COLUMN_TYPE_NUM_EXONS_TRANSCRIPT2;
    textAdd (columnHeaders,GFR_COLUMN_NAME_NUM_EXONS_TRANSCRIPT2);
  }
  else if (strEqual (type,GFR_COLUMN_NAME_EXON_COORDINATES_TRANSCRIPT2)) {
    bitSetOne (presentColumnTypes,GFR_COLUMN_TYPE_EXON_COORDINATES_TRANSCRIPT2);
    array (columnTypes,arrayMax (columnTypes),int) = GFR_COLUMN_TYPE_EXON_COORDINATES_TRANSCRIPT2;
    textAdd (columnHeaders,GFR_COLUMN_NAME_EXON_COORDINATES_TRANSCRIPT2);
  }
  else if (strEqual (type,GFR_COLUMN_NAME_START_TRANSCRIPT2)) {
    bitSetOne (presentColumnTypes,GFR_COLUMN_TYPE_START_TRANSCRIPT2);
    array (columnTypes,arrayMax (columnTypes),int) = GFR_COLUMN_TYPE_START_TRANSCRIPT2;
    textAdd (columnHeaders,GFR_COLUMN_NAME_START_TRANSCRIPT2);
  }
  else if (strEqual (type,GFR_COLUMN_NAME_END_TRANSCRIPT2)) {
    bitSetOne (presentColumnTypes,GFR_COLUMN_TYPE_END_TRANSCRIPT2);
    array (columnTypes,arrayMax (columnTypes),int) = GFR_COLUMN_TYPE_END_TRANSCRIPT2;
    textAdd (columnHeaders,GFR_COLUMN_NAME_END_TRANSCRIPT2);
  }
  else if (strEqual (type,GFR_COLUMN_NAME_GENE_SYMBOL_TRANSCRIPT2)) {
    bitSetOne (presentColumnTypes,GFR_COLUMN_TYPE_GENE_SYMBOL_TRANSCRIPT2);
    array (columnTypes,arrayMax (columnTypes),int) = GFR_COLUMN_TYPE_GENE_SYMBOL_TRANSCRIPT2;
    textAdd (columnHeaders,GFR_COLUMN_NAME_GENE_SYMBOL_TRANSCRIPT2);
  }
  else if (strEqual (type,GFR_COLUMN_NAME_DESCRIPTION_TRANSCRIPT2)) {
    bitSetOne (presentColumnTypes,GFR_COLUMN_TYPE_DESCRIPTION_TRANSCRIPT2);
    array (columnTypes,arrayMax (columnTypes),int) = GFR_COLUMN_TYPE_DESCRIPTION_TRANSCRIPT2;
    textAdd (columnHeaders,GFR_COLUMN_NAME_DESCRIPTION_TRANSCRIPT2);
  }
  else if (strEqual (type,GFR_COLUMN_NAME_INTER_READS)) {
    bitSetOne (presentColumnTypes,GFR_COLUMN_TYPE_INTER_READS);
    array (columnTypes,arrayMax (columnTypes),int) = GFR_COLUMN_TYPE_INTER_READS;
    textAdd (columnHeaders,GFR_COLUMN_NAME_INTER_READS);
  }
  else if (strEqual (type,GFR_COLUMN_NAME_PAIR_COUNT)) {
    bitSetOne (presentColumnTypes,GFR_COLUMN_TYPE_PAIR_COUNT);
    array (columnTypes,arrayMax (columnTypes),int) = GFR_COLUMN_TYPE_PAIR_COUNT;
    textAdd (columnHeaders,GFR_COLUMN_NAME_PAIR_COUNT);
  }
  else if (strEqual (type,GFR_COLUMN_NAME_ID)) {
    bitSetOne (presentColumnTypes,GFR_COLUMN_TYPE_ID);
    array (columnTypes,arrayMax (columnTypes),int) = GFR_COLUMN_TYPE_ID;
    textAdd (columnHeaders,GFR_COLUMN_NAME_ID);
  }
  else if (strEqual (type,GFR_COLUMN_NAME_READS_TRANSCRIPT1)) {
    bitSetOne (presentColumnTypes,GFR_COLUMN_TYPE_READS_TRANSCRIPT1);
    array (columnTypes,arrayMax (columnTypes),int) = GFR_COLUMN_TYPE_READS_TRANSCRIPT1;
    textAdd (columnHeaders,GFR_COLUMN_NAME_READS_TRANSCRIPT1);
  }
  else if (strEqual (type,GFR_COLUMN_NAME_READS_TRANSCRIPT2)) {
    bitSetOne (presentColumnTypes,GFR_COLUMN_TYPE_READS_TRANSCRIPT2);
    array (columnTypes,arrayMax (columnTypes),int) = GFR_COLUMN_TYPE_READS_TRANSCRIPT2;
    textAdd (columnHeaders,GFR_COLUMN_NAME_READS_TRANSCRIPT2);
  }
  else if (strEqual (type,GFR_COLUMN_NAME_SPER)) {
    bitSetOne (presentColumnTypes,GFR_COLUMN_TYPE_SPER);
    array (columnTypes,arrayMax (columnTypes),int) = GFR_COLUMN_TYPE_SPER;
    textAdd (columnHeaders,GFR_COLUMN_NAME_SPER);
  }  
  else if (strEqual (type,GFR_COLUMN_NAME_DASPER)) {
    bitSetOne (presentColumnTypes,GFR_COLUMN_TYPE_DASPER);
    array (columnTypes,arrayMax (columnTypes),int) = GFR_COLUMN_TYPE_DASPER;
    textAdd (columnHeaders,GFR_COLUMN_NAME_DASPER);
  }  
  else if (strEqual (type,GFR_COLUMN_NAME_RESPER)) {
    bitSetOne (presentColumnTypes,GFR_COLUMN_TYPE_RESPER);
    array (columnTypes,arrayMax (columnTypes),int) = GFR_COLUMN_TYPE_RESPER;
    textAdd (columnHeaders,GFR_COLUMN_NAME_RESPER);
  }
 else {
    die ("Unknown presentColumn: %s",type);
  }
}



int gfr_init (char *fileName) 
{
  int i;
  Texta tokens;
  lsGfr = ls_createFromFile (fileName);
  char* firstLine = ls_nextLine( lsGfr );
  if( firstLine==NULL) return 0;
  columnTypes = arrayCreate (20,int);
  columnHeaders = textCreate (20);
  presentColumnTypes = bitAlloc (100);  
  headerLine = hlr_strdup ( firstLine );
  tokens = textFieldtokP (headerLine,"\t");
  for (i = 0; i < arrayMax (tokens); i++) {
    gfr_addColumnType (textItem (tokens,i));
  }
  return 1;
}



void gfr_addNewColumnType (char* columnName)
{
  int i;

  i = 0;
  while (i < arrayMax (columnHeaders)) {
    if (strEqual (textItem (columnHeaders,i),columnName)) {
      break;
    } 
    i++;
  }
  if (i == arrayMax (columnHeaders)) {
    gfr_addColumnType (columnName);
  }
}



void gfr_deInit (void) 
{
  if (lsGfr != NULL) {
    ls_destroy (lsGfr);
  }
  arrayDestroy (columnTypes);
  textDestroy (columnHeaders);
  bitFree (&presentColumnTypes);
  hlr_free (headerLine);
}



static void gfr_freeEntry (GfrEntry* currEntry) 
{
  if (currEntry == NULL) {
    return;
  }
  if (bitReadOne (presentColumnTypes,GFR_COLUMN_TYPE_FUSION_TYPE)) {
    hlr_free (currEntry->fusionType);
  }
  if (bitReadOne (presentColumnTypes,GFR_COLUMN_TYPE_NAME_TRANSCRIPT1)) {
    hlr_free (currEntry->nameTranscript1);
  }
  if (bitReadOne (presentColumnTypes,GFR_COLUMN_TYPE_CHROMOSOME_TRANSCRIPT1)) {
    hlr_free (currEntry->chromosomeTranscript1);
  }
  if (bitReadOne (presentColumnTypes,GFR_COLUMN_TYPE_GENE_SYMBOL_TRANSCRIPT1)) {
    hlr_free (currEntry->geneSymbolTranscript1);
  }
  if (bitReadOne (presentColumnTypes,GFR_COLUMN_TYPE_EXON_COORDINATES_TRANSCRIPT1)) {
    arrayDestroy (currEntry->exonCoordinatesTranscript1);
  }
  if (bitReadOne (presentColumnTypes,GFR_COLUMN_TYPE_DESCRIPTION_TRANSCRIPT1)) {
    hlr_free (currEntry->descriptionTranscript1);
  }
  if (bitReadOne (presentColumnTypes,GFR_COLUMN_TYPE_NAME_TRANSCRIPT2)) {
    hlr_free (currEntry->nameTranscript2);
  }
  if (bitReadOne (presentColumnTypes,GFR_COLUMN_TYPE_CHROMOSOME_TRANSCRIPT2)) {
    hlr_free (currEntry->chromosomeTranscript2);
  }
  if (bitReadOne (presentColumnTypes,GFR_COLUMN_TYPE_GENE_SYMBOL_TRANSCRIPT2)) {
    hlr_free (currEntry->geneSymbolTranscript2);
  }
  if (bitReadOne (presentColumnTypes,GFR_COLUMN_TYPE_EXON_COORDINATES_TRANSCRIPT2)) {
    arrayDestroy (currEntry->exonCoordinatesTranscript2);
  }
  if (bitReadOne (presentColumnTypes,GFR_COLUMN_TYPE_DESCRIPTION_TRANSCRIPT2)) {
    hlr_free (currEntry->descriptionTranscript2);
  }
  if (bitReadOne (presentColumnTypes,GFR_COLUMN_TYPE_ID)) {
    hlr_free (currEntry->id);
  }
  if (bitReadOne (presentColumnTypes,GFR_COLUMN_TYPE_INTER_READS)) {
    arrayDestroy (currEntry->interReads);
  }
  if (bitReadOne (presentColumnTypes,GFR_COLUMN_TYPE_PAIR_COUNT)) {
    arrayDestroy (currEntry->pairCounts);
  }
  if (bitReadOne (presentColumnTypes,GFR_COLUMN_TYPE_READS_TRANSCRIPT1)) {
    textDestroy (currEntry->readsTranscript1);
  }
  if (bitReadOne (presentColumnTypes,GFR_COLUMN_TYPE_READS_TRANSCRIPT2)) {
    textDestroy (currEntry->readsTranscript2);
  }
  freeMem (currEntry);
  currEntry = NULL;
}



static GfrEntry* gfr_processNextEntry (int freeMemory) 
{
  static GfrEntry *currEntry = NULL;
  char *line,*token,*pos;
  WordIter w;
  int index,columnType;
  GfrPairCount *currGPC;
  GfrInterRead *currGIR;
  Texta tokens,items;
  int i;
  GfrExonCoordinate *currEC;
 
  if (!ls_isEof (lsGfr)) {
    while (line = ls_nextLine (lsGfr)) {
      if (line[0] == '\0') {
	continue;
      }
      if (freeMemory) {
        gfr_freeEntry (currEntry);
      }
      AllocVar (currEntry);
      index = 0;
      w = wordIterCreate (line,"\t",0);
      while (token = wordNext (w)) {
	columnType = arru (columnTypes,index,int);
	if (columnType == GFR_COLUMN_TYPE_NUM_INTER) {
	  currEntry->numInter = atoi (token);
	}
	else if (columnType == GFR_COLUMN_TYPE_INTER_MEAN_AB) {
	  currEntry->interMeanAB = atof (token);
	}
        else if (columnType == GFR_COLUMN_TYPE_INTER_MEAN_BA) {
	  currEntry->interMeanBA = atof (token);
	}
	else if (columnType == GFR_COLUMN_TYPE_PVALUE_AB) {
          currEntry->pValueAB = atof (token);
	}
        else if (columnType == GFR_COLUMN_TYPE_PVALUE_BA) {
          currEntry->pValueBA = atof (token);
	}
        else if (columnType == GFR_COLUMN_TYPE_NUM_INTRA1) {
          currEntry->numIntra1 = atoi (token);
	}
        else if (columnType == GFR_COLUMN_TYPE_NUM_INTRA2) {
          currEntry->numIntra2 = atoi (token);
	}
        else if (columnType == GFR_COLUMN_TYPE_FUSION_TYPE) {
          currEntry->fusionType = hlr_strdup (token);
	}
        else if (columnType == GFR_COLUMN_TYPE_NAME_TRANSCRIPT1) {
          currEntry->nameTranscript1 = hlr_strdup (token);
	}
        else if (columnType == GFR_COLUMN_TYPE_CHROMOSOME_TRANSCRIPT1) {
          currEntry->chromosomeTranscript1 = hlr_strdup (token);
	}
        else if (columnType == GFR_COLUMN_TYPE_STRAND_TRANSCRIPT1) {
          currEntry->strandTranscript1 = token[0];
	}
        else if (columnType == GFR_COLUMN_TYPE_NUM_EXONS_TRANSCRIPT1) {
          currEntry->numExonsTranscript1 = atoi (token);
	}
        else if (columnType == GFR_COLUMN_TYPE_EXON_COORDINATES_TRANSCRIPT1) {
          tokens = textFieldtok (token,"|");
          currEntry->exonCoordinatesTranscript1 = arrayCreate (100,GfrExonCoordinate);
          for (i = 0; i < arrayMax (tokens); i++) {
            currEC = arrayp (currEntry->exonCoordinatesTranscript1,arrayMax (currEntry->exonCoordinatesTranscript1),GfrExonCoordinate);
            pos = strchr (textItem (tokens,i),',');
            *pos = '\0';
            currEC->start = atoi (textItem (tokens,i));
            currEC->end = atoi (pos + 1);
          }
          textDestroy (tokens);
	}
        else if (columnType == GFR_COLUMN_TYPE_START_TRANSCRIPT1) {
          currEntry->startTranscript1 = atoi (token);
	}
        else if (columnType == GFR_COLUMN_TYPE_END_TRANSCRIPT1) {
	  currEntry->endTranscript1 = atoi (token);
	}
        else if (columnType == GFR_COLUMN_TYPE_GENE_SYMBOL_TRANSCRIPT1) {
          currEntry->geneSymbolTranscript1 = hlr_strdup (token);
	}
        else if (columnType == GFR_COLUMN_TYPE_DESCRIPTION_TRANSCRIPT1) {
          currEntry->descriptionTranscript1 = hlr_strdup (token);
	}
        else if (columnType == GFR_COLUMN_TYPE_NAME_TRANSCRIPT2) {
	  currEntry->nameTranscript2 = hlr_strdup (token);
	}
        else if (columnType == GFR_COLUMN_TYPE_CHROMOSOME_TRANSCRIPT2) {
          currEntry->chromosomeTranscript2 = hlr_strdup (token);
	}
        else if (columnType == GFR_COLUMN_TYPE_STRAND_TRANSCRIPT2) {
          currEntry->strandTranscript2 = token[0];
	}
        else if (columnType == GFR_COLUMN_TYPE_NUM_EXONS_TRANSCRIPT2) {
          currEntry->numExonsTranscript2 = atoi (token);
	}
        else if (columnType == GFR_COLUMN_TYPE_EXON_COORDINATES_TRANSCRIPT2) {
          tokens = textFieldtok (token,"|");
          currEntry->exonCoordinatesTranscript2 = arrayCreate (100,GfrExonCoordinate);
          for (i = 0; i < arrayMax (tokens); i++) {
            currEC = arrayp (currEntry->exonCoordinatesTranscript2,arrayMax (currEntry->exonCoordinatesTranscript2),GfrExonCoordinate);
            pos = strchr (textItem (tokens,i),',');
            *pos = '\0';
            currEC->start = atoi (textItem (tokens,i));
            currEC->end = atoi (pos + 1);
          }
          textDestroy (tokens);
	}
        else if (columnType == GFR_COLUMN_TYPE_START_TRANSCRIPT2) {
          currEntry->startTranscript2 = atoi (token);
	}
        else if (columnType == GFR_COLUMN_TYPE_END_TRANSCRIPT2) {
          currEntry->endTranscript2 = atoi (token);
	}
        else if (columnType == GFR_COLUMN_TYPE_GENE_SYMBOL_TRANSCRIPT2) {
          currEntry->geneSymbolTranscript2 = hlr_strdup (token);
	} 
        else if (columnType == GFR_COLUMN_TYPE_DESCRIPTION_TRANSCRIPT2) {
          currEntry->descriptionTranscript2 = hlr_strdup (token);
	}
        else if (columnType == GFR_COLUMN_TYPE_INTER_READS) {
          tokens = textFieldtok (token,"|");
          currEntry->interReads = arrayCreate (100,GfrInterRead);
          for (i = 0; i < arrayMax (tokens); i++) {
            if (textItem (tokens,i)[0] == '\0') {
              continue;
            }
            currGIR = arrayp (currEntry->interReads,arrayMax (currEntry->interReads),GfrInterRead);
            items = textFieldtok (textItem (tokens,i),",");	    
	    if( arrayMax( items ) > 6  ) {	      
	      currGIR->pairType = atoi (textItem (items,0));
              currGIR->number1 = atoi (textItem (items,1));
              currGIR->number2 = atoi (textItem (items,2));
	      currGIR->readStart1 = atoi (textItem (items,3));
	      currGIR->readEnd1 = atoi (textItem (items,4));
	      currGIR->readStart2 = atoi (textItem (items,5));
	      currGIR->readEnd2 = atoi (textItem (items,6));
	    } else {
	      currGIR->pairType = GFR_PAIR_TYPE_EXONIC_EXONIC;
	      currGIR->number1 = atoi (textItem (items,0));
	      currGIR->readStart1 = atoi (textItem (items,1));
	      currGIR->readEnd1 = atoi (textItem (items,2));
	      currGIR->number2 = atoi (textItem (items,3));
	      currGIR->readStart2 = atoi (textItem (items,4));
	      currGIR->readEnd2 = atoi (textItem (items,5));
	    }
	    currGIR->flag = 0;
            textDestroy (items);
          }
          textDestroy (tokens);
        }
  	else if (columnType == GFR_COLUMN_TYPE_PAIR_COUNT) {
          tokens = textFieldtok (token,"|");
          currEntry->pairCounts = arrayCreate (100,GfrPairCount);
          for (i = 0; i < arrayMax (tokens); i++) {
            currGPC = arrayp (currEntry->pairCounts,arrayMax (currEntry->pairCounts),GfrPairCount);
            items = textFieldtok (textItem (tokens,i),",");
	    if( arrayMax( items ) > 3  ) {
		  currGPC->pairType = atoi (textItem (items,0));
		  currGPC->count = atoi (textItem (items,1));
		  currGPC->number1 = atoi (textItem (items,2));
		  currGPC->number2 = atoi (textItem (items,3));
	    } else { 
	      currGPC->pairType = GFR_PAIR_TYPE_EXONIC_EXONIC;
	      currGPC->count = atoi (textItem (items,2));
	      currGPC->number1 = atoi (textItem (items,0));
	      currGPC->number2 = atoi (textItem (items,1));
	    }
            textDestroy (items);
	    }
          textDestroy (tokens);
        }
        else if (columnType == GFR_COLUMN_TYPE_ID) {
	  currEntry->id = hlr_strdup (token);
	}
        else if (columnType == GFR_COLUMN_TYPE_READS_TRANSCRIPT1) {
          tokens = textFieldtok (token,"|");
          currEntry->readsTranscript1 = textCreate (100);
          for (i = 0; i < arrayMax (tokens); i++) {
            textAdd (currEntry->readsTranscript1,textItem (tokens,i));
          }
          textDestroy (tokens);
        }
        else if (columnType == GFR_COLUMN_TYPE_READS_TRANSCRIPT2) {
          tokens = textFieldtok (token,"|");
          currEntry->readsTranscript2 = textCreate (100);
          for (i = 0; i < arrayMax (tokens); i++) {
            textAdd (currEntry->readsTranscript2,textItem (tokens,i));
          }
          textDestroy (tokens);
        }
       	else if (columnType == GFR_COLUMN_TYPE_SPER) {
	  currEntry->SPER = atof(token);
	}
	else if (columnType == GFR_COLUMN_TYPE_DASPER) {
	  currEntry->DASPER = atof(token);
	}
	else if (columnType == GFR_COLUMN_TYPE_RESPER) {
	  currEntry->RESPER = atof(token);
	}
        else {
	  die ("Unknown columnType: %d",columnType);
	}
	index++;
      }
      wordIterDestroy (w);
      return currEntry;
    }
  }
  if (freeMemory) {
    gfr_freeEntry (currEntry);
  }
  currEntry = NULL;
  return currEntry;
}



GfrEntry* gfr_nextEntry (void) 
{
  return gfr_processNextEntry (1); 
}



Array gfr_parse (void) 
{
  Array gfrEntries;
  GfrEntry *currEntry;

  gfrEntries = arrayCreate (100000,GfrEntry);
  while (currEntry = gfr_processNextEntry (0)) {
    array (gfrEntries,arrayMax (gfrEntries),GfrEntry) = *currEntry;
  }
  return gfrEntries;
}



static void gfr_addTab (Stringa buffer, int *first) 
{
  if (*first == 1) {
    *first = 0;
    return;
  }
  stringCatChar (buffer,'\t');
}



char* gfr_writeHeader (void)
{
  static Stringa buffer = NULL;
  int i;

  stringCreateClear (buffer,100);
  for (i = 0; i < arrayMax (columnHeaders); i++) {
    stringAppendf (buffer,"%s%s",textItem (columnHeaders,i), 
		   i < arrayMax (columnHeaders) - 1 ? "\t" : "");
  }
  return string (buffer);
}



char* gfr_writeGfrEntry (GfrEntry *currEntry)
{
  static Stringa buffer = NULL;
  int first;
  int i,j;
  int columnType;
  GfrPairCount *currGPC;
  GfrInterRead *currGIR;
  GfrExonCoordinate *currEC;

  stringCreateClear (buffer,100);
  first = 1;
  for (i = 0; i < arrayMax (columnTypes); i++) {
    columnType = arru (columnTypes,i,int);
    if (bitReadOne (presentColumnTypes,GFR_COLUMN_TYPE_NUM_INTER) && columnType == GFR_COLUMN_TYPE_NUM_INTER) {
      gfr_addTab (buffer,&first);
      stringAppendf (buffer,"%d",currEntry->numInter);
    }
    if (bitReadOne (presentColumnTypes,GFR_COLUMN_TYPE_INTER_MEAN_AB) && columnType == GFR_COLUMN_TYPE_INTER_MEAN_AB) {
      gfr_addTab (buffer,&first);
      stringAppendf (buffer,"%.2f",currEntry->interMeanAB);
    }
    if (bitReadOne (presentColumnTypes,GFR_COLUMN_TYPE_INTER_MEAN_BA) && columnType == GFR_COLUMN_TYPE_INTER_MEAN_BA) {
      gfr_addTab (buffer,&first);
      stringAppendf (buffer,"%.2f",currEntry->interMeanBA);
    }
    if (bitReadOne (presentColumnTypes,GFR_COLUMN_TYPE_PVALUE_AB) && columnType == GFR_COLUMN_TYPE_PVALUE_AB) {
      gfr_addTab (buffer,&first);
      stringAppendf (buffer,"%.5f",currEntry->pValueAB);
    }
    if (bitReadOne (presentColumnTypes,GFR_COLUMN_TYPE_PVALUE_BA) && columnType == GFR_COLUMN_TYPE_PVALUE_BA) {
      gfr_addTab (buffer,&first);
      stringAppendf (buffer,"%.5f",currEntry->pValueBA);
    }
    if (bitReadOne (presentColumnTypes,GFR_COLUMN_TYPE_NUM_INTRA1) && columnType == GFR_COLUMN_TYPE_NUM_INTRA1) {
      gfr_addTab (buffer,&first);
      stringAppendf (buffer,"%d",currEntry->numIntra1);
    }
    if (bitReadOne (presentColumnTypes,GFR_COLUMN_TYPE_NUM_INTRA2) && columnType == GFR_COLUMN_TYPE_NUM_INTRA2) {
      gfr_addTab (buffer,&first);
      stringAppendf (buffer,"%d",currEntry->numIntra2);
    }
    if (bitReadOne (presentColumnTypes,GFR_COLUMN_TYPE_FUSION_TYPE) && columnType == GFR_COLUMN_TYPE_FUSION_TYPE) {
      gfr_addTab (buffer,&first);
      stringAppendf (buffer,"%s",currEntry->fusionType);
    }
    if (bitReadOne (presentColumnTypes,GFR_COLUMN_TYPE_NAME_TRANSCRIPT1) && columnType == GFR_COLUMN_TYPE_NAME_TRANSCRIPT1) {
      gfr_addTab (buffer,&first);
      stringAppendf (buffer,"%s",currEntry->nameTranscript1);
    }
    if (bitReadOne (presentColumnTypes,GFR_COLUMN_TYPE_NUM_EXONS_TRANSCRIPT1) && columnType == GFR_COLUMN_TYPE_NUM_EXONS_TRANSCRIPT1) {
      gfr_addTab (buffer,&first);
      stringAppendf (buffer,"%d",currEntry->numExonsTranscript1);
    }
    if (bitReadOne (presentColumnTypes,GFR_COLUMN_TYPE_EXON_COORDINATES_TRANSCRIPT1) && columnType == GFR_COLUMN_TYPE_EXON_COORDINATES_TRANSCRIPT1) {
      gfr_addTab (buffer,&first);
      for (j = 0; j < arrayMax (currEntry->exonCoordinatesTranscript1); j++) {
        currEC = arrp (currEntry->exonCoordinatesTranscript1,j,GfrExonCoordinate);
        stringAppendf (buffer,"%d,%d%s",currEC->start,currEC->end,j < arrayMax (currEntry->exonCoordinatesTranscript1) - 1 ? "|" : "");
      }
    }
    if (bitReadOne (presentColumnTypes,GFR_COLUMN_TYPE_CHROMOSOME_TRANSCRIPT1) && columnType == GFR_COLUMN_TYPE_CHROMOSOME_TRANSCRIPT1) {
      gfr_addTab (buffer,&first);
      stringAppendf (buffer,"%s",currEntry->chromosomeTranscript1);
    }
    if (bitReadOne (presentColumnTypes,GFR_COLUMN_TYPE_STRAND_TRANSCRIPT1) && columnType == GFR_COLUMN_TYPE_STRAND_TRANSCRIPT1) {
      gfr_addTab (buffer,&first);
      stringAppendf (buffer,"%c",currEntry->strandTranscript1);
    }
    if (bitReadOne (presentColumnTypes,GFR_COLUMN_TYPE_START_TRANSCRIPT1) && columnType == GFR_COLUMN_TYPE_START_TRANSCRIPT1) {
      gfr_addTab (buffer,&first);
      stringAppendf (buffer,"%d",currEntry->startTranscript1);
    }
    if (bitReadOne (presentColumnTypes,GFR_COLUMN_TYPE_END_TRANSCRIPT1) && columnType == GFR_COLUMN_TYPE_END_TRANSCRIPT1) {
      gfr_addTab (buffer,&first);
      stringAppendf (buffer,"%d",currEntry->endTranscript1);
    }
    if (bitReadOne (presentColumnTypes,GFR_COLUMN_TYPE_GENE_SYMBOL_TRANSCRIPT1) && columnType == GFR_COLUMN_TYPE_GENE_SYMBOL_TRANSCRIPT1) {
      gfr_addTab (buffer,&first);
      stringAppendf (buffer,"%s",currEntry->geneSymbolTranscript1);
    }
    if (bitReadOne (presentColumnTypes,GFR_COLUMN_TYPE_DESCRIPTION_TRANSCRIPT1) && columnType == GFR_COLUMN_TYPE_DESCRIPTION_TRANSCRIPT1) {
      gfr_addTab (buffer,&first);
      stringAppendf (buffer,"%s",currEntry->descriptionTranscript1);
    }
    if (bitReadOne (presentColumnTypes,GFR_COLUMN_TYPE_NAME_TRANSCRIPT2) && columnType == GFR_COLUMN_TYPE_NAME_TRANSCRIPT2) {
      gfr_addTab (buffer,&first);
      stringAppendf (buffer,"%s",currEntry->nameTranscript2);
    }
    if (bitReadOne (presentColumnTypes,GFR_COLUMN_TYPE_NUM_EXONS_TRANSCRIPT2) && columnType == GFR_COLUMN_TYPE_NUM_EXONS_TRANSCRIPT2) {
      gfr_addTab (buffer,&first);
      stringAppendf (buffer,"%d",currEntry->numExonsTranscript2);
    }
    if (bitReadOne (presentColumnTypes,GFR_COLUMN_TYPE_EXON_COORDINATES_TRANSCRIPT2) && columnType == GFR_COLUMN_TYPE_EXON_COORDINATES_TRANSCRIPT2) {
      gfr_addTab (buffer,&first);
      for (j = 0; j < arrayMax (currEntry->exonCoordinatesTranscript2); j++) {
        currEC = arrp (currEntry->exonCoordinatesTranscript2,j,GfrExonCoordinate);
        stringAppendf (buffer,"%d,%d%s",currEC->start,currEC->end,j < arrayMax (currEntry->exonCoordinatesTranscript2) - 1 ? "|" : "");
      }
    }
    if (bitReadOne (presentColumnTypes,GFR_COLUMN_TYPE_CHROMOSOME_TRANSCRIPT2) && columnType == GFR_COLUMN_TYPE_CHROMOSOME_TRANSCRIPT2) {
      gfr_addTab (buffer,&first);
      stringAppendf (buffer,"%s",currEntry->chromosomeTranscript2);
    }
    if (bitReadOne (presentColumnTypes,GFR_COLUMN_TYPE_STRAND_TRANSCRIPT2) && columnType == GFR_COLUMN_TYPE_STRAND_TRANSCRIPT2) {
      gfr_addTab (buffer,&first);
      stringAppendf (buffer,"%c",currEntry->strandTranscript2);
    }
    if (bitReadOne (presentColumnTypes,GFR_COLUMN_TYPE_START_TRANSCRIPT2) && columnType == GFR_COLUMN_TYPE_START_TRANSCRIPT2) {
      gfr_addTab (buffer,&first);
      stringAppendf (buffer,"%d",currEntry->startTranscript2);
    }
    if (bitReadOne (presentColumnTypes,GFR_COLUMN_TYPE_END_TRANSCRIPT2) && columnType == GFR_COLUMN_TYPE_END_TRANSCRIPT2) {
      gfr_addTab (buffer,&first);
      stringAppendf (buffer,"%d",currEntry->endTranscript2);
    }
    if (bitReadOne (presentColumnTypes,GFR_COLUMN_TYPE_GENE_SYMBOL_TRANSCRIPT2) && columnType == GFR_COLUMN_TYPE_GENE_SYMBOL_TRANSCRIPT2) {
      gfr_addTab (buffer,&first);
      stringAppendf (buffer,"%s",currEntry->geneSymbolTranscript2);
    }
    if (bitReadOne (presentColumnTypes,GFR_COLUMN_TYPE_DESCRIPTION_TRANSCRIPT2) && columnType == GFR_COLUMN_TYPE_DESCRIPTION_TRANSCRIPT2) {
      gfr_addTab (buffer,&first);
      stringAppendf (buffer,"%s",currEntry->descriptionTranscript2);
    }
    if (bitReadOne (presentColumnTypes,GFR_COLUMN_TYPE_PAIR_COUNT) && columnType == GFR_COLUMN_TYPE_PAIR_COUNT) {
      gfr_addTab (buffer,&first);
      for (j = 0; j < arrayMax (currEntry->pairCounts); j++) {
        currGPC = arrp (currEntry->pairCounts,j,GfrPairCount);
        stringAppendf (buffer,"%d,%d,%d,%d%s",currGPC->pairType,currGPC->count,currGPC->number1,currGPC->number2,
                       j < arrayMax (currEntry->pairCounts) - 1 ? "|" : "");
      }
    }
    if (bitReadOne (presentColumnTypes,GFR_COLUMN_TYPE_INTER_READS) && columnType == GFR_COLUMN_TYPE_INTER_READS) {
      gfr_addTab (buffer,&first);
      for (j = 0; j < arrayMax (currEntry->interReads); j++) {
        currGIR = arrp (currEntry->interReads,j,GfrInterRead); 
        if (currGIR->flag == 0) {
          stringAppendf (buffer,"%d,%d,%d,%d,%d,%d,%d%s",
                         currGIR->pairType,currGIR->number1,currGIR->number2,
                         currGIR->readStart1,currGIR->readEnd1,
                         currGIR->readStart2,currGIR->readEnd2, 
                         j < arrayMax (currEntry->interReads) - 1 ? "|" : "");
        }
      }
    }
    if (bitReadOne (presentColumnTypes,GFR_COLUMN_TYPE_ID) && columnType == GFR_COLUMN_TYPE_ID) {
      gfr_addTab (buffer,&first);
      stringAppendf (buffer,"%s",currEntry->id);
    }
    if (bitReadOne (presentColumnTypes,GFR_COLUMN_TYPE_READS_TRANSCRIPT1) && columnType == GFR_COLUMN_TYPE_READS_TRANSCRIPT1) {
      gfr_addTab (buffer,&first);
      for (j = 0; j < arrayMax (currEntry->readsTranscript1); j++) {
        stringAppendf (buffer,"%s%s",textItem (currEntry->readsTranscript1,j),
                       j < arrayMax (currEntry->readsTranscript1) - 1 ? "|" : "");
      }
    }
    if (bitReadOne (presentColumnTypes,GFR_COLUMN_TYPE_READS_TRANSCRIPT2) && columnType == GFR_COLUMN_TYPE_READS_TRANSCRIPT2) {
      gfr_addTab (buffer,&first);
      for (j = 0; j < arrayMax (currEntry->readsTranscript2); j++) {
        stringAppendf (buffer,"%s%s",textItem (currEntry->readsTranscript2,j),
                       j < arrayMax (currEntry->readsTranscript2) - 1 ? "|" : "");
      }
    }
    if (bitReadOne (presentColumnTypes,GFR_COLUMN_TYPE_SPER) && columnType == GFR_COLUMN_TYPE_SPER) {
      gfr_addTab (buffer,&first);
      stringAppendf (buffer,"%f",currEntry->SPER);
    }
    if (bitReadOne (presentColumnTypes,GFR_COLUMN_TYPE_DASPER) && columnType == GFR_COLUMN_TYPE_DASPER) {
      gfr_addTab (buffer,&first);
      stringAppendf (buffer,"%f",currEntry->DASPER);
    }
    if (bitReadOne (presentColumnTypes,GFR_COLUMN_TYPE_RESPER) && columnType == GFR_COLUMN_TYPE_RESPER) {
      gfr_addTab (buffer,&first);
      stringAppendf (buffer,"%f",currEntry->RESPER);
    }
  }
  return string (buffer);
}

   
