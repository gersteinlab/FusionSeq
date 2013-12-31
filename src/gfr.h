#ifndef DEF_GFR_H
#define DEF_GFR_H



#define GFR_COLUMN_TYPE_NUM_INTER 1
#define GFR_COLUMN_TYPE_INTER_MEAN_AB 2
#define GFR_COLUMN_TYPE_INTER_MEAN_BA 3
#define GFR_COLUMN_TYPE_PVALUE_AB 4
#define GFR_COLUMN_TYPE_PVALUE_BA 5
#define GFR_COLUMN_TYPE_NUM_INTRA1 6
#define GFR_COLUMN_TYPE_NUM_INTRA2 7
#define GFR_COLUMN_TYPE_FUSION_TYPE 8
#define GFR_COLUMN_TYPE_NAME_TRANSCRIPT1 9
#define GFR_COLUMN_TYPE_NUM_EXONS_TRANSCRIPT1 10
#define GFR_COLUMN_TYPE_EXON_COORDINATES_TRANSCRIPT1 11
#define GFR_COLUMN_TYPE_CHROMOSOME_TRANSCRIPT1 12
#define GFR_COLUMN_TYPE_STRAND_TRANSCRIPT1 13
#define GFR_COLUMN_TYPE_START_TRANSCRIPT1 14
#define GFR_COLUMN_TYPE_END_TRANSCRIPT1 15
#define GFR_COLUMN_TYPE_NAME_TRANSCRIPT2 16
#define GFR_COLUMN_TYPE_NUM_EXONS_TRANSCRIPT2 17
#define GFR_COLUMN_TYPE_EXON_COORDINATES_TRANSCRIPT2 18
#define GFR_COLUMN_TYPE_CHROMOSOME_TRANSCRIPT2 19
#define GFR_COLUMN_TYPE_STRAND_TRANSCRIPT2 20
#define GFR_COLUMN_TYPE_START_TRANSCRIPT2 21
#define GFR_COLUMN_TYPE_END_TRANSCRIPT2 22
#define GFR_COLUMN_TYPE_PAIR_COUNT 23
#define GFR_COLUMN_TYPE_INTER_READS 24
#define GFR_COLUMN_TYPE_ID 25
#define GFR_COLUMN_TYPE_READS_TRANSCRIPT1 26
#define GFR_COLUMN_TYPE_READS_TRANSCRIPT2 27
#define GFR_COLUMN_TYPE_GENE_SYMBOL_TRANSCRIPT1 28
#define GFR_COLUMN_TYPE_GENE_SYMBOL_TRANSCRIPT2 29
#define GFR_COLUMN_TYPE_DESCRIPTION_TRANSCRIPT1 30
#define GFR_COLUMN_TYPE_DESCRIPTION_TRANSCRIPT2 31
#define GFR_COLUMN_TYPE_SPER 32
#define GFR_COLUMN_TYPE_DASPER 33
#define GFR_COLUMN_TYPE_RESPER 34


#define GFR_COLUMN_NAME_NUM_INTER "numInter"
#define GFR_COLUMN_NAME_INTER_MEAN_AB "interMeanAB"
#define GFR_COLUMN_NAME_INTER_MEAN_BA "interMeanBA"
#define GFR_COLUMN_NAME_PVALUE_AB "pValueAB"
#define GFR_COLUMN_NAME_PVALUE_BA "pValueBA"
#define GFR_COLUMN_NAME_NUM_INTRA1 "numIntra1"
#define GFR_COLUMN_NAME_NUM_INTRA2 "numIntra2"
#define GFR_COLUMN_NAME_FUSION_TYPE "fusionType"
#define GFR_COLUMN_NAME_NAME_TRANSCRIPT1 "nameTranscript1"
#define GFR_COLUMN_NAME_NUM_EXONS_TRANSCRIPT1 "numExonsTranscript1"
#define GFR_COLUMN_NAME_EXON_COORDINATES_TRANSCRIPT1 "exonCoordinatesTranscript1"
#define GFR_COLUMN_NAME_CHROMOSOME_TRANSCRIPT1 "chromosomeTranscript1"
#define GFR_COLUMN_NAME_STRAND_TRANSCRIPT1 "strandTranscript1"
#define GFR_COLUMN_NAME_START_TRANSCRIPT1 "startTranscript1"
#define GFR_COLUMN_NAME_END_TRANSCRIPT1 "endTranscript1"
#define GFR_COLUMN_NAME_NAME_TRANSCRIPT2 "nameTranscript2"
#define GFR_COLUMN_NAME_NUM_EXONS_TRANSCRIPT2 "numExonsTranscript2"
#define GFR_COLUMN_NAME_EXON_COORDINATES_TRANSCRIPT2 "exonCoordinatesTranscript2"
#define GFR_COLUMN_NAME_CHROMOSOME_TRANSCRIPT2 "chromosomeTranscript2"
#define GFR_COLUMN_NAME_STRAND_TRANSCRIPT2 "strandTranscript2"
#define GFR_COLUMN_NAME_START_TRANSCRIPT2 "startTranscript2"
#define GFR_COLUMN_NAME_END_TRANSCRIPT2 "endTranscript2"
#define GFR_COLUMN_NAME_PAIR_COUNT "pairCount"
#define GFR_COLUMN_NAME_INTER_READS "interReads"
#define GFR_COLUMN_NAME_ID "id"
#define GFR_COLUMN_NAME_READS_TRANSCRIPT1 "readsTranscript1"
#define GFR_COLUMN_NAME_READS_TRANSCRIPT2 "readsTranscript2"
#define GFR_COLUMN_NAME_GENE_SYMBOL_TRANSCRIPT1 "geneSymbolTranscript1"
#define GFR_COLUMN_NAME_GENE_SYMBOL_TRANSCRIPT2 "geneSymbolTranscript2"
#define GFR_COLUMN_NAME_DESCRIPTION_TRANSCRIPT1 "descriptionTranscript1"
#define GFR_COLUMN_NAME_DESCRIPTION_TRANSCRIPT2 "descriptionTranscript2"
#define GFR_COLUMN_NAME_SPER "SPER"
#define GFR_COLUMN_NAME_DASPER "DASPER"
#define GFR_COLUMN_NAME_RESPER "RESPER"



#define GFR_PAIR_TYPE_EXONIC_EXONIC 1
#define GFR_PAIR_TYPE_EXONIC_INTRONIC 2
#define GFR_PAIR_TYPE_EXONIC_JUNCTION 3
#define GFR_PAIR_TYPE_INTRONIC_EXONIC 4
#define GFR_PAIR_TYPE_INTRONIC_INTRONIC 5
#define GFR_PAIR_TYPE_INTRONIC_JUNCTION 6
#define GFR_PAIR_TYPE_JUNCTION_JUNCTION 7
#define GFR_PAIR_TYPE_JUNCTION_EXONIC 8
#define GFR_PAIR_TYPE_JUNCTION_INTRONIC 9



typedef struct {
  int pairType;
  int number1;
  int number2;
  int count;
} GfrPairCount;



typedef struct {
  int readStart1;
  int readStart2;
  int readEnd1;
  int readEnd2;
  int pairType;
  int number1;
  int number2;
  int flag;
} GfrInterRead;



typedef struct {
  int start;
  int end;
} GfrExonCoordinate;



typedef struct {
  int numInter;
  double interMeanAB;
  double interMeanBA;
  double pValueAB;
  double pValueBA;
  int numIntra1;
  int numIntra2;
  char *fusionType;
  char *nameTranscript1;
  char *chromosomeTranscript1;
  char strandTranscript1;
  int numExonsTranscript1;
  Array exonCoordinatesTranscript1;
  int startTranscript1;
  int endTranscript1;
  char *geneSymbolTranscript1;
  char *descriptionTranscript1;
  char *nameTranscript2;
  char *chromosomeTranscript2;
  char strandTranscript2;
  int numExonsTranscript2;
  Array exonCoordinatesTranscript2;
  int startTranscript2;
  int endTranscript2;
  char *geneSymbolTranscript2;
  char *descriptionTranscript2;
  Array interReads;
  Array pairCounts;
  char *id;
  Texta readsTranscript1;
  Texta readsTranscript2;
  double SPER;
  double DASPER;
  double RESPER;
} GfrEntry;



extern int gfr_init (char* fileName);
extern void gfr_addNewColumnType (char* columnName);
extern void gfr_deInit (void);
extern GfrEntry* gfr_nextEntry (void);
extern Array gfr_parse (void);
extern char* gfr_writeHeader (void);
extern char* gfr_writeGfrEntry (GfrEntry *currEntry);



#endif
