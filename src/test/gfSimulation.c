// &&TAG*

#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <unistd.h>
#include <sys/types.h>

#include <bios/log.h>
#include <bios/format.h>
#include <bios/intervalFind.h>
#include <bios/linestream.h>
#include <bios/fasta.h>



#define TWO_BIT_FILE "/home1/lh372/DATA/blat/hg18.2bit"



typedef struct {
  char *name;
  double value;
} Gene;



typedef struct {
  Gene *gA;
  Interval *iA;
  Gene *gB;
  Interval *iB;
  int flag;
} Fusion;



static Array getValues (char *fileName) 
{
  Array genes;
  char *line;
  LineStream ls;
  Gene *currG;
  char *pos;

  genes = arrayCreate (10000,Gene); 
  ls = ls_createFromFile (fileName);
  while (line = ls_nextLine (ls)) {
    pos = strchr (line,'\t');
    *pos = '\0';
    currG = arrayp (genes,arrayMax (genes),Gene);
    currG->name = hlr_strdup (line);
    currG->value = atof (pos + 1);
  }
  ls_destroy (ls);
  return genes;
}



static int sortIntervals (Interval *a, Interval *b)
{
  return strcmp (a->name,b->name);
}



static int sortGenes (Gene *a, Gene *b)
{
  if (a->value < b->value) {
    return 1;
  }
  if (a->value > b->value) {
    return -1;
  }
  return 0;
}



static Array generateFusions (Array genes, int numFusionsToSimulate)
{
  Array fusions;
  Fusion *currF;
  Gene *a,*b;
  int i;

  fusions = arrayCreate (numFusionsToSimulate,Fusion);
  while (arrayMax (fusions) < numFusionsToSimulate) {
    a = arrp (genes,rand () % arrayMax (genes),Gene);
    b = arrp (genes,rand () % arrayMax (genes),Gene);
    if (a == b) {
      continue;
    }
    i = 0; 
    while (i < arrayMax (fusions)) {
      currF = arrp (fusions,i,Fusion);
      if ((currF->gA == a && currF->gB == b) || 
          (currF->gB == a && currF->gA == b)) {
        break;
      }
      i++;
    }
    if (i != arrayMax (fusions)) {
      continue;
    }
    currF = arrayp (fusions,arrayMax (fusions),Fusion);
    currF->gA = a;
    currF->gB = b;
  }
  return fusions;
}



static SubInterval* chooseExon (Interval *interval, int readLength)
{
  SubInterval *currSI;
  int index;
  
 while (1) {
    index = rand () % arrayMax (interval->subIntervals);
    currSI = arrp (interval->subIntervals,index,SubInterval);
    if ((currSI->end - currSI->start) > readLength) {
      break;
    }
  }
  return currSI;
}



static void getReads (Texta reads, Interval *currInterval, int start, int end, int readLength, int numReads, double fusionFactor) 
{
  int range;
  int i;
  int offset;
  static Stringa buffer = NULL;
  int totalNumReads;

  totalNumReads = numReads * fusionFactor;
  stringCreateClear (buffer,100);
  range = end - start - readLength;
  for (i = 0; i < totalNumReads; i++) {
    offset = rand () % range;
    stringPrintf (buffer,"%s:%c:%d:%d:%d:%d",currInterval->chromosome,currInterval->strand,start + offset,start + offset + readLength - 1,1,readLength);
    textAdd (reads,string (buffer));
  }
}



static void simulateReads (Array fusions, int readLength, double millionsOfMappedReads, double fusionFactor, FILE *fp1, FILE *fp2)
{
  Fusion *currF;
  int i,j;
  SubInterval *exonPtrA,*exonPtrB;
  int numReadsA,numReadsB;
  Texta readsA,readsB;
  int minNumReads;
  Stringa targetsFile;
  Stringa buffer;
  FILE *fp;
  Texta tokens;
  Array targetSeqs;
  int index;
  int numFusions;
 
  readsA = textCreate (1000);
  readsB = textCreate (1000);
  numFusions = 0;
  for (i = 0; i < arrayMax (fusions); i++) {
    currF = arrp (fusions,i,Fusion);
    if (currF->flag == 1) {
      continue;
    }
    exonPtrA = chooseExon (currF->iA,readLength);
    exonPtrB = chooseExon (currF->iB,readLength); 
    if (arrayMax (intervalFind_getOverlappingIntervals (currF->iA->chromosome,exonPtrA->start,exonPtrA->end)) > 1) {
      warn ("Non-unique exon: %s, %s:%d-%d",currF->iA->name,currF->iA->chromosome,exonPtrA->start,exonPtrA->end);
      continue;
    }
    if (arrayMax (intervalFind_getOverlappingIntervals (currF->iB->chromosome,exonPtrB->start,exonPtrB->end)) > 1) {
      warn ("Non-unique exon: %s, %s:%d-%d",currF->iB->name,currF->iB->chromosome,exonPtrB->start,exonPtrB->end);
      continue;
    }
    numFusions++;
    numReadsA = currF->gA->value * millionsOfMappedReads * (exonPtrA->end - exonPtrA->start + 1) / 1000;
    numReadsB = currF->gB->value * millionsOfMappedReads * (exonPtrB->end - exonPtrB->start + 1) / 1000;
    minNumReads = MIN (numReadsA,numReadsB);
    getReads (readsA,currF->iA,exonPtrA->start,exonPtrA->end,readLength,minNumReads,fusionFactor);
    getReads (readsB,currF->iB,exonPtrB->start,exonPtrB->end,readLength,minNumReads,fusionFactor);
    fprintf (fp2,"#####\n");
    fprintf (fp2,"%s, %f, %d\n",currF->gA->name,currF->gA->value,numReadsA);
    fprintf (fp2,"%s, %f, %d\n",currF->gB->name,currF->gB->value,numReadsB);
    fprintf (fp2,"Num reads simulated: %d \n",(int)(minNumReads * fusionFactor));
   }
  warn ("Done simulating reads...");
  buffer = stringCreate (100);
  targetsFile = stringCreate (100);
  stringPrintf (targetsFile,"targets_%d.txt",getpid ());
  if (!( fp = fopen (string (targetsFile),"w")) ){
    die ("Unable to open target file: %s",string (targetsFile));
  }
  if (arrayMax (readsA) != arrayMax (readsB)) {
    die ("Expected the same number of reads");
  }
  for (i = 0; i < arrayMax (readsA); i++) {
    tokens = textFieldtokP (textItem (readsA,i),":");
    fprintf (fp,"%s:%s-%s\n",textItem (tokens,0),textItem (tokens,2),textItem (tokens,3));
    textDestroy (tokens);
    tokens = textFieldtokP (textItem (readsB,i),":");
    fprintf (fp,"%s:%s-%s\n",textItem (tokens,0),textItem (tokens,2),textItem (tokens,3));
    textDestroy (tokens);
  }
  fclose (fp);
  stringPrintf (buffer,"twoBitToFa %s stdout -noMask -seqList=%s",TWO_BIT_FILE,string (targetsFile));
  fasta_initFromPipe (string (buffer));
  targetSeqs = fasta_readAllSequences (0);
  fasta_deInit ();
  warn ("Done obtaining read sequences...");
  stringPrintf (buffer,"rm -rf %s",string (targetsFile));
  hlr_system (string (buffer),0);
  fprintf (fp1,"AlignmentBlocks\tSequence\n");
  index = 0;
  for (j = 0; j < arrayMax (readsA); j++) {
    fprintf (fp1,"%s|%s\t%s|%s\n",textItem (readsA,j),textItem (readsB,j),arrp (targetSeqs,index,Seq)->sequence,arrp (targetSeqs,index + 1,Seq)->sequence);
    index = index + 2;
  }
  warn ("Done writing output...");
  textDestroy (readsA);
  textDestroy (readsB);
  stringDestroy (buffer);
  stringDestroy (targetsFile);
  warn ("Number of fusions simulated: %d",numFusions);
}



int main (int argc, char *argv[])
{
  Array intervals;
  Interval *currI;
  Array genes;
  int i,j;
  Gene *currG;
  double expressionCutoff,millionsOfMappedReads,fusionFactor;
  int readLength,numFusionsToSimulate;
  Array fusions;
  Fusion *currF;
  Stringa buffer;
  FILE *fp1,*fp2;

  if (argc != 9) {
    usage ("%s <file.interval> <file.expression> <expressionCutoff> <readLength> <millionsOfMappedReads> <fusionFactor> <numFusionsToSimulate> <outputPrefix>",argv[0]);
  }
  srand (time (0));
  intervalFind_addIntervalsToSearchSpace (argv[1],0);
  intervals = intervalFind_parseFile (argv[1],0);
  arraySort (intervals,(ARRAYORDERF)sortIntervals);
  genes = getValues (argv[2]); 
  arraySort (genes,(ARRAYORDERF)sortGenes);
  expressionCutoff = atof (argv[3]);
  readLength = atoi (argv[4]);
  millionsOfMappedReads = atof (argv[5]);
  fusionFactor = atof (argv[6]);
  numFusionsToSimulate = atoi (argv[7]);
  buffer = stringCreate (100);
  stringPrintf (buffer,"%s.mrf",argv[8]);
  fp1 = fopen (string (buffer),"w");
  stringPrintf (buffer,"%s.simMeta",argv[8]);
  fp2 = fopen (string (buffer),"w");
  if (fp1 == NULL || fp2 == NULL) {
    die ("Unable to open output files!");
  }
  fprintf (fp2,"expressionCutoff: %f\n",expressionCutoff);
  fprintf (fp2,"readLength: %d\n",readLength);
  fprintf (fp2,"millionsOfMappedReads: %f\n",millionsOfMappedReads);
  fprintf (fp2,"fusionFactor: %f\n",fusionFactor);
  i = 0; 
  while (i < arrayMax (genes)) {
    currG = arrp (genes,i,Gene);
    if (currG->value < expressionCutoff) {
      break;
    }
    i++;
  }
  arraySetMax (genes,(i - 1) >= 0 ? (i - 1) : 0);
  warn ("Number of genes above threshold: %d",arrayMax (genes));
  if (arrayMax (genes) < numFusionsToSimulate) {
    die ("Not enough gene candidates to simulate %d fusions. Lower expressionCutoff!",numFusionsToSimulate);
  }
  fusions = generateFusions (genes,numFusionsToSimulate);
  warn ("Done generating fusions...");
  for (i = 0; i < arrayMax (fusions); i++) {
    currF = arrp (fusions,i,Fusion);
    currF->iA = NULL;
    currF->iB = NULL;
    for (j = 0; j < arrayMax (intervals); j++) {
      currI = arrp (intervals,j,Interval);
      if (strEqual (currF->gA->name,currI->name)) {
        currF->iA = currI;
      }
      if (strEqual (currF->gB->name,currI->name)) {
        currF->iB = currI;
      }
    }
    if (currF->iA == NULL || currF->iB == NULL) {
      die ("Unable to find coordinates: %s %s",currF->gA->name,currF->gB->name);
    }
    if ((currF->iA->end - currF->iA->start) < readLength) {
      currF->flag = 1;
      warn ("Short gene: %s",currF->iA->name);
    }
    if ((currF->iB->end - currF->iB->start) < readLength) {
      currF->flag = 1;
      warn ("Short gene: %s",currF->iB->name);
    }
  }
  warn ("Starting read simulation...");
  simulateReads (fusions,readLength,millionsOfMappedReads,fusionFactor,fp1,fp2);
  fclose (fp1);
  fclose (fp2);
  return 0;
}



