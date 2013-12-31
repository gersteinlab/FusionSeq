#include <stdlib.h>

#include <bios/log.h>
#include <bios/format.h>
#include <bios/html.h>
#include <bios/htmlLinker.h>
#include <bios/confp.h>

#include "gfr.h"

#define GFR_PAIR_NAME_EXONIC_EXONIC "exon-exon"
#define GFR_PAIR_NAME_EXONIC_INTRONIC "exon-intron"
#define GFR_PAIR_NAME_EXONIC_JUNCTION "exon-boundary"
#define GFR_PAIR_NAME_INTRONIC_EXONIC "intron-exon"
#define GFR_PAIR_NAME_INTRONIC_INTRONIC "intron-intron"
#define GFR_PAIR_NAME_INTRONIC_JUNCTION "intron-boundary"
#define GFR_PAIR_NAME_JUNCTION_JUNCTION "boundary-boundary"
#define GFR_PAIR_NAME_JUNCTION_EXONIC "boundary-exon"
#define GFR_PAIR_NAME_JUNCTION_INTRONIC "boundary-intron"

static config *Conf = NULL;

static char* processString (char *str)
{
  static Stringa buffer = NULL;
  Texta tokens;
  int i;

  if (!strchr (str,'|')) {
    return str;
  }
  stringCreateClear (buffer,100);
  tokens = textStrtokP (str,"|");
  for (i = 0; i < arrayMax (tokens); i++) {
    stringAppendf (buffer,"%s%s",textItem (tokens,i),i < arrayMax (tokens) - 1 ? "<br>" : "");
  }
  return string (buffer);
}  

char* getPairTypeName(	int pairType ) {
  switch ( pairType ) {
  case GFR_PAIR_TYPE_EXONIC_EXONIC:
    return GFR_PAIR_NAME_EXONIC_EXONIC;
  case GFR_PAIR_TYPE_EXONIC_INTRONIC:
    return GFR_PAIR_NAME_EXONIC_INTRONIC;
  case GFR_PAIR_TYPE_EXONIC_JUNCTION:
    return GFR_PAIR_NAME_EXONIC_JUNCTION;
  case GFR_PAIR_TYPE_INTRONIC_EXONIC:
    return GFR_PAIR_NAME_INTRONIC_EXONIC;
  case GFR_PAIR_TYPE_INTRONIC_INTRONIC:
    return GFR_PAIR_NAME_INTRONIC_INTRONIC;
  case GFR_PAIR_TYPE_INTRONIC_JUNCTION:
    return GFR_PAIR_NAME_INTRONIC_JUNCTION;
  case GFR_PAIR_TYPE_JUNCTION_JUNCTION:
    return GFR_PAIR_NAME_JUNCTION_JUNCTION;
  case GFR_PAIR_TYPE_JUNCTION_EXONIC:
    return GFR_PAIR_NAME_JUNCTION_EXONIC;    
  case GFR_PAIR_TYPE_JUNCTION_INTRONIC:
    return GFR_PAIR_NAME_JUNCTION_INTRONIC;
  }
  return "NULL";
}
char* getEntryNumber( int number, int pairType, int readNum ) {
  int flag=0;
  Stringa str=stringCreate(10);
  switch ( pairType ) {
  case GFR_PAIR_TYPE_EXONIC_EXONIC:    
  case GFR_PAIR_TYPE_EXONIC_INTRONIC:
  case GFR_PAIR_TYPE_INTRONIC_EXONIC:   
  case GFR_PAIR_TYPE_INTRONIC_INTRONIC:
    break;
  case GFR_PAIR_TYPE_INTRONIC_JUNCTION:
  case GFR_PAIR_TYPE_EXONIC_JUNCTION:
    if( readNum==2 ) flag=1;      
    break;
  case GFR_PAIR_TYPE_JUNCTION_EXONIC:
  case GFR_PAIR_TYPE_JUNCTION_INTRONIC:
    if( readNum==1 ) flag=1;
    break;
  case GFR_PAIR_TYPE_JUNCTION_JUNCTION:
    flag=1;
  }
  if( flag==1 ) {
    int rem = number % 2;
    if( rem == 0 )
      stringPrintf(str, "%d right", number/2);
    else 
      stringPrintf(str, "left %d", number/2+1);
  } else {
    stringPrintf(str, "%d", number);
  }
  return string(str);
}

int main (int argc, char *argv[]) 
{
  FILE* ftmp = NULL;
  
  if ((Conf = confp_open(getenv("FUSIONSEQ_CONFPATH"))) == NULL)
    return EXIT_FAILURE;
  
  cgiInit();
  cgiHeader("text/html");

  if (argc == 3) {
    GfrEntry *currGE;
    Stringa buffer;
    GfrPairCount *currGEPC;
    GfrInterRead *currGIR;
    int i;

    puts ("<html>");
    puts ("<head>");
    html_printGenericStyleSheet (12);
    puts ("<title>geneFusions Details</title>\n");
    puts ("</head>");
    puts ("<body>");
    buffer = stringCreate (100);
    stringPrintf (buffer, "%s/%s.gfr", confp_get(Conf, "WEB_DATA_DIR"),argv[1]);    
    gfr_init (string (buffer));
    while (currGE = gfr_nextEntry ()){
      fflush( stdout );
      if (!strEqual (currGE->id,argv[2])) {
        continue;
      }
      printf ("<h1>Detailed summary for potential gene fusion candidate</h1><br>");
      puts ("<table border=0 cellpadding=10>");
      puts ("<tr align=left valign=top>");
      puts ("<td width=400>");
      puts ("<h2>Summary information</h2><br>");
      printf ("<b>Identifier</b>: %s<br><br>\n",currGE->id);
      printf ("<b>Number of inter paired-end reads</b>: %d<br><br>\n",currGE->numInter);
      printf ("<b>Type</b>: %s<br><br>\n",currGE->fusionType);     
      
      stringPrintf(buffer, "%s/GFF/%s.gff", confp_get(Conf, "WEB_DATA_DIR"),currGE->id);       
      ftmp = fopen( string(buffer), "r" ); // displaying this only if data are present
      if (ftmp) {
	 printf("<b>Connected Reads</b>: <a href=%s&hgt.customText=%s/GFF/%s.gff target=blank>UCSC connectivity graph</a><br>\n",
              	htmlLinker_generateLinkToGenomeBrowserAtUCSC ("hg18","vertebrate","human",
			currGE->chromosomeTranscript1,
			currGE->startTranscript1 - atoi(confp_get(Conf, "UCSC_GENOME_BROWSER_FLANKING_REGION")),
			currGE->endTranscript2 + atoi(confp_get(Conf, "UCSC_GENOME_BROWSER_FLANKING_REGION"))),
              confp_get(Conf, "WEB_DATA_LINK"),currGE->id); 
	 fclose( ftmp );
      }
   
      puts ("</td>");
      puts ("<td>");
      puts ("<h2>Transcript connectivity graph</h2>");
      printf ("<img src=%s/IMAGES/%s.jpg alt=geneFusionImage>\n", confp_get(Conf, "WEB_DATA_LINK"), currGE->id);
      puts ("</td>");
      puts ("<td>");
      puts ("<h2>Transcript connectivity table</h2><br>");
      puts ("<table border=0>");
      puts ("<tr align=left>");
      puts ("<th width=200>Pair Type</th>");
      puts ("<th width=200>Entry transcript 1</th>");
      puts ("<th width=200>Entry transcript 2</th>");
      puts ("<th width=200>Counts</th>");
      puts ("</tr>");
      fflush( stdout );
      for (i = 0; i < arrayMax (currGE->pairCounts); i++) {
        currGEPC = arrp (currGE->pairCounts,i,GfrPairCount);	
	printf ("<tr><td>%s</td><td>%s</td><td>%s</td><td>%d</td></tr>\n", 
		getPairTypeName(currGEPC->pairType), 
		getEntryNumber(currGEPC->number1, currGEPC->pairType, 1),
		getEntryNumber(currGEPC->number2, currGEPC->pairType, 2),
		currGEPC->count);
      }
      puts ("</table>");
      puts ("</td>");
      puts ("</tr>");
      puts ("</table>");
      puts ("<br>");

      puts ("<h2>Transcript information</h2><br>");
      puts ("<table border=1 cellpadding=10 width=\"80%\">");
      puts ("<tr align=left>");
      puts ("<th width=\"20%\"></th>");
      puts ("<th><font color='blue'>Transcript 1</font></th>");
      puts ("<th><font color='orange'>Transcript 2</font></th>");
      puts ("</tr>");
      puts ("<tr align=left>");
      puts ("<td width=\"20%\"><b>Gene symbol(s)</b></td>");
      printf ("<td width=\"30%%\"><font color='blue'>%s</font></td>\n",processString (currGE->geneSymbolTranscript1));
      printf ("<td width=\"30%%\"><font color='orange'>%s</font></td>\n",processString (currGE->geneSymbolTranscript2));
      puts ("</tr>");
      puts ("<tr align=left>");
      puts ("<td width=\"20%\"><b>Coordinates</b></td>");
      printf ("<td width=\"30%%\">%s:%d-%d</td>\n",currGE->chromosomeTranscript1,currGE->startTranscript1,currGE->endTranscript1);
      printf ("<td width=\"30%%\">%s:%d-%d</td>\n",currGE->chromosomeTranscript2,currGE->startTranscript2,currGE->endTranscript2);
      puts ("</tr>");
      puts ("<tr align=left>");
      puts ("<td width=\"20%\"><b>Strand</b></td>");
      printf ("<td width=\"30%%\">%c</td>\n",currGE->strandTranscript1);
      printf ("<td width=\"30%%\">%c</td>\n",currGE->strandTranscript2);
      puts ("</tr>");
      puts ("<tr align=left>");
      puts ("<td width=\"20%\"><b>Gene description(s)</b></td>");
      printf ("<td width=\"30%%\">%s</td>\n",processString (currGE->descriptionTranscript1));
      printf ("<td width=\"30%%\">%s</td>\n",processString (currGE->descriptionTranscript2));
      puts ("</tr>");
      puts ("<tr align=left>");
      puts ("<td width=\"20%\"><b>Number of exons</b></td>");
      printf ("<td width=\"30%%\">%d</td>\n",currGE->numExonsTranscript1);
      printf ("<td width=\"30%%\">%d</td>\n",currGE->numExonsTranscript2);
      puts ("</tr>");
      puts ("<tr align=left>");
      puts ("<td width=\"20%\"><b>Number of intra paired-end reads</b></td>");
      printf ("<td width=\"30%%\">%d</td>\n",currGE->numIntra1);
      printf ("<td width=\"30%%\">%d</td>\n",currGE->numIntra2);
      puts ("</tr>");
      puts ("<tr align=left>");
      puts ("<td width=\"20%\"><b>Links</b></td>");
      printf ("<td width=\"30%%\">[<a href=%s&hgt.customText=%s/BED/%s_1.bed target=blank>UCSC genome browser</a>]&nbsp;&nbsp;&nbsp;[<a href=%s/FASTA/%s_1.fasta>FASTA file</a>]<br></td>\n",
              htmlLinker_generateLinkToGenomeBrowserAtUCSC ("hg18","vertebrate","human",
		      currGE->chromosomeTranscript1,
		      currGE->startTranscript1 - atoi(confp_get(Conf, "UCSC_GENOME_BROWSER_FLANKING_REGION")),
		      currGE->endTranscript1 + atoi(confp_get(Conf, "UCSC_GENOME_BROWSER_FLANKING_REGION"))),
              confp_get(Conf, "WEB_DATA_LINK"),
	      currGE->id,
	      confp_get(Conf, "WEB_DATA_LINK"),
	      currGE->id); 
      printf ("<td width=\"30%%\">[<a href=%s&hgt.customText=%s/BED/%s_2.bed target=blank>UCSC genome browser</a>]&nbsp;&nbsp;&nbsp;[<a href=%s/FASTA/%s_2.fasta>FASTA file</a>]<br></td></tr>\n",
              htmlLinker_generateLinkToGenomeBrowserAtUCSC ("hg18","vertebrate","human",
		      currGE->chromosomeTranscript2,
		      currGE->startTranscript2 - atoi(confp_get(Conf, "UCSC_GENOME_BROWSER_FLANKING_REGION")),
		      currGE->endTranscript2 + atoi(confp_get(Conf, "UCSC_GENOME_BROWSER_FLANKING_REGION"))),
              confp_get(Conf, "WEB_DATA_LINK"),
	      currGE->id,
	      confp_get(Conf, "WEB_DATA_LINK"),
	      currGE->id); 
      
      puts ("<tr align=left>");
      puts ("<td width=\"20%\"><b>Expression</b></td>"); 

      stringPrintf(buffer, "%s/BGRS/%s_%s.bgr.gz", 
		   confp_get(Conf, "WEB_DATA_DIR"),
		   argv[1],
		   currGE->chromosomeTranscript1);  
      ftmp = fopen( string(buffer), "r" ); // displaying this only if data are present
      puts("<td width=\"30%\">");
      if( ftmp ) {
	printf ("[<a href=%s&hgt.customText=%s/BGRS/%s_%s.bgr.gz target=blank>Expression %s</a>]",
		htmlLinker_generateLinkToGenomeBrowserAtUCSC ("hg18","vertebrate","human",
			currGE->chromosomeTranscript1,
			currGE->startTranscript1 - atoi(confp_get(Conf, "UCSC_GENOME_BROWSER_FLANKING_REGION")),
			currGE->endTranscript1 + atoi(confp_get(Conf, "UCSC_GENOME_BROWSER_FLANKING_REGION"))),
		confp_get(Conf, "WEB_DATA_LINK"),
		argv[1],
		currGE->chromosomeTranscript1,
		currGE->chromosomeTranscript1); 
	fclose(ftmp);
      }
      puts("</td>");

      stringPrintf(buffer, "%s/BGRS/%s_%s.bgr.gz", confp_get(Conf, "WEB_DATA_DIR"),argv[1],currGE->chromosomeTranscript2); 
      ftmp = fopen( string(buffer), "r" ); // displaying this only if data are present
      puts("<td width=\"30%\">");
      if( ftmp ) {
	printf ("[<a href=%s&hgt.customText=%s/BGRS/%s_%s.bgr.gz target=blank>Expression %s</a>]",
		htmlLinker_generateLinkToGenomeBrowserAtUCSC ("hg18","vertebrate","human",
			currGE->chromosomeTranscript2,
			currGE->startTranscript2 - atoi(confp_get(Conf, "UCSC_GENOME_BROWSER_FLANKING_REGION")),
			currGE->endTranscript2 + atoi(confp_get(Conf, "UCSC_GENOME_BROWSER_FLANKING_REGION"))),
		confp_get(Conf, "WEB_DATA_LINK"),
		argv[1],
		currGE->chromosomeTranscript2,
		currGE->chromosomeTranscript2); 
	fclose(ftmp);
      } 
      puts("</td>");
      puts("</tr>");
      puts ("</table><br><br>");
      
      puts ("<h2>Breakpoint analysis</h2><br>");
      puts ("<table border=1 width=\"80%\" cellpadding=10><thead><tr><th>Orientation</th><th>Alignments</th><th colspan=2>Breakpoints</th></tr></thead><tbody>");
      puts ("<tr><td>Orientation AB</td>");
	if (currGE->strandTranscript1=='+') {
	  currGE->strandTranscript2=='+' ? stringPrintf(buffer, "AB_trans1F_trans2F") : stringPrintf(buffer, "AB_trans1F_trans2R");
	} else if( currGE->strandTranscript1 == '-') {
	  currGE->strandTranscript2=='+' ? stringPrintf(buffer, "AB_trans1R_trans2F") : stringPrintf(buffer, "AB_trans1R_trans2R");
	} else {
	  die("Strand informatation is not correct (transcript 1): %c", currGE->strandTranscript1);
	}
	printf ("<td align=center><a href=%s/ALIGNMENTS/%s_AB_breakPointAlignments.txt><img src=%s/IMAGES/%s.png></img>&nbsp;AB</a></td>", 
		confp_get(Conf, "WEB_DATA_LINK"), 
		currGE->id, 
		confp_get(Conf, "WEB_DATA_LINK"), 
		string(buffer)); 
	printf ("<td align=center><a href=%s&hgt.customText=%s/WIGS/%s_AB_breakPointsTranscript1.wig target=blank>Breakpoints transcript 1 UCSC Genome Browser</a></td>", 
		htmlLinker_generateLinkToGenomeBrowserAtUCSC ("hg18","vertebrate","human",
			currGE->chromosomeTranscript1,
			currGE->startTranscript1 - atoi(confp_get(Conf, "UCSC_GENOME_BROWSER_FLANKING_REGION")),
			currGE->endTranscript1 + atoi(confp_get(Conf, "UCSC_GENOME_BROWSER_FLANKING_REGION"))),
		confp_get(Conf, "WEB_DATA_LINK"),
		currGE->id);
	printf ("<td align=center><a href=%s&hgt.customText=%s/WIGS/%s_AB_breakPointsTranscript2.wig target=blank>Breakpoints transcript 2 UCSC Genome Browser</a></td></tr>", 
		htmlLinker_generateLinkToGenomeBrowserAtUCSC ("hg18","vertebrate","human",
			currGE->chromosomeTranscript2,
			currGE->startTranscript2 - atoi(confp_get(Conf, "UCSC_GENOME_BROWSER_FLANKING_REGION")),
			currGE->endTranscript2 + atoi(confp_get(Conf, "UCSC_GENOME_BROWSER_FLANKING_REGION"))),
		confp_get(Conf, "WEB_DATA_LINK"),
		currGE->id);	
      fflush(stdout);
      puts   ("<tr><td>Orientation BA</td>");  
      if (currGE->strandTranscript1 == '+') {
	currGE->strandTranscript2=='+' ? stringPrintf(buffer, "BA_trans1F_trans2F") : stringPrintf(buffer, "BA_trans1F_trans2R");
      } else if( currGE->strandTranscript1 == '-') {
	currGE->strandTranscript2=='+' ? stringPrintf(buffer, "BA_trans1R_trans2F") : stringPrintf(buffer, "BA_trans1R_trans2R");
      } else {
	die("Strand informatation is not correct (transcript2): %c", currGE->strandTranscript2);
	}	
      printf ("<td align=center><a href=%s/ALIGNMENTS/%s_BA_breakPointAlignments.txt><img src=%s/IMAGES/%s.png></img>&nbsp;BA</a></td>",
	      confp_get(Conf, "WEB_DATA_LINK"),
	      currGE->id, 
	      confp_get(Conf, "WEB_DATA_LINK"),
	      string(buffer));
      printf ("<td align=center><a href=%s&hgt.customText=%s/WIGS/%s_BA_breakPointsTranscript2.wig target=blank>Breakpoints transcript 2 UCSC Genome Browser</a></td>",	
	      htmlLinker_generateLinkToGenomeBrowserAtUCSC ("hg18","vertebrate","human",
		      currGE->chromosomeTranscript2,
		      currGE->startTranscript2 - atoi(confp_get(Conf, "UCSC_GENOME_BROWSER_FLANKING_REGION")),
	      	      currGE->endTranscript2 + atoi(confp_get(Conf, "UCSC_GENOME_BROWSER_FLANKING_REGION"))),
	      confp_get(Conf, "WEB_DATA_LINK"),
	      currGE->id);
      printf ("<td align=center><a href=%s&hgt.customText=%s/WIGS/%s_BA_breakPointsTranscript1.wig target=blank>Breakpoints transcript 1 UCSC Genome Browser</a></td></tr>", 
	      htmlLinker_generateLinkToGenomeBrowserAtUCSC ("hg18","vertebrate","human",
		      currGE->chromosomeTranscript1,
		      currGE->startTranscript1 - atoi(confp_get(Conf, "UCSC_GENOME_BROWSER_FLANKING_REGION")),
		      currGE->endTranscript1 + atoi(confp_get(Conf, "UCSC_GENOME_BROWSER_FLANKING_REGION"))),
	      confp_get(Conf, "WEB_DATA_LINK"),
	      currGE->id);       

      puts ("</tbody></table>");
      puts ("<br><br><br>");
      fflush(stdout);
    
    
      puts ("<h2>Read coordinates</h2><br>");
      puts ("<table border=0>");
      puts ("<tr align=left>");
      puts ("<th width=\"10%\">Pair Type</th>");
      puts ("<th width=\"10%\">Entry Transcript 1</th>");
      puts ("<th width=\"10%\">Read start transcript 1</th>");
      puts ("<th width=\"10%\">Read end transcript 1</th>");
      puts ("<th width=\"10%\">Entry Transcript 2</th>");
      puts ("<th width=\"10%\">Read start transcript 2</th>");
      puts ("<th width=\"10%\">Read end transcript 2</th>");
      puts ("</tr>");     
      for (i = 0; i < arrayMax (currGE->interReads); i++) {
	currGIR = arrp (currGE->interReads,i,GfrInterRead);
	printf ("<tr><td>%s</td><td>%s</td><td>%d</td><td>%d</td><td>%s</td><td>%d</td><td>%d</td></tr>\n",
		getPairTypeName(currGIR->pairType), 
		getEntryNumber(currGIR->number1, currGIR->pairType, 1),
		currGIR->readStart1,currGIR->readEnd1,
		getEntryNumber(currGIR->number2,currGIR->pairType, 2),
		currGIR->readStart2,
		currGIR->readEnd2);
      }
      puts ("</table><br><br><br>");
      puts ("</body>");
      puts ("</html>");
    fflush (stdout);
    }
  }
  confp_close(Conf);
  
  return EXIT_SUCCESS;
}
