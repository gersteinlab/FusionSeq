#include <stdlib.h>

#include <bios/confp.h>
#include <bios/log.h>
#include <bios/format.h>
#include <bios/html.h>
#include <bios/htmlLinker.h>
#include <bios/linestream.h>
#include <bios/stringUtil.h>

#include "gfr.h"


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



static void generateOutput (char* prefix, char* typeSelected, int minNum)
{
  GfrEntry *currGE;
  Stringa buffer;
  char *pos;

  puts ("<html>");
  puts ("<head>");
  puts ("<title>Results - Gene Fusions</title>");
  html_printGenericStyleSheet (12);
  puts ("</head>");
  puts ("<body>");
  if (prefix[0] == '\0') {
    die ("Invalid prefix");
  }
  printf ("<h1>Results - %s</h1><br><br><br>",prefix);

  buffer = stringCreate(50);
  //Chromosome expression, if present
  LineStream ls;
  char* chrSignal=NULL;  
  stringPrintf(buffer, "ls -1 %s/BGRS/%s_chr*.bgr.gz 2> /dev/null", 
	       confp_get(Conf, "WEB_DATA_DIR"), 
	       prefix);
  ls = ls_createFromPipe(string(buffer));
  int countCol = 0;
  puts ("Expression signal: &nbsp;");
  fflush(stdout);
  while( chrSignal = ls_nextLine(ls)) {
        
	char* chrTmp = stringBetween( prefix, ".bgr.gz", chrSignal );
	chrTmp++;      
	printf ("[<a href=%s&hgt.customText=%s/BGRS/%s_%s.bgr.gz target='blank'>%s</a>]&nbsp;",
		htmlLinker_generateLinkToGenomeBrowserAtUCSC("hg18","vertebrate","human", chrTmp, 
			confp_get(Conf, "UCSC_GENOME_BROWSER_FLANKING_REGION"), 
			50000000 + confp_get(Conf, "UCSC_GENOME_BROWSER_FLANKING_REGION")),
			confp_get(Conf, "WEB_DATA_LINK"), 
			prefix, 
			chrTmp, 
			chrTmp); 
	if (countCol > 10) {
	  puts( "<BR>" );
	  countCol=0;
	}
	countCol++;
  }
  if( countCol==0) puts( "No data available yet" );
  ls_destroy(ls);
  puts ("<br><br>");
  puts ("For a definition of SPER, DASPER and RESPER see <a href=http://rnaseq.gersteinlab.org/fusionseq/>FusionSeq</a>");
  puts ("<br><br>");
  puts ("<br><table border=0 width=100% align=center cellpadding=10>");
  puts ("<tr align=left>");
  puts ("<th>SPER</th>");
  puts ("<th>DASPER</th>");
  puts ("<th>RESPER</th>");
  puts ("<th>Number of inter paired-end reads</th>");
  puts ("<th>Type</th>");
  puts ("<th>Genomic coordinates</th>");
  puts ("<th>Gene symbol</th>");
  puts ("<th>Description</th>");
  puts ("<th>Genomic coordinates</th>");
  puts ("<th>Gene symbol</th>");
  puts ("<th>Description</th>");
  puts ("<th></th>");
  puts ("</tr>");
  fflush(stdout);

  stringPrintf (buffer,"%s/%s.gfr", confp_get(Conf, "WEB_DATA_DIR"), prefix);
  gfr_init (string (buffer));
  int countElements = 0;
  while (currGE = gfr_nextEntry ()) {
    if (currGE->numInter < minNum) {
      continue;
    }
    if (strEqual (typeSelected,"all") || strEqual (currGE->fusionType,typeSelected) || 
	( strEqual(currGE->fusionType,"cis") && strEqual( typeSelected,"same") ) ||
	( strEqual(currGE->fusionType,"read-through") && strEqual( typeSelected,"same") ) ) {
      if (pos = strchr (currGE->descriptionTranscript1,'|')) {
        *pos = '\0';
      }
      if (pos = strchr (currGE->descriptionTranscript2,'|')) {
        *pos = '\0';
      }
      puts ("<tr>");
      printf ("<td align=left>%1.3f</td>\n",currGE->SPER);
      printf ("<td align=left>%1.3f</td>\n",currGE->DASPER);
      printf ("<td align=left>%1.3f</td>\n",currGE->RESPER);
      printf ("<td align=left>%d</td>\n",currGE->numInter);
      printf ("<td align=left>%s</td>\n",currGE->fusionType);
      printf ("<td align=left><a href=%s target=blank>%s:%d-%d</a></td>\n",
              htmlLinker_generateLinkToGenomeBrowserAtUCSC ("hg18","vertebrate","human",
			currGE->chromosomeTranscript1,
			currGE->startTranscript1 - atoi(confp_get(Conf, "UCSC_GENOME_BROWSER_FLANKING_REGION")),
			currGE->endTranscript1 + atoi(confp_get(Conf, "UCSC_GENOME_BROWSER_FLANKING_REGION"))),
     	      currGE->chromosomeTranscript1,currGE->startTranscript1,currGE->endTranscript1);
      printf ("<td align=left>%s</td>\n",processString (currGE->geneSymbolTranscript1));
      printf ("<td align=left>%s</td>\n",currGE->descriptionTranscript1);
      printf ("<td align=left><a href=%s target=blank>%s:%d-%d</a></td>\n",
              htmlLinker_generateLinkToGenomeBrowserAtUCSC ("hg18","vertebrate","human",
		     	currGE->chromosomeTranscript2,
			currGE->startTranscript2 - atoi(confp_get(Conf, "UCSC_GENOME_BROWSER_FLANKING_REGION")),
			currGE->endTranscript2 + atoi(confp_get(Conf, "UCSC_GENOME_BROWSER_FLANKING_REGION"))),
              currGE->chromosomeTranscript2,currGE->startTranscript2,currGE->endTranscript2);
      printf ("<td align=left>%s</td>\n",processString (currGE->geneSymbolTranscript2));
      printf ("<td align=left>%s</td>\n",currGE->descriptionTranscript2);
      printf ("<td align=left><a href=%s/showDetails_cgi?%s+%s>Details</a></td>\n", confp_get(Conf, "WEB_URL_CGI"), prefix,currGE->id);
      puts ("</tr>");
      countElements++;
    }
  }
  gfr_deInit ();
  stringDestroy (buffer);
  puts ("</table><br><br>");
  if( countElements == 0) puts("No fusion candidates can be found satisfying all specified criteria.");
  puts ("</body>");
  puts ("</html>");
  fflush (stdout);
}



int main (int argc, char *argv[]) 
{
  char *queryString;

  if ((Conf = confp_open(getenv("FUSIONSEQ_CONFPATH"))) == NULL)
    return EXIT_FAILURE;

  cgiInit();
  cgiHeader("text/html");
  queryString = cgiGet2Post();
  if (queryString[0] == '\0') {
    puts ("<html>");
    puts ("<head>");
    html_printGenericStyleSheet (12);
    puts ("<title>geneFusions</title>\n");
    puts ("</head>");
    puts ("<body>");
    puts ("<h1>Identification of potential gene fusions using paired-end reads</h1><br><br>");
    printf ("<form action=%s/geneFusions_cgi method=get>", confp_get(Conf, "WEB_URL_CGI"));
    puts ("<b>Data prefix</b>:&nbsp;");
    puts ("<input type=text name=prefix>");
    puts ("<br><br><br>");
    puts ("<b>Minimum number of paired-end reads connecting two genes</b>:&nbsp;");
    puts ("<select name=minNum>");
    puts ("<option value=2>2");
    puts ("<option value=3>3");
    puts ("<option value=5 selected>5");
    puts ("<option value=10>10");
    puts ("</select>");
    puts ("<br><br><br>");
    puts ("<b>Type of gene fusion</b>:&nbsp;");
    puts ("<select name=type>");
    puts ("<option value=read-through>Read-through events");
    puts ("<option value=cis>Cis events");
    puts ("<option value=intra>Intra-chromosomal events");
    puts ("<option value=same>Genes on the same chromosome");
    puts ("<option value=inter>Genes on different chromosomes");
    puts ("<option value=all selected>All potential gene fusions");
    puts ("</select>");
    puts ("<br><br><br>");
    puts ("<input type=submit value=Submit>");
    puts ("<input type=reset value=Reset>");
    puts ("</form>");
    puts ("</body>");
    puts ("</html>");
    fflush (stdout);
  }
  else {
    int first;
    Stringa item = stringCreate (20);
    Stringa value = stringCreate (20);
    char *iPtr,*vPtr,*prefix,*type;
    int minNum;

    first = 1;
    cgiGetInit ();
    while (cgiGetNextPair (&first,item,value)) {
      iPtr = string (item);
      vPtr = string (value);
      if (strEqual (iPtr,"prefix")) {
	prefix = hlr_strdup (vPtr);
      }
      if (strEqual (iPtr,"type")) {
	type = hlr_strdup (vPtr);
      }
      if (strEqual (iPtr,"minNum")) {
	minNum = atoi (vPtr);
      }
    }
    generateOutput (prefix,type,minNum);
  }
  confp_close(Conf);

  return EXIT_SUCCESS;
}
