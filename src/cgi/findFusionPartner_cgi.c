#include <stdlib.h>

#include <bios/confp.h>
#include <bios/log.h>
#include <bios/format.h>
#include <bios/html.h>
#include <bios/htmlLinker.h>
#include <bios/linestream.h>
#include <bios/stringUtil.h>

#include "gfr.h"

#define MATCH_TYPE_EXACT 1
#define MATCH_TYPE_SUBSTRING 2



typedef struct {
  Stringa buffer;
  int numInter;
} Entry;

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
  textDestroy (tokens);
  return string (buffer);
}



static int sortEntries (Entry *a, Entry *b) 
{
  return b->numInter - a->numInter;
}



static void addEntry (Array entries, GfrEntry *currGE, char *sample)
{
  Entry *currEntry;
  char *pos;

  currEntry = arrayp (entries,arrayMax (entries),Entry);
  currEntry->buffer = stringCreate (1000);
  currEntry->numInter = currGE->numInter;
  if (pos = strchr (currGE->descriptionTranscript1,'|')) {
    *pos = '\0';
  }
  if (pos = strchr (currGE->descriptionTranscript2,'|')) {
    *pos = '\0';
  }
  stringAppendf (currEntry->buffer,"<tr>");
  stringAppendf (currEntry->buffer,"<td align=left><a href=%s/geneFusions_cgi?prefix=%s&minNum=5&type=all>%s</a></td>\n",
		 confp_get(Conf, "WEB_URL_CGI"),
		 sample,sample);
  stringAppendf (currEntry->buffer,"<td align=left>%1.3f</td>\n",currGE->SPER);
  stringAppendf (currEntry->buffer,"<td align=left>%1.3f</td>\n",currGE->DASPER);
  stringAppendf (currEntry->buffer,"<td align=left>%1.3f</td>\n",currGE->RESPER);
  stringAppendf (currEntry->buffer,"<td align=left>%d</td>\n",currGE->numInter);
  stringAppendf (currEntry->buffer,"<td align=left>%s</td>\n",currGE->fusionType);
  stringAppendf (currEntry->buffer,"<td align=left><a href=%s target=blank>%s:%d-%d</a></td>\n",
                 htmlLinker_generateLinkToGenomeBrowserAtUCSC ("hg18","vertebrate","human",
			 currGE->chromosomeTranscript1,
			 currGE->startTranscript1 - atoi(confp_get(Conf, "UCSC_GENOME_BROWSER_FLANKING_REGION")),
			 currGE->endTranscript1 + atoi(confp_get(Conf, "UCSC_GENOME_BROWSER_FLANKING_REGION"))),
                 currGE->chromosomeTranscript1,
		 currGE->startTranscript1,
		 currGE->endTranscript1);
  stringAppendf (currEntry->buffer,"<td align=left>%s</td>\n",processString (currGE->geneSymbolTranscript1));
  stringAppendf (currEntry->buffer,"<td align=left>%s</td>\n",currGE->descriptionTranscript1);
  stringAppendf (currEntry->buffer,"<td align=left><a href=%s target=blank>%s:%d-%d</a></td>\n",
                 htmlLinker_generateLinkToGenomeBrowserAtUCSC ("hg18","vertebrate","human",
			 currGE->chromosomeTranscript2,
			 currGE->startTranscript2 - atoi(confp_get(Conf, "UCSC_GENOME_BROWSER_FLANKING_REGION")),
			 currGE->endTranscript2 + atoi(confp_get(Conf, "UCSC_GENOME_BROWSER_FLANKING_REGION"))),
                 currGE->chromosomeTranscript2,
		 currGE->startTranscript2,
		 currGE->endTranscript2);
  stringAppendf (currEntry->buffer,"<td align=left>%s</td>\n",processString (currGE->geneSymbolTranscript2));
  stringAppendf (currEntry->buffer,"<td align=left>%s</td>\n",currGE->descriptionTranscript2);
  stringAppendf (currEntry->buffer,"<td align=left><a href=%s/showDetails_cgi?%s+%s>Details</a></td>\n", 
		 confp_get(Conf, "WEB_URL_CGI"),
		 sample,currGE->id);
  stringAppendf (currEntry->buffer,"</tr>\n");
}



static void generateOutput (char* geneName, int matchType)
{
  GfrEntry *currGE;
  Stringa cmd;
  char *pos;
  LineStream ls;
  char *line;
  char *prefix = NULL;
  char *sample = NULL;
  Array entries;
  Entry *currEntry;
  int i;
  Texta tokens;
  int found;

  puts ("<html>");
  puts ("<head>");
  puts ("<title>Find Fusion Partner</title>");
  html_printGenericStyleSheet (12);
  puts ("</head>");
  puts ("<body>");
  if (geneName[0] == '\0') {
    die ("Invalid geneName");
  }
  toupperStr (geneName);
  entries = arrayCreate (100,Entry);
  cmd = stringCreate(50); 
  stringPrintf (cmd, "ls -1 %s/*.gfr", confp_get(Conf, "WEB_DATA_DIR"));
  ls = ls_createFromPipe (string (cmd));
  while (line = ls_nextLine (ls)) {
    if (gfr_init (line) == 0) 
      continue;
    strReplace (&prefix,line);
    pos = strrchr (prefix,'.');
    *pos = '\0';
    pos = strrchr (prefix,'/');
    strReplace (&sample,pos + 1);
    while (currGE = gfr_nextEntry ()){
      if (matchType == MATCH_TYPE_EXACT) {
        tokens = textStrtokP (currGE->geneSymbolTranscript1,"|");
        i = 0;
        found = 0;
        while (i < arrayMax (tokens)) {
          if (strCaseEqual (textItem (tokens,i),geneName)) {
            addEntry (entries,currGE,sample);
            found = 1;
            break;
          }
          i++;
        }
        textDestroy (tokens);
        if (found == 0) {
          tokens = textStrtokP (currGE->geneSymbolTranscript2,"|");
          i = 0;
          while (i < arrayMax (tokens)) {
            if (strCaseEqual (textItem (tokens,i),geneName)) {
              addEntry (entries,currGE,sample);
              break;
            }
            i++;
          }
          textDestroy (tokens);
        }
      }
      else if (matchType == MATCH_TYPE_SUBSTRING) {
        if (strCaseStr (currGE->geneSymbolTranscript1,geneName) ||
            strCaseStr (currGE->geneSymbolTranscript2,geneName)) {
          addEntry (entries,currGE,sample);
        }
      }
      else {
        die ("Unknown matchType: %d",matchType);
      }
    }
    gfr_deInit ();
  }
  ls_destroy (ls);
  if (arrayMax (entries) == 0) {
    puts ("No fusion candidates can be found satisfying all specified criteria.");
    puts ("</body>");
    puts ("</html>");
    return;
  }
  arraySort (entries,(ARRAYORDERF)sortEntries);
  printf ("<h1>Results for %s</h1><br><br><br>",geneName);
  puts ("For a definition of SPER, DASPER and RESPER see <a href=http://rnaseq.gersteinlab.org/fusionseq/>FusionSeq</a>");
  puts ("<br><br>");
  puts ("<br><table border=0 width=100% align=center cellpadding=10>");
  puts ("<tr align=left>");
  puts ("<th>Sample</th>");
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
  for (i = 0; i < arrayMax (entries); i++) {
    currEntry = arrp (entries,i,Entry);
    puts (string (currEntry->buffer));
  }
  puts ("</table><br><br>");
  puts ("</body>");
  puts ("</html>");
  fflush (stdout);
  stringDestroy (cmd);
}



int main (int argc, char *argv[]) 
{
  if ((Conf = confp_open(getenv("FUSIONSEQ_CONFPATH"))) == NULL)
    return EXIT_FAILURE;

  cgiInit();
  cgiHeader("text/html");
  if (argc < 2) {
    puts ("<html>");
    puts ("<head>");
    html_printGenericStyleSheet (12);
    puts ("<title>findFusionPartner</title>\n");
    puts ("</head>");
    puts ("<body>");
    puts ("<h1>Find gene fusion partners</h1><br><br>");
    printf ("<form action=%s/findFusionPartner_cgi?process method=post>", confp_get(Conf, "WEB_URL_CGI"));
    puts ("<b>Gene name</b>:&nbsp;");
    puts ("<input type=text name=geneName>");
    puts ("<br><br><br>");
    printf ("<input type=radio name=matchType value=%d checked>Exact match&nbsp;&nbsp;\n",MATCH_TYPE_EXACT);
    printf ("<input type=radio name=matchType value=%d>Match substring\n",MATCH_TYPE_SUBSTRING);
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
    char *iPtr,*vPtr,*geneName;
    int matchType;

    first = 1;
    cgiGetInit ();
    while (cgiGetNextPair (&first,item,value)) {
      iPtr = string (item);
      vPtr = string (value);
      if (strEqual (iPtr,"geneName")) {
	geneName = hlr_strdup (vPtr);
      }
      if (strEqual (iPtr,"matchType")) {
        matchType = atoi (vPtr);
      }
    }
    generateOutput (geneName,matchType);
  }

  confp_close(Conf);
  return EXIT_SUCCESS;
}
