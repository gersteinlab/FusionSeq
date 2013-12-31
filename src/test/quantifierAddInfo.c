#include <stdlib.h>

#include <bios/confp.h>
#include <bios/log.h>
#include <bios/format.h>
#include <bios/linestream.h>

#include "util.h"


int main (int argc, char *argv[])
{
	Array kgXrefs;
	Stringa buffer;
	LineStream ls;
	int count=0;
	char* geneSymbolTranscript;
	char* descriptionTranscript;
	char* line;
	char* exonID = NULL;

	config *conf;

	if ((conf = confp_open(getenv("FUSIONSEQ_CONFPATH"))) == NULL)
		return EXIT_FAILURE;

	buffer = stringCreate (100);

	stringPrintf (buffer,"%s/%s",
		      confp_get(conf, "ANNOTATION_DIR"),
		      confp_get(conf, "KNOWN_GENE_XREF_FILENAME"));
	kgXrefs = util_readKnownGeneXrefs (string (buffer));
	arraySort (kgXrefs,(ARRAYORDERF)sortKgXrefsByTranscriptName);
	stringDestroy (buffer);

	//  gfr_init ("-");
	//  ls = ls_createFromFile("-");
  
	while (line = ls_nextLine(ls)) {
		char *lineP = hlr_strdup(line);
		WordIter w = wordIterCreate( line, "\t", 0);
		char *nameTranscript = wordNext( w );
		char *p = rindex(nameTranscript, '_');
		if (p) {
			exonID = hlr_strdup( p+1 );
			*p='\0';
		}
		transcript2geneSymbolAndGeneDescription(kgXrefs,
							nameTranscript,
							&geneSymbolTranscript,
							&descriptionTranscript);
		if (exonID) {
			printf("%s_%s\t%s\t%s\t%s", 
				nameTranscript, 
				exonID,
				geneSymbolTranscript, 
				exonID, 
				descriptionTranscript);
			hlr_free(exonID);
		} else {
			printf("%s\t%s\t1\t%s", 
				nameTranscript, 
				geneSymbolTranscript, 
				descriptionTranscript);
		}
		printf("%s\n", lineP+strlen(nameTranscript));
		count++;
		hlr_free(lineP);
		wordIterDestroy(w);
	}
	ls_destroy (ls);
	warn ("%s_numGfrEntries: %d",argv[0],count);
	confp_close(conf);

	return EXIT_SUCCESS;
}

