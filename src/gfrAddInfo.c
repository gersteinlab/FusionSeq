#include <stdlib.h>

#include <bios/confp.h>
#include <bios/log.h>
#include <bios/format.h>

#include "util.h"
#include "gfr.h"



int main (int argc, char *argv[])
{
	GfrEntry *currGE;
	Array kgXrefs;
	Stringa buffer;
	int count;

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

	count = 0;
	gfr_init ("-");
	gfr_addNewColumnType (GFR_COLUMN_NAME_GENE_SYMBOL_TRANSCRIPT1);
	gfr_addNewColumnType (GFR_COLUMN_NAME_GENE_SYMBOL_TRANSCRIPT2);
	gfr_addNewColumnType (GFR_COLUMN_NAME_DESCRIPTION_TRANSCRIPT1);
	gfr_addNewColumnType (GFR_COLUMN_NAME_DESCRIPTION_TRANSCRIPT2);
	puts (gfr_writeHeader ());
	
	while (currGE = gfr_nextEntry ()){
		transcript2geneSymbolAndGeneDescription (kgXrefs,currGE->nameTranscript1,&currGE->geneSymbolTranscript1,&currGE->descriptionTranscript1);
		transcript2geneSymbolAndGeneDescription (kgXrefs,currGE->nameTranscript2,&currGE->geneSymbolTranscript2,&currGE->descriptionTranscript2);
		puts (gfr_writeGfrEntry (currGE));
		count++;
	}
	gfr_deInit ();
	warn ("%s_numGfrEntries: %d",argv[0],count);
	confp_close(conf);

	return EXIT_SUCCESS;
}

