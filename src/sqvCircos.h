#ifndef CCCIRCOS_H
#define CCCIRCOS_H

#define HEATMP_COLORS "vlgrey,grey,lgreen,green,lyellow,yellow,vlorange,orange,red"
#define HEATMP_WIDTH 50

#define DATTRK_WIDTH 0.04
#define DATTRK_SPACER 0.02
#define TILE_LAYERS 4

void conf_printHeader (FILE *fp,
		               char *circos_dir,
		               char *data_dir,
		               char *sdata_dir,
		               Chrdata_t *chromosomes,
		               char *prefix,
		               Locus locus,
		               int scale);
void conf_printUnits (FILE *fp, Array regions, Chrdata_t *chromosomes, int scale);
void conf_printDataTracks (FILE *fp,
		                   char *prefix,
		                   Locus locus,
		                   char *sdata_dir,
		                   char *data_dir,
		                   float *rpos,
		                   Array regions,
		                   Chrdata_t *chromosomes,
		                   SVCfg_t *settings);
void conf_printDatGenes (FILE *fp, char *sdata_dir, Array regions, Chrdata_t *chromosomes, float *rpos);
void conf_printDatExons (FILE *fp, char *sdata_dir, Array regions, Chrdata_t *chromosomes, float *rpos);
void conf_printFilterTars (FILE *fp, char *data_dir, char *prefix, Locus locus, float *rpos);
void conf_printLinks (FILE *fp, char *data_dir, float *rpos, char *prefix, Locus locus, int readlim);
void conf_printFooter (FILE *fp);

void incl_getGeneHlightFile (FILE *fp, Array regions, char *sdata_dir);
void incl_getExonHlightFile (FILE *fp, Array regions, char *sdata_dir);

#endif
