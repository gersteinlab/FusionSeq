#include <stdio.h>
#include <stdlib.h>

#include <bios/array.h>
#include <bios/format.h>
#include <bios/linestream.h>
#include <bios/log.h>

#include "sqvUtil.h"
#include "sqvCircos.h"

void conf_printHeader (FILE *fp,
		               char *circos_dir,
		               char *data_dir,
		               char *sdata_dir,
		               Chrdata_t *chromosomes,
		               char *prefix,
		               Locus locus,
		               int scale)
{
  fprintf (fp, "<colors>\n<<include %s/etc/colors.conf>>\n</colors>\n\n", circos_dir);
  fprintf (fp, "<fonts>\n<<include %s/etc/fonts.conf>>\n</fonts>\n\n", circos_dir);

  fprintf (fp, "<<include %s/test/ideogram.conf>>\n<<include %s/ticks%i.conf>>\n\n", circos_dir, sdata_dir, scale);

  if (chromosomes[0].instances == 1) {
    if (chromosomes[0].show == 23) {
      fprintf (fp, "karyotype = %s/karyotype.humanX.txt\n\n<image>\n", sdata_dir);
    }
    else if (chromosomes[0].show == 24) {
      fprintf (fp, "karyotype = %s/karyotype.humanY.txt\n\n<image>\n", sdata_dir);
    }
    else {
      fprintf (fp, "karyotype = %s/karyotype.human%i.txt\n\n<image>\n", sdata_dir, chromosomes[0].show);
    }
  }
  else {
    // Include if necessary
    fprintf (fp, "karyotype = %s/karyotype.human.txt\n\n<image>\n", sdata_dir);
  }

  fprintf (fp, "dir = %s/tmp\nfile = %s_%s_%d_%d.png\n", data_dir, prefix, locus.chromosome, locus.start, locus.end);
  fprintf (fp, "svg = yes\npng = yes\n");
  fprintf (fp, "radius = 1500p\n");
  fprintf (fp, "background     = white\n");
  fprintf (fp, "angle_offset   = -90\n</image>\n\n");
}

void conf_printUnits (FILE *fp, Array regions, Chrdata_t *chromosomes, int scale)
{
  int i;
  SRegion_t *tmp;
  char *chrname;
  int sflag = 1;
  Stringa chrshow = stringCreate (50);

  tmp = arrayp (regions, 0, SRegion_t);
  if (tmp->chromosome == 0) {
    fprintf (fp, "chromosomes_display_default = yes\n");
    fprintf (fp, "chromosomes_units = 1000000\n");
  }
  else {
    stringPrintf (chrshow, "chromosomes = ");

    fprintf (fp, "chromosomes_units = %i\n", scale);
    fprintf (fp, "chromosomes_display_default = no\n");

    for (i = 0; i < arrayMax (regions); i++) {
      tmp = arrayp (regions, i, SRegion_t);
      if (chromosomes[tmp->chromosome].show && tmp->show) {
        chrname = getHchrname (tmp->chromosome);

        if (sflag == 0) {
          stringAppendf (chrshow, ";");
        }
        stringAppendf (chrshow, "%s:%i-%i", chrname, tmp->start/scale, tmp->end/scale);
        sflag = 0;

        free (chrname);
      }
    }
    fprintf (fp, "%s\n", string (chrshow));
  }

  stringDestroy (chrshow);
}



void conf_printDataTracks (FILE *fp,
		                   char *prefix,
		                   Locus locus,
		                   char *sdata_dir,
		                   char *data_dir,
		                   float *rpos,
		                   Array regions,
		                   Chrdata_t *chromosomes,
		                   SVCfg_t *settings)
{
  if (settings->tracks.expr || settings->tracks.genes || settings->tracks.exons || settings->rfilter.showtars) {
    fprintf (fp, "<plots>\nz = 5\nlayers_overflow = hide\n");
    if (settings->tracks.genes) {
      conf_printDatGenes (fp, sdata_dir, regions, chromosomes, rpos);
      if (settings->tracks.exons) {
        conf_printDatExons (fp, sdata_dir, regions, chromosomes, rpos);
      }
    }
    if (settings->rfilter.showtars) {
      conf_printFilterTars (fp, data_dir, prefix, locus, rpos);
    }
    fprintf (fp, "</plots>\n");
  }
}

void conf_printDatGenes (FILE *fp, char *sdata_dir, Array regions, Chrdata_t *chromosomes, float *rpos)
{
  float npos = *rpos - DATTRK_WIDTH;

  fprintf (fp, "<plot>\nshow = yes\ntype = tile\n");

  incl_getGeneHlightFile (fp, regions, sdata_dir);

  fprintf (fp, "orientation = out\nlayers = %i\nmargin = 0.02u\n\n", TILE_LAYERS);
  fprintf (fp, "thickness = 20\npadding = 5\ncolor = blue\n\n");
  fprintf (fp, "stroke_thickness = 1\nstroke_color     = blue\n");
  fprintf (fp, "r0    = %1.2fr\nr1    = %1.2fr\n\n", npos, *rpos);

  fprintf (fp, "</plot>\n");

  npos -= DATTRK_SPACER;
  *rpos = npos;
}

void conf_printDatExons (FILE *fp, char *sdata_dir, Array regions, Chrdata_t *chromosomes, float *rpos)
{
  float npos = *rpos - DATTRK_WIDTH;

  fprintf (fp, "<plot>\nshow = yes\ntype = tile\n");

  incl_getExonHlightFile (fp, regions, sdata_dir);

  fprintf (fp, "orientation = in\nlayers = %i\nmargin = 0.02u\n\n", TILE_LAYERS);
  fprintf (fp, "thickness = 20\npadding = 5\n\n");
  fprintf (fp, "stroke_thickness = 1\nstroke_color     = dgreen\ncolor = green\n\n");

  fprintf (fp, "r0    = %1.2fr\nr1    = %1.2fr\n", npos, *rpos);

  fprintf (fp, "</plot>\n");

  npos -= DATTRK_SPACER;
  *rpos = npos;
}

void conf_printFilterTars (FILE *fp, char *data_dir, char *prefix, Locus locus, float *rpos)
{
  float npos = *rpos - DATTRK_WIDTH;

  fprintf (fp, "<plot>\nshow = yes\ntype = tile\n");
  fprintf (fp, "file = %s/tmp/%s_%s_%d_%d.hlight.txt\n\n",
  		  data_dir, prefix, locus.chromosome, locus.start, locus.end );
  fprintf (fp, "orientation = in\nlayers = %i\nmargin = 0.02u\n\n", TILE_LAYERS);
  fprintf (fp, "thickness = 20\npadding = 5\n\n");
  fprintf (fp, "stroke_thickness = 1\nstroke_color     = dred\ncolor = red\n\n");

  fprintf (fp, "r0    = %1.2fr\nr1    = %1.2fr\n", npos, *rpos);

  fprintf (fp, "</plot>\n");

  npos -= DATTRK_SPACER;
  *rpos = npos;
}

void conf_printLinks (FILE *fp, char *data_dir, float *rpos, char *prefix, Locus locus, int readlim)
{
  fprintf (fp, "<links>\n\nz = 0\nradius = %1.2fr\nbezier_radius = 0.2r\n\n", *rpos);

  fprintf (fp, "<link data>\nshow = yes\ncolor = vvdgrey\nthickness = 2\n");
  if (readlim) {
    fprintf (fp, "record_limit = %i\n", readlim);
  }
  fprintf (fp, "file = %s/tmp/%s_%s_%d_%d_s.data\n</link>\n\n",
		  data_dir, prefix, locus.chromosome, locus.start, locus.end );
  fprintf (fp, "</links>\n\n");
}

void conf_printFooter (FILE *fp)
{
  fprintf (fp, "anglestep = 0.5\nminslicestep = 10\nbeziersamples = 40\ndebug = no\nwarnings = no\nimagemap = no\n\n");
  fprintf (fp, "units_ok     = bupr\nunits_nounit = n\n\n");
}

void incl_getGeneHlightFile (FILE *fp, Array regions, char *sdata_dir)
{
  LineStream src;
  FILE *out;
  char *line;
  Texta entry;
  int i, astart, aend;

  Stringa buffer = stringCreate (50);

  stringPrintf (buffer, "%s/tmp/genes.hlight_s.txt", sdata_dir);
  if (!(out = fopen (string (buffer), "w"))) {
    fprintf (stderr, "Cannot open genes.hlight_s.txt\n");
    return;
  }

  SRegion_t *tmp;
  tmp = arrayp (regions, 0, SRegion_t);

  if (tmp->chromosome == 0) {
    fprintf (fp, "file = %s/genes.hlight.txt\n", sdata_dir);
  }
  else {
    for (i = 0; i < arrayMax (regions); i++) {
      tmp = arrayp (regions, i, SRegion_t);

      if (tmp->chromosome == 23) {
        stringPrintf (buffer, "%s/X/genes.hlight.txt", sdata_dir);
      }
      else if (tmp->chromosome == 24) {
        stringPrintf (buffer, "%s/Y/genes.hlight.txt", sdata_dir);
      }
      else {
        stringPrintf (buffer, "%s/%i/genes.hlight.txt", sdata_dir, tmp->chromosome);
      }

      if ((src = ls_createFromFile (string (buffer))) == NULL) {
        warn( "Cannot open genes.hlight.txt %s\n", string(buffer));
        return;
      }

      while ((line = ls_nextLine (src)) != NULL) {
        entry = textFieldtokP (line, " ");

        astart = atoi (textItem (entry, 1));
        aend   = atoi (textItem (entry, 2));

        if ((astart >= tmp->start && astart <= tmp->end) ||
            (aend >= tmp->start && aend <= tmp->end)) {
          fprintf (out, "%s\n", line);
        }
        textDestroy (entry);
      }
    }

    fprintf (fp, "file = %s/tmp/genes.hlight_s.txt\n", sdata_dir);
  }

  stringDestroy (buffer);
  ls_destroy (src);
  fclose (out);
}

void incl_getExonHlightFile (FILE *fp, Array regions, char *sdata_dir)
{
  LineStream src;
  FILE *out;
  char *line;
  Texta entry;
  int i, astart, aend;

  Stringa buffer = stringCreate (50);

  stringPrintf (buffer, "%s/tmp/exons.hlight_s.txt", sdata_dir);
  if (!(out = fopen (string (buffer), "w"))) {
	fprintf (stderr, "Cannot open exons.hlight_s.txt\n");
	return;
  }

  SRegion_t *tmp;
  tmp = arrayp (regions, 0, SRegion_t);

  if (tmp->chromosome == 0) {
	fprintf (fp, "file = %s/exons.hlight.txt\n", sdata_dir);
  }
  else {
	for (i = 0; i < arrayMax (regions); i++) {
	  tmp = arrayp (regions, i, SRegion_t);

	  if (tmp->chromosome == 23) {
		stringPrintf (buffer, "%s/X/exons.hlight.txt", sdata_dir);
	  }
	  else if (tmp->chromosome == 24) {
		stringPrintf (buffer, "%s/Y/exons.hlight.txt", sdata_dir);
	  }
	  else {
		stringPrintf (buffer, "%s/%i/exons.hlight.txt", sdata_dir, tmp->chromosome);
	  }

	  if ((src = ls_createFromFile (string (buffer))) == NULL) {
		fprintf (stderr, "Cannot open exons.hlight.txt\n");
		return;
	  }

	  while ((line = ls_nextLine (src)) != NULL) {
		entry = textFieldtokP (line, " ");

		astart = atoi (textItem (entry, 1));
		aend   = atoi (textItem (entry, 2));

		if ((astart >= tmp->start && astart <= tmp->end) ||
			(aend >= tmp->start && aend <= tmp->end)) {
		  fprintf (out, "%s\n", line);
		}
		textDestroy (entry);
	  }
	}

	fprintf (fp, "file = %s/tmp/exons.hlight_s.txt\n", sdata_dir);
  }

  stringDestroy (buffer);
  ls_destroy (src);
  fclose (out);
}
