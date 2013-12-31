#include <stdio.h>
#include <bios/array.h>
#include <bios/html.h>
#include <bios/format.h>
#include "sqvUtil.h"
#include "sqvWeb.h"

void filterChrRegions(Array *regions, Chrdata_t **chromosomes)
{
  Array tmpArr = arrayCreate(1, SRegion_t);
  SRegion_t *tmp, *tmp2;
  int i, j = 0;

  for (i = 0; i < arrayMax(*regions); i++) {
    tmp = arrayp(*regions, i, SRegion_t);
    if (tmp->chromosome != -1 && tmp->instance != -1) {
      if (tmp->mend == -1) {
        tmp->mend = chrsize(tmp->chromosome);
      }
      if (tmp->end == -1) {
        tmp->end = chrsize(tmp->chromosome);
      }
      if (tmp->end > tmp->mend) {
        tmp->end = tmp->mend;
      }
      if (tmp->start < tmp->mstart) {
        tmp->start = tmp->mstart;
      }

      tmp2 = arrayp(tmpArr, j, SRegion_t);
      tmp2->chromosome = tmp->chromosome;
      tmp2->instance   = tmp->instance;
      tmp2->show       = tmp->show;
      tmp2->start      = tmp->start;
      tmp2->end        = tmp->end;
      tmp2->mstart     = tmp->mstart;
      tmp2->mend       = tmp->mend;

      (*chromosomes)[tmp->chromosome].instances++;
      j++;
    }
  }

  for (i = 1; i <= 24; i++) {
    if ((*chromosomes)[i].show) {
      if ((*chromosomes)[i].instances == 0) {
        (*chromosomes)[i].show = 0;
      }
    }
  }

  arraySort(tmpArr, (ARRAYORDERF) &SRegionCmp);
  arrayDestroy(*regions);

  *regions = tmpArr;
}



void web_printSingleLocusForm(char *web_url, char *web_pub_url)
{
  puts ("<html>");
  puts ("<head>");
  puts ("<title>SeqViz: Single Locus Analysis</title>\n");
  html_printGenericStyleSheet (12);
  html_printAdditionalCSS(web_pub_url);
  html_linkjQuery(web_pub_url);
  puts ("</head>");
  puts ("<body>");
  puts ("<h1>SeqViz: Identification of all reads which align to one region</h1><br><br>");
  printf ("<form action=\"%s/seqViz_cgi\" method=\"get\">\n", web_url);
  puts ("<b>Data prefix:</b>&nbsp;");
  puts ("<input type=\"text\" name=\"prefix\">");
  puts ("<br><br><br>");
  puts ("<b>Locus:</b>&nbsp;");
  puts ("<input type=\"text\" name=\"location\">");
  puts ("Example: chr21:2449394-3499493");
  puts ("<br><br><br>");
  puts ("<b>Limit number of reads:</b>&nbsp;");
  puts ("<input type=\"text\" name=\"readlim\">");
  puts ("<br><br><br>");
  puts ("<b>Minimum insert size:</b>&nbsp;");
  puts ("<input type=\"text\" name=\"minspan\">");
  puts ("<br><br><br>");
  puts ("<input type=\"submit\" value=\"Submit\">");
  puts ("<input type=\"reset\" value=\"Reset\">");
  puts ("</form>");
  puts ("</body>");
  puts ("</html>");
}

void web_printHTMLHead(Array regions, Chrdata_t *chromosomes, char *web_pub_url)
{
  int i;
  SRegion_t *tmp;

  puts ("<html>");
  puts ("<head>");
  puts ("<title>SeqViz - Single Location</title>");
  html_printGenericStyleSheet (10);
  html_printAdditionalCSS(web_pub_url);
  html_linkjQuery(web_pub_url);

  tmp = arrayp(regions, 0, SRegion_t);
  if (tmp->chromosome >= 1 && tmp->chromosome <= 24) {
    puts ("<script type=\"text/javascript\">");
    puts ("$(function() {");

    for (i = 0; i < arrayMax(regions); i++) {
      tmp = arrayp(regions, i, SRegion_t);
      web_printSliderJS (tmp);
    }

    puts ("});");
    //    puts ("\n\nvar imgLoaded=0;\n\nfunction loadImage() {\n\talert(\"loadImage\");\n\timgCircos=document.getElementById(\"circosImg\");\nalert(imgCircos);\n\tif( !isNull(imgCircos) ) { \n\t\tvar oldImg=imgCircos.src;\n\t\timgCircos.src=\"\";\n\t\timgCircos.src=oldImg;}\nelse {\n}\n\nfunction checkImage() {\n\talert(\"checkImage:\"+imgLoaded);\n\n\tif( imgLoaded===0 ) {\n\t\tloadImage();\nalert(\"called loadImage, now setTimeout\");\n\t\tsetTimeout(checkImage,1000);\n\t} else {\n\t\talert(\"img loaded\");\n\t}\n}\n\nfunction setImageLoaded() { \n\talert(\"setImageLoaded\");\n\timgLoaded=1;\n}\n\nwindow.onload=checkImage;\n\n");
    puts ("</script>");
  }
  puts ("</head>");
  puts ("<body>");
  puts ("<div class=\"container\">");
}



void web_printSliderJS (SRegion_t *region)
{
  char *schrname;
  char *mchrname;

  schrname = getSchrname (region->chromosome);
  mchrname = getMchrname (region->chromosome);

  printf ("\t$(\"#slider-%s-%i\").slider({\n", schrname, region->instance);
  puts ("\t\trange: true,");
  printf ("\t\tmin: %i,\n", region->mstart);
  printf ("\t\tmax: %i,\n", region->mend);
  printf ("\t\tvalues: [%i, %i],\n", region->start, region->end);
  puts ("\t\tslide: function(event, ui) {");
  printf ("\t\t\t$(\"#%s-%i-start\").val(ui.values[0]);\n", schrname, region->instance);
  printf ("\t\t\t$(\"#%s-%i-end\").val(ui.values[1]);\n", schrname, region->instance);
  puts ("\t\t}");
  puts ("\t});");
  printf ("\t$(\"#%s-%i-start\").val($(\"#slider-%s-%i\").slider(\"values\", 0));\n", schrname, region->instance, schrname, region->instance);
  printf ("\t$(\"#%s-%i-end\").val($(\"#slider-%s-%i\").slider(\"values\", 1));\n\n", schrname, region->instance, schrname, region->instance);

  printf ("\t$(\"#slider-%s-%i\").slider({\n", mchrname, region->instance);
  puts ("\t\trange: true,");
  printf ("\t\tmin: %i,\n", 0);
  printf ("\t\tmax: %i,\n", chrsize(region->chromosome));
  printf ("\t\tvalues: [%i, %i],\n", region->mstart, region->mend);
  puts ("\t\tslide: function(event, ui) {");
  printf ("\t\t\t$(\"#%s-%i-start\").val(ui.values[0]);\n", mchrname, region->instance);
  printf ("\t\t\t$(\"#%s-%i-end\").val(ui.values[1]);\n", mchrname, region->instance);
  puts ("\t\t}");
  puts ("\t});");
  printf ("\t$(\"#%s-%i-start\").val($(\"#slider-%s-%i\").slider(\"values\", 0));\n", mchrname, region->instance, mchrname, region->instance);
  printf ("\t$(\"#%s-%i-end\").val($(\"#slider-%s-%i\").slider(\"values\", 1));\n\n", mchrname, region->instance, mchrname, region->instance);

  free (schrname);
  free (mchrname);
}

void web_printPageHeader(char *prefix, char *location)
{
  if (prefix[0] == '\0') {
    fprintf (stderr, "Invalid prefix\n");
    exit (1);
  }
  printf ("<h1>SeqViz - %s - %s</h1>\n", prefix, location);
}


void web_printSidebar (char *web_url,
                       char *prefix,
                       char *location,
                       char *cirImgUrl,
                       char *cirSvgUrl,
                       Array regions,
                       Chrdata_t *chromosomes,
                       SVCfg_t *settings)
{
  int i;
  char *chrname;
  char *schrname;
  Stringa url;

  //--------------------------------------------------------
  // Display options block
  //--------------------------------------------------------
  putsn (1, "<div class=\"sidebar\">");
  putsn (2, "<div class=\"m_wrapper\">");
  putsn (3, "<div class=\"m_header\">Ideogram</div>");
  putsn (3, "<div class=\"m_bodyBlank\">");
  printf ("\t\t\t\t<form action=\"%s/seqViz_cgi\" method=\"get\">\n", web_url);
  putsn (5, "<table cellpadding=\"0\" cellspacing=\"0\">");
  putsn (6, "<tr>");
  putsn (7, "<td class=\"f_label\">Prefix</td>");
  putsn (7, "<td class=\"f_field\">");
  printf ("\t\t\t\t\t\t\t<input type=\"text\" name=\"prefix\" value=\"%s\" />\n", prefix);
  putsn (7, "</td>");
  putsn (6, "</tr>");
  putsn (6, "<tr>");
  putsn (7, "<td class=\"f_label\">Location</td>");
  putsn (7, "<td class=\"f_field\">");
  printf ("\t\t\t\t\t\t\t<input type=\"text\" name=\"location\" value=\"%s\" />\n", location);
  putsn (7, "</td>");
  putsn (6, "</tr>");
  putsn (6, "<tr>");
  putsn (7, "<td class=\"f_label\">Segmentation filter</td>");
  putsn (7, "<td class=\"f_field\">");

  if (settings->rfilter.on) {
    putsn (8, "<span class=\"l_no\">Off&nbsp;<input type=\"radio\" name=\"rfilter\" value=\"0\" /></span>");
    putsn (8, "<span class=\"l_yes\"><input type=\"radio\" name=\"rfilter\" value=\"1\" checked=\"checked\" />&nbsp;On</span>");
  }
  else {
    putsn (8, "<span class=\"l_no\">Off&nbsp;<input type=\"radio\" name=\"rfilter\" value=\"0\" checked=\"checked\" /></span>");
    putsn (8, "<span class=\"l_yes\"><input type=\"radio\" name=\"rfilter\" value=\"1\" />&nbsp;On</span>");
  }
  putsn (7, "</td>");
  putsn (6, "</tr>");

  if (settings->rfilter.on) {
    putsn (6, "<tr>");
    putsn (7, "<td class=\"f_cas\" colspan=\"2\">");
    putsn (8, "<fieldset>");
    putsn (9, "<legend>Filter params</legend>");
    putsn (9, "<table cellpading=\"0\" cellspacing=\"0\">");
    putsn (10, "<tr>");
    putsn (11, "<td class=\"f_label\">Threshold</td>");
    puts ("\t\t\t\t\t\t\t\t\t\t\t");
    printf ("<td class=\"f_field\"><input type=\"text\" name=\"rthold\" value=\"%i\" size=\"15\" /></td>\n", settings->rfilter.thold);
    putsn (10, "</tr>");
    putsn (10, "<tr>");
    putsn (11, "<td class=\"f_label\">Maxgap</td>");
    puts ("\t\t\t\t\t\t\t\t\t\t\t");
    printf ("<td class=\"f_field\"><input type=\"text\" name=\"rmaxgap\" value=\"%i\" size=\"15\" /></td>\n", settings->rfilter.maxgap);
    putsn (10, "</tr>");
    putsn (10, "<tr>");
    putsn (11, "<td class=\"f_label\">Minrun</td>");
    puts ("\t\t\t\t\t\t\t\t\t\t\t");
    printf ("<td class=\"f_field\"><input type=\"text\" name=\"rminrun\" value=\"%i\" size=\"15\" /></td>\n", settings->rfilter.minrun);
    putsn (10, "</tr>");
    putsn (10, "<tr>");
    putsn (11, "<td class=\"f_label\">Tars</td>");
    putsn (11, "<td class=\"f_field\">");
    if (settings->rfilter.showtars) {
      putsn (8, "<span class=\"l_no\">Hide&nbsp;<input type=\"radio\" name=\"rshowtars\" value=\"0\" /></span>");
      putsn (8, "<span class=\"l_yes\"><input type=\"radio\" name=\"rshowtars\" value=\"1\" checked=\"checked\" />&nbsp;Show</span>");
    } else {
      putsn (8, "<span class=\"l_no\">Hide&nbsp;<input type=\"radio\" name=\"rshowtars\" value=\"0\" checked=\"checked\" /></span>");
      putsn (8, "<span class=\"l_yes\"><input type=\"radio\" name=\"rshowtars\" value=\"1\" />&nbsp;Show</span>");
    }
    putsn (11, "</td>");
    putsn (10, "</tr>");
    putsn (9, "</table>");
    putsn (8, "</fieldset>");
    putsn (7, "</td>");
    putsn (6, "</tr>");
  }

  putsn (6, "<tr>");
  putsn (7, "<td class=\"f_label\">Limit reads (0 for no limit)</td>");
  putsn (7, "<td class=\"f_field\">");
  printf ("\t\t\t\t\t\t\t\t<input type=\"text\" name=\"readlim\", value=\"%i\" />", settings->readlim);
  putsn (7, "</td>");
  putsn (6, "</tr>");

  putsn (6, "<tr>");
  putsn (7, "<td class=\"f_label\">Min insert size</td>");
  putsn (7, "<td class=\"f_field\">");
  printf ("\t\t\t\t\t\t\t\t<input type=\"text\" name=\"minspan\", value=\"%i\" />", settings->minspan);
  putsn (7, "</td>");
  putsn (6, "</tr>");

  putsn (6, "<tr>");
  putsn (7, "<td class=\"f_label\" colspan=\"2\">");
  if (settings->tracks.genes) {
    putsn (8, "<input type=\"hidden\" name=\"genes\" value=\"1\" />\n");
    if (settings->tracks.exons) {
      putsn (8, "<input type=\"hidden\" name=\"exons\" value=\"1\" />\n");
    }
  }
  if (settings->tracks.expr) {
    putsn (8, "<input type=\"hidden\" name=\"expr\" value=\"on\" />\n");
  }

  web_printChrHidden (chromosomes, 8);
  web_printRegionHidden (regions, 8);

  putsn (8, "<input type=\"submit\" value=\"Update\" />");
  putsn (7, "</td>");
  putsn (6, "</tr>");

  putsn (6, "<tr>");
  putsn (7, "<td class=\"f_label\" colspan=\"2\">");
  printf ("\t\t\t\t\t\t\t\t<a href=\"%s\" target=\"_blank\">Download PNG</a><br />\n", cirImgUrl);
  printf ("\t\t\t\t\t\t\t\t<a href=\"%s\" target=\"_blank\">Download SVG</a>\n", cirSvgUrl);
  putsn (7, "</td>");
  putsn (6, "</tr>");
  putsn (5, "</table>");
  putsn (4, "</form>");
  putsn (3, "</div>");
  putsn (2, "</div>");

  //--------------------------------------------------------
  // Print chromosomes form
  //--------------------------------------------------------
  putsn (2, "<div class=\"m_wrapper\">");
  putsn (3, "<div class=\"m_header\">Chromosomes</div>");
  putsn (3, "<div class=\"m_bodyBlank\">");
  printf ("\t\t\t\t<form action=\"%s/seqViz_cgi\" method=\"get\">\n", web_url);
  putsn (5, "<table cellpadding=\"0\" cellspacing=\"0\">");
  putsn (6, "<tr>");
  putsn (7, "<td>");
  printf ("\t\t\t\t\t\t\t<input type=\"hidden\" name=\"prefix\" value=\"%s\" />\n", prefix);
  printf ("\t\t\t\t\t\t\t<input type=\"hidden\" name=\"location\" value=\"%s\" />\n", location);
  putsn (8, "<table cellpadding=\"0\" cellspacing=\"0\">");
  for (i = 1; i < 13; i++) {
    chrname = getHchrname (i);
    url = stringCreate (50);
    schrname = getSchrname (i);

    stringPrintf (url, "%s", url_makeBase (prefix, location));
    stringAppendf (url, "%s", url_makeChrs (chromosomes));
    if (!chromosomes[i].show) {
    	stringAppendf (url, "&%s=on", chrname);
    }
    stringAppendf (url, "%s", url_makeRegions (regions, 0, 0));
    stringAppendf (url, "&%s-%i-show=on", schrname, chromosomes[i].instances+1);
    stringAppendf (url, "%s", url_makeOpts (settings));

    putsn (9, "<tr>");
    putsn (10, "<td class=\"f_field\">");
    if (chromosomes[i].instances > 0) {
      if (chromosomes[i].show) {
        printf ("\t\t\t\t\t\t\t\t\t\t\t<input type=\"checkbox\" name=\"%s\" checked=\"checked\" />\n", chrname);
      }
      else {
    	  printf ("\t\t\t\t\t\t\t\t\t\t\t<input type=\"checkbox\" name=\"%s\" />\n", chrname);
      }
    }
    putsn (10, "</td>");
    putsn (10, "<td class=\"f_label\">");
    printf ("\t\t\t\t\t\t\t\t\t\t\t%s&nbsp;<a href=\"seqViz_cgi%s\" class=\"l_go\">+</a>\n", chrname, string(url));
    putsn (10, "</td>");
    putsn (9, "</tr>");

    free (schrname);
    stringDestroy (url);
    free (chrname);
  }

  putsn (8, "</table>");
  putsn (7, "</td>");
  putsn (7, "<td>");
  putsn (8, "<table cellpadding=\"0\" cellspacing=\"0\">");
  for (i = 13; i <= 24; i++) {
	chrname = getHchrname (i);
	url = stringCreate (50);
	schrname = getSchrname (i);

	 stringPrintf (url, "%s", url_makeBase (prefix, location));
	 stringAppendf (url, "%s", url_makeChrs (chromosomes));
	 if (!chromosomes[i].show) {
	   stringAppendf (url, "&%s=on", chrname);
	 }
	 stringAppendf (url, "%s", url_makeRegions (regions, 0, 0));
	 stringAppendf (url, "&%s-%i-show=on", schrname, chromosomes[i].instances+1);
	 stringAppendf (url, "%s", url_makeOpts (settings));

	putsn (9, "<tr>");
	putsn (10, "<td class=\"f_field\">");
	if (chromosomes[i].instances > 0) {
	  if (chromosomes[i].show) {
		printf ("\t\t\t\t\t\t\t\t\t\t\t<input type=\"checkbox\" name=\"%s\" checked=\"checked\" />\n", chrname);
	  }
	  else {
		  printf ("\t\t\t\t\t\t\t\t\t\t\t<input type=\"checkbox\" name=\"%s\" />\n", chrname);
	  }
	}
	putsn (10, "</td>");
	putsn (10, "<td class=\"f_label\">");
	printf ("\t\t\t\t\t\t\t\t\t\t\t%s&nbsp;<a href=\"seqViz_cgi%s\" class=\"l_go\">+</a>\n", chrname, string(url));
	putsn (10, "</td>");
	putsn (9, "</tr>");

	free (schrname);
    stringDestroy (url);
    free (chrname);
  }
  putsn (8, "</table>");
  putsn (7, "</td>");
  putsn (6, "</tr>");
  putsn (6, "<tr>");
  putsn (7, "<td class=\"f_field\" colspan=\"2\">");
  if (settings->tracks.genes) {
    putsn (8, "<input type=\"hidden\" name=\"genes\" value=\"1\" />\n");
    if (settings->tracks.exons) {
      putsn (8, "<input type=\"hidden\" name=\"exons\" value=\"1\" />\n");
    }
  }
  if (settings->tracks.expr) {
    putsn (8, "<input type=\"hidden\" name=\"expr\" value=\"on\" />\n");
  }
  if (settings->rfilter.on) {
    putsn (8, "<input type=\"hidden\" name=\"rfilter\" value=\"1\" />");
    if (settings->rfilter.thold)
      printf ("\t\t\t\t\t\t\t\t<input type=\"hidden\" name=\"rthold\" value=\"%i\" />", settings->rfilter.thold);
    if (settings->rfilter.maxgap)
      printf ("\t\t\t\t\t\t\t\t<input type=\"hidden\" name=\"rmaxgap\" value=\"%i\" />", settings->rfilter.maxgap);
    if (settings->rfilter.minrun)
      printf ("\t\t\t\t\t\t\t\t<input type=\"hidden\" name=\"rminrun\" value=\"%i\" />", settings->rfilter.minrun);
    if (settings->rfilter.showtars)
      printf ("\t\t\t\t\t\t\t\t<input type=\"hidden\" name=\"rshowtars\" value=\"%i\" />", settings->rfilter.showtars);
  }
  if (settings->readlim) {
    printf ("\t\t\t\t\t\t\t\t<input type=\"hidden\" name=\"readlim\", value=\"%i\" />", settings->readlim);
  }
  if (settings->minspan) {
    printf ("\t\t\t\t\t\t\t\t<input type=\"hidden\" name=\"minspan\", value=\"%i\" />", settings->minspan);
  }

  web_printRegionHidden (regions, 8);
  putsn (8, "<input type=\"submit\" value=\"Update\">");
  putsn (7, "</td>");
  putsn (6, "</tr>");
  putsn (5, "</table>");
  putsn (4, "</form>");
  putsn (3, "</div>");
  putsn (2, "</div>");
  putsn (1, "</div>");

}

void web_printBody (char *web_url,
                    char *prefix,
                    char *location,
                    Array regions,
                    Chrdata_t *chromosomes,
                    SVCfg_t *settings)
{
  SRegion_t *tmp;
  char *schrname;
  char *mchrname;
  char *chrname;
  Stringa url;
  int i;

  putsn (1, "<br class=\"clear\" />");
  printf ("\t<form action=\"%s/seqViz_cgi\" method=\"get\">\n", web_url);
  printf ("\t\t<input type=\"hidden\" name=\"prefix\" value=\"%s\" />\n", prefix);
  printf ("\t\t<input type=\"hidden\" name=\"location\" value=\"%s\" />\n", location);

  //--------------------------------------------------------
  // Print regions form
  //--------------------------------------------------------
  putsn (2, "<div class=\"m_wrapper\">");
  printf ("\t\t\t<div class=\"m_header\">Added regions&nbsp;<a href=\"seqViz_cgi?prefix=%s&location=%s%s\" class=\"l_warn\">- Remove all</a></div>\n",
		  prefix, location, url_makeOpts (settings));

  tmp = arrayp (regions, 0, SRegion_t);
  if (tmp->chromosome >= 1 && tmp->chromosome <= 24) {
    // The loop
    for (i = 0; i < arrayMax (regions); i++) {
      tmp = arrayp (regions, i, SRegion_t);
      schrname = getSchrname (tmp->chromosome);
      mchrname = getMchrname (tmp->chromosome);
      chrname  = getHchrname (tmp->chromosome);
      url = stringCreate (50);

      putsn (3, "<div class=\"m_body\">");
      putsn (4, "<fieldset>");
      putsn (5, "<legend>");
      if (tmp->show) {
        printf ("\t\t\t\t\t\t<input type=\"checkbox\" name=\"%s-%i-show\" checked=\"checked\" />\n", schrname, tmp->instance);
      }
      else {
        printf ("\t\t\t\t\t\t<input type=\"checkbox\" name=\"%s-%i-show\" />\n", schrname, tmp->instance);
      }
      printf ("\t\t\t\t\t\t<label for=\"%s-%i-show\">%s</label>\n", schrname, tmp->instance, chrname);
      putsn (5, "</legend>");
      putsn (5, "<table cellpadding=\"0\" cellspacing=\"0\">");
      putsn (6, "<tr>");
      putsn (7, "<td rowspan=\"2\" class=\"f_label\">Zoomed Region</td>");
      putsn (7, "<td class=\"f_field\">");
      printf ("\t\t\t\t\t\t\t\t<label for=\"%s-%i-start\">Start:</label>\n", schrname, tmp->instance);
      printf ("\t\t\t\t\t\t\t\t<input type=\"text\" name=\"%s-%i-start\" id=\"%s-%i-start\" style=\"border:0; color:#f6931f; font-weight:bold;\" />\n",
    		schrname, tmp->instance, schrname, tmp->instance);
      printf ("\t\t\t\t\t\t\t\t<label for=\"%s-%i-end\">End:</label>\n", schrname, tmp->instance);
      printf ("\t\t\t\t\t\t\t\t<input type=\"text\" name=\"%s-%i-end\" id=\"%s-%i-end\" style=\"border:0; color:#f6931f; font-weight:bold;\" />\n",
        	schrname, tmp->instance, schrname, tmp->instance);
      putsn (7, "</td>");
      putsn (6, "</tr>");
      putsn (6, "<tr>");
      putsn (7, "<td class=\"f_field slide\">");
      putsn (8, "<br class=\"clear\" />");
      printf ("\t\t\t\t\t\t\t\t<div id=\"slider-%s-%i\"></div>\n", schrname, tmp->instance);
      putsn (8, "<br class=\"clear\" />");
      putsn (7, "</td>");
      putsn (6, "</tr>");
      putsn (6, "<tr>");
      putsn (7, "<td rowspan=\"2\" class=\"f_label\">Entire Chromosome</td>");
      putsn (7, "<td class=\"f_field\">");
      printf ("\t\t\t\t\t\t\t\t<label for=\"%s-%i-start\">Start:</label>\n", mchrname, tmp->instance);
      printf ("\t\t\t\t\t\t\t\t<input type=\"text\" name=\"%s-%i-start\" id=\"%s-%i-start\" style=\"border:0; color:#f6931f; font-weight:bold;\" />\n",
        	mchrname, tmp->instance, mchrname, tmp->instance);
      printf ("\t\t\t\t\t\t\t\t<label for=\"%s-%i-end\">End:</label>\n", mchrname, tmp->instance);
      printf ("\t\t\t\t\t\t\t\t<input type=\"text\" name=\"%s-%i-end\" id=\"%s-%i-end\" style=\"border:0; color:#f6931f; font-weight:bold;\" />\n",
        	mchrname, tmp->instance, mchrname, tmp->instance);
      putsn (7, "</td>");
      putsn (6, "</tr>");
      putsn (6, "<tr>");
      putsn (7, "<td class=\"f_field slide\">");
      putsn (8, "<br class=\"clear\" />");
      printf ("\t\t\t\t\t\t\t\t<div id=\"slider-%s-%i\"></div>\n", mchrname, tmp->instance);
      putsn (8, "<br class=\"clear\" />");
      putsn (7, "</td>");
      putsn (6, "</tr>");
      putsn (6, "<tr>");
      putsn (6, "<td class=\"f_label\">Options</td>");
      putsn (6, "<td class=\"f_field\">");

      stringPrintf (url, "%s%s%s%s",
    		        url_makeBase (prefix, location),
    		        url_makeChrs (chromosomes),
    		        url_makeRegions (regions, tmp->chromosome, tmp->instance),
      		        url_makeOpts (settings));

      printf ("<a href=\"seqViz_cgi%s\" class=\"l_warn\">- Remove</a>&nbsp;\n", string (url));
      putsn (7, "</td>");
      putsn (6, "</tr>");
      putsn (5, "</table>");
      putsn (4, "</fieldset>");
      putsn (3, "</div>");

      free (chrname);
      free (schrname);
      free (mchrname);
      stringDestroy (url);
    }
  }
  putsn (3, "<div class=\"m_body\">");
  // TODO: hidden inputs if split to 2 forms?
  web_printChrHidden (chromosomes, 4);
  if (settings->rfilter.on) {
    putsn (8, "<input type=\"hidden\" name=\"rfilter\" value=\"1\" />");
    if (settings->rfilter.thold)
	  printf ("\t\t\t\t\t\t\t\t<input type=\"hidden\" name=\"rthold\" value=\"%i\" />", settings->rfilter.thold);
	if (settings->rfilter.maxgap)
	  printf ("\t\t\t\t\t\t\t\t<input type=\"hidden\" name=\"rmaxgap\" value=\"%i\" />", settings->rfilter.maxgap);
	if (settings->rfilter.minrun)
	  printf ("\t\t\t\t\t\t\t\t<input type=\"hidden\" name=\"rminrun\" value=\"%i\" />", settings->rfilter.minrun);
	if (settings->rfilter.minrun)
	  printf ("\t\t\t\t\t\t\t\t<input type=\"hidden\" name=\"rminrun\" value=\"%i\" />", settings->rfilter.minrun);
  }
  if (settings->readlim) {
    printf ("\t\t\t\t\t\t\t\t<input type=\"hidden\" name=\"readlim\", value=\"%i\" />", settings->readlim);
  }
  if (settings->minspan) {
    printf ("\t\t\t\t\t\t\t\t<input type=\"hidden\" name=\"minspan\", value=\"%i\" />", settings->minspan);
  }
  putsn (4, "<input type=\"submit\" value=\"Update\" />");
  putsn (3, "</div>");
  putsn (2, "</div>");

  //--------------------------------------------------------
  // Print data tracks form
  //--------------------------------------------------------
  putsn (2, "<div class=\"m_wrapper\">");
  putsn (3, "<div class=\"m_header\">Data Tracks</div>");
  putsn (3, "<div class=\"m_bodyBlank\">");
  putsn (4, "<table cellpadding=\"0\" cellspacing=\"0\">");
  putsn (5, "<tr>");
  if (settings->tracks.expr) {
    putsn (6, "<td class=\"f_field\"><input type=\"checkbox\" name=\"expr\" checked=\"checked\" /></td>");
  }
  else {
    putsn (6, "<td class=\"f_field\"><input type=\"checkbox\" name=\"expr\" /></td>");
  }
  putsn (6, "<td class=\"f_label\" colspan=\"2\"><label for=\"expr\">Gene Expression</label></td>");
  putsn (5, "</tr>");
  putsn (5, "<tr>");
  if (settings->tracks.genes) {
    putsn (6, "<td class=\"f_field\"><input type=\"checkbox\" name=\"genes\" checked=\"checked\" /></td>");
    putsn (6, "<td class=\"f_label\"><label for=\"genes\">Genes</label></td>");
    putsn (6, "<td class=\"f_label\">");
    if (settings->tracks.exons) {
      putsn (7, "<span class=\"l_no\">Hide exons&nbsp;<input type=\"radio\" name=\"exons\" value=\"0\" /></span>");
      putsn (7, "<span class=\"l_yes\"><input type=\"radio\" name=\"exons\" value=\"1\" checked=\"checked\" />&nbsp;Show exons</span>");
    }
    else {
      putsn (7, "<span class=\"l_no\">Hide exons&nbsp;<input type=\"radio\" name=\"exons\" value=\"0\" checked=\"checked\" /></span>");
      putsn (7, "<span class=\"l_yes\"><input type=\"radio\" name=\"exons\" value=\"1\" />&nbsp;Show exons</span>");
    }
    putsn (6, "</td>");
  }
  else {
    putsn (6, "<td class=\"f_field\"><input type=\"checkbox\" name=\"genes\" /></td>");
    putsn (6, "<td class=\"f_label\" colspan=\"2\"><label for=\"gexons\">Genes</label></td>");
  }
  putsn (5, "</tr>");
  putsn (5, "<tr>");
  putsn (6, "<td class=\"f_field\"><input type=\"checkbox\" name=\"isoforms\" disabled=\"disabled\" /></td>");
  putsn (6, "<td class=\"f_label\" colspan=\"2\">");
  putsn (7, "<fieldset>");
  putsn (8, "<legend>Isoforms</legend>");
  putsn (8, "<table cellspacing=\"0\" cellpadding=\"0\">");
  putsn (9, "<tr>");
  putsn (10, "<td class=\"f_field\">Coming Soon</td>");
  putsn (9, "</tr>");
  putsn (8, "</table>");
  putsn (7, "</fieldset>");
  putsn (6, "</td>");
  putsn (5, "</tr>");
  putsn (4, "</table>");
  putsn (3, "</div>");
  putsn (3, "<div class=\"m_body\">");
  putsn (4, "<input type=\"submit\" value=\"Update\" />");
  putsn (3, "</div>");
  putsn (2, "</div>");
  putsn (1, "</form>");
  puts ("</div>");
  puts ("</body>");
  puts ("</html>");
}

void web_printChrHidden (Chrdata_t *chromosomes, int ntabs)
{
  char *chrname;
  int i, j;

  for (i = 1; i <= 24; i++) {
    if (chromosomes[i].show) {
      chrname = getHchrname (i);

      for (j = 0; j < ntabs; j++) {
        printf ("\t");
      }
      printf ("<input type=\"hidden\" name=\"%s\" value=\"on\" />\n", chrname);
      free (chrname);
    }
  }
}

void web_printRegionHidden (Array regions, int ntabs)
{
  char *mchrname;
  char *schrname;
  Stringa tabs = stringCreate(ntabs+1);
  SRegion_t *tmp;
  int i, j;

  stringPrintf (tabs, "");

  for (j = 0; j < ntabs; j++) {
    stringAppendf (tabs, "\t");
  }

  tmp = arrayp(regions, 0, SRegion_t);
  if (tmp->chromosome >= 1 && tmp->chromosome <= 24) {
    for (i = 0; i < arrayMax(regions); i++) {
      tmp = arrayp(regions, i, SRegion_t);
      mchrname = getMchrname (tmp->chromosome);
      schrname = getSchrname (tmp->chromosome);

      if (tmp->show) {
        printf ("%s<input type=\"hidden\" name=\"%s-%i-show\" value=\"on\" />\n",
            string (tabs), schrname, tmp->instance);
      }
      printf ("%s<input type=\"hidden\" name=\"%s-%i-start\" value=\"%i\" />\n",
          string (tabs), schrname, tmp->instance, tmp->start);
      printf ("%s<input type=\"hidden\" name=\"%s-%i-end\" value=\"%i\" />\n",
          string (tabs), schrname, tmp->instance, tmp->end);
      printf ("%s<input type=\"hidden\" name=\"%s-%i-start\" value=\"%i\" />\n",
          string (tabs), mchrname, tmp->instance, tmp->mstart);
      printf ("%s<input type=\"hidden\" name=\"%s-%i-end\" value=\"%i\" />\n",
          string (tabs), mchrname, tmp->instance, tmp->mend);

      free (mchrname);
      free (schrname);
    }
  }

  stringDestroy (tabs);
}

char *url_makeBase (char *prefix, char *location)
{
  char *url = strappend (strdup ("?prefix="), prefix);
  char *loc = strappend (strdup ("&location="), location);
  url = strappend (url, loc);

  free (loc);
  return url;
}

char *url_makeOpts (SVCfg_t *settings)
{
  Stringa url = stringCreate (50);
  char *u_url;

  stringPrintf (url, "");

  if (settings->tracks.expr) {
    stringAppendf (url, "&expr=on");
  }
  if (settings->tracks.genes) {
	stringAppendf (url, "&genes=on");
    if (settings->tracks.exons) {
      stringAppendf (url, "&exons=1");
    }
  }
  if (settings->rfilter.on) {
    stringAppendf (url, "&rfilter=1");
    if (settings->rfilter.thold)
      stringAppendf (url, "&rthold=%i", settings->rfilter.thold);
    if (settings->rfilter.maxgap)
      stringAppendf (url, "&rmaxgap=%i", settings->rfilter.maxgap);
    if (settings->rfilter.minrun)
      stringAppendf (url, "&rminrun=%i", settings->rfilter.minrun);
    if (settings->rfilter.showtars)
      stringAppendf (url, "&rshowtars=%i", settings->rfilter.showtars);
  }
  if (settings->readlim) {
    stringAppendf (url, "&readlim=%i", settings->readlim);
  }
  if (settings->minspan) {
    stringAppendf (url, "&minspan=%i", settings->minspan);
  }

  u_url = strdup (string (url));
  stringDestroy (url);

  return u_url;
}

char *url_makeChrs (Chrdata_t *chromosomes)
{
  Stringa buffer = stringCreate (50);
  char *chrname;
  char *url;
  int i;

  stringPrintf (buffer, "");

  for (i = 1; i <= 24; i++) {
    if (chromosomes[i].show) {
      chrname = getHchrname (i);
      stringAppendf (buffer, "&%s=on", chrname);
      free (chrname);
    }
  }

  url = strdup (string (buffer));
  stringDestroy (buffer);
  return url;
}

char *url_makeRegions (Array regions, int exchr, int exinst)
{
  Stringa buffer = stringCreate (50);
  SRegion_t *tmp;
  char *schrname;
  char *mchrname;
  char *url;
  int i;

  stringPrintf (buffer, "");

  for (i = 0; i < arrayMax(regions); i++) {
    tmp = arrayp (regions, i, SRegion_t);

    if (tmp->chromosome != exchr || tmp->instance != exinst) {
      schrname = getSchrname (tmp->chromosome);
      mchrname = getMchrname (tmp->chromosome);

      if (tmp->show) {
        stringAppendf (buffer, "&%s-%i-show=on", schrname, tmp->instance);
      }
      stringAppendf (buffer, "&%s-%i-start=%i", schrname, tmp->instance, tmp->start);
      stringAppendf (buffer, "&%s-%i-end=%i", schrname, tmp->instance, tmp->end);
      stringAppendf (buffer, "&%s-%i-start=%i", mchrname, tmp->instance, tmp->mstart);
      stringAppendf (buffer, "&%s-%i-end=%i", mchrname, tmp->instance, tmp->mend);

      free (schrname);
      free (mchrname);
    }
  }

  url = strdup (string (buffer));
  stringDestroy (buffer);
  return url;
}



void html_printAdditionalCSS(char *web_pub_url)
{
  printf ("<link rel=\"stylesheet\" href=\"%s/css/style.css\" type=\"text/css\" media=\"all\" />\n", web_pub_url);
}

void html_linkjQuery(char *web_pub_url)
{
  printf ("<link rel=\"stylesheet\" href=\"%s/css/start/jquery-ui-1.7.2.custom.css\" type=\"text/css\" media=\"all\" />\n", web_pub_url);
  printf ("<script src=\"%s/js/jquery-1.3.2.min.js\" type=\"text/javascript\"></script>\n", web_pub_url);
  printf ("<script src=\"%s/js/jquery-ui-1.7.2.custom.min.js\" type=\"text/javascript\"></script>\n", web_pub_url);
}
