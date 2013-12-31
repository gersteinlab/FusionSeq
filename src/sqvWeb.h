#ifndef CCWEB_H
#define CCWEB_H

void filterChrRegions(Array *regions, Chrdata_t **chromosomes);

void web_printSingleLocusForm(char *web_url, char *web_pub_url);
void web_printHTMLHead(Array regions, Chrdata_t *chromosomes, char *web_pub_url);
void web_printSliderJS (SRegion_t *region);
void web_printPageHeader(char *prefix, char *location);
void web_printChrHidden (Chrdata_t *chromosomes, int ntabs);
void web_printRegionHidden (Array regions, int ntabs);
void web_printSidebar (char *web_url,
                       char *prefix,
                       char *location,
                       char *cirImgUrl,
                       char *cirSvgUrl,
                       Array regions,
                       Chrdata_t *chromosomes,
                       SVCfg_t *settings);
void web_printBody (char *web_url,
                    char *prefix,
                    char *location,
                    Array regions,
                    Chrdata_t *chromosomes,
                    SVCfg_t *settings);

char *url_makeBase (char *prefix, char *location);
char *url_makeOpts (SVCfg_t *settings);
char *url_makeChrs (Chrdata_t *chromosomes);
char *url_makeRegions (Array regions, int exchr, int exinst);

void html_printAdditionalCSS(char *web_pub_url);
void html_linkjQuery(char *web_pub_url);

#endif
