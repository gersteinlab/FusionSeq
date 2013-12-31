#ifndef CCUTIL_H
#define CCUTIL_H

#define ILEN 32

#define nmin(x, y) ((x<y)?x:y)
#define nmax(x, y) ((x>y)?x:y)

typedef struct {
  struct {
    int expr;
    int exons;
    int genes;
  } tracks;
  struct {
    int on;
    int thold;
    int maxgap;
    int minrun;
    int showtars;
  } rfilter;
  int readlim;
  int minspan;
} SVCfg_t;

typedef struct {
  int chromosome;
  int instance;
  int show;
  int start;
  int end;
  int mstart;
  int mend;
} SRegion_t;

typedef struct {
  int show;
  int instances;
} Chrdata_t;

/**
 * MappedRead
 */
typedef struct {
  char* chromosome;
  int start;
  int end;
} Locus;

typedef struct {
  Locus read1;
  Locus read2;
  int id;
} PEreads;


void swapInt (int *a, int *b);
unsigned int randomN (int n);

int inRegion (int istart, int iend, int rstart, int rend);
int inRegions (Array regions, int chromosome, int start, int end);

int SRegionCmp(SRegion_t *region1, SRegion_t *region2);
int PEReadSizeCmp (PEreads *per1, PEreads *per2);
int sortPosition( PEreads *per1, PEreads *per2);

char *getMchrname (int chromosome);
char *getSchrname (int chromosome);
char *getHchrname (int chromosome);

int getScale (Array regions);
int getChrnum (char *chrname);
int chrsize (int chromosome);

// String utils
int putsn (int ntabs, const char *str);
char *striappend (char *str, int value, int base);
char *strappend (char *destination, const char *source);
int antoi (const char *str, size_t n);
char *itoa (int value, char *str, int base);
void reverse (char *str);

#endif
