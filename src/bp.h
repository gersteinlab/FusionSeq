#ifndef DEF_BP_H
#define DEF_BP_H



typedef struct {
  int offset;
  char *read;
} BreakPointRead;



typedef struct {
  char *tileCoordinate1;
  char *tileCoordinate2;
  Array breakPointReads; // of type BreakPointRead
} BreakPoint;



extern void bp_init (const char* fileName);
extern void bp_deInit (void);
extern Array bp_getBreakPoints (void);
extern char* bp_writeBreakPoint (BreakPoint *currBP);



#endif
