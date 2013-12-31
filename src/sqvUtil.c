#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <bios/format.h>
#include <bios/array.h>

#include "sqvUtil.h"
#include "sqvWeb.h"

/**
 * swap two integers
 */
void swapInt ( int *a, int *b)
{
  int tmp = *a;
  *a = *b;
  *b = tmp;
};

int inRegion (int istart, int iend, int rstart, int rend)
{
  if ((istart >= rstart && istart <= rend) &&
      (iend >= rstart && iend <= rend)) {
    return 1;
  }
  return 0;
}

unsigned int randomN (int n)
{
  return random() / ((2 << 30) - 1) * n;
}


int inRegions (Array regions, int chromosome, int start, int end)
{
  SRegion_t *tmp;
  int i;

  tmp = arrayp (regions, 0, SRegion_t);
  if (tmp->chromosome == 0) {
    return 1;
  }
  else {
    for (i = 0; i < arrayMax (regions); i++) {
      tmp = arrayp (regions, i, SRegion_t);
      if (tmp->chromosome == chromosome && tmp->show) {
    	if ((start >= tmp->start && start <= tmp->end) &&
            (end >= tmp->start && end <= tmp->end)) {
          return 1;
    	}
      }
    }
    return 0;
  }
}

int SRegionCmp (SRegion_t *region1, SRegion_t *region2)
{
  if (region1->chromosome == region2->chromosome) {
    return region1->instance - region2->instance;
  } else {
    return region1->chromosome - region2->chromosome;
  }
}

int sortPosition( PEreads *per1, PEreads *per2) {
  int chrnum11, chrnum12, chrnum21, chrnum22;
  int start1, start2, chr1, chr2;
  // first sort the PEreads individually -- they should be already sorted
  if (!strcmp (per1->read1.chromosome+3, "X")) {
    chrnum11 = 23;
  } else {
    if ( !strcmp (per1->read1.chromosome+3, "Y")) {
      chrnum11 = 24; 
    } else {
      chrnum11 = atoi (per1->read1.chromosome+3 );
    }
  }
  if (!strcmp (per1->read2.chromosome+3, "X")) {
    chrnum12 = 23;
  } else {
    if ( !strcmp (per1->read2.chromosome+3, "Y")) {
      chrnum12 = 24; 
    } else {
      chrnum12 = atoi (per1->read2.chromosome+3 );
    }
  }
  if( chrnum11 == chrnum12 ) {
    start1 = MIN( per1->read1.start,  per1->read2.start);
    chr1 = chrnum11;
  } else { 
    chr1 = MIN (chrnum11, chrnum12 );
    if( chrnum11 < chrnum12 )
      start1 = per1->read1.start;
    else
      start1 = per1->read2.start;
  }
  // end pair 1

  if (!strcmp (per2->read1.chromosome+3, "X")) {
    chrnum21 = 23;
  } else {
    if ( !strcmp (per2->read1.chromosome+3, "Y")) {
      chrnum21 = 24; 
    } else {
      chrnum21 = atoi (per2->read1.chromosome+3 );
    }
  }
  if (!strcmp (per2->read2.chromosome+3, "X")) {
    chrnum22 = 23;
  } else {
    if ( !strcmp (per2->read2.chromosome+3, "Y")) {
      chrnum22 = 24; 
    } else {
      chrnum22 = atoi (per2->read2.chromosome+3 );
    }
  }
  if( chrnum21 == chrnum22 ) {
    start2 = MIN( per2->read2.start,  per2->read2.start);
    chr2 = chrnum21;
  } else { 
    chr2 = MIN (chrnum21, chrnum22 );
    if( chrnum21 < chrnum22 )
      start2 = per2->read1.start;
    else
      start2 = per2->read2.start;
  }
  // end pair2
  if( chr1 == chr2 )
    return start1 - start2;
  else 
    return chr1 - chr2;
}

int PEReadSizeCmp (PEreads *per1, PEreads *per2)
{
  int chrnum11, chrnum12, chrnum21, chrnum22;
  int instsize1, instsize2;

  if (strcmp (per1->read1.chromosome+3, per1->read2.chromosome+3) ||
	  strcmp (per2->read1.chromosome+3, per2->read2.chromosome+3)) {
    if (!strcmp (per1->read1.chromosome+3, "X")) {
      chrnum11 = 23;
    } else if (!strcmp (per1->read1.chromosome+3, "Y")) {
      chrnum11 = 24;
    } else {
      chrnum11 = atoi (per1->read1.chromosome+3);
    }

    if (!strcmp (per1->read2.chromosome+3, "X")) {
      chrnum12 = 23;
    } else if (!strcmp (per1->read2.chromosome+3, "Y")) {
      chrnum12 = 24;
    } else {
      chrnum12 = atoi (per1->read2.chromosome+3);
    }

    if (!strcmp (per2->read1.chromosome+3, "X")) {
	  chrnum21 = 23;
	} else if (!strcmp (per2->read1.chromosome+3, "Y")) {
	  chrnum21 = 24;
	} else {
	  chrnum21 = atoi (per2->read1.chromosome+3);
	}

	if (!strcmp (per2->read2.chromosome+3, "X")) {
	  chrnum22 = 23;
	} else if (!strcmp (per2->read2.chromosome+3, "Y")) {
	  chrnum22 = 24;
	} else {
	  chrnum22 = atoi (per2->read2.chromosome+3);
	}

	instsize1 = chrnum12 - chrnum11;
	instsize2 = chrnum22 - chrnum21;
  }
  else {
    instsize1 = per1->read1.end - per1->read2.start;
    instsize2 = per2->read1.end - per2->read2.start;
  }

  if (instsize1 < 0)
    instsize1 = 0 - instsize1;
  if (instsize2 < 0)
    instsize2 = 0 - instsize2;

  return instsize2 - instsize1;
}

char *getMchrname (int chromosome)
{
  if (chromosome < 10) {
    return striappend ("mhs0", chromosome, 10);
  }
  else if (chromosome >= 10 && chromosome < 23) {
    return striappend ("mhs", chromosome, 10);
  }
  else if (chromosome == 23) {
    return strdup ("mhs0X");
  }
  else if (chromosome == 24) {
    return strdup ("mhs0Y");
  }
  else
  {
    return NULL;
  }
}

char *getSchrname (int chromosome)
{
  if (chromosome < 10) {
    return striappend ("shs0", chromosome, 10);
  }
  else if (chromosome >= 10 && chromosome < 23) {
    return striappend ("shs", chromosome, 10);
  }
  else if (chromosome == 23) {
    return strdup ("shs0X");
  }
  else if (chromosome == 24) {
    return strdup ("shs0Y");
  }
  else {
    return NULL;
  }
}

char *getHchrname (int chromosome)
{
  if (chromosome <= 22) {
    return striappend ("hs", chromosome, 10);
   }
  else if (chromosome == 23) {
    return strdup ("hsX");
   }
  else if (chromosome == 24) {
    return strdup ("hsY");
  }
  else {
    return NULL;
  }
}

int getScale (Array regions)
{
  int rlength = 0;
  SRegion_t *tmp;
  int i;

  tmp = arrayp(regions, 0, SRegion_t);
  if (tmp->chromosome >= 1 && tmp->chromosome <= 24) {
    for (i = 0; i < arrayMax(regions); i++) {
      tmp = arrayp(regions, i, SRegion_t);
      rlength += tmp->end - tmp->start;
    }
  }

  if (rlength > 500000000 || rlength == 0) {
    // scale 0
    return 1000000;
  }
  else if (rlength <= 500000000 && rlength > 50000000) {
    // scale 1
    return 100000;
  }
  else if (rlength <= 50000000 && rlength > 5000000) {
    // scale 2
    return 10000;
  }
  else if (rlength <= 5000000 && rlength > 500000) {
    // scale 3
    return 1000;
  }
  else if (rlength <= 500000 && rlength > 50000) {
    // scale 4
    return 100;
  }
  else if (rlength <= 50000 && rlength > 5000) {
    // scale 5
    return 10;
  }
  else {
    // scale 6
    return 1;
  }
}

int getChrnum (char *chrname)
{
  if (!strcmp (chrname, "X")) {
    return 23;
  }
  else if (!strcmp (chrname, "Y")) {
    return 24;
  }
  else {
    return atoi (chrname);
  }
}

// X as 23, Y as 24
int chrsize(int chromosome)
{
  switch (chromosome) {
  case 1:
    return 247249719;
    break;
  case 2:
    return 242951149;
    break;
  case 3:
    return 199501827;
    break;
  case 4:
    return 191273063;
    break;
  case 5:
    return 180857866;
    break;
  case 6:
    return 170899992;
    break;
  case 7:
    return 158821424;
    break;
  case 8:
    return 146274826;
    break;
  case 9:
    return 140273252;
    break;
  case 10:
    return 135374737;
    break;
  case 11:
    return 134452384;
    break;
  case 12:
    return 132349534;
    break;
  case 13:
    return 114142980;
    break;
  case 14:
    return 106368585;
    break;
  case 15:
    return 100338915;
    break;
  case 16:
    return 88827254;
    break;
  case 17:
    return 78774742;
    break;
  case 18:
    return 76117153;
    break;
  case 19:
    return 63811651;
    break;
  case 20:
    return 62435964;
    break;
  case 21:
    return 46944323;
    break;
  case 22:
    return 49691432;
    break;
  case 23:
    return 154913754;
    break;
  case 24:
    return 57772954;
    break;
  default:
    return 0;
    break;
  }
}

int putsn (int ntabs, const char *str)
{
  int i;
  for (i = 0; i < ntabs; i++)
    putchar ('\t');
  return puts (str);
}

// Postcondition: str is freed
char *striappend(char *str, int value, int base)
{
  char *tmp = (char *) malloc((ILEN) * sizeof(*tmp));
  char *new;

  tmp = itoa(value, tmp, base);
  new = (char *) malloc((strlen(tmp) + strlen(str) + 1) * sizeof(*new));

  strcpy(new, str);
  strcat(new, tmp);

  free(tmp);
  return new;
}

char *strappend(char *destination, const char *source)
{
  destination = (char *) realloc(destination,
		  (strlen(destination) + strlen(source) + 1)*sizeof(*destination));

  strcat(destination, source);
  return destination;
}

char *itoa (int value, char *str, int base)
{
  int i = 0;
  int sign;

  if ((sign = value) < 0)
    value *= -1;

  do {
    str[i] = value % 10 + '0';
    i++;
  } while ((value /= 10) > 0);

  if (sign < 0) {
    str[i] = '\0';
    i++;
  }

  str[i] = '\0';
  reverse(str);
  return str;
}

int antoi(const char *str, size_t n)
{
  char *tmp = (char *) malloc((n+1) * sizeof(*str));
  int i;

  tmp = strncpy(tmp, str, n);
  tmp[n] = '\0';
  i = atoi(tmp);
  free(tmp);
  return i;
}

void reverse(char *str)
{
  char c;
  int i, j;
  for (i = 0, j = strlen(str)-1; i < j; i++, j--) {
    c = str[j];
    str[j] = str[i];
    str[i] = c;
  }
}
