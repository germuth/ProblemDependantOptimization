/*   cover.h
**
**   This file contains definitions common to all modules (or many of
**   them). This file should be included in *every* module.
**
*/

#ifndef _cover_
#define _cover_


#include <stdio.h>
#include <limits.h>
#include <cstdlib>
//#define maxv 40
#define maxv 70
/* binomial coefficients are tried to calculate up to binCoef[maxv][?]
 * (overflow is checked and the program is not halted
 */

//#define min(X,Y) ((X) < (Y) ? (X) : (Y))


/* initialize random number generator */
//Note that the random number generated is seeded with the current time which
//only has precision up to a second
//this means if two simulations are run within a second, they will run identically

//if windows
//#define randomize() srand((int)time(NULL))
#define randomize() srand((int)time(NULL))
#define setPRNGseed(X) srand(X)
#define rnd(num) (rand() % (num))
//#define rnd(num) 13
#define random01() ((double) rand() / 32767)
//if UNIX

//if windows
/* give the seed X for the random number generator */


/* gives an integer random number between 0 and num - 1, inclusive */

#define MAXRANDOM INT_MAX
/* gives a floating point random number in [0,1) */


typedef unsigned rankType;
/* this type should be able to contain binCoef[v][k] and binCoef[v][m] */

typedef unsigned short coveredType;
/* this type should be able to contain maxkSetCount. Keep it unsigned, if
   you want the definition of maxkSetCount to continue to work. */

#define maxkSetCount ((((unsigned long) 1) << (8 * sizeof(coveredType))) - 1)
/* maximum number of k-sets. If your type `char' is different from 8 bits
   or you have some other problems, replace this with the absolute number.
   See also `coveredType' above. */

typedef unsigned binCoefType;
/* this type should be able to contain binomial coefficients needed */

typedef unsigned char varietyType;
/* this type should be able to contain maxv + 1 (!) */

typedef int costType;
/* this type should be able to contain maximum cost */

typedef int costDType;
/* this type should be able to contain the difference between any two
   solutions */

extern int v, k, t, m, b; /* global parameters */
extern int pdoFlag;
extern int startFromFileFlag, greedyStartFlag;
extern float pdoK, pdoJ;
extern int pdoPrint, pdoPrintFreq, pdoMaxJDF;
extern int check;
extern int startB, finalB, bSearch;
extern costType finalCost;
extern char resultFileName[];
extern char logFileName[];
extern int testCount;
extern int bcount;
extern int endLimit;
extern int onTheFly;
extern int coverNumber;
extern int solX;
extern long unsigned memoryLimit;
extern int searchB;
extern float SBFact;
extern int pack;
extern int check;
extern int verbose;
extern unsigned int PRNGseed;

/* functions common to all modules */

extern int compareVarieties(varietyType *a, varietyType *b);
extern void coverError(int errnum);

/* error numbers */

#define BINCOEF_OVERFLOW       1
#define INTERNAL_OVERFLOW      2
#define INVALID_PARAMETERS     3
#define V_TOO_LARGE            4
#define MEM_ALLOC_ERROR        5
#define SEE_ABOVE_ERROR        6
#define B_TOO_LARGE            7
#define TOO_MUCH_SPACE         8
#define INTERNAL_ERROR         9
#define RANKTYPE_OVERFLOW     10
#define CHECK_SOLUTION_ERROR  11

#endif  /* #ifndef _cover_ */
