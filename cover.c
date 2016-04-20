/*   cover.c
**   This program tries to find covering designs using PDO.
**
**   Modified from Kari J. Nurmela at the Digital Systems Laboratory at
**   Helsinki University of Technology.
**   Modified By Aaron Germuth and Iliya Bluskov, University of Northern
**   British Columbia, 2016.
**
**   This program can be freely used whenever the following stipulations
**   are followed:
**
**   TODO write licence
*/

#include <stdio.h>
#include <stdlib.h>
//changed to help OSX compilation?
//#include <malloc/malloc.h>
#include <malloc.h>
#include <sys/time.h>
//not available on windows
//#include <sys/resource.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "cover.h"
#include "bincoef.h"
#include "tables.h"
#include "exp.h"
#include "pdo.h"

int v=6, k=5, t=3, m=4, b=4;
int pdoFlag = 1;
int pack = 0;
int check = 1;
int startB;
int finalB = 0;
int bSearch = 0;
int finalCost = 0;

int endLimit = 0;
int startFromFileFlag = 0;
int greedyStartFlag = 0;

float pdoK = 10;
int pdoJ = 2;
int pdoPrint = 2;
int pdoPrintFreq = 250;
int pdoMaxJDF = 5000000;

int onTheFly = 0;
int coverNumber = 1;
int solX = 0;
long unsigned memoryLimit = 0;

int testCount = 1;
int searchB = 0;
float SBFact = 0.95;

int verbose = 2;
unsigned int PRNGseed;

/*
** printSolution prints the current solution.
**
*/
void printSolution(FILE *fp) {
  int i,j;
  varietyType set[maxv + 1];
  varietyType *vptr;

  for(j = 0; j < b; j++) {
    printSubset(fp, kset[j], k);
    fprintf(fp, "\n");
  }
}


/*
** coverError prints an error message and exits the program
**
*/
//errno is a reserved function on windows environments....
void coverError(int errnum)
{
  static char *errmsg[] = {
    "No Error",
    "Binomial coefficient overflow",
    "Internal overflow",
    "Invalid parameters",
    "Parameter v too large, not enough space reserved",
    "Memory allocation error",
    "(See above)",
    "Parameter b is larger than maxkSetCount",
    "Space demands exceed the limit given by MemoryLimit option",
    "Internal error. Make a bug report.",
    "RankType is too small to contain the binomial coefficients needed.",
    "Cost change calculation gives wrong result.",
  };

  fprintf(stderr, "\n\nERROR: %s\n\n", errmsg[errnum]);
  exit(errnum);
}


/*
** newBAfterSuccess gradually decreases b (for option SearchB=1)
**
*/

static int newBAfterSuccess(int oldB)
{
  int lb;

  lb = (int) (SBFact * oldB + 0.5);
  if(lb == oldB)
    lb--;
  return lb;
}


/*
** searchBContinues calculates the new b, if needed (for option
** SearchB = 1)
**
*/


static int newSplitB(int b, int hi, int lo, int found)
{
  int hlp;

  if(hi - lo <= 1)
    return 0;
  bIs(lo + (hi - lo + 1) / 2);
  return 1;
}

static int searchBContinues(int found, int *hiB, int *loB)
{
  if(!searchB)
    return 0;
  if(*loB == -1) /* no failed yet */
    if(found) {
      *hiB = b;
      bIs(newBAfterSuccess(b));
      return 1;
    }
    else {
      *loB = b;
      return newSplitB(b, *hiB, *loB, found);
    }
  else {
    if(found)
      *hiB = b;
    else
      *loB = b;
    return newSplitB(b, *hiB, *loB, found);
  }
}


/*
** printParams prints the parameters given to the program.
**
*/
void printParams(FILE *fp)
{
  fprintf(fp, "Design parameters:\n"
	  "------------------\n"
	  "t - (v,m,k,l) = %d - (%d,%d,%d,%d)\nb = %d\n\n",
	  t, v, m, k, coverNumber, b);
  fprintf(fp, "Optimization parameters:\n"
	  "------------------------\n");
  fprintf(fp, "EndLimit      = %d\n"
    "FinalB        = %d\n"
    "PDO-Search    = %d\n"
    "PDO-K         = %d\n"
    "PDO-J         = %d\n"
    "PDO-Print     = %d\n"
    "PDO-Print-Freq= %d\n"
    "PDO-Max-JDF   = %d\n"
    "greedyStart   = %d\n"
    "startFromFile = %d\n"
    "OntheFly      = %d\n"
    "Packing       = %d\n"
    "SolX          = %d\n"
    "MemoryLimit   = %lu\n"
    "check         = %d\n\n", endLimit, finalB, pdoFlag, pdoK, pdoJ, pdoPrint,
        pdoPrintFreq, pdoMaxJDF, greedyStartFlag, startFromFileFlag, onTheFly,
        pack, solX, memoryLimit, check);
  fflush(fp);
}


/*
** compareVarieties is needed for qsort int randomNeighbor()
**
*/
int compareVarieties(varietyType *a, varietyType *b) {
  if(*a < *b)
    return -1;
  else {
    if(*a > *b)
      return 1;
    else
      return 0;
  }
}

int main(int argc, char **argv) {
    costType retVal;
    int j, i, count, bcounter;
    int iterSum;
    // struct rusage before, after;
    costType costSum, costSquareSum;
    float CPU, CPUsum;
    int solFound = 0;
    char *logName, *resultName;
    FILE *logFp, *resFp;
    int hiB = -1, loB = -1;
    float costSD = 0.;

    randomize();

    parseArguments(argc, argv);

    printf("\n" "cover-PDO v1.0 - find covering designs using PDO\n"
        "============================================================\n\n");
    calculateBinCoefs();   /* compute tables for binomial coefficients */
    calculate_exps();      /* and approximate exponentiation           */

    //neighbour and cover tables
    computeTables(t, k, m, v);       /* compute tables for this design */

    time_t rawtime;
    struct tm * timeinfo;
    time ( &rawtime );
    timeinfo = localtime ( &rawtime );
    printf ( "Started at %s\n", asctime (timeinfo) );

    startB = b;
    while((!pack && (b >= finalB) || (pack && b <= finalB))){
        if(pack){
            printf("Searching for a (%d,%d,%d,%d,%d) packing in %d blocks. (v,k,m,t,lambda)\n",
                v,k,m,t,coverNumber,b);
        }else{
            printf("Searching for a (%d,%d,%d,%d,%d) covering in %d blocks. (v,k,m,t,lambda)\n",
                v,k,m,t,coverNumber,b);
        }
        if(pack){
            asprintf(&logName, "./solutions/P(%d,%d,%d,%d,%d) - %d.log", v,k,m,t,coverNumber,b);
        }else{
            asprintf(&logName, "./solutions/C(%d,%d,%d,%d,%d) - %d.log", v,k,m,t,coverNumber,b);
        }
        logFp = fopen(logName, "w");
        if(!logFp) {
            fprintf(stderr, "Can't open log file %s.\n", logName);
            coverError(SEE_ABOVE_ERROR);
        }
        if(verbose && !pdoFlag){
            printParams(stdout);
        }
        printParams(logFp);
        fprintf(logFp, "\nRuns:\n-----\n");

        finalCost = pdo();
        if(finalCost <= endLimit) {
            solFound = 1;
            sortSolution();
        }
            //getrusage(RUSAGE_SELF, &after);
            // CPU = after.ru_utime.tv_sec + after.ru_utime.tv_usec / 1000000.0 -
            // (before.ru_utime.tv_sec + before.ru_utime.tv_usec / 1000000.0);
            if(verbose){
               printf("Result:\n" "-------\n"
                "EndLimit      = %d\n\n", endLimit);
            }
            fprintf(logFp, "cost          = %d\n", finalCost);
            fflush(logFp);
            if(verbose){
                if(finalCost <= endLimit) {
                    printf("Solution:\n" "---------\n");
                } else {
                    printf("EndLimit was not reached.\n\n");
                    if(verbose >= 2) {
                        printf("Inadequate solution:\n" "--------------------\n");
                        printSolution(stdout);
                    }
                }
            }
            if(finalCost <= endLimit) {
                if(verbose) {
                    printSolution(stdout);
                }
                if(pack){
                    asprintf(&resultName, "./solutions/P(%d,%d,%d,%d,%d) - %d.res", v,k,m,t,coverNumber,b);
                }else{
                    asprintf(&resultName, "./solutions/C(%d,%d,%d,%d,%d) - %d.res", v,k,m,t,coverNumber,b);
                }
                resFp = fopen(resultName, "w");
                if(!resFp) {
                    fprintf(stderr, "Can't open file %s.\n", resultName);
                    coverError(SEE_ABOVE_ERROR);
                }
                printSolution(resFp);
                fclose(resFp);
            }
            if(check) {
                int checks = checkSolution();
                if(checks != finalCost){
                    printf("PDO reported finalCost was %d. Manual Check returned %d\n", finalCost, checks);
                    coverError(CHECK_SOLUTION_ERROR);
                } else if(verbose){
                    printf("Final cost check OK.\n\n");
                }
            } else if(verbose){
                printf("\n");
            }

            if(bSearch){
                if(pack){
                    b++;
                }else{
                    b--;
                }
            }else{
                break;
            }
            fflush(stdout);
            fflush(logFp);
    }

    printf("Done \n");
    // fprintf(logFp, "\nStatistics:\n-----------\n"
    //     "CPU-time      = %.2f\n", CPUsum);

    freeTables();
    return !solFound; /* returns 0 if a solution was found */
}
