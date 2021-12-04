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
#include "setoper.h"
#include "solcheck.h"
#include "arg.h"
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
float pdoJ = 3;
int pdoPrint = 2;
int pdoPrintFreq = 500;
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
** printParams prints the parameters given to the program.
**
*/
void printParams(FILE *fp){
  fprintf(fp, "Design parameters:\n"
	  "------------------\n"
	  "t - (v,m,k,l) = %d - (%d,%d,%d,%d)\nb = %d\n\n",
	  t, v, m, k, coverNumber, b);
  fprintf(fp, "Optimization parameters:\n"
	  "------------------------\n");
  fprintf(fp, "EndLimit      = %d\n"
    "bSearch       = %d\n"
    "FinalB        = %d\n"
    "PDO-Search    = %d\n"
    "PDO-K         = %f\n"
    "PDO-J         = %f\n"
    "PDO-Print     = %d\n"
    "PDO-Print-Freq= %d\n"
    "PDO-Max-JDF   = %d\n"
    "greedyStart   = %d\n"
    "startFromFile = %d\n"
    "OntheFly      = %d\n"
    "Packing       = %d\n"
    "SolX          = %d\n"
    "MemoryLimit   = %lu\n"
    "check         = %d\n\n", endLimit, bSearch, finalB, pdoFlag, pdoK, pdoJ, pdoPrint,
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

char *logName, *resultName;
FILE *logFp, *resFp;

void printProblemParams(){
    if(pack){
        printf("Searching for a (%d,%d,%d,%d,%d) packing in %d blocks. (v,k,m,t,lambda)\n\n",
        v,k,m,t,coverNumber,b);
    }else{
        printf("Searching for a (%d,%d,%d,%d,%d) covering in %d blocks. (v,k,m,t,lambda)\n\n",
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
}

void printProblemOutput(){
    if(verbose){
        printf("Result:\n" "-------\n"
        "EndLimit      = %d\n\n", endLimit);
    }
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
}

int main(int argc, char **argv) {
    costType retVal;
    int j, i, count, bcounter;
    int solFound = 0;

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

    if(finalB == 0 && pack){
        finalB = 9999999; //finalB needs to be arbitrarly large with packings
    }

    startB = b;
    while((!pack && (b >= finalB)) || (pack && (b <= finalB))){
        //Creating Files and Output Text
        printProblemParams();

        finalCost = pdo();
        if(finalCost <= endLimit) {
            solFound = 1;
            sortSolution();
        }
            //getrusage(RUSAGE_SELF, &after);
            // CPU = after.ru_utime.tv_sec + after.ru_utime.tv_usec / 1000000.0 -
            // (before.ru_utime.tv_sec + before.ru_utime.tv_usec / 1000000.0);
        printProblemOutput();

        if(check) { //double check that solution is valid
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

        //adjust block count if requested
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

    freeTables();
    printf("Program Terminated Successfully \n");

    return !solFound; /* returns 0 if a solution was found */
}

