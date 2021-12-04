/*   tables.c
**
**   This file contains the functions needed when computing and handling
**   the neighborhood and coverings data tables.
**
*/

#include <cstdlib>
#include <cstring>
#include <cmath>
#include "cover.h"
#include "bincoef.h"
#include "setoper.h"
#include "tables.h"
#include <algorithm>

rankType *kset = NULL; //this is the solution
int neighborLen, coverLen, coveredLen;
rankType *neighbors, *coverings;
coveredType *covered;
costType *costs;
costDType *costds;


/*
** allocateMemory() allocates the memory for the tables with given v,k,m,t
** and calculates neighborLen, coverLen, coveredLen
**
*/

void allocateMemory(void)
{
  unsigned elemCountNeighbors, elemCountCoverings, elemCountCovered;
  unsigned long tmp;
  int i;

  if(verbose)
    printf("Memory allocation:\n"
	   "------------------\n");

  neighborLen = k * (v - k);

  if(overflowBinCoef(v, k) || overflowBinCoef(v, m))
    coverError(BINCOEF_OVERFLOW);
  if(binCoef[v][k] != (rankType) binCoef[v][k] ||
     binCoef[v][m] != (rankType) binCoef[v][m])
    coverError(RANKTYPE_OVERFLOW);
  if(v > maxv)
    coverError(V_TOO_LARGE);
  if(t > k || t > m || k > v || m > v || b <= 0 || t <= 0)
    coverError(INVALID_PARAMETERS);

  tmp = coverLen = 0;
  for(i = 0; i <= fmin(k - t, m - t); i++) {
    if(overflowBinCoef(k, t + i) || overflowBinCoef(v - k, m - t - i))
      coverError(BINCOEF_OVERFLOW);
    coverLen += binCoef[k][t + i] * binCoef[v - k][m - t - i];
    if(coverLen < tmp)
      coverError(INTERNAL_OVERFLOW);
    tmp = (unsigned long) coverLen;
  }
  coverLen++; /* sentinels after cover sets of each k-set */
  if(coverLen < tmp)
    coverError(INTERNAL_OVERFLOW);

  if(!onTheFly) {
    elemCountNeighbors = neighborLen * binCoef[v][k];
    if(abs(elemCountNeighbors -
	   (unsigned)((float) neighborLen * (float) binCoef[v][k])) > 1000)
      coverError(INTERNAL_OVERFLOW);
    if(verbose)
      printf("neighbors:%11u   elems\n", elemCountNeighbors);
    elemCountCoverings = coverLen * binCoef[v][k];
    if(abs(elemCountCoverings -
	   (unsigned)((float) coverLen * (float) binCoef[v][k])) > 1000)
      coverError(INTERNAL_OVERFLOW);
    if(verbose)
      printf("coverings:%11u   elems\n", elemCountCoverings);
  }
  else {
    elemCountNeighbors = 0;
    elemCountCoverings = coverLen * 2;
    if(elemCountCoverings <= coverLen)
      coverError(INTERNAL_OVERFLOW);
    if(verbose)
      printf("coverings:%11u   elems\n", elemCountCoverings);
  }

  elemCountCovered = coveredLen = binCoef[v][m];
  if(verbose)
    printf("covered:  %11u   elems\n\n", elemCountCovered);

  /* are the space demands too much? */
  if(memoryLimit > 0 &&
     ((float) elemCountNeighbors * sizeof(rankType) +
      (float) elemCountCoverings * sizeof(rankType) +
      (float) elemCountCovered * sizeof(coveredType)
      > memoryLimit))
    coverError(TOO_MUCH_SPACE);
  /* if not, try to get the memory */
  else {
    if(elemCountNeighbors)
      neighbors = (rankType *) calloc(elemCountNeighbors, sizeof(rankType));
    coverings = (rankType *) calloc(elemCountCoverings, sizeof(rankType));
    covered = (coveredType *) calloc(elemCountCovered, sizeof(coveredType));
  }

  /* was the memory allocation OK? */
  if(!covered || !coverings || (!neighbors && !onTheFly))
    coverError(MEM_ALLOC_ERROR);
}


/*
** calculateNeighbors computes the neighbor ranks for each rank of a k-set.
**
*/
void calculateNeighbors(void) {
  rankType r;
  varietyType subset[maxv + 1], csubset[maxv + 1];
  varietyType subsubset[maxv + 1], subcsubset[maxv + 1], mergeset[maxv + 1];
  varietyType *ssptr, *scptr, *mptr;
  rankType *nptr;
  int i;

  nptr = neighbors;
  getFirstSubset(subset, k);
  for(r = 0; r < (rankType) binCoef[v][k]; r++) {         //for each v choose k - choose a block
    makeComplement(subset, csubset, v);
    getFirstSubset(subsubset, k - 1); /* optimoinnin varaa */
    do {                                                  //for each k choose k-1 - choose which points to keep
      getFirstSubset(subcsubset, 1); /* optimoinnin varaa */
      do {                                                //for each v-k choose 1 - choose what NEW point to add
        ssptr = subsubset;
        scptr = subcsubset;
        mptr = mergeset;
        //I believe these are sentinel values
        subsubset[k - 1] = (varietyType) k;
        subcsubset[1] = (varietyType) v - k;
        //need to add keepers (subsubset) and the new point (subcsubset) together into one set
        //but we want them in order, so merge them
        for(i = 0; i < k; i++){
          if(subset[(int) *ssptr] < csubset[(int) *scptr]){
            *mptr++ = subset[(int) *ssptr++];
          } else {
            *mptr++ = csubset[(int) *scptr++];
          }
        }
        subsubset[k - 1] = (varietyType) maxv + 1; /* sentinel */
        subcsubset[1] = (varietyType) maxv + 1; /* sentinel */
        *mptr = maxv + 1; /* sentinel */
        *nptr++ = rankSubset(mergeset, k);
      } while(getNextSubset(subcsubset, 1, v - k));
    } while(getNextSubset(subsubset, k - 1, k));
    getNextSubset(subset, k, v);
  }
}


/*
** The ranks of the m-sets covered by a certain k-set are stored in
** increasing order in the table "coverings". The algorithm doesn't
** produce them in this order, so they must be sorted. We use a stan-
** dard qsort, and compareRanks() is the compare function provided for
** the sort.
**
*/

int compareRanks(rankType *a, rankType *b)
{
  if(*a < *b)
    return -1;
  else {
    if(*a > *b)
      return 1;
    else
      return 0;
  }
}


/*
** calculateOneCoverings() is the most complex function in this file. It
** computes the ranks of the m-sets covered by a given k-set.
** `subset' is the k-set being processed. `csubset' is
** the complement of this k-set. `ti' is the index that goes from t to
** min(k,m), that is there must be at least t varieties in the inter-
** section of the k-set and the m-set. `subsubset' is a subset of
** `subset'. `subcsubset' is subset of csubset. When `subset', `ti',
** `subsubset' and `subcsubset' get all the possible values, all that
** remains to do is to merge `subsubset' and `subcsubset' to form the
** m-set that is covered by `subset' (k-set) in `ti' varieties. The
** ranks of covered m-sets are stored in `buf'.
**
*/
void calculateOneCovering(rankType kRank, rankType *buf)
{
  static varietyType subset[maxv + 1], csubset[maxv + 1];
  static varietyType subsubset[maxv + 1];
  static varietyType subcsubset[maxv + 1], mergeset[maxv + 1];
  static varietyType *ssptr, *scptr, *mptr;
  static rankType *coverptr;
  static int i, ti;

  coverptr = buf;
  unrankSubset(kRank, subset, k);                 //subset is curr block
  subset[k] = maxv + 1; /* sentinel */
  makeComplement(subset, csubset, v);             //csubset is points not in curr block
  for(ti = t; ti <= fmin(k, m); ti++) {
    getFirstSubset(subsubset, ti);                //t-set that this block covers
    do {
      getFirstSubset(subcsubset, m - ti);         //other points not in block we can add to t-set to make m-set
      do {
        ssptr = subsubset;
        scptr = subcsubset;
        mptr = mergeset;
        subsubset[ti] = (varietyType) k;
        subcsubset[m - ti] = (varietyType) v - k;
        for(i = 0; i < m; i++){
          if(subset[(int) *ssptr] < csubset[(int) *scptr]) {
            *mptr++ = subset[(int) *ssptr++];
          } else {
            *mptr++ = csubset[(int) *scptr++];
          }
        }
        subsubset[ti] = (varietyType) (maxv + 1); /* sentinel */
        subcsubset[m - ti] = (varietyType) (maxv + 1); /* sentinel */
        *mptr = (varietyType) (maxv + 1); /* sentinel */
        *coverptr++ = rankSubset(mergeset, m);
      } while(getNextSubset(subcsubset, m - ti, v - k));
    } while(getNextSubset(subsubset, ti, k));
  }
  *coverptr = binCoef[v][m]; /* sentinel */
  std::sort(buf, kset+coverLen - 1);
}


/*
** calculateCoverings() computes the ranks of m-sets covered by all
** possible k-sets into table `coverings'.
**
*/
void calculateCoverings(void) {
  rankType r;

  for(r = 0; r < (rankType) binCoef[v][k]; r++) //for each possible block
    calculateOneCovering(r, coverings + ((int) r * coverLen));
}


/*
** freeTables() frees all dynamically allocated tables so that
** computeTables() can be called again.
**
*/

void freeTables(void) {
  if(!onTheFly) {
    free((void *) neighbors);
  }
  free((void *) coverings);
  free((void *) covered);
  if(kset)
    free((void *) kset);
}


/*
** computeTables() allocates memory for the tables and calculates them.
** It also sets the global variables t, k, m, v.
*/

void computeTables(int tl, int kl, int ml, int vl)
{
  t = tl;
  k = kl;
  m = ml;
  v = vl;
  allocateMemory();
  if(!onTheFly) {
    calculateNeighbors();
    calculateCoverings();
  }
}


/*
** While t,k,m,l are set by computeTables(), b is set separately, because
** it doesn't have an effect on the tables except table `kset'. We can do
** several simulated annealing processes with different values of b
** without computing the tables again.
**
*/
void bIs(int bl)
{
  b = bl;        /* b is the number of k-sets */
  if(b > maxkSetCount)
    coverError(B_TOO_LARGE);
  if(kset) {
    if(!(kset = (rankType *) realloc((char *) kset, b * sizeof(rankType))) ||
       !(costs = (costType *) realloc((char *) costs, (b + 1) * sizeof(costType))) ||
       !(costds = (costDType *) realloc((char *) costds, (b + 1) * sizeof(costDType))))
      coverError(MEM_ALLOC_ERROR);
  }
  else
    if(!(kset = (rankType *) malloc(b * sizeof(rankType))) ||
       !(costs = (costType *) malloc((b + 1) * sizeof(costType))) ||
       !(costds = (costDType *) malloc((b + 1) * sizeof(costDType))))
      coverError(MEM_ALLOC_ERROR);
}


/*
** `sortSolution()' sorts the current solution in the lexicographical order.
**
*/

void sortSolution(void)
{
  std::sort(kset, kset+b);
}
