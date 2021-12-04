#include <math.h>
#include <cstdlib>
#include <algorithm>
#include "cover.h"
#include "bincoef.h"
#include "tables.h"
#include "exp.h"
#include "pdo.h"
#include "setoper.h"
//always accept sideways moves
//moves at the same level are counted together towards the number of moves needed to find better solution
//initialize jdf to 10*neighbourhood-size
//update jdf = ceil[ n*jdf + newestAmountOFMoves / n+1 ]
//try 2*jdf amount of moves

//the jumpdownfunction used in the rest of the program is jdf (integers)
//however, in order ot update the jdf, the floating point values are always saved
//and used for all calculations. After a calculation the value is ceiling-ed
//to update the jdf
//this prevents accumulating ceilings from slowly raising the jump down function
int *jdf;
double *jdf_precise;
int *jdc;
int currCost;
int counter;
int maxCost = 0;
int printCount = 0;
rankType *oldKSet;
FILE *fp;

/*
** `setNumber' is the index to the table `kset'. `setNumber' indicates
** the index of the k-set in the proposed change. The current solution
** can be changed to the proposed next solution by assigning
** kset[setNumber] = nextS. `stored', `storedPtr', `currSto' and
**`nextSto' are for on-the-fly annealing. `costs[x]' holds the difference
** of costs associated with a m-set covered x times and x+1 times.
*/
static int setNumber;
static rankType nextS;
static rankType stored[2];
static int currSto, nextSto;
rankType *storedPtr[2];

/*
** Populates costs array
** Contains the punishment for how many times it was covered
** For example costs[0] is how much cost is added to a solution if you cover
** a single m-set 0 times
** costs[3] is if you cover it 3 times
** if lambda is 1 we would expect [1,0,0,0,0,0,0,0]
**  for a covering design (there is no cost for over-covering)
*/
void calculateCosts(void) {
  int i;

  if(pack) /* packing design */
    for(i = 0; i <= b; i++)
      if(i < coverNumber)
	costs[i] = (costType) 0;
      else
	costs[i] = (costDType) (i - coverNumber);
  else     /* covering design */
    for(i = 0; i <= b; i++)
      if(i < coverNumber)
	costs[i] = (costType) (coverNumber - i);
      else
	costs[i] = 0;
  for(i = 0; i < b; i++)
    costds[i] = costs[i] - costs[i + 1];
}

/*
** randomNeighbor() computes the rank of a random neigbor of the k-set with
** rank `curr'.
*/
rankType randomNeighbor(rankType curr)
{
  varietyType subset[maxv + 1];
  varietyType csubset[maxv + 1];

  unrankSubset(curr, subset, k);
  makeComplement(subset, csubset, v);
  subset[rnd(k)] = csubset[rnd(v-k)];
  std::sort(subset, subset+k);
  return rankSubset(subset, k);
}

typedef struct {
    rankType first;
    char second;
} pairdata_t;

bool comparePair(const pairdata_t& a, const pairdata_t& b)
{
  return (a.first < b.first);
}

/*
** computeNeighbor() calculates the cost difference between the current
** solution and a random neighbor of the current solution. It employs
** the sentinels at the end of the covered sets of each k-set. (That's
** why they were put there in calculateNeighbors().)
**
*/

static costType computeNeighbor(void)
{
  costType costDelta = 0;

  int i;
  rankType currS;
  rankType *currPtr, *nextPtr;
  const long limit = (coverLen-1)*2-1;

  setNumber = rnd(b);

  currS = kset[setNumber];
  nextS = onTheFly ? randomNeighbor(currS) :
    neighbors[currS * neighborLen + rnd(neighborLen)];

  if(onTheFly) {
    currPtr = NULL;
    for(i = 0; i < 2; i++)
      if(currS == stored[i]) {
          currPtr = storedPtr[i];
          currSto = i;
      }
    nextPtr = NULL;
    for(i = 0; i < 2; i++)
      if(nextS == stored[i]) {
          nextPtr = storedPtr[i];
          nextSto = i;
      }
    for(i = 0; !currPtr; i++)
      if(nextSto != i) {
          storedPtr[i] = currPtr = coverings + i * coverLen;
          currSto = i;
          stored[i] = currS;
          calculateOneCovering(currS, currPtr);
      }
    for(i = 0; !nextPtr; i++)
      if(currSto != i) {
          storedPtr[i] = nextPtr = coverings + i * coverLen;
          nextSto = i;
          stored[i] = nextS;
          calculateOneCovering(nextS, nextPtr);
      }
  }
  else {
    currPtr = coverings + currS * coverLen;
    nextPtr = coverings + nextS * coverLen;
  }

  for(i = 0; i < (coverLen - 1) * 2; i++)
    if(*currPtr == *nextPtr) {
      currPtr++;
      nextPtr++;
      i++;
    }
    else if(*currPtr < *nextPtr) {
      costDelta += costds[covered[*currPtr++] - 1];
    }
    else
      costDelta -= costds[covered[*nextPtr++]];

  return costDelta;
}


/*
** acceptNeighbor() changes the current solution to the latest solution
** computed.
*/

static void acceptNeighbor(void)
{
  int i;
  rankType currS;
  rankType *coveringsPtrA;
  rankType *coveringsPtrB;

  currS = kset[setNumber];
  if(onTheFly)
    coveringsPtrA = coverings + currSto * coverLen;
  else
    coveringsPtrA = coverings + currS * coverLen;

  if(onTheFly)
    coveringsPtrB = coverings + nextSto * coverLen;
  else
    coveringsPtrB = coverings + nextS * coverLen;
  #pragma clang loop unroll_count(100) vectorize_width(100) interleave_count(100)
  for(i = 0; i < coverLen - 1; i++) {
    covered[coveringsPtrA[i]]--;
    covered[coveringsPtrB[i]]++;
  }
  kset[setNumber] = nextS;
}

rankType getNextBlockFromFile(){
    int j;
    //TODO add some sort of verification here and print error if file sucks
    varietyType *block = (varietyType *)calloc(k+1, sizeof(int));
    for(j = 0; j < k; j++){
        int point;
        fscanf(fp, "%d", &point);
        block[j] = point;
    }
    block[k] = maxv + 1;
    // printf("%d %d %d %d %d %d\n", block[0],block[1],block[2],block[3],block[4], block[5]);
    return rankSubset(block, k);
    // printf("%d\n", kset[i]);
}

rankType getNextGreedyBlock(){
    int i2, j;
    //look at first m-set not covered
    //construct partial block that covers it
    //if block has leftover spaces, attempt to cover next m-set

    int idx = 0;
    varietyType *block = (varietyType *)calloc(k + 1, sizeof(int));

    for(j = 0; j < binCoef[v][m]; j++){
        if(k < idx+t){
            //not enough room to cover another one
            break;
        }
        if(!covered[j]){
            //this m-set not covered
            //take first t elements of m-set, and add to block
            varietyType *mset = (varietyType *)calloc(m + 1, sizeof(int));
            unrankSubset(covered[j], mset, m);
            for(i2 = 0; i2 < t; i2++){
                block[i2] = mset[i2 + idx];
            }
            idx += t;
        }
    }
    //fill in block with other stuff
    varietyType *com = (varietyType *)calloc(v - idx, sizeof(int));
    makeComplement(block, com, v);
    for(; idx < k; idx++){
        block[idx] = com[rnd(v-idx)];//TODO hope it doesn't hit same number twice
    }

    std::sort(block, block+k);
    //block done
    return rankSubset(block, k);
}

void printProgress(void){
    if(pdoPrint == 0){
        return;
    }
    int indices[] = {0,1,2,3,4,5,6,7,8,9};
    if(pdoPrint == 2){
        //print first 10 costs that have been found rather than just first 10
        int i = 0;
        int curr = 0;
        for(; i < maxCost && curr < 10; i++){
            if(jdc[i] > 1){
                indices[curr] = i;
                curr++;
            }
        }
        if(curr != 10){
            std::sort(indices, indices+10);
        }
    }else if(pdoPrint == 3){ //show 10 levels above current
        int i = 0;
        for(; i < 10; i++){
            indices[i] = currCost + i;
        }
    }
    printf("                                                                                                             ");
    printf("\rlevel\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n",
    indices[0], indices[1], indices[2], indices[3],
    indices[4], indices[5], indices[6], indices[7],
    indices[8], indices[9]);
    printf("                                                                                                             ");
    printf("\rjdc\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n",
    jdc[indices[0]], jdc[indices[1]], jdc[indices[2]], jdc[indices[3]],
    jdc[indices[4]], jdc[indices[5]], jdc[indices[6]], jdc[indices[7]],
    jdc[indices[8]], jdc[indices[9]]);
    printf("                                                                                                             ");
    printf("\rjdf\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d",
    jdf[indices[0]], jdf[indices[1]], jdf[indices[2]], jdf[indices[3]],
    jdf[indices[4]], jdf[indices[5]], jdf[indices[6]], jdf[indices[7]],
    jdf[indices[8]], jdf[indices[9]], currCost);
    printf("\033[F");
    printf("\033[F");
    // printf("\033[1A");
    // printf("\033[1A");
    // printf("\033[1A");
}

static void updateJDF(void){
    /*
    double ans = jdc[currCost] * jdf[currCost] + counter;
    ans /= (double)(jdc[currCost] + 1);
    */

    //Jan de Heer, more efficient way to compute
    //and avoids overflows from initial multiplication
    double ans = counter - jdf_precise[currCost];
    ans /= (double)(jdc[currCost] + 1);
    ans += (double)(jdf_precise[currCost]);

    jdf_precise[currCost] = ans;
    jdf[currCost] = (int)ceil(ans);

    //only try up to a million times?
    if(pdoMaxJDF && jdf[currCost] > pdoMaxJDF){
        jdf[currCost] = pdoMaxJDF;
    }
}

static costType initSolution(void) {
    int i, j, count, maxCount, l;
    costType initCost;
    coveredType *ptr;
    rankType *coveringsPtr;

    fp = fopen("startFile.txt", "r");
    for(i = 0; i < b; i++) {
        if(startFromFileFlag && !bSearch){
            kset[i] = getNextBlockFromFile();
        }else if(greedyStartFlag && !bSearch){
            kset[i] = getNextGreedyBlock();
        } else if( bSearch == 2 &&  //if taking from last results
                startB != b &&      //and this isn't the first one (no one to take frome)
                !(pack && i == b-1)){//and this isn't the last block (new packing has one new extra block)
            kset[i] = oldKSet[i]; //copy it over from previous solution
        } else {
            kset[i] = rnd(binCoef[v][k]); //generate randomly
        }

        if(onTheFly) {
            calculateOneCovering(kset[i], coverings);
            coveringsPtr = coverings;
        } else {
            coveringsPtr = coverings + (int) kset[i] * coverLen;
        }
        for(j = 0; j < coverLen - 1; j++) {
            covered[coveringsPtr[j]]++;
        }
    }
    if(bSearch == 2 && (startB != b) && oldKSet){
        //TODO free is crashing program for some reason
        //crashes occasionally with v=17 k=6 m=4 t=4 b=105 bFinal=114 bSearch=2 pack=1
        if(!pack){
            free((void *) oldKSet);
        }
    }
    fclose(fp);
    for(i = 0, initCost = (costType) 0, ptr = covered; i < coveredLen; i++, ptr++)
        initCost += costs[*ptr];

    return initCost;
}


costType pdo() {
    int i;

    if(b > maxkSetCount){
        coverError(B_TOO_LARGE);
    }
    //blank covering data
    for(i = 0; i < binCoef[v][m]; i++){
        covered[i] = 0;
    }

    //need to allocate kset, costs, costds
    if(bSearch == 2 && startB != b){ //allocate new storage if we are copying over
        oldKSet = kset;
        // printf("%d %d\n", oldKSet[0], kset[0]);
        if(!(kset = (rankType *) malloc(b * sizeof(rankType)))){
            coverError(MEM_ALLOC_ERROR);
        }
        // printf("%d %d\n", oldKSet[0], kset[0]);
    }else{ //reallocate if we are overwriting
        free((void *) kset);
        free((void *) costs);
        free((void *) costds);
        if(!(kset = (rankType *) malloc(b * sizeof(rankType))) ||
                !(costs = (costType *) malloc((b + 1) * sizeof(costType))) ||
                !(costds = (costDType *) malloc((b + 1) * sizeof(costDType))))
        // if(!(kset = (rankType *) realloc((char *) kset, b * sizeof(rankType))) ||
        //         !(costs = (costType *) realloc((char *) costs, (b + 1) * sizeof(costType))) ||
        //         !(costds = (costDType *) realloc((char *) costds, (b + 1) * sizeof(costDType))))
            coverError(MEM_ALLOC_ERROR);
    }
    //costs[lambda] holds costs for covering a single m-set lambda times
    maxCost = 0;
    calculateCosts();
    for(i = 0; i < b; i++){
        if(costs[i] > maxCost){
            maxCost = costs[i];
        }
    }
    //if the most we can be pushished for one m-set is max-cost, then the highest
    //possible cost is number of m-sets * max pushishment per m-set
    maxCost *= binCoef[v][m];

    //Jump Down Function: cost -> moves needed to find good neighbour
    jdf = (int *)calloc(maxCost, sizeof(int));
    jdf_precise = (double *)calloc(maxCost, sizeof(double));
    for(i = 0; i < maxCost; i++){
        jdf_precise[i] = pdoK * b * neighborLen;
        jdf[i] = (int)jdf_precise[i];
    }
    //Jump Down Count: jdc[i] = number of times we have jumped down from a solution at level i
    jdc = (int *)calloc(maxCost, sizeof(int));
    for(i = 0; i < maxCost; i++){
        jdc[i] = 1;
    }

    counter = 0;
    currCost = initSolution();

    //get random neighbour and its cost
    int costDifference = computeNeighbor();

    printProgress();

    while( currCost > endLimit){
        counter++;
        if( costDifference > 0) { //if( neighborCost > currCost){
            if( counter > pdoJ * jdf[currCost] ){
                //failed to find jump down from this level, accept bad neighbour anyway
                acceptNeighbor();
                counter = 0;
                currCost += costDifference;
            }
        }else if(costDifference == 0) { //if( neighborCost == currCost ){
            acceptNeighbor();
            //counter keeps counting
        }else { //found better solution
            updateJDF();
            jdc[currCost]++;//increment number of times we have jumped down
            acceptNeighbor();
            if(pdoPrint && !printCount--){
                printProgress();
                printCount = pdoPrintFreq;
            }
            counter = 0;
            currCost += costDifference;
            printf("currCost: %15ld\n\n", currCost);
            printf("\033[F");
            printf("\033[F");
        }
        //try random neighbour
        costDifference = computeNeighbor();

    }
    if(pdoPrint){
        printf("\n\n\n\n");
    }

    return currCost;//found solution
}
