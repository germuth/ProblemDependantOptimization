/*   arg.c
**
**   This file contains functions for parsing the arguments from the
**   command line.
**
*/


#include <cstdio>
#include <cstring>
#include <cmath>
#include "cover.h"
#include "arg.h"
#include "tables.h"


struct varstruct {
  char *varname;
  char *format;
  void *varaddr;
} vars[] =
{
  {"k", "%d", &k},
  {"t", "%d", &t},
  {"m", "%d", &m},
  {"v", "%d", &v},
  {"l", "%d", &coverNumber}, /* lambda is coverNumber in the program! */
  {"lambda", "%d", &coverNumber}, /* lambda is coverNumber in the program! */
  {"b", "%d", &b},
  {"finalB", "%d", &finalB},
  {"FB", "%d", &finalB},
  {"bFinal", "%d", &finalB},
  {"BF", "%d", &finalB},
  {"bSearch", "%d", &bSearch},
  {"BS", "%d", &bSearch},
  {"endLimit", "%d", &endLimit},
  {"EL", "%d", &endLimit},
  {"pdo", "%d", &pdoFlag},
  {"pdoK", "%f", &pdoK},
  {"pdoJ", "%f", &pdoJ},
  {"pdoPrint", "%d", &pdoPrint},
  {"pdoPrintFreq", "%d", &pdoPrintFreq},
  {"pdoMaxJDF", "%d", &pdoMaxJDF},
  {"startFromFile", "%d", &startFromFileFlag},
  {"greedyStart", "%d", &greedyStartFlag},
  {"ontheFly", "%d", &onTheFly},
  {"OF", "%d", &onTheFly},
  {"pack", "%d", &pack},
  {"P", "%d", &pack},
  {"solX", "%d", &solX},
  {"SX", "%d", &solX},
  {"memoryLimit", "%lu", &memoryLimit},
  {"ML", "%lu", &memoryLimit},
  {"check", "%d", &check},
  {"SearchB", "%d", &searchB},
  {"SB", "%d", &searchB},
  {"SBFact", "%f", &SBFact},
  {"SBF", "%f", &SBFact},
  {"PRNGseed", "%u", &PRNGseed},
  {"PRNG", "%u", &PRNGseed},
};

#define OPT_COUNT (sizeof(vars) / sizeof(struct varstruct))


static void errmsg(char *msg)
{
  int i;
  static char format[40];

  fprintf(stderr, "Invalid option: %s\nValid options are:\n\n"
	  "       variable   format   default\n"
	  "       --------   ------   -------\n", msg);
  for(i = 0; i < OPT_COUNT; i++) {
    fprintf(stderr, "%15s   %6s   ", vars[i].varname, vars[i].format);
    strcpy(format, vars[i].format);
    strcat(format, "\n");
    if(strchr(vars[i].format, 'd'))
      fprintf(stderr, format, *((int *)vars[i].varaddr));
    else if(strchr(vars[i].format, 'f'))
      fprintf(stderr, "%4.2f\n", *((float *)vars[i].varaddr));
    else if(strchr(vars[i].format, 's'))
      fprintf(stderr, format, (char *)vars[i].varaddr);
    else if(strstr(vars[i].format, "lu"))
      fprintf(stderr, format, *((long unsigned *)vars[i].varaddr));
    /* add any types you wish */
  }
  coverError(SEE_ABOVE_ERROR);
}


#define TMP_BUF_SIZE 50

void parseArguments(int argc, char **argv)
{
  int i, j, eq;
  char name[TMP_BUF_SIZE];
  char value[TMP_BUF_SIZE];
  char *tmp;
  int found;
  int m_defined = 0, t_defined = 0;

  for(i = 1; i < argc; i++) {
    if((tmp = strchr(argv[i], '=')) == NULL)
      errmsg(argv[i]);
    eq = tmp - argv[i] + 1;
    strncpy(name, argv[i], fmin(eq - 1, TMP_BUF_SIZE - 1));
    name[(int)fmin(eq - 1, TMP_BUF_SIZE - 1)] = '\0';
    strncpy(value, argv[i] + eq, (int)fmin(strlen(argv[i]) - eq, TMP_BUF_SIZE - 1));
    value[(int)fmin(strlen(argv[i]) - eq, TMP_BUF_SIZE - 1)] = '\0';
    found = 0;
    for(j = 0; j < OPT_COUNT; j++)
      if(!strcmp(vars[j].varname, name)) {
	found++;
	if(sscanf(value, vars[j].format, vars[j].varaddr) != 1)
	  errmsg(argv[i]);
	if(!strcmp(name, "m") && !t_defined) { /* default: t=m if either one */
	  t = m;                               /* is not given */
	  m_defined = 1;
	}
	if(!strcmp(name, "t") && !m_defined) {
	  m = t;
	  t_defined = 1;
	}
	if(!strcmp(name, "PRNGseed") || !strcmp(name, "PRNG"))
	  setPRNGseed(PRNGseed);
	break;
      }
    if(!found)
      errmsg(argv[i]);
  }
}
