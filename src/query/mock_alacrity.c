#include "public/adios_query.h"
   
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "adios_read.h"



int64_t adios_query_alac_estimate_method(ADIOS_QUERY* q, int timeStep) 
{}
int64_t adios_query_alac_evaluate_method(ADIOS_QUERY* q, int timeStep, uint64_t _maxResult) 
{}
int  adios_query_alac_get_selection_method(ADIOS_QUERY* q, 
					   uint64_t batchSize, 
					   ADIOS_SELECTION* outputBoundry, 
					   ADIOS_SELECTION** result)
{}

int  adios_query_alac_free_method(ADIOS_QUERY* query) 
{}
void adios_query_alac_clean_method()
{}



void adios_query_alac_init_method() 
{
  const char* conffile = 0;
#ifdef DEBUG
  int msglvl = 200;
#else
  int msglvl = 0;
#endif

  printf("     [moock has initialized with msglvl = %d]\n", msglvl);
}
