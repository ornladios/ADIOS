#include <stdio.h>
#include <stdlib.h>
#include <string.h>
    
#include "adios_query_hooks.h"
#include "public/adios_query.h"

#define ASSIGN_FNS(a,b) \
  (*t) [b].adios_query_init_method_fn          = adios_query_##a##_init_method; \
  (*t) [b].adios_query_free_method_fn          = adios_query_##a##_free_method; \
  (*t) [b].adios_query_estimate_method_fn      = adios_query_##a##_estimate_method;  	\
  (*t) [b].adios_query_get_selection_method_fn = adios_query_##a##_get_selection_method; \
  (*t) [b].adios_query_clean_method_fn         = adios_query_##a##_clean_method; 


void adios_query_hooks_init(struct adios_query_hooks_struct ** t, enum ADIOS_QUERY_TOOL tool )
{ 
	static int has_init_called[ADIOS_QUERY_TOOL_COUNT] = {0};
	if (tool < 0 || tool >= ADIOS_QUERY_TOOL_COUNT) {
		fprintf(stderr, "adios_query_hooks_init(): tool ID %d is an invalid tool identifier\n", (int)tool);
		exit(EXIT_FAILURE);
	}
	if (has_init_called[tool]) {
		return;
	}
	has_init_called[tool] = 1;

  fflush(stdout);
 
  *t = (struct adios_query_hooks_struct *) calloc (ADIOS_QUERY_TOOL_COUNT, sizeof (struct adios_query_hooks_struct));

  switch (tool) {
#ifdef ALACRITY
  case ADIOS_QUERY_TOOL_ALACRITY:
	  ASSIGN_FNS(alac, ADIOS_QUERY_TOOL_ALACRITY);
	  break;
#endif
#ifdef FASTBIT
  case ADIOS_QUERY_TOOL_FASTBIT:
	  ASSIGN_FNS(fastbit, ADIOS_QUERY_TOOL_FASTBIT);
	  break;
#endif
  default:
	  printf("unknown query tool type %d\n", (int)tool);
	  exit(EXIT_FAILURE);
  }
}

#undef ASSIGN_FNS


