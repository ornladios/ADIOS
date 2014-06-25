#include <stdio.h>
#include <stdlib.h>
#include <string.h>
    
#include "adios_query_hooks.h"
#include "adios_query.h"

#define ASSIGN_FNS(a,b) \
  (*t) [b].adios_query_init_method_fn          = adios_query_##a##_init_method; \
  (*t) [b].adios_query_free_method_fn          = adios_query_##a##_free_method; \
  (*t) [b].adios_query_estimate_method_fn      = adios_query_##a##_estimate_method;  	\
  (*t) [b].adios_query_get_selection_method_fn = adios_query_##a##_get_selection_method; \
  (*t) [b].adios_query_clean_method_fn         = adios_query_##a##_clean_method; 


void adios_query_hooks_init(struct adios_query_hooks_struct ** t, enum ADIOS_QUERY_TOOL tool )
{ 
  static int has_init_called = 0;

  if (has_init_called) {
      return;
  }
  
  fflush(stdout);

 
  *t = (struct adios_query_hooks_struct *) calloc (ADIOS_QUERY_TOOL_COUNT, sizeof (struct adios_query_hooks_struct));

  //#ifdef _USE_FASTBIT
  // default initiation

  /*if (tool ==  ADIOS_QUERY_TOOL_FASTBIT){
	  ASSIGN_FNS(fastbit, ADIOS_QUERY_TOOL_FASTBIT);

  }else */if ( tool == ADIOS_QUERY_TOOL_ALACRITY){
	  ASSIGN_FNS(alac, ADIOS_QUERY_TOOL_ALACRITY);
  }else {
	  printf("unknown query tool type \n");
	  exit(EXIT_FAILURE);
  }

  has_init_called = 1;

  //#endif

}

#undef ASSIGN_FNS


