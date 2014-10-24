#include <stdio.h>
#include <stdlib.h>
#include <string.h>
    
#include "adios_query_hooks.h"
#include "public/adios_query.h"

#define ASSIGN_FNS(a,b) \
  (*t) [b].method_name = #b; /* stringify the enum constant's name */ \
  (*t) [b].adios_query_free_fn          = adios_query_##a##_free; \
  (*t) [b].adios_query_estimate_fn      = adios_query_##a##_estimate; \
  (*t) [b].adios_query_can_evaluate_fn  = adios_query_##a##_can_evaluate;  \
  (*t) [b].adios_query_evaluate_fn      = adios_query_##a##_evaluate; \
  (*t) [b].adios_query_finalize_fn      = adios_query_##a##_finalize; 


void adios_query_hooks_init (struct adios_query_hooks_struct ** t)
{ 
    static int did_init = 0;
    if (did_init) {
        return;
    }
    did_init = 1;

    fflush(stdout);
 
    *t = (struct adios_query_hooks_struct *) calloc (ADIOS_QUERY_METHOD_COUNT, sizeof (struct adios_query_hooks_struct));

    int i=0;
    for (i=0; i<ADIOS_QUERY_METHOD_COUNT; i++) {
      (*t) [i].adios_query_free_fn = 0;
      (*t) [i].adios_query_estimate_fn = 0;
      (*t) [i].adios_query_can_evaluate_fn = 0;
      (*t) [i].adios_query_evaluate_fn = 0;
      (*t) [i].adios_query_finalize_fn = 0;
    }

#ifdef ALACRITY
    ASSIGN_FNS(alac, ADIOS_QUERY_METHOD_ALACRITY);
#endif
#ifdef FASTBIT
    ASSIGN_FNS(fastbit, ADIOS_QUERY_METHOD_FASTBIT);
#endif
}

#undef ASSIGN_FNS


