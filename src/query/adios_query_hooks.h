#ifndef ADIOS_QUERY_HOOKS_H
#define ADIOS_QUERY_HOOKS_H

#include "adios_query.h"
#include "adios_read.h"

#define FORWARD_DECLARE(a) \
  int64_t adios_query_##a##_estimate_method(ADIOS_QUERY* q);	\
  int     adios_query_##a##_get_selection_method(ADIOS_QUERY* q, uint64_t batchSize, ADIOS_SELECTION* outputBoundry, ADIOS_SELECTION** result); \
  int     adios_query_##a##_init_method ();					\
  int     adios_query_##a##_free_method (ADIOS_QUERY* q);			\
  int     adios_query_##a##_clean_method();


// default query tool is fastbit
FORWARD_DECLARE(fastbit)
FORWARD_DECLARE(alacrity)

typedef int      (* ADIOS_QUERY_FREE_METHOD_FN) (ADIOS_QUERY* q);
typedef int      (* ADIOS_QUERY_INIT_METHOD_FN) ();
typedef int      (* ADIOS_QUERY_CLEAN_METHOD_FN) ();
typedef int      (* ADIOS_QUERY_GET_SELECTION_METHOD_FN) (ADIOS_QUERY* q,  uint64_t batchSize, ADIOS_SELECTION* o, ADIOS_SELECTION** result);
typedef int64_t  (* ADIOS_QUERY_ESTIMATE_METHOD_FN) (ADIOS_QUERY* q);
//typedef int64_t  (* ADIOS_QUERY_EVALUATE_METHOD_FN) (ADIOS_QUERY* q, int timeStep, uint64_t maxResult);

struct adios_query_hooks_struct
{
  ADIOS_QUERY_INIT_METHOD_FN         adios_query_init_method_fn;
  //ADIOS_QUERY_EVALUATE_METHOD_FN     adios_query_evaluate_method_fn;
  ADIOS_QUERY_GET_SELECTION_METHOD_FN     adios_query_get_selection_method_fn;
  ADIOS_QUERY_CLEAN_METHOD_FN        adios_query_clean_method_fn;
  ADIOS_QUERY_FREE_METHOD_FN         adios_query_free_method_fn;
  ADIOS_QUERY_ESTIMATE_METHOD_FN     adios_query_estimate_method_fn;
};
  
void adios_query_hooks_init(struct adios_query_hooks_struct ** t);

#undef FORWARD_DECLARE
#endif
