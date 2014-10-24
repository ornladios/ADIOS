#ifndef ADIOS_QUERY_HOOKS_H
#define ADIOS_QUERY_HOOKS_H

#include "public/adios_query.h"
#include "public/adios_read.h"

#define FORWARD_DECLARE(a) \
  int     adios_query_##a##_can_evaluate(ADIOS_QUERY* q); \
  int64_t adios_query_##a##_estimate(ADIOS_QUERY* q, int timeStep);			\
  int     adios_query_##a##_evaluate(ADIOS_QUERY* q, int timeStep, uint64_t batchSize, ADIOS_SELECTION* outputBoundry, ADIOS_SELECTION** result); \
  int     adios_query_##a##_free(ADIOS_QUERY* q); \
  int     adios_query_##a##_finalize();

FORWARD_DECLARE(fastbit)
FORWARD_DECLARE(alac)

typedef int      (* ADIOS_QUERY_FREE_FN) (ADIOS_QUERY* q);
typedef int      (* ADIOS_QUERY_FINALIZE_FN) ();
typedef int      (* ADIOS_QUERY_EVALUATE_FN) (ADIOS_QUERY* q, int timeStep, uint64_t batchSize, ADIOS_SELECTION* o, ADIOS_SELECTION** result);
typedef int64_t  (* ADIOS_QUERY_ESTIMATE_FN) (ADIOS_QUERY* q, int timeStep);
typedef int  (* ADIOS_QUERY_CAN_EVALUATE_FN) (ADIOS_QUERY* q);

struct adios_query_hooks_struct
{
  const char *                  method_name;
  ADIOS_QUERY_EVALUATE_FN       adios_query_evaluate_fn;
  ADIOS_QUERY_FINALIZE_FN       adios_query_finalize_fn;
  ADIOS_QUERY_FREE_FN           adios_query_free_fn;
  ADIOS_QUERY_ESTIMATE_FN       adios_query_estimate_fn;
  ADIOS_QUERY_CAN_EVALUATE_FN   adios_query_can_evaluate_fn;
};
  
void adios_query_hooks_init (struct adios_query_hooks_struct ** t);

#undef FORWARD_DECLARE
#endif
