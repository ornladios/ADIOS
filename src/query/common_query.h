#ifndef __COMMON_QUERY_H__
#define __COMMON_QUERY_H__

#ifdef __cplusplus
extern "C" {
#endif

#include "public/adios_query.h"
#include "core/futils.h"

#define NO_EVAL_BEFORE -1

// called from Read init only; it does NOT call methods' init
void common_query_init();

//void common_query_set_method(enum ADIOS_QUERY_METHOD method);

int common_query_is_method_available(enum ADIOS_QUERY_METHOD method);
void common_query_set_method (ADIOS_QUERY* q, enum ADIOS_QUERY_METHOD method);
int adios_get_actual_timestep(ADIOS_QUERY* q, int timeStep);

ADIOS_QUERY* common_query_create(ADIOS_FILE* f, 				 
				 ADIOS_SELECTION* queryBoundry,
				 const char* varName,
				 enum ADIOS_PREDICATE_MODE op,
				 const char* value); 
					

ADIOS_QUERY* common_query_combine(ADIOS_QUERY* q1, 
				   enum ADIOS_CLAUSE_OP_MODE operator,		    
				   ADIOS_QUERY* q2);

int64_t common_query_estimate(ADIOS_QUERY* q, int timestep);

//void common_query_set_timestep(int timeStep);

int common_query_evaluate(ADIOS_QUERY* q, 
			  ADIOS_SELECTION* outputBoundry,
			  int timestep,
			  uint64_t batchSize, 
			  ADIOS_SELECTION** result);

void common_query_free(ADIOS_QUERY* q);

// called from Read finalize only; 
// this function then calls all query methods' finalize
void common_query_finalize();
  
#ifdef __cplusplus
}
#endif

#endif /* __COMMON_QUERY_H__ */
