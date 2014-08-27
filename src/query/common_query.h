#ifndef __COMMON_QUERY_H__
#define __COMMON_QUERY_H__

#ifdef __cplusplus
extern "C" {
#endif

#include "public/adios_query.h"

#define NO_EVAL_BEFORE -1

void common_query_init(enum ADIOS_QUERY_TOOL tool);

ADIOS_QUERY* common_query_create(ADIOS_FILE* f, 				 
				 const char* varName,
				 ADIOS_SELECTION* queryBoundry,
				 enum ADIOS_PREDICATE_MODE op,
				 const char* value); 
					

ADIOS_QUERY* common_query_combine(ADIOS_QUERY* q1, 
				   enum ADIOS_CLAUSE_OP_MODE operator,		    
				   ADIOS_QUERY* q2);

int64_t common_query_estimate(ADIOS_QUERY* q);

void common_query_set_timestep(int timeStep);

int common_query_get_selection(ADIOS_QUERY* q, 
				uint64_t batchSize, 
				ADIOS_SELECTION* outputBoundry,
				ADIOS_SELECTION** result);

void common_query_free(ADIOS_QUERY* q);
void common_query_clean();
  
#ifdef __cplusplus
}
#endif

#endif /* __ADIOS_QUERY_H__ */
