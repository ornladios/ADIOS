
#include <adios_read.h>
#include "common_query.h"

void adios_query_init(enum ADIOS_QUERY_TOOL tool) 
{
  common_query_init(tool);
}

ADIOS_QUERY* adios_query_create(ADIOS_FILE* f, 
				uint64_t *start,
				uint64_t *count,
				const char* varName,
				enum ADIOS_PREDICATE_MODE op,
				const char* value)
{
  return common_query_create(f, start, count, varName, op, value);
}
					

ADIOS_QUERY* adios_query_combine(ADIOS_QUERY* q1, 
				 enum ADIOS_CLAUSE_OP_MODE operator,		    
				 ADIOS_QUERY* q2)
{
  return common_query_combine(q1, operator, q2);
}

int64_t adios_query_estimate(ADIOS_QUERY* q, 
			     int timeStep)
{
  return common_query_estimate(q, timeStep);
}

 

int64_t adios_query_evaluate(ADIOS_QUERY* q, 
			     int timeStep,
			     uint64_t maxResult)
{
  return common_query_evaluate(q, timeStep, maxResult);
}

void adios_query_get_selection(ADIOS_QUERY* q, 
			       int timeStep, 
			       int batchSize, // limited by maxResult
			       ADIOS_SELECTION* result)
{
  common_query_get_selection(q, timeStep, batchSize, result);
}

void adios_query_free(ADIOS_QUERY* q)
{
  common_query_free(q);
}

void adios_query_clean()
{
  common_query_clean();
}

