#include "public/adios_read.h"
#include "common_query.h"

int adios_query_is_method_available(enum ADIOS_QUERY_METHOD method) {
	return common_query_is_method_available(method);
}

ADIOS_QUERY* adios_query_create(ADIOS_FILE* f, 
				const char* varName,
				ADIOS_SELECTION* queryBoundary,
				enum ADIOS_PREDICATE_MODE op,
				const char* value)
{
  return common_query_create(f, varName, queryBoundary, op, value);
}
					

ADIOS_QUERY* adios_query_combine(ADIOS_QUERY* q1, 
				 enum ADIOS_CLAUSE_OP_MODE operator,		    
				 ADIOS_QUERY* q2)
{
  return common_query_combine(q1, operator, q2);
}

void adios_query_set_method (ADIOS_QUERY* q, enum ADIOS_QUERY_METHOD method) 
{
     common_query_set_method (q, method);
}

int64_t adios_query_estimate(ADIOS_QUERY* q, int timestep)
{
  return common_query_estimate(q, timestep);
}

 
/* //obsolete
void adios_query_set_timestep(int timeStep)
{
  return common_query_set_timestep(timeStep);
}
*/
int  adios_query_evaluate(ADIOS_QUERY* q, 
			  int timeStep, 
			  uint64_t batchSize, // limited by maxResult
			  ADIOS_SELECTION* outputBoundary,
			  ADIOS_SELECTION** queryResult)
{
  return common_query_evaluate(q, timeStep,  batchSize, outputBoundary, queryResult);
}

void adios_query_free(ADIOS_QUERY* q)
{
  common_query_free(q);
}


