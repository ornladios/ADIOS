#include "public/adios_read.h"
#include "common_query.h"

void adios_query_init(enum ADIOS_QUERY_TOOL tool) 
{
  common_query_init(tool);
}

ADIOS_QUERY* adios_query_create(ADIOS_FILE* f, 
				const char* varName,
				ADIOS_SELECTION* queryBoundry,
				enum ADIOS_PREDICATE_MODE op,
				const char* value)
{
  return common_query_create(f, varName, queryBoundry, op, value);
}
					

ADIOS_QUERY* adios_query_combine(ADIOS_QUERY* q1, 
				 enum ADIOS_CLAUSE_OP_MODE operator,		    
				 ADIOS_QUERY* q2)
{
  return common_query_combine(q1, operator, q2);
}

int64_t adios_query_estimate(ADIOS_QUERY* q)
{
  return common_query_estimate(q);
}

 
void adios_query_set_timestep(int timeStep)
{
  return common_query_set_timestep(timeStep);
}

int  adios_query_get_selection(ADIOS_QUERY* q, 
			       //const char* varName,
			       //int timeStep, 
			       uint64_t batchSize, // limited by maxResult
			       ADIOS_SELECTION* outputBoundry,
			       ADIOS_SELECTION** queryResult)
{
  common_query_get_selection(q,  batchSize, outputBoundry, queryResult);
}

void adios_query_free(ADIOS_QUERY* q)
{
  common_query_free(q);
}

void adios_query_clean()
{
  common_query_clean();
}

