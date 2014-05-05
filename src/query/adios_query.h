#ifndef __ADIOS_QUERY_H__
#define __ADIOS_QUERY_H__

#ifdef __cplusplus
extern "C" {
#endif

#include <adios_read.h>

#define ADIOS_QUERY_TOOL_COUNT  2

enum ADIOS_QUERY_TOOL 
{
        ADIOS_QUERY_TOOL_FASTBIT = 0,
        ADIOS_QUERY_TOOL_OTHER = 1
};
    

enum ADIOS_PREDICATE_MODE 
{
	ADIOS_LT = 0,
	ADIOS_LTEQ = 1,
	ADIOS_GT = 2,
	ADIOS_GTEQ = 3,
	ADIOS_EQ = 4,
	ADIOS_NE = 5
};

enum ADIOS_CLAUSE_OP_MODE 
{
	ADIOS_QUERY_OP_AND = 0,
	ADIOS_QUERY_OP_OR  = 1
};

typedef struct {
    char* condition;
    void* dataSelection;

  // keeping start/count to map 1d results from fastbit to N-d
    uint64_t* _start;
    uint64_t* _count;  

    ADIOS_SELECTION* _sel;
    void* dataSlice;

    ADIOS_VARINFO* var;
    ADIOS_FILE* f;

    enum ADIOS_PREDICATE_MODE op;
    char* _value;
    uint64_t dataSize;

    void* left;
    void* right;
    enum ADIOS_CLAUSE_OP_MODE leftRightOp;
} ADIOS_QUERY;
   


/* functions */
void adios_query_init();

ADIOS_QUERY* adios_query_create(ADIOS_FILE* f, 
				uint64_t *start,
				uint64_t *count,
				const char* varName,
				enum ADIOS_PREDICATE_MODE op,
				const char* value); //this value needs to be &int &double etc, not a string!
					

ADIOS_QUERY* adios_query_combine(ADIOS_QUERY* q1, 
				 enum ADIOS_CLAUSE_OP_MODE operator,		    
				 ADIOS_QUERY* q2);

int64_t adios_query_estimate(ADIOS_QUERY* q, 
			     int timeStep);

int64_t adios_query_evaluate(ADIOS_QUERY* q, 
			     int timeStep,
			     uint64_t maxResult);

void adios_query_get_selection(ADIOS_QUERY* q, 
			       int timeStep, 
			       int batchSize, // limited by maxResult
			       ADIOS_SELECTION* result);

void adios_query_free(ADIOS_QUERY* q);

void adios_query_clean();
  
#ifdef __cplusplus
}
#endif

#endif /* __ADIOS_QUERY_H__ */
