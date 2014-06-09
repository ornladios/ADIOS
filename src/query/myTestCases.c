#include "adios_query.h"

void multiBoundBox(ADIOS_FILE* f)
{
  printf("\n=============== testing multiple bound box  ===========\n");
  uint64_t start1[] = {0, 0, 0};
  uint64_t count1[] = {256, 1,32};

  ADIOS_SELECTION* box1 = adios_selection_boundingbox(3, start1, count1);

  uint64_t start2[] = {0, 1, 0};
  uint64_t count2[] = {256, 1,32};
  ADIOS_SELECTION* box2 = adios_selection_boundingbox(3, start2, count2);

  const char* varName1 = "/Timestep_0/cells/X";
  enum ADIOS_PREDICATE_MODE op1 = ADIOS_GT;
  const char* value1 = "0.96874";

  const char* varName2 = "/Timestep_0/cells/X";
  enum ADIOS_PREDICATE_MODE op2 = ADIOS_LT;
  const char* value2 = "0.96876";

  ADIOS_QUERY* q1 = adios_query_create(f, varName1, box1, op1, value1);
  ADIOS_QUERY* q2 = adios_query_create(f, varName2, box2, op2, value2);

  // if box has a different shape, e.g. different count values, then combine() returns error
  ADIOS_QUERY* q = adios_query_combine(q1, ADIOS_QUERY_OP_AND, q2);

  if (q != NULL) {
    int timestep = 0;
    uint64_t max = 10000;
    
    //int64_t hitSize = adios_query_evaluate(q, timestep, max);
    int64_t batchSize = 50;
    
    // box3 is the same shape as other boxes
    // if it is in different shape from box1/box2, then error
    while (1) {
      ADIOS_SELECTION* currBatch = NULL;
      //ADIOS_SELECTION* box3 = adios_selection_boundingbox(3, start3, count3);
      int hasMore =  adios_query_get_selection(q, batchSize, box1, &currBatch);
      adios_selection_delete(currBatch);
      
      if (hasMore == 0) {
	break;
      }
    }
    
    fastbit_selection_free(q->_queryInternal);
    adios_query_free(q);
  }

  adios_query_free(q1);
  adios_query_free(q2);

  adios_selection_delete(box1);
  adios_selection_delete(box2);
}

void defaultBoundBox(ADIOS_FILE* f) 
{
  printf("\n=============== testing default bound box (no box specified) for all ===========\n");
  ADIOS_SELECTION* box = 0;

  const char* varName1 = "/Timestep_0/cells/X";
  enum ADIOS_PREDICATE_MODE op1 = ADIOS_LT;
  const char* value1 = "0.96874";

  const char* varName2 = "/Timestep_0/cells/Y";
  enum ADIOS_PREDICATE_MODE op2 = ADIOS_GT;
  const char* value2 = "0.96874";

  ADIOS_QUERY* q1 = adios_query_create(f, varName1, box, op1, value2);
  ADIOS_QUERY* q2 = adios_query_create(f, varName2, box, op2, value2);

  ADIOS_QUERY* q = adios_query_combine(q1, ADIOS_QUERY_OP_AND, q2);

  if (q != NULL) {
    int timestep = 0;
    uint64_t max = 10000;
    //int64_t hitSize = adios_query_evaluate(q, timestep, max);
    
    int64_t batchSize = 50;
    while (1) {
      ADIOS_SELECTION* currBatch = NULL;
      int hasMore =  adios_query_get_selection(q, batchSize, box, &currBatch);
      adios_selection_delete(currBatch);
      
      if (hasMore == 0) {
	break;
      }
    }

    fastbit_selection_free(q->_queryInternal);
    adios_query_free(q);
  }

  adios_query_free(q1);
  adios_query_free(q2);

  adios_selection_delete(box);

}

void onePointList(ADIOS_FILE* f)
{
 printf("\n=============== testing a point array for all ===========\n");

  int numPoints = 6;

  //uint64_t points[] = {200,0,31, 199,0,21,  201,0,31, 198,0,21, 197,0,21, 196,0,21, 202,0,31,  203,0,31,  204,0,31,  205,0,31,  206,0,31, 207,0,31, 195,0,21, 208,0,31,209,0,31};
  uint64_t points[] = {1,31,10,    1,30,10,   1,31,11,     3,31,10,  1,31,28, 200,0,31};  
  ADIOS_SELECTION* box = adios_selection_points(3, numPoints, points);

  const char* varName1 = "/Timestep_0/cells/X";
  enum ADIOS_PREDICATE_MODE op1 = ADIOS_LT;
  const char* value1 = "0.96874";

  const char* varName2 = "/Timestep_0/cells/Y";
  enum ADIOS_PREDICATE_MODE op2 = ADIOS_GT;
  const char* value2 = "0.96874";

  ADIOS_QUERY* q1 = adios_query_create(f, varName1, box, op1, value2);
  ADIOS_QUERY* q2 = adios_query_create(f, varName2, box, op2, value2);

  ADIOS_QUERY* q = adios_query_combine(q1, ADIOS_QUERY_OP_AND, q2);

  if (q!= NULL) {
    int timestep = 0;
    uint64_t max = 10000;
    //int64_t hitSize = adios_query_evaluate(q, timestep, max);
    
    int64_t batchSize = 50;
    while (1) {
      ADIOS_SELECTION* currBatch = NULL;
      int hasMore =  adios_query_get_selection(q, batchSize, box, &currBatch);
      adios_selection_delete(currBatch);
      
      if (hasMore == 0) {
	break;
      }
    }
    
    fastbit_selection_free(q->_queryInternal);
    adios_query_free(q);
  }

  adios_query_free(q1);
  adios_query_free(q2);

  adios_selection_delete(box);
}

void oneBoundBoxForAllVar(ADIOS_FILE* f)
{
  printf("\n=============== testing one bound box for all ===========\n");
  uint64_t start[] = {1, 10, 10};
  uint64_t count[] = {1, 22, 22};

  ADIOS_SELECTION* box = adios_selection_boundingbox(3, start, count);

  const char* varName1 = "/Timestep_0/cells/X";
  enum ADIOS_PREDICATE_MODE op1 = ADIOS_LT;
  const char* value1 = "0.96874";

  const char* varName2 = "/Timestep_0/cells/Y";
  enum ADIOS_PREDICATE_MODE op2 = ADIOS_GT;
  const char* value2 = "0.96874";

  ADIOS_QUERY* q1 = adios_query_create(f, varName1, box, op1, value2);
  ADIOS_QUERY* q2 = adios_query_create(f, varName2, box, op2, value2);

  ADIOS_QUERY* q = adios_query_combine(q1, ADIOS_QUERY_OP_AND, q2);

  if (q!= NULL) {
    int timestep = 0;
    uint64_t max = 10000;
    //int64_t hitSize = adios_query_evaluate(q, timestep, max);
    
    int64_t batchSize = 50;
    
    int i = 0;
    printf("times steps for variable is: %d \n",q1->_var->nsteps);
    for (i=0; i<q1->_var->nsteps; i++) {
      adios_query_set_timestep(i);
      
      while (1) {
	ADIOS_SELECTION* currBatch = NULL;
	int hasMore =  adios_query_get_selection(q, batchSize, box, &currBatch);
	adios_selection_delete(currBatch);
	
	if (hasMore == 0) {
	  break;
	}
      }
      
    }

    fastbit_selection_free(q->_queryInternal);
    adios_query_free(q);
  }
  adios_query_free(q1);
  adios_query_free(q2);

  adios_selection_delete(box);
}


int main (int argc, char ** argv) 
{
    if (argc < 2) {
        usage(argv[0]);
	return 1;
    }

    adios_query_init(ADIOS_QUERY_TOOL_FASTBIT);

    ADIOS_FILE * f;
    MPI_Comm    comm_dummy = 0;  /* MPI_Comm is defined through adios_read.h */

    f = adios_read_open_file (argv[1], ADIOS_READ_METHOD_BP, comm_dummy);
    if (f == NULL) {
        printf ("::%s\n", adios_errmsg());
	return -1;
    }
    
    defaultBoundBox(f);
    multiBoundBox(f);
    
    oneBoundBoxForAllVar(f); 
    onePointList(f);

    adios_query_clean();
    adios_read_close(f);
    return 1;
    
}
