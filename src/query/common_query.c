#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <inttypes.h>
#include <assert.h>

#include "common_query.h"
#include "adios_query_hooks.h"
#include "public/adios_error.h"
#include "core/common_read.h"
#include "core/adios_logger.h"

static struct adios_query_hooks_struct * query_hooks = 0;


ADIOS_SELECTION* getAdiosDefaultBoundingBox(ADIOS_VARINFO* v) 
{
  if (v->ndim == 0) {
    return NULL;
  }
  /*uint64_t start[v->ndim];
  uint64_t count[v->ndim];

  int i=0;
  for (i=0; i<v->ndim; i++) {
    start[i] = 0;
    count[i] = v->dims[i];
  }
  */
  uint64_t* start = malloc(v->ndim * sizeof(uint64_t));
  uint64_t* count = malloc(v->ndim * sizeof(uint64_t));

  int i=0;
                                                                                                                                                                         
  for (i=0; i<v->ndim; i++) {
    start[i] = 0;
    count[i] = v->dims[i];
  }   

  ADIOS_SELECTION* result =  common_read_selection_boundingbox(v->ndim, start, count);
  return result;
}

void common_query_init()
{
    static int did_init = 0;
    if (!did_init) {
        adios_query_hooks_init(&query_hooks);
        did_init = 1;
    }
}

void common_query_finalize()
{
    static int did_finalize = 0;
    if (did_finalize) {
        enum ADIOS_QUERY_METHOD m;
        for (m=0; m < ADIOS_QUERY_METHOD_COUNT; m++) {
            query_hooks[m].adios_query_finalize_fn();
        }
        free(query_hooks);
        did_finalize = 0;
    }
}


void common_query_set_method (ADIOS_QUERY* q, enum ADIOS_QUERY_METHOD method) 
{
    q->method = method;
}

// Choose a query method which can work on this query
static enum ADIOS_QUERY_METHOD get_method (ADIOS_QUERY* q) 
{
    enum ADIOS_QUERY_METHOD m;
    if (q->method != ADIOS_QUERY_METHOD_UNKNOWN) {
        // Was set by user manually
        return q->method;
    }
    // Look for a method that can evaluate this query
    for (m=0; m < ADIOS_QUERY_METHOD_COUNT; m++) {      
        // without checking whether *evaluate_fn is defined, 
        // it causes crash when idx is not used for fastbit. (i.e. m=0, returns 0, m=1, crashes at "found = nullpoiint(q)"
        if (query_hooks[m].adios_query_can_evaluate_fn == NULL) {
	   continue;
	}
        int found = query_hooks[m].adios_query_can_evaluate_fn(q);
        if (found) {
            return m;
        }
    }
    // return default that always works
    q->method = ADIOS_QUERY_METHOD_FASTBIT;
    return ADIOS_QUERY_METHOD_FASTBIT;
}

void common_query_free(ADIOS_QUERY* q)
{
  if (q->_deleteSelectionWhenFreed) {
    common_read_selection_delete(q->_sel);
  }
  if (q->method < ADIOS_QUERY_METHOD_COUNT) {
      query_hooks[q->method].adios_query_free_fn(q);
  }
}

static int getTotalByteSize (ADIOS_FILE* f, ADIOS_VARINFO* v, ADIOS_SELECTION* sel, 
                             uint64_t* total_byte_size, uint64_t* dataSize)
{
  *total_byte_size = common_read_type_size (v->type, v->value);    
  *dataSize = 1; 

  if (sel == 0) {
    uint64_t s = 0;
    for (s=0; s<v->ndim; s++) {
         *total_byte_size *=v->dims[s];
         *dataSize *= v->dims[s];
         //	log_debug(" dim %" PRIu64 "default count %" PRIu64 "\n", s, v->dims[s]);
    }
    return 0;
  }

  switch (sel->type) {
  case  ADIOS_SELECTION_BOUNDINGBOX:    
    {
      const ADIOS_SELECTION_BOUNDINGBOX_STRUCT *bb = &(sel->u.bb);
      uint64_t* count = bb->count;            
      uint64_t* start = bb->start;            

      int s=0;

      for (s=0; s<v->ndim; s++) {
	   if (start[s]+count[s] > v->dims[s]) {
	     log_error(" Invalid bounding box start %" PRIu64 " + count %" PRIu64 " exceeds dim size: %" PRIu64 "\n", start[s], count[s], v->dims[s]);
	     return -1;
	   }
	   *total_byte_size *=count[s];
	   *dataSize *= count[s];
//	   log_debug(" dim %" PRIu64 "count %" PRIu64 " \n", s, count[s]);
      }
      
//	   log_debug("\tThe data size is = %" PRIu64 " \n", *dataSize);
      break;
    }
  case ADIOS_SELECTION_POINTS:
    {
      const ADIOS_SELECTION_POINTS_STRUCT *pts = &(sel->u.points);
      *total_byte_size *= pts->npoints;
      *dataSize = pts->npoints;
      break;
    }
  case ADIOS_SELECTION_WRITEBLOCK:
    {
      const ADIOS_SELECTION_WRITEBLOCK_STRUCT *wb = &(sel->u.block);

      common_read_inq_var_blockinfo(f, v);
      int i=0;
      int min = v->nblocks[0];
      int absBlockCounter = wb->index;
      for (i=0; i<v->nsteps; i++) 
	{
	  int nBlocksAtStep = v->nblocks[i];	  
	  if (nBlocksAtStep < min) {
	     min = nBlocksAtStep;
	  }
	  log_debug("\t\t   currstep=%d nblocks=%d\n", i, nBlocksAtStep);
	  if (i < gCurrentTimeStep) {
	    absBlockCounter += nBlocksAtStep;
	  }
	}

      if (wb->index > min) {
	  log_error("Error: Unable to handle this block index %d over all the timesteps. Stop.\n", wb->index);
	  return -1;
      }

      int j=0;
      for (j=0; j<v->ndim; j++)
	{
          *total_byte_size *= v->blockinfo[absBlockCounter].count[j];
          *dataSize *= v->blockinfo[absBlockCounter].count[j];
	}

      log_debug("\t\t   block %d, abs id:%d, bytes: %" PRIu64 ", size =  %" PRIu64 " \n", wb->index, absBlockCounter, *total_byte_size, *dataSize);

      break;
    }
  default:
    break;
  }
  return 0;
}

static void initialize(ADIOS_QUERY* result)
{
  result->_onTimeStep = -1; // no data recorded
  result->_maxResultDesired = 0; // init
  result->_lastRead = 0; // init
  result->_hasParent = 0;
  result->_deleteSelectionWhenFreed = 0;
  result->method = ADIOS_QUERY_METHOD_UNKNOWN;
}

ADIOS_QUERY* common_query_create(ADIOS_FILE* f, 
        const char* varName,
        ADIOS_SELECTION* queryBoundary,
        enum ADIOS_PREDICATE_MODE op,
        const char* value)
{
    if (query_hooks == NULL) {
        log_error("ADIOS Query Library Error: Query environment is not initialized.\n");
        exit(EXIT_FAILURE);
    }

    if (queryBoundary != NULL) {
        if ((queryBoundary->type != ADIOS_SELECTION_BOUNDINGBOX)
                && (queryBoundary->type != ADIOS_SELECTION_POINTS)
                && (queryBoundary->type != ADIOS_SELECTION_WRITEBLOCK)) 
        {
            adios_error (err_unsupported_selection, 
                    "Query create: selection type is not supported in queries. "
                    "Choose either boundingbox, points or writeblock selection\n");	       
            return NULL;
        }
    }

    if (value == NULL) {
        adios_error (err_invalid_query_value, "Query create: NULL for value is provided.\n");
        return NULL;
    }
    if (f == NULL) {
        adios_error (err_invalid_file_pointer, "Query create: NULL for input file is provided.\n");
        return NULL;
    }

    // get data from bp file
    ADIOS_VARINFO* v = common_read_inq_var(f, varName);
    if (v == NULL) {
        adios_error (err_invalid_varname, "Query create: Invalid variable '%s':\n%s", 
                     varName, adios_get_last_errmsg());
        return NULL;
    }

    uint64_t total_byte_size, dataSize;

    int defaultBoundaryUsed = 0;
    if (queryBoundary == NULL) {
#ifdef ALACRITY
        queryBoundary = getAdiosDefaultBoundingBox(v);
        defaultBoundaryUsed = 1;
#endif
    }

    if (getTotalByteSize(f, v, queryBoundary, &total_byte_size, &dataSize) < 0) {
        if (defaultBoundaryUsed) {
            common_read_selection_delete(queryBoundary);
        }
        exit(EXIT_FAILURE);
    }   


    //
    // create selection string for fastbit
    //
    ADIOS_QUERY* result = (ADIOS_QUERY*)calloc(1, sizeof(ADIOS_QUERY));
    result->_condition = malloc(strlen(varName)+strlen(value)+ 10); // 10 is enough for op and spaces 
    if (op == ADIOS_LT) {
        sprintf(result->_condition, "(%s < %s)", varName, value);
    } else if (op == ADIOS_LTEQ) {
        sprintf(result->_condition, "(%s <= %s)", varName, value);
    } else if (op == ADIOS_GT) {
        sprintf(result->_condition, "(%s > %s)", varName, value);
    } else if (op == ADIOS_GTEQ) {
        sprintf(result->_condition, "(%s >= %s)", varName, value);
    } else if (op == ADIOS_EQ) {
        sprintf(result->_condition, "(%s = %s)", varName, value);
    } else {
        sprintf(result->_condition, "(%s != %s)", varName, value);
    }

    (result->_condition)[strlen(result->_condition)] = 0;

    initialize(result);

    result->_var = v;
    result->_f = f;

    result->_dataSlice = malloc(total_byte_size);
    result->_rawDataSize = dataSize;

    result->_sel = queryBoundary;
    result->_deleteSelectionWhenFreed = defaultBoundaryUsed;

    result->_op = op;
    result->_value = strdup(value);

    //initialize two pointers
    result->_left = NULL;
    result->_right = NULL;

    return result;
}
					
static int isSelectionCompatible(ADIOS_SELECTION* first, ADIOS_SELECTION* second)			  
{
  if ((first == NULL) || (second == NULL)) {
    return 1;
  }

  switch (first->type) {
  case  ADIOS_SELECTION_BOUNDINGBOX:    
    if (second->type != ADIOS_SELECTION_BOUNDINGBOX) {
        log_error("Error! Not supported: comparing bounding box to another type \n");
	return 0;
    }
    
    return 1;
  case ADIOS_SELECTION_POINTS:
    if (second->type != ADIOS_SELECTION_POINTS) {
        log_error("Error! Not supported: comparing adios points to another type \n");
	return 0;
    }
    const ADIOS_SELECTION_POINTS_STRUCT *pt1 = &(first->u.points);
    const ADIOS_SELECTION_POINTS_STRUCT *pt2 = &(second->u.points);
    
    if (pt1 -> npoints != pt2->npoints) {
      return 0;
    }
    return 1;
  case ADIOS_SELECTION_WRITEBLOCK:
    if (second->type != ADIOS_SELECTION_WRITEBLOCK) {
        log_error("Error! Not supported: comparing adios blocks to another type \n");
	return 0;
    }      
    return 1;
  default:
    return 1;
  }    
}

/* static uint64_t getVariableSize(ADIOS_VARINFO* v) 
{
  uint64_t dataSize = 1;
  int s;
  for (s=0; s<v->ndim; s++) {
    dataSize *= v->dims[s];
  }
  return dataSize;
}*/

//
// return 1 if yes.
//
static int isCompatible(ADIOS_QUERY* q1, ADIOS_QUERY* q2) {
  if (q1->_rawDataSize != q2->_rawDataSize) {
    log_error("Error! Not supported: combining query with different sizes!\n");
    return 0;
  }
  
  if ((q1->_sel != NULL) && (q2->_sel != NULL)) {
    return isSelectionCompatible(q1->_sel, q2->_sel);
  } 

  // all other cases, as long as data sizes match, fastbit can work on it.
  return 1;
}

ADIOS_QUERY* common_query_combine(ADIOS_QUERY* q1, 
				  enum ADIOS_CLAUSE_OP_MODE operator,		    
				  ADIOS_QUERY* q2)
{
    // combine selection sel3 = q1.fastbitSelection & q2.fastbitSelection
    //create a new query (q1.cond :op: q2.cond, sel3);
    //ADIOS_QUERY* result = (ADIOS_QUERY*)malloc(sizeof(ADIOS_QUERY));

    if ((q1 == NULL) || (q2 == NULL)) {
        log_error("Error: detected NULL query when combining.\n");
        adios_error (err_incompatible_queries, "Query combine: NULL passed as query.\n");
        return NULL;
    }

    if (isCompatible(q1, q2) != 1) {
        adios_error (err_incompatible_queries, 
                "Query combine: the two queries' selections are not compatible.\n");
        return NULL;
    }

    ADIOS_QUERY* result = (ADIOS_QUERY*)calloc(1, sizeof(ADIOS_QUERY));
    result->_condition = malloc(strlen(q1->_condition)+strlen(q2->_condition)+10);

    if (operator == ADIOS_QUERY_OP_AND) {
        sprintf(result->_condition, "(%s and %s)", q1->_condition, q2->_condition);
    } else {
        sprintf(result->_condition, "(%s or %s)", q1->_condition, q2->_condition);
    }
    result->_condition[strlen(result->_condition)]=0; 

    q1->_hasParent = 1;
    q2->_hasParent = 1;
    result->_left = q1;
    result->_right = q2;
    result->_leftToRightOp = operator;

    initialize(result);
    return result;
}

int64_t common_query_estimate(ADIOS_QUERY* q)
{
    enum ADIOS_QUERY_METHOD m = get_method (q);
    return query_hooks[m].adios_query_estimate_fn(q);
}

void common_query_set_timestep(int timeStep)
{  
    gCurrentTimeStep = timeStep;
}


static void updateBlockSize(const ADIOS_SELECTION_WRITEBLOCK_STRUCT* wb, ADIOS_QUERY* leaf) 
{
    int i=0;
    int serializedBlockNum = 0;

    ADIOS_VARINFO* v = leaf->_var;
    uint64_t total_byte_size = common_read_type_size (v->type, v->value); ;
    uint64_t dataSize=0;

    for (i=0; i<gCurrentTimeStep; i++) 
    {
        int nBlocksAtStep = v->nblocks[i];  
        serializedBlockNum += nBlocksAtStep;     
    }

    int j=0;
    for (j=0; j<v->ndim; j++)
    {
        total_byte_size *= v->blockinfo[serializedBlockNum].count[j];
        dataSize *= v->blockinfo[serializedBlockNum].count[j];
    }

    log_debug("\t\t   block %d (linedup as %d), bytes: %" PRIu64 ", size = %" PRIu64 " \n", 
            wb->index, serializedBlockNum, total_byte_size, dataSize);

    if (dataSize != leaf->_rawDataSize) {
        log_debug("\t\t reallocate dataSlice due to block size change\n");

        leaf->_rawDataSize = dataSize;
        if (leaf->_dataSlice != NULL) {
            free(leaf->_dataSlice);
        }
        leaf->_dataSlice = malloc(total_byte_size);
    }
}

static int updateBlockSizeIfNeeded(ADIOS_QUERY* q) 
{
    // leaf query
    if ((q->_left == NULL) && (q->_right == NULL)) 
    {
        if (q->_sel == NULL) {
            log_error("No selections detected. \n");
            return -1;
        }

        if (q->_var == NULL) {
            log_error("No variable recorded. \n");
            return -1;
        }

        if (gCurrentTimeStep > q->_var->nsteps) {
            log_error("The given timestep %d exceeds variable (id %d)'s nsteps. \n", gCurrentTimeStep, q->_var->varid);
            return -1;
        }

        if (q->_sel->type != ADIOS_SELECTION_WRITEBLOCK) {
            return 0;
        }
        const ADIOS_SELECTION_WRITEBLOCK_STRUCT *wb = &(q->_sel->u.block);
        updateBlockSize(wb, q);      
        return 1;
    }

    int result = 0;
    if (q->_left != NULL) {
        int leftUpdate = updateBlockSizeIfNeeded(q->_left);
        if (leftUpdate < 0) {
            return -1;
        }
        result += leftUpdate;
    } 
    if (q->_right != NULL) {
        int rightUpdate = updateBlockSizeIfNeeded(q->_right);
        if (rightUpdate < 0) {
            return -1;
        }
        result += rightUpdate;
    }

    return result;
}

static int checkCompatibility(ADIOS_QUERY* q) 
{
    if ((q->_left != NULL) && (q->_right != NULL)) {
        return isCompatible(q->_left, q->_right); 
    }
    return 1; // ok, no need to check  
}

static ADIOS_VARBLOCK * computePGBounds(ADIOS_QUERY *q, int wbindex, int timestep, int *out_ndim) 
{
    if (!q->_left && !q->_right) {
        // In this case, we have reached a leaf query node, so directly
        // retrieve the varblock from the varinfo
        assert(q->_var);

        // Read the blockinfo if not already present
        if (!q->_var->blockinfo) {
            adios_read_set_data_view(q->_f, LOGICAL_DATA_VIEW);
            common_read_inq_var_blockinfo(q->_f, q->_var);
        }

        // Note: adios_get_absolute_writeblock_index ensures that timestep and wbindex
        // are both in bounds, signalling an adios_error if not. However, there will be
        // no variable name cited in the error, so perhaps better error handling would
        // be desirable in the future
        const int abs_wbindex = adios_get_absolute_writeblock_index(q->_var, wbindex, timestep);

        // Finally, return ndim and the varblock
        *out_ndim = q->_var->ndim;
        return &q->_var->blockinfo[abs_wbindex];
    } else if (!q->_left || !q->_right) {
        // In this case, we have only one subtree, so just return the
        // ndim and varblock from that subtree directly, since there's
        // nothing to compare against

        ADIOS_QUERY *present_subtree = q->_left ? (ADIOS_QUERY*)q->_left : (ADIOS_QUERY*)q->_right;
        return computePGBounds(present_subtree, wbindex, timestep, out_ndim);
    } else {
        // In this final case, we have two subtrees, and we must compare
        // the resultant varblock from each one to ensure they are equal
        // before returning

        ADIOS_QUERY *left = (ADIOS_QUERY *)q->_left;
        ADIOS_QUERY *right = (ADIOS_QUERY *)q->_right;

        // Next, retrieve the ndim and varblock for each subtree
        int left_ndim, right_ndim;
        ADIOS_VARBLOCK *left_vb = computePGBounds(left, wbindex, timestep, &left_ndim);
        ADIOS_VARBLOCK *right_vb = computePGBounds(right, wbindex, timestep, &right_ndim);

        // If either subtree returns an invalid (NULL) varblock, fail immediately
        if (!left_vb || !right_vb) {
            return NULL;
        }

        // Check that the ndims are equal, failing if not
        int ndim;
        if (left_ndim != right_ndim) {
            return NULL;
        } else {
            ndim = left_ndim;
        }

        // Check the start/count coordinate in each dimension for equality,
        // failing if any coordinate is not equal between the subtrees
        int i;
        for (i = 0; i < ndim; i++) {
            if (left_vb->start[i] != right_vb->start[i] ||
                    left_vb->count[i] != right_vb->count[i]) {
                return NULL;
            }
        }

        // Finally, we have ensured that both subtrees yield valid and equal
        // varblocks, so return the common ndim and varblock (arbitrarily use
        // left_vb, since right and left equal)
        *out_ndim = ndim;
        return left_vb;
    }
}

static ADIOS_SELECTION * convertWriteblockToBoundingBox(
        ADIOS_QUERY *q, ADIOS_SELECTION_WRITEBLOCK_STRUCT *wb, int timestep) 
{
    assert(!wb->is_absolute_index && !wb->is_sub_pg_selection); // The user should not be using the internal ADIOS writeblock flags

    int pg_ndim;
    ADIOS_VARBLOCK *pg_bounds = computePGBounds(q, wb->index, timestep, &pg_ndim);
    if (!pg_bounds)
    	return NULL;

    ADIOS_SELECTION *bb = common_read_selection_boundingbox(
                        pg_ndim, pg_bounds->start, pg_bounds->count);
    return bb;
}

int common_query_get_selection(ADIOS_QUERY* q, 
			       //const char* varName,
			       //int timeStep, 
			        uint64_t batchSize, // limited by maxResult
				ADIOS_SELECTION* outputBoundary, 
				ADIOS_SELECTION** result)
{
    if ((q->_onTimeStep >= 0) && (q->_onTimeStep != gCurrentTimeStep)) {
        int updateResult = updateBlockSizeIfNeeded(q);
        if (updateResult < 0) {
            log_error("Error with this timestep %d. Can not proceed. \n", gCurrentTimeStep);
            return -1;
        }
        if (updateResult > 0) { // blocks were updated. check compatitibity
            if (checkCompatibility(q) <= 0) {
                return -1;
            }
        }
    }

    int freeOutputBoundary = 0;
    if (outputBoundary->type == ADIOS_SELECTION_WRITEBLOCK) {
        outputBoundary = convertWriteblockToBoundingBox(q, &outputBoundary->u.block, gCurrentTimeStep);
        if (!outputBoundary) {
        	adios_error(err_invalid_argument,
        			    "Attempt to use writeblock output selection on a query where not "
        			    "all variables participating have the same varblock bounding box "
        			    "at that writeblock index (index = %d)\n",
        			    outputBoundary->u.block.index);
        	return -1;
        }
        freeOutputBoundary = 1;
    }

    enum ADIOS_QUERY_METHOD m = get_method (q);

    int retval = query_hooks[m].adios_query_get_selection_fn(
                               q,  batchSize, outputBoundary, result);

    if (freeOutputBoundary) common_read_selection_delete(outputBoundary);
    return retval;
}


