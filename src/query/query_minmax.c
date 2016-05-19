/*
 * query_minmax.c
 *
 *  Created on: Jan 2016
 *      Author: Norbert Podhorszki
 */
#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <float.h>
#include "public/adios_error.h"
#include "public/adios_read_ext.h"
#include "public/adios_query.h"
#include "public/adios_selection.h"
#include "core/common_read.h"
#include "core/adios_logger.h"
#include "common_query.h"
#include "query_utils.h"
#include "config.h"  // HAVE_STRTOLD


typedef struct {
    int nblocks;
    char *blocks;  // 0-1 boolean flag for each writeblock, 1=matches query
    int is_outputBoundary_set; // did we set outputBoundary
    ADIOS_SELECTION *outputBoundary; // remember output selection from first eval call (for one step)
    ADIOS_SELECTION *rightmostsel; // rightmost leaf's selection saved at top, from can_evaluate()

    int current_blockid; // end of last evaluation (remember it to be able 
                         // to continue in consecutive evaluate calls)
} MINMAX_INTERNAL;

static void create_internal (ADIOS_QUERY *q)
{
    if (q->queryInternal == NULL) {
        MINMAX_INTERNAL* internal = calloc (1, sizeof(MINMAX_INTERNAL));
        q->queryInternal = (void *)internal;
    }
}

static void free_internal (ADIOS_QUERY *q)
{
    if (q->queryInternal != NULL) {
        MINMAX_INTERNAL* qi = (MINMAX_INTERNAL*) q->queryInternal;
        if (qi->blocks)
            free (qi->blocks);
        free (qi);
        q->queryInternal = NULL;
    }
}

#define INTERNAL(q) ((MINMAX_INTERNAL*) (q->queryInternal))

static void internal_alloc_blocks (ADIOS_QUERY*q, int nblocks)
{
    assert (q->queryInternal); 
    MINMAX_INTERNAL* qi = (MINMAX_INTERNAL*)(q->queryInternal);
    qi->nblocks = nblocks;
    qi->blocks = calloc (nblocks, sizeof(char));
    assert (qi->blocks);
}




#define COMPARE_VALUES(v1,op,v2) {  \
    switch (op) {                   \
        case ADIOS_LT:              \
            return (v1 < v2);       \
            break;                  \
        case ADIOS_LTEQ:            \
            return (v1 <= v2);      \
            break;                  \
        case ADIOS_GT:              \
            return (v1 > v2);       \
            break;                  \
        case ADIOS_GTEQ:            \
            return (v1 >= v2);      \
            break;                  \
        case ADIOS_EQ:              \
            return (v1 == v2);      \
            break;                  \
        case ADIOS_NE:              \
            return (v1 != v2);      \
            break;                  \
    }                               \
}

#if HAVE_STRTOLD
#  define  LONGDOUBLE long double
#  define  STRTOLONGDOUBLE(x,y) strtold(x,y)
#else
#  define  LONGDOUBLE double
#  define  STRTOLONGDOUBLE(x,y) strtod(x,y)
#endif

/* Compare two values with 'op', where one value comes in as string, the other is
   hidden behind a void* pointer and its type depends on the adios 'type'
   returns 1 if the comparison is true, 0 otherwise
 */
static int compare_values (void *v_pred, enum ADIOS_PREDICATE_MODE op, void *v_void, enum ADIOS_DATATYPES vartype)
{
    signed long long v1_int, v2_int;
    unsigned long long v1_uint, v2_uint;
    double v1_real, v2_real;
    LONGDOUBLE v1_ld, v2_ld;

    switch (vartype) 
    {
        case adios_unsigned_byte:
            v2_uint = (unsigned long long) *((unsigned char *) v_void);
            break;
        case adios_byte:
            v2_int = (signed long long) *((signed char *) v_void);
            break;
        case adios_unsigned_short:
            v2_uint = (unsigned long long) *((unsigned short *) v_void);
            break;
        case adios_short:
            v2_int = (signed long long) *((signed short *) v_void);
            break;
        case adios_unsigned_integer:
            v2_uint = (unsigned long long) *((unsigned int *) v_void);
            break;
        case adios_integer:
            v2_int = (signed long long) *((signed int *) v_void);
            break;
        case adios_unsigned_long:
            v2_uint = *((unsigned long long *) v_void);
            break;
        case adios_long:
            v2_int =  *((signed long long *) v_void);
            break;
        case adios_real:
            v2_real = (double) *((float *) v_void);
            break;
        case adios_double:
            v2_real = (double) *((double *) v_void);
            break;
        case adios_long_double:
            v2_ld = *((LONGDOUBLE *) v_void);
            break;

        case adios_complex:
            //fprintf(outf,(f ? format : "(%g,i%g)"), ((float *) data)[2*item], ((float *) data)[2*item+1]);
        case adios_double_complex:
            //fprintf(outf,(f ? format : "(%g,i%g)" ), ((double *) data)[2*item], ((double *) data)[2*item+1]);
        default:
            return 0;
            break;
    } // end switch

    switch (vartype) 
    {
        case adios_unsigned_byte:
        case adios_unsigned_short:
        case adios_unsigned_integer:
        case adios_unsigned_long:
            v1_uint = *(unsigned long long *) v_pred;
            COMPARE_VALUES (v1_uint, op, v2_uint)
            break;

        case adios_byte:
        case adios_short:
        case adios_integer:
        case adios_long:
            v1_int = *(signed long long *) v_pred;
            COMPARE_VALUES (v1_int, op, v2_int)
            break;

        case adios_real:
        case adios_double:
            v1_real = *(double *) v_pred;
            COMPARE_VALUES (v1_real, op, v2_real)
            break;

        case adios_long_double:
            v1_ld = *(LONGDOUBLE *) v_pred;
            COMPARE_VALUES (v1_ld, op, v2_ld)
            break;

        case adios_complex:
        case adios_double_complex:
        default:
            return 0;
            break;
    } // end switch
    return 0;
}

static void * string_to_value (char *v_str, enum ADIOS_DATATYPES vartype)
{
    void * retval;
    static signed long long v_int;
    static unsigned long long v_uint;
    static double v_real;
    LONGDOUBLE v_ld;
    switch (vartype) 
    {
        case adios_unsigned_byte:
        case adios_unsigned_short:
        case adios_unsigned_integer:
        case adios_unsigned_long:
            v_uint = strtoll (v_str, NULL, 10);
            retval = &v_uint;
            break;

        case adios_byte:
        case adios_short:
        case adios_integer:
        case adios_long:
            v_int = strtoll (v_str, NULL, 10);
            retval = &v_int;
            break;

        case adios_real:
        case adios_double:
            v_real = strtod (v_str,NULL);
            retval = &v_real;
            break;

        case adios_long_double:
            v_ld = STRTOLONGDOUBLE(v_str,NULL);
            retval = &v_ld;
            break;

        case adios_complex:
        case adios_double_complex:
        default:
            return 0;
            break;
    } // end switch
    return retval;
}

/*
 * evaluate a single query item (Variable PredicateOP Value) 
 * In: blocks array flag has 1s which writeblocks have to be checked
 *     *sel is the selection used in other part of the query tree 
 * Return the number of matches
 */

static int minmax_evaluate_node (ADIOS_QUERY* q, int timestep, int nblocks, char * blocks, ADIOS_SELECTION **sel, bool estimate) 
{
    // LEAF NODE: evaluate this
    int nmatches = 0;
    int i;
   
    // varinfo contains blocks for many timesteps, we need the index where the current timestep starts
    int block_start_idx = 0; 
    for (i = 0; i < timestep; i++) {
        // FIXME: this may be incorrect for variables that are not written at every timestep into the file
        block_start_idx += q->varinfo->nblocks[i];
    }
    
    // sanity check
    assert (q->varinfo);
    assert (q->varinfo->blockinfo);
    assert (q->varinfo->statistics);
    assert (q->varinfo->statistics->blocks);
    assert (nblocks == q->varinfo->nblocks[timestep]);

    // speed up for special case when selection is a single writeblock
    int loop_start = 0;
    int loop_end = nblocks;
    if (q->sel && q->sel->type == ADIOS_SELECTION_WRITEBLOCK) 
    {
        int index = (q->sel->u.block.is_absolute_index ? 
                     q->sel->u.block.index - block_start_idx :
                     q->sel->u.block.index);
        assert (index > 0);
        assert (index < nblocks);
        memset (blocks, 0, nblocks);
        blocks [ index ] = 1;
        loop_start = index;
        loop_end = index+1;
    }

    void * pred_val = string_to_value (q->predicateValue, q->varinfo->type);

    for (i=loop_start; i < loop_end; i++) 
    {    
        if (blocks[i] && q->sel && *sel != q->sel)
        {
            // if selection is different from previously used selection (in the query tree), then we
            // have to do boundary checks again here
            if (q->sel->type == ADIOS_SELECTION_BOUNDINGBOX) 
            { 
                if (q->varinfo->global) 
                {
                    ADIOS_SELECTION_BOUNDINGBOX_STRUCT bb = q->sel->u.bb;
                    uint64_t *pg_start = q->varinfo->blockinfo[i+block_start_idx].start;
                    uint64_t *pg_count = q->varinfo->blockinfo[i+block_start_idx].count;
                    int k;
                    for (k = 0; k < bb.ndim; k++) {
                        if (bb.start[k]+bb.count[k] <= pg_start[k] ||
                            pg_start[k]+pg_count[k] <= bb.start[k] )
                            blocks[i] = 0;
                    }
                } 
                else
                { // what do we do for non-global arrays?
                }
            }
            /* see speed up above 
            else if (q->sel->type == ADIOS_SELECTION_WRITEBLOCK) 
            {
                int index = (q->sel->u.block.is_absolute_index ? 
                        q->sel->u.block.index - block_start_idx :
                        q->sel->u.block.index);
                if (i != index)
                    blocks[i] = 0;
            }*/
        }

        if (blocks[i])  // block is still in boundary
        {
            // check the formula finally
            switch (q->predicateOp) 
            {
                case ADIOS_LT: 
                    //blocks[i] = (v > q->varinfo->statistics->blocks->mins[i+block_start_idx]);
                    blocks[i] = compare_values (pred_val, ADIOS_GT, q->varinfo->statistics->blocks->mins[i+block_start_idx], q->varinfo->type);
                    break;
                case ADIOS_LTEQ:
                    //blocks[i] = (v >= q->varinfo->statistics->blocks->mins[i+block_start_idx]);
                    blocks[i] = compare_values (pred_val, ADIOS_GTEQ, q->varinfo->statistics->blocks->mins[i+block_start_idx], q->varinfo->type);
                    break;
                case ADIOS_GT:
                    //blocks[i] = (v < q->varinfo->statistics->blocks->maxs[i+block_start_idx]);
                    blocks[i] = compare_values (pred_val, ADIOS_LT, q->varinfo->statistics->blocks->maxs[i+block_start_idx], q->varinfo->type);
                    break;
                case ADIOS_GTEQ:
                    //blocks[i] = (v <= q->varinfo->statistics->blocks->maxs[i+block_start_idx]);
                    blocks[i] = compare_values (pred_val, ADIOS_LTEQ, q->varinfo->statistics->blocks->maxs[i+block_start_idx], q->varinfo->type);
                    break;
                case ADIOS_EQ:
                    // we MAY have a match in block if the predicate value falls inside of the min..max range
                    //blocks[i] = (v >= q->varinfo->statistics->blocks->mins[i+block_start_idx] &&
                    //             v <= q->varinfo->statistics->blocks->maxs[i+block_start_idx]); 
                    blocks[i] = compare_values (pred_val, ADIOS_GTEQ, q->varinfo->statistics->blocks->mins[i+block_start_idx], q->varinfo->type);
                    blocks[i] = compare_values (pred_val, ADIOS_LTEQ, q->varinfo->statistics->blocks->maxs[i+block_start_idx], q->varinfo->type);
                    break;
                case ADIOS_NE:
                    // we only know for sure that the block is not a match if all elements
                    // are the same (min=max) and the predicate value is that same value
                    //blocks[i] = !(v == q->varinfo->statistics->blocks->mins[i+block_start_idx] &&
                    //              v == q->varinfo->statistics->blocks->maxs[i+block_start_idx]); 
                    blocks[i] = !(compare_values (pred_val, ADIOS_EQ, q->varinfo->statistics->blocks->mins[i+block_start_idx], q->varinfo->type) &&
                                  compare_values (pred_val, ADIOS_EQ, q->varinfo->statistics->blocks->maxs[i+block_start_idx], q->varinfo->type));
                    break;
            }
        }

        if (blocks[i])  // block is still matching after evaluation
        {
            nmatches++;
        }
    }


    // update selection going out and up the tree
    *sel = q->sel;
    return nmatches;
}

/*
 * Evaluate the (sub)query on blocks that are still in play (blocks_in)
 * Return new blocks flag array that still match after evaluating the subquery.
 * Also speed up a bit by passing the selection from the previous query node, and if it's
 * the same as for the current query node then boundary check is skipped.
 * At top level, it should be called with a full-1 blocks flag array, as the routine filters out
 * non matches. 
 * At top level, it should be called with a <pointer to a NULL ADIOS_SELECTION*> so that *sel==NULL.
 * At top level return, the number of matches is returned and blocks contain the flags sporadically.
 */
static int minmax_process_rec(ADIOS_QUERY* q, int timestep, int nblocks, char * blocks, ADIOS_SELECTION **sel, bool estimate) 
{
    int nmatches = 0;
    if (q->left == NULL && q->right == NULL) {
        //LEAF NODE: evaluate this
        nmatches = minmax_evaluate_node (q, timestep, nblocks, blocks, sel, estimate);
        return nmatches;
    }

    // combine nodes: evaluate subqueries
    char *rightblocks = blocks;
    int rn, ln;
    if (q->left) {
        ln = minmax_process_rec((ADIOS_QUERY*) q->left, timestep, nblocks, blocks, sel, estimate);
    } else {
        ln = nblocks; // fake value to pass the condition in the right side (AND & 0 skips right side)
    }

    if (q->right) {

        if (q->combineOp == ADIOS_QUERY_OP_OR) {
            // OR operation needs a separate flag array for the right side subquery
            rightblocks = malloc (nblocks * sizeof(char));
            memset (rightblocks, 1, nblocks);
        } else {
            // AND operation simply passes the result flag array from left to right
            rightblocks = blocks;
        }

        if ( q->combineOp != ADIOS_QUERY_OP_AND || ln > 0) {
            rn = minmax_process_rec((ADIOS_QUERY*) q->right, timestep, nblocks, rightblocks, sel, estimate);
        } else {
            rn = 0; // skip evaluating right side since left produced already zero results
        }

        if (q->combineOp == ADIOS_QUERY_OP_OR) {
            // OR operation needs to merge the two arrays with OR 
            int i;
            nmatches = 0;
            for (i = 0; i < nblocks; i++) {
                blocks[i] |= rightblocks[i];
                if (blocks[i])
                    nmatches++;
            }
            free (rightblocks);
        } else {
            nmatches = rn;
        }
    } else {
        nmatches = ln;
    }

    return nmatches;
}


/*
 * Evaluate the expression tree on each writeblock that intersects with the query's selection(s)
 * Uses a flag array to signal which blocks are still in play during evaluation.
 * Each node in the expression tree is visited once and it evaluates the subquery on all blocks still playing.
 * Return the number of matches  (and queryInternal->blocks will have the flags for matching blocks)
 */
static int minmax_process(ADIOS_QUERY* q, int timestep, bool estimate) 
{
    /* At this point, it is ensured that every subquery refers to the same number of writeblocks */
    int nblocks = ((MINMAX_INTERNAL*)(q->queryInternal))->nblocks;
    char *blocks = ((MINMAX_INTERNAL*)(q->queryInternal))->blocks;

    // set every block as match originally
    memset (blocks, 1, nblocks); 

    ADIOS_SELECTION *nullsel = NULL;
    int nmatches = minmax_process_rec(q, timestep, nblocks, blocks, &nullsel, estimate); 
    
    return nmatches;
}


static ADIOS_SELECTION * build_results (ADIOS_QUERY *q, uint64_t retrieval_size, ADIOS_SELECTION* outputBoundry)
{
    int nblocks = INTERNAL(q)->nblocks;
    char *blocks = INTERNAL(q)->blocks;

    /* make the result list of selections from the matching block IDs */
    ADIOS_SELECTION *result = (ADIOS_SELECTION *) calloc (retrieval_size, sizeof(ADIOS_SELECTION));
    ADIOS_SELECTION *r = result;

    int n = retrieval_size;
    int i = INTERNAL(q)->current_blockid;
    assert (i < nblocks);
    for (; i < nblocks; i++)
    {
        if (blocks[i])
        {
            // create a selection in result[j] 
            r->type = ADIOS_SELECTION_WRITEBLOCK;
            r->u.block.index = i;
            r->u.block.is_absolute_index = 1;
            r++;
            n--;
        }
        if (n==0) 
            break;
    }
    assert (i <= nblocks);
    INTERNAL(q)->current_blockid = i;
    return result;    
}


static int selections_are_minmax_compatible (ADIOS_SELECTION *sel1, ADIOS_SELECTION *sel2)
{
    if (sel1==sel2)
        return 1;
    if (!sel1 && sel2)
        return 0;
    if (sel1 && !sel2)
        return 0;
    if (sel1->type == ADIOS_SELECTION_BOUNDINGBOX)
    {
        // the number of dimensions and the sizes should match
        if (sel1->u.bb.ndim != sel2->u.bb.ndim)
            return 0;
        int i=0;
        for (i = 0; i < sel1->u.bb.ndim; i++) {
            if (sel1->u.bb.count[i] != sel2->u.bb.count[i]) {
                return 0;
            }
        }
    }
    else if (sel1->type == ADIOS_SELECTION_WRITEBLOCK) 
    {
        // the writeblock ID must be the same
        if (sel1->u.block.index != sel2->u.block.index)
            return 0;
    }
    else 
    {
        return 0;
    }
    return 1;
}

/* Determine if the query can be evaluated by traversing the tree. 
   As a side effect, all statistics and blockinfo will be retrieved for all variables
*/
static int can_evaluate(ADIOS_QUERY* q, int timestep, ADIOS_SELECTION **sel, int *nblocks)
{
    int supported = 0;
    if (!q->left && !q->right) {
        // If this is a query leaf node, we support MINMAX query iff
        // - the selection is bounding box or writeblock or NULL
        // - the variable has min/max statistics
        if (!q->sel || 
            q->sel->type == ADIOS_SELECTION_BOUNDINGBOX || 
            q->sel->type == ADIOS_SELECTION_WRITEBLOCK
           ) 
        {
            if (!q->varinfo)
                q->varinfo = common_read_inq_var (q->file, q->varName); // get per block statistics
            if (q->varinfo && !q->varinfo->statistics)
                common_read_inq_var_stat (q->file, q->varinfo, 0, 1); // get per block statistics
            if (q->varinfo && !q->varinfo->blockinfo)
                common_read_inq_var_blockinfo (q->file, q->varinfo); // get per block dimensions
            if (q->varinfo  && q->varinfo->statistics  && q->varinfo->statistics->blocks && q->varinfo->blockinfo) {
                supported = 1;
                if (q->sel && q->sel->type == ADIOS_SELECTION_BOUNDINGBOX && q->sel->u.bb.ndim != q->varinfo->ndim) {
                    supported = 0;
                }
                if (q->varinfo->type == adios_complex || q->varinfo->type == adios_double_complex ||
                    q->varinfo->type == adios_string || q->varinfo->type == adios_string_array) {
                    supported = 0;
                }
            }
            *nblocks = q->varinfo->nblocks[timestep];
        } else {
            *nblocks = 0;
        }
        *sel = q->sel;
    } else {
        // Else, this is an internal node, and we support MINMAX query if
        // - ALL descendent nodes support MINMAX
        // - have the same number of writeblocks
        // - the selections are compatible
        int lsupported=1, rsupported=1;
        int l_nblocks, r_nblocks;
        ADIOS_SELECTION *lsel, *rsel;
        if (q->left) {
            lsupported = can_evaluate((ADIOS_QUERY *)q->left, timestep, &lsel, &l_nblocks);
            *nblocks = l_nblocks;
            *sel = lsel;
        }
        if (q->right) {
            rsupported = can_evaluate((ADIOS_QUERY *)q->right, timestep, &rsel, &r_nblocks);
            *nblocks = r_nblocks;
            *sel = rsel;
        }
        supported = lsupported && rsupported;
        if (supported && q->left && q->right) {
            supported = l_nblocks==r_nblocks;
            if (supported)
                supported = selections_are_minmax_compatible (lsel, rsel);
        }
    }
    return supported;
}

// Do the evaluation first time for this timestep
// Return the total number of results available, -1 on error
static int do_evaluate_now (ADIOS_QUERY *q, int timestep)
{
    // run the can_evaluate routine to get the statistics and block info
    ADIOS_SELECTION *qsel;
    int nblocks;
    int supported = can_evaluate (q, timestep, &qsel, &nblocks);
    if (!supported) {
        adios_error (err_incompatible_queries, 
                "%s: the query is not compatible with the minmax query method\n", __func__);
        return -1;
    }

    free_internal (q);
    create_internal (q);
    internal_alloc_blocks (q, nblocks);
    INTERNAL(q)->current_blockid = 0;
    INTERNAL(q)->rightmostsel = qsel;
    q->resultsReadSoFar = 0;
    INTERNAL(q)->is_outputBoundary_set = 0;

    // evaluate query for ALL blocks, fill q->queryInternal->blocks bool array 
    q->maxResultsDesired =  minmax_process(q, timestep, false);
    return q->maxResultsDesired;
} 

/*====================================================================================*/
/*                                  Public functions                                  
*/

int adios_query_minmax_can_evaluate(ADIOS_QUERY* q)
{
    // we can evaluate only iff
    // - every query item is on a compatible selection
    // - that selection is NULL, bounding box or writeblock
    // - the variable in each query item has min/max statistics and
    // - the number of writeblocks of each variable is the same
    int supported = 0;
    ADIOS_SELECTION *sel;
    int nblocks;
    supported = can_evaluate (q, 0, &sel, &nblocks);
    return supported;
}


int64_t adios_query_minmax_estimate(ADIOS_QUERY* q, int timestep) 
{
    const int absoluteTimestep = adios_get_actual_timestep(q, timestep);
    // timestep is always 0 for streaming; the absolute timestep for files
    // absoluteTimestep makes it possible to realize we have a new step
    // in a stream here

    int retval = do_evaluate_now (q, timestep);
    if (retval > -1) {
        // this is treated as the first call to evaluate the query for a new timestep
        // so no need to evaluate again when the evaluate function is called for the same timestep
        q->onTimeStep = absoluteTimestep;
    }
    return retval;
}


void adios_query_minmax_evaluate(ADIOS_QUERY* q,
                   int timestep,
                   uint64_t batchSize, 
                   ADIOS_SELECTION* outputBoundry,
                   ADIOS_QUERY_RESULT * queryResult)
{
#ifdef BREAKDOWN
    double tStart = 0, tEnd = 0;
    tStart = dclock();
#endif
    const int absoluteTimestep = adios_get_actual_timestep(q, timestep);
    // timestep is always 0 for streaming; the absolute timestep for files
    // absoluteTimestep makes it possible to realize we have a new step
    // in a stream here


    if (q->onTimeStep != absoluteTimestep) 
    { 
        // this is the first call to evaluate the query for a new timestep
        /*
        if (qsel != outputBoundry) {
            adios_error (err_incompatible_queries, 
                    "%s: the output boundary selection must be the same as for the query itself "
                    "when using the minmax query method\n", __func__);
            queryResult->status = ADIOS_QUERY_RESULT_ERROR;
            return;
        }
         */

        int nselections = do_evaluate_now (q, timestep);
        if (nselections == -1) {
            queryResult->status = ADIOS_QUERY_RESULT_ERROR;
            return;
        }
        q->onTimeStep = absoluteTimestep;
        INTERNAL(q)->outputBoundary = outputBoundry;
        INTERNAL(q)->is_outputBoundary_set = 1;
    } 
    else 
    { 
        assert (q->queryInternal);
        if (!INTERNAL(q)->is_outputBoundary_set)
        {
            INTERNAL(q)->outputBoundary = outputBoundry;
        }
        else if (((MINMAX_INTERNAL*)(q->queryInternal))->outputBoundary != outputBoundry)
        {
            adios_error (err_incompatible_queries,
                    "%s: follow-up query evaluation calls must use the same outputBoundary selection"
                    "as the first evaluation call\n", __func__);
            queryResult->status = ADIOS_QUERY_RESULT_ERROR;
            return;
        }
    }

    /* FIXME: This check ought to be done before evaluating the query at the first call, but
     * internal->rightmostsel is only set after can_evaluate() was executed in do_evaluate_now().
     */
    if (!selections_are_minmax_compatible (INTERNAL(q)->rightmostsel, INTERNAL(q)->outputBoundary)) {
        adios_error (err_incompatible_queries,
                "%s: the outputBoundary selection is not compatible with the "
                "selections used in the query conditions\n", __func__);
        return;
    }

    // calculate how many results we will return at this time
    uint64_t retrievalSize = q->maxResultsDesired - q->resultsReadSoFar;
    if (retrievalSize <= 0) {
        queryResult->nselections = 0;
        queryResult->selections = NULL;
        queryResult->status = ADIOS_QUERY_NO_MORE_RESULTS;
        return;
    }
    if (retrievalSize > batchSize) {
        retrievalSize = batchSize;
    }

    queryResult->selections = build_results (q, retrievalSize, outputBoundry);
    queryResult->nselections = retrievalSize;
    queryResult->npoints = 0;

    q->resultsReadSoFar += retrievalSize;

#ifdef BREAKDOWN
    tEnd= dclock();
    printf("time [minmax plugin] : %f \n", tEnd - tStart);
#endif

    if (q->resultsReadSoFar < q->maxResultsDesired)
        queryResult->status = ADIOS_QUERY_HAS_MORE_RESULTS;
    else
        queryResult->status = ADIOS_QUERY_NO_MORE_RESULTS;
    return;
}


int adios_query_minmax_free(ADIOS_QUERY* query) {

    if (query == NULL) 
        return 0;
    free_internal (query);
    /*
       Do not free the tree in a bottom-to-up manner
       because every query piece is supposed to be freed
       by the user one by one
    */
    /*
    if (query->left == NULL && query->right == NULL) {
        return 1;
    } else if  (query->right) {
        return adios_query_minmax_free(query->right);
    } else if (query->left) {
        return adios_query_minmax_free(query->left);
    }
    */
    return 1;
}

void adios_query_minmax_finalize() { /* there is nothing to finalize */ }


