#include "config.h"
#include <string.h>
#include <unistd.h>
#include <stdint.h>
#include <stdlib.h>

#define __INCLUDED_FROM_FORTRAN_API__
#include "public/adios_read_v2.h"
#include "common_query.h"
#include "core/futils.h"


#ifdef __cplusplus
extern "C"  /* prevent C++ name mangling */
#endif

#ifdef BUILD_WITH_CMAKE
  #include "FC.h"
#endif



int FC_FUNC_(adios_query_is_method_available_f2c,ADIOS_QUERY_IS_METHOD_AVAILABLE_F2C) (int *method)
{
    return common_query_is_method_available((enum ADIOS_QUERY_METHOD)*method);
}

void FC_FUNC_(adios_query_create,ADIOS_QUERY_CREATE) (
        int64_t     * fp, 
        int64_t     * queryBoundary,
        const char  * varName,
        int         * op,
        const char  * value, 
        int64_t     * q,
        int           varName_size,
        int           value_size
        )
{
    char * buf1 = 0;
    char * buf2 = 0;
    ADIOS_QUERY *query;

    buf1 = futils_fstr_to_cstr (varName, varName_size);
    buf2 = futils_fstr_to_cstr (value, value_size);
    if (buf1 != 0 && buf2 != 0) {
        query = common_query_create( (ADIOS_FILE*) *fp, 
                                     (ADIOS_SELECTION*) *queryBoundary, 
                                     buf1, 
                                     (enum ADIOS_PREDICATE_MODE) *op, 
                                     buf2);
        free (buf1);
        free (buf2);
        *q = (int64_t)query;
    } else {
        *q = 0;
    }
}
					

void FC_FUNC_(adios_query_combine,ADIOS_QUERY_COMBINE) (
        int64_t * q1, 
        int     * op, 
        int64_t * q2, 
        int64_t * q)
{
    ADIOS_QUERY *query = common_query_combine( (ADIOS_QUERY*)*q1, 
                                               (enum ADIOS_CLAUSE_OP_MODE)*op, 
                                               (ADIOS_QUERY*)*q2);
    *q = (int64_t)query;
}

void FC_FUNC_(adios_query_set_method,ADIOS_QUERY_SET_METHOD) (int64_t * q, int * method) 
{
    common_query_set_method ((ADIOS_QUERY*)*q, (enum ADIOS_QUERY_METHOD)*method);
}

int64_t FC_FUNC_(adios_query_estimate,ADIOS_QUERY_ESTIMATE) (int64_t * q, int * timestep)
{
    return common_query_estimate((ADIOS_QUERY*)*q, *timestep);
}

 
void FC_FUNC_(adios_query_evaluate,ADIOS_QUERY_EVALUATE) (
        int64_t  * q, 
        int64_t  * sel_outputboundary,
        int      * timestep, 
        uint64_t * batchsize, 
        int64_t  * sel_result,
        int      * err
        )
{
    ADIOS_SELECTION * result;
    *err = common_query_evaluate( (ADIOS_QUERY*)*q, 
                                  (ADIOS_SELECTION*) *sel_outputboundary,
                                  *timestep, 
                                  *batchsize, 
                                  &result);
    if (!*err) {
        *sel_result = (int64_t)result;
    } else {
        *sel_result = 0;
    }
}

void FC_FUNC_(adios_query_free,ADIOS_QUERY_FREE) (int64_t* q)
{
    common_query_free((ADIOS_QUERY*)*q);
}


