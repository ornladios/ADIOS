#ifndef __QUERY_UTILS_H__
#define __QUERY_UTILS_H__

#ifdef __cplusplus
extern "C" {
#endif

#include "public/adios_query.h"

/* helper functions for all query methods */


/* Return global writeblock ID (among all timesteps) of 
   a writeblock relative to a given timestep */
int query_utils_getGlobalWriteBlockId(int idxRelativeToTimeStep, int timeStep, ADIOS_VARINFO* v);


/* Return 1 if path exists on the file system, 0 otherwise */
int query_utils_file_exists (char * path);

//current timestamp
double dclock(void);

#ifdef __cplusplus
}
#endif

#endif /* __QUERY_UTILS_H__ */
