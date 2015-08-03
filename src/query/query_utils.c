#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <errno.h>
#include <inttypes.h>
#include <assert.h>

#include "core/common_read.h"
#include "core/adios_logger.h"

int query_utils_getGlobalWriteBlockId(int idxRelativeToTimeStep, int timeStep, ADIOS_VARINFO* v) 
{
    int absBlockCounter = idxRelativeToTimeStep;
    int i=0;
    for (i=0; i<v->nsteps; i++)
    {
        int nBlocksAtStep = v->nblocks[i];
        if (i < timeStep) {
            absBlockCounter += nBlocksAtStep;
        }
    }

    return absBlockCounter;
}


int query_utils_file_exists (char * path)
{
    struct stat sb;
    int i = stat ( path, &sb );
    if ( i == 0 )
        /* File found */
        return 1;
    return 0;
}

// records the current time-stamp
double dclock(void) {
	struct timeval tv;
	gettimeofday(&tv, 0);

	return (double) tv.tv_sec + (double) tv.tv_usec * 1e-6;
}
