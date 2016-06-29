#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <limits.h>
#include <errno.h>
#include <assert.h>
#include <ctype.h>

#include "config.h"
#include "core/util_mpi.h"


static int unique (uint32_t * nids, int size)
{
    int i, j, k;
    uint32_t temp;

    // sort the nids first
    for (i = 1; i < size; i++)
    {
        for (j = 0; j < size - i; j++)
        {
            if (nids[j] > nids[j + 1])
            {
                temp = nids[j];
                nids[j] = nids[j + 1];
                nids[j + 1] = temp;
            }
        }
    }

    // remove duplicates
    i = 0;
    k = 0;
    while (i < size)
    {
        nids[k] = nids[i];

        j = i + 1;
        while (j < size && nids[i] == nids[j])
        {
            j++;
        }
   
        if (j < size)
        {
            k++;
            i = j;            
        }
        else
        {
            break;
        }
    }

    return k + 1;
}

static uint32_t nid_atoi ()
{
    int name_len;
    char * nid_str, * str_buf = malloc (MPI_MAX_PROCESSOR_NAME);
    uint32_t nid;

    MPI_Get_processor_name (str_buf, &name_len);
    nid_str = str_buf;
    while (*nid_str != '\0' && (!isdigit (*nid_str) || *nid_str == '0'))
    {
        nid_str++;
    }

    if (*nid_str == '\0')
    {
        // report an error
    }

    nid = atoi (nid_str);
    free (str_buf);

    return nid;
}

// This helper routine returns a vector of unique NID's.
// It is caller's responsiblity to free nids afterwards.
int get_unique_nids (MPI_Comm comm, uint32_t ** nids)
{
    int size;
    uint32_t my_nid;

    my_nid = nid_atoi ();
    MPI_Comm_size (comm, &size);
    * nids = (uint32_t *) malloc (size * 4);
    assert (* nids);

    MPI_Allgather (&my_nid, 1, MPI_INT,
                   *nids, 1, MPI_INT,
                   comm);
    return unique (*nids, size);
}


