/* ADIOS C Example: read vars from a BP file
 *
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mpi.h"
#include "adios_read.h"

int main (int argc, char ** argv) 
{
    char        filename [256];
    int         rank, size, i, j;
    MPI_Comm    comm = MPI_COMM_WORLD;
    enum ADIOS_DATATYPES attr_type;
    void * data = NULL;
    uint64_t * start = NULL, * count = NULL, bytes_read = 0;

    MPI_Init (&argc, &argv);
    MPI_Comm_rank (comm, &rank);
    MPI_Comm_size (comm, &size);

    ADIOS_FILE * f = adios_fopen ("restart.bp", comm);
    if (f == NULL)
    {
        printf ("%s\n", adios_errmsg());
        return -1;
    }

    ADIOS_GROUP * g = adios_gopen (f, "my_group");
    if (g == NULL)
    {
        printf ("%s\n", adios_errmsg());
        return -1;
    }

    for (i = 0; i < g->vars_count; i++)
    {
        ADIOS_VARINFO * v = adios_inq_var (g, g->var_namelist[i]);

        uint64_t total_size = adios_type_size (v->type, v->value);

        for (j = 0; j < v->ndim; j++)
           total_size *= v->dims[j];

        data = malloc (total_size);
        start = (uint64_t *) malloc (v->ndim * sizeof (uint64_t) + 8);
        count = (uint64_t *) malloc (v->ndim * sizeof (uint64_t) + 8);
        if (data == NULL || start == NULL || count == NULL)
        {
            fprintf (stderr, "malloc failed.\n");
            return -1;
        }

        for (j = 0; j < v->ndim; j++)
        {
            start[j] = 0;
            count[j] = v->dims[j];   
        }
       
        bytes_read = adios_read_var (g, g->var_namelist[i], start, count, data);

        if (v->ndim == 0)
        {
            if (v->type == adios_integer)
                printf ("%s: %d\n", g->var_namelist[i], * (int *)data);
            else if (v->type == adios_double)
                printf ("%s: %e\n", g->var_namelist[i], * (double *)data);
            else if (v->type == adios_string)
            {
                printf ("%s: %s\n", g->var_namelist[i], data);
                printf ("total size = %llu\n", total_size);
            }

        }
        else if (v->ndim == 1)
        {
        }

        free (count);
        free (start);
        free (data);
    }

    adios_gclose (g);
    adios_fclose (f);

    MPI_Barrier (comm);

    MPI_Finalize ();
    return 0;
}
