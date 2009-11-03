/* ADIOS C Example: read attributes from a BP file
 *
*/
#include <stdio.h>
#include <string.h>
#include "mpi.h"
#include "adios_read.h"

int main (int argc, char ** argv) 
{
    char        filename [256];
    int         rank, size, i;
    int         NX = 10, level = 1; 
    double      t[NX], mean = 0;
    MPI_Comm    comm = MPI_COMM_WORLD;
    enum ADIOS_DATATYPES attr_type;
    int attr_size;
    void * data = NULL;

    MPI_Init (&argc, &argv);
    MPI_Comm_rank (comm, &rank);
    MPI_Comm_size (comm, &size);

    ADIOS_FILE * f = adios_fopen ("restart.bp", comm);
    if (f == NULL)
    {
        printf ("%s\n", adios_errmsg());
        return -1;
    }

    ADIOS_GROUP * g = adios_gopen (f, "temperature");
    if (g == NULL)
    {
        printf ("%s\n", adios_errmsg());
        return -1;
    }

    for (i = 0; i < g->attrs_count; i++)
    {

        adios_get_attr (g, g->attr_namelist[i], &attr_type, &attr_size, &data);

        if (attr_type == adios_integer)
        {
            printf ("Attr: %s = %d\n", g->attr_namelist[i], * (int *)data);
        }
        else if (attr_type == adios_double)
        {
            printf ("Attr: %s = %e\n", g->attr_namelist[i], * (double *)data);
        }
        else if (attr_type == adios_string)
        {
            printf ("Attr: %s = %s\n", g->attr_namelist[i], (char *)data);
        }
        free (data);
        data = 0;
    }

    adios_gclose (g);
    adios_fclose (f);

    MPI_Barrier (comm);

    MPI_Finalize ();
    return 0;
}
