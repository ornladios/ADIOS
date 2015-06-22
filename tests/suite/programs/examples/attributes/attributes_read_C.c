/* 
 * ADIOS is freely available under the terms of the BSD license described
 * in the COPYING file in the top level directory of this source distribution.
 *
 * Copyright (c) 2008 - 2009.  UT-BATTELLE, LLC. All rights reserved.
 */

/* ADIOS C Example: read attributes from a BP file
 *
 * This is possible only with the generic read API.
 * so the GREAD stuff and the xml file is not used here.
 *
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mpi.h"
#include "adios_read.h"
#include "core/adios_logger.h"

int main (int argc, char ** argv) 
{
    int         rank, size, i;
    MPI_Comm    comm = MPI_COMM_WORLD;
    enum ADIOS_DATATYPES attr_type;
    enum ADIOS_READ_METHOD method = ADIOS_READ_METHOD_BP;
    int attr_size;
    void * data = NULL;

    MPI_Init (&argc, &argv);
    MPI_Comm_rank (comm, &rank);
    MPI_Comm_size (comm, &size);

    adios_read_init_method (method, comm, "verbose=3");
    adios_logger_open ("log_read_C", rank);
    ADIOS_FILE * f = adios_read_open ("attributes_C.bp", method, comm, ADIOS_LOCKMODE_NONE, 0.0);
    if (f == NULL)
    {
        log_error ("%s\n", adios_errmsg());
        return -1;
    }

    for (i = 0; i < f->nattrs; i++)
    {

        adios_get_attr (f, f->attr_namelist[i], &attr_type, &attr_size, &data);

        log_test("rank %d: attr: %s %s = ", rank, adios_type_to_string(attr_type), f->attr_namelist[i]);
        int type_size = adios_type_size (attr_type, data);
        int nelems = attr_size / type_size;
        int k;
        char *p = (char*)data;
        for (k=0; k<nelems; k++) 
        {
            if (k>0) log_test(", ");
            switch (attr_type)  
            {
                case adios_integer:
                    log_test ("%d", *(int *)p);
                    break;
                case adios_double:
                    log_test ("%e", *(double *)p);
                    break;
                case adios_string:
                    log_test ("\"%s\"", (char *)p);
                    break;
                case adios_string_array:
                    log_test ("\"%s\"", *(char **)p);
                    break;
                default:
                    log_test ("??????\n");
            }
            p=p+type_size;
        }
        log_test("\n");
        free (data);
        data = 0;
    }

    adios_read_close (f);
    MPI_Barrier (comm);
    adios_read_finalize_method (ADIOS_READ_METHOD_BP);
    adios_logger_close();
    MPI_Finalize ();
    return 0;
}
