/* 
 * ADIOS is freely available under the terms of the BSD license described
 * in the COPYING file in the top level directory of this source distribution.
 *
 * Copyright (c) 2008 - 2009.  UT-BATTELLE, LLC. All rights reserved.
 */

/* ADIOS C Example: 
 *  read all variables and attributes from 
 *    all groups in a BP file
 *
 * This is a sequential program.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "adios_read.h"

const char * value_to_string (enum ADIOS_DATATYPES type, void * data, int idx);

int main (int argc, char ** argv) 
{
    int         i, j, k,l;
    MPI_Comm    comm_dummy = 0;  /* MPI_Comm is defined through adios_read.h */

    if (argc < 2) {
        printf("Usage: %s <BP-file>\n", argv[0]);
        return 1;
    }

    adios_read_init_method (ADIOS_READ_METHOD_BP, comm_dummy, "show_hidden_attrs");
    ADIOS_FILE * f;
    f = adios_read_open_file (argv[1], ADIOS_READ_METHOD_BP, comm_dummy);
    if (f == NULL) {
        printf ("%s\n", adios_errmsg());
        return -1;
    }

    /* For all variables */
    printf("  Variables=%d:\n", f->nvars);
    for (i = 0; i < f->nvars; i++) {
        ADIOS_VARINFO * v = adios_inq_var_byid (f, i);
        adios_inq_var_stat (f, v, 0, 1);
        adios_inq_var_blockinfo (f, v);

        uint64_t total_size = adios_type_size (v->type, v->value);
        for (j = 0; j < v->ndim; j++)
            total_size *= v->dims[j];

        printf("    %-9s  %s", adios_type_to_string(v->type), f->var_namelist[i]);
        if (v->ndim == 0) {
            /* Scalars do not need to be read in, we get it from the metadata
               when using adios_inq_var */
            printf(" = %s\n", value_to_string(v->type, v->value, 0));
        } else {
            /* Arrays, print min/max statistics*/
            printf("[%lld",v->dims[0]);
            for (j = 1; j < v->ndim; j++)
                printf(", %lld",v->dims[j]);
            //printf("] = \n");

            if (v->type == adios_integer)
                printf("] : min=%d  max=%d\n", 
                        (*(int*)v->statistics->min), (*(int*)v->statistics->max));
            else if (v->type == adios_double)
                printf("] : min=%lg  max=%lg\n", 
                        (*(double*)v->statistics->min), (*(double*)v->statistics->max));

            /* Print block info */
            for (l=0; l<v->nsteps; l++) {
                printf("        step %3d: \n", l);
                for (j=0; j<v->nblocks[l]; j++) {
                    printf("          block %3d: [", j);
                    for (k=0; k<v->ndim; k++) {
                        printf("%3lld:%3lld", v->blockinfo[j].start[k],
                                v->blockinfo[j].start[k]+v->blockinfo[j].count[k]-1);
                        if (k<v->ndim-1)
                            printf(", ");
                    }
                    printf("]\n");
                }
            }

        }

        adios_free_varinfo (v);
    } /* variables */

    /* For all attributes */
    printf("  Attributes=%d:\n", f->nattrs);
    for (i = 0; i < f->nattrs; i++) {
        enum ADIOS_DATATYPES atype;
        int  asize;
        void *adata;
        adios_get_attr_byid (f, i, &atype, &asize, &adata);
        printf("    %-9s  %s = %s\n", adios_type_to_string(atype), 
                f->attr_namelist[i], value_to_string(atype, adata, 0));
        free(adata);
    } /* attributes */

    adios_read_close (f);

    return 0;
}


const char * value_to_string (enum ADIOS_DATATYPES type, void * data, int idx)
{
    static char s [100];
    s [0] = 0;


    switch (type)
    {
        case adios_unsigned_byte:
            sprintf (s, "%u", ((uint8_t *) data)[idx]);
            break;

        case adios_byte:
            sprintf (s, "%d", ((int8_t *) data)[idx]);
            break;

        case adios_short:
            sprintf (s, "%hd", ((int16_t *) data)[idx]);
            break;

        case adios_unsigned_short:
            sprintf (s, "%hu", ((uint16_t *) data)[idx]);
            break;

        case adios_integer:
            sprintf (s, "%d", ((int32_t *) data)[idx]);
            break;

        case adios_unsigned_integer:
            sprintf (s, "%u", ((uint32_t *) data)[idx]);
            break;

        case adios_long:
            sprintf (s, "%lld", ((int64_t *) data)[idx]);
            break;

        case adios_unsigned_long:
            sprintf (s, "%llu", ((uint64_t *) data)[idx]);
            break;

        case adios_real:
            sprintf (s, "%g", ((float *) data)[idx]);
            break;

        case adios_double:
            sprintf (s, "%lg", ((double *) data)[idx]);
            break;

        case adios_long_double:
            sprintf (s, "%Lg", ((long double *) data)[idx]);
            break;

        case adios_string:
            return (char*) ((char *)data+idx);
            break;

        case adios_complex:
            sprintf (s, "(%g, %g)", 
                    ((float *) data)[2*idx], ((float *) data)[2*idx+1]);
            break;

        case adios_double_complex:
            sprintf (s, "(%lg, %lg)", 
                    ((double *) data)[2*idx], ((double *) data)[2*idx+1]);
            break;

        default:
            sprintf (s, "unknown");
            break;

    }

    return s;
}
