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
    int         i, j, k, l, t;
    MPI_Comm    comm_dummy = 0;  /* MPI_Comm is defined through adios_read.h */
    void      * data = NULL;
    uint64_t    start[] = {0,0,0,0,0,0,0,0,0,0};
    uint64_t    count[10];
    ADIOS_SELECTION *sel;

    if (argc < 2) {
        printf("Usage: %s <BP-file>\n", argv[0]);
        return 1;
    }

    ADIOS_FILE * f;
    //int step;
    //for (step=0; step < 2; step++) {
        f = adios_read_open_file (argv[1], ADIOS_READ_METHOD_BP, comm_dummy);
        if (f == NULL) {
            printf ("%s\n", adios_errmsg());
            return -1;
        }

        /* For all variables */
        printf("  Variables=%d:\n", f->nvars);
        for (i = 0; i < f->nvars; i++) {
            ADIOS_VARINFO * v = adios_inq_var_byid (f, i);
            adios_inq_var_stat (f, v, 0, 0);

            uint64_t total_size = adios_type_size (v->type, v->value);
            for (j = 0; j < v->ndim; j++)
                total_size *= v->dims[j];

            printf("    %-9s  %s", adios_type_to_string(v->type), f->var_namelist[i]);
            if (v->ndim == 0) {
                /* Scalars do not need to be read in, we get it from the metadata
                   when using adios_inq_var */
                printf(" = %s\n", value_to_string(v->type, v->value, 0));
            } else {
                /* Arrays have to be read in from the file */
                if (v->nsteps > 1) {
                    printf(" %d*",v->nsteps);
                }
                printf("[%" PRIu64,v->dims[0]);
                for (j = 1; j < v->ndim; j++)
                    printf(", %" PRIu64,v->dims[j]);
                //printf("] = \n");
                
                if (v->type == adios_integer)
                    printf("] = min=%d  max=%d\n", (*(int*)v->statistics->min), (*(int*)v->statistics->max));
                else if (v->type == adios_double)
                    printf("] = min=%lg  max=%lg\n", (*(double*)v->statistics->min), (*(double*)v->statistics->max));
                
                if (total_size > 1024*1024*1024) {
                    printf("        // too big, do not read in\n");
                } else {
                    data = malloc (total_size);
                    if (data == NULL) {
                        fprintf (stderr, "malloc failed.\n");
                        return -1;
                    }

                    for (j = 0; j < v->ndim; j++) 
                        count[j] = v->dims[j];   

                    for (t=0; t<v->nsteps; t++) {
                        sel = adios_selection_boundingbox (v->ndim, start, count);
                        adios_schedule_read_byid (f, sel, i, t, 1, data);
                        adios_perform_reads (f, 1);

                        printf("      Step %d:\n", t);
                        if (adios_errno) {
                            printf ("%s\n", adios_errmsg());
                        } else if (total_size > 1024*1024) {
                            printf ("Too big to print\n");
                        } else if (v->ndim == 1) {
                            printf ("        [");
                            for (j = 0; j < v->dims[0]; j++) 
                                printf("%s ", value_to_string(v->type, data, j));
                            printf ("]\n");
                        } else if (v->ndim == 2) {
                            for (j = 0; j < v->dims[0]; j++) {
                                printf ("        row %d: [", j);
                                for (k = 0; k < v->dims[1]; k++) 
                                    printf("%s ", value_to_string(v->type, data, j*v->dims[1] + k));
                                printf ("]\n");
                            }
                        } else if (v->ndim == 3) {
                            for (j = 0; j < v->dims[0]; j++) {
                                printf ("      block %d: \n", j);
                                for (k = 0; k < v->dims[1]; k++) {
                                    printf ("        row %d: [", k);
                                    for (l = 0; l < v->dims[2]; l++) {
                                        // NCSU ALACRITY-ADIOS - Fixed bug, k*v->dims[1] changed to  k*v->dims[2]
                                        printf("%s ", value_to_string(v->type, data, j*v->dims[1]*v->dims[2] + k*v->dims[2] + l));
                                    }
                                    printf ("]\n");
                                }
                                printf ("\n");
                            }
                        } else {
                            printf ("    cannot print arrays with >3 dimensions\n");
                        }
                    }
                    free (data);
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
            int type_size = adios_type_size (atype, adata);
            int nelems = asize / type_size;
            printf("    %-9s  %s = ", adios_type_to_string(atype), f->attr_namelist[i]);
            char *p = (char*)adata;
            if (nelems>1) printf("{");
            for (j=0; j<nelems; j++) {
                if (j>0) printf(", ");
                printf ("%s", value_to_string(atype, p, 0));
                p += type_size;
            }
            if (nelems>1) printf("}");
            printf("\n");
            free(adata);
        } /* attributes */

        adios_read_close (f);

    //} /* loop 'step' */
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
            sprintf (s, "%" PRId64, ((int64_t *) data)[idx]);
            break;

        case adios_unsigned_long:
            sprintf (s, "%" PRIu64, ((uint64_t *) data)[idx]);
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

        case adios_string_array:
            return (char*) *((char **)data+idx);
            break;

        case adios_complex:
            sprintf (s, "(%g, %g)", 
                    ((float *) data)[2*idx], ((float *) data)[2*idx+1]);
            break;

        case adios_double_complex:
            sprintf (s, "(%lg, %lg)", 
                    ((double *) data)[2*idx], ((double *) data)[2*idx+1]);
            break;
        
        case adios_unknown:
            break;
    }

    return s;
}
