/* 
 * ADIOS is freely available under the terms of the BSD license described
 * in the COPYING file in the top level directory of this source distribution.
 *
 * Copyright (c) 2008 - 2009.  UT-BATTELLE, LLC. All rights reserved.
 */

/* ADIOS C Example: read global arrays from a BP file
 *
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mpi.h"
#include "adios_read.h"

int main (int argc, char ** argv) 
{
    int         rank, size, i, j;
    MPI_Comm    comm = MPI_COMM_WORLD;
    enum ADIOS_READ_METHOD method = ADIOS_READ_METHOD_BP;

    MPI_Init (&argc, &argv);
    MPI_Comm_rank (comm, &rank);
    MPI_Comm_size (comm, &size);

    adios_read_init_method (method, comm, "verbose=3");
    ADIOS_FILE * f = adios_read_open_file("adios_stat.bp", method, comm);
    if (f == NULL)
    {
        fprintf (stderr, "%s\n", adios_errmsg());
        return -1;
    }

    ADIOS_VARINFO * v = adios_inq_var (f, "temperature");
    if (v) {
        /* get statistics, for each individual time-step */ 
        adios_inq_var_stat (f, v, 1, 0);


        printf("Global MIN of temperature: %lf\n", * (double *) v->statistics->min);
        printf("Global MAX of temperature: %lf\n", * (double *) v->statistics->max);
        printf("Global AVG of temperature: %lf\n", * (double *) v->statistics->avg);
        printf("Global STD DEV of temperature: %lf\n", * (double *) v->statistics->std_dev);

        printf("\n");
        printf("---------------------------------------------------------------------------\n");
        if (v->statistics->steps) {
            printf("MIN\t\tMAX\t\tAVG\t\tSTD_DEV\t\tHISTOGRAM\n");
            for(i = 0; i < v->nsteps; i++)
            {
                if (v->statistics->steps->mins[i]) 
                    printf("%lf\t", * (double *) v->statistics->steps->mins[i]);
                else
                    printf("--\t\t");
                if (v->statistics->steps->maxs[i]) 
                    printf("%lf\t", * (double *) v->statistics->steps->maxs[i]);
                else
                    printf("--\t\t");
                if (v->statistics->steps->avgs[i]) 
                    printf("%lf\t", * (double *) v->statistics->steps->avgs[i]);
                else
                    printf("--\t\t");
                if (v->statistics->steps->std_devs[i]) 
                    printf("%lf\t", * (double *) v->statistics->steps->std_devs[i]);
                else
                    printf("--\t\t");
                if (v->statistics->histogram) {
                    for(j = 0; j <= v->statistics->histogram->num_breaks; j++) {
                        printf("%d ", v->statistics->histogram->frequencies[i][j]);
                    }
                    printf("\n");
                } else {
                    printf("--\t\t");
                }
                printf("\n");
            }	
        } else {
            printf ("Per step statistics is missing\n"); 
        }
        printf("---------------------------------------------------------------------------\n");
        printf("\n");

        if (v->statistics->histogram) {
            printf("Break points:\t\t\t");
            for(j = 0; j < v->statistics->histogram->num_breaks; j++)
                printf("%6.2lf\t", v->statistics->histogram->breaks[j]);

            printf("\n");

            printf("Frequencies:\t\t\t");
            for(j = 0; j <= v->statistics->histogram->num_breaks; j++)
                printf("%6d\t", v->statistics->histogram->gfrequencies[j]);
        }

        printf("\n\n");

#if 0 
        printf ("Auto covariance of MIN values of temperature over time: %lf\n", adios_stat_cov (v, v, "min", 0, 12, 0));
        printf ("Auto correlation of MAX values of temperature over time, with lag 2 units: %lf\n", adios_stat_cor (v, v, "max", 0, 8, 2));
        printf("\n\n");
#endif 

        adios_free_varinfo (v);

    } else {
        fprintf (stderr, "ERROR: Cannot inquire statistics of variable 'temperature': %s\n", adios_errmsg());
    }

    v = adios_inq_var (f, "complex");
    if (v) {
        adios_inq_var_stat (f, v, 1, 0);
        double *C = v->statistics->min;
        printf("Global Minimum of variable complex - Magnitude: %lf\n", C[0]);
        printf("Global Minimum of variable complex - Real part: %lf\n", C[1]);
        printf("Global Minimum of variable complex - Imaginary part: %lfi\n", C[2]);

        double ** Cmin;
        Cmin = (double **) v->statistics->steps->mins;

        printf("\nMagnitude\t\tReal\t\t\tImaginary\n");
        for (j = 0; j < v->nsteps; j++)
        {
            printf ("%lf\t\t%lf\t\t%lf\n", Cmin[j][0], Cmin[j][1], Cmin[j][2]);
        }
        printf("\n");
        adios_free_varinfo (v);

    } else {
        fprintf (stderr, "ERROR: Cannot inquire statistics of variable 'complex': %s\n", adios_errmsg());
    }

    adios_read_close (f);
    MPI_Barrier (comm);
    adios_read_finalize_method (method);
    MPI_Finalize ();
    return 0;
}
