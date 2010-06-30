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
    char        filename [256];
    int         rank, size, i, j, k;
    MPI_Comm    comm = MPI_COMM_WORLD;
    void * data = NULL;
    uint64_t start[3], count[3], bytes_read = 0;

    MPI_Init (&argc, &argv);
    MPI_Comm_rank (comm, &rank);
    MPI_Comm_size (comm, &size);

    ADIOS_FILE * f = adios_fopen ("adios_stat.bp", comm);
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

    ADIOS_VARINFO * v = adios_inq_var (g, "temperature");

    printf("Global MIN of temperature: %lf\n", * (double *) v->gmin);
    printf("Global MAX of temperature: %lf\n", * (double *) v->gmax);
    printf("Global AVG of temperature: %lf\n", * (double *) v->gavg);
    printf("Global STD DEV of temperature: %lf\n", * (double *) v->gstd_dev);

	printf("\n");
	printf("---------------------------------------------------------------------------\n");
	printf("MIN\t\tMAX\t\tAVG\t\tSTD_DEV\t\tHISTOGRAM\n");
	for(i = 0; v->ndim >= 0 && (i < v->dims[0]); i ++)
	{
		printf("%lf\t", * (double *) v->mins[i]);
		printf("%lf\t", * (double *) v->maxs[i]);
		printf("%lf\t", * (double *) v->avgs[i]);
		printf("%lf\t", * (double *) v->std_devs[i]);
		for(j = 0; j <= v->hist->num_breaks; j++)
			printf("%d ", v->hist->frequenciess[i][j]);
		printf("\n");
		printf("\n");
	}	
	printf("---------------------------------------------------------------------------\n");
	printf("\n");

	printf("Break points:\t\t\t");
	for(j = 0; j < v->hist->num_breaks; j++)
		printf("%6.2lf\t", v->hist->breaks[j]);

	printf("\n");

	printf("Frequencies:\t\t\t");
	for(j = 0; j <= v->hist->num_breaks; j++)
		printf("%6d\t", v->hist->gfrequencies[j]);

	printf("\n");
	printf ("\nAuto covariance of MIN values of temperature over time: %lf\n", adios_stat_cov (v, v, "min", 0, 12, 0));
	printf ("Auto correlationof MAX values of temperature over time, with lag 2 units: %lf\n", adios_stat_cor (v, v, "max", 0, 8, 2));

	printf("\n\n");
	v = adios_inq_var (g, "complex");
    double *C = v->gmin;
    printf("Global Minimum of variable complex - Magnitude: %lf\n", C[0]);
    printf("Global Minimum of variable complex - Real part: %lf\n", C[1]);
    printf("Global Minimum of variable complex - Imaginary part: %lfi\n", C[2]);

	double ** Cmin;
	Cmin = (double **) v->mins;

	printf("\nMagnitude\t\tReal\t\t\tImaginary\n");
	for (j = 0; v->ndim >= 0 &&  (j < v->dims[0]); j ++)
	{
		printf ("%lf\t\t%lf\t\t%lf\n", Cmin[j][0], Cmin[j][1], Cmin[j][2]);
	}
	printf("\n");

    free (data);
    adios_free_varinfo (v);

    adios_gclose (g);
    adios_fclose (f);

    MPI_Barrier (comm);

    MPI_Finalize ();
    return 0;
}
