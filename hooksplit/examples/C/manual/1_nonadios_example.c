/* 
 * ADIOS is freely available under the terms of the BSD license described
 * in the COPYING file in the top level directory of this source distribution.
 *
 * Copyright (c) 2008 - 2009.  UT-BATTELLE, LLC. All rights reserved.
 */

/*
 * Starting example described in the User's manual
 *
 * The initial, non-ADIOS example
 *
 */
#include <stdio.h>
#include "mpi.h"
int main (int argc, char ** argv) 
{
    char     filename [256];
    int      rank;
    int      NX = 10;
    double   t[NX];
    FILE   * fp;
    int      i;

    MPI_Init (&argc, &argv);
    MPI_Comm_rank (MPI_COMM_WORLD, &rank);
    sprintf (filename, "restart_%5.5d.dat", rank);
    
    for (i=0; i<NX; i++)
        t[i] = rank*NX + i;

    fp = fopen (filename, "w");
    fwrite ( &NX, sizeof(int), 1, fp);
    fwrite (t,  sizeof(double), NX, fp);
    fclose (fp);

    MPI_Finalize ();
    return 0;
}

