/* 
 * ADIOS is freely available under the terms of the BSD license described
 * in the COPYING file in the top level directory of this source distribution.
 *
 * Copyright (c) 2008 - 2009.  UT-BATTELLE, LLC. All rights reserved.
 */

/* ADIOS C Example: write variables along with an unstructured mesh. 
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <fcntl.h>
#include <errno.h>
#include "mpi.h"
#include "public/adios.h"
#include "public/adios_types.h"

//will work with 12 cores, which are arranged by npx=4, npy=3 (4x3)
char   npx_str[256];       // # of procs in x dim (string value)
char   npy_str[256];       // # of procs in y dim (string value)
int    npx;                // # of procs in x direction
int    npy;                // # of procs in y direction
int    nproc;               // # of total procs

void printUsage(char *prgname)
{
    printf("Usage: mpirun -np <N> %s <nx> <ny>\n"
           "    <nx> <ny>  2D decomposition values in each dimension of an 2D array\n"
           "         The product of these number must be equal the number of processes\n"
           "         e.g. for N=12 you may use  4 3\n"
        ,prgname);
}


int processArgs(int argc, char ** argv)
{
    if (argc < 3) {
        printUsage (argv[0]);
        return 1;
    }

    strncpy(npx_str, argv[1], sizeof(npx_str));
    strncpy(npy_str, argv[2], sizeof(npy_str));

    npx = atoi(npx_str);
    npy = atoi(npy_str);

    if (npx*npy != nproc) {
        printf ("ERROR: Product of decomposition numbers in X and Y dimension %d != number of processes %d\n", npx*npy, nproc);
        printUsage(argv[0]);
        return 1;
    }

    return 0;
}


int main (int argc, char ** argv) 
{
    MPI_Comm    comm = MPI_COMM_WORLD;
    int         npoints, num_cells;
    int         rank;
    int         ndx, ndy;             // size of array per processor
    double      *N;                   // node centered variable
    double      *C;                   // cell centered variable
    double      *points;              //X,Y coordinate
    int         *cells;

    // Offsets and sizes
    int         offs_x, offs_y;       //offset in x and y direction
    int         nx_local, ny_local;   //local address
    int         nx_global, ny_global; //global address
    int         posx, posy;           // position index in the array
    int         i,j;

    int64_t     m_adios_group;
    uint64_t    adios_groupsize, adios_totalsize;
    int64_t     adios_handle;

    MPI_Init (&argc, &argv);
    MPI_Comm_rank (comm, &rank);
    MPI_Comm_size (comm, &nproc);

    if (processArgs(argc, argv)) {
        return 1;
    }

    //will work with each core writing ndx = 65, ndy = 129, (65*4,129*3) global
    ndx = 4;
    ndy = 3;

    npoints = ndx * ndy * npx * npy;
    num_cells = (ndx * npx - 1) * (ndy * npy -1) * 2;

    //2D array with block,block decomposition
    posx = rank%npx;           // 1st dim
    posy = rank/npx;           // 2nd dim
    offs_x = posx * ndx;
    offs_y = posy * ndy;
    nx_local = ndx;
    ny_local = ndy;
    nx_global = npx * ndx;
    ny_global = npy * ndy;

    // local mesh + data
    N = malloc (ndx * ndy * sizeof(double));
    points = malloc (2 * ndx * ndy * sizeof(double));

    if( posx == npx-1 && posy < npy-1 )            //only need to extend in Y direction
    {
        C = malloc (sizeof(double) * (ndx-1)*ndy*2);
        cells = malloc(sizeof(int) * ((ndx-1)*ndy*2*3));
    }
    else if( posy == npy-1 && posx < npx-1 )      ////only need to extend in X direction
    {
        C = malloc (sizeof(double) * ndx*(ndy-1)*2);
        cells = malloc(sizeof(int) * ndx*(ndy-1)*2*3);
//        printf("rank %d cells size is %d\n", rank, ndx*(ndy-1)*2*3);
    }
    else if( posx == npx-1 && posy == npy-1 )    //do not need to extend in any direction
    {
        C = malloc (sizeof(double) * (ndx-1)*(ndy-1)*2);
        cells = malloc(sizeof(int) * (ndx-1)*(ndy-1)*2*3);
//        printf("rank %d cells size is %d\n", rank, ndx*(ndy-1)*2*3);
    }
    else               // if( posx < npx-1 && posy < npy-1 )   //need to extend in both X and Y direction
    {
        C = malloc (sizeof(double) * ndx*ndy*2);
        cells = malloc(sizeof(int) * ndx*ndy*2*3);
    }

    // generate local data
    int lp = ndx * ndy; // number of points in this local array
    int op = ndx * ndy * rank; // offset in global points array 
    for( i = 0; i < ndx; i++ )
        for( j = 0; j < ndy; j++)
        {
            points[(i*ndy + j)*2] = offs_x + posy*ndx + i*ndx/ndx + (double)ndx*j/ndy;
            points[(i*ndy + j)*2+1] = offs_y + ndy*j/ndy;
        }

    for( i = 0; i < lp; i++ )
        N[i] = 1.0*rank;

    // cells
    int lc;
    int oc;

    if( posx == npx-1 && posy < npy-1 )
    {
        lc = (ndx-1)*ndy*2;
        oc = posx*ndx*ndy*2 + posy*(ndx*npx-1)*ndy*2;
        for (i = 0; i < ndx-1; i++)
        {
            for (j = 0; j < ndy; j++)
            {
                int p = i*ndy+j;
                if( i<ndx-1 && j<ndy-1 )
                {
                    cells[6*p+0] = op+p; cells[6*p+1] = op+p+ndy; cells[6*p+2] = op+p+ndy+1;
                    cells[6*p+3] = op+p; cells[6*p+4] =op+p+ndy+1; cells[6*p+5] = op+p+1;
                }
                else   //extend in Y direction only
                {
                    cells[6*p+0] = op+p; cells[6*p+1] = op+p+ndy; cells[6*p+2] = op+nx_global*ndy+(i+1)*ndy;
                    cells[6*p+3] = op+p; cells[6*p+4] = op+nx_global*ndy+(i+1)*ndy; cells[6*p+5] = op+nx_global*ndy+i*ndy;
                }
            }
        }
    }
    else if( posy == npy-1 && posx < npx-1 )
    {
        lc = ndx*(ndy-1)*2;
        oc = posy*(ndx*npx-1)*ndy*2 + posx*ndx*(ndy-1)*2;
        for (i = 0; i < ndx; i++)
            for (j = 0; j < ndy-1; j++)
            {
                int p = i*(ndy-1)+j;
                int p1 = i*ndy+j;
                if( i<ndx-1 && j<ndy-1 )
                {
                    cells[6*p+0] = op+p1; cells[6*p+1] = op+p1+ndy; cells[6*p+2] = op+p1+ndy+1;
                    cells[6*p+3] = op+p1; cells[6*p+4] =op+p1+ndy+1; cells[6*p+5] = op+p1+1;
                }
                else   //extend in x direction only
                {
                    cells[6*p+0] = op+p1; cells[6*p+1] = op+ndx*ndy+j; cells[6*p+2] = op+ndx*ndy+j+1;
                    cells[6*p+3] = op+p1; cells[6*p+4] = op+ndx*ndy+j+1; cells[6*p+5] = op+p1+1;
                }
            }
    }
    else if( posx == npx-1 && posy == npy-1 )
    {
        lc = (ndx-1)*(ndy-1)*2;
        oc = posy*(ndx*npx-1)*ndy*2 + posx*ndx*(ndy-1)*2;
        for (i = 0; i < ndx-1; i++)
            for (j = 0; j < ndy-1; j++)
            {
                int p = i*(ndy-1)+j;
                int p1 = i*ndy+j;
                cells[6*p+0] = op+p1; cells[6*p+1] = op+p1+ndy; cells[6*p+2] = op+p1+ndy+1;
                cells[6*p+3] = op+p1; cells[6*p+4] =op+p1+ndy+1; cells[6*p+5] = op+p1+1;
            }
    }
    else
    {
        lc = ndx*ndy*2;
        oc = posx*ndx*ndy*2 + posy*(ndx*npx-1)*ndy*2;
        for (i = 0; i < ndx; i++)
            for (j = 0; j < ndy; j++)
            {
                int p = i*ndy+j;
                if( i<ndx-1 && j<ndy-1 )
                {
                    cells[6*p+0] = op+p; cells[6*p+1] = op+p+ndy; cells[6*p+2] = op+p+ndy+1;
                    cells[6*p+3] = op+p; cells[6*p+4] =op+p+ndy+1; cells[6*p+5] = op+p+1;
                }
                else if( i==ndx-1 && j<ndy-1 )
                {
                    cells[6*p+0] = op+p; cells[6*p+1] = op+ndx*ndy+j; cells[6*p+2] = op+ndx*ndy+j+1;
                    cells[6*p+3] = op+p; cells[6*p+4] = op+ndx*ndy+j+1; cells[6*p+5] = op+p+1;
                }
                else if( i<ndx-1 && j==ndy-1 )
                {
                    cells[6*p+0] = op+p; cells[6*p+1] = op+p+ndy; cells[6*p+2] = op+nx_global*ndy+(i+1)*ndy;
                    cells[6*p+3] = op+p; cells[6*p+4] = op+nx_global*ndy+(i+1)*ndy; cells[6*p+5] = op+nx_global*ndy+i*ndy;
                }
                else // inter corner point
                {
                    cells[6*p+0] = op+p; cells[6*p+1] = op+ndx*ndy+j; cells[6*p+2] = op+nx_global*ndy+ndx*ndy;
                    cells[6*p+3] = op+p; cells[6*p+4] = op+nx_global*ndy+ndx*ndy; cells[6*p+5] = op+nx_global*ndy+i*ndy;
                }
            }
    }

    for (i=0; i<lc; i++)
        C[i] = 1.0*rank;    
 
    char * schema_version = "1.1";
    //char * dimemsions = "nx_global,ny_global";
 
	adios_init_noxml (comm);
    adios_set_max_buffer_size (50);

    adios_declare_group (&m_adios_group, "tri2d", "", adios_flag_yes);
    adios_select_method (m_adios_group, "MPI", "", "");

    adios_define_var (m_adios_group, "nx_global"
			,"", adios_integer
			,0, 0, 0);
    adios_define_var (m_adios_group, "ny_global"
            ,"", adios_integer
            ,0, 0, 0);
    adios_define_var (m_adios_group, "nproc"
                ,"", adios_integer                
                ,0, 0, 0);
    adios_define_var (m_adios_group, "npoints"
                ,"", adios_integer
                ,0, 0, 0);
    adios_define_var (m_adios_group, "num_cells"
                ,"", adios_integer
                ,0, 0, 0);
    adios_define_var (m_adios_group, "offs_x"
                ,"", adios_integer
                ,0, 0, 0);
    adios_define_var (m_adios_group, "offs_y"
                ,"", adios_integer
                ,0, 0, 0);
    adios_define_var (m_adios_group, "nx_local"
                ,"", adios_integer
                ,0, 0, 0);
    adios_define_var (m_adios_group, "ny_local"
                ,"", adios_integer
                ,0, 0, 0);
    adios_define_var (m_adios_group, "op"
                ,"", adios_integer
                ,0, 0, 0);
    adios_define_var (m_adios_group, "lp"
                ,"", adios_integer
                ,0, 0, 0);
    adios_define_var (m_adios_group, "oc"
                ,"", adios_integer
                ,0, 0, 0);
    adios_define_var (m_adios_group, "lc"
                ,"", adios_integer
                ,0, 0, 0);
    adios_define_var (m_adios_group, "points"
                    ,"", adios_double
                    ,"lp,2", "npoints,2", "op,0");
    adios_define_var (m_adios_group, "cells"
                    ,"", adios_integer
                    ,"lc,3", "num_cells,3", "oc,0");
    adios_define_var (m_adios_group, "N"
                    ,"", adios_double
                    ,"lp", "npoints", "op");
    adios_define_var (m_adios_group, "C"
                    ,"", adios_double
                    ,"lc", "num_cells", "oc");

    adios_define_attribute (m_adios_group, "description", "/nproc", adios_string, "Number of writers", "");
    adios_define_attribute (m_adios_group, "description", "/npoints", adios_string, "Number of points", "");
    adios_define_attribute (m_adios_group, "description", "/num_cells", adios_string, "Number of triangles", "");

    adios_define_schema_version (m_adios_group, schema_version);
    adios_define_mesh_unstructured ("points", "cells", "num_cells", "triangle", 0, "2", m_adios_group, "trimesh");
    adios_define_mesh_timevarying ("no", m_adios_group, "trimesh");

    adios_define_var_mesh (m_adios_group, "N", "trimesh");
    adios_define_var_centering (m_adios_group, "N", "point");
    adios_define_attribute (m_adios_group, "description", "/N", adios_string, "Node centered data", "");
    adios_define_var_mesh (m_adios_group, "C", "trimesh");
    adios_define_var_centering (m_adios_group, "C", "cell");
    adios_define_attribute (m_adios_group, "description", "/C", adios_string, "Cell centered data", "");

    adios_open (&adios_handle, "tri2d", "tri2d_noxml.bp", "w", comm);

    adios_groupsize = 13*sizeof(int) \
    + sizeof(double) * (lp) \
    + sizeof(double) * (lc) \
    + sizeof(double) * (lp) * (2) \
    + sizeof(int) * (lc) * (3);

    adios_group_size (adios_handle, adios_groupsize, &adios_totalsize);
    adios_write (adios_handle, "nproc", &nproc);
    adios_write (adios_handle, "npoints", &npoints);
    adios_write (adios_handle, "num_cells", &num_cells);
    adios_write (adios_handle, "nx_global", &nx_global);
    adios_write (adios_handle, "ny_global", &ny_global);
    adios_write (adios_handle, "offs_x", &offs_x);
    adios_write (adios_handle, "offs_y", &offs_y);
    adios_write (adios_handle, "nx_local", &nx_local);
    adios_write (adios_handle, "ny_local", &ny_local);
    adios_write (adios_handle, "lp", &lp);
    adios_write (adios_handle, "op", &op);
    adios_write (adios_handle, "lc", &lc);
    adios_write (adios_handle, "oc", &oc);
    adios_write (adios_handle, "N", N);
    adios_write (adios_handle, "C", C);
    adios_write (adios_handle, "points", points);
    adios_write (adios_handle, "cells", cells);

    adios_close (adios_handle);

    MPI_Barrier (comm);

    free (N);
    free (points);
    free (C);
    free (cells);

	adios_finalize (rank);

	MPI_Finalize ();
	return 0;
}
