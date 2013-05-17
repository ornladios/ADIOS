#include "skel_xml_output.h"
#include <stdio.h>
#include <stdlib.h>
#include "mpi.h"


int common_skel_write_coarse_xml_data (double open_time, double write_time, double close_time, double total_time);


int FC_FUNC_(skel_write_coarse_xml_data_f, SKEL_WRITE_COARSE_XML_DATA_F) (double open_time, double write_time, double close_time, double total_time)
{
    return common_skel_write_coarse_xml_data (open_time, write_time, close_time, total_time);
}


int skel_write_coarse_xml_data (double open_time, double write_time, double close_time, double total_time)
{
    return common_skel_write_coarse_xml_data (open_time, write_time, close_time, total_time);
}


/*
 * Used by a generated skeletal program to write the times gathered by the skeletal to the xml file.
 *
 * See also: adios_timing_write_xml_common() in src/core/adios_timing.c
 *
 */
int common_skel_write_coarse_xml_data (double open_time, double write_time, double close_time, double total_time)
{
    char* filename = "gts_skel_time.xml";

    int size, rank, i;
    MPI_Comm_size (MPI_COMM_WORLD, &size);
    MPI_Comm_rank (MPI_COMM_WORLD, &rank);

    //double close_time = 1.67f * rank;
    double * close_times = NULL;

    if (rank == 0)
    {
        close_times = (double*) malloc (sizeof (double) * size);
    }

    // Collect timing on proc 0

    MPI_Gather (
        &close_time,
        1,  // sendcount
        MPI_DOUBLE, // sendtype
        close_times,
        1, // recvcount
        MPI_DOUBLE, // recvtype
        0, // root
        MPI_COMM_WORLD
    );

    // Write it out (proc 0 only)
    if (rank == 0)
    {
        FILE* f = fopen (filename, "a");

        char* key_str = "ad_close";

        fprintf (f, "<skel_result>\n");
        fprintf (f, "  <adios_timing cores='%i' keys='%s'>\n", size, key_str);

        for (i = 0; i < size; i++)
        {
            fprintf (f, "    <proc id='%i' vals='%f' />\n", i, close_times[i]);
        }

        fprintf (f, "  </adios_timing>\n");
        fprintf (f, "</skel_result>\n");

        free (close_times);
    }
    return 0;
}
