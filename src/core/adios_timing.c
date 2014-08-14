/*
 * ADIOS is freely available under the terms of the BSD license described
 * in the COPYING file in the top level directory of this source distribution.
 *
 * Copyright (c) 2008 - 2009.  UT-BATTELLE, LLC. All rights reserved.
 */


#include "public/adios.h"
#include "public/adios_error.h"
#include "core/adios_internals.h"
#include "core/adios_timing.h"
#include "core/adios_logger.h"
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>

#include "public/adios_mpi.h"
//#include "mpi.h"


/*
 * Dump the timing information to a file.
 * Called both from C and Fortran API's (adios.c and adiosf.c)
*/
void adios_timing_write_xml_common (int64_t fd_p, const char* filename)
{
#ifdef SKEL_TIMING

    struct adios_file_struct * fd = (struct adios_file_struct *) fd_p;
    if (!fd)
    {
        adios_error (err_invalid_file_pointer,
                     "Invalid handle passed to adios_get_timing_name\n");
        return;
    }

    if (!fd->timing_obj)
    {
        log_error ("No timing info available, file not written\n");
        // No timing info, don't write anything.
        return;
    }

    int size, rank, i, global_event_count, count_to_send;
    int * counts;
    int * displs;
    struct adios_timing_event_struct* events;
    MPI_Datatype event_type;
    MPI_Comm_size (MPI_COMM_WORLD, &size);
    MPI_Comm_rank (MPI_COMM_WORLD, &rank);

    if (rank == 0)
    {
        counts = (int*) malloc (sizeof (int) * size);
    }

    // Collect all of the events on proc 0
    // First, per proc event counts

    count_to_send = (fd->timing_obj->event_count > ADIOS_TIMING_MAX_EVENTS) ?
                      ADIOS_TIMING_MAX_EVENTS : fd->timing_obj->event_count;


    MPI_Gather (
        &count_to_send, // sendbuf
        1,              // sendcount
        MPI_INT,        // sendtype
        counts,         // recvbuf
        1,           // recvcount
        MPI_INT,        // recvtype
        0,              // root
        MPI_COMM_WORLD  // comm
    );

    if (rank == 0)
    {

        displs = (int*) malloc (sizeof (int) * size);
        displs[0] = 0;
        global_event_count = counts[0];

        for (i = 1; i < size; i++)
        {
            displs[i] = displs[i-1] + counts[i-1];
            global_event_count += counts[i];
        }

        events = (struct adios_timing_event_struct*) malloc (
            sizeof (struct adios_timing_event_struct) * global_event_count);
    }

    // structure of the adios_timing_event_struct (int, int, double)
    int blocklens[]  = {2,1};
    int disps[]      = {0,8};
    MPI_Datatype types[] = {MPI_INT,MPI_DOUBLE};

    MPI_Type_create_struct (
        2, // count
        blocklens, // array_of_blocklengths
        disps, // array_of_displacements
        types, // array_of_types
        &event_type
    );
    MPI_Type_commit (&event_type);

    // Now the events
    MPI_Gatherv (
        &fd->timing_obj->events, // sendbuf
        count_to_send, // sendcount
        event_type, // sendtype
        events, //recvbuf
        counts, // recvcounts
        displs, // displacements
        event_type, // recvtype
        0, // root
        MPI_COMM_WORLD // comm
    );

    // Write the events to a file
    if (rank == 0)
    {
        FILE* f = fopen (filename, "w");
        int event_rank;

        i = 0;
        for (event_rank = 0; event_rank < size; event_rank++)
        {
            for ( ; i < displs[event_rank] + counts[event_rank]; i++) 
            {
                fprintf (f, "%i,%i%s,%f\n", event_rank, events[i].type, events[i].is_start?"S":"E", events[i].time);
            }
        }

        fclose(f);
    }


    if (rank == 0)
    {
        if (counts)
            free (counts);
    }

#else
    log_warn ("Timing information is not currently available.\n"
              "To use the Skel timing functions, you must enable them when building ADIOS.\n"
              "Use --enable-skel-timing during the configuration step.\n");
#endif


// The old XML timing output
#if 0
    int size, rank, i, p;
    MPI_Comm_size (MPI_COMM_WORLD, &size);
    MPI_Comm_rank (MPI_COMM_WORLD, &rank);

    double * internal_times = NULL;
    double * user_times = NULL;

    // Allocate space for aggregation
    if (rank == 0)
    {
        internal_times = (double*) malloc (sizeof (double) *
                             fd->timing_obj->internal_count * size);
        user_times = (double*) malloc (sizeof (double) *
                             fd->timing_obj->user_count * size);
    }

    // Aggregate timing info on rank 0
    // Handle internal counts and user counts separately

    MPI_Gather (
        fd->timing_obj->times + ADIOS_TIMING_MAX_USER_TIMERS,
        fd->timing_obj->internal_count,  // sendcount
        MPI_DOUBLE, // sendtype
        internal_times,
        fd->timing_obj->internal_count, // recvcount
        MPI_DOUBLE, // recvtype
        0, // root
        MPI_COMM_WORLD
    );

    MPI_Gather (
        fd->timing_obj->times,  // sendbuf
        fd->timing_obj->user_count,  // sendcount
        MPI_DOUBLE, // sendtype
        user_times,
        fd->timing_obj->user_count, // recvcount
        MPI_DOUBLE, // recvtype
        0, // root
        MPI_COMM_WORLD
    );


    // Now write all timing info from rank 0

    if (rank == 0)
    {
        FILE* f = fopen (filename, "w");

        // Rank 0 starts the xml and includes the text labels
        fprintf (f, "<skel_result><adios_timing cores='%i' keys='", size);
        for (i = 0; i < fd->timing_obj->user_count; i++)
        {
            if (fd->timing_obj->names[i])
            {
                fprintf (f, "%s, ", fd->timing_obj->names[i]);
            }
            else
            {
                fprintf (f, "user%i, ", i);
            }
        }
        for (i = 0; i < fd->timing_obj->internal_count; i++)
        {
            if (fd->timing_obj->names[ADIOS_TIMING_MAX_USER_TIMERS + i])
            {
                fprintf (f, "%s", fd->timing_obj->names[ADIOS_TIMING_MAX_USER_TIMERS + i]);
            }
            else
            {
                fprintf (f, "internal%i", i);
            }
            if (i != fd->timing_obj->internal_count - 1) // Skip trailing comma
            {
                fprintf (f, ", ");
            }
        }
        fprintf (f, "' "); // Close the keys attribute

        // Assume there is only one method in play
        fprintf (f, "method='%s' ", fd->group->methods->method->method);

        struct timeval tv;
        gettimeofday (&tv, NULL);
        double time = tv.tv_sec+(tv.tv_usec/1000000.0);
        fprintf (f, "start_time='%f' ", time);

        fprintf (f, ">\n"); // Close the adios_timing element


// Use the aggregated values here

        for (p = 0; p < size; p++)
        {
            // This part should be the same for all procs
            fprintf (f, "<proc id='%i' vals='", p);
            for (i = 0; i < fd->timing_obj->user_count; i++)
            {
                fprintf (f, "%f, ", user_times[p*fd->timing_obj->user_count+i]);
                //fprintf (f, "%f, ", fd->timing_obj->times[i]);
            }
            for (i = 0; i < fd->timing_obj->internal_count; i++)
            {
                fprintf (f, "%f", internal_times[p*fd->timing_obj->internal_count+i]);
                //fprintf (f, "%f", fd->timing_obj->times[ADIOS_TIMING_MAX_USER_TIMERS + i]);
                if (i != fd->timing_obj->internal_count - 1)
                {
                    fprintf (f, ", ");
                }
            }
            fprintf (f, "' />\n");
        }


        fprintf (f, "</adios_timing></skel_result>\n");
        fclose (f);

    }
#endif

}



//Build the internal functions only when timing is enabled.
#ifdef SKEL_TIMING
int adios_get_timing_internal_count (int64_t fd_p, int64_t * tc)
{
    struct adios_file_struct * fd = (struct adios_file_struct *) fd_p;
    if (!fd)
    {
        adios_error (err_invalid_file_pointer,
                     "Invalid handle passed to adios_get_timing_count\n");

        return 1;
    }

    if (! fd->timing_obj)
    {
        *tc = 0;
    }
    else
    {
        *tc = fd->timing_obj->internal_count;
    }

    return 0;
}


int adios_get_timing_name (int64_t fd_p, int64_t index, char* name)
{
    struct adios_file_struct * fd = (struct adios_file_struct *) fd_p;
    if (!fd)
    {
        adios_error (err_invalid_file_pointer,
                     "Invalid handle passed to adios_get_timing_name\n");

        return 1;
    }

    strcpy (name, fd->timing_obj->names[index]);
    //*name = fd->timing_obj->names[index];

    return 0;
}


int adios_get_timing_value (int64_t fd_p, int64_t index, double* value)
{
    struct adios_file_struct * fd = (struct adios_file_struct *) fd_p;
    if (!fd)
    {
        adios_error (err_invalid_file_pointer,
                     "Invalid handle passed to adios_get_timing_value\n");

        return 1;
    }

    *value = fd->timing_obj->times[index];

    return 0;
}


void adios_timing_go (struct adios_timing_struct * ts, int64_t index)
{
    // Grab the time
    double now = MPI_Wtime();

    // Do accounting for time summary
    ts->times[index] -= now;

    // Log the event
    struct adios_timing_event_struct * new_event =
        &(ts->events[ts->event_count % ADIOS_TIMING_MAX_EVENTS]);
    new_event->type = index;
    new_event->is_start = 1;
    new_event->time = now;
    ts->event_count++;

}


void adios_timing_stop (struct adios_timing_struct * ts, int64_t index)
{
    // Grab the time
    double now = MPI_Wtime();

    // Do accounting for time summary
    ts->times[index] += now;

    // Log the event
    struct adios_timing_event_struct * new_event =
        &(ts->events[ts->event_count % ADIOS_TIMING_MAX_EVENTS]);

    new_event->type = index;
    new_event->is_start = 0;
    new_event->time = now;
    ts->event_count++;
}


struct adios_timing_struct *  adios_timing_create (int timer_count, char** timer_names)
{
    int i;
    struct adios_timing_struct * ts = (struct adios_timing_struct *)
                                       malloc (sizeof (struct adios_timing_struct) );

    ts->internal_count = timer_count;
    ts->user_count = 0;
    ts->names = (char**) malloc ( (ADIOS_TIMING_MAX_USER_TIMERS + timer_count) * sizeof (char*) );
    ts->times = (double*) malloc ( (ADIOS_TIMING_MAX_USER_TIMERS + timer_count) * sizeof (double) );
    ts->event_count = 0;


    // Clear all timers
    memset(ts->times, 0, (ADIOS_TIMING_MAX_USER_TIMERS + timer_count) * sizeof (double) );
    memset(ts->names, 0, (ADIOS_TIMING_MAX_USER_TIMERS + timer_count) * sizeof (char*) );

    for (i = 0; i < timer_count; i++)
    {
        ts->names[ADIOS_TIMING_MAX_USER_TIMERS + i] = (char*) malloc (strlen(timer_names[i]) + 1 * sizeof (char) );
        strcpy (ts->names[ADIOS_TIMING_MAX_USER_TIMERS + i], timer_names[i]);
    }

    return ts;
}



void adios_timing_destroy (struct adios_timing_struct * timing_obj)
{
    if (timing_obj)
    {

        if (timing_obj->times)
        {
            free (timing_obj->times);
        }
        free (timing_obj);
    }
}

#endif // ifdef SKEL_TIMING


