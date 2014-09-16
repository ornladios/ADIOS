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
#if defined SKEL_TIMING && !defined _NOMPI //No timing information on single process

    struct adios_file_struct * fd = (struct adios_file_struct *) fd_p;
    if (!fd)
    {
        adios_error (err_invalid_file_pointer,
                     "Invalid handle passed to adios_get_timing_name\n");
        return;
    }

    if (!fd->group || !fd->group->prev_timing_obj)
    {
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

    count_to_send = (fd->group->prev_timing_obj->event_count > ADIOS_TIMING_MAX_EVENTS) ?
                      ADIOS_TIMING_MAX_EVENTS : fd->group->prev_timing_obj->event_count;


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
    MPI_Aint disps[]      = {0,2*sizeof(int)};
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
        &fd->group->prev_timing_obj->events, // sendbuf
        count_to_send, // sendcount
        event_type, // sendtype
        events, //recvbuf
        counts, // recvcounts
        displs, // displacements
        event_type, // recvtype
        0, // root
        MPI_COMM_WORLD // comm
    );

    // Gather the write sizes
    int *write_sizes = NULL;
    if (rank == 0)
    {
        write_sizes = (int*) malloc (sizeof(int) * size);
    }

    MPI_Gather (
        &fd->write_size_bytes, //sendbuf
        1, //sendcount
        MPI_INT, //sendtype
        write_sizes, //recvbuf
        1, //recvcount
        MPI_INT, //recvtype
        0, //root
        MPI_COMM_WORLD //comm
    );

    // Write the events to a file
    if (rank == 0)
    {
        FILE* f = fopen (filename, "w");
        int event_rank;

        for (i = 0; i < size; i++)
        {
            fprintf (f, "'%i'%i\n", i, write_sizes[i]);
        }

        // Write the labels
        for (i = 0; i < fd->group->prev_timing_obj->internal_count; i++)
        {
            fprintf (f, ":%i:%s\n", ADIOS_TIMING_MAX_USER_TIMERS + i,
                     fd->group->prev_timing_obj->names[ADIOS_TIMING_MAX_USER_TIMERS + i]); 
        }

        // Now the event data
        i = 0;
        for (event_rank = 0; event_rank < size; event_rank++)
        {
            for ( ; i < displs[event_rank] + counts[event_rank]; i++) 
            {
                fprintf (f, "%i,%i%s,%f\n", event_rank, events[i].type,
                         events[i].is_start?"S":"E", events[i].time);
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

    if (! fd->group->prev_timing_obj)
    {
        *tc = 0;
    }
    else
    {
        *tc = fd->group->prev_timing_obj->internal_count;
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

    strcpy (name, fd->group->prev_timing_obj->names[index]);
    //*name = fd->group->prev_timing_obj->names[index];

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

    *value = fd->group->prev_timing_obj->times[index];

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
    ts->names = (char**) malloc ( (ADIOS_TIMING_MAX_USER_TIMERS + timer_count) * sizeof (char*) );
    ts->times = (double*) malloc ( (ADIOS_TIMING_MAX_USER_TIMERS + timer_count) * sizeof (double) );

    adios_clear_timers (ts);

    for (i = 0; i < timer_count; i++)
    {
        ts->names[ADIOS_TIMING_MAX_USER_TIMERS + i] = (char*) malloc (strlen(timer_names[i]) + 1 * sizeof (char) );
        strcpy (ts->names[ADIOS_TIMING_MAX_USER_TIMERS + i], timer_names[i]);
    }

    return ts;
}

void adios_clear_timers (struct adios_timing_struct * ts)
{
    ts->user_count = 0;
    ts->event_count = 0;

    // Clear all timers
    memset(ts->times, 0, (ADIOS_TIMING_MAX_USER_TIMERS + ts->internal_count) * sizeof (double) );
    memset(ts->names, 0, (ADIOS_TIMING_MAX_USER_TIMERS + ts->internal_count) * sizeof (char*) );
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


