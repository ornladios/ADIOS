/* 
 * ADIOS is freely available under the terms of the BSD license described
 * in the COPYING file in the top level directory of this source distribution.
 *
 * Copyright (c) 2008 - 2009.  UT-BATTELLE, LLC. All rights reserved.
 */

/* 
   ADIOS Staging server main function.
*/
#define __MAIN__
#include "globals.h"
#include "mem.h"
#include "options.h"
#if HAVE_PTHREAD
#include "worker_thread.h"
#endif
#include <string.h> // memset

/* Two functions some other files should define */
void transport_thread_main (void *arg);
void transport_finalize (void);

int main (int argc, char ** argv) 
{
    memset(&gd, 0, sizeof(struct global_data_struct));
    gd.mpi_comm = MPI_COMM_WORLD;
    gd.application_finished = false;
    gd.nc_current = -1;
    gd.terminate = false;

    MPI_Init (&argc, &argv);
    MPI_Comm_rank (gd.mpi_comm, &gd.mpi_rank);
    MPI_Comm_size (gd.mpi_comm, &gd.mpi_size);
    /*
    int x = pthread_mutex_init(&gd.mpi_mutex, NULL);
    if (x != 0) {
        fprintf(stderr, "ERROR: mpi mutex initialization failed.");
        exit(1);
    }
    */

    /* set options to default then get user options */
    user_max_block_size = MAX_BLOCK_SIZE;
    user_max_client_pulls = MAX_CLIENT_PULLS;
    user_max_memory_allowed = 0;
    logfile_separate_ranks = false;
    verbose = 0;
    logpath = NULL;
    int exitcode = options_process_args(argc, argv);
    if (exitcode) {
        if (exitcode==255) exitcode=0; // e.g. --help, exits with 0
        goto terminate;
    }

    ardma_logger_init(logpath, verbose, (logfile_separate_ranks ? gd.mpi_rank : -1)); 

    gd.wtq = queue_init(QUEUE_UNBOUNDED);
    gd.twq = queue_init(QUEUE_UNBOUNDED);
    if (!gd.wtq || !gd.twq) {
        log_error("rank %d: ERROR: could not create transport-worker queues\n", gd.mpi_rank);
        MPI_Abort(gd.mpi_comm, 1);
    }

    worker_global_init();

#if HAVE_PTHREAD
    /* Start worker thread.  It may spawn new threads on its own */
    pthread_attr_t attr;
    pthread_attr_init(&attr);
    pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
    int rc = pthread_create(&gd.worker_thread, &attr, worker_thread_main, NULL); 
    if (rc) {
        log_error("ERROR: rank %d: Cannot create worker thread, err code = %d\n", gd.mpi_rank, rc);
        MPI_Abort(gd.mpi_comm, 1);
    }
    pthread_attr_destroy(&attr);
#endif

    /* The main thread becomes the transport thread. The worker thread will be called as 
       function if no threading is used*/
    transport_thread_main(NULL);

#if HAVE_PTHREAD
    /* function above exits when application disappears, so do the worker thread */
    void *status;
    rc = pthread_join(gd.worker_thread, &status);
    if (rc) {
        log_error("ERROR: rank %d: Cannot join worker thread, err code = %d\n", gd.mpi_rank, rc);
    } else {
        log_debug("rank %d: Joined worker thread.\n", gd.mpi_rank);
    }
#endif

    queue_destroy(gd.wtq);
    queue_destroy(gd.twq);
    transport_finalize();
    ardma_logger_finalize();

terminate:
    //pthread_mutex_lock(&gd.mpi_mutex);
    log_debug("rank %d: call MPI_Barrier.\n", gd.mpi_rank);
    MPI_Barrier (gd.mpi_comm);
    log_debug("rank %d: call MPI_Finalize.\n", gd.mpi_rank);
    MPI_Finalize ();
    //pthread_mutex_unlock(&gd.mpi_mutex);
    //pthread_mutex_destroy(&gd.mpi_mutex);
    log_debug("rank %d: exit.\n", gd.mpi_rank);
    return exitcode;
}
