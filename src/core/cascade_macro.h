/* Code snippets to run some code cascaded by one after the other MPI processes */
/* Usage example:

    int         rank, size;
    MPI_Comm    comm = MPI_COMM_WORLD;
    ...
    MPI_Init (&argc, &argv);
    MPI_Comm_rank (comm, &rank);
    MPI_Comm_size (comm, &size);
    ...

    CASCADE_BEGIN(comm,rank,size)
    printf ("rank %d: doing some work alone...\n", rank);
    ...
    CASCADE_END(comm,rank,size)
*/
#ifndef ADIOS_CASCADE_MACRO_H
#define ADIOS_CASCADE_MACRO_H


#define CASCADE_BEGIN(comm,rank,nproc) {            \
    int cascade_token=1;                            \
    MPI_Status status;                              \
    MPI_Request request;                            \
    if (rank > 0) {                                 \
        MPI_Recv (&cascade_token, 1, MPI_INT,       \
                  rank-1, rank-1, comm, &status);   \
    } 

/* Note that the opening bracket is not closed, only at the end of the next macro */


#define CASCADE_END(comm,rank,nproc)                \
    if (rank<nproc-1) {                             \
        MPI_Isend (&cascade_token, 1, MPI_INT,      \
                   rank+1, rank, comm, &request);   \
    }                                               \
    MPI_Barrier(comm);                              \
} /* closes the opening bracket of CASCADE_START */



#endif /*ADIOS_CASCADE_MACRO_H*/
