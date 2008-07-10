typedef int MPI_COMM;
typedef int MPI_DATATYPE;

#define MPI_SUCCESS 0
#define MPI_FAILURE 1

#define MPI_COMM_WORLD 1

#define MPI_LOGICAL  0
#define MPI_INTEGER  1
#define MPI_REAL     2
#define MPI_DOUBLE_PRECISION  3
#define MPI_2DOUBLE_PRECISION  4
#define MPI_CHARACTER 5
#define MPI_BYTE 5

/*
 * you might not think you'd need this
 * many types, but in fact you do; PARAMESH
 * will routinelly generated hundreds of thousands
 * of types.
 */
#define MPI_MAXTYPES 100000

typedef int MPI_OP;
typedef int MPI_REQUEST;
typedef int MPI_STATUS;

#define MPI_MinLoc  1
#define MPI_MaxLoc  2
#define MPI_MIN     3  
#define MPI_MAX     4
#define MPI_SUM     5
#define MPI_LOR     6

#define MPI_ANY_SOURCE 999
#define MPI_ANY_TAG    999
#define MPI_PROC_NULL  999

typedef int MPI_TAG;
