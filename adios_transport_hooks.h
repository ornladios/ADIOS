// this is defined in the lint program to get empty implementations
#ifdef ADIOS_EMPTY_TRANSPORTS
#define FORWARD_DECLARE(a) \
void adios_##a##_init (const char * parameters \
                      ,struct adios_method_struct * method \
                      ) {} \
void adios_##a##_open (struct adios_file_struct * fd \
                      ,struct adios_method_struct * method \
                      ) {} \
void adios_##a##_write (struct adios_file_struct * fd \
                       ,struct adios_var_struct * v \
                       ,void * data \
                       ,struct adios_method_struct * method \
                       ) {} \
void adios_##a##_get_write_buffer (struct adios_file_struct * fd \
                                  ,struct adios_var_struct * v \
                                  ,unsigned long long * size \
                                  ,void ** buffer \
                                  ,struct adios_method_struct * method \
                                  ) {} \
void adios_##a##_read (struct adios_file_struct * fd \
                      ,struct adios_var_struct * v \
                      ,void * buffer \
                      ,struct adios_method_struct * method \
                      ) {} \
void adios_##a##_close (struct adios_file_struct * fd \
                       ,struct adios_method_struct * method \
                       ) {} \
void adios_##a##_finalize (int mype, struct adios_method_struct * method) {} \
void adios_##a##_end_iteration (struct adios_method_struct * method) {} \
void adios_##a##_start_calculation (struct adios_method_struct * method) {} \
void adios_##a##_stop_calculation (struct adios_method_struct * method) {}
#else
#define FORWARD_DECLARE(a) \
void adios_##a##_init (const char * parameters \
                      ,struct adios_method_struct * method \
                      ); \
void adios_##a##_open (struct adios_file_struct * fd \
                      ,struct adios_method_struct * method \
                      ); \
void adios_##a##_write (struct adios_file_struct * fd \
                       ,struct adios_var_struct * v \
                       ,void * data \
                       ,struct adios_method_struct * method \
                       ); \
void adios_##a##_get_write_buffer (struct adios_file_struct * fd \
                                  ,struct adios_var_struct * v \
                                  ,unsigned long long * size \
                                  ,void ** buffer \
                                  ,struct adios_method_struct * method \
                                  ); \
void adios_##a##_read (struct adios_file_struct * fd \
                      ,struct adios_var_struct * v \
                      ,void * buffer \
                      ,struct adios_method_struct * method \
                      ); \
void adios_##a##_close (struct adios_file_struct * fd \
                       ,struct adios_method_struct * method \
                       ); \
void adios_##a##_finalize (int mype, struct adios_method_struct * method); \
void adios_##a##_end_iteration (struct adios_method_struct * method); \
void adios_##a##_start_calculation (struct adios_method_struct * method); \
void adios_##a##_stop_calculation (struct adios_method_struct * method);
#endif

#define MATCH_STRING_TO_METHOD(b,d) if (!strcmp (buf,b)) {*method=d;return 1;}

#define ASSIGN_FNS(a,b) \
(*t) [b].adios_init_fn = adios_##a##_init; \
(*t) [b].adios_open_fn = adios_##a##_open; \
(*t) [b].adios_write_fn = adios_##a##_write; \
(*t) [b].adios_get_write_buffer_fn = adios_##a##_get_write_buffer; \
(*t) [b].adios_read_fn = adios_##a##_read; \
(*t) [b].adios_close_fn = adios_##a##_close; \
(*t) [b].adios_finalize_fn = adios_##a##_finalize; \
(*t) [b].adios_end_iteration_fn = adios_##a##_end_iteration; \
(*t) [b].adios_start_calculation_fn = adios_##a##_start_calculation; \
(*t) [b].adios_stop_calculation_fn = adios_##a##_stop_calculation;

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
//// SETUP YOUR NEW TRANSPORT METHODS BELOW (FOLLOW THE PATTERN):          ////
//// 1. Add an entry to the ADIOS_IO_METHOD updating the ADIOS_METHOD_COUNT////
//// 2. Add a FOWARD_DECLARE line (assuming standard naming)               ////
//// 3. Add an entry to ADIOS_PARSE_METHOD_SETUP for the string and ID     ////
//// 4. Add an entry to ADIOS_INIT_TRANSPORTS_SETUP for name to ID         ////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

struct adios_method_struct;
struct adios_file_struct;
struct adios_var_struct;

// the list of the methods that have been integrated
enum ADIOS_IO_METHOD {ADIOS_METHOD_UNKNOWN     = -2
                     ,ADIOS_METHOD_NULL        = -1
                     ,ADIOS_METHOD_MPI         = 0
                     ,ADIOS_METHOD_DATATAP     = 1
                     ,ADIOS_METHOD_POSIX       = 2
                     ,ADIOS_METHOD_DART        = 3
                     ,ADIOS_METHOD_VTK         = 4
                     ,ADIOS_METHOD_POSIX_ASCII = 5
                     ,ADIOS_METHOD_MPI_CIO     = 6
                     ,ADIOS_METHOD_COUNT       = 7
                     };

// forward declare the functions (or dummies for internals use)
FORWARD_DECLARE(mpi)
FORWARD_DECLARE(mpi_cio)
FORWARD_DECLARE(posix)
FORWARD_DECLARE(vtk)
FORWARD_DECLARE(posix_ascii)
#if USE_PORTALS
FORWARD_DECLARE(datatap)
FORWARD_DECLARE(dart)
#endif

// add the string<->ID mapping here (also add ID in adios_internals.h)
#if USE_PORTALS
#define ADIOS_PARSE_METHOD_SETUP \
    MATCH_STRING_TO_METHOD("MPI",ADIOS_METHOD_MPI)                 \
    MATCH_STRING_TO_METHOD("DATATAP",ADIOS_METHOD_DATATAP)         \
    MATCH_STRING_TO_METHOD("PBIO",ADIOS_METHOD_DATATAP)            \
    MATCH_STRING_TO_METHOD("POSIX",ADIOS_METHOD_POSIX)             \
    MATCH_STRING_TO_METHOD("FB",ADIOS_METHOD_POSIX)                \
    MATCH_STRING_TO_METHOD("DART",ADIOS_METHOD_DART)               \
    MATCH_STRING_TO_METHOD("VTK",ADIOS_METHOD_VTK)                 \
    MATCH_STRING_TO_METHOD("POSIX_ASCII",ADIOS_METHOD_POSIX_ASCII) \
    MATCH_STRING_TO_METHOD("MPI_CIO",ADIOS_METHOD_MPI_CIO) \
    MATCH_STRING_TO_METHOD("NULL",ADIOS_METHOD_NULL)
#else
#define ADIOS_PARSE_METHOD_SETUP \
    MATCH_STRING_TO_METHOD("MPI",ADIOS_METHOD_MPI)                 \
    MATCH_STRING_TO_METHOD("POSIX",ADIOS_METHOD_POSIX)             \
    MATCH_STRING_TO_METHOD("FB",ADIOS_METHOD_POSIX)                \
    MATCH_STRING_TO_METHOD("VTK",ADIOS_METHOD_VTK)                 \
    MATCH_STRING_TO_METHOD("POSIX_ASCII",ADIOS_METHOD_POSIX_ASCII) \
    MATCH_STRING_TO_METHOD("NULL",ADIOS_METHOD_NULL)
#endif

// add the initialization of the functions for the calls here
#if USE_PORTALS
#define ADIOS_INIT_TRANSPORTS_SETUP \
    ASSIGN_FNS(mpi,ADIOS_METHOD_MPI)                 \
    ASSIGN_FNS(mpi_cio,ADIOS_METHOD_MPI_CIO)         \
    ASSIGN_FNS(posix,ADIOS_METHOD_POSIX)             \
    ASSIGN_FNS(datatap,ADIOS_METHOD_DATATAP)         \
    ASSIGN_FNS(dart,ADIOS_METHOD_DART)               \
    ASSIGN_FNS(vtk,ADIOS_METHOD_VTK)                 \
    ASSIGN_FNS(posix_ascii,ADIOS_METHOD_POSIX_ASCII)
#else
#define ADIOS_INIT_TRANSPORTS_SETUP \
    ASSIGN_FNS(mpi,ADIOS_METHOD_MPI)                 \
    ASSIGN_FNS(mpi_cio,ADIOS_METHOD_MPI_CIO)         \
    ASSIGN_FNS(posix,ADIOS_METHOD_POSIX)             \
    ASSIGN_FNS(vtk,ADIOS_METHOD_VTK)                 \
    ASSIGN_FNS(posix_ascii,ADIOS_METHOD_POSIX_ASCII)
#endif
