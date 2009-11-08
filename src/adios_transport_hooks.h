#ifndef ADIOS_TRANSPORT_HOOKS_H
#define ADIOS_TRANSPORT_HOOKS_H

#include "config.h"
#include <stdint.h>
#include <string.h>

// this is defined in the lint program to get empty implementations
#ifdef ADIOS_EMPTY_TRANSPORTS
#define FORWARD_DECLARE(a) \
void adios_##a##_init (const char * parameters \
                      ,struct adios_method_struct * method \
                      ) {} \
int adios_##a##_open (struct adios_file_struct * fd \
                     ,struct adios_method_struct * method, void * comm \
                     ) {return 0;} \
enum ADIOS_FLAG adios_##a##_should_buffer (struct adios_file_struct * fd \
                                          ,struct adios_method_struct * method \
                                          ) {return 0;} \
void adios_##a##_write (struct adios_file_struct * fd \
                       ,struct adios_var_struct * v \
                       ,void * data \
                       ,struct adios_method_struct * method \
                       ) {} \
void adios_##a##_get_write_buffer (struct adios_file_struct * fd \
                                  ,struct adios_var_struct * v \
                                  ,uint64_t * size \
                                  ,void ** buffer \
                                  ,struct adios_method_struct * method \
                                  ) {} \
void adios_##a##_read (struct adios_file_struct * fd \
                      ,struct adios_var_struct * v \
                      ,void * buffer \
                      ,uint64_t buffer_size \
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
int adios_##a##_open (struct adios_file_struct * fd \
                     ,struct adios_method_struct * method, void * comm \
                     ); \
enum ADIOS_FLAG adios_##a##_should_buffer (struct adios_file_struct * fd \
                                          ,struct adios_method_struct * method \
                                          ); \
void adios_##a##_write (struct adios_file_struct * fd \
                       ,struct adios_var_struct * v \
                       ,void * data \
                       ,struct adios_method_struct * method \
                       ); \
void adios_##a##_get_write_buffer (struct adios_file_struct * fd \
                                  ,struct adios_var_struct * v \
                                  ,uint64_t * size \
                                  ,void ** buffer \
                                  ,struct adios_method_struct * method \
                                  ); \
void adios_##a##_read (struct adios_file_struct * fd \
                      ,struct adios_var_struct * v \
                      ,void * buffer \
                      ,uint64_t buffer_size \
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
// VTK and POSIX_ASCII are placeholders reserved for future use
enum ADIOS_IO_METHOD {ADIOS_METHOD_UNKNOWN     = -2
                     ,ADIOS_METHOD_NULL        = -1
                     ,ADIOS_METHOD_MPI         = 0
                     ,ADIOS_METHOD_DATATAP     = 1
                     ,ADIOS_METHOD_POSIX       = 2
                     ,ADIOS_METHOD_DART        = 3
                     ,ADIOS_METHOD_VTK         = 4
                     ,ADIOS_METHOD_POSIX_ASCII = 5
                     ,ADIOS_METHOD_MPI_CIO     = 6
                     ,ADIOS_METHOD_PHDF5       = 7
                     ,ADIOS_METHOD_PROVENANCE  = 8
                     ,ADIOS_METHOD_MPI_STRIPE  = 9
                     ,ADIOS_METHOD_MPI_STRIPE2 = 10
                     ,ADIOS_METHOD_MPI_STAGGER = 11
                     ,ADIOS_METHOD_MPI_AGG     = 12
                     ,ADIOS_METHOD_ADAPTIVE    = 13
                     ,ADIOS_METHOD_COUNT       = 14
                     };

// forward declare the functions (or dummies for internals use)
FORWARD_DECLARE(mpi)
FORWARD_DECLARE(mpi_stripe2)
FORWARD_DECLARE(mpi_cio)
FORWARD_DECLARE(mpi_stripe)
FORWARD_DECLARE(mpi_stagger)
FORWARD_DECLARE(mpi_aggregate)
FORWARD_DECLARE(datatap)
FORWARD_DECLARE(posix)
FORWARD_DECLARE(provenance)
FORWARD_DECLARE(phdf5)
FORWARD_DECLARE(adaptive)

#if USE_PORTALS
FORWARD_DECLARE(dart)
#endif

#endif
