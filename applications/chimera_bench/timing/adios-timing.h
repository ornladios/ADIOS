#ifndef ADIOS_TIMING_H
#define ADIOS_TIMING_H 1

#include "mpi.h"

#define MAX_IO_COUNT 20
#define MAX_PRINT_BUF_SIZE 1*1024*1024
#define MAX_CYCLE_COUNT 400

/*
 * initialize profiling 
 *
 */
int init_prof(char *prof_file_name
             ,int max_io_count
             ,int max_cycle_count
			 ,MPI_Comm comm
			 ,int comm_size
			 ,int my_rank
			 );

/*
 * Fortran interface for init_prof()
 *
 */
int init_prof_(char *prof_file_name
              ,int *max_io_count
              ,int *max_cycle_count
              ,MPI_Fint *comm
              ,int *comm_size
              ,int *my_rank
              ,int *pfile_name_size
              );

/*
 * record start time of a simulation cycle
 *
 * Fortran interface 
 */
void cycle_start_(int *cycle);

/*
 * record end time of a simulation cycle
 *
 * Fortran interface 
 */
void cycle_end_(int *cycle);

/*
 * record start time of open
 *
 * Fortran interface 
 */
void open_start_(int *cycle);

/*
 * record end time of open
 *
 * Fortran interface 
 */
void open_end_(int *cycle);

/*
 * record start time of write
 *
 * Fortran interface 
 */
void write_start_(int *cycle);

/*
 * record end time of write
 *
 * Fortran interface 
 */
void write_end_(int *cycle);

/*
 * record start time of close
 *
 * Fortran interface 
 */
void close_start_(int *cycle);

/*
 * record end time of close
 *
 * Fortran interface 
 */
void close_end_(int *cycle);

/*
 * report profling results 
 *
 * Fortran interface 
 */
int finalize_prof_();

/*
 * get local timestamp
 */
double get_timestamp(); 

/* 
 * perform allreduces necessary to gather statistics on a single double 
 * value recorded across a comm.
 * 
 * code borrowed from mpi-tile-io package
 */
void get_time_dist(MPI_Comm comm
                  ,double *time
				  ,double *min
				  ,double *max
				  ,double *mean
				  ,double *var
				  );

/*
 * initialize profiling 
 *
 * used inside ADIOS so we can record io timing info for each every adios group
 */
int init_prof_internal(char *prof_file_name, struct adios_group_list_struct *group_list);

/*
 * initialize profiling 
 *
 * Fortran interface
 */
int init_prof_all_(char *prof_file_name, int prof_file_name_size);

/*
 * record open start time for specified group
 *
 * Fortran interface
 */
void open_start_for_group_(long long *gp_prof_handle, char *group_name, int *cycle, int *gp_name_size);

/*
 * record open end time for specified group
 *
 * Fortran interface
 */
void open_end_for_group_(long long *gp_prof_handle, int *cycle);

/*
 * record write start time for specified group
 *
 * Fortran interface
 */
void write_start_for_group_(long long *gp_prof_handle, int *cycle);

/*
 * record write end time for specified group
 *
 * Fortran interface
 */
void write_end_for_group_(long long *gp_prof_handle, int *cycle);

/*
 * record close start time for specified group
 *
 * Fortran interface
 */
void close_start_for_group_(long long *gp_prof_handle, int *cycle);

/*
 * record close end time for specified group
 *
 * Fortran interface
 */
void close_end_for_group_(long long *gp_prof_handle, int *cycle);

/*
 * Report timing info for all groups
 *
 * Fortran interface  
 */
int finalize_prof_all_();

#endif

