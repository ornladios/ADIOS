#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <sys/stat.h>
#include <math.h>
#include "mpi.h"

#include "adios.h"
#include "adios_transport_hooks.h"
#include "adios_internals.h"

#include "adios-timing.h"

/*
 * the following variables have been deprecated in
 * ADIOS; so just set them to use default MPI values 
 */
MPI_Comm adios_mpi_comm_world = MPI_COMM_WORLD;
MPI_Comm adios_mpi_comm_self = MPI_COMM_SELF;

typedef struct _cycle_timing_info {
    int cycle;
    double cycle_start_time;
    double cycle_end_time;
    double cycle_time;
} cycle_timing_info;

typedef struct _timing_info {
    int cycle;
    double open_start_time;
    double open_end_time;
    double open_time;
    double write_start_time;
    double write_end_time;
    double write_time;
    double close_start_time;
    double close_end_time;
    double close_time;
} timing_info;

typedef struct _group_prof_info {
    struct adios_group_struct *adios_group_id;
    MPI_Comm communicator;
    int comm_size;
    int myid;
    int max_io_count;
    int io_count;
    timing_info *timing; 
} group_prof_info;

typedef struct _io_prof_setting {
    char *prof_file_name;
    MPI_File prof_file_ptr;  
    int split_result;
    group_prof_info *gp_prof_info;
    int num_groups;
} io_prof_setting;

static io_prof_setting profile;

/*
 * simplified profiling facility:
 * suitable for profiling code involving only one adios group and with nothing
 * between successive adios_write statements
 *
 */
typedef struct _simple_io_prof_setting {
    MPI_Comm communicator;
    int comm_size;
    int myid;
    int max_io_count;
    int io_count;
    timing_info *timing; 
    int max_cycle_count;
    int cycle_count;
    cycle_timing_info *cycle_timing; 
    char *prof_file_name;
    FILE *prof_file;  
} simple_io_prof_setting;

static simple_io_prof_setting simple_profile;

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
             )
{

    simple_profile.communicator = comm;
    simple_profile.comm_size = comm_size;
    simple_profile.myid = my_rank;
    simple_profile.max_io_count = max_io_count;
    simple_profile.io_count = 0;
    simple_profile.timing = (timing_info *) malloc(sizeof(timing_info) * max_io_count);
    simple_profile.max_cycle_count = max_cycle_count;
    simple_profile.cycle_count = 0;
    simple_profile.cycle_timing = (cycle_timing_info *) malloc(sizeof(cycle_timing_info) * max_cycle_count);
   
    if(!prof_file_name || !strcmp(prof_file_name, "")) {
        fprintf(stderr, "ADIOS Profiling: using file \"./log\" to save profiling results.\n");    
        simple_profile.prof_file_name = "./log";          
    }
    else {
        simple_profile.prof_file_name = prof_file_name;
    }

    return 0;
}

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
              )
{
    return init_prof(prof_file_name, *max_io_count, *max_cycle_count, MPI_Comm_f2c(*comm), *comm_size, *my_rank);
}

/*
 * record start time of a simulation cycle
 *
 * Fortran interface 
 */
void cycle_start_(int *cycle)
{
    if(simple_profile.cycle_count < simple_profile.max_cycle_count) {
        simple_profile.cycle_timing[simple_profile.cycle_count].cycle = *cycle;
        simple_profile.cycle_timing[simple_profile.cycle_count].cycle_start_time = get_timestamp();
    }
    else {
        fprintf(stderr, "ADIOS Profiling: cycle count exceeds max limit (%d)\n", simple_profile.max_cycle_count);
        exit(-1);
    }
}

/*
 * record end time of a simulation cycle
 *
 * Fortran interface 
 */
void cycle_end_(int *cycle)
{
    simple_profile.cycle_timing[simple_profile.cycle_count].cycle = *cycle;
    simple_profile.cycle_timing[simple_profile.cycle_count].cycle_end_time = get_timestamp();
    simple_profile.cycle_timing[simple_profile.cycle_count].cycle_time = simple_profile.cycle_timing[simple_profile.cycle_count].cycle_end_time 
                                              - simple_profile.cycle_timing[simple_profile.cycle_count].cycle_start_time;
    simple_profile.cycle_count ++;
}

/*
 * record start time of open
 *
 * Fortran interface 
 */
void open_start_(int *cycle)
{
    if(simple_profile.io_count < simple_profile.max_io_count) {
        simple_profile.timing[simple_profile.io_count].cycle = *cycle;
        simple_profile.timing[simple_profile.io_count].open_start_time = get_timestamp();
    }
    else {
        fprintf(stderr, "ADIOS Profiling: io count exceeds max limit (%d)\n", simple_profile.max_io_count);
        exit(-1);
    }
}

/*
 * record end time of open
 *
 * Fortran interface 
 */
void open_end_(int *cycle)
{
    simple_profile.timing[simple_profile.io_count].cycle = *cycle;
    simple_profile.timing[simple_profile.io_count].open_end_time = get_timestamp();
    simple_profile.timing[simple_profile.io_count].open_time = simple_profile.timing[simple_profile.io_count].open_end_time 
                                              - simple_profile.timing[simple_profile.io_count].open_start_time;
}

/*
 * record start time of write
 *
 * Fortran interface 
 */
void write_start_(int *cycle)
{
    simple_profile.timing[simple_profile.io_count].cycle = *cycle;
    simple_profile.timing[simple_profile.io_count].write_start_time = get_timestamp();
}

/*
 * record end time of write
 *
 * Fortran interface 
 */
void write_end_(int *cycle)
{
    simple_profile.timing[simple_profile.io_count].cycle = *cycle;
    simple_profile.timing[simple_profile.io_count].write_end_time = get_timestamp();
    simple_profile.timing[simple_profile.io_count].write_time = simple_profile.timing[simple_profile.io_count].write_end_time 
                                              - simple_profile.timing[simple_profile.io_count].write_start_time;
}

/*
 * record start time of close
 *
 * Fortran interface 
 */
void close_start_(int *cycle)
{
    simple_profile.timing[simple_profile.io_count].cycle = *cycle;
    simple_profile.timing[simple_profile.io_count].close_start_time = get_timestamp();
}

/*
 * record end time of close
 *
 * Fortran interface 
 */
void close_end_(int *cycle)
{
    simple_profile.timing[simple_profile.io_count].cycle = *cycle;
    simple_profile.timing[simple_profile.io_count].close_end_time = get_timestamp();
    simple_profile.timing[simple_profile.io_count].close_time = simple_profile.timing[simple_profile.io_count].close_end_time 
                                              - simple_profile.timing[simple_profile.io_count].close_start_time;
    simple_profile.io_count ++;
}

/*
 * report profling results 
 *
 * Fortran interface 
 */
int finalize_prof_()
{
    int i;
    int rc;

    if(simple_profile.myid == 0) {
        simple_profile.prof_file = fopen(simple_profile.prof_file_name, "a+");
        if(!simple_profile.prof_file) {
            return -1;
        } 
        fprintf(simple_profile.prof_file, "I/O Timing results\n");
        fprintf(simple_profile.prof_file, "Operations    : \t %-18s \t %-18s \t %-18s \t %-18s\n",
            "min", "max", "mean", "var");
    }
    
	// IO timing info
    for(i = 0; i < simple_profile.io_count; i ++) {
        // gather timing info to proc 0
        double open_time_min, open_time_max, open_time_mean, open_time_var;
        double open_stime_min, open_stime_max, open_stime_mean, open_stime_var;
        double open_etime_min, open_etime_max, open_etime_mean, open_etime_var;
        double write_time_min, write_time_max, write_time_mean, write_time_var;
        double write_stime_min, write_stime_max, write_stime_mean, write_stime_var;
        double write_etime_min, write_etime_max, write_etime_mean, write_etime_var;
        double close_time_min, close_time_max, close_time_mean, close_time_var;
        double close_stime_min, close_stime_max, close_stime_mean, close_stime_var;
        double close_etime_min, close_etime_max, close_etime_mean, close_etime_var;
        double total_io_time, total_io_time_min, total_io_time_max, total_io_time_mean, total_io_time_var;

        get_time_dist(simple_profile.communicator
                      ,&simple_profile.timing[i].open_time
                      ,&open_time_min
                      ,&open_time_max
                      ,&open_time_mean
                      ,&open_time_var
                      );

        get_time_dist(simple_profile.communicator
                      ,&simple_profile.timing[i].open_start_time
                      ,&open_stime_min
                      ,&open_stime_max
                      ,&open_stime_mean
                      ,&open_stime_var
                      );

        get_time_dist(simple_profile.communicator
                      ,&simple_profile.timing[i].open_end_time
                      ,&open_etime_min
                      ,&open_etime_max
                      ,&open_etime_mean
                      ,&open_etime_var
                      );

        get_time_dist(simple_profile.communicator
                      ,&simple_profile.timing[i].write_time
                      ,&write_time_min
                      ,&write_time_max
                      ,&write_time_mean
                      ,&write_time_var
                      );

        get_time_dist(simple_profile.communicator
                      ,&simple_profile.timing[i].write_start_time
                      ,&write_stime_min
                      ,&write_stime_max
                      ,&write_stime_mean
                      ,&write_stime_var
                      );

        get_time_dist(simple_profile.communicator
                      ,&simple_profile.timing[i].write_end_time
                      ,&write_etime_min
                      ,&write_etime_max
                      ,&write_etime_mean
                      ,&write_etime_var
                      );

        get_time_dist(simple_profile.communicator
                      ,&simple_profile.timing[i].close_time
                      ,&close_time_min
                      ,&close_time_max
                      ,&close_time_mean
                      ,&close_time_var
                      );    

        get_time_dist(simple_profile.communicator
                      ,&simple_profile.timing[i].close_start_time
                      ,&close_stime_min
                      ,&close_stime_max
                      ,&close_stime_mean
                      ,&close_stime_var
                      );

        get_time_dist(simple_profile.communicator
                      ,&simple_profile.timing[i].close_end_time
                      ,&close_etime_min
                      ,&close_etime_max
                      ,&close_etime_mean
                      ,&close_etime_var
                      );

        total_io_time = simple_profile.timing[i].close_end_time - simple_profile.timing[i].open_start_time; 

        get_time_dist(simple_profile.communicator
                      ,&total_io_time
                      ,&total_io_time_min
                      ,&total_io_time_max
                      ,&total_io_time_mean
                      ,&total_io_time_var
                      ); 

        if(simple_profile.myid == 0) {
            fprintf(simple_profile.prof_file, "cycle no \t %d\n", simple_profile.timing[i].cycle);
            fprintf(simple_profile.prof_file, "io count \t %d\n", i);

            fprintf(simple_profile.prof_file, "# Open        : \t %-18f \t %-18f \t %-18f \t %-18f\n",
                    open_time_min, open_time_max, open_time_mean, open_time_var);
            fprintf(simple_profile.prof_file, "# Open start  : \t %-18f \t %-18f \t %-18f \t %-18f\n",
                    open_stime_min, open_stime_max, open_stime_mean, open_stime_var);
            fprintf(simple_profile.prof_file, "# Open end    : \t %-18f \t %-18f \t %-18f \t %-18f\n",
                    open_etime_min, open_etime_max, open_etime_mean, open_etime_var);

            fprintf(simple_profile.prof_file, "# Write       : \t %-18f \t %-18f \t %-18f \t %-18f\n",
                    write_time_min, write_time_max, write_time_mean, write_time_var);
            fprintf(simple_profile.prof_file, "# Write start : \t %-18f \t %-18f \t %-18f \t %-18f\n",
                    write_stime_min, write_stime_max, write_stime_mean, write_stime_var);
            fprintf(simple_profile.prof_file, "# Write end   : \t %-18f \t %-18f \t %-18f \t %-18f\n",
                    write_etime_min, write_etime_max, write_etime_mean, write_etime_var);

            fprintf(simple_profile.prof_file, "# Close       : \t %-18f \t %-18f \t %-18f \t %-18f\n",
                    close_time_min, close_time_max, close_time_mean, close_time_var);
            fprintf(simple_profile.prof_file, "# Close start : \t %-18f \t %-18f \t %-18f \t %-18f\n",
                    close_stime_min, close_stime_max, close_stime_mean, close_stime_var);
            fprintf(simple_profile.prof_file, "# Close end   : \t %-18f \t %-18f \t %-18f \t %-18f\n",
                    close_etime_min, close_etime_max, close_etime_mean, close_etime_var);

            fprintf(simple_profile.prof_file, "# Total       : \t %-18f \t %-18f \t %-18f \t %-18f\n",
                    total_io_time_min, total_io_time_max, total_io_time_mean, total_io_time_var);
        }
    }

    // simulation cycle timing info
    if(simple_profile.myid == 0) {        
		fprintf(simple_profile.prof_file, "\nTotal Simulation Timing results\n");
        fprintf(simple_profile.prof_file, "Cycle No.         : \t %-18s \t %-18s \t %-18s \t %-18s\n",
            "min", "max", "mean", "var");
    }
    
    for(i = 0; i < simple_profile.cycle_count; i ++) {
        // gather timing info to proc 0
        double cycle_time_min, cycle_time_max, cycle_time_mean, cycle_time_var;
        double cycle_stime_min, cycle_stime_max, cycle_stime_mean, cycle_stime_var;
        double cycle_etime_min, cycle_etime_max, cycle_etime_mean, cycle_etime_var;

        get_time_dist(simple_profile.communicator
                      ,&simple_profile.cycle_timing[i].cycle_time
                      ,&cycle_time_min
                      ,&cycle_time_max
                      ,&cycle_time_mean
                      ,&cycle_time_var
                      );

        get_time_dist(simple_profile.communicator
                      ,&simple_profile.cycle_timing[i].cycle_start_time
                      ,&cycle_stime_min
                      ,&cycle_stime_max
                      ,&cycle_stime_mean
                      ,&cycle_stime_var
                      );

        get_time_dist(simple_profile.communicator
                      ,&simple_profile.cycle_timing[i].cycle_end_time
                      ,&cycle_etime_min
                      ,&cycle_etime_max
                      ,&cycle_etime_mean
                      ,&cycle_etime_var
                      );

        if(simple_profile.myid == 0) {
            fprintf(simple_profile.prof_file, "cycle no \t %d\n", simple_profile.cycle_timing[i].cycle);

            fprintf(simple_profile.prof_file, "# Total Time  : \t %-18f \t %-18f \t %-18f \t %-18f\n",
                    cycle_time_min, cycle_time_max, cycle_time_mean, cycle_time_var);
            fprintf(simple_profile.prof_file, "# Start Time  : \t %-18f \t %-18f \t %-18f \t %-18f\n",
                    cycle_stime_min, cycle_stime_max, cycle_stime_mean, cycle_stime_var);
            fprintf(simple_profile.prof_file, "# End Time    : \t %-18f \t %-18f \t %-18f \t %-18f\n",
                    cycle_etime_min, cycle_etime_max, cycle_etime_mean, cycle_etime_var);
        }
    }

    if(simple_profile.myid == 0) {
        if(simple_profile.prof_file) {
            rc = fflush(simple_profile.prof_file);
            if(!rc) return rc;
            rc = fclose(simple_profile.prof_file); 
            if(!rc) return rc;
        }
    }
    return 0;
}

/*
 * get local timestamp
 */
double get_timestamp()
{
#ifdef TIMER_MPI_WTIME 
    return MPI_Wtime();
#else   
    struct timeval timestamp;
    gettimeofday(&timestamp, NULL);

    double realtime = timestamp.tv_sec + timestamp.tv_usec / 1.0e6;
    return realtime;
#endif
}

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
                  )
{
    int nr_procs;
    double sq_time_part, sum, sum_sq;
    
    MPI_Comm_size(comm, &nr_procs);

    MPI_Allreduce(time, max, 1, MPI_DOUBLE, MPI_MAX, comm);
    MPI_Allreduce(time, min, 1, MPI_DOUBLE, MPI_MIN, comm);
    MPI_Allreduce(time, &sum, 1, MPI_DOUBLE, MPI_SUM, comm);

    *mean = sum / nr_procs;

    /* calculate our part of the variance */
    sq_time_part = *time - *mean;
    sq_time_part = sq_time_part * sq_time_part;
    MPI_Allreduce(&sq_time_part, &sum_sq, 1, MPI_DOUBLE, MPI_SUM, comm);

    if (nr_procs > 1) {
        *var = sqrt( sum_sq / ((double) nr_procs - 1.0) ); // sample standard deviation
    }
    else {
        *var = 0.0;
    }

    return;
}

/*
 * initialize profiling 
 *
 * used inside ADIOS so we can record io timing info for each every adios group
 */
int init_prof_internal(char *prof_file_name, struct adios_group_list_struct *group_list)
{
    if(!prof_file_name || !strcmp(prof_file_name, "")) {
        profile.split_result = 1;
        fprintf(stderr, "ADIOS Profiling: log file not specified; write timing result"
                        " to seperate files for each group.\n");    
        profile.prof_file_name = NULL;          
    }
    else {
        profile.split_result = 0;
        profile.prof_file_name = prof_file_name;
    }

    profile.num_groups = 0;
    struct adios_group_list_struct *gp_ptr = group_list;  
    while(gp_ptr) {
        profile.num_groups ++;
        gp_ptr = gp_ptr->next;
    }

    profile.gp_prof_info = (group_prof_info *) malloc(profile.num_groups * sizeof(group_prof_info));

    gp_ptr = group_list;  
    int i = 0;
    while(i < profile.num_groups) {
        profile.gp_prof_info[i].adios_group_id = gp_ptr->group;
        profile.gp_prof_info[i].communicator = MPI_COMM_NULL;
        profile.gp_prof_info[i].comm_size = 0;
        profile.gp_prof_info[i].myid = 0;

        // TODO: hardcode max_io_count now; later change to use link list or set by env variable
        profile.gp_prof_info[i].max_io_count = MAX_IO_COUNT; 
        profile.gp_prof_info[i].io_count = 0;

        // NOTE: 
        // timing_info is initialized with all bits set to 0. 
        // Therefore, timing[i].open_start_time == 0 means this process didn't do IO at that cycle
        profile.gp_prof_info[i].timing = (timing_info *) calloc(profile.gp_prof_info[i].max_io_count, sizeof(timing_info));         

        gp_ptr = gp_ptr->next;
        i ++;
    }

    return 0;
}

/*
 * initialize profiling 
 *
 * Fortran interface
 */
int init_prof_all_(char *prof_file_name, int prof_file_name_size) {
    struct adios_group_list_struct * g = adios_get_groups();
    init_prof_internal(prof_file_name, g);
}

/*
 * find profiling struct by group id
 */
group_prof_info *get_prof_by_group_name(char *group_name)
{
    int i = 0;
    for(; i < profile.num_groups; i ++) {
//fprintf(stderr,"------------ADIOS Profiling: group: %s : %s\n", profile.gp_prof_info[i].adios_group_id->name, group_name);
        if(!strcmp(profile.gp_prof_info[i].adios_group_id->name, group_name)) {
            return &profile.gp_prof_info[i];
        }
    }

    return NULL; 
}

/*
 * record open start time for specified group
 *
 * Fortran interface
 */
void open_start_for_group_(long long *gp_prof_handle, char *group_name, int *cycle, int *gp_name_size)
{
int rank;
MPI_Comm_rank(adios_mpi_comm_world, &rank);
fprintf(stderr, "------------ADIOS Profiling: open start group: %s rank: %d\n", group_name, rank);

    group_prof_info *gp_prof = get_prof_by_group_name(group_name);

    if(!gp_prof) {
        fprintf(stderr, "ADIOS Profiling: unkown group \"%s\"\n", group_name);
    }

    *gp_prof_handle = (long long)gp_prof;

    if(gp_prof->io_count < gp_prof->max_io_count) {
        gp_prof->timing[gp_prof->io_count].cycle = *cycle;
        gp_prof->timing[gp_prof->io_count].open_start_time = get_timestamp();
    }
    else {
        fprintf(stderr, "ADIOS Profiling: io count exceeds max limit (%d)\n", gp_prof->max_io_count);
        exit(-1);
    }
fprintf(stderr, "------------ADIOS Profiling: open start end group: %s rank: %d\n", group_name, rank);
}

/*
 * record open end time for specified group
 *
 * Fortran interface
 */
void open_end_for_group_(long long *gp_prof_handle, int *cycle)
{
    group_prof_info *gp_prof = (group_prof_info *)*gp_prof_handle;
    if(!gp_prof) {
        fprintf(stderr, "ADIOS Profiling: unkown group\n");
    }

    gp_prof->timing[gp_prof->io_count].cycle = *cycle;
    gp_prof->timing[gp_prof->io_count].open_end_time = get_timestamp();
    gp_prof->timing[gp_prof->io_count].open_time = gp_prof->timing[gp_prof->io_count].open_end_time 
                                              - gp_prof->timing[gp_prof->io_count].open_start_time;
}

/*
 * record write start time for specified group
 *
 * Fortran interface
 */
void write_start_for_group_(long long *gp_prof_handle, int *cycle)
{
    group_prof_info *gp_prof = (group_prof_info *)*gp_prof_handle;
    if(!gp_prof) {
        fprintf(stderr, "ADIOS Profiling: unkown group\n");
    }

    gp_prof->timing[gp_prof->io_count].cycle = *cycle;
    gp_prof->timing[gp_prof->io_count].write_start_time = get_timestamp();
}

/*
 * record write end time for specified group
 *
 * Fortran interface
 */
void write_end_for_group_(long long *gp_prof_handle, int *cycle)
{
    group_prof_info *gp_prof = (group_prof_info *)*gp_prof_handle;
    if(!gp_prof) {
        fprintf(stderr, "ADIOS Profiling: unkown group\n");
    }

    gp_prof->timing[gp_prof->io_count].cycle = *cycle;
    gp_prof->timing[gp_prof->io_count].write_end_time = get_timestamp();
    gp_prof->timing[gp_prof->io_count].write_time = gp_prof->timing[gp_prof->io_count].write_end_time 
                                              - gp_prof->timing[gp_prof->io_count].write_start_time;
}

MPI_Comm get_comm_by_group(struct adios_group_struct *group);

/*
 * record close start time for specified group
 *
 * Fortran interface
 */
void close_start_for_group_(long long *gp_prof_handle, int *cycle)
{
int rank;
MPI_Comm_rank(adios_mpi_comm_world, &rank);
fprintf(stderr, "------------ADIOS Profiling: close start group: rank: %d\n", rank);


    group_prof_info *gp_prof = (group_prof_info *)*gp_prof_handle;
    if(!gp_prof) {
        fprintf(stderr, "ADIOS Profiling: unkown group\n");
    }

    // Right now is the only chance to get group's communicator
    // because at the end of adios_close all elements in adios_var struct
    // will get cleaned
    gp_prof->communicator = get_comm_by_group(gp_prof->adios_group_id);

    gp_prof->timing[gp_prof->io_count].cycle = *cycle;
    gp_prof->timing[gp_prof->io_count].close_start_time = get_timestamp();
}

/*
 * record close end time for specified group
 *
 * Fortran interface
 */
void close_end_for_group_(long long *gp_prof_handle, int *cycle)
{
int rank;
MPI_Comm_rank(adios_mpi_comm_world, &rank);
fprintf(stderr, "------------ADIOS Profiling: close end group: rank: %d\n", rank);


    group_prof_info *gp_prof = (group_prof_info *)*gp_prof_handle;
    if(!gp_prof) {
        fprintf(stderr, "ADIOS Profiling: unkown group\n");
    }

    gp_prof->timing[gp_prof->io_count].cycle = *cycle;
    gp_prof->timing[gp_prof->io_count].close_end_time = get_timestamp();
    gp_prof->timing[gp_prof->io_count].close_time = gp_prof->timing[gp_prof->io_count].close_end_time 
                                              - gp_prof->timing[gp_prof->io_count].close_start_time;
    gp_prof->io_count ++;
}

/*
 * get communicator handle for specified group
 */
MPI_Comm get_comm_by_group(struct adios_group_struct *group)
{
    if (group->group_comm) {
        struct adios_var_struct *var = adios_find_var_by_name(group->vars, group->group_comm, group->all_unique_var_names);

        if (var) {
            if (var->data) {
                return *(MPI_Comm *)var->data;
            }
            else {
                fprintf (stderr, "ADIOS Profiling: coordination-communicator: %s not provided. "
                                 "Using adios_mpi_comm_world instead\n"
                        ,group->group_comm
                        );
                return adios_mpi_comm_world;
            }
        }
        else {
            fprintf (stderr, "ADIOS Profiling: coordination-communicator: %s not found in "
                             "adios-group.  Using adios_mpi_comm_world instead\n"
                    ,group->group_comm
                    );

            return adios_mpi_comm_world;
        }
    }
    else
    {
        return MPI_COMM_NULL;
    }
}

/*
 * Report timing info for all groups
 *
 * Fortran interface  
 */
int finalize_prof_all_()
{
    int i;
    int g;
    int rc;
    int global_rank;
    int global_size;
    MPI_Comm_size(adios_mpi_comm_world, &global_size);
    MPI_Comm_rank(adios_mpi_comm_world, &global_rank);
    char print_buf[MAX_PRINT_BUF_SIZE];
    char *cur = print_buf;
    int num_written = 0;
    int print_flag = 0;

    for(g = 0; g < profile.num_groups; g ++) {
        group_prof_info *gp_prof = &profile.gp_prof_info[g];

        if(gp_prof->communicator == MPI_COMM_NULL) {
            // if group communicator is not specified in config.xml, then assume
            // this processes did IO independently. This complies with adios_mpi
            // and adios_mpi_cio convention 
            gp_prof->communicator = adios_mpi_comm_self;
        }

        MPI_Comm_size(gp_prof->communicator, &gp_prof->comm_size);
        MPI_Comm_rank(gp_prof->communicator, &gp_prof->myid);      

        memset(print_buf, 0, MAX_PRINT_BUF_SIZE);
        cur = print_buf;
        print_flag = 0;
 
        if(global_rank == 0) {
            print_flag = 1;
            num_written = sprintf(cur, "Timing results for group: %s\n", gp_prof->adios_group_id->name);
            cur += num_written;
            num_written = sprintf(cur, "Operations    : \t %-18s \t %-18s \t %-18s \t %-18s\n",
                "min", "max", "mean", "var");
            cur += num_written;
        }
fprintf(stderr, "--------11111 timing result for group: %s : rank %d\n", gp_prof->adios_group_id->name, global_rank);


        for(i = 0; i < gp_prof->io_count; i ++) {
            // gather timing info to proc 0
            double open_time_min, open_time_max, open_time_mean, open_time_var;
            double open_stime_min, open_stime_max, open_stime_mean, open_stime_var;
            double open_etime_min, open_etime_max, open_etime_mean, open_etime_var;
            double write_time_min, write_time_max, write_time_mean, write_time_var;
            double write_stime_min, write_stime_max, write_stime_mean, write_stime_var;
            double write_etime_min, write_etime_max, write_etime_mean, write_etime_var;
            double close_time_min, close_time_max, close_time_mean, close_time_var;
            double close_stime_min, close_stime_max, close_stime_mean, close_stime_var;
            double close_etime_min, close_etime_max, close_etime_mean, close_etime_var;
            double total_io_time, total_io_time_min, total_io_time_max, total_io_time_mean, total_io_time_var;

            // since there may be only part of processes within the comm which actually did IO,
            // we need to first check how many processes did IO
            int did_io = (gp_prof->timing[i].open_start_time != 0)? 1 : 0;
            MPI_Comm comm_io;
            rc = MPI_Comm_split(gp_prof->communicator, did_io, gp_prof->myid, &comm_io);
            if(rc != MPI_SUCCESS) {
                int error_class; 
                MPI_Error_class(rc, &error_class);
                fprintf(stderr, "ADIOS Profiling: Error %d returned by MPI_Comm_split() for group %s\n", 
                    error_class, gp_prof->adios_group_id->name);
                return -1;
            }

            int comm_io_rank, comm_io_size;
            MPI_Comm_rank(comm_io, &comm_io_rank);
            MPI_Comm_size(comm_io, &comm_io_size);
            
            if(did_io) {
/*                if(comm_io_size == 1) {
                    open_time_min = gp_prof->timing[i].open_time;
                    open_time_max = gp_prof->timing[i].open_time;
                    open_time_mean = gp_prof->timing[i].open_time;
                    open_time_var = 0;
                    open_stime_min = gp_prof->timing[i].open_start_time;
                    open_stime_max = gp_prof->timing[i].open_start_time;
                    open_stime_mean = gp_prof->timing[i].open_start_time;
                    open_stime_var = 0;
                    open_etime_min = gp_prof->timing[i].open_end_time;
                    open_etime_max = gp_prof->timing[i].open_end_time;
                    open_etime_mean = gp_prof->timing[i].open_end_time;
                    open_etime_var = 0;
                    write_time_min = gp_prof->timing[i].write_time;
                    write_time_max = gp_prof->timing[i].write_time;
                    write_time_mean = gp_prof->timing[i].write_time;
                    write_time_var = 0;
                    write_stime_min = gp_prof->timing[i].write_start_time;
                    write_stime_max = gp_prof->timing[i].write_start_time;
                    write_stime_mean = gp_prof->timing[i].write_start_time;
                    write_stime_var = 0;
                    write_etime_min = gp_prof->timing[i].write_end_time;
                    write_etime_max = gp_prof->timing[i].write_end_time;
                    write_etime_mean = gp_prof->timing[i].write_end_time;
                    write_etime_var = 0;
                    close_time_min = gp_prof->timing[i].close_time;
                    close_time_max = gp_prof->timing[i].close_time;
                    close_time_mean = gp_prof->timing[i].close_time;
                    close_time_var = 0;
                    close_stime_min = gp_prof->timing[i].close_start_time;
                    close_stime_max = gp_prof->timing[i].close_start_time;
                    close_stime_mean = gp_prof->timing[i].close_start_time;
                    close_stime_var = 0;
                    close_etime_min = gp_prof->timing[i].close_end_time;
                    close_etime_max = gp_prof->timing[i].close_end_time;
                    close_etime_mean = gp_prof->timing[i].close_end_time;
                    close_etime_var = 0;
                    total_io_time = gp_prof->timing[i].close_end_time - gp_prof->timing[i].open_start_time;
                    total_io_time_min = total_io_time;
                    total_io_time_max = total_io_time;
                    total_io_time_mean = total_io_time;
                    total_io_time_var = 0;
                }
                else {
*/  
              get_time_dist(comm_io
                              ,&gp_prof->timing[i].open_time
                              ,&open_time_min
                              ,&open_time_max
                              ,&open_time_mean
                              ,&open_time_var
                              );

                get_time_dist(comm_io
                              ,&gp_prof->timing[i].open_start_time
                              ,&open_stime_min
                              ,&open_stime_max
                              ,&open_stime_mean
                              ,&open_stime_var
                              );

                get_time_dist(comm_io
                              ,&gp_prof->timing[i].open_end_time
                              ,&open_etime_min
                              ,&open_etime_max
                              ,&open_etime_mean
                              ,&open_etime_var
                              );

                get_time_dist(comm_io
                              ,&gp_prof->timing[i].write_time
                              ,&write_time_min
                              ,&write_time_max
                              ,&write_time_mean
                              ,&write_time_var
                              );

                get_time_dist(comm_io
                              ,&gp_prof->timing[i].write_start_time
                              ,&write_stime_min
                              ,&write_stime_max
                              ,&write_stime_mean
                              ,&write_stime_var
                              );

                get_time_dist(comm_io
                              ,&gp_prof->timing[i].write_end_time
                              ,&write_etime_min
                              ,&write_etime_max
                              ,&write_etime_mean
                              ,&write_etime_var
                              );

                get_time_dist(comm_io
                              ,&gp_prof->timing[i].close_time
                              ,&close_time_min
                              ,&close_time_max
                              ,&close_time_mean
                              ,&close_time_var
                              );    

                get_time_dist(comm_io
                              ,&gp_prof->timing[i].close_start_time
                              ,&close_stime_min
                              ,&close_stime_max
                              ,&close_stime_mean
                              ,&close_stime_var
                              );

                get_time_dist(comm_io
                              ,&gp_prof->timing[i].close_end_time
                              ,&close_etime_min
                              ,&close_etime_max
                              ,&close_etime_mean
                              ,&close_etime_var
                              );

                total_io_time = gp_prof->timing[i].close_end_time - gp_prof->timing[i].open_start_time; 

                get_time_dist(comm_io
                              ,&total_io_time
                              ,&total_io_time_min
                              ,&total_io_time_max
                              ,&total_io_time_mean
                              ,&total_io_time_var
                              ); 

//                }

                if(comm_io_rank == 0) {
                    print_flag = 1;
                    num_written = sprintf(cur, "reported by %d\n", global_rank);
                    cur += num_written;
                    num_written = sprintf(cur, "cycle no \t %d\n", gp_prof->timing[i].cycle);
                    cur += num_written;
                    num_written = sprintf(cur, "io count \t %d\n", i);
                    cur += num_written;

                    num_written = sprintf(cur, "# Open        : \t %-18f \t %-18f \t %-18f \t %-18f\n",
                            open_time_min, open_time_max, open_time_mean, open_time_var);
                    cur += num_written;
                    num_written = sprintf(cur, "# Open start  : \t %-18f \t %-18f \t %-18f \t %-18f\n",
                            open_stime_min, open_stime_max, open_stime_mean, open_stime_var);
                    cur += num_written;
                    num_written = sprintf(cur, "# Open end    : \t %-18f \t %-18f \t %-18f \t %-18f\n",
                            open_etime_min, open_etime_max, open_etime_mean, open_etime_var);
                    cur += num_written;

                    num_written = sprintf(cur, "# Write       : \t %-18f \t %-18f \t %-18f \t %-18f\n",
                            write_time_min, write_time_max, write_time_mean, write_time_var);
                    cur += num_written;
                    num_written = sprintf(cur, "# Write start : \t %-18f \t %-18f \t %-18f \t %-18f\n",
                            write_stime_min, write_stime_max, write_stime_mean, write_stime_var);
                    cur += num_written;
                    num_written = sprintf(cur, "# Write end   : \t %-18f \t %-18f \t %-18f \t %-18f\n",
                            write_etime_min, write_etime_max, write_etime_mean, write_etime_var);
                    cur += num_written;

                    num_written = sprintf(cur, "# Close       : \t %-18f \t %-18f \t %-18f \t %-18f\n",
                            close_time_min, close_time_max, close_time_mean, close_time_var);
                    cur += num_written;
                    num_written = sprintf(cur, "# Close start : \t %-18f \t %-18f \t %-18f \t %-18f\n",
                            close_stime_min, close_stime_max, close_stime_mean, close_stime_var);
                    cur += num_written;
                    num_written = sprintf(cur, "# Close end   : \t %-18f \t %-18f \t %-18f \t %-18f\n",
                            close_etime_min, close_etime_max, close_etime_mean, close_etime_var);
                    cur += num_written;

                    num_written = sprintf(cur, "# Total       : \t %-18f \t %-18f \t %-18f \t %-18f\n",
                            total_io_time_min, total_io_time_max, total_io_time_mean, total_io_time_var);
                    cur += num_written;
                }
            }
 
            MPI_Comm_free(&comm_io);
        }

fprintf(stderr, "--------22222 timing result for group: %s : rank %d\n", gp_prof->adios_group_id->name, global_rank);
        unsigned long long my_length = cur - print_buf;
        unsigned long long my_offset = 0;
        MPI_Offset write_offset;
        MPI_Request request_hd;

        // calculate offset
        MPI_Scan(&my_length, &my_offset, 1, MPI_LONG_LONG, MPI_SUM, adios_mpi_comm_world);
        write_offset = my_offset - my_length;

        if(profile.split_result == 1) {
            profile.prof_file_name = (char *)calloc(strlen(gp_prof->adios_group_id->name)+1+5, sizeof(char));    
            if(!profile.prof_file_name) {
                fprintf(stderr, "ADIOS Profiling: Cannot allocate memory\n");
                return -1;                
            } 
            strcat(profile.prof_file_name, gp_prof->adios_group_id->name);
            strcat(profile.prof_file_name, ".prof");
        }

    	// open log file using MPI IO API
        int err_class;
        rc = MPI_File_open(adios_mpi_comm_world
                          ,profile.prof_file_name
                          ,MPI_MODE_WRONLY | MPI_MODE_APPEND | MPI_MODE_CREATE
                          ,NULL, &profile.prof_file_ptr
                          );

        if (rc != MPI_SUCCESS) {
            fprintf (stderr, "ADIOS Profiling: Error opening file %s\n", profile.prof_file_name);
        }   

        // set file pointer to the right position
        rc = MPI_File_seek(profile.prof_file_ptr, write_offset, MPI_SEEK_SET);

fprintf(stderr, "--------33333 timing result for group: %s : rank %d\n", gp_prof->adios_group_id->name, global_rank);
        // write timing info to log file using explicit offset
        // the timing info will be ordered by global_rank
        MPI_Status status;
        if(print_flag) {
            rc = MPI_File_write(profile.prof_file_ptr
				                  ,print_buf
				                  ,my_length
				                  ,MPI_CHAR
				                  ,&status
				                  );

            if(rc != MPI_SUCCESS) {
                int error_class; 
                MPI_Error_class(rc, &error_class);
                fprintf(stderr, "ADIOS Profiling: Error %d returned by MPI_File_write() for group %s\n", 
                    error_class, gp_prof->adios_group_id->name);
                return -1;
            }
        }

fprintf(stderr, "--------44444 timing result for group: %s : rank %d\n", gp_prof->adios_group_id->name, global_rank);
        // flush and close log file
        rc = MPI_File_close (&profile.prof_file_ptr);
        if(rc != MPI_SUCCESS) {
            int error_class; 
            MPI_Error_class(rc, &error_class);
            fprintf(stderr, "ADIOS Profiling: Error %d returned by MPI_File_close() for group %s\n", 
                error_class, gp_prof->adios_group_id->name);
            return -1;
        }
        free(profile.prof_file_name);
fprintf(stderr, "--------55555 timing result for group: %s : rank %d\n", gp_prof->adios_group_id->name, global_rank);
    }

    // cleanup
    for(i = 0; i < profile.num_groups; i ++) {
        free(profile.gp_prof_info[i].timing);    
    }
    free(profile.gp_prof_info);

    return 0;
}
