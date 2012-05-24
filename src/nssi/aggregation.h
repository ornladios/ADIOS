/*
 * aggregation.h
 *
 *  Created on: Mar 9, 2009
 *      Author: thkorde
 */

#ifndef AGGREGATION_H_
#define AGGREGATION_H_

//#include "adios.h"
#include "adios_types.h"
//#include "adios_internals.h"

#include "adios_nssi_args.h"

struct aggregation_chunk_details_t {
    int  fd;
    char var_path[ADIOS_PATH_MAX];
    char var_name[ADIOS_PATH_MAX];
    int  ndims;

    void     *buf;          /* the data */
    uint64_t  len;          /* length of buf in bytes */
    int       num_elements; /* number of datatype elements in buf (len/atype_size) */

    enum ADIOS_DATATYPES atype; /* adios type of data in buf*/
    int                  atype_size;

    char    **offset_path;
    char    **offset_name;
    uint64_t *offset;     /* starting corner (eg. 0,0,0 is the origin of a cube) */
    char    **count_path;
    char    **count_name;
    uint64_t *count;      /* num elements in each dimension (eg. 3,3,3 is a cube of size 3) */
};
typedef struct aggregation_chunk_details_t aggregation_chunk_details_t;

int use_aggregation(const int fd);
int use_caching(const int fd);
int use_collective(const int fd);
int use_independent(const int fd);
int use_direct(const int fd);
void add_file(const int fd, const write_type write_type);
void add_chunk(aggregation_chunk_details_t *chunk);
void cleanup_aggregation_chunks(const int fd);
void cleanup_aggregation_chunks(const int fd, const char *var_name);
int try_aggregation(const int fd);
int try_aggregation(const int fd, const char *var_name);
int aggregate_data_ready_to_write(const int fd, const char *var_name);
int cache_data_ready_to_write(const int fd, const char *var_name);
aggregation_chunk_details_t **get_chunks(const int fd, int *chunk_count);
aggregation_chunk_details_t **get_chunks(const int fd, const char *var_name, int *chunk_count);
void print_chunk(aggregation_chunk_details_t *c);



#endif /* AGGREGATION_H_ */
