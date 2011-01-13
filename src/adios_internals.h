/* 
 * ADIOS is freely available under the terms of the BSD license described
 * in the COPYING file in the top level directory of this source distribution.
 *
 * Copyright (c) 2008 - 2009.  UT-BATTELLE, LLC. All rights reserved.
 */

#ifndef ADIOS_INTERNALS_H
#define ADIOS_INTERNALS_H

#include <stdint.h>
#include <stdlib.h>

// need the enum for the transports
#include "adios_transport_hooks.h"
#include "adios_bp_v1.h"

enum ADIOS_METHOD_MODE {adios_mode_write  = 1
                       ,adios_mode_read   = 2
                       ,adios_mode_update = 3 // not supported yet
                       ,adios_mode_append = 4
                       };

struct adios_dimension_struct;
struct adios_var_struct;

// NCSU - Generic data for statistics
struct adios_stat_struct 
{
	void * data;
};


struct adios_var_struct
{
    uint16_t id;
    uint16_t parent_id;

    char * name;
    char * path;
    enum ADIOS_DATATYPES type;
    struct adios_dimension_struct * dimensions;
    enum ADIOS_FLAG got_buffer;
    enum ADIOS_FLAG is_dim;   // if it is a dimension, we need to
                              // track for netCDF usage
    uint64_t write_offset;  // offset this var was written at  [for writes]

    enum ADIOS_FLAG free_data;    // primarily used for writing
    void * data;                  // primarily used for reading
    uint64_t data_size;           // primarily used for reading

	// NCSU - Adding stat related variables
	struct adios_stat_struct ** stats; // 2D array. Complex numbers can contain upto 3 parts
	uint32_t bitmap;

    struct adios_var_struct * next;
};

// NCSU - structure for histogram
struct adios_hist_struct
{   
    double min; //minimum value of histogram ** for when we use complex variables
    double max; //maximum value of histogram
    uint32_t num_breaks; //number of break points for the histogram
    uint32_t * frequencies; //array of frequencies for the histogram
    double * breaks; //breaks array for the histogram, output this to gnuplot
};

struct adios_attribute_struct
{
    uint16_t id;
    char * name;
    char * path;
    enum ADIOS_DATATYPES type;
    void * value;
    struct adios_var_struct * var;
    uint64_t write_offset;  // offset this var was written at  [for writes]

    struct adios_attribute_struct * next;
};

struct adios_method_struct
{
    enum ADIOS_IO_METHOD m;
    char * base_path;
    char * method;
    void * method_data;
    char * parameters;
    int iterations;
    int priority;
    struct adios_group_struct * group;
};

struct adios_method_list_struct
{
    struct adios_method_struct * method;
    struct adios_method_list_struct * next;
};

enum ADIOS_MESH_TYPE
{
     ADIOS_MESH_UNIFORM      = 1
    ,ADIOS_MESH_STRUCTURED   = 2
    ,ADIOS_MESH_RECTILINEAR  = 3
    ,ADIOS_MESH_UNSTRUCTURED = 4
};

struct adios_mesh_uniform_struct;
struct adios_mesh_rectilinear_struct;
struct adios_mesh_structured_struct;
struct adios_mesh_unstructured_struct;

struct adios_mesh_struct
{
    enum ADIOS_FLAG time_varying;
    enum ADIOS_MESH_TYPE type;
    union 
    {
        struct adios_mesh_uniform_struct * uniform;
        struct adios_mesh_rectilinear_struct * rectilinear;
        struct adios_mesh_structured_struct * structured;
        struct adios_mesh_unstructured_struct * unstructured;
    };
};

struct adios_group_struct
{
    uint16_t id;
    uint16_t member_count;
    uint64_t group_offset;

    char * name;
    int var_count;
    enum ADIOS_FLAG adios_host_language_fortran;
    enum ADIOS_FLAG all_unique_var_names;
    struct adios_var_struct * vars;
    struct adios_var_struct * vars_written;
    struct adios_attribute_struct * attributes;
    char * group_comm;
    char * group_by;
    char * time_index_name;
    uint32_t time_index;
    enum ADIOS_FLAG stats_on;
    uint32_t process_id;

    struct adios_method_list_struct * methods;

    struct adios_mesh_struct * mesh;
};

struct adios_group_list_struct
{
    struct adios_group_struct * group;

    struct adios_group_list_struct * next;
};

struct adios_file_struct
{
    char * name;
    uint32_t subfile_index; // needs to be set in ADIOS method if generates subfiles
    struct adios_group_struct * group;
    enum ADIOS_METHOD_MODE mode;
    uint64_t data_size;
    uint64_t write_size_bytes;

    enum ADIOS_FLAG shared_buffer;

    uint64_t pg_start_in_file; //  where this pg started in the file

    uint64_t base_offset;   // where writing last ocurred

    char * buffer;          // buffer we use for building the output
    uint64_t offset;        // current offset to write at
    uint64_t bytes_written; // largest offset into buffer written to
    uint64_t buffer_size;   // how big the buffer is currently

    uint64_t vars_start;    // offset for where to put the var/attr count
    uint16_t vars_written;  // count of vars/attrs to write
};

struct adios_dimension_item_struct
{
    uint64_t rank;
    uint16_t id;
    enum ADIOS_FLAG time_index;
};

struct adios_dimension_struct
{
    struct adios_dimension_item_struct dimension;
    struct adios_dimension_item_struct global_dimension;
    struct adios_dimension_item_struct local_offset;
    struct adios_dimension_struct * next;
};

typedef void (* ADIOS_INIT_FN) (const char * parameters
                               ,struct adios_method_struct * method
                               );
typedef int (* ADIOS_OPEN_FN) (struct adios_file_struct * fd
                              ,struct adios_method_struct * method, void * comm
                              );
typedef enum ADIOS_FLAG (* ADIOS_SHOULD_BUFFER_FN)
                                       (struct adios_file_struct * fd
                                       ,struct adios_method_struct * method
                                       );
typedef void (* ADIOS_WRITE_FN) (struct adios_file_struct * fd
                                ,struct adios_var_struct * v
                                ,void * data
                                ,struct adios_method_struct * method
                                );
typedef void (* ADIOS_GET_WRITE_BUFFER_FN) (struct adios_file_struct * fd
                                           ,struct adios_var_struct * v
                                           ,uint64_t * size
                                           ,void ** buffer
                                           ,struct adios_method_struct * method
                                           );
typedef void (* ADIOS_READ_FN) (struct adios_file_struct * fd
                               ,struct adios_var_struct * v
                               ,void * buffer
                               ,uint64_t buffer_size
                               ,struct adios_method_struct * method
                               );
typedef void (* ADIOS_CLOSE_FN) (struct adios_file_struct * fd
                                ,struct adios_method_struct * method
                                );
typedef void (* ADIOS_FINALIZE_FN) (int mype
                                   ,struct adios_method_struct * method
                                   );
typedef void (* ADIOS_END_ITERATION_FN) (struct adios_method_struct * method);
typedef void (* ADIOS_START_CALCULATION_FN)
                                        (struct adios_method_struct * method);
typedef void (* ADIOS_STOP_CALCULATION_FN)
                                        (struct adios_method_struct * method);

struct adios_transport_struct
{
    ADIOS_INIT_FN adios_init_fn;
    ADIOS_OPEN_FN adios_open_fn;
    ADIOS_SHOULD_BUFFER_FN adios_should_buffer_fn;
    ADIOS_WRITE_FN adios_write_fn;
    ADIOS_GET_WRITE_BUFFER_FN adios_get_write_buffer_fn;
    ADIOS_READ_FN adios_read_fn;
    ADIOS_CLOSE_FN adios_close_fn;
    ADIOS_FINALIZE_FN adios_finalize_fn;
    ADIOS_END_ITERATION_FN adios_end_iteration_fn;
    ADIOS_START_CALCULATION_FN adios_start_calculation_fn;
    ADIOS_STOP_CALCULATION_FN adios_stop_calculation_fn;
};

struct adios_buffer_part_entry
{
    void * buffer;
    size_t buffer_size;
};

struct adios_parse_buffer_struct
{
    struct adios_var_struct * vars;
    enum ADIOS_FLAG all_unique_var_names;
    uint64_t buffer_len;
    void * buffer;
};

////////////////////////
// mesh support data structures
////////////////////////
struct adios_mesh_item_struct
{
    double rank;
    struct adios_var_struct * var;
};

struct adios_mesh_item_list_struct
{
    struct adios_mesh_item_struct item;
    struct adios_mesh_item_list_struct * next;
};

struct adios_mesh_var_list_struct
{
    struct adios_var_struct * var;
    struct adios_mesh_var_list_struct * next;
};

struct adios_mesh_cell_list_struct
{
    enum ADIOS_FLAG cells_uniform;
    struct adios_mesh_item_struct count;
    struct adios_var_struct * data;
    struct adios_mesh_item_struct type;
};

struct adios_mesh_cell_list_list_struct
{
    struct adios_mesh_cell_list_struct cell_list;
    struct adios_mesh_cell_list_list_struct * next;
};

//////////////////////////////////////////////////////////////
// Main mesh structs
//////////////////////////////////////////////////////////////
struct adios_mesh_uniform_struct
{
    struct adios_mesh_item_list_struct * dimensions;
    struct adios_mesh_item_list_struct * origin;
    struct adios_mesh_item_list_struct * spacing;
};

struct adios_mesh_rectilinear_struct
{
    enum ADIOS_FLAG coordinates_single_var;
    struct adios_mesh_item_list_struct * dimensions;
    struct adios_mesh_var_list_struct * coordinates;
};

struct adios_mesh_structured_struct
{
    enum ADIOS_FLAG points_single_var;
    struct adios_mesh_item_struct * nspace;
    struct adios_mesh_item_list_struct * dimensions;
    struct adios_mesh_var_list_struct * points;
};

struct adios_mesh_unstructured_struct
{
    struct adios_mesh_item_struct * components;
    struct adios_mesh_item_struct * points_count;
    struct adios_var_struct * points;
    struct adios_mesh_cell_list_list_struct * cell_list;
};

//////////////////////////////////////////////////////////////////////////////
// Function Delcarations
//////////////////////////////////////////////////////////////////////////////
uint64_t adios_size_of_var (struct adios_var_struct * v, void * data);
uint64_t adios_size_of_attribute (struct adios_attribute_struct * a);

uint64_t adios_data_size (struct adios_group_struct * g);

struct adios_method_list_struct * adios_get_methods (void);
struct adios_group_list_struct * adios_get_groups (void);

struct adios_var_struct * adios_find_var_by_name (struct adios_var_struct * root
                                                 ,const char * name
                                                 ,enum ADIOS_FLAG unique_names
                                                 );
struct adios_var_struct * adios_find_var_by_id (struct adios_var_struct * root
                                               ,uint16_t id
                                               );
struct adios_attribute_struct * adios_find_attribute_by_id
                                         (struct adios_attribute_struct * root
                                         ,uint16_t id
                                         );

struct adios_attribute_struct * adios_find_attribute_var_by_name
                                       (struct adios_attribute_struct * root
                                       ,const char * name
                                       );

struct adios_attribute_struct * adios_find_attribute_var_by_id
                                       (struct adios_attribute_struct * root
                                       ,uint16_t id
                                       );

/*void adios_parse_buffer (struct adios_file_struct * fd, char * buffer
                        ,uint64_t len
                        );*/

int adios_parse_dimension (const char * dimension
                          ,const char * global_dimension
                          ,const char * local_offset
                          ,struct adios_group_struct * g
                          ,struct adios_dimension_struct * dim
                          );

/*void adios_extract_string (char ** out, const char * in, int size);*/

int adios_common_define_attribute (int64_t group, const char * name
                                  ,const char * path
                                  ,enum ADIOS_DATATYPES type
                                  ,const char * value
                                  ,const char * var
                                  );

void adios_append_method (struct adios_method_struct * method);

void adios_add_method_to_group (struct adios_method_list_struct ** root
                               ,struct adios_method_struct * method
                               );

void adios_append_group (struct adios_group_struct * group);

// is the var name unique
enum ADIOS_FLAG adios_append_var (struct adios_var_struct ** root
                                 ,struct adios_var_struct * var
                                 ,uint16_t id
                                 );

void adios_append_dimension (struct adios_dimension_struct ** root
                            ,struct adios_dimension_struct * dimension
                            );

void adios_append_attribute (struct adios_attribute_struct ** root
                            ,struct adios_attribute_struct * attribute
                            ,uint16_t id
                            );

int adios_common_declare_group (int64_t * id, const char * name
                               ,enum ADIOS_FLAG host_language_fortran
                               ,const char * coordination_comm
                               ,const char * coordination_var
                               ,const char * time_index
                               ,enum ADIOS_FLAG stats
                               );

int adios_common_define_var (int64_t group_id, const char * name
                            ,const char * path, enum ADIOS_DATATYPES type
                            ,const char * dimensions
                            ,const char * global_dimensions
                            ,const char * local_offsets
                            );

int adios_common_define_var_characteristcs  (struct adios_group_struct * g, const char * var_name
                            				,const char * bin_interval 
                            				,const char * bin_min
                            				,const char * bin_max
                            				,const char * bin_count
                            				);

void adios_common_get_group (int64_t * group_id, const char * name);
int adios_common_free_group (int64_t id);

// ADIOS file format functions

uint16_t adios_calc_var_overhead_v1 (struct adios_var_struct * v);
uint32_t adios_calc_attribute_overhead_v1 (struct adios_attribute_struct * a);
uint64_t adios_calc_overhead_v1 (struct adios_file_struct * fd);

int adios_write_version_v1 (char ** buffer
                           ,uint64_t * buffer_size
                           ,uint64_t * buffer_offset
                           );
int adios_write_process_group_header_v1 (struct adios_file_struct * fd
                                        ,uint64_t total_size
                                        );

// data is only there for sizing
uint64_t adios_write_var_header_v1 (struct adios_file_struct * fd
                                   ,struct adios_var_struct * v
                                   );
int adios_generate_var_characteristics_v1 (struct adios_file_struct * fd
                                          ,struct adios_var_struct * var
                                          );
uint16_t adios_write_var_characteristics_v1 (struct adios_file_struct * fd
                                            ,struct adios_var_struct * var
                                            );
int adios_write_var_payload_v1 (struct adios_file_struct * fd
                               ,struct adios_var_struct * var
                               );
int adios_write_attribute_v1 (struct adios_file_struct * fd
                             ,struct adios_attribute_struct * a
                             );
int adios_write_open_vars_v1 (struct adios_file_struct * fd);
int adios_write_close_vars_v1 (struct adios_file_struct * fd);
int adios_write_open_attributes_v1 (struct adios_file_struct * fd);
int adios_write_close_attributes_v1 (struct adios_file_struct * fd);
int adios_write_index_v1 (char ** buffer
                         ,uint64_t * buffer_size
                         ,uint64_t * buffer_offset
                         ,uint64_t index_start
                         ,struct adios_index_process_group_struct_v1 * pg_root
                         ,struct adios_index_var_struct_v1 * vars_root
                         ,struct adios_index_attribute_struct_v1 * attrs_root
                         );

void adios_build_index_v1 (struct adios_file_struct * fd
                       ,struct adios_index_process_group_struct_v1 ** pg_root
                       ,struct adios_index_var_struct_v1 ** vars_root
                       ,struct adios_index_attribute_struct_v1 ** attrs_root
                       );
void adios_merge_index_v1 (
                   struct adios_index_process_group_struct_v1 ** main_pg_root
                  ,struct adios_index_var_struct_v1 ** main_vars_root
                  ,struct adios_index_attribute_struct_v1 ** main_attrs_root
                  ,struct adios_index_process_group_struct_v1 * new_pg_root
                  ,struct adios_index_var_struct_v1 * new_vars_root
                  ,struct adios_index_attribute_struct_v1 * new_attrs_root
                  );
void adios_clear_index_v1 (struct adios_index_process_group_struct_v1 * pg_root
                          ,struct adios_index_var_struct_v1 * vars_root
                          ,struct adios_index_attribute_struct_v1 * attrs_root
                          );

uint64_t adios_get_type_size (enum ADIOS_DATATYPES type, void * var);
uint64_t adios_get_var_size (struct adios_var_struct * var
                            ,struct adios_group_struct * group, void * data
                            );

uint64_t adios_get_stat_size (void * data, enum ADIOS_DATATYPES type, enum ADIOS_STAT stat_id);
uint8_t adios_get_stat_set_count (enum ADIOS_DATATYPES type);

const char * adios_type_to_string_int (int type);
const char * adios_file_mode_to_string (int mode);

// the following are defined in adios_transport_hooks.c
void adios_init_transports (struct adios_transport_struct ** transports);
int adios_parse_method (const char * buf, enum ADIOS_IO_METHOD * method
                       ,int * requires_group_comm
                       );

/* some internal functions that adios_internals.c and adios_internals_mxml.c share */
int adios_int_is_var (const char * temp); // 1 == yes, 0 == no
int adios_int_is_num (char * temp); // 1 == yes, 0 == no

// queue code for adaptive message passing
#ifdef __cplusplus
extern "C" {
#endif

#if 1
/* LIST */
typedef struct ListElmt_
{
    void             * data;
    struct ListElmt_ * next;
} ListElmt;

typedef struct List_
{
    volatile int      size;
    int      (* match) (const void * key1, const void * key2);
    void     (* destroy) (void * data);
    ListElmt * head;
    ListElmt * tail;
} List;

void list_init (List * list, void (* destroy) (void * data));
void list_destroy (List * list);
int list_ins_next (List * list, ListElmt * element, const void * data);
int list_rem_next (List * list, ListElmt * element, void ** data);
#define list_size(list) ((list)->size)
#define list_head(list) ((list)->head)
#define list_tail(list) ((list)->tail)
#define list_is_head(list, element) ((element) == (list)->head ? 1 : 0)
#define list_is_tail(element) ((element)->next == NULL ? 1 : 0)
#define list_data(element) ((element)->data)
#define list_next(element) ((element)->next)

/* QUEUE */
typedef List Queue;

#define queue_init list_init
#define queue_destroy list_destroy
int queue_enqueue (Queue * queue, const void * data);
int queue_dequeue (Queue * queue, void ** data);
#define queue_peek(queue) ((queue)->head == NULL ? NULL : (queue)->head->data)
#define queue_size list_size

#endif /* 0 */

#ifdef __cplusplus
}
#endif

#endif
