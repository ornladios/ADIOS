#include <stdint.h>
#include "binpack-general.h"

#define STR_LEN 1000

enum ADIOS_METHOD_MODE {adios_mode_write  = 1
                       ,adios_mode_read   = 2
                       ,adios_mode_update = 3 // not supported yet
                       ,adios_mode_append = 4
                       };

struct adios_dimension_struct;
struct adios_var_struct;

struct adios_global_bounds_struct
{
    struct adios_dimension_struct * dimensions;
    struct adios_dimension_struct * offsets;
};

struct adios_var_struct
{
    uint32_t id;
    char * name;
    char * path;
    enum ADIOS_DATATYPES type;
    struct adios_dimension_struct * dimensions;
    struct adios_global_bounds_struct * global_bounds;
    enum ADIOS_FLAG copy_on_write;
    enum ADIOS_FLAG got_buffer;

    enum ADIOS_FLAG free_data;    // primarily used for writing
    void * data;                  // primarily used for reading
    unsigned long long data_size; // primarily used for reading

    struct adios_var_struct * next;
};

struct adios_attribute_struct
{
    char * name;
    char * path;

    // if var.data == 0, then it is a var.  Otherwise, that is the value.
    struct adios_var_struct var;

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
    uint32_t id;
    char * name;
    int var_count;
    enum ADIOS_FLAG adios_host_language_fortran;
    struct adios_var_struct * vars;
    struct adios_attribute_struct * attributes;
    const char * group_by;
    const char * group_comm;
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
    long long base_offset;
    long long offset;
    struct adios_group_struct * group;
    enum ADIOS_METHOD_MODE mode;
};

struct adios_dimension_item_struct
{
    uint64_t rank;
    struct adios_var_struct * var;
};

struct adios_dimension_struct
{
    struct adios_dimension_item_struct dimension;
    struct adios_dimension_struct * next;
};

typedef void (* ADIOS_INIT_FN) (const char * parameters
                               ,struct adios_method_struct * method
                               );
typedef void (* ADIOS_OPEN_FN) (struct adios_file_struct * fd
                               ,struct adios_method_struct * method
                               );
typedef void (* ADIOS_WRITE_FN) (struct adios_file_struct * fd
                                ,struct adios_var_struct * v
                                ,void * data
                                ,struct adios_method_struct * method
                                );
typedef void (* ADIOS_GET_WRITE_BUFFER_FN) (struct adios_file_struct * fd
                                           ,struct adios_var_struct * v
                                           ,unsigned long long * size
                                           ,void ** buffer
                                           ,struct adios_method_struct * method
                                           );
typedef void (* ADIOS_READ_FN) (struct adios_file_struct * fd
                               ,struct adios_var_struct * v
                               ,void * buffer
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
    long long buffer_len;
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
int adios_set_buffer_size (void);

unsigned long long adios_method_buffer_alloc (unsigned long long size);
int adios_method_buffer_free (unsigned long long size);

unsigned long long adios_size_of_var (struct adios_var_struct * v, void * data);
unsigned long long adios_size_of_attribute (struct adios_attribute_struct * a);

unsigned long long adios_data_size (struct adios_group_struct * g);

int adios_parse_config (const char * config);
struct adios_method_list_struct * adios_get_methods (void);
struct adios_group_list_struct * adios_get_groups (void);

struct adios_var_struct * adios_find_var_by_name (struct adios_var_struct * root
                                                 ,const char * name
                                                 );

struct adios_var_struct * adios_find_attribute_var_by_name
                                       (struct adios_attribute_struct * root
                                       ,const char * name
                                       );

void adios_pre_element_fetch (struct adios_bp_element_struct * element
                             ,void ** buffer, long long * buffer_size
                             ,void * private_data
                             );

void adios_post_element_fetch (struct adios_bp_element_struct * element
                              ,void * buffer, long long buffer_size
                              ,void * private_data
                              );

void adios_parse_buffer (struct adios_file_struct * fd, char * buffer
                        ,long long len
                        );

void adios_parse_dimension (char * dimension, struct adios_group_struct * g
                           ,struct adios_dimension_struct * dim
                           );

void adios_extract_string (char * out, const char * in, int size);

int adios_do_write_var (struct adios_var_struct * v
                       ,void * buf
                       ,unsigned long long buf_size
                       ,unsigned long long buf_start
                       ,unsigned long long * buf_end
                       );
int adios_do_write_attribute (struct adios_attribute_struct * a
                             ,void * buf
                             ,unsigned long long buf_size
                             ,unsigned long long buf_start
                             ,unsigned long long * buf_end
                             );

int adios_common_define_attribute (long long group, const char * name
                                  ,const char * path, const char * value
                                  ,const char * type, const char * var
                                  );

void adios_append_method (struct adios_method_struct * method);

void adios_add_method_to_group (struct adios_method_list_struct ** root
                               ,struct adios_method_struct * method
                               );

void adios_append_global_bounds (struct adios_global_bounds_struct * bounds);

void adios_append_group (struct adios_group_struct * group);

void adios_append_var (struct adios_var_struct ** root
                      ,struct adios_var_struct * var
                      );

void adios_append_dimension (struct adios_dimension_struct ** root
                            ,struct adios_dimension_struct * dimension
                            );

void adios_append_attribute (struct adios_attribute_struct ** root
                            ,struct adios_attribute_struct * attribute
                            );

void * adios_dupe_data (struct adios_var_struct * v, void * data);

void adios_dims_to_bp_dims (char * name
                           ,struct adios_dimension_struct * adios_dims
                           ,struct adios_global_bounds_struct * global_bounds
                           ,int * rank
                           ,struct adios_bp_dimension_struct * bp_dims
                           );

int adios_common_declare_group (long long * id, const char * name
                               ,const char * coordination_comm
                               ,const char * coordination_var
                               );

int adios_common_define_global_bounds (long long group_id
                                      ,const char * dimensions
                                      ,const char * offsets
                                      ,struct adios_global_bounds_struct ** b
                                      );

int adios_common_define_var (long long group_id, const char * name
                            ,const char * path, int type
                            ,int copy_on_write
                            ,const char * dimensions
                            ,struct adios_global_bounds_struct * global_bounds
                            );

int adios_common_select_method (int priority, const char * method
                               ,const char * parameters, const char * group 
                               ,const char * base_path, int iters
                               );

void adios_common_get_group (long long * group_id, const char * name);

// the following are defined in adios_transport_hooks.c
void adios_init_transports (struct adios_transport_struct ** transports);
int adios_parse_method (const char * buf, enum ADIOS_IO_METHOD * method);
