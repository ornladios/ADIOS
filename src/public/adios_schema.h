#ifndef ADIOS_SCHEMA_H
#define ADIOS_SCHEMA_H

//#include <stdbool.h>

enum ADIOS_MESH_TYPE
{
     ADIOS_MESH_UNIFORM      = 1
    ,ADIOS_MESH_STRUCTURED   = 2
    ,ADIOS_MESH_RECTILINEAR  = 3
    ,ADIOS_MESH_UNSTRUCTURED = 4
};

//mesh structure used by write method
// ADIOS Schema: modifying mesh struct
struct adios_mesh_struct
{
    // ADIOS Schema: adding mesh names
    // Groups can have multiple meshes
    char * name;
    enum ADIOS_FLAG time_varying;
    enum ADIOS_MESH_TYPE type;
/*    union
    {
        struct adios_mesh_uniform_struct * uniform;
        struct adios_mesh_rectilinear_struct * rectilinear;
        struct adios_mesh_structured_struct * structured;
        struct adios_mesh_unstructured_struct * unstructured;
    };*/
    struct adios_mesh_struct * next;
};


typedef struct
{
    int num_dimensions;
    uint64_t * dimensions;
    double * origins;
    double * spacings;
    double * maximums;
} MESH_UNIFORM;

typedef struct
{
    int use_single_var;        // 1 means coordinates-single-var,0 means coordinates-multi-var
    int num_dimensions;
    uint64_t * dimensions;
    char ** coordinates;       // name of the variable(s) containing the rectilinear spacing values
} MESH_RECTILINEAR;

typedef struct
{
    int use_single_var;        // 1 means points-single-var, 0 mean points-multi-var
    int num_dimensions;
    uint64_t * dimensions;
    int nspaces;
    char ** points;            // name of the variable(s) containing the point coordinates 
} MESH_STRUCTURED;

// ADIOS Schema: supported cell types
enum ADIOS_CELL_TYPE
{
     ADIOS_CELL_PT         = 1
    ,ADIOS_CELL_LINE       = 2
    ,ADIOS_CELL_TRI        = 3
    ,ADIOS_CELL_QUAD       = 4
    ,ADIOS_CELL_HEX        = 5
    ,ADIOS_CELL_PRI        = 6
    ,ADIOS_CELL_TET        = 7
    ,ADIOS_CELL_PYR        = 8
};

typedef struct
{
    int nspaces;
    uint64_t npoints;
    int nvar_points;           // how much vars for points-multi-var, 1 for points-single-var
    char ** points;
    int ncsets;
    uint64_t * ccounts;
    char ** cdata;
    enum ADIOS_CELL_TYPE * ctypes;
} MESH_UNSTRUCTURED;


typedef struct {   //type returned by adios_inq_mesh for read method
    int id;
    char * name;
    int time_varying;           //0 means not time-varying, 1 means time-varying
    enum ADIOS_MESH_TYPE type;
    union
    {
        MESH_UNIFORM * uniform;
        MESH_RECTILINEAR * rectilinear;
        MESH_STRUCTURED * structured;
        MESH_UNSTRUCTURED * unstructured;
    } ;
} ADIOS_MESH;


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

/*
//////////////////////////////////////////////////////////////
// Main mesh structs
//////////////////////////////////////////////////////////////
struct adios_mesh_uniform_struct
{
    struct adios_mesh_item_list_struct * dimensions;
    struct adios_mesh_item_list_struct * origin;
    struct adios_mesh_item_list_struct * spacing;
    // ADIOS Schema: adding option to provide origin and maximum
    // instead restricting users to origin and spacing
    struct adios_mesh_item_list_struct * maximum;
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
    // ADIOS Schema: adding single/multi points option
    // adding nspace to allow 2D mesh in 3D for example,
    // finally adding the concept of cellset/cellsetcount
    enum ADIOS_FLAG points_single_var;
    struct adios_mesh_item_struct * nspace;
    struct adios_mesh_var_list_struct * points;
    struct adios_mesh_item_struct * points_count;
    struct adios_mesh_cell_list_list_struct * cell_list;
    struct adios_mesh_item_struct * cell_set_count;
};
*/
#endif
