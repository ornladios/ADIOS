#ifndef _FLEXPATH_H
#define _FLEXPATH_H


#include "core/adios_logger.h"

#define READER_CONTACT_FILE "reader_info.txt"
#define WRITER_CONTACT_FILE "writer_info.txt"
#define READER_READY_FILE "reader_ready.txt"
#define WRITER_READY_FILE "writer_ready.txt"
#define FP_RANK_ATTR_NAME "fp_rank_num"
#define FP_DST_ATTR_NAME "fp_dst_rank"
#define FP_DIM_ATTR_NAME "fp_dim"
#define FP_NDIMS_ATTR_NAME "fp_ndims"

#define CLOSE_MSG 0
#define OPEN_MSG 1
#define ACK_MSG 2
#define INIT_MSG 3
#define EOS_MSG 4

#define perr(...) if(getenv("FP_DEBUG")) fprintf(stderr, __VA_ARGS__);

#define fp_log(LOG, ...)                             \
            if(getenv("FP_DEBUG")) {    \
                if(strcmp(getenv("FP_DEBUG"),"ALL")==0) {          \
                    fprintf(stderr, __VA_ARGS__);   \
                } else if(strcmp(getenv("FP_DEBUG"),LOG)==0) {     \
                    fprintf(stderr, __VA_ARGS__);   \
                }                                   \
            }

#define fp_write_log(LOG, ...)                                      \
            if(getenv("FP_DEBUG")) {                                \
                if(strcmp(getenv("FP_DEBUG"),"ALL")==0) {           \
                    fprintf(stderr, "%d %s:", flexpathWriteData.rank, LOG);   \
                    fprintf(stderr, __VA_ARGS__);                   \
                } else {                                            \
                    char* env_tok;                                  \
                    char* env = strdup(getenv("FP_DEBUG"));         \
                    env_tok = strtok(env, ",");                     \
                    while(env_tok) {                                \
                        if(strcmp(env_tok, LOG)==0) {               \
                    fprintf(stderr, "%d %s:", flexpathWriteData.rank, LOG);   \
                    fprintf(stderr, __VA_ARGS__);                   \
                        }                                           \
                        env_tok = strtok(NULL, ",");                \
                    }                                               \
                }                                                   \
            }
            

//adios_logger(4,1, __VA_ARGS__);

#define CONTACT_STR_LEN 50

typedef enum { FORMAT=0, DATA, EVGROUP } Flush_type;

/********** offset_struct **********/
/*
 * Contains the offset information for a variable for all writers.
 * offsets_per_rank is == ndims.
 */
typedef struct _offset_struct{
    int offsets_per_rank;
    int total_offsets;
    int * local_dimensions;
    int * local_offsets;
} offset_struct;

static FMField offset_struct_field_list[]=
{
    {"offsets_per_rank", "integer", sizeof(int), FMOffset(offset_struct*, offsets_per_rank)},
    {"total_offsets", "integer", sizeof(int), FMOffset(offset_struct*, total_offsets)},
    {"local_dimensions", "integer[total_offsets]", sizeof(int), FMOffset(offset_struct*, local_dimensions)},
    {"local_offsets", "integer[total_offsets]", sizeof(int), FMOffset(offset_struct*, local_offsets)},
    {NULL, NULL, 0, 0}
};

static FMStructDescRec offset_struct_format_list[] =
{
    {"offset_struct", offset_struct_field_list, sizeof(offset_struct), NULL},
    {NULL, NULL, 0, 0}
};

/********** global_var **********/
typedef struct _var {
    char * name;
    int noffset_structs;
    offset_struct * offsets;    
} global_var, *global_var_ptr;

static FMField global_var_field_list[]=
{
    {"name", "string", sizeof(char*), FMOffset(global_var_ptr, name)},
    {"noffset_structs", "integer", sizeof(int), FMOffset(global_var_ptr, noffset_structs)},
    {"offsets", "offset_struct[noffset_structs]", sizeof(offset_struct), FMOffset(global_var_ptr, offsets)},
    {NULL, NULL, 0, 0}
};

/********** evgroup **********/
typedef struct _evgroup {    
    int condition;
    int step;
    int num_vars;
    global_var* vars;
} evgroup, *evgroup_ptr;

static FMField evgroup_field_list[]=
{
    {"condition", "integer", sizeof(int), FMOffset(evgroup_ptr, condition)},
    {"step", "integer", sizeof(int), FMOffset(evgroup_ptr, step)},
    {"num_vars", "integer", sizeof(int), FMOffset(evgroup_ptr, num_vars)},
    {"vars", "global_var[num_vars]", sizeof(global_var), FMOffset(evgroup_ptr, vars)},
    {NULL, NULL, 0, 0}
};

static FMStructDescRec evgroup_format_list[] =
{   
    {"evgroup", evgroup_field_list, sizeof(evgroup), NULL},
    {"offset_struct", offset_struct_field_list, sizeof(offset_struct), NULL},
    {"global_var", global_var_field_list, sizeof(global_var), NULL},
    {NULL,NULL,0,NULL}
};

/********** op_msg **********/
typedef struct _op_msg
{
    int process_id;
    char * file_name;
    int type; //4 = end_of_stream, 3 = init, 2 = ack, 1 = open, 0 = close,
    int step;
    int condition;
} op_msg, *op_msg_ptr;

static FMField op_file_field_list[] =
{
    {"process_id", "integer", sizeof(int), FMOffset(op_msg_ptr, process_id)},
    {"file_name", "string", sizeof(char*), FMOffset(op_msg_ptr, file_name)},
    {"type", "integer", sizeof(int), FMOffset(op_msg_ptr, type)},
    {"step", "integer", sizeof(int), FMOffset(op_msg_ptr, step)},
    {"condition", "integer", sizeof(int), FMOffset(op_msg_ptr, condition)},
    {NULL, NULL, 0, 0}
};

static FMStructDescRec op_format_list[] =
{
    {"op_msg", op_file_field_list, sizeof(op_msg), NULL},
    {NULL, NULL, 0, NULL}
};

/********** Update step msg **********/
typedef struct _update_step_msg // one way message from writer to reader.
{
    int process_id; //mpi rank of the writer.
    int step;
    int finalized; // has the writer called finalize? 1 = yes.
} update_step_msg, *update_step_msg_ptr;

static FMField update_step_msg_field_list[] =
{
    {"process_id", "integer", sizeof(int), FMOffset(update_step_msg_ptr, process_id)},
    {"step" "integer", sizeof(int), FMOffset(update_step_msg_ptr, step)},
    {"finalized", "integer", sizeof(int), FMOffset(update_step_msg_ptr, finalized)},
    {NULL, NULL, 0, 0}
};

static FMStructDescRec update_step_msg_format_list[] =
{
    {"update_step_msg", update_step_msg_field_list, sizeof(update_step_msg), NULL},
    {NULL, NULL, 0, NULL}
};

/********** flush_msg **********/
typedef struct flush_msg_ {
    Flush_type type;
    int rank;
    int condition;
} Flush_msg, *Flush_msg_ptr;

static FMField flush_field_list[] =
{   
    {"type", "integer", sizeof(Flush_type), FMOffset(Flush_msg_ptr, type)},
    {"rank", "integer", sizeof(int), FMOffset(Flush_msg_ptr, rank)},
    {"condition", "integer", sizeof(int), FMOffset(Flush_msg_ptr, condition)},
    {NULL, NULL, 0, 0}
};
 
static FMStructDescRec flush_format_list[] =
{   
    {"flush", flush_field_list, sizeof(Flush_msg), NULL},
    {NULL,NULL,0,NULL}
};

/********** format_msg **********/
typedef struct format_msg_ {
    int id_len;
    int rep_id_len;
    char* format_id;
    char* rep_id;
    int condition;
} Format_msg, *Format_msg_ptr;

static FMField format_field_list[] =
{
    {"id_len", "integer", sizeof(int), FMOffset(Format_msg_ptr, id_len)},
    {"rep_id_len", "integer", sizeof(int), FMOffset(Format_msg_ptr, rep_id_len)},
    {"format_id", "char[id_len]", sizeof(char), FMOffset(Format_msg_ptr, format_id)},
    {"rep_id", "char[rep_id_len]", sizeof(char), FMOffset(Format_msg_ptr, rep_id)},
    {"condition", "integer", sizeof(int), FMOffset(Format_msg_ptr, condition)},
    {NULL, NULL, 0, 0}
};

static FMStructDescRec format_format_list[] =
{
    {"formatMsg", format_field_list, sizeof(Format_msg), NULL},
    {NULL,NULL,0,NULL}
};

/********** var_msg **********/
typedef struct var_msg_ {
    char* var_name;
    int rank;
    int condition;
} Var_msg, *Var_msg_ptr;

static FMField var_field_list[] =
{
    {"var_name", "string", sizeof(char*), FMOffset(Var_msg_ptr, var_name)},
    {"rank", "integer", sizeof(int), FMOffset(Var_msg_ptr, rank)},
    {NULL, NULL, 0, 0}
};

static FMStructDescRec var_format_list[] =
{
    {"varMsg", var_field_list, sizeof(Var_msg), NULL},
    {NULL, NULL, 0, NULL}
};

/********** data_format **********/
static FMStructDescRec data_format_list[] =
{
    {"anonymous", NULL, 0, NULL},
    {NULL, NULL, 0, NULL}
};

static char *getFixedName(char *name);


#endif
