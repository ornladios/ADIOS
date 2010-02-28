/*
 * ADIOS is freely available under the terms of the BSD license described
 * in the COPYING file in the top level directory of this source distribution.
 *
 * Copyright (c) 2008 - 2009.  UT-BATTELLE, LLC. All rights reserved.
 */

#include <unistd.h>
#include <fcntl.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdio.h>

#include "config.h"
#include "mpi.h"
#include "adios.h"
#include "adios_types.h"
#include "adios_bp_v1.h"
#include "adios_transport_hooks.h"
#include "adios_internals.h"
#ifdef HAVE_NSSI
#include "nssi_client.h"
#include "adios_nssi_args.h"
#include "adios_nssi_config.h"
#include "nssi_logger.h"
#endif

#include "io_timer.h"


#define NUM_GP 24
void adios_nssi_end_iteration(
        struct adios_method_struct * method)
{
}
void adios_nssi_start_calculation(
        struct adios_method_struct * method)
{
}
void adios_nssi_stop_calculation(
        struct adios_method_struct * method)
{
}
void adios_nssi_get_write_buffer(
        struct adios_file_struct *f,
        struct adios_var_struct *v,
        uint64_t *size,
        void **buffer,
        struct adios_method_struct *method)
{
}

#define HAVE_NSSI

#ifndef HAVE_NSSI
void adios_nssi_init(
        const char *parameters,
        struct adios_method_struct * method)
{
}
void adios_nssi_finalize(
        int mype,
        struct adios_method_struct * method)
{
}
enum ADIOS_FLAG adios_nssi_should_buffer(
        struct adios_file_struct *f,
        struct adios_method_struct *method)
{
    return adios_flag_unknown;
}
int adios_nssi_open(
        struct adios_file_struct *f,
        struct adios_method_struct *method,
        void * comm)
{
    return -1;
}
void adios_nssi_close(
        struct adios_file_struct *f,
        struct adios_method_struct *method)
{
}
void adios_nssi_write(
        struct adios_file_struct *f,
        struct adios_var_struct *v,
        void *data,
        struct adios_method_struct *method)
{
}
void adios_nssi_read(
        struct adios_file_struct *f,
        struct adios_var_struct *v,
        void *buffer,
        uint64_t buffersize,
        struct adios_method_struct *method)
{
}
#else

///////////////////////////
// Datatypes
///////////////////////////
struct adios_nssi_data_struct
{
    int      fd;
    MPI_Comm group_comm;
    int      rank;
    int      size;

    void * comm; // temporary until moved from should_buffer to open
};

/* Need a struct to encapsulate var offset info.
 */
struct var_offset {
    char      opath[ADIOS_PATH_MAX];
    char      oname[ADIOS_PATH_MAX];
    void     *ovalue;
    uint64_t  osize;
};

/* list of variable offsets */
List var_offset_list;

///////////////////////////
// Global Variables
///////////////////////////
static int adios_nssi_initialized = 0;

static int default_svc;
static MPI_Comm ncmpi_collective_op_comm;
static int global_rank=-1;
static int collective_op_rank=-1;
static nssi_service *svcs;
struct adios_nssi_config nssi_cfg;

//static log_level adios_nssi_debug_level;


///////////////////////////
// Function Declarations
///////////////////////////


///////////////////////////
// Function Definitions
///////////////////////////

static struct var_offset *var_offset_create(const char *path, const char *name, void *value, uint64_t size)
{
    struct var_offset *vo=malloc(sizeof(struct var_offset));

    strcpy(vo->opath, path);
    strcpy(vo->oname, name);
    vo->ovalue = value;
    vo->osize  = size;

    return(vo);
}
static void var_offset_free(void *vo)
{
    free(vo);
}
static int var_offset_equal(const struct var_offset *vo1, const struct var_offset *vo2)
{
    if ((strcmp(vo1->opath, vo2->opath) == 0) && (strcmp(vo1->oname, vo2->oname) == 0)) return TRUE;

    return FALSE;
}
static struct var_offset *var_offset_find(const char *path, const char *name)
{
    ListElmt *elmt;
    struct var_offset *vo;

    //printf("looking for opath(%s) oname(%s)\n", path, name);

    elmt = list_head(&var_offset_list);
    while(elmt) {
        vo = list_data(elmt);
        //printf("comparing to opath(%s) oname(%s)\n", vo->opath, vo->oname);
        if ((strcmp(path, vo->opath) == 0) && (strcmp(name, vo->oname) == 0)) {
            //printf("opath(%s) oname(%s) matches search\n", vo->opath, vo->oname);
            return vo;
        }
        elmt = list_next(elmt);
    }

    return NULL;
}
static void var_offset_printall(void)
{
    ListElmt *elmt;
    struct var_offset *vo;

    elmt = list_head(&var_offset_list);
    while(elmt) {
        vo = list_data(elmt);
        printf("opath(%s) oname(%s)\n", vo->opath, vo->oname);
        elmt = list_next(elmt);
    }

    return NULL;
}

static void parse_dimension_size(
        struct adios_group_struct *group,
        struct adios_var_struct *pvar_root,
        struct adios_attribute_struct *patt_root,
        struct adios_dimension_item_struct *dim,
        size_t *dimsize)
{
    struct adios_var_struct *var_linked = NULL;
    struct adios_attribute_struct *attr_linked;
    if (dim->id) {
        var_linked = adios_find_var_by_id (pvar_root , dim->id);
        if (!var_linked) {
            attr_linked = adios_find_attribute_by_id (patt_root, dim->id);
            if (!attr_linked->var) {
                switch (attr_linked->type) {
                case adios_unsigned_byte:
                    *dimsize = *(uint8_t *)attr_linked->value;
                    break;
                case adios_byte:
                    *dimsize = *(int8_t *)attr_linked->value;
                    break;
                case adios_unsigned_short:
                    *dimsize = *(uint16_t *)attr_linked->value;
                    break;
                case adios_short:
                    *dimsize = *(int16_t *)attr_linked->value;
                    break;
                case adios_unsigned_integer:
                    *dimsize = *(uint32_t *)attr_linked->value;
                    break;
                case adios_integer:
                    *dimsize = *(int32_t *)attr_linked->value;
                    break;
                case adios_unsigned_long:
                    *dimsize = *(uint64_t *)attr_linked->value;
                    break;
                case adios_long:
                    *dimsize = *(int64_t *)attr_linked->value;
                    break;
                default:
                    fprintf (stderr, "Invalid datatype for array dimension on "
                            "var %s: %s\n"
                            ,attr_linked->name
                            ,adios_type_to_string_int (var_linked->type)
                    );
                    break;
                }
            } else {
                var_linked = attr_linked->var;
            }
        }
        if (var_linked && var_linked->data) {
            *dimsize = *(int *)var_linked->data;
        }
    } else {
        if (dim->time_index == adios_flag_yes) {
            *dimsize = 1;
        } else {
            *dimsize = dim->rank;
        }
    }

    return;
}
static void parse_dimension_name(
        struct adios_group_struct *group,
        struct adios_var_struct *pvar_root,
        struct adios_attribute_struct *patt_root,
        struct adios_dimension_item_struct *dim,
        char *dimname)
{
    struct adios_var_struct *var_linked = NULL;
    struct adios_attribute_struct *attr_linked;
    if (dim->id) {
        var_linked = adios_find_var_by_id (pvar_root , dim->id);
        if (!var_linked) {
            attr_linked = adios_find_attribute_by_id (patt_root, dim->id);
            if (!attr_linked->var) {
//				strcpy(dimname, attr_linked->name);
                sprintf(dimname, "%s", attr_linked->name);
            } else {
                var_linked = attr_linked->var;
            }
        }
        if (var_linked && var_linked->name) {
//			strcpy(dimname, var_linked->name);
            sprintf(dimname, "%s", var_linked->name);
        }
    } else {
        if (dim->time_index == adios_flag_yes) {
//			strcpy(dimname, group->time_index_name);
            sprintf(dimname, "%s", group->time_index_name);
        } else {
            dimname[0] = '\0';
        }
    }

    return;
}


static int gen_offset_list(
        struct adios_group_struct *group,
        struct adios_var_struct *pvar_root,
        struct adios_attribute_struct *patt_root)
{
    struct adios_var_struct *v;
    struct adios_dimension_struct *dims;
    struct var_offset *vo;
    char offset_name[255];
    uint64_t *offset;

    v = pvar_root;
    while (v) {
        dims=v->dimensions;
        int loffs_idx=0;
        while (dims) {
            parse_dimension_name(group, pvar_root, patt_root, &dims->local_offset, offset_name);
            if (offset_name[0] == '\0') {
                sprintf(offset_name, "offset_%d", loffs_idx);
            }
            //printf("gen: offset_name(%s)\n", offset_name);
            offset=malloc(sizeof(uint64_t));
            parse_dimension_size(group, pvar_root, patt_root, &dims->local_offset, offset);
//			if (myrank==0) {
//				printf(":o(%d)", nc4_offsets[loffs_idx]);
//			}

            vo = var_offset_create("" /*v->path*/, offset_name, offset, sizeof(uint64_t));
            list_ins_next(&var_offset_list, list_head(&var_offset_list), vo);

            loffs_idx++;
            dims = dims->next;
        }
        v = v->next;
    }
}

static void create_offset_list_for_var(
        struct adios_write_args *args,
        struct adios_var_struct *v,
        struct adios_group_struct *group,
        struct adios_var_struct *pvar_root,
        struct adios_attribute_struct *patt_root)
{
    struct adios_dimension_struct *dims;
    struct var_offset *vo;
    char offset_name[255];
    uint64_t offset;

    args->offsets.offsets_len=0;
    args->offsets.offsets_val=NULL;

    if ((v) && (v->dimensions)) {
        int local_offset_count=0;
        dims=v->dimensions;
        while (dims) {
            local_offset_count++;
            dims = dims->next;
        }

        args->offsets.offsets_len=local_offset_count;
        args->offsets.offsets_val=malloc(local_offset_count * sizeof(struct adios_offset));

        dims=v->dimensions;
        int loffs_idx=0;
        while (dims) {
            parse_dimension_name(group, pvar_root, patt_root, &dims->local_offset, offset_name);
            if (offset_name[0] == '\0') {
                sprintf(offset_name, "offset_%d", loffs_idx);
            }
            //printf("create: offset_name(%s)\n", offset_name);
            args->offsets.offsets_val[loffs_idx].vpath=strdup("" /*v->path*/);
            args->offsets.offsets_val[loffs_idx].vname=strdup(offset_name);
//			if (myrank==0) {
//				printf(":o(%d)", nc4_offsets[loffs_idx]);
//			}

            loffs_idx++;
            dims = dims->next;
        }
    }
}

static int read_var(
        int fd,
        struct adios_var_struct *pvar,
        int myrank,
        int nproc)
{
    int return_code=0;
    int i, rc;

    adios_read_args args;
    adios_read_res  res;

    args.fd       = fd;
    args.max_read = pvar->data_size;
    args.vpath = strdup(pvar->path);
    args.vname = strdup(pvar->name);

    rc = nssi_call_rpc_sync(&svcs[default_svc],
            ADIOS_READ_OP,
            &args,
            pvar->data,
            pvar->data_size,
            &res);
    if (rc != NSSI_OK) {
        //log_error(adios_nssi_debug_level, "unable to call remote adios_read");
        return_code=-2;
    }

    free(args.vpath);
    free(args.vname);

    return return_code;
}
static int write_var(
        int fd,
        struct adios_group_struct *group,
        struct adios_var_struct *pvar_root,
        struct adios_attribute_struct *patt_root,
        struct adios_var_struct *pvar,
        uint64_t var_size,
        enum ADIOS_FLAG fortran_flag,
        int myrank,
        int nproc)
{
    int i;
    int rc;
    int return_code=0;

    adios_write_args args;
    adios_write_res  res;

//    var_offset_printall();

    args.fd    = fd;
    args.vsize = var_size;
    args.vpath = strdup(pvar->path);
    args.vname = strdup(pvar->name);
    args.is_scalar = FALSE;
    if (var_offset_find("" /*pvar->path*/, pvar->name)) {
        printf("var(%s) is an offset\n", pvar->name);
        args.is_offset = TRUE;
    } else {
        printf("var(%s) is NOT an offset\n", pvar->name);
        args.is_offset = FALSE;
    }
    args.writer_rank=myrank;
    args.offsets.offsets_len=0;
    args.offsets.offsets_val=NULL;
    if (pvar->dimensions) {
        create_offset_list_for_var(
                &args,
                pvar,
                group,
                group->vars,
                group->attributes);
    }

    MPI_Barrier(MPI_COMM_WORLD);
    rc = nssi_call_rpc_sync(&svcs[default_svc],
            ADIOS_WRITE_OP,
            &args,
            pvar->data,
            var_size,
            &res);
    if (rc != NSSI_OK) {
        //log_error(adios_nssi_debug_level, "unable to call remote adios_write");
        return_code=-2;
    }
    MPI_Barrier(MPI_COMM_WORLD);

    free(args.vpath);
    free(args.vname);
    free(args.offsets.offsets_val);

    return return_code;
}


static void adios_var_to_comm_nssi(
        enum ADIOS_FLAG host_language_fortran,
        void *data,
        MPI_Comm *comm)
{
    if (data) {
        int t = *(int *) data;
        if (host_language_fortran == adios_flag_yes) {
            *comm = MPI_Comm_f2c (t);
        } else {
            *comm = *(MPI_Comm *) data;
        }
    } else {
        fprintf (stderr, "coordination-communication not provided. "
                "Using MPI_COMM_WORLD instead\n");
        *comm = MPI_COMM_WORLD;
    }

    return;
}


void adios_nssi_init(
        const char *parameters,
        struct adios_method_struct *method)
{
    int rc=NSSI_OK;
    struct adios_nssi_data_struct *md = (struct adios_nssi_data_struct *)method->method_data;
    int verbose=3;

    if (!adios_nssi_initialized) {
        adios_nssi_initialized = 1;
    }
    method->method_data = malloc(sizeof(struct adios_nssi_data_struct));
    md = (struct adios_nssi_data_struct *)method->method_data;
    md->fd         = -1;
    md->rank       = -1;
    md->size       = 0;
    md->group_comm = MPI_COMM_NULL;

    logger_init((log_level)verbose, NULL);

#ifdef HAVE_PORTALS
    nssi_ptl_init(PTL_IFACE_CLIENT, getpid() + 1000);
    nssi_rpc_init(NSSI_RPC_PTL, NSSI_RPC_XDR);
#endif
#ifdef HAVE_INFINIBAND
    nssi_ib_init(NULL);
    rc = nssi_rpc_init(NSSI_RPC_IB, NSSI_RPC_XDR);
#endif

    /* Register the client operations */
    NSSI_REGISTER_CLIENT_STUB(ADIOS_OPEN_OP, adios_open_args, void, adios_open_res);
    NSSI_REGISTER_CLIENT_STUB(ADIOS_GROUP_SIZE_OP, adios_group_size_args, void, void);
    NSSI_REGISTER_CLIENT_STUB(ADIOS_READ_OP, adios_read_args, void, adios_read_res);
    NSSI_REGISTER_CLIENT_STUB(ADIOS_WRITE_OP, adios_write_args, void, adios_write_res);
    NSSI_REGISTER_CLIENT_STUB(ADIOS_END_ITER_OP, adios_end_iter_args, void, void);
    NSSI_REGISTER_CLIENT_STUB(ADIOS_START_CALC_OP, adios_start_calc_args, void, void);
    NSSI_REGISTER_CLIENT_STUB(ADIOS_STOP_CALC_OP, adios_stop_calc_args, void, void);
    NSSI_REGISTER_CLIENT_STUB(ADIOS_CLOSE_OP, adios_close_args, void, void);

    list_init(&var_offset_list, var_offset_free);

    parse_nssi_config(getenv("ADIOS_NSSI_CONFIG_FILE"), &nssi_cfg);

    return;
}


enum ADIOS_FLAG adios_nssi_should_buffer(
        struct adios_file_struct *f,
        struct adios_method_struct *method)
{
    int rc=NSSI_OK;

    struct adios_nssi_data_struct *md = (struct adios_nssi_data_struct *)method->method_data;

    adios_group_size_args args;

    args.fd = md->fd;
    args.data_size = f->write_size_bytes;

    rc = nssi_call_rpc_sync(&svcs[default_svc],
            ADIOS_GROUP_SIZE_OP,
            &args,
            NULL,
            0,
            NULL);

    return adios_flag_no;
}

int adios_nssi_open(
        struct adios_file_struct *f,
        struct adios_method_struct *method,
        void *comm)
{
    int rc=NSSI_OK;

    struct adios_nssi_data_struct * md = (struct adios_nssi_data_struct *)method->method_data;

    adios_open_args args;
    adios_open_res  res;

    if (md->fd != -1) {
        // file already open
        return adios_flag_no;
    }

    md->comm = comm;
    adios_var_to_comm_nssi(f->group->adios_host_language_fortran, md->comm, &md->group_comm);
    if (md->group_comm != MPI_COMM_NULL) {
        MPI_Comm_rank(md->group_comm, &md->rank);
        MPI_Comm_size(md->group_comm, &md->size);
    } else {
        md->group_comm=MPI_COMM_SELF;
    }
    f->group->process_id = md->rank;

    if (md->size <= nssi_cfg.num_servers) {
        default_svc = md->rank;
    } else {
        if ((md->size%nssi_cfg.num_servers) > 0) {
            default_svc = md->rank/((md->size/nssi_cfg.num_servers)+1);
        } else {
            default_svc = md->rank/(md->size/nssi_cfg.num_servers);
        }
    }

    /* create a new communicator for just those clients, who share a default service. */
    MPI_Comm_split(MPI_COMM_WORLD, default_svc, md->rank, &ncmpi_collective_op_comm);
    /* find my rank in the new communicator */
    MPI_Comm_rank(ncmpi_collective_op_comm, &collective_op_rank);

    //log_debug(adios_nssi_debug_level, "global_rank(%d) collective_op_rank(%d) default_service(%d)", md->rank, collective_op_rank, default_svc);

    svcs=(nssi_service *)calloc(nssi_cfg.num_servers, sizeof(nssi_service));
    /* !global_rank0 has a preferred server for data transfers.  connect to preferred server.
     * connect to other servers on-demand.
     */
    double GetSvcTime=MPI_Wtime();
    printf("get staging-service: default_svc(%d) nid(%lld) pid(%llu) hostname(%s) port(%d)\n",
            default_svc,
            nssi_cfg.nssi_server_ids[default_svc].nid,
            nssi_cfg.nssi_server_ids[default_svc].pid,
            nssi_cfg.nssi_server_ids[default_svc].hostname,
            nssi_cfg.nssi_server_ids[default_svc].port);
    rc = nssi_get_service(nssi_cfg.nssi_server_ids[default_svc], -1, &svcs[default_svc]);
    if (rc != NSSI_OK) {
        //log_error(adios_nssi_debug_level, "Couldn't connect to netcdf master: %s", nssi_err_str(rc));
        return;
    }


    gen_offset_list(
            f->group,
            f->group->vars,
            f->group->attributes);
//    var_offset_printall();

    args.fname = malloc(sizeof(char) * (strlen(method->base_path) + strlen(f->name) + 1));
    sprintf(args.fname, "%s%s", method->base_path, f->name);
    args.gname = strdup(method->group->name);
    switch(f->mode) {
        case adios_mode_read:
            args.mode = ADIOS_MODE_READ;
            break;
        case adios_mode_write:
            args.mode = ADIOS_MODE_WRITE;
            break;
        case adios_mode_append:
            args.mode = ADIOS_MODE_APPEND;
            break;
        case adios_mode_update:
            args.mode = ADIOS_MODE_UPDATE;
            break;
    }

    rc = nssi_call_rpc_sync(&svcs[default_svc],
            ADIOS_OPEN_OP,
            &args,
            NULL,
            0,
            &res);

    md->fd = res.fd;

    free(args.fname);
    free(args.gname);

    return 1;
}

void adios_nssi_write(
        struct adios_file_struct *f,
        struct adios_var_struct *v,
        void *data,
        struct adios_method_struct *method)
{
    struct adios_nssi_data_struct * md = (struct adios_nssi_data_struct *)method->method_data;
    static int first_write = 1;

    if (f->mode == adios_mode_write || f->mode == adios_mode_append) {
//		if (first_write == 1) {
//			write_header(fd, md);
//			first_write = 0;
//		}

        if (md->rank==0) {
//			fprintf(stderr, "-------------------------\n");
//			fprintf(stderr, "write var: %s start!\n", v->name);
        }
        uint64_t var_size = adios_get_var_size (v, f->group, data);
        printf("vname(%s) vsize(%ld)\n", v->name, var_size);
        write_var(md->fd,
                f->group,
                f->group->vars,
                f->group->attributes,
                v,
                var_size,
                f->group->adios_host_language_fortran,
                md->rank,
                md->size);
    } else {
        //fprintf(stderr, "entering unknown nc4 mode %d!\n", fd->mode);
    }
    if (md->rank==0) {
//		fprintf(stderr, "write var: %s end!\n", v->name);
        //fprintf(stderr, "-------------------------\n");
    }

    return;
}


void adios_nssi_read(
        struct adios_file_struct *f,
        struct adios_var_struct *v,
        void *buffer,
        uint64_t buffersize,
        struct adios_method_struct *method)
{
    struct adios_nssi_data_struct * md = (struct adios_nssi_data_struct *)method->method_data;

    if(f->mode == adios_mode_read) {
        v->data = buffer;
        v->data_size = buffersize;

        if (md->rank==0) {
//			fprintf(stderr, "-------------------------\n");
//			fprintf(stderr, "read var: %s! start\n", v->name);
        }
        read_var(md->fd,
                v,
                md->rank,
                md->size);
        if (md->rank==0) {
//			fprintf(stderr, "read var: %s! end\n", v->name);
            //fprintf(stderr, "-------------------------\n");
        }
    }

    return;
}

static void adios_nssi_do_read(
        struct adios_file_struct *f,
        struct adios_method_struct *method)
{
    // This function is not useful for nc4 since adios_read/write do real read/write
}

void adios_nssi_close(
        struct adios_file_struct *f,
        struct adios_method_struct *method)
{
    int rc=NSSI_OK;
    struct adios_nssi_data_struct * md = (struct adios_nssi_data_struct*)method->method_data;
    struct adios_attribute_struct * a = f->group->attributes;
    int myrank=md->rank;

    adios_close_args args;

    if (f->mode == adios_mode_read) {
        if (md->rank==0) {
            fprintf(stderr, "-------------------------\n");
            fprintf(stderr, "reading done, nc4 file is virtually closed;\n");
            fprintf(stderr, "-------------------------\n");
        }
    } else if (f->mode == adios_mode_write || f->mode == adios_mode_append) {
        //fprintf(stderr, "entering nc4 write attribute mode!\n");
        if (md->rank==0) {
            fprintf(stderr, "-------------------------\n");
            fprintf(stderr, "writing done, nc4 file is virtually closed;\n");
            fprintf(stderr, "-------------------------\n");
        }
    }

    MPI_Barrier(md->group_comm);
    if (md->rank == 0) {
        args.fd = md->fd;
        rc = nssi_call_rpc_sync(&svcs[default_svc],
                ADIOS_CLOSE_OP,
                &args,
                NULL,
                0,
                NULL);
        if (rc != NSSI_OK) {
            //log_error(adios_nssi_debug_level, "unable to call remote adios_read");
        }
    }
    MPI_Barrier(md->group_comm);

    return;
}

void adios_nssi_finalize(
        int mype,
        struct adios_method_struct *method)
{
    struct adios_nssi_data_struct * md = (struct adios_nssi_data_struct*)method->method_data;
    int myrank=md->rank;

    md->group_comm = MPI_COMM_NULL;
    md->fd = -1;
    md->rank = -1;
    md->size = 0;

    free_nssi_config(&nssi_cfg);

    if (adios_nssi_initialized)
        adios_nssi_initialized = 0;
}

#endif
