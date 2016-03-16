/*
 * ADIOS is freely available under the terms of the BSD license described
 * in the COPYING file in the top level directory of this source distribution.
 *
 * Copyright (c) 2008 - 2009.  UT-BATTELLE, LLC. All rights reserved.
 */


/**************************************************/
/* Read method for NSSI memory-to-memory coupling */
/**************************************************/

#include <stdlib.h>
#include <string.h>
//#include <errno.h>  /* ENOMEM */
#include "adios.h"
#include "bp_utils.h"
#include "bp_types.h"
#include "adios_types.h"
#include "adios_read.h"
#include "adios_read_hooks.h"
#include "adios_error.h"
#include "futils.h"
#include "globals.h"

#ifdef HAVE_NSSI
#include "nssi_client.h"
#include "adios_nssi_args.h"
#include "adios_nssi_config.h"
#include "nssi_logger.h"
#endif

#include "io_timer.h"


static uint64_t number_of_fopens = 0;  /* for versioning, works only if one file is fopened (in a loop) in the application */

struct adios_read_nssi_data_struct
{
    nssi_request start_calc_req;
    int          has_outstanding_req;
    int          default_svc_index;  /* service to use when there is no open file (eg. finalize) */

    char *fname;               // path of file
    int timestep;              // counting the access
    int disconnect_at_fclose;  // disconnect from NSSI in fclose()

    uint64_t fd;
    MPI_Comm group_comm;
    int      size;
    int      rank;

    int      svc_index;
    MPI_Comm collective_op_comm;
    int      collective_op_size;
    int      collective_op_rank;

    int8_t use_single_server;
};


///////////////////////////
// Global Variables
///////////////////////////
static int adios_nssi_initialized = 0;

static char *job_id=NULL;
static int global_rank=-1;
static int global_size=-1;
static nssi_service *svcs;
struct adios_nssi_config nssi_cfg;

//static log_level adios_nssi_debug_level;
static int DEBUG=0;




/* If init is used, we connect to NSSI here, otherwise we connect in fopen.
   If multiple fopen..fclose cycles are used, init/finalize must be used too to
   avoid multiple connection/disconnection in fopen/fclose.
*/
int adios_read_nssi_init (MPI_Comm comm)
{
    int rc=NSSI_OK;
    int verbose=5;
    char logfile[1024];
    int log_rank;
    struct adios_read_nssi_method_data_struct *private;

    if (!adios_nssi_initialized) {
        adios_nssi_initialized = 1;
    }

    MPI_Comm_rank(MPI_COMM_WORLD, &global_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &global_size);

    if (DEBUG>3) fprintf(stderr, "rank(%d) enter adios_read_nssi_init\n", global_rank);

#ifdef HAVE_PORTALS
    nssi_ptl_init(PTL_IFACE_CLIENT, getpid() + 1000);
    nssi_rpc_init(NSSI_RPC_PTL, NSSI_RPC_XDR);
#endif
#ifdef HAVE_INFINIBAND
    nssi_ib_init(NULL);
    rc = nssi_rpc_init(NSSI_RPC_IB, NSSI_RPC_XDR);
#endif

    /* Register the client operations */
//    NSSI_REGISTER_CLIENT_STUB(ADIOS_READ_INITIALIZE_OP,    adios_read_initialize_args,    void, void);
    NSSI_REGISTER_CLIENT_STUB(ADIOS_READ_FOPEN_OP,            adios_read_fopen_args,            void, adios_read_fopen_res);
    NSSI_REGISTER_CLIENT_STUB(ADIOS_READ_FCLOSE_OP,           adios_read_fclose_args,           void, void);
//    NSSI_REGISTER_CLIENT_STUB(ADIOS_READ_GOPEN_OP,         adios_read_gopen_args,         void, adios_read_gopen_res);
//    NSSI_REGISTER_CLIENT_STUB(ADIOS_READ_GOPEN_BYID_OP,    adios_read_gopen_byid_args,    void, adios_read_gopen_byid_res);
//    NSSI_REGISTER_CLIENT_STUB(ADIOS_READ_GCLOSE_OP,        adios_read_gclose_args,        void, adios_read_gclose_res);
//    NSSI_REGISTER_CLIENT_STUB(ADIOS_READ_GETATTR_OP,       adios_read_getattr_args,       void, adios_read_getattr_res);
//    NSSI_REGISTER_CLIENT_STUB(ADIOS_READ_GETATTR_BYID_OP,  adios_read_getattr_byid_args,  void, adios_read_getattr_byid_res);
//    NSSI_REGISTER_CLIENT_STUB(ADIOS_READ_INQVAR_OP,        adios_read_inqvar_args,        void, adios_read_inqvar_res);
//    NSSI_REGISTER_CLIENT_STUB(ADIOS_READ_INQVAR_BYID_OP,   adios_read_inqvar_byid_args,   void, adios_read_inqvar_byid_res);
    NSSI_REGISTER_CLIENT_STUB(ADIOS_READ_GET_VARTYPE_SIZE_OP, adios_read_get_vartype_size_args, void, adios_read_get_vartype_size_res);
    NSSI_REGISTER_CLIENT_STUB(ADIOS_READ_READ_VAR_OP,         adios_read_read_var_args,         void, adios_read_read_var_res);
//    NSSI_REGISTER_CLIENT_STUB(ADIOS_READ_READ_VAR_BYID_OP, adios_read_read_var_byid_args, void, adios_read_read_var_byid_res);
//    NSSI_REGISTER_CLIENT_STUB(ADIOS_READ_FINALIZE_OP,      adios_read_finalize_args,      void, void);

    parse_nssi_config(getenv("ADIOS_NSSI_CONFIG_FILE"), &nssi_cfg);

    return(0);
}

int adios_read_nssi_finalize ()
{
    int rc=NSSI_OK;
    int myrank;

    if (DEBUG>3) fprintf(stderr, "rank(%d) enter adios_read_nssi_finalize\n", global_rank);

    free_nssi_config(&nssi_cfg);

    if (adios_nssi_initialized)
        adios_nssi_initialized = 0;


//    // disconnect from NSSI only if we the reader is connected (the writer not anymore)
//    if (globals_adios_is_nssi_connected_from_reader() &&
//        !globals_adios_is_nssi_connected_from_both())
//    {
//        nssi_finalize();
//        DBG_PRINTF("-- %s: disconnected from NSSI\n", __func__);
//    }
//    globals_adios_set_nssi_disconnected_from_reader();
}

ADIOS_FILE * adios_read_nssi_fopen (const char * fname, MPI_Comm comm)
{
    int rc=NSSI_OK;
    ADIOS_FILE * fp;
    struct adios_read_nssi_data_struct * ds;
    int i;

    adios_read_fopen_args args;
    adios_read_fopen_res  res;

    adios_errno = 0;

    if (DEBUG>3) fprintf(stderr, "enter adios_read_nssi_fopen: fname=%s\n", fname);

    ds = (struct adios_read_nssi_data_struct *) malloc (sizeof(struct adios_read_nssi_data_struct));
    if (!ds) {
        adios_error (err_no_memory, "Cannot allocate memory for file info.");
        return(NULL);
    }

    /* fill out NSSI method specific struct */
    ds->size = 0;
    ds->rank = 0;
    ds->fname = strdup(fname);
    ds->timestep = 1; /*number_of_fopens+1;  /* Read data of separate versions from NSSI */

    if (DEBUG>3) fprintf(stderr, "global_rank(%d): enter adios_read_nssi_fopen (%s)\n", global_rank, fname);

    ds->group_comm = comm;
    if (DEBUG>3) fprintf(stderr, "global_rank(%d): adios_read_nssi_fopen: setup group_comm\n", global_rank);
    if (ds->group_comm != MPI_COMM_NULL) {
        if (DEBUG>3) fprintf(stderr, "global_rank(%d): adios_read_nssi_fopen: get rank and size\n", global_rank);
        MPI_Comm_rank(ds->group_comm, &ds->rank);
        MPI_Comm_size(ds->group_comm, &ds->size);
        if (DEBUG>3) fprintf(stderr, "global_rank(%d): adios_read_nssi_fopen: size(%d) rank(%d)\n", global_rank, ds->size, ds->rank);
    } else {
        ds->group_comm=MPI_COMM_SELF;
        MPI_Comm_rank(ds->group_comm, &ds->rank);
        MPI_Comm_size(ds->group_comm, &ds->size);
    }

    if (ds->size <= nssi_cfg.num_servers) {
        // there are fewer clients than servers.
        // assume file-per-process and use a single server for this file.
        ds->use_single_server=TRUE;
        if (ds->size < global_size) {
            // a subset of all clients is writing
            ds->svc_index = ((global_rank/ds->size)%nssi_cfg.num_servers);
        } else {
            ds->svc_index = 0;
        }
    } else {
        ds->use_single_server=FALSE;
        if ((ds->size%nssi_cfg.num_servers) > 0) {
            ds->svc_index = ds->rank/((ds->size/nssi_cfg.num_servers)+1);
        } else {
            ds->svc_index = ds->rank/(ds->size/nssi_cfg.num_servers);
        }
    }
    if (ds->default_svc_index == -1) {
        ds->default_svc_index=ds->svc_index;
    }

    /* create a new communicator for just those clients, who share a default service. */
    if (DEBUG>3) fprintf(stderr, "global_rank(%d): adios_read_nssi_fopen: before MPI_Comm_split\n", global_rank);
    MPI_Comm_split(ds->group_comm, ds->svc_index, ds->rank, &ds->collective_op_comm);
    if (DEBUG>3) fprintf(stderr, "global_rank(%d): adios_read_nssi_fopen: after MPI_Comm_split\n", global_rank);
    /* find my rank in the new communicator */
    if (DEBUG>3) fprintf(stderr, "global_rank(%d): adios_read_nssi_fopen: before MPI_Comm_size\n", global_rank);
    MPI_Comm_size(ds->collective_op_comm, &ds->collective_op_size);
    if (DEBUG>3) fprintf(stderr, "global_rank(%d): adios_read_nssi_fopen: before MPI_Comm_rank\n", global_rank);
    MPI_Comm_rank(ds->collective_op_comm, &ds->collective_op_rank);

    if (DEBUG>3) fprintf(stderr, "global_rank(%d) ds->rank(%d) ds->collective_op_rank(%d) default_service(%d)\n", global_rank, ds->rank, ds->collective_op_rank, ds->svc_index);

    svcs=(nssi_service *)calloc(nssi_cfg.num_servers, sizeof(nssi_service));
    /* !global_rank0 has a preferred server for data transfers.  connect to preferred server.
     * connect to other servers on-demand.
     */
    double GetSvcTime=MPI_Wtime();
    if (DEBUG>3) fprintf(stderr, "get staging-service: ds->svc_index(%d) nid(%lld) pid(%llu) hostname(%s) port(%d)\n",
            ds->svc_index,
            nssi_cfg.nssi_server_ids[ds->svc_index].nid,
            nssi_cfg.nssi_server_ids[ds->svc_index].pid,
            nssi_cfg.nssi_server_ids[ds->svc_index].hostname,
            nssi_cfg.nssi_server_ids[ds->svc_index].port);
    rc = nssi_get_service(nssi_cfg.nssi_server_ids[ds->svc_index], -1, &svcs[ds->svc_index]);
    if (rc != NSSI_OK) {
        fprintf(stderr, "NSSI ERROR: nssi_get_service failed\n");
        return(NULL);
    }

    if (job_id==NULL) {
        if (ds->rank==0) {
            job_id = getenv("PBS_JOBID");
            if (job_id == NULL) {
                fprintf(stderr, "adios_read_nssi_init: unable to determine job id.  defaulting id to \"UNKNOWN_JOB_ID\".\n");
                job_id = strdup("UNKNOWN_JOB_ID");
            } else {
                int len=strlen(job_id)+36+1;
                job_id=calloc(len,1);

                struct uuid_st;
                extern int uuid_create   (      struct uuid_st **_uuid);
                extern int uuid_destroy  (      struct uuid_st  *_uuid);
                extern int uuid_make     (      struct uuid_st  *_uuid, unsigned int _mode, ...);
                extern int uuid_export   (const struct uuid_st  *_uuid, unsigned int _fmt,       void **_data_ptr, size_t *_data_len);

                struct uuid_st *uuid;
                char *uuid_str=NULL;
                uuid_create(&uuid);
                uuid_make(uuid, 1);
                uuid_export(uuid, 1, &uuid_str, NULL);
                uuid_destroy(uuid);

                sprintf(job_id, "%s.%s", getenv("PBS_JOBID"), uuid_str);

                free(uuid_str);

                MPI_Bcast (job_id, len, MPI_BYTE, 0, ds->group_comm);
            }
        } else {
            job_id = getenv("PBS_JOBID");
            if (job_id == NULL) {
                fprintf(stderr, "adios_read_nssi_init: unable to determine job id.  defaulting id to \"UNKNOWN_JOB_ID\".\n");
                job_id = strdup("UNKNOWN_JOB_ID");
            } else {
                int len=strlen(job_id)+36+1;
                job_id=calloc(len,1);

                MPI_Bcast (job_id, len, MPI_BYTE, 0, ds->group_comm);
            }
        }
    }

    fprintf(stderr, "adios_read_nssi_fopen: job_id=%s\n", job_id);



    memset(&args, 0, sizeof(args));
    memset(&res,  0, sizeof(res));
    args.client_id = strdup(job_id);
    args.gname = strdup("");
    args.fname = strdup(ds->fname);
    args.requested_timestep = ds->timestep;

    Func_Timer("ADIOS_READ_FOPEN_OP",
            rc = nssi_call_rpc_sync(&svcs[ds->svc_index],
            ADIOS_READ_FOPEN_OP,
            &args,
            NULL,
            0,
            &res););
    free(args.client_id);
    free(args.gname);
    free(args.fname);
    if (rc != NSSI_OK) {
        fprintf(stderr, "NSSI ERROR: ADIOS_READ_FOPEN_OP failed\n");
        return(NULL);
    }

    if (res.fd==-1) {
        fprintf(stderr, "ADIOS_READ_FOPEN_OP failed: fd=%ld\n", res.fd);
        return(NULL);
    }
    ds->fd = res.fd;

    fp = (ADIOS_FILE *) malloc (sizeof (ADIOS_FILE));
    if (!fp) {
        adios_error (err_no_memory, "Cannot allocate memory for file info.");
        return(NULL);
    }

    /* fill out ADIOS_FILE struct */
    fp->fh = (uint64_t) ds;
    fp->groups_count = 1;
    fp->vars_count = 0;
    fp->attrs_count = 0;
    fp->tidx_start = 0;
    fp->ntimesteps = 1;
    fp->file_size = 0;
    fp->version = 1;
    fp->endianness = 0; /* FIXME: not always Little Endian. Does it matter? */
    a2s_alloc_namelist (&fp->group_namelist,fp->groups_count);
    for (i=0;i<fp->groups_count;i++) {
        if (!fp->group_namelist[i]) {
            adios_error (err_no_memory, "Could not allocate buffer for %d strings in adios_fopen()", fp->groups_count);
            adios_read_nssi_fclose(fp);
            return(NULL);
        }
        else  {
            strcpy(fp->group_namelist[i],"NSSI");
        }
    }

    number_of_fopens++;

    if (DEBUG>3) fprintf(stderr, "exit adios_read_nssi_fopen: fname=%s fd=%ld timestep=%ld\n",
        fname, ds->fd, ds->timestep);

    return(fp);
}

int adios_read_nssi_fclose (ADIOS_FILE *fp)
{
    int rc=NSSI_OK;
    int i,j;

    adios_read_fclose_args args;

    struct adios_read_nssi_data_struct * ds = (struct adios_read_nssi_data_struct *) fp->fh;

    if (DEBUG>3) fprintf(stderr, "enter adios_read_nssi_fclose: fname=%s\n", ds->fname);

    adios_errno = 0;

    MPI_Barrier(ds->collective_op_comm);

    if (ds->collective_op_rank == 0) {
        memset(&args, 0, sizeof(args));
        args.client_id     = strdup(job_id);
        args.fd            = ds->fd;
        args.fname         = strdup(ds->fname);
        args.open_timestep = ds->timestep;

        Func_Timer("ADIOS_READ_FCLOSE_OP",
                rc = nssi_call_rpc_sync(&svcs[ds->svc_index],
                        ADIOS_READ_FCLOSE_OP,
                        &args,
                        NULL,
                        0,
                        NULL););
        free(args.client_id);
        free(args.fname);
        if (rc != NSSI_OK) {
            fprintf(stderr, "NSSI ERROR: ADIOS_READ_FCLOSE_OP failed\n");
        }
    }

    a2s_free_namelist ((fp->group_namelist),fp->groups_count);
    if (ds->fname) { free(ds->fname); ds->fname = 0; }
    free(ds);
    free(fp);

    if (DEBUG>3) fprintf(stderr, "exit adios_read_nssi_fclose: fname=%s\n", ds->fname);

    return 0;
}


int adios_read_nssi_get_dimension_order (const ADIOS_FILE *fp)
{
    return 0;
}

/* This function can be called if user places
   the wrong sequences of dims for a var
*/
void adios_read_nssi_reset_dimension_order (ADIOS_FILE *fp, int is_fortran)
{
    /* unimplemented */
}


ADIOS_GROUP * adios_read_nssi_gopen (ADIOS_FILE *fp, const char * grpname)
{
    /* NSSI has no groups, so any grpname is accepted and the same empty stuff is returned */
    return adios_read_nssi_gopen_byid(fp, 0);
}

ADIOS_GROUP * adios_read_nssi_gopen_byid (ADIOS_FILE *fp, int grpid)
{
    struct adios_read_nssi_data_struct * ds = (struct adios_read_nssi_data_struct *) fp->fh;
    ADIOS_GROUP * gp;

    /* NSSI has no groups, so any grpid is accepted and the same empty stuff is returned */

    adios_errno = 0;
    gp = (ADIOS_GROUP *) malloc(sizeof(ADIOS_GROUP));
    if (!gp) {
        adios_error (err_no_memory, "Could not allocate memory for group info");
        return NULL;
    }

    /* fill out ADIOS_GROUP struct */
    gp->grpid = grpid;
    gp->gh = (uint64_t) 0;
    gp->fp = fp;
    gp->vars_count = 0;
    gp->attrs_count = 0;
    gp->var_namelist = 0;
    gp->attr_namelist = 0;

    return gp;
}

int adios_read_nssi_gclose (ADIOS_GROUP *gp)
{
    struct adios_read_nssi_data_struct * ds = (struct adios_read_nssi_data_struct *) gp->fp->fh;

    adios_errno = 0;

    a2s_free_namelist ((gp->var_namelist),gp->vars_count);
    a2s_free_namelist ((gp->attr_namelist),gp->attrs_count);
    free(gp);
    return 0;
}



int adios_read_nssi_get_attr (ADIOS_GROUP * gp, const char * attrname, enum ADIOS_DATATYPES * type,
                    int * size, void ** data)
{
    /* NSSI does not support attributes */
    adios_error (err_invalid_attrname, "NSSI read method does not support attributes!");
    *size = 0;
    *type = adios_unknown;
    *data = 0;
    return adios_errno;
}

int adios_read_nssi_get_attr_byid (ADIOS_GROUP * gp, int attrid,
                    enum ADIOS_DATATYPES * type, int * size, void ** data)
{
    /* NSSI does not support attributes */
    adios_error (err_invalid_attrid, "NSSI read method does not support attributes!");
    *size = 0;
    *type = adios_unknown;
    *data = 0;
    return adios_errno;
}


ADIOS_VARINFO * adios_read_nssi_inq_var (ADIOS_GROUP *gp, const char * varname)
{
    /* NSSI has no inquiry capability, report somthing dummy */
    return adios_read_nssi_inq_var_byid(gp, 0);
}

ADIOS_VARINFO * adios_read_nssi_inq_var_byid (ADIOS_GROUP *gp, int varid)
{
    struct adios_read_nssi_data_struct * ds = (struct adios_read_nssi_data_struct *) gp->fp->fh;
    ADIOS_VARINFO * vi;
    int i,k;

    adios_errno = 0;
    vi = (ADIOS_VARINFO *) malloc(sizeof(ADIOS_VARINFO));
    if (!vi) {
        adios_error (err_no_memory, "Could not allocate memory for variable info.");
        return NULL;
    }

    /* NSSI has no inquiry capability, report somthing dummy */
    vi->varid = varid;
    vi->type = adios_unknown;
    vi->ndim = 0;
    vi->dims = NULL;
    vi->timedim = -1;
    vi->value = NULL;
    vi->gmin = NULL;
    vi->gmax = NULL;
    vi->mins = NULL;
    vi->maxs = NULL;

    return vi;
}

void adios_read_nssi_free_varinfo (ADIOS_VARINFO *vp)
{
    if (vp) {
        if (vp->dims)   free(vp->dims);
        if (vp->value)  free(vp->value);
        if (vp->gmin && vp->gmin != vp->value)   free(vp->gmin);
        if (vp->gmax && vp->gmax != vp->value)   free(vp->gmax);
        if (vp->mins)   free(vp->mins);
        if (vp->maxs)   free(vp->maxs);
        free(vp);
    }
}

//static int adios_read_nssi_get (const char * varname, enum ADIOS_DATATYPES vartype,
//                                struct adios_read_nssi_data_struct * ds,
//                                int * offset, int * readsize, void * data)
//{
//
//    struct obj_data *od;
//    int elemsize = common_read_type_size(vartype, NULL);
//    int err;
//
//    DBG_PRINTF("-- %s, rank %d: get data: varname=%s version=%d, lb=(%d,%d,%d) ub=(%d,%d,%d)}\n",
//        __func__, ds->rank, varname, ds->timestep, offset[1], offset[0], offset[2],
//        offset[1]+readsize[1]-1, offset[0]+readsize[0]-1, offset[2]+readsize[2]-1);
//
//    err =  NSSI_get (varname, ds->timestep, elemsize,
//                     offset[1], offset[0], offset[2],
//                     offset[1]+readsize[1]-1,
//                     offset[0]+readsize[0]-1,
//                     offset[2]+readsize[2]-1,
//                     data
//                    );
//    /*if (err == -ENOMEM) {
//        adios_error (err_no_memory, "Not enough memory for NSSI to perform NSSI_get()");
//        return -err_no_memory;
//    }
//    else*/ if (err) {
//        adios_error (err_corrupted_variable, "NSSI failed to read variable %s.", varname);
//        return -err_corrupted_variable;
//    }
//
//    return 0;
//}

int64_t adios_read_nssi_read_var (ADIOS_GROUP * gp, const char * varname,
                        const uint64_t * start, const uint64_t * count,
                        void * data)
{
    int rc=NSSI_OK;
    int64_t total_size;
    struct adios_read_nssi_data_struct * ds = (struct adios_read_nssi_data_struct *) gp->fp->fh;
    enum ADIOS_DATATYPES vartype;
    int elemsize;
    int err;
    int i;

    adios_read_get_vartype_size_args vts_args;
    adios_read_get_vartype_size_res  vts_res;
    adios_read_read_var_args args;
    adios_read_read_var_res  res;

    memset(&vts_args, 0, sizeof(vts_args));
    memset(&vts_res,  0, sizeof(vts_res));
    vts_args.fd            = ds->fd;
    vts_args.open_timestep = ds->timestep;
    vts_args.client_id     = strdup(job_id);
    vts_args.vpath         = strdup("");
    vts_args.vname         = strdup(varname);
    Func_Timer("ADIOS_READ_GET_VARTYPE_SIZE_OP",
            rc = nssi_call_rpc_sync(&svcs[ds->svc_index],
                    ADIOS_READ_GET_VARTYPE_SIZE_OP,
                    &vts_args,
                    NULL,
                    0,
                    &vts_res););
    free(vts_args.vpath);
    free(vts_args.vname);
    if (rc != NSSI_OK) {
        fprintf(stderr, "NSSI ERROR: ADIOS_READ_GET_VARTYPE_SIZE_OP failed\n");
        return(-1);
    }

    memset(&args, 0, sizeof(args));
    memset(&res,  0, sizeof(res));
    args.fd=ds->fd;
    args.open_timestep = ds->timestep;
    args.client_id = strdup(job_id);
    args.vpath=strdup("");
    args.vname=strdup(varname);
    total_size = 1;
    for (i=0; i<3; i++) {
        args.offsets[i] = (int) start[i];
        args.counts[i]  = (int) count[i];
        total_size     *= count[i];
    }
    total_size *= vts_res.vartype_size;
    args.max_read=total_size;

    if (DEBUG>3) fprintf(stderr, "-- %s, rank %d: get data: varname=%s offsets=(%d,%d,%d) counts=(%d,%d,%d) total_size=%ld\n",
        __func__, ds->rank, varname,
        args.offsets[0], args.offsets[1], args.offsets[2],
        args.counts[0], args.counts[1], args.counts[2],
        total_size);

    Func_Timer("ADIOS_READ_READ_VAR_OP",
            rc = nssi_call_rpc_sync(&svcs[ds->svc_index],
                    ADIOS_READ_READ_VAR_OP,
                    &args,
                    data,
                    total_size,
                    &res););
    free(args.vpath);
    free(args.vname);
    if (rc != NSSI_OK) {
        fprintf(stderr, "NSSI ERROR: ADIOS_READ_READ_VAR_OP failed\n");
        return(-1);
    }

    return(res.bytes_read);
}

int64_t adios_read_nssi_read_var_byid (ADIOS_GROUP    * gp,
                             int              varid,
                             const uint64_t  * start,
                             const uint64_t  * count,
                             void           * data)
{
    adios_error (err_invalid_varid, "NSSI does not know variable indicies, only variable names can be used.");
    return -err_invalid_varid;
}


int64_t adios_read_nssi_read_local_var (ADIOS_GROUP * gp, const char * varname,
                                      int vidx, const uint64_t * start,
                                      const uint64_t * count, void * data)
{  
    adios_error (err_operation_not_supported, "adios_read_local_var() is not supported with NSSI method.");
    return adios_errno;
}
