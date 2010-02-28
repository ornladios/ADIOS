/**  @file main.c
 *
 *   @brief Driver for the LWFS name server.
 *
 *   @author Ron Oldfield (raoldfi\@sandia.gov).
 *   $Revision: 1264 $.
 *   $Date: 2007-02-27 15:30:26 -0700 (Tue, 27 Feb 2007) $.
 */


#ifndef __STDC_CONSTANT_MACROS
#define __STDC_CONSTANT_MACROS
#endif

#include "config.h"

#include <unistd.h>
#include <string.h>
#include <stdlib.h>
#include <sys/mman.h>

#include "adios.h"

#ifdef HAVE_NSSI
#ifdef HAVE_PORTALS
#include "nssi_ptls.h"
#endif
#ifdef HAVE_INFINIBAND
#include "nssi_ib.h"
#endif
#include "nssi_server.h"
#include "nssi_logger.h"

#include "adios_nssi_args.h"
#include "adios_nssi_config.h"
#endif

#include "io_timer.h"

#include <mpi.h>
#include <algorithm>
#include <map>

using namespace std;



#ifdef __LIBCATAMOUNT__
#define ntohs(value) 0
#endif


/* Need a struct to encapsulate open file info.
 */
struct open_file {
    char fname[ADIOS_PATH_MAX];
//    open_file(const struct lwfs_ib_connection *c) {
//        qp_num=c->msg_qp_num;
//    }
    open_file(const char *name) {
        strcpy(fname, name);
    }
};
/* Need a comparison operator to pass into the open_file_map
 */
struct open_file_lt
{
    bool operator()(const struct open_file &of1, const struct open_file &of2) const
    {
//        log_debug(rpc_debug_level, "cqp1.qp_num == %u", cqp1.qp_num);
//        log_debug(rpc_debug_level, "cqp2.qp_num == %u", cqp2.qp_num);

        if (strcmp(of1.fname, of2.fname) <0) return TRUE;

        return FALSE;
    }
};

/* Map of open files */
static map<struct open_file, int64_t, open_file_lt> open_file_map;
typedef map<struct open_file, int64_t, open_file_lt>::iterator open_file_map_iterator_t;
typedef pair<struct open_file, int64_t> open_file_map_pair_t;
static pthread_mutex_t open_file_map_mutex=PTHREAD_RECURSIVE_MUTEX_INITIALIZER_NP;
static pthread_cond_t  open_file_map_cond =PTHREAD_COND_INITIALIZER;

/* Need a struct to encapsulate open file info.
 */
struct offset {
    int64_t orank;
    char opath[ADIOS_PATH_MAX];
    char oname[ADIOS_PATH_MAX];
//    offset(const struct lwfs_ib_connection *c) {
//        qp_num=c->msg_qp_num;
//    }
    offset(const int64_t rank, const char *path, const char *name) {
        orank=rank;
        strcpy(opath, path);
        strcpy(oname, name);
    }
};
/* Need a comparison operator to pass into the offset_map
 */
struct offset_lt
{
    bool operator()(const struct offset &o1, const struct offset &o2) const
    {
        if (o1.orank < o2.orank) return TRUE;
        if ((o1.orank == o2.orank) && (strcmp(o1.opath, o2.opath) < 0)) return TRUE;
        if ((o1.orank == o2.orank) && (strcmp(o1.opath, o2.opath) == 0) && (strcmp(o1.oname, o2.oname) < 0)) return TRUE;

        return FALSE;
    }
};

/* Map of open files */
static map<struct offset, void *, offset_lt> offset_map;
typedef map<struct offset, void *, offset_lt>::iterator offset_map_iterator_t;
typedef pair<struct offset, void *> offset_map_pair_t;
static pthread_mutex_t offset_map_mutex=PTHREAD_RECURSIVE_MUTEX_INITIALIZER_NP;
static pthread_cond_t  offset_map_cond =PTHREAD_COND_INITIALIZER;

/* -------------------- PRIVATE FUNCTIONS ---------- */
void open_file_add(char *fname, int64_t fd)
{
    open_file_map[fname]=fd;
}

int64_t open_file_get(char *fname)
{
    open_file of(fname);
    int64_t fd=-1;

    open_file_map_iterator_t iter=open_file_map.find(of);
    if (iter != open_file_map.end()) {
        fd=iter->second;
    }

    return(fd);
}

//void open_file_del(lwfs_ib_connection *conn)
//{
//    conn_map_iterator_t iter;
//    conn_qp cqp(conn->msg_qp_num);
//
//    log_debug(rpc_debug_level, "begin");
//    log_debug(rpc_debug_level, "deleting connection with qp_num=%d", conn->msg_qp_num);
//    conn_map.erase(cqp);
//    log_debug(rpc_debug_level, "end");
//}
void open_file_del(char *fname)
{
    open_file_map_iterator_t iter;
    open_file of(fname);

    open_file_map.erase(of);
}


void offset_add(int64_t rank, char *path, char *name, void *data)
{
    offset o(rank, path, name);
    printf("adding offset rank(%ld) opath(%s) oname(%s)\n", rank, path, name);
    offset_map[o]=data;
}

void *offset_get(int64_t rank, char *path, char *name)
{
    offset o(rank, path, name);
    void *data=NULL;

    printf("looking for offset rank(%ld) opath(%s) oname(%s)\n", rank, path, name);

    offset_map_iterator_t iter=offset_map.find(o);
    if (iter != offset_map.end()) {
        data=iter->second;
    }

    return(data);
}

//void offset_del(lwfs_ib_connection *conn)
//{
//    conn_map_iterator_t iter;
//    conn_qp cqp(conn->msg_qp_num);
//
//    log_debug(rpc_debug_level, "begin");
//    log_debug(rpc_debug_level, "deleting connection with qp_num=%d", conn->msg_qp_num);
//    conn_map.erase(cqp);
//    log_debug(rpc_debug_level, "end");
//}
void offset_del(int64_t rank, char *path, char *name)
{
    offset_map_iterator_t iter;
    offset o(rank, path, name);

    offset_map.erase(o);
}


/**
 * The next 3 utility functions are lifted from IOR.
 */
/******************************************************************************/
/*
 * Extract key/value pair from hint string.
 */

void
ExtractHint(char * settingVal,
            char * valueVal,
            char * hintString)
{
    char * settingPtr,
         * valuePtr,
         * tmpPtr1,
         * tmpPtr2;

    settingPtr = (char *)strtok(hintString, "=");
    valuePtr = (char *)strtok(NULL, " \t\r\n");
    tmpPtr1 = settingPtr;
    tmpPtr2 = (char *)strstr(settingPtr, "MPIIO_HINT__");
    if (tmpPtr1 == tmpPtr2) {
        settingPtr += strlen("MPIIO_HINT__");
    }
    strcpy(settingVal, settingPtr);
    strcpy(valueVal, valuePtr);
} /* ExtractHint() */


/******************************************************************************/
/*
 * Set hints for MPIIO, HDF5, or NCMPI.
 */
#define MAX_HINT_STR 1024
void
SetHints(MPI_Info * mpiHints, char * hintsFileName)
{
    char           hintString[MAX_HINT_STR],
                   settingVal[MAX_HINT_STR],
                   valueVal[MAX_HINT_STR];
    extern char ** environ;
    int            i;
    FILE         * fd;

    /*
     * This routine checks for hints from the environment and/or from the
     * hints files.  The hints are of the form:
     * 'MPIIO_HINT__<hint>=<value>', <hint> is the full name of the hint
     * to be set, and <value> is the hint value.
     * E.g., 'setenv MPIIO_HINT__panfs_concurrent_write 1'
     * or 'MPIIO_HINT__panfs_concurrent_write=1' in the hints file.
     */
    MPI_Info_create(mpiHints);

    /* get hints from environment */
    for (i = 0; environ[i] != NULL; i++) {
        /* if this is an IOR_HINT, pass the hint to the info object */
        if (strncmp(environ[i], "MPIIO_HINT", strlen("MPIIO_HINT")) == 0) {
            strcpy(hintString, environ[i]);
            ExtractHint(settingVal, valueVal, hintString);
            MPI_Info_set(*mpiHints, settingVal, valueVal);
        }
    }

    /* get hints from hints file */
    if (strcmp(hintsFileName, "") != 0) {

        /* open the hint file */
        fd = fopen(hintsFileName, "r");
        if (fd == NULL) {
            printf("cannot open hints file\n");
        } else {
            /* iterate over hints file */
            while(fgets(hintString, MAX_HINT_STR, fd) != NULL) {
                if (strncmp(hintString, "MPIIO_HINT", strlen("MPIIO_HINT")) == 0) {
                    ExtractHint(settingVal, valueVal, hintString);
                    MPI_Info_set(*mpiHints, settingVal, valueVal);
                }
            }
            /* close the hints files */
            if (fclose(fd) != 0) printf("cannot close hints file\n");
        }
    }
} /* SetHints() */


/******************************************************************************/
/*
 * Show all hints (key/value pairs) in an MPI_Info object.
 */

void ShowHints(MPI_Info * mpiHints)
{
    char key[MPI_MAX_INFO_VAL],
         value[MPI_MAX_INFO_VAL];
    int  flag,
         i,
         nkeys;

    MPI_Info_get_nkeys(*mpiHints, &nkeys);

    for (i = 0; i < nkeys; i++) {
        MPI_Info_get_nthkey(*mpiHints, i, key);
        MPI_Info_get(*mpiHints, key, MPI_MAX_INFO_VAL-1, value, &flag);
        printf("mpiHint[%d]: %s = %s\n", i, key, value);
    }
} /* ShowHints() */

/* -------------------- SERVER-SIDE STUBS ---------- */


/**
 * @brief Open a netcdf dataset.
 *
 * Open an ADIOS dataset.
 */
int nssi_staging_open_stub(
        const unsigned long request_id,
        const nssi_remote_pid *caller,
        const adios_open_args *args,
        const nssi_rma *data_addr,
        const nssi_rma *res_addr)
{
    int rc = 0;
    char omode[2];
    adios_open_res res;  /* this is what we send back to the client */
    MPI_Info mpiHints = MPI_INFO_NULL;

    int64_t fd;

    memset(&res, 0, sizeof(res));

    printf("calling nssi_staging_open_stub(%s, %d)\n", args->fname, args->mode);

//    SetHints(&mpiHints, "");
//    ShowHints(&mpiHints);

    fd = open_file_get(args->fname);
    if (fd == -1) {
        omode[0]='\0';
        omode[1]='\0';
        switch(args->mode) {
        case ADIOS_MODE_READ:
            omode[0]='r';
            break;
        case ADIOS_MODE_WRITE:
            omode[0]='w';
            break;
        case ADIOS_MODE_APPEND:
            omode[0]='a';
            break;
        case ADIOS_MODE_UPDATE:
            omode[0]='u';
            break;
        default:
            break;
        }

        double callTime;
        Start_Timer(callTime);
        printf("start adios_open\n");
        if ((rc = adios_open(&fd, args->gname, args->fname, omode, NULL)) != 0) {
            printf("Error opening file \"%s\": %d\n", args->fname, rc);
            goto cleanup;
        }
        printf("end adios_open\n");
        Stop_Timer("adios_open", callTime);

        open_file_add(args->fname, fd);
    }

    res.fd=fd;

cleanup:
    /* send the ncid and return code back to client */
    rc = nssi_send_result(caller, request_id, rc, &res, res_addr);

    return rc;
}

int nssi_staging_group_size_stub(
        const unsigned long request_id,
        const nssi_remote_pid *caller,
        const adios_group_size_args *args,
        const nssi_rma *data_addr,
        const nssi_rma *res_addr)
{
    int rc = 0;
    uint64_t total_size=0;

    printf("calling nssi_staging_group_size_stub(%ld)\n", args->fd);

    double callTime;
    Start_Timer(callTime);
    rc = adios_group_size(args->fd, args->data_size, &total_size);
    if (rc != 0) {
        printf("adios_group_size failed: %d\n", rc);
    }
    Stop_Timer("adios_group_size", callTime);

    /* send result to client */
    rc = nssi_send_result(caller, request_id, rc, NULL, res_addr);

    return rc;
}

int nssi_staging_close_stub(
        const unsigned long request_id,
        const nssi_remote_pid *caller,
        const adios_close_args *args,
        const nssi_rma *data_addr,
        const nssi_rma *res_addr)
{
    int rc = 0;

    printf("calling nssi_staging_close_stub(%ld)\n", args->fd);

    double callTime;
    Start_Timer(callTime);
    adios_close(args->fd);
    Stop_Timer("adios_close", callTime);

    /* send result to client */
    rc = nssi_send_result(caller, request_id, rc, NULL, res_addr);

    return rc;
}

int nssi_staging_read_stub(
        const unsigned long request_id,
        const nssi_remote_pid *caller,
        const adios_read_args *args,
        const nssi_rma *data_addr,
        const nssi_rma *res_addr)
{
    int rc = 0;
    adios_read_res res;
    char vpathname[ADIOS_PATH_MAX];
    int pathlen;
    char *v=NULL;

    printf("calling nssi_staging_read_stub(%ld)\n", args->fd);

//    pathlen=strlen(args->vpath);
//    if (args->vpath[pathlen-1]=='/') {
//        sprintf(vpathname, "%s%s", args->vpath, args->vname);
//    } else {
//        sprintf(vpathname, "%s/%s", args->vpath, args->vname);
//    }

    v=(char *)calloc(args->max_read, 1);

    double callTime;
    Start_Timer(callTime);
    adios_read(args->fd, args->vname, v, args->max_read);
    Stop_Timer("adios_read", callTime);

    res.bytes_read=args->max_read;

    rc = nssi_put_data(caller, v, res.bytes_read, data_addr, -1);
    if (rc != NSSI_OK) {
        printf("Could not put var data on client\n");
        goto cleanup;
    }

cleanup:
    /* send result to client */
    rc = nssi_send_result(caller, request_id, rc, &res, res_addr);

    return rc;
}

int nssi_staging_write_stub(
        const unsigned long request_id,
        const nssi_remote_pid *caller,
        const adios_write_args *args,
        const nssi_rma *data_addr,
        const nssi_rma *res_addr)
{
    int rc = 0;
    adios_write_res res;
    char vpathname[ADIOS_PATH_MAX];
    int pathlen;
    char *v=NULL;
    int i=0;

    printf("calling nssi_staging_write_stub(fd=%ld, vsize=%ld)\n", args->fd, args->vsize);

//    pathlen=strlen(args->vpath);
//    if (args->vpath[pathlen-1]=='/') {
//        sprintf(vpathname, "%s%s", args->vpath, args->vname);
//    } else {
//        sprintf(vpathname, "%s/%s", args->vpath, args->vname);
//    }

    v=(char *)malloc(args->vsize);

    rc = nssi_get_data(caller, v, args->vsize, data_addr);
    if (rc != NSSI_OK) {
        printf("Could not get var data on client\n");
        goto cleanup;
    }

    printf("vname(%s) vsize(%ld) is_offset(%d) is_scalar(%d) rank(%ld)\n", args->vname, args->vsize, args->is_offset, args->is_scalar, args->writer_rank);

    if (args->is_offset) {
        // put offset in list
        void *old_data=offset_get(args->writer_rank, "" /*args->vpath*/, args->vname);
        if (old_data != NULL) {
            offset_del(args->writer_rank, "" /*args->vpath*/, args->vname);
            free(old_data);
        }
        offset_add(args->writer_rank, "" /*args->vpath*/, args->vname, v);
    }

    if (!args->is_scalar) {
        // write all offsets for clients rank to update adios internals
        for(i=0;i<args->offsets.offsets_len;i++) {
            void *odata=offset_get(args->writer_rank, args->offsets.offsets_val[i].vpath, args->offsets.offsets_val[i].vname);
            if (odata != NULL) {
                printf("updating oname(%s)\n", args->offsets.offsets_val[i].vname);
                adios_write(args->fd, args->offsets.offsets_val[i].vname, odata);
            }
        }
    }

    double callTime;
    Start_Timer(callTime);
    adios_write(args->fd, args->vname, v);
    Stop_Timer("adios_write", callTime);

    res.bytes_written=args->vsize;


cleanup:
    /* send result to client */
    rc = nssi_send_result(caller, request_id, rc, &res, res_addr);

    return rc;
}

int nssi_staging_start_calc_stub(
        const unsigned long request_id,
        const nssi_remote_pid *caller,
        const adios_start_calc_args *args,
        const nssi_rma *data_addr,
        const nssi_rma *res_addr)
{
    int rc = 0;

    printf("calling nssi_staging_start_calc_stub(%ld)\n", args->fd);


    /* send result to client */
    rc = nssi_send_result(caller, request_id, rc, NULL, res_addr);

    return rc;
}

int nssi_staging_stop_calc_stub(
        const unsigned long request_id,
        const nssi_remote_pid *caller,
        const adios_stop_calc_args *args,
        const nssi_rma *data_addr,
        const nssi_rma *res_addr)
{
    int rc = 0;

    printf("calling nssi_staging_stop_calc_stub(%ld)\n", args->fd);


    /* send result to client */
    rc = nssi_send_result(caller, request_id, rc, NULL, res_addr);

    return rc;
}

int nssi_staging_end_iter_stub(
        const unsigned long request_id,
        const nssi_remote_pid *caller,
        const adios_end_iter_args *args,
        const nssi_rma *data_addr,
        const nssi_rma *res_addr)
{
    int rc = 0;

    printf("calling nssi_staging_end_iter_stub(%ld)\n", args->fd);


    /* send result to client */
    rc = nssi_send_result(caller, request_id, rc, NULL, res_addr);

    return rc;
}

/* -------- END SERVER-SIDE STUBS -------------- */

int nssi_staging_server_init(const char *adios_config_file)
{
    int rc=NSSI_OK;

    printf("start adios_init(%s)\n", adios_config_file);
    rc = adios_init(adios_config_file);
    if (rc != 1) {
        printf("adios_init() failed: %d\n", rc);
        return(-1);
    }
    printf("end adios_init(%s)\n", adios_config_file);

    /* register server stubs */
    NSSI_REGISTER_SERVER_STUB(ADIOS_OPEN_OP, nssi_staging_open_stub, adios_open_args, adios_open_res);
    NSSI_REGISTER_SERVER_STUB(ADIOS_GROUP_SIZE_OP, nssi_staging_group_size_stub, adios_group_size_args, void);
    NSSI_REGISTER_SERVER_STUB(ADIOS_READ_OP, nssi_staging_read_stub, adios_read_args, adios_read_res);
    NSSI_REGISTER_SERVER_STUB(ADIOS_WRITE_OP, nssi_staging_write_stub, adios_write_args, adios_write_res);
    NSSI_REGISTER_SERVER_STUB(ADIOS_END_ITER_OP, nssi_staging_end_iter_stub, adios_end_iter_args, void);
    NSSI_REGISTER_SERVER_STUB(ADIOS_START_CALC_OP, nssi_staging_start_calc_stub, adios_start_calc_args, void);
    NSSI_REGISTER_SERVER_STUB(ADIOS_STOP_CALC_OP, nssi_staging_stop_calc_stub, adios_stop_calc_args, void);
    NSSI_REGISTER_SERVER_STUB(ADIOS_CLOSE_OP, nssi_staging_close_stub, adios_close_args, void);

    return 0;
}



static void generate_contact_info(nssi_remote_pid *myid)
{
    nssi_remote_pid *all_pids=NULL;
    int rank, np;
    char contact_path[1024];

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    printf("rank (%d)\n", rank);

    if (rank==0) {
        MPI_Comm_size(MPI_COMM_WORLD, &np);
        all_pids=(nssi_remote_pid *)malloc(np*sizeof(nssi_remote_pid));
    }
    MPI_Gather(myid, sizeof(nssi_remote_pid), MPI_BYTE,
               all_pids, sizeof(nssi_remote_pid), MPI_BYTE,
               0, MPI_COMM_WORLD);
    if (rank==0) {
        char *contact_file=getenv("ADIOS_NSSI_CONTACT_INFO");
        if (contact_file==NULL) {
            printf("ADIOS_NSSI_CONTACT_INFO env var is undefined.\n");
            free(all_pids);
            return;
        }
        sprintf(contact_path, "%s.%04d", contact_file, rank);
        printf("creating contact file (%s)\n", contact_path);
        FILE *f=fopen(contact_path, "w");
        for (int i=0;i<np;i++) {
            fprintf(f, "%u@%u@%s@%u\n",
                    all_pids[i].nid, all_pids[i].pid,
                    all_pids[i].hostname, (unsigned int)ntohs(all_pids[i].port));
        }
//        fprintf(f, "%u@%u@%s@%u\n",
//                myid->nid, myid->pid,
//                myid->hostname, (unsigned int)ntohs(myid->port));
        fclose(f);
        free(all_pids);
    }
}


/**
 * @brief The LWFS xfer-server.
 */
int main(int argc, char **argv)
{
    int rc = NSSI_OK;

    nssi_service nssi_svc;
//    log_level debug_level;
//    char logfile[1024];
    int rank, np;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &np);

    /* options that can be overriden by the command-line */
    bool daemon_flag = false;
    int verbose = 3;  /* default debug_level */
    int num_threads = 0;
    int server_pid = 128;   /* process ID of the server */
    int server_port = 7728; /* TCP port of the server */

    memset(&nssi_svc, 0, sizeof(nssi_service));

//    /* initialize and enable logging */
//    if (args_info.logfile_arg != NULL) {
//        sprintf(logfile, "%s.%04d", args_info.logfile_arg, rank);
//        logger_init((log_level)args_info.verbose_arg, logfile);
//    } else {
//        logger_init((log_level)args_info.verbose_arg, NULL);
//    }
//    netcdf_debug_level=(log_level)(log_level)args_info.verbose_arg;
//    debug_level = (log_level)args_info.verbose_arg;

    logger_init((log_level)verbose, NULL);

    if (daemon_flag) {
        nssi_daemonize();
    }

#ifdef HAVE_PORTALS
    nssi_ptl_init(PTL_IFACE_SERVER, server_pid);
    rc = nssi_rpc_init(NSSI_RPC_PTL, NSSI_RPC_XDR);
    if (rc != NSSI_OK) {
        printf("could not init rpc: %s\n",
                nssi_err_str(rc));
        return rc;
    }
    nssi_remote_pid myid;
    memset(&myid, 0, sizeof(nssi_remote_pid));
    nssi_ptl_get_id(&myid);
    generate_contact_info(&myid);
#endif
#ifdef HAVE_INFINIBAND
    memset(&nssi_svc.req_addr.match_id, 0, sizeof(nssi_remote_pid));
    strcpy(nssi_svc.req_addr.match_id.hostname, args_info.server_addr_arg);
    nssi_svc.req_addr.match_id.port = args_info.server_port_arg;

    nssi_ib_init(&nssi_svc);
    rc = nssi_rpc_init(NSSI_RPC_IB, NSSI_RPC_XDR);
    if (rc != NSSI_OK) {
        printf("could not init rpc: %s\n",
                nssi_err_str(rc));
        return rc;
    }
    generate_contact_info(&nssi_svc.req_addr.match_id);
#endif

    printf("Initialize staging service\n");

    /* initialize the lwfs service */
    rc = nssi_service_init(0, NSSI_SHORT_REQUEST_SIZE, &nssi_svc);
    if (rc != NSSI_OK) {
        printf("could not init nssi_svc: %s\n",
                nssi_err_str(rc));
        return -1;
    }

    /* initialize staging service */
    rc = nssi_staging_server_init(argv[1]);

    /* start processing requests */
    nssi_svc.max_reqs = -1;
    rc = nssi_service_start(&nssi_svc, num_threads);
    if (rc != NSSI_OK) {
        printf("exited nssi_svc: %s\n",
                nssi_err_str(rc));
    }

    /* shutdown the nssi_svc */
    printf("shutting down service library\n");
    nssi_service_fini(&nssi_svc);

    nssi_rpc_fini();

    MPI_Finalize();

    return rc;
}
