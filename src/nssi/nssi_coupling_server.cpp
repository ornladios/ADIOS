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
#include <list>

using namespace std;



#ifdef __LIBCATAMOUNT__
#define ntohs(value) 0
#endif


static int global_rank=-1;
static int DEBUG=0;


struct var_details_t {
    char var_path[ADIOS_PATH_MAX];
    char var_name[ADIOS_PATH_MAX];
    int  ndims;
    int8_t is_scalar;

    uint64_t  cache_offset; /* offset into the timestep cache_buffer where the var data can be found */
    char     *cache_ptr;    /* address of var data (== cache_buffer+cache_offset) */
    uint64_t  len;          /* var data length in bytes */
    uint64_t  num_elements; /* number of datatype elements in buf (len/atype_size) */

    enum ADIOS_DATATYPES atype; /* adios type of data in buf*/
    int32_t              atype_size;

    char    **offset_path;
    char    **offset_name;
    uint64_t *offset;     /* starting corner (eg. 0,0,0 is the origin of a cube) */
    char    **count_path;
    char    **count_name;
    uint64_t *count;      /* num elements in each dimension (eg. 3,3,3 is a cube of size 3) */
    char    **global_path;
    char    **global_name;
    uint64_t *global;      /* num elements in each dimension (eg. 3,3,3 is a cube of size 3) */

    var_details_t() {
        var_path[0]='\0';
        var_name[0]='\0';
        is_scalar=TRUE;
        ndims=0;
        cache_offset=0;
        cache_ptr=NULL;
        len=0;
        num_elements=0;
        atype=adios_unknown;
        atype_size=0;
        offset_path=NULL;
        offset_name=NULL;
        offset=NULL;
        count_path=NULL;
        count_name=NULL;
        count=NULL;
        global_path=NULL;
        global_name=NULL;
        global=NULL;
    }
};
struct var_details_lt
{
    bool operator()(const char* vn1, const char* vn2) const
    {
        if (DEBUG>2) printf("var_details_lt(vn1=%s ; vn2=%s)\n",
                vn1, vn2);

        if (strcmp(vn1, vn2) <0) return TRUE;

        return FALSE;
    }
};

/* Map of variable details */
typedef map<char*, var_details_t*, var_details_lt> var_details_map_t;
typedef map<char*, var_details_t*, var_details_lt>::iterator var_details_map_iterator_t;
typedef pair<char*, var_details_t*> var_details_map_pair_t;
static pthread_mutex_t var_details_map_mutex=PTHREAD_RECURSIVE_MUTEX_INITIALIZER_NP;
static pthread_cond_t  var_details_map_cond =PTHREAD_COND_INITIALIZER;

typedef list<char*> var_details_list_t;
typedef list<char*>::iterator var_details_list_iterator_t;


struct file_timestep_t {
    int32_t  timestep;
    int8_t   is_complete;
    int8_t   is_locked;
    uint64_t group_size;
    char    *cache_buffer;
    int64_t  cache_size;
    int64_t  bytes_left;

    var_details_map_t vars_map;
    var_details_list_t vars_list;

    file_timestep_t(int32_t ts) {
        timestep=ts;
        is_complete=FALSE;
        is_locked=FALSE;
        group_size=0;
        cache_buffer=NULL;
        cache_size=0;
    }
};

/* list of timesteps */
typedef list<file_timestep_t*> file_timestep_list_t;
typedef list<file_timestep_t*>::iterator file_timestep_list_iterator_t;
static pthread_mutex_t file_timestep_list_mutex=PTHREAD_RECURSIVE_MUTEX_INITIALIZER_NP;
static pthread_cond_t  file_timestep_list_cond =PTHREAD_COND_INITIALIZER;


/* Need a struct to encapsulate open file info.
 */
struct open_file {
    char    writer_id[ADIOS_PATH_MAX];
    char    reader_id[ADIOS_PATH_MAX];
    char    gname[ADIOS_PATH_MAX];
    char    ofname[ADIOS_PATH_MAX];
    char    omode[2];
    int64_t ofdesc;
    int8_t  is_open;
    int64_t current_timestep;
    int64_t last_written_timestep;
    file_timestep_list_t timesteps;

    open_file(const char *name) {
        strcpy(ofname, name);
        ofdesc=-1;
        current_timestep=0;
        last_written_timestep=0;
        is_open=FALSE;
    }
    open_file(const char *name, const int64_t desc) {
        strcpy(ofname, name);
        ofdesc=desc;
        current_timestep=0;
        last_written_timestep=0;
        is_open=FALSE;
    }
};
struct open_file_lt
{
    bool operator()(const char* ofn1, const char* ofn2) const
    {
        if (DEBUG>2) printf("open_file_lt(ofn1=%s ; ofn2=%s)\n",
                ofn1, ofn2);

        if (strcmp(ofn1, ofn2) <0) return TRUE;

        return FALSE;
    }
};

/* Map of open files */
static map<char*, open_file*, open_file_lt> open_file_map;
typedef map<char*, open_file*, open_file_lt>::iterator open_file_map_iterator_t;
typedef pair<char*, open_file*> open_file_map_pair_t;


static map<int64_t, open_file*> open_filedesc_map;
typedef map<int64_t, open_file*>::iterator open_filedesc_map_iterator_t;
typedef pair<int64_t, open_file*> open_filedesc_map_pair_t;


static pthread_mutex_t open_file_map_mutex=PTHREAD_RECURSIVE_MUTEX_INITIALIZER_NP;
static pthread_cond_t  open_file_map_cond =PTHREAD_COND_INITIALIZER;



/* -------------------- PRIVATE FUNCTIONS ---------- */

void open_file_printall(char *prefix)
{
    open_file *of=NULL;
    open_file_map_iterator_t iter=open_file_map.begin();
    while (iter != open_file_map.end()) {
        char *key=iter->first;
        of=iter->second;
        if (DEBUG>5) printf("%s: myrank(%d): open_file(key=%s ; of=%p ; of->ofname=%s ; of->ofdesc=%ld)\n",
                prefix, global_rank, key, of, of->ofname, of->ofdesc);

        iter++;
    }
}

open_file *open_file_add(char *fname)
{
    open_file *of=NULL;

    if (DEBUG>2) printf("open_file_add(fname=%s)\n",
            fname);

    char *key=strdup(fname);

    of = new open_file(fname);
    of->ofdesc = (int64_t)of;
    open_file_map[key]=of;
    open_filedesc_map[of->ofdesc]=of;

    open_file_printall("open_file_add");

    return(of);
}
open_file *open_file_get(char *fname)
{
    open_file *of=NULL;

    open_file_printall("open_file_get(by fname)");

    open_file_map_iterator_t iter=open_file_map.find(fname);
    if (iter != open_file_map.end()) {
        of=iter->second;
    }

    return(of);
}
open_file *open_file_get(int64_t fd)
{
    open_file *of=NULL;

    open_file_printall("open_file_get(by fd)");

    open_filedesc_map_iterator_t iter=open_filedesc_map.find(fd);
    if (iter != open_filedesc_map.end()) {
        of=iter->second;
    }

    return(of);
}
void open_file_del(char *fname)
{
    open_file *of=open_file_get(fname);

    if (of!=NULL) {
        open_filedesc_map.erase(of->ofdesc);
    }

    open_file_map.erase(fname);
}

file_timestep_t *timestep_get_current(open_file *of)
{
    file_timestep_t *ts=NULL;
    file_timestep_list_iterator_t iter=of->timesteps.begin();
    for(;iter != of->timesteps.end();iter++) {
        if ((*iter)->timestep == of->current_timestep) {
            ts=*iter;
            break;
        }
    }
    return(ts);
}
file_timestep_t *timestep_get(open_file *of, uint64_t ts_num)
{
    file_timestep_t *ts=NULL;
    file_timestep_list_iterator_t iter=of->timesteps.begin();
    for(;iter != of->timesteps.end();iter++) {
        if ((*iter)->timestep == ts_num) {
            ts=*iter;
            break;
        }
    }
    return(ts);
}

void var_details_add(struct open_file *of, var_details_t *vd)
{
    file_timestep_t *ts=NULL;
    ts = timestep_get_current(of);

    char *key=strdup(vd->var_name);

    ts->vars_map[key]=vd;
    ts->vars_list.push_back(key);
}
var_details_t *var_details_get(struct open_file *of, char *vname)
{
    file_timestep_t *ts=NULL;
    var_details_t *vd=NULL;

    ts = timestep_get_current(of);
    var_details_map_iterator_t iter=ts->vars_map.find(vname);
    if (iter != ts->vars_map.end()) {
        vd=iter->second;
    }

    return(vd);
}
var_details_t *var_details_get(struct open_file *of, uint64_t ts_num, char *vname)
{
    file_timestep_t *ts=NULL;
    var_details_t *vd=NULL;

    ts = timestep_get(of, ts_num);
    var_details_map_iterator_t iter=ts->vars_map.find(vname);
    if (iter != ts->vars_map.end()) {
        vd=iter->second;
    }

    return(vd);
}
void var_details_del(struct open_file *of, char *vname)
{
    file_timestep_t *ts=NULL;
    ts = timestep_get_current(of);
    ts->vars_map.erase(vname);
}



int grank, gsize;

MPI_Comm comm_self=MPI_COMM_SELF;
MPI_Comm comm_world=MPI_COMM_WORLD;


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


int32_t getTypeSize(
        enum ADIOS_DATATYPES type,
        void *val)
{
    switch (type)
    {
    case adios_byte:
    case adios_unsigned_byte:
        return 1;

    case adios_string:
        return strlen ((char *) val);

    case adios_short:
    case adios_unsigned_short:
        return 2;

    case adios_integer:
    case adios_unsigned_integer:
        return 4;

    case adios_real:
        return 4;

    case adios_long:
    case adios_unsigned_long:
        return 8;

    case adios_double:
        return 8;

    case adios_long_double:
        return 16;

    case adios_complex:
        return 2 * 4;

    case adios_double_complex:
        return 2 * 8;

    default:
        return -1;
    }
}


#ifdef __cplusplus
extern "C" {
#endif

extern int adios_nssi_filter_is_anon_dim(
        int fd,
        const char *dimname);
extern void adios_nssi_filter_set_anon_dim(
        int fd,
        const char *dimname,
        const uint64_t dimvalue);

#ifdef __cplusplus
}
#endif
int write_cache(const char *client_id, int8_t free_memory)
{
    int rc=0;

    int timestamps_written=0;

    if (DEBUG>2) printf("myrank(%d): enter write_cache(client_id=%s)\n", grank, client_id);


    open_file_map_iterator_t of_iter=open_file_map.begin();
    for(;of_iter != open_file_map.end();of_iter++) {
        open_file *of=of_iter->second;

        if ((client_id != NULL) && (strcmp(client_id, of->writer_id)) != 0) {
            fprintf(stdout, "flusher != writer - file(%s) client_id(%s) writer_id(%s)\n", of->ofname, client_id, of->writer_id);
            continue;
        }

        int64_t fd;
        if (DEBUG>3) printf("start adios_open\n");
        if (DEBUG>3) printf("adios_open: using MPI_COMM_WORLD\n");
        Func_Timer("adios_open", rc = adios_open(&fd, of->gname, of->ofname, of->omode, comm_world););
        if (rc != 0) {
            printf("Error opening file \"%s\": %d\n", of->ofname, rc);
        }
        if (DEBUG>3) printf("end adios_open\n");

        file_timestep_list_iterator_t ts_iter=of->timesteps.begin();
        while(ts_iter != of->timesteps.end()) {
            var_details_t *vd=NULL;

            file_timestep_t *ts=*ts_iter;

            fprintf(stdout, "timestep - file(%s) timestep(%ld) is_locked(%d) is_complete(%d)\n",
                    of->ofname, ts->timestep, ts->is_locked, ts->is_complete);
            if (ts->is_locked == TRUE) {
                /* a reader has this timestep locked */
                fprintf(stdout, "timestep is locked - file(%s) timestep(%ld) reader_id(%s)\n", of->ofname, ts->timestep, of->reader_id);
//                break;
            }
            if (ts->is_complete == FALSE) {
                /* the writer has not finished this timestep */
                fprintf(stdout, "timestep is NOT complete - file(%s) timestep(%ld)\n", of->ofname, ts->timestep);
                break;
            }

            uint64_t total_size=0;
            Func_Timer("adios_group_size", rc = adios_group_size(fd, ts->group_size, &total_size););
            if (rc != 0) {
                printf("adios_group_size failed: %d\n", rc);
            }

            var_details_list_iterator_t vd_list_iter=ts->vars_list.begin();
            for(;vd_list_iter != ts->vars_list.end();vd_list_iter++) {
                var_details_map_iterator_t vd_iter=ts->vars_map.find(*vd_list_iter);
                if (vd_iter == ts->vars_map.end()) {
                    fprintf(stdout, "couldn't find vname(%s) in the vd map\n", *vd_list_iter);
                    continue;
                }

                vd = vd_iter->second;

                if (!vd->is_scalar) {
                    for(int i=0;i<vd->ndims;i++) {
                        vd->offset[i]=0;
                        if (DEBUG>3) printf("writing offset myrank(%d) vd(%d) vpath(%s) vname(%s) opath(%s) oname(%s) odata(%lu)\n",
                                grank, i, vd->var_path, vd->var_name, vd->offset_path[i], vd->offset_name[i], vd->offset[i]);
                        if (adios_nssi_filter_is_anon_dim(fd, vd->offset_name[i]) == TRUE) {
                            if (DEBUG>2) printf("server_rank(%d) writing anon dim vname(%s)\n", global_rank, vd->offset_name[i]);
                            adios_nssi_filter_set_anon_dim(fd, vd->offset_name[i], vd->offset[i]);
                        } else {
                            if (DEBUG>2) printf("server_rank(%d) writing aggregated offset vname(%s)\n", global_rank, vd->offset_name[i]);
                            Func_Timer("adios_set_path_var", adios_set_path_var(fd, vd->offset_path[i], vd->offset_name[i]););
                            Func_Timer("adios_write", adios_write(fd, vd->offset_name[i], &(vd->offset[i])););
                        }
                    }
                    for(int i=0;i<vd->ndims;i++) {
                        if (DEBUG>3) printf("writing count myrank(%d) vd(%d) vpath(%s) vname(%s) dpath(%s) dname(%s) ddata(%lu)\n",
                                grank, i, vd->var_path, vd->var_name, vd->count_path[i], vd->count_name[i], vd->global[i]);
                        if (adios_nssi_filter_is_anon_dim(fd, vd->count_name[i]) == TRUE) {
                            if (DEBUG>2) printf("server_rank(%d) writing anon dim vname(%s)\n", global_rank, vd->count_name[i]);
                            adios_nssi_filter_set_anon_dim(fd, vd->count_name[i], vd->global[i]);
                        } else {
                            if (DEBUG>2) printf("server_rank(%d) writing aggregated dim vname(%s)\n", global_rank, vd->global_name[i]);
                            Func_Timer("adios_set_path_var", adios_set_path_var(fd, vd->count_path[i], vd->count_name[i]););
                            Func_Timer("adios_write", adios_write(fd, vd->count_name[i], &(vd->global[i])););
                        }
                    }
                }

                if (DEBUG>3) printf("writing array myrank(%d) fd(%ld) vname(%s) vdata(%lu)\n", grank, fd, vd->var_name, *(uint64_t*)vd->cache_ptr);
                if (DEBUG>2) printf("server_rank(%d) writing aggregated array vname(%s)\n", global_rank, vd->var_name);
                Func_Timer("adios_set_path_var", adios_set_path_var(fd, vd->var_path, vd->var_name););
                Func_Timer("adios_write", adios_write(fd, vd->var_name, vd->cache_ptr););
            }

            if ((free_memory==TRUE) && (ts->is_locked == FALSE)) {
                /* release resources associated with the timestep */
                var_details_map_iterator_t vd_iter=ts->vars_map.begin();
                while(vd_iter != ts->vars_map.end()) {
                    char *key=vd_iter->first;

                    vd = vd_iter->second;
                    for(int i=0;i<vd->ndims;i++) {
                        if (vd->global_path != NULL) free(vd->global_path[i]);
                        if (vd->global_name != NULL) free(vd->global_name[i]);
                        if (vd->count_path != NULL) free(vd->count_path[i]);
                        if (vd->count_name != NULL) free(vd->count_name[i]);
                        if (vd->offset_path != NULL) free(vd->offset_path[i]);
                        if (vd->offset_name != NULL) free(vd->offset_name[i]);
                    }
                    if (vd->global_path != NULL) free(vd->global_path);
                    if (vd->global_name!= NULL) free(vd->global_name);
                    if (vd->global != NULL) free(vd->global);
                    if (vd->count_path != NULL) free(vd->count_path);
                    if (vd->count_name != NULL) free(vd->count_name);
                    if (vd->count != NULL) free(vd->count);
                    if (vd->offset_path != NULL) free(vd->offset_path);
                    if (vd->offset_name != NULL) free(vd->offset_name);
                    if (vd->offset != NULL) free(vd->offset);

                    ts->vars_map.erase(vd_iter++);

                    free(key);
                }
                if (ts->cache_buffer != NULL) free(ts->cache_buffer);

                ts->vars_list.clear();

                delete ts;
            }

            timestamps_written++;
            of->last_written_timestep=ts->timestep;

            if ((free_memory==TRUE) && (ts->is_locked == FALSE)) {
                of->timesteps.erase(ts_iter++);
            } else {
                ts_iter++;
            }
        }

        Func_Timer("adios_close", adios_close(fd););
    }

    if (DEBUG>2) printf("myrank(%d): exit write_cache(client_id=%s)\n", grank, client_id);

    return(timestamps_written);
}

/* -------------------- SERVER-SIDE STUBS ---------- */


/**
 * @brief Open a netcdf dataset.
 *
 * Open an ADIOS dataset.
 */
int nssi_coupling_open_stub(
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
    int64_t fd=-1;
    open_file *of=NULL;

    memset(&res, 0, sizeof(res));

    if (DEBUG>2) printf("myrank(%d): enter nssi_coupling_open_stub(%s, %d)\n", grank, args->fname, args->mode);

//    SetHints(&mpiHints, "");
//    ShowHints(&mpiHints);

    of = open_file_get(args->fname);
    if (DEBUG>3) printf("myrank(%d): nssi_coupling_open_stub(%s, %d)\n", grank, args->fname, args->mode);
    if (of == NULL) {
        of=open_file_add(args->fname);

        strcpy(of->writer_id, args->client_id);
        strcpy(of->gname, args->gname);
    }
    if (of->is_open == FALSE) {
        of->omode[0]='\0';
        of->omode[1]='\0';
        switch(args->mode) {
        case ADIOS_MODE_READ:
            of->omode[0]='r';
            break;
        case ADIOS_MODE_WRITE:
            of->omode[0]='w';
            break;
        case ADIOS_MODE_APPEND:
            of->omode[0]='a';
            break;
        case ADIOS_MODE_UPDATE:
            of->omode[0]='u';
            break;
        default:
            break;
        }
        of->current_timestep++;
        file_timestep_t *ts=new file_timestep_t(of->current_timestep);
        of->timesteps.push_back(ts);

        of->is_open=TRUE;
    }

    if (DEBUG>2) printf("myrank(%d): exit nssi_coupling_open_stub(%s, %d): fd=%ld, fd=%p\n", grank, args->fname, args->mode, of->ofdesc, of->ofdesc);

    res.fd=of->ofdesc;

cleanup:
    /* send the ncid and return code back to client */
    rc = nssi_send_result(caller, request_id, rc, &res, res_addr);

    return rc;
}

int flush_unlocked_timesteps(void)
{
    return(write_cache(NULL, TRUE));
}

int nssi_coupling_group_size_stub(
        const unsigned long request_id,
        const nssi_remote_pid *caller,
        const adios_group_size_args *args,
        const nssi_rma *data_addr,
        const nssi_rma *res_addr)
{
    int rc = 0;
    uint64_t total_size=0;
    open_file *of=NULL;

    if (DEBUG>2) printf("myrank(%d): enter nssi_coupling_group_size_stub(fd=%ld, pg_size=%ld)\n", grank, args->fd, args->data_size);

    of=open_file_get(args->fd);
    if (of == NULL) {
        printf("nssi_coupling_group_size_stub failed: file not open(fd=%ld)\n", args->fd);
        rc=-1;
    } else {
        file_timestep_t *current_ts=timestep_get_current(of);
        if (current_ts == NULL) {
            printf("nssi_coupling_group_size_stub failed: could not get current timestep(fd=%ld)\n", args->fd);
            rc=-1;
        } else {
            if (DEBUG>2) printf("myrank(%d): creating timestep cache (of=%p ; ts=%p ; timestep=%d)\n",
                    grank, of, current_ts, current_ts->timestep);

            current_ts->group_size=args->data_size;
            current_ts->cache_size=args->data_size;
            current_ts->bytes_left=args->data_size;
            current_ts->cache_buffer=(char *)malloc(args->data_size);
            while(current_ts->cache_buffer == NULL) {
                // insufficient memory to start a new timestep.  check if there are unlocked timesteps to flush.
                if (flush_unlocked_timesteps() == FALSE) {
                    // there are no timesteps to flush.  fail.
                    rc=-1;
                    break;
                }
                // retry
                current_ts->cache_buffer=(char *)malloc(args->data_size);
            }
            if (DEBUG>2) printf("myrank(%d): timestep cache (of=%p ; ts=%p ; timestep=%d ; cache_size=%ld ; buffer=%p)\n",
                    grank, of, current_ts, current_ts->timestep, current_ts->cache_size, current_ts->cache_buffer);
        }
    }

    /* send result to client */
    rc = nssi_send_result(caller, request_id, rc, NULL, res_addr);

    if (DEBUG>2) printf("myrank(%d): exit nssi_coupling_group_size_stub(%ld)\n", grank, args->fd);

    return rc;
}

int nssi_coupling_close_stub(
        const unsigned long request_id,
        const nssi_remote_pid *caller,
        const adios_close_args *args,
        const nssi_rma *data_addr,
        const nssi_rma *res_addr)
{
    int rc = 0;
    open_file *of=NULL;

    if (DEBUG>2) printf("myrank(%d): enter nssi_coupling_close_stub(%ld, %s)\n", grank, args->fd, args->fname);

    of=open_file_get(args->fd);
    if (of == NULL) {
        printf("nssi_coupling_close_stub failed: file not open(fd=%ld)\n", args->fd);
        rc=-1;
    } else {
        file_timestep_t *current_ts=timestep_get_current(of);
        if (current_ts == NULL) {
            printf("nssi_coupling_close_stub failed: could not get current timestep(fd=%ld)\n", args->fd);
            rc=-1;
        } else {
            current_ts->is_complete=TRUE;
        }

        of->is_open=FALSE;
    }

    /* send result to client */
    rc = nssi_send_result(caller, request_id, rc, NULL, res_addr);

    if (DEBUG>2) printf("myrank(%d): exit nssi_coupling_close_stub(%ld, %s)\n", grank, args->fd, args->fname);

    return rc;
}

int nssi_coupling_read_stub(
        const unsigned long request_id,
        const nssi_remote_pid *caller,
        const adios_read_args *args,
        const nssi_rma *data_addr,
        const nssi_rma *res_addr)
{
    int rc = 0;

    if (DEBUG>2) printf("myrank(%d): enter nssi_coupling_read_stub(%ld)\n", grank, args->fd);

    // ENOTSUP - use read API

    if (DEBUG>2) printf("myrank(%d): exit nssi_coupling_read_stub(%ld)\n", grank, args->fd);

    return rc;
}

static void recursive_copy_chunk_in(char *src,
                                    var_details_t *vd,
                                    uint64_t src_offset,
                                    uint64_t dst_offset,
                                    uint64_t *src_index,
                                    uint64_t *dst_index,
                                    int current_dim)
{
    uint64_t my_src_offset=0;
    uint64_t my_dst_offset=0;

    if (current_dim < vd->ndims-1) {
        for (int i=0;i<vd->count[current_dim];i++) {
            my_src_offset = src_index[current_dim];
            my_dst_offset = dst_index[current_dim];
//            if (DEBUG > 3) printf("src_offset[%d](%lu)\n", current_dim, vd->offset[current_dim]);
            my_dst_offset += (vd->offset[current_dim] * vd->atype_size);
            for (int j=current_dim+1;j<vd->ndims;j++) {
                my_src_offset *= vd->count[j];
                my_dst_offset *= vd->global[j];
            }

            src_index[current_dim+1]=0;
            dst_index[current_dim+1]=0;
            recursive_copy_chunk_in(src, vd, src_offset+my_src_offset, dst_offset+my_dst_offset,
                                    src_index, dst_index, current_dim+1);
            src_index[current_dim] += vd->atype_size;
            dst_index[current_dim] += vd->atype_size;
        }
    } else {
        dst_offset += (vd->offset[current_dim] * vd->atype_size);
        memcpy(vd->cache_ptr + dst_offset,
               src + src_offset,
               vd->count[current_dim]*vd->atype_size);
    }
}

static void copy_chunk_in(char *src,
                          var_details_t *vd)
{
    uint64_t *src_index=(uint64_t *)calloc(vd->ndims, sizeof(uint64_t));
    uint64_t *dst_index=(uint64_t*)calloc(vd->ndims, sizeof(uint64_t));
    uint64_t src_offset=0;
    uint64_t dst_offset=0;
    int current_dim=0;

    memset(src_index, 0, vd->ndims*sizeof(uint64_t));
    memset(dst_index, 0, vd->ndims*sizeof(uint64_t));
    recursive_copy_chunk_in(src, vd, src_offset, dst_offset, src_index, dst_index, current_dim);

    free(src_index);
    free(dst_index);
}

static void recursive_copy_chunk_out(char *dst,
                                     var_details_t *vd,
                                     uint64_t dst_offset,
                                     uint64_t src_offset,
                                     uint64_t *dst_index,
                                     uint64_t *src_index,
                                     int current_dim)
{
    uint64_t my_dst_offset=0;
    uint64_t my_src_offset=0;

    if (current_dim < vd->ndims-1) {
        for (int i=0;i<vd->count[current_dim];i++) {
            my_dst_offset = dst_index[current_dim];
            my_src_offset = src_index[current_dim];
//            if (DEBUG > 3) printf("dst_offset[%d](%lu)\n", current_dim, vd->offset[current_dim]);
            my_src_offset += (vd->offset[current_dim] * vd->atype_size);
            for (int j=current_dim+1;j<vd->ndims;j++) {
                my_dst_offset *= vd->count[j];
                my_src_offset *= vd->global[j];
            }

            dst_index[current_dim+1]=0;
            src_index[current_dim+1]=0;
            recursive_copy_chunk_out(dst, vd, dst_offset+my_dst_offset, src_offset+my_src_offset,
                                     dst_index, src_index, current_dim+1);
            dst_index[current_dim] += vd->atype_size;
            src_index[current_dim] += vd->atype_size;
        }
    } else {
        src_offset += (vd->offset[current_dim] * vd->atype_size);
        memcpy(dst + dst_offset,
               vd->cache_ptr + src_offset,
               vd->count[current_dim]*vd->atype_size);
    }
}

static void copy_chunk_out(char *dst,
                           var_details_t *vd)
{
    uint64_t *dst_index=(uint64_t *)calloc(vd->ndims, sizeof(uint64_t));
    uint64_t *src_index=(uint64_t*)calloc(vd->ndims, sizeof(uint64_t));
    uint64_t dst_offset=0;
    uint64_t src_offset=0;
    int current_dim=0;

    memset(dst_index, 0, vd->ndims*sizeof(uint64_t));
    memset(src_index, 0, vd->ndims*sizeof(uint64_t));
    recursive_copy_chunk_out(dst, vd, dst_offset, src_offset, dst_index, src_index, current_dim);

    free(dst_index);
    free(src_index);
}

int nssi_coupling_write_stub(
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

    if (DEBUG>2) printf("myrank(%d): enter nssi_coupling_write_stub(fd=%ld, vsize=%ld)\n", grank, args->fd, args->vsize);

    open_file *of=open_file_get(args->fd);
    file_timestep_t *ts=timestep_get_current(of);

    if (DEBUG>2) printf("myrank(%d): timestep cache (of=%p ; ts=%p ; timestep=%d ; cache_size=%ld ; bytes_left=%ld ; buffer=%p)\n",
            grank, of, ts, ts->timestep, ts->cache_size, ts->bytes_left, ts->cache_buffer);

    v=(char *)malloc(args->vsize);

    Func_Timer("nssi_get_data", rc = nssi_get_data(caller, v, args->vsize, data_addr););
    if (rc != NSSI_OK) {
        printf("Could not get var data on client\n");
        goto cleanup;
    }

    if (DEBUG>3) printf("server_rank(%d) vname(%s) vsize(%ld) is_scalar(%d) writer_rank(%ld)\n",
            global_rank, args->vname, args->vsize, args->is_scalar, args->writer_rank);

    if (!args->is_scalar) {
        if (DEBUG>2) printf("server_rank(%d) caching non-scalar vname(%s) vsize(%ld)\n", global_rank, args->vname, args->vsize);
        if (DEBUG>3) printf("allocated v(%p), len(%ld)\n", v, args->vsize);

        uint64_t write_offset=1;
        uint64_t write_size=1;
        var_details_t *vd=NULL;

        vd = var_details_get(of, args->vname);
        if (vd == NULL) {
            if (DEBUG>3) printf("************\nvar_details NOT found for %s\n************\n", args->vname);

            vd = new var_details_t;

            vd->is_scalar=FALSE;

//            vd->fd = args->fd;
            strcpy(vd->var_path, args->vpath);
            strcpy(vd->var_name, args->vname);
            vd->ndims = args->offsets.offsets_len;

            vd->num_elements = 1;
            vd->global_path = (char **)calloc(args->gdims.gdims_len, sizeof(char *));
            vd->global_name = (char **)calloc(args->gdims.gdims_len, sizeof(char *));
            vd->global  = (uint64_t *)calloc(args->gdims.gdims_len, sizeof(uint64_t));;
            for (int i=0;i<args->gdims.gdims_len;i++) {
                vd->global_path[i] = strdup(args->gdims.gdims_val[i].vpath);
                vd->global_name[i] = strdup(args->gdims.gdims_val[i].vname);
                vd->global[i] = args->gdims.gdims_val[i].vdata;

                vd->num_elements *= args->gdims.gdims_val[i].vdata;
            }

            vd->offset_path = (char **)calloc(args->offsets.offsets_len, sizeof(char *));
            vd->offset_name = (char **)calloc(args->offsets.offsets_len, sizeof(char *));
            vd->offset = (uint64_t *)calloc(args->offsets.offsets_len, sizeof(uint64_t));
            for (int i=0;i<args->offsets.offsets_len;i++) {
                vd->offset_path[i] = strdup(args->offsets.offsets_val[i].vpath);
                vd->offset_name[i] = strdup(args->offsets.offsets_val[i].vname);
                vd->offset[i] = args->offsets.offsets_val[i].vdata;
            }

            vd->count_path = (char **)calloc(args->ldims.ldims_len, sizeof(char *));
            vd->count_name = (char **)calloc(args->ldims.ldims_len, sizeof(char *));
            vd->count  = (uint64_t *)calloc(args->ldims.ldims_len, sizeof(uint64_t));;
            for (int i=0;i<args->ldims.ldims_len;i++) {
                vd->count_path[i] = strdup(args->ldims.ldims_val[i].vpath);
                vd->count_name[i] = strdup(args->ldims.ldims_val[i].vname);
                vd->count[i] = args->ldims.ldims_val[i].vdata;
            }

            vd->atype      = (enum ADIOS_DATATYPES)args->atype;
            vd->atype_size = getTypeSize(vd->atype, v);
            vd->len        = vd->num_elements * vd->atype_size;

            vd->cache_offset = ts->cache_size-ts->bytes_left;
            vd->cache_ptr    = ts->cache_buffer+vd->cache_offset;

            ts->bytes_left  -= vd->len;

            var_details_add(of, vd);

            if (DEBUG > 3) printf("********** ts->cache_size(%lu)-ts->bytes_left(%lu)==vd->cache_offset(%lu)\n",
                     ts->cache_size, ts->bytes_left, vd->cache_offset);
            if (DEBUG > 3) printf("********** vd->num_elements(%lu)*vd->atype_size(%u)==vd->len(%lu)\n",
                     vd->num_elements, vd->atype_size, vd->len);
        }

        if (DEBUG > 3) printf("********** ts->cache_size(%lu) ;  vd->cache_offset(%lu) ; vd->cache_buffer(%p) ; vd->cache_ptr(%p)\n",
                 ts->cache_size, vd->cache_offset, ts->cache_buffer, vd->cache_ptr);

        for (int i=0;i<args->offsets.offsets_len;i++) {
            vd->offset[i] = args->offsets.offsets_val[i].vdata;
            write_offset *= args->offsets.offsets_val[i].vdata;

            if (DEBUG > 3) printf("********** var_name(%s) ;  vd->offset[%d](%lu)\n",
                     vd->var_name, i, vd->offset[i]);
        }
        write_offset *= vd->atype_size;

        for (int i=0;i<args->ldims.ldims_len;i++) {
            vd->count[i] = args->ldims.ldims_val[i].vdata;
            write_size *= args->ldims.ldims_val[i].vdata;

            if (DEBUG > 3) printf("********** var_name(%s) ;  vd->count[%d](%lu)\n",
                     vd->var_name, i, vd->count[i]);
        }
        write_size *= vd->atype_size;

        /* sanity checks */
        if (write_size != args->vsize) {
            printf("write size conflict: write_size(%lu) != args->vsize(%lu)\n",
                    write_size, args->vsize);
        }
        if (write_offset+write_size > vd->len) {
            printf("variable overflow: write_offset(%lu)+write_size(%lu)=(%lu) > vd->len(%lu)\n",
                    write_offset, write_size, write_offset+write_size, vd->len);
        }
        if (vd->cache_offset+write_size > ts->cache_size) {
            printf("cache overflow: vd->cache_offset(%lu)+write_size(%lu)=(%lu) > ts->cache_size(%lu)\n",
                    vd->cache_offset, write_size, vd->cache_offset+write_size, ts->cache_size);
            abort();
        }

        if (DEBUG > 3) printf("********** write_size(%lu) ;  write_offset(%lu)\n",
                 write_size, write_offset);

        copy_chunk_in(v, vd);

    } else {
        if (DEBUG>2) printf("server_rank(%d) writing scalar vname(%s) vsize(%ld)\n", global_rank, args->vname, args->vsize);

        uint64_t write_offset=1;
        uint64_t write_size=1;
        var_details_t *vd=NULL;

        vd = var_details_get(of, args->vname);
        if (vd == NULL) {
            if (DEBUG>3) printf("************\nvar_details NOT found for %s\n************\n", args->vname);

            vd = new var_details_t;

            vd->is_scalar=TRUE;

//            vd->fd = args->fd;
            strcpy(vd->var_path, args->vpath);
            strcpy(vd->var_name, args->vname);

            vd->ndims = 0;
            vd->num_elements = 1;

            vd->atype      = (enum ADIOS_DATATYPES)args->atype;
            vd->atype_size = getTypeSize(vd->atype, v);
            vd->len        = vd->num_elements * vd->atype_size;

            vd->cache_offset = ts->cache_size-ts->bytes_left;
            vd->cache_ptr    = ts->cache_buffer+vd->cache_offset;

            ts->bytes_left  -= vd->len;

            var_details_add(of, vd);
        }

        write_offset = 0;
        write_size   = vd->atype_size;

        /* sanity checks */
        if (write_size != args->vsize) {
            printf("write size conflict: write_size(%lu) != args->vsize(%lu)\n",
                    write_size, args->vsize);
        }
        if (write_offset+write_size > vd->len) {
            printf("variable overflow: write_offset(%lu)+write_size(%lu)=(%lu) > vd->len(%lu)\n",
                    write_offset, write_size, write_offset+write_size, vd->len);
        }
        if (vd->cache_offset+write_size > ts->cache_size) {
            printf("cache overflow: vd->cache_offset(%lu)+write_size(%lu)=(%lu) > ts->cache_size(%lu)\n",
                    vd->cache_offset, write_size, vd->cache_offset+write_size, ts->cache_size);
        }

        memcpy(vd->cache_ptr+write_offset, v, write_size);

    }

    free(v);

    res.bytes_written=args->vsize;

cleanup:

    /* send result to client */
    rc = nssi_send_result(caller, request_id, rc, &res, res_addr);

    if (DEBUG>2) printf("myrank(%d): exit nssi_coupling_write_stub(fd=%ld, vsize=%ld)\n", grank, args->fd, args->vsize);

    return rc;
}

int nssi_coupling_start_calc_stub(
        const unsigned long request_id,
        const nssi_remote_pid *caller,
        const adios_start_calc_args *args,
        const nssi_rma *data_addr,
        const nssi_rma *res_addr)
{
    int rc = 0;

    if (DEBUG>2) printf("myrank(%d): enter nssi_coupling_start_calc_stub(%ld)\n", grank, args->fd);

    /* send result to client */
    rc = nssi_send_result(caller, request_id, rc, NULL, res_addr);

    if (DEBUG>2) printf("myrank(%d): exit nssi_coupling_start_calc_stub(%ld)\n", grank, args->fd);

    return rc;
}

int nssi_coupling_stop_calc_stub(
        const unsigned long request_id,
        const nssi_remote_pid *caller,
        const adios_stop_calc_args *args,
        const nssi_rma *data_addr,
        const nssi_rma *res_addr)
{
    int rc = 0;

    if (DEBUG>2) printf("myrank(%d): enter nssi_coupling_stop_calc_stub(%ld)\n", grank, args->fd);


    /* send result to client */
    rc = nssi_send_result(caller, request_id, rc, NULL, res_addr);

    if (DEBUG>2) printf("myrank(%d): exit nssi_coupling_stop_calc_stub(%ld)\n", grank, args->fd);

    return rc;
}

int nssi_coupling_end_iter_stub(
        const unsigned long request_id,
        const nssi_remote_pid *caller,
        const adios_end_iter_args *args,
        const nssi_rma *data_addr,
        const nssi_rma *res_addr)
{
    int rc = 0;

    if (DEBUG>2) printf("myrank(%d): enter nssi_coupling_end_iter_stub(%ld)\n", grank, args->fd);


    /* send result to client */
    rc = nssi_send_result(caller, request_id, rc, NULL, res_addr);

    if (DEBUG>2) printf("myrank(%d): exit nssi_coupling_end_iter_stub(%ld)\n", grank, args->fd);

    return rc;
}

int nssi_coupling_finalize_stub(
        const unsigned long request_id,
        const nssi_remote_pid *caller,
        const adios_finalize_args *args,
        const nssi_rma *data_addr,
        const nssi_rma *res_addr)
{
    int rc = 0;

    if (DEBUG>2) printf("myrank(%d): enter nssi_coupling_finalize_stub(%s)\n", grank, args->client_id);


//    write_cache(args->client_id);


    /* send result to client */
    rc = nssi_send_result(caller, request_id, rc, NULL, res_addr);

    if (DEBUG>2) printf("myrank(%d): exit nssi_coupling_finalize_stub(%s)\n", grank, args->client_id);

    return rc;
}

/**
 * @brief Open a netcdf dataset.
 *
 * Open an ADIOS dataset.
 */
int nssi_coupling_read_fopen_stub(
        const unsigned long request_id,
        const nssi_remote_pid *caller,
        const adios_read_fopen_args *args,
        const nssi_rma *data_addr,
        const nssi_rma *res_addr)
{
    int rc = 0;
    char omode[2];
    adios_open_res res;  /* this is what we send back to the client */
    MPI_Info mpiHints = MPI_INFO_NULL;
    int64_t fd=-1;
    open_file *of=NULL;

    memset(&res, 0, sizeof(res));

    if (DEBUG>2) printf("myrank(%d): enter nssi_coupling_read_fopen_stub(%s, %d)\n", grank, args->fname, args->requested_timestep);

//    SetHints(&mpiHints, "");
//    ShowHints(&mpiHints);

    of = open_file_get(args->fname);
    if (DEBUG>3) printf("myrank(%d): nssi_coupling_read_fopen_stub(%s, %d)\n", grank, args->fname, args->requested_timestep);
    if (of == NULL) {
        res.fd=-1; /* the requested file is not currently in core.  out of core not supported.  return error. */
        goto cleanup;
    }

    file_timestep_t *ts=timestep_get(of, args->requested_timestep);
    if (ts == NULL) {
        res.fd=-1; /* the requested timestep is not currently in core.  out of core not supported.  return error. */
        goto cleanup;
    } else {
        if ((ts->is_complete==TRUE) && (ts->is_locked==FALSE)) {
            strcpy(of->reader_id, args->client_id);
            ts->is_locked=TRUE;
            if (DEBUG>3) printf("nssi_coupling_read_fopen_stub(%s, %d): timestep is complete and locked\n", grank, args->fname, args->requested_timestep, of->ofdesc, of->ofdesc);
        }
    }

    res.fd=of->ofdesc;

cleanup:
    if (DEBUG>2) printf("myrank(%d): exit nssi_coupling_read_fopen_stub(%s, %d): fd=%ld\n", grank, args->fname, args->requested_timestep, res.fd);

    /* send the ncid and return code back to client */
    rc = nssi_send_result(caller, request_id, rc, &res, res_addr);

    return rc;
}

int nssi_coupling_read_fclose_stub(
        const unsigned long request_id,
        const nssi_remote_pid *caller,
        const adios_read_fclose_args *args,
        const nssi_rma *data_addr,
        const nssi_rma *res_addr)
{
    int rc = 0;
    open_file *of=NULL;

    if (DEBUG>2) printf("myrank(%d): enter nssi_coupling_read_fclose_stub(%ld, %s)\n", grank, args->fd, args->fname);

    of=open_file_get(args->fd);
    if (of == NULL) {
        printf("nssi_coupling_read_fclose_stub failed: file not open(fd=%ld)\n", args->fd);
        rc=-1;
    } else {
        file_timestep_t *ts=timestep_get(of, args->open_timestep);
        if (ts == NULL) {
            printf("nssi_coupling_read_fclose_stub failed: could not find timestep(fd=%ld, ts_num=%lu)\n", args->fd, args->open_timestep);
            rc=-1;
        } else {
            if(strcmp(of->reader_id, args->client_id) == 0) {
                of->reader_id[0]='\0';
                ts->is_locked=FALSE;
            } else {
                printf("nssi_coupling_read_fclose_stub failed: timestep(fd=%ld, ts_num=%lu) fopen by different client (me=%s,opener=%s)\n",
                        args->fd, args->open_timestep, args->client_id, of->reader_id);
                rc=-1;
            }
        }
    }

    /* send result to client */
    rc = nssi_send_result(caller, request_id, rc, NULL, res_addr);

    if (DEBUG>2) printf("myrank(%d): exit nssi_coupling_read_fclose_stub(%ld, %s)\n", grank, args->fd, args->fname);

    return rc;
}

int nssi_coupling_read_get_vartype_size_stub(
        const unsigned long request_id,
        const nssi_remote_pid *caller,
        const adios_read_get_vartype_size_args *args,
        const nssi_rma *data_addr,
        const nssi_rma *res_addr)
{
    int rc = 0;
    open_file *of=NULL;
    adios_read_get_vartype_size_res res;

    if (DEBUG>2) printf("myrank(%d): enter nssi_coupling_read_get_vartype_size_stub(%ld)\n", grank, args->fd);

    of=open_file_get(args->fd);
    if (of == NULL) {
        printf("nssi_coupling_read_get_vartype_size_stub failed: file not open(fd=%ld)\n", args->fd);
        rc=-1;
    } else {
        file_timestep_t *ts=timestep_get(of, args->open_timestep);
        if (ts == NULL) {
            printf("nssi_coupling_read_get_vartype_size_stub failed: could not find timestep(fd=%ld, ts_num=%lu)\n",
                    args->fd, args->open_timestep);
            rc=-1;
        } else {
            if ((strcmp(of->reader_id, args->client_id) == 0) && (ts->is_locked==TRUE)) {
                var_details_t *vd=var_details_get(of, args->open_timestep, args->vname);
                if (vd == NULL) {
                    printf("nssi_coupling_read_get_vartype_size_stub failed: could not find variable(fd=%ld, ts_num=%lu, varname=%s)\n",
                            args->fd, args->open_timestep, args->vname);
                    rc=-1;
                } else {
                    res.vartype_size=vd->atype_size;
                }
            } else {
                printf("nssi_coupling_read_get_vartype_size_stub failed: timestep(fd=%ld, ts_num=%lu) fopen by different client (me=%s,opener=%s)\n",
                        args->fd, args->open_timestep, args->client_id, of->reader_id);
                rc=-1;
            }
        }
    }

    /* send result to client */
    rc = nssi_send_result(caller, request_id, rc, &res, res_addr);

    if (DEBUG>2) printf("myrank(%d): exit nssi_coupling_read_get_vartype_size_stub(fd=%ld res.vartype_size=%lu)\n",
            grank, args->fd, res.vartype_size);

    return rc;
}

int nssi_coupling_read_read_var_stub(
        const unsigned long request_id,
        const nssi_remote_pid *caller,
        const adios_read_read_var_args *args,
        const nssi_rma *data_addr,
        const nssi_rma *res_addr)
{
    int rc = 0;
    open_file *of=NULL;
    adios_read_read_var_res res;

    char *v=NULL;

    res.bytes_read=-1;

    if (DEBUG>2) printf("myrank(%d): enter nssi_coupling_read_read_var_stub(%ld)\n", grank, args->fd);

    of=open_file_get(args->fd);
    if (of == NULL) {
        printf("nssi_coupling_read_read_var_stub failed: file not open(fd=%ld)\n", args->fd);
        rc=-1;
    } else {
        file_timestep_t *ts=timestep_get(of, args->open_timestep);
        if (ts == NULL) {
            printf("nssi_coupling_read_read_var_stub failed: could not find timestep(fd=%ld, ts_num=%lu)\n", args->fd, args->open_timestep);
            rc=-1;
        } else {
            if ((strcmp(of->reader_id, args->client_id) == 0) && (ts->is_locked==TRUE)) {
                uint64_t read_size=1;
                var_details_t *vd=var_details_get(of, args->open_timestep, args->vname);

                for (int i=0;i<vd->ndims;i++) {
                    vd->offset[i] = args->offsets[i];

                    if (DEBUG > 3) printf("********** var_name(%s) ;  vd->offset[%d](%lu)\n",
                             vd->var_name, i, vd->offset[i]);
                }
                for (int i=0;i<vd->ndims;i++) {
                    vd->count[i] = args->counts[i];
                    read_size *= args->counts[i];

                    if (DEBUG > 3) printf("********** var_name(%s) ;  vd->count[%d](%lu) ;  vd->globals[%d](%lu)\n",
                             vd->var_name, i, vd->count[i], i, vd->global[i]);
                }
                int8_t is_scalar=FALSE;
                if (read_size==1) {
                    is_scalar=TRUE;
                }
                read_size *= vd->atype_size;

                /* sanity checks */
                if (read_size != args->max_read) {
                    printf("read size conflict: read_size(%lu) != args->max_read(%lu)\n",
                            read_size, args->max_read);
                }
                if (vd->cache_offset+read_size > ts->cache_size) {
                    printf("cache overflow: vd->cache_offset(%lu)+read_size(%lu)=(%lu) > ts->cache_size(%lu)\n",
                            vd->cache_offset, read_size, vd->cache_offset+read_size, ts->cache_size);
                    abort();
                }

                if (DEBUG > 3) printf("********** read_size(%lu)\n", read_size);

                v=(char*)malloc(read_size);
                if (is_scalar==TRUE) {
                    memcpy(v, vd->cache_ptr, read_size);
                } else {
                    copy_chunk_out(v, vd);
                }

                Func_Timer("nssi_put_data", rc = nssi_put_data(caller, v, read_size, data_addr, -1););
                if (rc != NSSI_OK) {
                    printf("Could not put var data on client\n");
                    goto cleanup;
                }

                res.bytes_read=read_size;

            } else {
                printf("nssi_coupling_read_read_var_stub failed: timestep(fd=%ld, ts_num=%lu) fopen by different client (me=%s,opener=%s)\n",
                        args->fd, args->open_timestep, args->client_id, of->reader_id);
                rc=-1;
            }
        }
    }

cleanup:
    /* send result to client */
    rc = nssi_send_result(caller, request_id, rc, &res, res_addr);

    if (v) free(v);

    if (DEBUG>2) printf("myrank(%d): exit nssi_coupling_read_read_var_stub(%ld)\n", grank, args->fd);

    return rc;
}

/* -------- END SERVER-SIDE STUBS -------------- */

int nssi_coupling_server_init(const char *adios_config_file)
{
    int rc=NSSI_OK;

    if (DEBUG>3) printf("start adios_init(%s)\n", adios_config_file);
    rc = adios_init(adios_config_file);
    if (rc != 1) {
        printf("adios_init() failed: %d\n", rc);
        return(-1);
    }
    if (DEBUG>3) printf("end adios_init(%s)\n", adios_config_file);


    /* register server stubs */
    NSSI_REGISTER_SERVER_STUB(ADIOS_OPEN_OP,       nssi_coupling_open_stub,       adios_open_args,       adios_open_res);
    NSSI_REGISTER_SERVER_STUB(ADIOS_GROUP_SIZE_OP, nssi_coupling_group_size_stub, adios_group_size_args, void);
    NSSI_REGISTER_SERVER_STUB(ADIOS_READ_OP,       nssi_coupling_read_stub,       adios_read_args,       adios_read_res);
    NSSI_REGISTER_SERVER_STUB(ADIOS_WRITE_OP,      nssi_coupling_write_stub,      adios_write_args,      adios_write_res);
    NSSI_REGISTER_SERVER_STUB(ADIOS_END_ITER_OP,   nssi_coupling_end_iter_stub,   adios_end_iter_args,   void);
    NSSI_REGISTER_SERVER_STUB(ADIOS_START_CALC_OP, nssi_coupling_start_calc_stub, adios_start_calc_args, void);
    NSSI_REGISTER_SERVER_STUB(ADIOS_STOP_CALC_OP,  nssi_coupling_stop_calc_stub,  adios_stop_calc_args,  void);
    NSSI_REGISTER_SERVER_STUB(ADIOS_CLOSE_OP,      nssi_coupling_close_stub,      adios_close_args,      void);
    NSSI_REGISTER_SERVER_STUB(ADIOS_FINALIZE_OP,   nssi_coupling_finalize_stub,   adios_finalize_args,   void);

    NSSI_REGISTER_SERVER_STUB(ADIOS_READ_FOPEN_OP,            nssi_coupling_read_fopen_stub,            adios_read_fopen_args,            adios_read_fopen_res);
    NSSI_REGISTER_SERVER_STUB(ADIOS_READ_FCLOSE_OP,           nssi_coupling_read_fclose_stub,           adios_read_fclose_args,           void);
    NSSI_REGISTER_SERVER_STUB(ADIOS_READ_GET_VARTYPE_SIZE_OP, nssi_coupling_read_get_vartype_size_stub, adios_read_get_vartype_size_args, adios_read_get_vartype_size_res);
    NSSI_REGISTER_SERVER_STUB(ADIOS_READ_READ_VAR_OP,         nssi_coupling_read_read_var_stub,         adios_read_read_var_args,         adios_read_read_var_res);

    return 0;
}



static void generate_contact_info(nssi_remote_pid *myid)
{
    nssi_remote_pid *all_pids=NULL;
    int rank, np;
    char contact_path[1024];

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    //printf("rank (%d)\n", rank);

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
//        sprintf(contact_path, "%s.%04d", contact_file, rank);
        sprintf(contact_path, "%s.tmp", contact_file);
        if (DEBUG>3) printf("creating temp contact file (%s)\n", contact_path);
        FILE *f=fopen(contact_path, "w");
        if (f==NULL) {
            perror("fopen");
        }
        for (int i=0;i<np;i++) {
            fprintf(f, "%u@%u@%s@%u\n",
                    all_pids[i].nid, all_pids[i].pid,
                    all_pids[i].hostname, (unsigned int)ntohs(all_pids[i].port));
        }
//        fprintf(f, "%u@%u@%s@%u\n",
//                myid->nid, myid->pid,
//                myid->hostname, (unsigned int)ntohs(myid->port));
        fclose(f);
        if (DEBUG>3) printf("renaming temp contact file (%s) to contact file (%s)\n", contact_path, contact_file);
        rename(contact_path, contact_file);
        free(all_pids);
    }
    MPI_Barrier(MPI_COMM_WORLD);
}


/**
 * @brief The LWFS xfer-server.
 */
int main(int argc, char **argv)
{
    int rc = NSSI_OK;

    nssi_service nssi_svc;
//    log_level debug_level;
    char logfile[1024];
    int rank, np;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &global_rank);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &np);
    grank=rank;
    gsize=np;

    /* options that can be overriden by the command-line */
    bool daemon_flag = false;
    int verbose = 5;  /* default debug_level */
    int num_threads = 0;
    int server_pid = 128;   /* process ID of the server */
    int server_port = 7728; /* TCP port of the server */

    memset(&nssi_svc, 0, sizeof(nssi_service));

    /* initialize and enable logging */
//    if (args_info.logfile_arg != NULL) {
//        sprintf(logfile, "%s.%04d", "nssi_coupling_server.log", rank);
//        logger_init((log_level)verbose, logfile);
//    } else {
//        logger_init((log_level)verbose, NULL);
//    }
//    netcdf_debug_level=(log_level)(log_level)args_info.verbose_arg;
//    debug_level = (log_level)args_info.verbose_arg;

//    logger_init((log_level)verbose, NULL);

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

    if (DEBUG>3) printf("Initialize staging service\n");

    /* initialize the lwfs service */
    rc = nssi_service_init(0, NSSI_SHORT_REQUEST_SIZE, &nssi_svc);
    if (rc != NSSI_OK) {
        printf("could not init nssi_svc: %s\n",
                nssi_err_str(rc));
        return -1;
    }

    /* initialize staging service */
    rc = nssi_coupling_server_init(argv[1]);

    /* start processing requests */
    nssi_svc.max_reqs = -1;
    rc = nssi_service_start(&nssi_svc, num_threads);
    if (rc != NSSI_OK) {
        printf("exited nssi_svc: %s\n",
                nssi_err_str(rc));
    }

    flush_unlocked_timesteps();

    adios_finalize(rank);

    /* shutdown the nssi_svc */
    if (DEBUG>3) printf("shutting down service library\n");
    nssi_service_fini(&nssi_svc);

    nssi_rpc_fini();

    MPI_Finalize();

    return rc;
}
