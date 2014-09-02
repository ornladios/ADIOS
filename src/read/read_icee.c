/*
  read_icee.c       
  Goal: to create evpath io connection layer in conjunction with 
  write/adios_icee.c
*/
// system libraries
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/queue.h>
#include <sys/socket.h>
#include <sys/times.h>
#include <netinet/in.h>
#include <sys/time.h>
#include <sys/uio.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <ctype.h>
#include <pthread.h>
#include <unistd.h>
#include <assert.h>

// evpath libraries
#include <ffs.h>
#include <atl.h>
//#include <gen_thread.h>
#include <evpath.h>

#ifdef __MACH__
#include <mach/clock.h>
#include <mach/mach.h>
#endif

// local libraries
#include "config.h"
#include "public/adios.h"
#include "public/adios_types.h"
#include "public/adios_read_v2.h"
#include "core/adios_read_hooks.h"
#include "core/adios_logger.h"
#include "public/adios_error.h"

// conditional libraries
#ifdef DMALLOC
#include "dmalloc.h"
#endif

#define MALLOC0(type_t, var) var = malloc(sizeof(type_t)); memset(var, '\0', sizeof(type_t));
#define MYALLOC(var, nbyte)       \
    if (nbyte) {                  \
        var = malloc(nbyte);      \
        assert((var != NULL) && "malloc failure");      \
        memset(var, '\0', nbyte); \
    } else                        \
        var = NULL;

#define MYFREE(var) free(var); var = NULL;

#define MYMIN(a,b)               \
    ({ __typeof__ (a) _a = (a);  \
        __typeof__ (b) _b = (b); \
        _a < _b ? _a : _b; })

///////////////////////////
// Global Variables
///////////////////////////
#include "adios_icee.h"

#define DUMP(fmt, ...) fprintf(stderr, ">>> "fmt"\n", ## __VA_ARGS__); 

typedef struct icee_llist
{
    void*             item;
    struct icee_llist* next;
} icee_llist_t, *icee_llist_ptr_t;

int (*icee_llist_fp)(const void*, const void*);

icee_llist_ptr_t
icee_llist_create(void* item)
{
    icee_llist_ptr_t node = malloc(sizeof(icee_llist_t));
    node->item = item;
    node->next = NULL;
    
    return node;
}

icee_llist_ptr_t
icee_llist_append(const icee_llist_ptr_t root, void* item)
{
    assert(root != NULL);
        
    icee_llist_ptr_t p = root;
    
    while (p->next != NULL)
        p = p->next;

    icee_llist_ptr_t node = malloc(sizeof(icee_llist_t));
    node->item = item;
    node->next = NULL;
    
    p->next = node;

    return node;
}

icee_llist_ptr_t
icee_llist_search(const icee_llist_ptr_t root, int (*fp)(const void*, const void*), const void* arg)
{
    icee_llist_ptr_t p = root;
    
    while (p != NULL)
    {
        if ((*fp)(p->item, arg)) break;
        p = p->next;
    }

    return p;
}

int icee_fileinfo_checkname(const void* item, const void* fname)
{
    icee_fileinfo_rec_ptr_t fp = (icee_fileinfo_rec_ptr_t) item;
    if (strcmp(fp->fname, (char *)fname) == 0)
        return 1;
    else
        return 0;
}

icee_llist_ptr_t icee_filelist = NULL;

void icee_fileinfo_append(const icee_fileinfo_rec_ptr_t root, icee_fileinfo_rec_ptr_t fp)
{
    assert(root != NULL);
    
    icee_fileinfo_rec_ptr_t p = root;

    while (p->next != NULL)
        p = p->next;

    p->next = fp;
}

void icee_fileinfo_copy(icee_fileinfo_rec_ptr_t dest, const icee_fileinfo_rec_ptr_t src)
{
    memcpy(dest, src, sizeof(icee_fileinfo_rec_t));

    dest->fname = strdup(src->fname);
    dest->next = NULL;
}

void icee_fileinfo_free(icee_fileinfo_rec_ptr_t fp)
{
    icee_varinfo_rec_ptr_t vp = fp->varinfo;

    while (vp != NULL)
    {
        free(vp->varname);
        free(vp->gdims);
        free(vp->ldims);
        free(vp->offsets);
        free(vp->data);

        icee_varinfo_rec_ptr_t prev = vp;
        vp = vp->next;

        free(prev);
    }

    free(fp);
    fp = NULL;    
}

icee_varinfo_rec_ptr_t
icee_varinfo_search_byname(icee_varinfo_rec_ptr_t head, const char* name)
{
    icee_varinfo_rec_ptr_t vp = head;

    while (vp != NULL)
    {
        if (strcmp(vp->varname, name) == 0)
            break;
        vp = vp->next;
    }

    return vp;
}

void icee_varinfo_copy(icee_varinfo_rec_ptr_t dest, const icee_varinfo_rec_ptr_t src)
{
    memcpy(dest, src, sizeof(icee_varinfo_rec_t));

    dest->varname = strdup(src->varname);

    uint64_t dimsize = src->ndims * sizeof(uint64_t);
    dest->gdims = malloc(dimsize);
    dest->ldims = malloc(dimsize);
    dest->offsets = malloc(dimsize);

    memcpy(dest->gdims, src->gdims, dimsize);
    memcpy(dest->ldims, src->ldims, dimsize);
    memcpy(dest->offsets, src->offsets, dimsize);
    
    dest->data = malloc(src->varlen);
    memcpy(dest->data, src->data, src->varlen);

    dest->next = NULL;
}

double dsum(const int len, const double* data)
{
    double s = 0.0;
    int i;
    for (i=0; i<len; i++)
        s += data[i];

    return s;
}

void icee_data_print(const int type, const uint64_t varlen, const char* data)
{
    fprintf(stderr, "%10s : %p\n", "*data", data);
    switch (type)
    {
    case 2: // int
        if (data)
        {
            fprintf(stderr, "%10s : %d,%d,%d,...\n", "data", 
                    ((int*)data)[0], ((int*)data)[1], ((int*)data)[2]);
            fprintf(stderr, "%10s : %g\n", "sum", dsum(varlen/8, (double*)data));
        }
        break;
    case 6: // double
        if (data)
        {
            fprintf(stderr, "%10s : %g,%g,%g,...\n", "data", 
                    ((double*)data)[0], ((double*)data)[1], ((double*)data)[2]);
            fprintf(stderr, "%10s : %g\n", "sum", dsum(varlen/8, (double*)data));
        }
        break;
    }
}

void icee_dims_print(const char* name, const int ndims, const uint64_t *dims)
{
    
    switch (ndims)
    {
    case 0:
        fprintf(stderr, "%10s : none\n", name);
        break;
    case 1:
        fprintf(stderr, "%10s : %llu\n", name, dims[0]);
        break;
    case 2:
        fprintf(stderr, "%10s : %llu,%llu\n", 
                name, dims[0], dims[1]);
        break;
    case 3:
        fprintf(stderr, "%10s : %llu,%llu,%llu\n", 
                name, dims[0], dims[1], dims[2]);
        break;
    default:
        fprintf(stderr, "%10s : %llu,%llu,%llu,...\n", 
                name, dims[0], dims[1], dims[2]);
        break;
    }
}

void icee_varinfo_print(const icee_varinfo_rec_ptr_t vp)
{
    int i;

    fprintf(stderr, "===== varinfo (%p) =====\n", vp);

    if (vp)
    {
        fprintf(stderr, "%10s : %s\n", "varname", vp->varname);
        fprintf(stderr, "%10s : %d\n", "varid", vp->varid);
        fprintf(stderr, "%10s : %d\n", "type", vp->type);
        fprintf(stderr, "%10s : %d\n", "typesize", vp->typesize);
        fprintf(stderr, "%10s : %d\n", "ndims", vp->ndims);
        icee_dims_print("gdims", vp->ndims, vp->gdims);
        icee_dims_print("ldims", vp->ndims, vp->ldims);
        icee_dims_print("offsets", vp->ndims, vp->offsets);
        fprintf(stderr, "%10s : %llu\n", "varlen", vp->varlen);
        icee_data_print(vp->type, vp->varlen, vp->data);
    }
    else
    {
        fprintf(stderr, "varinfo is invalid\n");
    }
}

static int
icee_fileinfo_handler(CManager cm, void *vevent, void *client_data, attr_list attrs)
{
    icee_fileinfo_rec_ptr_t event = vevent;
    log_debug("%s (%s)\n", __FUNCTION__, event->fname);

    icee_fileinfo_rec_ptr_t lfp = malloc(sizeof(icee_fileinfo_rec_t));
    icee_fileinfo_copy(lfp, event);

    icee_varinfo_rec_ptr_t eventvp = event->varinfo;
    icee_varinfo_rec_ptr_t *lvpp = &lfp->varinfo;

    while (eventvp != NULL)
    {
        *lvpp = malloc(sizeof(icee_varinfo_rec_t));
        icee_varinfo_copy(*lvpp, eventvp);
        if (adios_verbose_level > 3) DUMP("id,name = %d,%s", (*lvpp)->varid, (*lvpp)->varname);
        
        lvpp = &(*lvpp)->next;
        eventvp = eventvp->next;
    }

    if (icee_filelist == NULL)
        icee_filelist = icee_llist_create((void *)lfp);
    else
    {
        icee_llist_ptr_t head;
        head = icee_llist_search(icee_filelist, icee_fileinfo_checkname, (const void*) lfp->fname);

        if (head == NULL)
            icee_llist_append(icee_filelist, (void*) lfp);
        else
            icee_fileinfo_append((icee_fileinfo_rec_ptr_t)head->item, lfp);
    }

    return 1;
}

static int adios_read_icee_initialized = 0;

CManager cm;
EVstone stone;
EVstone remote_stone;
EVstone stone_r;

/********** Core ADIOS Read functions. **********/

/*
 * Gathers basic MPI information; sets up reader CM.
 */
int
adios_read_icee_init_method (MPI_Comm comm, PairStruct* params)
{   
    log_debug ("%s\n", __FUNCTION__);

    int remote_port = 59999;
    char *remote_host = "localhost";
    int cm_port = 59997;
    char *cm_host = "localhost";
    char *cm_attr = NULL;
    attr_list contact_list;

    PairStruct * p = params;

    while (p)
    {
        if (!strcasecmp (p->name, "cm_attr"))
        {
            cm_attr = p->value;
        }
        else if (!strcasecmp (p->name, "cm_host"))
        {
            cm_host = p->value;
        }
        else if (!strcasecmp (p->name, "cm_port"))
        {
            cm_port = atoi(p->value);
        }
        else if (!strcasecmp (p->name, "remote_host"))
        {
            remote_host = p->value;
        }
        else if (!strcasecmp (p->name, "remote_port"))
        {
            remote_port = atoi(p->value);
        }
        p = p->next;
    }

    if (cm_attr)
    {
        contact_list = attr_list_from_string(cm_attr);
    }
    else
    {
        contact_list = create_attr_list();
        add_int_attr(contact_list, attr_atom_from_string("IP_PORT"), cm_port);
        add_string_attr(contact_list, attr_atom_from_string("IP_HOST"), cm_host);
    }

    //log_info ("cm_attr : %s\n", cm_attr);
    log_info ("cm_host : %s\n", cm_host);
    log_info ("cm_port : %d\n", cm_port);

    if (!adios_read_icee_initialized)
    {
        cm = CManager_create();

        // Listen first
        {
            //cm = CManager_create();
            attr_list contact_list_r;
            contact_list_r = create_attr_list();
            add_int_attr(contact_list_r, attr_atom_from_string("IP_PORT"), cm_port);
            if (CMlisten_specific(cm, contact_list_r) == 0) 
            {
                fprintf(stderr, "error: unable to initialize connection manager.\n");
                exit(-1);
            }

            log_debug("Contact list \"%s\"\n", attr_list_to_string(contact_list_r));

            stone_r = EValloc_stone(cm);
            EVassoc_terminal_action(cm, stone_r, icee_fileinfo_format_list, icee_fileinfo_handler, NULL);

            if (!CMfork_comm_thread(cm)) 
            {
                printf("Fork of communication thread failed, exiting\n");
                exit(-1);
            }
        }


        // Send client info
        stone = EValloc_stone(cm);

        contact_list = create_attr_list();
        add_int_attr(contact_list, attr_atom_from_string("IP_PORT"), remote_port);
        add_string_attr(contact_list, attr_atom_from_string("IP_HOST"), remote_host);

        EVaction evaction = EVassoc_bridge_action(cm, stone, contact_list, remote_stone);
        if (evaction == -1)
        {
            fprintf(stderr, "No connection. Exit.\n");
            exit(1);
        }

        EVsource source = EVcreate_submit_handle(cm, stone, icee_clientinfo_format_list);
        
        icee_clientinfo_rec_t info;
        info.client_host = cm_host;
        info.client_port = cm_port;
        
        EVsubmit(source, &info, NULL);
        //CManager_close(cm);
        

        /*
        cm = CManager_create();
        contact_list = create_attr_list();
        add_int_attr(contact_list, attr_atom_from_string("IP_PORT"), cm_port);

        if (CMlisten_specific(cm, contact_list) == 0) 
        {
            fprintf(stderr, "error: unable to initialize connection manager.\n");
            exit(-1);
        }

        log_debug("Contact list \"%s\"\n", attr_list_to_string(contact_list));

        stone = EValloc_stone(cm);
        EVassoc_terminal_action(cm, stone, icee_fileinfo_format_list, icee_fileinfo_handler, NULL);

        //CMrun_network(cm);
        //CMlisten(cm);
        if (!CMfork_comm_thread(cm)) 
        {
            printf("Fork of communication thread failed, exiting\n");
            exit(-1);
        }
        */
        
        adios_read_icee_initialized = 1;
    }

    return 0;
}

ADIOS_FILE*
adios_read_icee_open_file(const char * fname, MPI_Comm comm)
{
    adios_error (err_operation_not_supported,
                 "ICEE staging method does not support file mode for reading. "
                 "Use adios_read_open() to open a staged dataset.\n");
    return NULL;
}

ADIOS_FILE*
adios_read_icee_open(const char * fname,
                     MPI_Comm comm,
                     enum ADIOS_LOCKMODE lock_mode,
                     float timeout_sec)
{
    log_debug("%s\n", __FUNCTION__);

    icee_llist_ptr_t head = NULL;

    float waited_sec = 0.0;    
    while (waited_sec <= timeout_sec)
    {
        head = icee_llist_search(icee_filelist, icee_fileinfo_checkname, 
                                 (const void*) fname);
        if (head != NULL)
            break;

        usleep(0.1*1E6);
        waited_sec += 0.1;
    }
    
    if (head == NULL)
    {
        log_error ("No data yet\n");
        return NULL;
    }
    
    icee_fileinfo_rec_ptr_t fp = (icee_fileinfo_rec_ptr_t)head->item;
    
    ADIOS_FILE *adiosfile = malloc(sizeof(ADIOS_FILE));

    adiosfile->fh = (uint64_t)fp;

    adiosfile->nvars = fp->nvars;
	adiosfile->var_namelist = malloc(fp->nvars * sizeof(char *));

    icee_varinfo_rec_ptr_t vp = fp->varinfo;

    int i;
    for (i = 0; i < fp->nvars; i++)
    {
        adiosfile->var_namelist[i] = strdup(vp->varname);
        vp = vp->next;
    }

    adiosfile->nattrs = 0;
    adiosfile->attr_namelist = NULL;

    adiosfile->current_step = 0;
    adiosfile->last_step = 1;

    adiosfile->path = strdup(fname);

    adiosfile->version = -1;
    adiosfile->file_size = 0;
    adios_errno = err_no_error;        

    return adiosfile;
}

int 
adios_read_icee_finalize_method ()
{
    log_debug("%s\n", __FUNCTION__);

    if (adios_read_icee_initialized)
    {
        CManager_close(cm);
        adios_read_icee_initialized = 0;
    }

    return 0;
}

void 
adios_read_icee_release_step(ADIOS_FILE *adiosfile) 
{
    log_debug("%s\n", __FUNCTION__);
}

int 
adios_read_icee_advance_step(ADIOS_FILE *adiosfile, int last, float timeout_sec) 
{
    log_debug("%s\n", __FUNCTION__);
    adios_errno = 0;

    icee_fileinfo_rec_ptr_t fp = (icee_fileinfo_rec_ptr_t) adiosfile->fh;
    icee_fileinfo_rec_ptr_t next = NULL;

    float waited_sec = 0.0;
    while (waited_sec <= timeout_sec)
    {
        next = fp->next;

        if (next != NULL)
            break;

        usleep(0.1*1E6);
        waited_sec += 0.1;
    }
    
    if (next != NULL)
    {
        icee_llist_ptr_t head = NULL;
        head = icee_llist_search(icee_filelist, icee_fileinfo_checkname, 
                                 (const void*) fp->fname);
        assert(head != NULL);

        head->item = next;
        adiosfile->fh = (uint64_t) next;

        icee_fileinfo_free(fp);
    }
    else
        adios_error (err_step_notready, 
                     "No more data yet\n");

    return adios_errno;
}

int 
adios_read_icee_close(ADIOS_FILE * fp)
{
    log_debug("%s\n", __FUNCTION__);

    return 0;
}

ADIOS_FILE *
adios_read_icee_fopen(const char *fname, MPI_Comm comm) 
{
    log_error("No support yet: %s\n", __FUNCTION__);
    return NULL;
}

int 
adios_read_icee_is_var_timed(const ADIOS_FILE* fp, int varid) 
{  
    log_error("No support yet: %s\n", __FUNCTION__);
    return 0; 
}

void 
adios_read_icee_get_groupinfo(const ADIOS_FILE *fp, int *ngroups, char ***group_namelist, uint32_t **nvars_per_group, uint32_t **nattrs_per_group) 
{
    log_debug("%s\n", __FUNCTION__);

    icee_fileinfo_rec_ptr_t p = (icee_fileinfo_rec_ptr_t) fp->fh;

    if (p)
    {
        *ngroups = 1;
        *group_namelist = (char **) malloc (sizeof (char*));
        *group_namelist[0] = strdup ("noname");
    }
    
    return;
}

int 
adios_read_icee_check_reads(const ADIOS_FILE* fp, ADIOS_VARCHUNK** chunk) 
{ 
    log_error("No support yet: %s\n", __FUNCTION__);
    return 0; 
}

int adios_read_icee_perform_reads(const ADIOS_FILE *adiosfile, int blocking)
{
    log_debug("%s\n", __FUNCTION__);
    return 0;
}

int
adios_read_icee_inq_var_blockinfo(const ADIOS_FILE* fp,
				      ADIOS_VARINFO* varinfo)
{
    log_error("No support yet: %s\n", __FUNCTION__);
    return 0; 
}

int
adios_read_icee_inq_var_stat(const ADIOS_FILE* fp,
				 ADIOS_VARINFO* varinfo,
				 int per_step_stat,
				 int per_block_stat)
{
    log_error("No support yet: %s\n", __FUNCTION__);
    return 0; 
}


int 
adios_read_icee_schedule_read_byid(const ADIOS_FILE *adiosfile,
				       const ADIOS_SELECTION *sel,
				       int varid,
				       int from_steps,
				       int nsteps,
				       void *data)
{   
    icee_fileinfo_rec_ptr_t fp = (icee_fileinfo_rec_ptr_t) adiosfile->fh;
    log_debug("%s (%d:%s)\n", __FUNCTION__, varid, fp->fname);
    assert(varid < fp->nvars);

    if(nsteps != 1){
        adios_error (err_invalid_timestep,
                     "Only one step can be read from a stream at a time. "
                     "You requested % steps in adios_schedule_read()\n", 
                     nsteps);
        return err_invalid_timestep;
    }
    
    icee_varinfo_rec_ptr_t vp = NULL;
    vp = icee_varinfo_search_byname(fp->varinfo, adiosfile->var_namelist[varid]);
    if (adios_verbose_level > 3) icee_varinfo_print(vp);

    if(!vp){
        adios_error(err_invalid_varid,
                    "Invalid variable id: %d\n",
                    varid);
        return adios_errno;
    }

    if (sel==0)
        memcpy(data, vp->data, vp->varlen);
    else
        switch(sel->type)
        {
        case ADIOS_SELECTION_WRITEBLOCK:
        {
            //DUMP("fp->rank: %d", fp->rank);
            //DUMP("u.block.index: %d", sel->u.block.index);

            if (fp->rank != sel->u.block.index)
                adios_error(err_unspecified,
                            "Block id missmatch. "
                            "Not yet supported by ICEE\n");

            // Possible memory overwrite
            memcpy(data, vp->data, vp->varlen);
            break;
        }
        case ADIOS_SELECTION_BOUNDINGBOX:
        {
            if (vp->ndims != sel->u.bb.ndim)
                adios_error(err_invalid_dimension,
                            "Dimension mismatch\n");

            int i;
            for (i = 0; i < vp->ndims; i++)
            {
                //DUMP("g,l,o = %llu,%llu,%llu", vp->gdims[i], vp->ldims[i], vp->offsets[i]);
                //DUMP("start,count = %llu,%llu", sel->u.bb.start[i], sel->u.bb.count[i]);
                if (sel->u.bb.start[i] != vp->offsets[i])
                    adios_error(err_expected_read_size_mismatch,
                                "Requested data is out of the local size. "
                                "Not yet supported by ICEE\n");

                if (sel->u.bb.start[i] + sel->u.bb.count[i] != vp->offsets[i] + vp->ldims[i])
                    adios_error(err_expected_read_size_mismatch,
                                "Requested data is out of the local size. "
                                "Not yet supported by ICEE\n");
            }
        
            memcpy(data, vp->data, vp->varlen);
            break;
        }
        case ADIOS_SELECTION_AUTO:
        {
            // Possible memory overwrite
            memcpy(data, vp->data, vp->varlen);
            break;
        }
        case ADIOS_SELECTION_POINTS:
        {
            adios_error(err_operation_not_supported,
                        "ADIOS_SELECTION_POINTS not yet supported by ICEE.\n");
            break;
        }
        }

    return adios_errno;
}

int 
adios_read_icee_schedule_read(const ADIOS_FILE *adiosfile,
			const ADIOS_SELECTION * sel,
			const char * varname,
			int from_steps,
			int nsteps,
			void * data)
{
    log_error("No support yet: %s\n", __FUNCTION__);
    return 0;
}

int 
adios_read_icee_get_attr (int *gp, const char *attrname,
                                 enum ADIOS_DATATYPES *type,
                                 int *size, void **data)
{
    log_error("No support yet: %s\n", __FUNCTION__);
    return 0;
}

int 
adios_read_icee_get_attr_byid (const ADIOS_FILE *adiosfile, int attrid,
				   enum ADIOS_DATATYPES *type,
				   int *size, void **data)
{
    log_error("No support yet: %s\n", __FUNCTION__);
    return 0;
}

ADIOS_VARINFO* 
adios_read_icee_inq_var(const ADIOS_FILE * adiosfile, const char* varname)
{
    log_debug("%s (%s)\n", __FUNCTION__, varname);

    return NULL;
}

ADIOS_VARINFO* 
adios_read_icee_inq_var_byid (const ADIOS_FILE * adiosfile, int varid)
{
    log_debug("%s (%d)\n", __FUNCTION__, varid);

    icee_fileinfo_rec_ptr_t fp = (icee_fileinfo_rec_ptr_t) adiosfile->fh;
    assert(varid < fp->nvars);

    ADIOS_VARINFO *a = calloc(1, sizeof(ADIOS_VARINFO));
    
    icee_varinfo_rec_ptr_t vp = NULL;
    vp = icee_varinfo_search_byname(fp->varinfo, adiosfile->var_namelist[varid]);
    //icee_varinfo_print(vp);

    if (vp)
    {
        a->varid = vp->varid;
        a->type = vp->type;
        a->ndim = vp->ndims;

        uint64_t dimsize = vp->ndims * sizeof(uint64_t);
        a->dims = malloc(dimsize);
        memcpy(a->dims, vp->ldims, dimsize);

        if (vp->ndims == 0)
        {
            a->value = malloc(vp->typesize);
            memcpy(a->value, vp->data, vp->typesize);
        }
        else
            a->value = NULL;
    }

    return a;
}

void 
adios_read_icee_free_varinfo (ADIOS_VARINFO *adiosvar)
{
    log_debug("%s\n", __FUNCTION__);
    free(adiosvar);
    return;
}


ADIOS_TRANSINFO* 
adios_read_icee_inq_var_transinfo(const ADIOS_FILE *gp, 
                                  const ADIOS_VARINFO *vi)
{    
    log_debug("%s\n", __FUNCTION__);
    ADIOS_TRANSINFO *trans = malloc(sizeof(ADIOS_TRANSINFO));
    memset(trans, 0, sizeof(ADIOS_TRANSINFO));
    trans->transform_type = adios_transform_none;
    return trans;
}


int 
adios_read_icee_inq_var_trans_blockinfo(const ADIOS_FILE *gp, 
                                        const ADIOS_VARINFO *vi, 
                                        ADIOS_TRANSINFO *ti)
{
    log_error("No support yet: %s\n", __FUNCTION__);
    return 0;
}

void 
adios_read_icee_reset_dimension_order (const ADIOS_FILE *adiosfile, 
                                       int is_fortran)
{
    log_error("No support yet: %s\n", __FUNCTION__);
    return;
}
