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
#include "core/transforms/adios_transforms_common.h"
//#include "core/adios_icee.h"

// conditional libraries
#ifdef DMALLOC
#include "dmalloc.h"
#endif

#define MALLOC0(type_t, var) var = malloc(sizeof(type_t)); memset(var, '\0', sizeof(type_t));
#define MYALLOC(var, nbyte)                         \
    if (nbyte) {                                    \
        var = malloc(nbyte);                        \
        assert((var != NULL) && "malloc failure");  \
        memset(var, '\0', nbyte);                   \
    } else                                          \
        var = NULL;

#define MYFREE(var) free(var); var = NULL;

#define MYMIN(a,b)                              \
    ({ __typeof__ (a) _a = (a);                 \
        __typeof__ (b) _b = (b);                \
        _a < _b ? _a : _b; })

#define MYMAX(a,b)                              \
    ({ __typeof__ (a) _a = (a);                 \
        __typeof__ (b) _b = (b);                \
        _a > _b ? _a : _b; })

///////////////////////////
// Global Variables
///////////////////////////
#include "core/globals.h"

#define DUMP(fmt, ...) fprintf(stderr, ">>> "fmt"\n", ## __VA_ARGS__); 


/* Data Structure
   + fileinfo 1.1 - fileinfo 1.2 - ...
   |    + var 1.1      + var 1.2 
   |    + var 2.1      + var 2.2
   |    + ...          + ...
   + fileinfo 2.1 - fileinfo 2.2 - ...

   x.y (x=filename, y=version)
*/

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

void
icee_llist_map(const icee_llist_ptr_t root, void (*fp)(const void*))
{
    icee_llist_ptr_t p = root;
    
    while (p != NULL)
    {
        (*fp)(p->item);
        p = p->next;
    }
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
pthread_mutex_t fileinfo_lock;

void icee_fileinfo_append(const icee_fileinfo_rec_ptr_t root, icee_fileinfo_rec_ptr_t fp)
{
    assert(root != NULL);
    
    icee_fileinfo_rec_ptr_t p = root;

    while ((p->timestep != fp->timestep) && (p->next != NULL))
    {
        p = p->next; 
    }

    if (p->timestep == fp->timestep)
    {
        icee_varinfo_rec_ptr_t v = p->varinfo;
        
        while (v->next != NULL)
            v = v->next;

        v->next = fp->varinfo;
        p->merge_count++;

        free(fp);
    }
    else
    {
        p->next = fp;
    }

    /*
      while (p->next != NULL)
      {
      p = p->next;
      }

      p->next = fp;
    */
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

void icee_sel_bb_print(const ADIOS_SELECTION *sel)
{
    fprintf(stderr, "===== selection (%p) =====\n", sel);

    if (sel)
    {
        switch(sel->type)
        {
        case ADIOS_SELECTION_WRITEBLOCK:
            fprintf(stderr, "%10s : %s\n", "type", "writeblock");
            break;
        case ADIOS_SELECTION_BOUNDINGBOX:
            fprintf(stderr, "%10s : %s\n", "type", "boundingbox");
            fprintf(stderr, "%10s : %d\n", "ndims", sel->u.bb.ndim);
            icee_dims_print("start", sel->u.bb.ndim, sel->u.bb.start);
            icee_dims_print("count", sel->u.bb.ndim, sel->u.bb.count);
            break;
        case ADIOS_SELECTION_AUTO:
            fprintf(stderr, "%10s : %s\n", "type", "auto");
            break;
        case ADIOS_SELECTION_POINTS:
            fprintf(stderr, "%10s : %s\n", "type", "points");
            break;
        default:
            fprintf(stderr, "%10s : %s\n", "type", "undefined");
            break;
        }
    }
    else
    {
        fprintf(stderr, "selection is invalid\n");
    }
}

static int
icee_fileinfo_handler(CManager cm, void *vevent, void *client_data, attr_list attrs)
{
    icee_fileinfo_rec_ptr_t event = vevent;
    log_debug("%s (%s)\n", __FUNCTION__, event->fname);

    icee_fileinfo_rec_ptr_t lfp = malloc(sizeof(icee_fileinfo_rec_t));
    icee_fileinfo_copy(lfp, event);
    lfp->merge_count = 1;

    icee_varinfo_rec_ptr_t eventvp = event->varinfo;
    icee_varinfo_rec_ptr_t *lvpp = &lfp->varinfo;

    while (eventvp != NULL)
    {
        *lvpp = malloc(sizeof(icee_varinfo_rec_t));
        icee_varinfo_copy(*lvpp, eventvp);
        if (adios_verbose_level > 5) DUMP("id,name = %d,%s", (*lvpp)->varid, (*lvpp)->varname);
        
        lvpp = &(*lvpp)->next;
        eventvp = eventvp->next;
    }

    pthread_mutex_lock (&fileinfo_lock);
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
    pthread_mutex_unlock (&fileinfo_lock);

    if (adios_verbose_level > 5)
        icee_llist_map(icee_filelist, icee_fileinfo_print);

    return 1;
}

void on_icee_passivecheckin_reply (CManager cm, CMConnection conn, icee_passivecheckin_rec_t *m)
{
    log_debug("%s\n", __FUNCTION__);
}

void on_icee_fileinfo_recv (CManager cm, CMConnection conn, icee_fileinfo_rec_t *m)
{
    log_debug("%s\n", __FUNCTION__);
    icee_fileinfo_handler(cm, (void *)m, NULL, NULL);
}

static int adios_read_icee_initialized = 0;

//CManager icee_read_cm;
CManager icee_read_cm[ICEE_MAX_PARALLEL];
static int icee_read_num_parallel = 1;

static int is_read_cm_passive = 0;
static int num_remote_server = 0;
CManager pcm[ICEE_MAX_REMOTE];

/* Row-ordered matrix representatin */
typedef struct {
    int      typesize;
    int      ndims;
    uint64_t dims[10];
    uint64_t accumdims[10];
    void*    data;
} icee_matrix_t;

/* View (subset) representation */
typedef struct {
    uint64_t vdims[10];
    uint64_t offsets[10];
    int      leastcontiguousdim;
    icee_matrix_t* mat;
} icee_matrix_view_t;

void 
mat_init (icee_matrix_t *m, 
          int typesize,
          int ndims,
          const uint64_t *dims,
          void *data)
{
    assert(ndims <= 10);
    m->typesize = typesize;
    m->ndims = ndims;
    memcpy(m->dims, dims, ndims * sizeof(uint64_t));

    int i;
    uint64_t p = 1;
    for (i=ndims-1; i>=0; i--)
    {
        m->accumdims[i] = p;
        p *= dims[i];
    }

    m->data = data;
}

void 
view_init (icee_matrix_view_t *v,
           icee_matrix_t *m,
           const int64_t *vdims,
           const uint64_t *offsets)
{
    v->mat = m;
    memcpy(v->vdims, vdims, m->ndims * sizeof(uint64_t));
    memcpy(v->offsets, offsets, m->ndims * sizeof(uint64_t));

    int i;
    for (i=m->ndims-1; i>=0; i--)
    {
        v->leastcontiguousdim = i;
        if ((offsets[i] != 0) || (vdims[i] != m->dims[i]))
            break;
    }
}

/* Copy data between two views. Dimension and size should match */
void
view_copy (icee_matrix_view_t *dest, icee_matrix_view_t *src)
{
    assert(dest->mat->ndims == src->mat->ndims);
    
    int i;
    for (i=0; i<dest->mat->ndims; i++)
        assert(dest->vdims[i] == src->vdims[i]);

    // Contiguous merging
    if ((dest->leastcontiguousdim == 0) && (src->leastcontiguousdim == 0))
    {
        int s, d;
        d = dest->offsets[0] * dest->mat->accumdims[0];
        s = src->offsets[0] * src->mat->accumdims[0];
        memcpy(dest->mat->data + d * dest->mat->typesize, 
               src->mat->data + s * dest->mat->typesize, 
               dest->vdims[0] * dest->mat->accumdims[0] * dest->mat->typesize);

        return;
    }
    
    // Non-contiguous merging
    switch (dest->mat->ndims)
    {
    case 1:
    {
        int s, d;
        d = dest->offsets[0];
        s = src->offsets[0];
        memcpy(dest->mat->data + d * dest->mat->typesize, 
               src->mat->data + s * dest->mat->typesize, 
               dest->vdims[0] * dest->mat->typesize);
        break;
    }
    case 2:
    {
        int i, s, d;
        for (i=0; i<dest->vdims[0]; i++)
        {
            d = (i + dest->offsets[0]) * dest->mat->accumdims[0] 
                + dest->offsets[1];
            s = (i + src->offsets[0]) * src->mat->accumdims[0] 
                + src->offsets[1];
            memcpy(dest->mat->data + d * dest->mat->typesize, 
                   src->mat->data + s * dest->mat->typesize, 
                   dest->vdims[1] * dest->mat->typesize);
        }
        break;
    }
    case 3:
    {
        int i, j, s, d;
        for (i=0; i<dest->vdims[0]; i++)
        {
            for (j=0; j<dest->vdims[1]; j++)
            {
                d = (i + dest->offsets[0]) * dest->mat->accumdims[0] 
                    + (j + dest->offsets[1]) * dest->mat->accumdims[1]
                    + dest->offsets[2];
                s = (i + src->offsets[0]) * src->mat->accumdims[0] 
                    + (j + src->offsets[1]) * src->mat->accumdims[1] 
                    + src->offsets[2];
                memcpy(dest->mat->data + d * dest->mat->typesize, 
                       src->mat->data + s * dest->mat->typesize, 
                       dest->vdims[2] * dest->mat->typesize);
            }
        }
        break;
    }
    default:
        adios_error(err_expected_read_size_mismatch,
                    "The variable dimension is out of the range. ",
                    "Not yet supported by ICEE\n");
        break;
    }
}

void 
set_contact_list (attr_list contact_list, 
                  const icee_transport_t icee_transport, 
                  char *host, int port)
{
    switch (icee_transport)
    {
    case ENET:
        add_string_attr(contact_list, 
                        attr_atom_from_string("CM_TRANSPORT"), 
                        "enet");
        add_string_attr(contact_list, 
                        attr_atom_from_string("CM_ENET_HOST"), 
                        host);
        add_int_attr(contact_list, 
                     attr_atom_from_string("CM_ENET_PORT"), 
                     port);
        break;
    case NNTI:
        add_string_attr(contact_list, 
                        attr_atom_from_string("CM_TRANSPORT"), 
                        "nnti");
        add_string_attr(contact_list, 
                        attr_atom_from_string("IP_HOST"), 
                        host);
        add_int_attr(contact_list, 
                     attr_atom_from_string("NNTI_PORT"), 
                     port);
        add_string_attr(contact_list, 
                        attr_atom_from_string("CM_NNTI_TRANSPORT"), 
                        "ib");
        break;
    case IB:
        add_string_attr(contact_list, 
                        attr_atom_from_string("CM_TRANSPORT"), 
                        "ib");
        add_string_attr(contact_list, 
                        attr_atom_from_string("IP_HOST"), 
                        host);
        add_int_attr(contact_list, 
                     attr_atom_from_string("IP_PORT"), 
                     port);
        break;
    default:
        add_string_attr(contact_list, 
                        attr_atom_from_string("IP_HOST"), 
                        host);
        add_int_attr(contact_list, 
                     attr_atom_from_string("IP_PORT"), 
                     port);
        break;
            
    }
}

/********** Core ADIOS Read functions. **********/

int
adios_read_icee_init_method (MPI_Comm comm, PairStruct* params)
{   
    log_debug ("%s\n", __FUNCTION__);

    int cm_port = 59997;
    char *cm_host = "localhost";
    int cm_remote_port = 59999;
    char *cm_remote_host = "localhost";
    char *cm_attr = NULL;
    //attr_list contact_list;
    icee_transport_t icee_transport_init = TCP;
    icee_transport_t icee_transport = TCP;

    icee_contactinfo_rec_t *remote_contact = NULL;
    int i;

    int use_single_remote_server = 1;
    char *remote_list_str = NULL;
    char *attr_list_str = NULL;

    int use_native_contact = 0;

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
        else if (!strcasecmp (p->name, "cm_remote_host"))
        {
            cm_remote_host = p->value;
        }
        else if (!strcasecmp (p->name, "cm_remote_port"))
        {
            cm_remote_port = atoi(p->value);
        }
        else if (!strcasecmp (p->name, "remote_list"))
        {
            use_single_remote_server = 0;
            if (p->value)
                remote_list_str = strdup(p->value);
        }
        else if (!strcasecmp (p->name, "attr_list"))
        {
            use_single_remote_server = 0;
            if (p->value)
                attr_list_str = strdup(p->value);
        }
        else if (!strcasecmp (p->name, "transport"))
        {
            if (strcasecmp(p->value, "TCP") == 0)
                icee_transport = TCP;
            else if (strcasecmp(p->value, "ENET") == 0)
                icee_transport = ENET;
            else if (strcasecmp(p->value, "NNTI") == 0)
                icee_transport = NNTI;
            else if (strcasecmp(p->value, "IB") == 0)
                icee_transport = IB;
            else
                log_error ("No support: %s\n", p->value);
        }
        else if (!strcasecmp (p->name, "transport_init"))
        {
            if (strcasecmp(p->value, "TCP") == 0)
                icee_transport_init = TCP;
            else if (strcasecmp(p->value, "ENET") == 0)
                icee_transport_init = ENET;
            else if (strcasecmp(p->value, "NNTI") == 0)
                icee_transport_init = NNTI;
            else if (strcasecmp(p->value, "IB") == 0)
                icee_transport_init = IB;
            else
                log_error ("No support: %s\n", p->value);
        }
        else if (!strcasecmp (p->name, "num_parallel"))
        {
            icee_read_num_parallel = atoi(p->value);
        }
        else if (!strcasecmp (p->name, "is_passive"))
        {
            is_read_cm_passive = atoi(p->value);
        }
        else if (!strcasecmp (p->name, "use_native_contact"))
        {
            use_native_contact = atoi(p->value);
        }

        p = p->next;
    }

    pthread_mutex_init(&fileinfo_lock, NULL);

    if (use_single_remote_server)
    {
        num_remote_server = 1;

        attr_list contact_list = create_attr_list();
        set_contact_list(contact_list, icee_transport_init, cm_remote_host, cm_remote_port);

        icee_contactinfo_rec_t *p;
        p = malloc(sizeof(icee_contactinfo_rec_t));
        char *contact_string = attr_list_to_string(contact_list);
        p->contact_string = contact_string;
        p->stone_id = 0; // we assume. it can be wrong.
        p->next = NULL;

        remote_contact = p;
    }
    else
    {
        num_remote_server = 0;

        icee_contactinfo_rec_t *p;
        icee_contactinfo_rec_t *prev;

        char* token = strtok(remote_list_str, ",");
        while (token)
        {
            char host[256];
            int port = 0;

            if (token[0] == ':')
            {
                strcpy(host, cm_remote_host);
                port = atoi(token+1);
            }
            else
            {
                char *pch = strchr(token, ':');
                if (pch != NULL)
                {
                    strncpy(host, token, pch - token);
                    host[pch-token] = '\0';
                    port = atoi(pch+1);
                }
                else
                {
                    int len = strlen(token);
                    strncpy(host, token, len);
                    assert(len < 256);
                    host[len] = '\0';
                    port = cm_remote_port;
                }
            }
            log_debug("Remote server list: (%d) %s:%d\n", num_remote_server, host, port);

            p = malloc(sizeof(icee_contactinfo_rec_t));
            attr_list contact_list;
            
            contact_list = create_attr_list();
            set_contact_list(contact_list, icee_transport_init, host, port);
            p->contact_string = attr_list_to_string(contact_list);
            p->stone_id = 0; // we assume. it can be wrong.
            p->next = NULL;

            if (num_remote_server == 0)
                remote_contact = p;
            else
                prev->next = p;

            prev = p;
            num_remote_server++;
            token = strtok(NULL, ",");
        }
    }

    if (attr_list_str != NULL)
    {
        num_remote_server = 0;

        icee_contactinfo_rec_t *p;
        icee_contactinfo_rec_t *prev;

        char* token = strtok(attr_list_str, ",");
        while (token)
        {
            int remote_stone = 0;
            char string_list[256];

            sscanf(token, "%d:%s", &remote_stone, &string_list[0]);

            p = malloc(sizeof(icee_contactinfo_rec_t));
            attr_list contact_list;
            
            p->stone_id = remote_stone;
            p->contact_string = strdup(string_list);
            p->next = NULL;

            if (num_remote_server == 0)
                remote_contact = p;
            else
                prev->next = p;

            prev = p;
            num_remote_server++;
            token = strtok(NULL, ",");
        }
    }

    if (icee_read_num_parallel > ICEE_MAX_PARALLEL)
    {
        icee_read_num_parallel = ICEE_MAX_PARALLEL;
        log_info ("Max. number of threads is set to %d\n", icee_read_num_parallel);
    }

    log_debug ("transport : %s\n", icee_transport_name[icee_transport]);

    /*
      log_info ("cm_host : %s\n", cm_host);
      log_info ("cm_port : %d\n", cm_port);

      for (i = 0; i < num_remote_server; i++) 
      {
      log_info ("remote_list : %s:%d\n", remote_server[i].client_host, remote_server[i].client_port);
      }
    */

    if (!adios_read_icee_initialized)
    {
        if (is_read_cm_passive)
        {
            icee_contactinfo_rec_t *prev;
            for (i = 0; i < num_remote_server; i++)
            {
                attr_list contact_list;
                icee_contactinfo_rec_t *p = (i == 0)? remote_contact : prev->next;

                pcm[i] = CManager_create();

                if (!CMfork_comm_thread(pcm[i]))
                    printf("Fork of communication thread[%d] failed.\n", i);

                contact_list = attr_list_from_string(p->contact_string);

                log_debug("Passive remote contact: \"%s\"\n", attr_list_to_string(contact_list));
                if (adios_verbose_level > 5) dump_attr_list(contact_list);

                /*
                  attr_list contact_list  = create_attr_list();
                  add_string_attr(contact_list, attr_atom_from_string("IP_HOST"),
                  remote_server[i].client_host);
                  add_int_attr(contact_list, attr_atom_from_string("IP_PORT"),
                  remote_server[i].client_port);
                */

                CMConnection conn = CMinitiate_conn(pcm[i], contact_list);

                int n = 0;
                while (conn == NULL)
                {
                    log_error ("Passive connection failed (%d). Try again ...\n", i);
                    dump_attr_list(contact_list);

                    sleep(2);
                    conn = CMinitiate_conn(pcm[i], contact_list);
                    
                    if (n > 5) break;
                    n++;
                }

                if (conn == NULL)
                {
                    log_error ("Initializing passive connection failed (%d)\n", i);
                    
                }

                CMFormat fm_checkin, fm_fileinfo;
                fm_checkin = CMregister_format(pcm[i], icee_passivecheckin_format_list);

                CMregister_handler(fm_checkin, icee_passivecheckin_reply_handler, on_icee_passivecheckin_reply);

                fm_fileinfo = CMregister_format(pcm[i], icee_fileinfo_format_list);
                CMregister_handler(fm_fileinfo, icee_fileinfo_recv_handler, on_icee_fileinfo_recv);

                icee_passivecheckin_rec_t m;

                int condition = CMCondition_get(pcm[i], conn);
                CMCondition_set_client_data(pcm[i], condition, NULL);
                m.condition = condition;

                if (CMwrite(conn, fm_checkin, (void*)&m) != 1)
                    log_error ("Passive check-in failed (%d)\n", i);

                prev = p;
            }
            log_debug("Passive connection established");
            goto done;
        }

        EVstone   stone[ICEE_MAX_PARALLEL], remote_stone;
        EVsource  source;
        attr_list contact[ICEE_MAX_PARALLEL];

        icee_contactinfo_rec_t contact_msg[ICEE_MAX_PARALLEL];

        for (i=0; i<icee_read_num_parallel; i++)
        {
            icee_read_cm[i] = CManager_create();

            contact[i] = create_attr_list();
            set_contact_list(contact[i], icee_transport, cm_host, cm_port+i);

            if (CMlisten_specific(icee_read_cm[i], contact[i]) == 0)
                printf("Error: unable to initialize connection manager[%d].\n", i);
            
            if (!CMfork_comm_thread(icee_read_cm[i])) 
                printf("Fork of communication thread[%d] failed.\n", i);

            stone[i] = EValloc_stone(icee_read_cm[i]);
            if (adios_verbose_level > 5) 
            {
                log_debug("Reader contact: \"%d:%s\"\n", stone[i], attr_list_to_string(CMget_contact_list(icee_read_cm[i])));
                dump_attr_list(CMget_contact_list(icee_read_cm[i]));
            }

            EVassoc_terminal_action(icee_read_cm[i], stone[i], icee_fileinfo_format_list, icee_fileinfo_handler, NULL);

            contact_msg[i].stone_id = stone[i];
            attr_list contact_list;
            if (use_native_contact)
                contact_list = CMget_contact_list(icee_read_cm[i]);
            else
                contact_list = contact[i];
                
            contact_msg[i].contact_string = attr_list_to_string(contact_list);
            contact_msg[i].next = NULL;

            if (i>0)
                contact_msg[i-1].next = &contact_msg[i];
        }
        
        EVstone split_stone;
        EVaction split_action;

        split_stone = EValloc_stone(icee_read_cm[0]);
        split_action = EVassoc_split_action(icee_read_cm[0], split_stone, NULL);
        icee_contactinfo_rec_t *prev;
        for (i = 0; i < num_remote_server; i++) 
        {
            attr_list contact_list;
            EVstone remote_stone, output_stone;
            output_stone = EValloc_stone(icee_read_cm[0]);
            icee_contactinfo_rec_t *p = (i == 0)? remote_contact : prev->next;
            
            remote_stone = p->stone_id;
            contact_list = attr_list_from_string(p->contact_string);

            EVaction action;
            action = EVassoc_bridge_action(icee_read_cm[0], output_stone, contact_list, remote_stone);

            int n = 0;
            while (action == -1)
            {
                log_error ("Connection failed (%d). Try again ...\n", i);
                dump_attr_list(contact_list);
                
                sleep(2);
                action = EVassoc_bridge_action(icee_read_cm[0], output_stone, contact_list, remote_stone);
                
                if (n > 5) break;
                n++;
            }

            EVaction_add_split_target(icee_read_cm[0], split_stone, split_action, output_stone);

            prev = p;

            log_debug("Remote contact: \"%d:%s\"\n", remote_stone, attr_list_to_string(contact_list));
            if (adios_verbose_level > 5) dump_attr_list(contact_list);

        }

        source = EVcreate_submit_handle(icee_read_cm[0], split_stone, icee_contactinfo_format_list);
        
        //if (adios_verbose_level > 5) icee_contactinfo_print(contact_msg);
        
        EVsubmit(source, contact_msg, NULL);

    done:
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
    
    while (fp->merge_count != fp->nchunks)
    {
        log_debug("Waiting the rest of blocks (%d/%d)\n", fp->merge_count, fp->nchunks);
        
        usleep(0.1*1E6);
    }

    ADIOS_FILE *adiosfile = malloc(sizeof(ADIOS_FILE));

    adiosfile->fh = (uint64_t)fp;

#if 0
    int hashsize = 10;
    qhashtbl_t *tbl = qhashtbl(hashsize);
    
    icee_varinfo_rec_ptr_t vp = fp->varinfo;
    while (vp != NULL)
    {
        tbl->put(tbl, vp->varname, NULL);
        vp = vp->next;
    }

    adiosfile->nvars = tbl->size(tbl);
    adiosfile->var_namelist = malloc(tbl->num * sizeof(char *));

    int i;
    size_t idx = 0;
    for (i=0; i<hashsize; i++)
    {
        qhnobj_t *qhobj = tbl->slots[i].head;
        while (qhobj != NULL)
        {
            //DUMP("%lu:%s", idx, qhobj->key);
            adiosfile->var_namelist[idx] = strdup(qhobj->key);
            qhobj->value = (void*) idx;
            idx++;
            qhobj = qhobj->next;
        }
    }

    /*
    // Update back to varinfo
    vp = fp->varinfo;
    while (vp != NULL)
    {
    vp->varid = (int) tbl->get(tbl, vp->varname);
    vp = vp->next;
    }
    */

    tbl->free(tbl);
#endif

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
        int i;
        if (is_read_cm_passive == 1)
            for (i = 0; i < num_remote_server; i++) 
                CManager_close(pcm[i]);
        else
            for (i=0; i<icee_read_num_parallel; i++)
                CManager_close(icee_read_cm[i]);

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
    int i;
    icee_fileinfo_rec_ptr_t fp = (icee_fileinfo_rec_ptr_t) adiosfile->fh;
    log_debug("%s (%d:%s)\n", __FUNCTION__, varid, fp->fname);
    //assert((varid < fp->nvars) || (fp->nvars == 0));

    if (nsteps != 1)
    {
        adios_error (err_invalid_timestep,
                     "Only one step can be read from a stream at a time. "
                     "You requested % steps in adios_schedule_read()\n", 
                     nsteps);
        return err_invalid_timestep;
    }
    
    icee_varinfo_rec_ptr_t vp = NULL;
    vp = icee_varinfo_search_byname(fp->varinfo, adiosfile->var_namelist[varid]);
    if (adios_verbose_level > 5) icee_varinfo_print(vp);

    if (!vp)
    {
        adios_error(err_invalid_varid, "Invalid variable id: %d\n", varid);
        return adios_errno;
    }

    while (fp->merge_count != fp->nchunks)
    {
        log_debug("Waiting the rest of blocks (%d/%d)\n", fp->merge_count, fp->nchunks);
        
        usleep(0.1*1E6);
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

            if (fp->comm_rank != sel->u.block.index)
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

            log_debug("Merging operation (total nvars: %d).\n", fp->nchunks);
            if (adios_verbose_level > 5) icee_sel_bb_print(sel);

            while (vp != NULL)
            {
                icee_matrix_t m_sel;
                icee_matrix_t m_var;
                icee_matrix_view_t v_sel;
                icee_matrix_view_t v_var;
                uint64_t start[10];
                int64_t count[10]; // should be signed to check validity
                uint64_t s_offsets[10], v_offsets[10];
                int i;

                if (adios_verbose_level > 5) icee_varinfo_print(vp);
                    
                mat_init(&m_sel, vp->typesize, vp->ndims, sel->u.bb.count, data);
                mat_init(&m_var, vp->typesize, vp->ndims, vp->ldims, vp->data);

                for (i=0; i<vp->ndims; i++)
                    start[i] = MYMAX(sel->u.bb.start[i], vp->offsets[i]);

                for (i=0; i<vp->ndims; i++)
                {
                    count[i] = 
                        MYMIN(sel->u.bb.start[i]+sel->u.bb.count[i],
                              vp->offsets[i]+vp->ldims[i]) - start[i];
                }
                    
                for (i=0; i<vp->ndims; i++)
                {
                    if (count[i] <= 0)
                    {
                        log_debug("No ROI. Skip\n");
                        goto next;
                    }
                }

                for (i=0; i<vp->ndims; i++)
                    s_offsets[i] = start[i] - sel->u.bb.start[i];

                for (i=0; i<vp->ndims; i++)
                    v_offsets[i] = start[i] - vp->offsets[i];

                view_init (&v_sel, &m_sel, count, s_offsets);
                view_init (&v_var, &m_var, count, v_offsets);
                view_copy (&v_sel, &v_var);
                    
            next:
                vp = icee_varinfo_search_byname(vp->next, adiosfile->var_namelist[varid]);
            }
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
    //assert((varid < fp->nvars) || (fp->nvars == 0));

    ADIOS_VARINFO *a = calloc(1, sizeof(ADIOS_VARINFO));
    
    icee_varinfo_rec_ptr_t vp = NULL;
    vp = icee_varinfo_search_byname(fp->varinfo, adiosfile->var_namelist[varid]);
    //icee_varinfo_print(vp);

    if (vp)
    {
        a->varid = vp->varid;
        a->type = vp->type;
        a->ndim = vp->ndims;
        a->nsteps = 1;
        
        if (vp->ndims == 0)
        {
            a->value = malloc(vp->typesize);
            memcpy(a->value, vp->data, vp->typesize);
        }
        else
        {
            uint64_t dimsize = vp->ndims * sizeof(uint64_t);
            a->dims = malloc(dimsize);
            memcpy(a->dims, vp->gdims, dimsize);
            a->global = 1;
            
            if (a->dims[0] == 0)
            {
                a->global = 0;
                memcpy(a->dims, vp->ldims, dimsize);
            }
            
            a->value = NULL;
        }
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

int
adios_read_icee_get_dimension_order (const ADIOS_FILE *adiosfile)
{
    log_error("No support yet: %s\n", __FUNCTION__);
    return 0; // Stub method: always return false (C order)
}

void
adios_read_icee_reset_dimension_order (const ADIOS_FILE *adiosfile, 
                                       int is_fortran)
{
    log_error("No support yet: %s\n", __FUNCTION__);
    return;
}
