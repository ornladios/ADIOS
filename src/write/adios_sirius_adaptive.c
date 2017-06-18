#include <unistd.h>
#include <fcntl.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <errno.h>
#include <assert.h>
#include <sys/stat.h>
#include <float.h>
// xml parser
#include <mxml.h>
#include <glib.h>
#include <zfp.h>
#include <omp.h>

// see if we have MPI or other tools
#include "config.h"

#include "public/adios.h"
#include "public/adios_types.h"
#include "public/adios_error.h"
#include "core/adios_transport_hooks.h"
#include "core/adios_internals.h"
#include "core/adios_internals_mxml.h"
#include "core/adios_logger.h"
#include "core/common_adios.h"
#include "core/util.h"

#include "pqueue.h"

#include "mpi.h"

#define MAXLEVEL 10
#define MAX_NODE_DEGREE 100
//#define DUMP_FEATURE
//#define TEST_REDUCTION

static char *io_method[MAXLEVEL]; //the IO methods for data output for each level
static char *io_parameters[MAXLEVEL]; //the IO method parameters
static char *io_paths[MAXLEVEL]; //the IO method output paths (prefix to filename)
static int nlevels=2; // Number of levels

double threshold = 0.1;
double compr_tolerance = 0.1;
int save_delta = 1;
int compress_delta = 1;
int dec_ratio = 10;

GHashTable ** nodes_ght = 0;

typedef struct node_t
{
	pqueue_pri_t pri;
	int    val;
	size_t pos;
} node_t;

typedef struct edge_cost_t
{
    node_t * pq_node;
    double cost;
} edge_cost_t;

typedef enum {NONE = 0, ABSOLUTE = 1, GRAD} op_type;

op_type thresh_type = NONE;

int
decompress (double * array, int nx, double tolerance,
            double * array_compressed, size_t array_size_compressed)
{
    zfp_type type;     /* array scalar type */
    zfp_field* field;  /* array meta data */
    zfp_stream* zfp;   /* compressed stream */
    bitstream * stream; /* bit stream to write to or read from */
    size_t zfpsize;    /* byte size of compressed stream */

    /* allocate meta data for the 3D array a[nz][ny][nx] */
    type = zfp_type_double;
    field = zfp_field_1d (array, type, nx);

    /* allocate meta data for a compressed stream */
    zfp = zfp_stream_open (NULL);
    zfp_stream_set_accuracy (zfp, tolerance, type);

    /* associate bit stream with allocated buffer */
    stream = stream_open (array_compressed, array_size_compressed);
    zfp_stream_set_bit_stream (zfp, stream);
    zfp_stream_rewind (zfp);

    assert (zfp_decompress(zfp, field));

    /* clean up */
    zfp_field_free (field);
    zfp_stream_close (zfp);
    stream_close (stream);

    return 0;
}

int
compress (double * array, int nx, double tolerance,
          double ** array_compressed)
{
    zfp_type type;     /* array scalar type */
    zfp_field* field;  /* array meta data */
    zfp_stream* zfp;   /* compressed stream */
    void* buffer;      /* storage for compressed stream */
    size_t bufsize;    /* byte size of compressed buffer */
    bitstream * stream; /* bit stream to write to or read from */
    size_t zfpsize;    /* byte size of compressed stream */

    /* allocate meta data for the 3D array a[nz][ny][nx] */
    type = zfp_type_double;
    field = zfp_field_1d (array, type, nx);

    /* allocate meta data for a compressed stream */
    zfp = zfp_stream_open (NULL);

    /* set compression mode and parameters via one of three functions */
    /*  zfp_stream_set_rate(zfp, rate, type, 3, 0); */
    /*  zfp_stream_set_precision(zfp, precision, type); */
    zfp_stream_set_accuracy (zfp, tolerance, type);

    /* allocate buffer for compressed data */
    bufsize = zfp_stream_maximum_size (zfp, field);
    buffer = malloc (bufsize);
    assert (buffer);

    /* associate bit stream with allocated buffer */
    stream = stream_open (buffer, bufsize);
    zfp_stream_set_bit_stream (zfp, stream);
    zfp_stream_rewind (zfp);

    /* compress array and output compressed stream */
    zfpsize = zfp_compress (zfp, field);
    assert (zfpsize);

    /* clean up */
    zfp_field_free (field);
    zfp_stream_close (zfp);
    stream_close (stream);

    * array_compressed = (double *) buffer;

    return zfpsize;
}

static int
cmp_pri(pqueue_pri_t next, pqueue_pri_t curr)
{
	return (next >= curr);
}


static pqueue_pri_t
get_pri(void *a)
{
	return ((node_t *) a)->pri;
}


static void
set_pri(void *a, pqueue_pri_t pri)
{
	((node_t *) a)->pri = pri;
}


static size_t
get_pos(void *a)
{
	return ((node_t *) a)->pos;
}


static void
set_pos(void *a, size_t pos)
{
	((node_t *) a)->pos = pos;
}

#define left(i)   ((i) << 1)
#define right(i)  (((i) << 1) + 1)
#define parent(i) ((i) >> 1)


pqueue_t *
pqueue_init(size_t n,
            pqueue_cmp_pri_f cmppri,
            pqueue_get_pri_f getpri,
            pqueue_set_pri_f setpri,
            pqueue_get_pos_f getpos,
            pqueue_set_pos_f setpos)
{
    pqueue_t *q;

    if (!(q = malloc(sizeof(pqueue_t))))
        return NULL;

    /* Need to allocate n+1 elements since element 0 isn't used. */
    if (!(q->d = malloc((n + 1) * sizeof(void *)))) {
        free(q);
        return NULL;
    }

    q->size = 1;
    q->avail = q->step = (n+1);  /* see comment above about n+1 */
    q->cmppri = cmppri;
    q->setpri = setpri;
    q->getpri = getpri;
    q->getpos = getpos;
    q->setpos = setpos;

    return q;
}


void
pqueue_free(pqueue_t *q)
{
    free(q->d);
    free(q);
}


size_t
pqueue_size(pqueue_t *q)
{
    /* queue element 0 exists but doesn't count since it isn't used. */
    return (q->size - 1);
}


static void
bubble_up(pqueue_t *q, size_t i)
{
    size_t parent_node;
    void *moving_node = q->d[i];
    pqueue_pri_t moving_pri = q->getpri(moving_node);

    for (parent_node = parent(i);
         ((i > 1) && q->cmppri(q->getpri(q->d[parent_node]), moving_pri));
         i = parent_node, parent_node = parent(i))
    {
        q->d[i] = q->d[parent_node];
        q->setpos(q->d[i], i);
    }

    q->d[i] = moving_node;
    q->setpos(moving_node, i);
}


static size_t
maxchild(pqueue_t *q, size_t i)
{
    size_t child_node = left(i);

    if (child_node >= q->size)
        return 0;

    if ((child_node+1) < q->size &&
        q->cmppri(q->getpri(q->d[child_node]), q->getpri(q->d[child_node+1])))
        child_node++; /* use right child instead of left */

    return child_node;
}


static void
percolate_down(pqueue_t *q, size_t i)
{
    size_t child_node;
    void *moving_node = q->d[i];
    pqueue_pri_t moving_pri = q->getpri(moving_node);

    while ((child_node = maxchild(q, i)) &&
           q->cmppri(moving_pri, q->getpri(q->d[child_node])))
    {
        q->d[i] = q->d[child_node];
        q->setpos(q->d[i], i);
        i = child_node;
    }

    q->d[i] = moving_node;
    q->setpos(moving_node, i);
}


int
pqueue_insert(pqueue_t *q, void *d)
{
    void *tmp;
    size_t i;
    size_t newsize;

    if (!q) return 1;

    /* allocate more memory if necessary */
    if (q->size >= q->avail) {
        newsize = q->size + q->step;
        if (!(tmp = realloc(q->d, sizeof(void *) * newsize)))
            return 1;
        q->d = tmp;
        q->avail = newsize;
    }

    /* insert item */
    i = q->size++;
    q->d[i] = d;
    bubble_up(q, i);

    return 0;
}


void
pqueue_change_priority(pqueue_t *q,
                       pqueue_pri_t new_pri,
                       void *d)
{
    size_t posn;
    pqueue_pri_t old_pri = q->getpri(d);

    q->setpri(d, new_pri);
    posn = q->getpos(d);
    if (q->cmppri(old_pri, new_pri))
        bubble_up(q, posn);
    else
        percolate_down(q, posn);
}


int
pqueue_remove(pqueue_t *q, void *d)
{
    size_t posn = q->getpos(d);
    q->d[posn] = q->d[--q->size];
    if (q->cmppri(q->getpri(d), q->getpri(q->d[posn])))
        bubble_up(q, posn);
    else
        percolate_down(q, posn);

    return 0;
}


void *
pqueue_pop(pqueue_t *q)
{
    void *head;

    if (!q || q->size == 1)
        return NULL;

    head = q->d[1];
    q->d[1] = q->d[--q->size];
    percolate_down(q, 1);

    return head;
}


void *
pqueue_peek(pqueue_t *q)
{
    void *d;
    if (!q || q->size == 1)
        return NULL;
    d = q->d[1];
    return d;
}


void
pqueue_dump(pqueue_t *q,
            FILE *out,
            pqueue_print_entry_f print)
{
    int i;

    fprintf(stdout,"posn\tleft\tright\tparent\tmaxchild\t...\n");
    for (i = 1; i < q->size ;i++) {
        fprintf(stdout,
                "%d\t%d\t%d\t%d\t%ul\t",
                i,
                left(i), right(i), parent(i),
                (unsigned int)maxchild(q, i));
        print(out, q->d[i]);
    }
}

#if 0
static void
set_pos(void *d, size_t val)
{
    /* do nothing */
}


static void
set_pri(void *d, pqueue_pri_t pri)
{
    /* do nothing */
}
#endif

void
pqueue_print(pqueue_t *q,
             FILE *out,
             pqueue_print_entry_f print)
{
    pqueue_t *dup;
	void *e;

    dup = pqueue_init(q->size,
                      q->cmppri, q->getpri, set_pri,
                      q->getpos, set_pos);
    dup->size = q->size;
    dup->avail = q->avail;
    dup->step = q->step;

    memcpy(dup->d, q->d, (q->size * sizeof(void *)));

    while ((e = pqueue_pop(dup)))
		print(out, e);

    pqueue_free(dup);
}


static int
subtree_is_valid(pqueue_t *q, int pos)
{
    if (left(pos) < q->size) {
        /* has a left child */
        if (q->cmppri(q->getpri(q->d[pos]), q->getpri(q->d[left(pos)])))
            return 0;
        if (!subtree_is_valid(q, left(pos)))
            return 0;
    }
    if (right(pos) < q->size) {
        /* has a right child */
        if (q->cmppri(q->getpri(q->d[pos]), q->getpri(q->d[right(pos)])))
            return 0;
        if (!subtree_is_valid(q, right(pos)))
            return 0;
    }
    return 1;
}


int
pqueue_is_valid(pqueue_t *q)
{
    return subtree_is_valid(q, 1);
}

struct var_struct
{
    char * name;
    char * path;
    enum ADIOS_DATATYPES type;
    enum ADIOS_FLAG multidim;
    char * global_dimensions;
    char * local_dimensions;
    char * local_offsets;
    void * data;
    uint64_t size; // in bytes

    struct var_struct *prev;
    struct var_struct *next;
};

struct level_struct
{
    int64_t fd;                        // ADIOS file descriptor to this level's output
    char *filename;                    // full path to this level's output
    char *grp_name;                    // each level has its own group name and group structure
    int64_t grp;
    int varcnt;                        // number of variables going into this level
    struct var_struct *vars;      // last inserted variable into this level
    struct var_struct *vars_head; // starting of variables in this level
    uint64_t totalsize;  // size of variables in this level
    pthread_t thread;
};

struct adios_sa_data_struct
{
    int64_t fpr;
    MPI_Comm group_comm;
    int rank;
    int size;
    void *comm;
    struct adios_bp_buffer_struct_v1 b;
    struct adios_group_struct * group;
    char * file_mode;

    struct level_struct level[MAXLEVEL];
};



// temporary solution for compiling error
static int declare_group (int64_t * id, const char * name
                          ,const char * time_index
                          ,enum ADIOS_FLAG stats
                         )
{
    int ret = adios_common_declare_group (id, name, adios_flag_no
                                      ,""
                                      ,""
                                      ,time_index
                                      ,adios_flag_no
                                      );
    if (ret == 1)
    {
        struct adios_group_struct * g = (struct adios_group_struct *) *id;
        g->all_unique_var_names = adios_flag_no;
    }
    
    return ret;
}

// temporary solution for compiling error
static int select_method (int64_t group, const char * method
                         ,const char * parameters
                         ,const char * base_path
                         )
{
    return adios_common_select_method_by_group_id (0 
                                                  ,method
                                                  ,parameters
                                                  ,group
                                                  ,base_path
                                                  ,0
                                                  );
}

static void define_iogroups (struct adios_method_struct * method)
{
    int len, l;
    struct adios_sa_data_struct * md = (struct adios_sa_data_struct *)
                             method->method_data;
    
    for (l = 0; l < nlevels; l++)
    {
        len = 5 + strlen (method->group->name); //new groupname= tg_groupname
        md->level[l].grp_name = (char *)malloc (len);
        memset (md->level[l].grp_name, 0x00, len);
        sprintf (md->level[l].grp_name, "%s_L%d",method->group->name, l);
        declare_group (&(md->level[l].grp), md->level[l].grp_name, "", adios_flag_yes);
        select_method (md->level[l].grp, io_method[l], io_parameters[l],"");
    }
}

static int convert_file_mode(enum ADIOS_METHOD_MODE mode, char * file_mode)
{
    switch (mode)
    {
        case adios_mode_read:
            strcpy (file_mode,"r");
            break;

        case adios_mode_write:
            strcpy (file_mode,"w");
            break;

        case adios_mode_append:
            strcpy (file_mode,"a");
            break;

        case adios_mode_update:
            strcpy (file_mode,"u");
            break;
        default:
            fprintf (stderr, "adios_open: unknown file mode: %s\n", file_mode);
            return -1;
            break;
    }

    return 0;
}


static void init_output_parameters(const PairStruct *params)
{
    const PairStruct *p = params;
    nlevels = 0;
    int level_params = 0;
    int level_paths = 0;

    while (p) {
        if (!strcasecmp (p->name, "method")) {
            errno = 0;
            io_method[nlevels] = strdup (p->value);
            if (!errno) {
                log_debug ("method %d set to %s for SIRIUS method\n", nlevels, io_method[nlevels]);
            } else {
                log_error ("Invalid 'method' parameter given to the SIRIUS method: '%s'\n", p->value);
                io_method[nlevels] = NULL;
            }
            nlevels++;
        } else if (!strcasecmp (p->name, "parameters")) {
            errno = 0;
            if(p->value)
                io_parameters[level_params] = strdup (p->value);
            else
                io_parameters[level_params] = strdup (" ");
            if (!errno) {
                log_debug ("parameters %d set to %s for SIRIUS method\n", level_params, io_parameters[level_params]);
            } else {
                log_error ("Invalid 'parameters' parameter given to the SIRIUS"
                           "method: '%s'\n", p->value);
                io_parameters[level_params] = NULL;
            }
            level_params++;
        } else if (!strcasecmp (p->name, "path")) {
            errno = 0;            
            io_paths[level_paths] = strdup (p->value);
            if (!errno) {
                log_debug ("path %d set to %s for SIRIUS method\n", level_paths, io_parameters[level_paths]);
            } else {
                log_error ("Invalid 'path' parameter given to the SIRIUS"
                           "method: '%s'\n", p->value);
                io_paths[level_paths] = NULL;
            }
            level_paths++;
        } else if (!strcasecmp (p->name, "thresh_type"))
        {
            if (!strcasecmp (p->value, "absolute"))
            {
                thresh_type = ABSOLUTE;
            }
            else if (!strcasecmp (p->value, "grad"))
            {
                thresh_type = GRAD;
            }
        } else if (!strcasecmp (p->name, "thresh"))
        {
            threshold = atof (p->value);
        } else if (!strcasecmp (p->name, "save-delta"))
        {
            save_delta = atoi (p->value);
        } else if (!strcasecmp (p->name, "compress-delta"))
        {
            compress_delta = atoi (p->value);

        } else if (!strcasecmp (p->name, "compression-tolerance"))
        {
            compr_tolerance = atof (p->value);
        } else if (!strcasecmp (p->name, "decimation-ratio"))
        {
            dec_ratio = atoi (p->value);
        } else {
            log_error ("Parameter name %s is not recognized by the SIRIUS "
                       "method\n", p->name);
        }

        p = p->next;
    }
    assert(nlevels==level_params);
    assert(nlevels==level_paths);
}


void adios_sirius_adaptive_init(const PairStruct * parameters,
                       struct adios_method_struct * method)
{
    struct adios_sa_data_struct * md = (struct adios_sa_data_struct *)
        method->method_data;

    method->method_data = malloc (sizeof (struct adios_sa_data_struct));
    md = (struct adios_sa_data_struct *) method->method_data;

    init_output_parameters(parameters);
}


static void init_method_parameters(struct adios_sa_data_struct * md)
{
    int l;
    for(l=0; l < nlevels; l++)
    {
        md->level[l].varcnt=0;
        md->level[l].vars=NULL;
        md->level[l].vars_head=NULL;
        md->level[l].fd = 0;
        md->level[l].filename = NULL;
        md->level[l].grp_name = NULL;
        md->level[l].totalsize = 0;
    }
}


int adios_sirius_adaptive_open (struct adios_file_struct * fd
                               ,struct adios_method_struct * method
                               ,MPI_Comm comm
                               )
{

    struct adios_sa_data_struct 
        * md = (struct adios_sa_data_struct *) method->method_data;
    char mode[2];
    int l;

    switch (fd->mode)
    {
        case adios_mode_read:
        {
            adios_error (err_invalid_file_mode
                        ,"SIRIUS_ADAPTIVE method: "
                         "Read mode is not supported.\n"
                        );
            return -1;
        }

        case adios_mode_append:
        case adios_mode_update:
        case adios_mode_write:
        {
            md->group_comm = comm;
            if (md->group_comm != MPI_COMM_NULL)
            {
                MPI_Comm_rank (md->group_comm, &md->rank);
                MPI_Comm_size (md->group_comm, &md->size);
            }

            fd->group->process_id = md->rank;

            init_method_parameters(md);

            define_iogroups(method);

            for (l = 0; l < nlevels; l++)
            {
                //check if the directory exists and create it if it doesn't
                struct stat sb;
                if((stat(io_paths[l], &sb) != 0) || !S_ISDIR(sb.st_mode))
                {
                    //directory doesn't exist
                    //FIXME: there is a case where something already exists but
                    //isn't a directory. Hard to imagine though so I am ignoring
                    //it for the time being
                    mkdir (io_paths[l], 0700);
                }
                md->level[l].filename = malloc (strlen(io_paths[l]) + strlen(fd->name) + 1);
                sprintf (md->level[l].filename, "%s/%s", io_paths[l], fd->name);
                convert_file_mode (fd->mode, mode);
               
                // Now call the transport
                common_adios_open (&(md->level[l].fd)
                                  ,md->level[l].grp_name
                                  ,md->level[l].filename
                                  ,mode
                                  ,comm
                                  );
            }

            break;
        }

        default:
        {
            adios_error (err_invalid_file_mode
                        ,"SIRIUS_ADAPTIVE method: "
                         "Unknown file mode requested: %d\n"
                        ,fd->mode
                        );

            return adios_flag_no;
        }
    }

    return 1;
}

enum BUFFERING_STRATEGY
adios_sirius_adaptive_should_buffer (struct adios_file_struct * fd
                                    ,struct adios_method_struct * method
                                    )
{
    //this method handles its own buffering
    return no_buffering;
}


//initial variable structure
static void init_var_struct (struct var_struct * var)
{
    var->name = NULL;
    var->path = NULL;
    var->type = adios_unknown;
    var->next = NULL;
    var->global_dimensions = (char *) calloc (128, sizeof(char));
    var->local_dimensions = (char *) calloc (128, sizeof(char));
    var->local_offsets = (char *) calloc (128, sizeof(char));
    var->size = 0;
}


static int do_write (int64_t fd_p, const char * name, void * var)
{
    struct adios_file_struct * fd = (struct adios_file_struct *) fd_p;

    if (!fd)
    {
        adios_error (err_invalid_file_pointer, "Invalid handle passed to adios_write\n");
        return 1;
    }

    struct adios_var_struct * v = fd->group->vars;
    struct adios_method_list_struct * m = fd->group->methods;

    if (m && m->next == NULL && m->method->m == ADIOS_METHOD_NULL)
    {
        // nothing to do so just return
        return 0;
    }

    v = adios_find_var_by_name (fd->group, name);

    if (!v)
    {
        adios_error (err_invalid_varname, "Bad var name (ignored) in SIRIUS adios_write(): '%s'\n", name);
        return 1;
    }

    common_adios_write_byid (fd, v, var);

    return 0;
}


static enum ADIOS_ERRCODES
alloc_var_struct (struct adios_sa_data_struct * md, int level)
{
    struct var_struct
        * var = (struct var_struct *) malloc (sizeof(struct var_struct));
    if (!var)
    {
        adios_error (err_no_memory, "No memory to allocate"
                    "yet another var in SIRIUS_ADAPTIVE method\n"
                    );

        return err_no_memory;
    }

    var->prev = md->level[level].vars;
    var->next = NULL;
    if (md->level[level].varcnt == 0)
    {
        //assign the header of the variable list
        md->level[level].vars_head = var;
    }

    md->level[level].vars = var;

    // initialize the variable structure
    init_var_struct (md->level[level].vars);

    return err_no_error;
}

static uint64_t get_var_dimensions (struct adios_var_struct * v, int ndims, uint64_t *gdims, uint64_t *ldims, uint64_t *offsets)
{
    struct adios_dimension_struct * d = v->dimensions;
    int dims_count = 0;
    uint64_t nelems = 1;
    while (d)
    {
        uint64_t dim = 0;
        //local dimension
        dim = adios_get_dim_value (&d->dimension);
        ldims[dims_count]=dim;

        //global dimension
        dim = adios_get_dim_value (&d->global_dimension);
        gdims[dims_count]=dim;

        //local offsets
        dim = adios_get_dim_value (&d->local_offset);
        offsets[dims_count]=dim;

        nelems *= ldims[dims_count];
        dims_count++;
        d=d->next;
    }
    return nelems;
}

static char * print_dimensions (int ndims, uint64_t *values)
{
    char * s = calloc (ndims*16, sizeof(char));
    int i = 0;
    for (i=0; i < ndims; i++)
    {
        if (i==0)
            sprintf(s, "%" PRIu64, values[i]);
        else
            sprintf(s, "%s,%" PRIu64, s, values[i]);
    }
    return s;
}

void get_coord (uint64_t element, int ndims, uint64_t * ldims, uint64_t * coord)
{
    int i;

    for (i = ndims - 1; i > -1; i--)
    {
        coord[i] = element % ldims[i];
        element = element / ldims[i];
    }
}

uint64_t get_linearized_index(int ndims, uint64_t * ldims, uint64_t * coord)
{
    int i;
    uint64_t index = 0;

    for (i = 0; i < ndims; i++)
    {
        index *= ldims[i];
        index = index + coord[i];
    }

    return index;
}

double get_value_by_coord (void * data, 
                           int ndims, 
                           uint64_t * ldims, 
                           uint64_t * coord)
{
    uint64_t idx = get_linearized_index (ndims, ldims, coord);
    return *((double *) data + idx);
}

int insert_node (double * newz, double * newr, double * newfield, int * size,
                  double z, double r, double field)
{
    int found;

    found = 0;
    for (int node = 0; node < *size; node++)
    {
        if (z == newz[node] && r == newr[node])
        {
            found = 1;
            return node;
        }
    }

    if (!found)
    {   
        newz[*size] = z;
        newr[*size] = r;
        newfield[*size] = field;

        (*size)++;

        return (*size) - 1;
    }
}

int find_mincost (int ** conn, edge_cost_t ** cost_matrix, 
                  int nvertices, double * r, double * z)
{
    int min_idx;
    double min_cost = DBL_MAX;

    for (int i = 0; i < nvertices; i++)
    {
        int j = 0;

        while (j < MAX_NODE_DEGREE && conn[i][j] != -1)
        {
            if (cost_matrix[i][j].cost <= min_cost)
            {
                min_cost = cost_matrix[i][j].cost;
                min_idx = i * MAX_NODE_DEGREE + j;
            }

            j++;
        }

    }

    return min_idx;
}


pqueue_t * build_pq (int ** conn, edge_cost_t *** cost_matrix_p,
                   int nvertices, double * r, double * z)
{
    edge_cost_t ** cost_matrix = (edge_cost_t **) malloc (nvertices * 8);
    assert (cost_matrix);

    * cost_matrix_p = cost_matrix;

    for (int i = 0; i < nvertices; i++)
    {
        cost_matrix[i] = (edge_cost_t *) malloc (MAX_NODE_DEGREE * sizeof (edge_cost_t));
        assert (cost_matrix[i]);

        for (int j = 0; j < MAX_NODE_DEGREE; j++)
        {
            cost_matrix[i][j].cost = DBL_MAX;
            cost_matrix[i][j].pq_node = 0;
        }
    }

    pqueue_t * pq = pqueue_init (nvertices * MAX_NODE_DEGREE, 
                                 cmp_pri, get_pri, set_pri, 
                                 get_pos, set_pos);

    for (int i = 0; i < nvertices; i++)
    {
        int n1 = i;
        int j = 0;

        while (j < MAX_NODE_DEGREE && conn[i][j] != -1)
        {
            int n2 = conn[i][j];
            cost_matrix[i][j].cost = sqrt(pow (r[n1] - r[n2], 2) + pow (z[n1] - z[n2], 2));

            node_t * pqn = (node_t *) malloc (sizeof (node_t));
            pqn->pri = cost_matrix[i][j].cost;
            pqn->val = i * MAX_NODE_DEGREE + j;

            cost_matrix[i][j].pq_node = pqn;

            pqueue_insert(pq, pqn);

            j++;
        }

    }

    return pq;
}

void sort2 (int * n1, int * n2)
{
    int t;

    if (* n1 > * n2)
    {
        t = * n1;
        * n1 = * n2;
        * n2 = t;
    }
}

void sort3 (int * n1, int * n2, int * n3)
{
    int t;

    if (* n1 > * n2)
    {
        t = * n1;
        * n1 = * n2;
        * n2 = t;
    }

    if (* n2 > * n3)
    {
        t = * n2;
        * n2 = * n3;
        * n3 = t;
    }

    if (* n1 > * n2)
    {
        t = * n1;
        * n1 = * n2;
        * n2 = t;
    }

}

void insert_hash (int nvertices, int node, int k, int v)
{
    assert (node >= 0 && node < nvertices);

    int * pkey = malloc (4);
    * pkey = k;

    int * pval = malloc (4);
    * pval = v;

    assert (TRUE == g_hash_table_insert (nodes_ght[node], pkey, pval));
}

void replace_hash (int nvertices, int node, int k, int v)
{
    assert (node >= 0 && node < nvertices);

    int * pkey = malloc (4);
    * pkey = k;

    int * pval = malloc (4);
    * pval = v;

    if (FALSE == g_hash_table_replace (nodes_ght[node], pkey, pval));
}

void remove_hash (int nvertices, int node, int k)
{
    assert (node >= 0 && node < nvertices);

    assert (TRUE == g_hash_table_remove (nodes_ght[node], &k));
}

void insert_triangle (int nvertices, int n1, int n2, int n3, int ** conn)
{
    int j = 0, found;

    found = 0;
    j = 0;
    while (j < MAX_NODE_DEGREE && conn[n1][j] != -1)
    {
        if (conn[n1][j] == n2)
        {
            found = 1;
        }

        j++;
    }

    if (j == MAX_NODE_DEGREE)
    {
        printf ("Reaching max MAX_NODE_DEGREE.\n");
    }
    else
    {
        if (!found)
        {
           conn[n1][j] = n2;
           insert_hash (nvertices, n2, n1, n1 * MAX_NODE_DEGREE + j);
        }
    }

    found = 0;
    j = 0;
    while (j < MAX_NODE_DEGREE && conn[n1][j] != -1)
    {
        if (conn[n1][j] == n3)
        {
            found = 1;
        }
        
        j++;
    }

    if (j == MAX_NODE_DEGREE)
    {
        printf ("Reaching max MAX_NODE_DEGREE.\n");
    }
    else
    {
        if (!found)
        {
           conn[n1][j] = n3;
           insert_hash (nvertices, n3, n1, n1 * MAX_NODE_DEGREE + j);
        }
    }

    found = 0;
    j = 0;
    while (j < MAX_NODE_DEGREE && conn[n2][j] != -1)
    {
        if (conn[n2][j] == n3)
        {
            found = 1;
        }

        j++;
    }

    if (j == MAX_NODE_DEGREE)
    {
        printf ("Reaching max MAX_NODE_DEGREE.\n");
    }
    else
    {
        if (!found)
        {
           conn[n2][j] = n3;
           insert_hash (nvertices, n3, n2, n2 * MAX_NODE_DEGREE + j);
        }
    }

}


int new_node(int ** conn, int v1, int n)
{
    int i = 0;

    while (i < MAX_NODE_DEGREE && conn[v1][i] != -1)
    {
        if (conn[v1][i] == n)
            return 0;

        i++;
    }

    return 1;
}

int ** build_conn (int nvertices, int * mesh, int nmesh)
{
    int ** conn = (int **) malloc (nvertices * 8);
    nodes_ght = (GHashTable **) malloc (nvertices * 8); 

    for (int i = 0; i < nvertices; i++)
    {
        conn[i] = (int *) malloc (MAX_NODE_DEGREE * 4);
        for (int j = 0; j < MAX_NODE_DEGREE; j++)
        {
            conn[i][j] = -1;
        }

        nodes_ght[i] = g_hash_table_new_full (g_int_hash, g_int_equal, free, free);
    }

    for (int i = 0; i < nmesh; i++)
    {
        int n1 = * (mesh + i * 3); 
        int n2 = * (mesh + i * 3 + 1);
        int n3 = * (mesh + i * 3 + 2);

        assert (n1 < nvertices && n2 < nvertices || n3 < nvertices);
   
        sort3 (&n1, &n2, &n3);
        insert_triangle (nvertices, n1, n2, n3, conn);
    }

    return conn;
}

int intersect (int ** conn, int n1, int n2, int * n3_list)
{
    int i, j, c;
    assert (n3_list);

    i = 0;
    c = 0;
    while (i < MAX_NODE_DEGREE && conn[n1][i] != -1)
    {
        j = 0;
        while (j < MAX_NODE_DEGREE && conn[n2][j] != -1)
        {
            if (conn[n1][i] == conn[n2][j])
            {
                n3_list[c++] = conn[n1][i];
            }

            j++;
        }

        i++;
    }

    return c;
}

void prep_mesh (int ** conn, int nvertices, int nvertices_new)
                
{
    for (int i = 0; i < nvertices; i++)
    {
        int j = 0; 
        while (j < MAX_NODE_DEGREE && conn[i][j] != -1)
        {
            int n1 = i, n2 = conn[i][j];
            int k = 0;

            while (k < MAX_NODE_DEGREE && conn[n2][k] != -1)
            {
                if (conn[n2][k] == n1) break;

                k++;
            }

            if (k < MAX_NODE_DEGREE && conn[n2][k] == -1)
            {
                conn[n2][k] = n1;
            }

            j++;
        }
    }
}

int to_offset (int * nodes_cut, int nnodes_cut, int my_node_id)
{
    int i;

    for (i = 0; i < nnodes_cut; i++)
    {
        if (my_node_id < nodes_cut[i])
        {
            return i;
        }
    }

    return nnodes_cut;
}

void rebuild_conn (int ** conn, int nvertices, int nvertices_new, 
                   int * nodes_cut)
{
    for (int i = 0; i < nvertices; i++)
    {
        int j = 0;

        while (j < MAX_NODE_DEGREE && conn[i][j] != -1)
        {
            conn[i][j] -= to_offset (nodes_cut, nvertices - nvertices_new,
                                     conn[i][j]);
            j++;
        }
    }
}
#if 0
uint64_t * alloc_ht (int size)
{
    uint64_t * ht = (uint64_t *) malloc (size * 8);

    assert (ht);

    for (int i = 0; i < size; i++)
    {
        ht[i] = 0;
    }

    return ht;
}


int ht_insert (int n1, int n2, int n3, uint64_t * ht, int unit, int size)
{
    assert (n1 < n2 && n2 < n3);

    uint64_t v = n1 * unit * unit + n2 * unit + n3;
    int i = 0;

    while (i < size && ht[i] != 0)
    {
        if (ht[i++] == v) return 0;
    }

//    printf ("i = %d, size = %d\n", i, size);
    assert (i < size * 10);

    ht[i] = v;

    return 1;
}
#endif

int build_mesh (int ** conn, int nvertices, int nvertices_new,
                int nmesh, int * nodes_cut, int ** mesh_new)
{
//    rebuild_conn (conn, nvertices, nvertices_new, nodes_cut);

    int * mesh = malloc (nmesh * 3 * 4);
    assert (mesh);

    int * n3_list = 0, len = 0, lastcell = 0;

#define MAX_COMMON_NODES 50
    n3_list = malloc (MAX_COMMON_NODES * 4);
    assert (n3_list);

    GHashTable * ght = g_hash_table_new_full (g_str_hash,g_str_equal, free, free);

    for (int i = 0; i < nvertices; i++)
    {
        int j = 0, k = 0;

        while (j < MAX_NODE_DEGREE && conn[i][j] != -1)
        {
            int n1 = i, n2 = conn[i][j];

            len = intersect (conn, n1, n2, n3_list);
            if (len > 0)
            {
                for (k = 0; k < len; k++)
                {
                    int n3 = n3_list[k];
                    char temp_str[128];
                    int t1 = n1, t2 = n2, t3 = n3;

                    sort3 (&t1, &t2, &t3);

                    assert (t1 < t2 && t2 < t3);

                    sprintf (temp_str, "%d,%d,%d", t1, t2, t3);
                    gchar * key_str = g_strdup (temp_str);

                    if (g_hash_table_insert (ght, key_str, 0) == TRUE)
                    {
                        * (mesh + lastcell * 3) = n1 - to_offset (nodes_cut, nvertices - nvertices_new, n1);
;
                        * (mesh + lastcell * 3 + 1) = n2 - to_offset (nodes_cut, nvertices - nvertices_new, n2);

                        * (mesh + lastcell * 3 + 2) = n3 - to_offset (nodes_cut, nvertices - nvertices_new, n3);
; 
                        lastcell++; 
                    }
                }
            }

            j++;
        }
    }

    g_hash_table_destroy (ght);

    free (n3_list);
    n3_list = 0;

    * mesh_new = mesh;

    return lastcell;
}

int * build_nodes_cut_list (int ** conn, int nvertices, int nvertices_new)
{
    int * nodes_cut = malloc ((nvertices - nvertices_new) * 4);
    assert (nodes_cut);

    int next = 0;

    for (int i = 0; i < nvertices; i++)
    {
        if (conn[i][0] == -1)
        {
            nodes_cut[next++] = i;
        }
    }

    assert (nvertices - nvertices_new == next);

    return nodes_cut;
}

void free_nodes_cut_list (int * nodes_cut)
{
    assert (nodes_cut);
    free (nodes_cut);
}

void build_field (int ** conn, int nvertices, int nvertices_new, int * nodes_cut,
                  double * r, double * z, double * field,
                  double * r_new, double * z_new, double * field_new)
{
    int i, prev, off;

    prev = 0;
    off = 0;
    for (i = 0; i < nvertices - nvertices_new + 1; i++)
    {
        int elems_to_cp = (i < nvertices - nvertices_new ? 
                           nodes_cut[i] - prev : nvertices - prev);
        
        memcpy (r_new, r + off, elems_to_cp * 8); 
        memcpy (z_new, z + off, elems_to_cp * 8); 
        memcpy (field_new, field + off, elems_to_cp * 8);

        r_new += elems_to_cp;
        z_new += elems_to_cp;
        field_new += elems_to_cp;

        off += elems_to_cp + 1;
        prev = nodes_cut[i] + 1;
    }
}

int get_node_degree (int ** conn, int n)
{
    int i = 0;

    while (i < MAX_NODE_DEGREE && conn[n][i] != -1)
    {
        i++;
    }

    return i;
}

void update_field (int v1, int v2, double * r, double * z, double * field)
{
    r[v1] = (r[v1] + r[v2]) / 2;
    z[v1] = (z[v1] + z[v2]) / 2;
    field[v1] = (field[v1] + field[v2]) / 2;
}

void update_cost (int ** conn, edge_cost_t ** cost_matrix,
                  pqueue_t * pq, int n1, int j,
                  double * r, double * z)
{
    int n2 = conn[n1][j];

    if (n2 == -1)
    {
        cost_matrix[n1][j].cost = DBL_MAX;

        if (cost_matrix[n1][j].pq_node)
        {
            // to be on the safe side, first change the priority.
            // This is a fix for a bug from pqueue.
            pqueue_change_priority (pq, cost_matrix[n1][j].cost,
                                    cost_matrix[n1][j].pq_node);

            assert (!pqueue_remove(pq, cost_matrix[n1][j].pq_node));
            free (cost_matrix[n1][j].pq_node);
            cost_matrix[n1][j].pq_node = 0;
        }
    }
    else
    {
        cost_matrix[n1][j].cost = sqrt(pow (r[n1] - r[n2], 2) + pow (z[n1] - z[n2], 2));

        if (cost_matrix[n1][j].pq_node)
        {
            pqueue_change_priority (pq, cost_matrix[n1][j].cost, 
                                    cost_matrix[n1][j].pq_node);

        }
        else
        {
            node_t * pqn = (node_t *) malloc (sizeof (node_t));
            pqn->pri = cost_matrix[n1][j].cost;
            pqn->val = n1 * MAX_NODE_DEGREE + j;
            cost_matrix[n1][j].pq_node = pqn;
            
            pqueue_insert(pq, pqn);
        }

    }
}

void free_conn (int ** conn, int nvertices)
{
    for (int i = 0; i < nvertices; i++)
    {
        free (conn[i]);
    }

    free (conn);
}

void free_cost_matrix (edge_cost_t ** cost_matrix, int nvertices)
{
    for (int i = 0; i < nvertices; i++)
    {
        free (cost_matrix[i]);
    }

    free (cost_matrix);
}

double sign (double p1x, double p1y,
             double p2x, double p2y,
             double p3x, double p3y)
{
    return (p1x - p3x) * (p2y - p3y) - (p2x - p3x) * (p1y - p3y);
}

int PointInTriangle (double px, double py,
                     double v1x, double v1y,
                     double v2x, double v2y,
                     double v3x, double v3y)
{
    int b1, b2, b3;

    b1 = sign(px, py, v1x, v1y, v2x, v2y) <= 0.0f;
    b2 = sign(px, py, v2x, v2y, v3x, v3y) <= 0.0f;
    b3 = sign(px, py, v3x, v3y, v1x, v1y) <= 0.0f;

    return ((b1 == b2) && (b2 == b3));
}

void find_clst (double * r, double * z, int node,
                double * r_reduced, double * z_reduced,
                int n1, int n2, int n3,
                int * n_clst)
{
    double d1 = sqrt (pow (r[node] - r_reduced[n1], 2)
                    + pow (z[node] - z_reduced[n1], 2));
    double d2 = sqrt (pow (r[node] - r_reduced[n2], 2)
                    + pow (z[node] - z_reduced[n2], 2));
    double d3 = sqrt (pow (r[node] - r_reduced[n3], 2)
                    + pow (z[node] - z_reduced[n3], 2));

    if (d1 <= d2)
    {
        if (d1 <= d3)
        {
            * n_clst = n1;
        }
        else
        {
            * n_clst = n3;
        }
    }
    else
    {
        if (d2 <= d3)
        {
            * n_clst = n2;
        }
        else
        {
            * n_clst = n3;
        }
    }
}

void test_delta (double * r, double * z, double * field,
                int nvertices, int * mesh, int nmesh,
                double * r_reduced, double * z_reduced,
                double * field_reduced, int nvertices_new,
                int * mesh_reduced, int nmesh_new,
                double * field_delta, double ** pfield_full
               )
{

    double * field_full = (double *) malloc (nvertices * 8);
    assert (field_full);

    for (int i = 0; i < nvertices; i++)
    {
        for (int m = 0; m < nmesh_new; m++)
        {
            int n1 = * (mesh_reduced + m * 3);
            int n2 = * (mesh_reduced + m * 3 + 1);
            int n3 = * (mesh_reduced + m * 3 + 2);

            int in = PointInTriangle (r[i], z[i],
                                      r_reduced[n1], z_reduced[n1],
                                      r_reduced[n2], z_reduced[n2],
                                      r_reduced[n3], z_reduced[n3]);
#if 0
            find_clst (r, z, i
                       r_reduced, z_reduced,
                       n1, n2, n3,
                       &n_clst);
#endif
            if (in)
            {
                double estimate = (field_reduced[n1] + field_reduced[n2] + field_reduced[n3]) / 3.0;
                field_full[i] = estimate + field_delta[i];
                break;
            }
            else if (m == nmesh_new - 1)
            {
//                double estimate  = field_reduce[n_clst];
                field_full[i] = 0.0;
            }
        }
    }

    * pfield_full = field_full;
}

int partition (double * x, int * idx, int l, int r)
{
    double t;
    double pivot = x[l];
    int i = l; 
    int j = r + 1;
    int t_id;
		
    while(1)
    {
        do ++i; while (x[i] <= pivot && i <= r);
   	do --j; while (x[j] > pivot);
   	if(i >= j) break;
   	t = x[i]; x[i] = x[j]; x[j] = t;
        t_id = idx[i]; idx[i] = idx[j]; idx[j] = t_id;
    }

    t = x[l]; x[l] = x[j]; x[j] = t;
    t_id = idx[l]; idx[l] = idx[j]; idx[j] = t_id;

    return j;
}

void quickSort (double * a, int * idx, int l, int r)
{
    int j;

    if (l < r)
    {
        j = partition (a, idx, l, r);
        quickSort (a, idx, l, j-1);
        quickSort (a, idx, j+1, r);
    }
}

void sort_by_x1 (double * x, int nvertices,
                 double ** px_sorted, int ** px_idx)
{
    double * x_sorted = (double *) malloc (nvertices * 8);
    assert (x_sorted);

    memcpy (x_sorted, x, nvertices * 8);

    int * x_idx = (int *) malloc (nvertices * 4);
    assert (x_idx);

    for (int i = 0; i < nvertices; i++)
    {
        x_idx[i] = i;
    }

    quickSort (x_sorted, x_idx, 0, nvertices - 1);

    * px_sorted = x_sorted;
    * px_idx = x_idx;
}

void sort_by_x (double * x, int nvertices,
                double ** px_sorted, int ** px_idx)
{

    double * x_sorted = (double *) malloc (nvertices * 8);
    assert (x_sorted);

    memcpy (x_sorted, x, nvertices * 8);

    int * x_idx = (int *) malloc (nvertices * 4);
    assert (x_idx);

    for (int i = 0; i < nvertices; i++)
    {
        x_idx[i] = i;
    }

    for (int i = 0; i < nvertices - 1; i++)
    {
        for (int j = 0; j < nvertices - i - 1; j++)
        {
            if (x_sorted[j] > x_sorted[j + 1])
            {
                double temp_v = x_sorted[j];
                int temp_idx = x_idx[j];

                x_sorted[j] = x_sorted[j + 1];
                x_idx[j] = x_idx[j + 1];

                x_sorted[j + 1] = temp_v;
                x_idx[j + 1] = temp_idx;
            }
        } 
    }

    * px_sorted = x_sorted;
    * px_idx = x_idx;
}

void sort_by_rz (double * r, double * z, int nvertices,
                 double ** pr_sorted, double ** pz_sorted,
                 int ** pr_idx, int ** pz_idx)
{
double time1 = MPI_Wtime ();
    sort_by_x1 (r, nvertices,
               pr_sorted, pr_idx
              );

double time2 = MPI_Wtime ();
    sort_by_x1 (z, nvertices,
               pz_sorted, pz_idx
              );
double time3 = MPI_Wtime ();
//printf ("sort time: %f, %f\n", time2 - time1, time3 - time2);
}

double min (double x, double y, double z)
{
    if (x <= y)
    {
        if (x <= z) 
            return x;
        else
            return z;
    }
    else
    {
        if (y <= z)
            return y;
         else
            return z;
    }
}

double max (double x, double y, double z)
{
    if (x <= y)
    {
        if (y <= z)
            return z;
        else
            return y;
    }
    else
    {
        if (x <= z)
            return z;
         else
            return x;
    }
}

int find_low (double * sorted, int nvertices, double min)
{
    int i = 0;

    while (i < nvertices)
    {
        if (min <= sorted[i]) break;
        i++;
    }

    assert (i < nvertices);

    return i;
}

int find_low1 (double * sorted, int nvertices, double min)
{
    int low = 0, high = nvertices - 1, mid;
   
    while (1)
    {
        mid = (low + high) / 2;

        if (mid == 0) return 0;

        if (sorted[mid] == min) return mid;
        if (sorted[mid - 1] == min) return mid - 1;

        if (sorted[mid] > min && sorted[mid - 1] < min) return mid;

        if (sorted[mid] > min && sorted[mid - 1] > min)
        {
            high = mid;
        }
        else if (sorted[mid] < min && sorted[mid - 1] < min)
        {
            low = mid;
        }
    }
}

int find_high (double * sorted, int nvertices, double max)
{
    int i = nvertices - 1;

    while (i >= 0)
    {
        if (max >= sorted[i]) break;
        i--;
    }

    assert (i >= 0);

    return i;
}

int find_high1 (double * sorted, int nvertices, double max)
{
    int low = 0, high = nvertices - 1, mid;

    while (1)
    {
        mid = (low + high) / 2;

        if (sorted[mid] == max) return mid;
        if (sorted[mid + 1] == max) return mid + 1;

        if (sorted[mid] < max && sorted[mid + 1] > max) return mid;

        if (sorted[mid] < max && sorted[mid + 1] < max)
        {
            low = mid;
        }
        else if (sorted[mid] > max && sorted[mid + 1] > max)
        {
            high = mid;
        }
    }
}

int find_possible_nodes (double min_r, double max_r, 
                         double min_z, double max_z,
                         double * r_sorted, double * z_sorted,
                         int nvertices, int * r_idx, int * z_idx, 
                         int ** plist)
{
    int r_low, r_high, z_low, z_high;

    r_low = find_low1 (r_sorted, nvertices, min_r);
    r_high = find_high1 (r_sorted, nvertices, max_r);
    z_low = find_low1 (z_sorted, nvertices, min_z);
    z_high = find_high1 (z_sorted, nvertices, max_z);

    int * list = (int *) malloc ( (r_high - r_low + 1) * 4);
    int next_node = 0;

    for (int i = r_low; i <= r_high; i++)
    {
        int node_id = r_idx[i];

        int found = 0;
        for (int j = z_low; j <= z_high; j++)
        {
            if (z_idx[j] == node_id)
            {
                found = 1;
                break;
            }
        }

        if (found) list[next_node++] = node_id;
    }

    * plist = list;

    return next_node;
}

void get_delta1 (double * r, double * z, double * field,
                 int nvertices, int * mesh, int nmesh,
                 double * r_reduced, double * z_reduced,
                 double * field_reduced, int nvertices_new,
                 int * mesh_reduced, int nmesh_new,
                 double ** pfield_delta
                )
{
    double * delta = (double *) malloc (nvertices * 8);
    assert (delta);
    for (int i = 0; i < nvertices; i++)
    {
        delta[i] = 0.0;
    }

    double * r_sorted = 0, * z_sorted = 0;
    int * r_idx = 0, * z_idx = 0;

double time1 = MPI_Wtime();
    sort_by_rz (r, z, nvertices,
                &r_sorted, &z_sorted,
                &r_idx, &z_idx);

double time2 = MPI_Wtime();
#pragma omp parallel for
    for (int m = 0; m < nmesh_new; m++)
    {
        int n1 = * (mesh_reduced + m * 3);
        int n2 = * (mesh_reduced + m * 3 + 1);
        int n3 = * (mesh_reduced + m * 3 + 2);

        double min_r = min (r_reduced[n1], r_reduced[n2], r_reduced[n3]);
        double min_z = min (z_reduced[n1], z_reduced[n2], z_reduced[n3]);
        double max_r = max (r_reduced[n1], r_reduced[n2], r_reduced[n3]);
        double max_z = max (z_reduced[n1], z_reduced[n2], z_reduced[n3]);

        int * plist = 0;
double time_find0 = MPI_Wtime ();
        int nlist = find_possible_nodes (min_r, max_r, min_z, max_z,
                                         r_sorted, z_sorted, nvertices, 
                                         r_idx, z_idx, &plist);
double time_find1 = MPI_Wtime ();
        for (int i = 0; i < nlist; i++) 
        {
            int in = PointInTriangle (r[plist[i]], z[plist[i]],
                                      r_reduced[n1], z_reduced[n1],
                                      r_reduced[n2], z_reduced[n2],
                                      r_reduced[n3], z_reduced[n3]);
            if (in)
            {
                double estimate = (field_reduced[n1] + field_reduced[n2] + field_reduced[n3]) / 3.0;
                delta[plist[i]] = field[plist[i]] - estimate;
            }
//          double estimate  = field_reduce[n_clst];
//            delta[plist[i]] = 0.0;
        }
    }

double time3 = MPI_Wtime ();
printf ("get delta time = %f\n", time3 - time1);

    * pfield_delta = delta;
}

void get_delta (double * r, double * z, double * field,
                int nvertices, int * mesh, int nmesh,
                double * r_reduced, double * z_reduced,
                double * field_reduced, int nvertices_new,
                int * mesh_reduced, int nmesh_new,
                double ** pfield_delta
               )
{
    double start_time = MPI_Wtime ();
    double * delta = (double *) malloc (nvertices * 8);
    assert (delta);
    int tid;

#pragma omp parallel for
    for (int i = 0; i < nvertices; i++)
    {
//        tid = omp_get_thread_num();
//        printf("from thread = %d\n", tid);

        delta[i] = - DBL_MAX;

        for (int m = 0; m < nmesh_new; m++)
        {
            int n1 = * (mesh_reduced + m * 3);
            int n2 = * (mesh_reduced + m * 3 + 1);
            int n3 = * (mesh_reduced + m * 3 + 2);

            int in = PointInTriangle (r[i], z[i],
                                      r_reduced[n1], z_reduced[n1],
                                      r_reduced[n2], z_reduced[n2],
                                      r_reduced[n3], z_reduced[n3]);

            if (in)
            {
                double estimate = (field_reduced[n1] + field_reduced[n2]
                                 + field_reduced[n3]) / 3.0;
                delta[i] = field[i] - estimate;
                break;
            }
            else if (m == nmesh_new - 1)
            {
                delta[i] = 0.0;
            }
        }
    }

    * pfield_delta = delta;
    double end_time = MPI_Wtime ();
    printf ("get delta time = %f\n", end_time - start_time);
}

double calc_area (double * r, double * z, double * data,
                  int nvertices, int * mesh, int nmesh
                 )
{
    int ntaggedCells = 0;
    double surface_size = 0.0;

    for (int m = 0; m < nmesh; m++)
    {
        int n1 = * (mesh + m * 3);
        int n2 = * (mesh + m * 3 + 1);
        int n3 = * (mesh + m * 3 + 2);

        double avg_mag = (data[n1] + data[n2] + data[n3]) / 3;

        if (avg_mag > threshold)
        {
            ntaggedCells++;

            double a = sqrt (pow (r[n1] - r[n2], 2) + pow (z[n1] - z[n2], 2));
            double b = sqrt (pow (r[n1] - r[n3], 2) + pow (z[n1] - z[n3], 2));
            double c = sqrt (pow (r[n2] - r[n3], 2) + pow (z[n2] - z[n3], 2));
            double p = (a + b + c) / 2;

            surface_size += sqrt (p * (p - a) * (p - b) * (p - c));
        }
    }  // loop through the node connectivity array

    //printf ("ntaggedCells = %d, new surface = %f\n", ntaggedCells, surface_size);

}

void decimate (double * rorg, double * zorg, double * fieldorg, 
               int nvertices, int * mesh, int nmesh,
               double ** r_reduced, double ** z_reduced, 
               double ** field_reduced, int * nvertices_new,
               int ** mesh_reduced, int * nmesh_new
              )
{
    double * r, * z, * field;
    double * r_new, * z_new, * field_new;
    int * mesh_new;
    int vertices_cut = 0, min_idx, pq_v1;
    edge_cost_t ** cost_matrix;

    r = (double *) malloc (nvertices * 8);
    z = (double *) malloc (nvertices * 8);
    field = (double *) malloc (nvertices * 8);
    assert (r && z && field);

    memcpy (r, rorg, nvertices * 8);
    memcpy (z, zorg, nvertices * 8);
    memcpy (field, fieldorg, nvertices * 8);
    
    int ** conn = build_conn (nvertices, mesh, nmesh);
    pqueue_t * pq = build_pq (conn, &cost_matrix, nvertices, r, z);
#if 0
    int min_idx = find_mincost (conn, cost_matrix, nvertices, r, z);
#endif
    node_t * pq_min = pqueue_pop (pq);
    min_idx = pq_min->val;
    pq_v1 = min_idx / MAX_NODE_DEGREE;

    assert (pq_v1 >=0 && pq_v1 < nvertices);

    free (cost_matrix[pq_v1][min_idx % MAX_NODE_DEGREE].pq_node);
    cost_matrix[pq_v1][min_idx % MAX_NODE_DEGREE].pq_node = 0;

double t0 = MPI_Wtime();
    while ((double)vertices_cut / (double)nvertices < (1.0 - 1.0 / dec_ratio))
//    while (vertices_cut < 10000)
    {
        int v1 = min_idx / MAX_NODE_DEGREE;
        assert (v1 >=0 && v1 < nvertices);

        int v2 = conn[v1][min_idx % MAX_NODE_DEGREE];
        assert (v2 >=0 && v2 < nvertices);

        sort2 (&v1, &v2);

        update_field (v1, v2, r, z, field);

        int i = 0, j = 0, m = 0, k = 0;
        while (j < MAX_NODE_DEGREE && conn[v1][j] != -1 
            && conn[v1][j] != v2)
        {
            update_cost (conn, cost_matrix, pq, v1, j, r, z);
            j++;
        }

        assert (j < MAX_NODE_DEGREE && conn[v1][j] == v2);

        // 1. To remove (v1, v2)
        // 2. keep the node ID of v1, remove node ID of v2
        // 
        // 3. change v1's position, and value.
        if (j == MAX_NODE_DEGREE - 1)
        {
            conn[v1][j] = -1;
            assert (0);
        }
        else
        {
            remove_hash (nvertices, v2, v1);

            while (j < MAX_NODE_DEGREE - 1 && conn[v1][j] != -1)
            {
                conn[v1][j] = conn[v1][j + 1];

                update_cost (conn, cost_matrix, pq, v1, j, r, z);

                if (conn[v1][j] != -1)
                {
                    replace_hash (nvertices, conn[v1][j], v1, v1 * MAX_NODE_DEGREE + j);
                    j++;
                }
            }

            assert (j < MAX_NODE_DEGREE - 1);
        }

        i = 0;
        m = 0;
        while ((j + m) < MAX_NODE_DEGREE && conn[v2][i] != -1)
        {
            if (new_node(conn, v1, conn[v2][i]))
            {
                conn[v1][j + m] = conn[v2][i];

                update_cost (conn, cost_matrix, pq, v1, j + m, r, z);
                insert_hash (nvertices, conn[v1][j + m], v1, v1 * MAX_NODE_DEGREE + j + m);
                m++;
            }

            remove_hash (nvertices, conn[v2][i], v2);
            conn[v2][i] = -1;
            update_cost (conn, cost_matrix, pq, v2, i, r, z);

            i++;
        }

        assert (j + m < MAX_NODE_DEGREE);
#if 1
        GHashTableIter iter;
        int * key, * value;

        g_hash_table_iter_init (&iter, nodes_ght[v1]);
        while (g_hash_table_iter_next (&iter, (gpointer *)&key, (gpointer *)&value))
        {
            i = * key;
            j = (* value) % MAX_NODE_DEGREE;

            update_cost (conn, cost_matrix, pq, i, j, r, z);
        }

        g_hash_table_iter_init (&iter, nodes_ght[v2]);

        int v2_ht_size = g_hash_table_size (nodes_ght[v2]);
        int toff;
        int * temp_table = malloc (v2_ht_size * 2 * 4);

        assert  (temp_table);

        toff = 0;
        while (g_hash_table_iter_next (&iter, (gpointer *)&key, (gpointer *)&value))
        {
            * (temp_table + toff * 2) = * key;
            * (temp_table + toff * 2 + 1) = * value;

            toff++;
        }

        toff = 0;
        while (toff < v2_ht_size)
        {
            i = * (temp_table + toff * 2);
            j = * (temp_table + toff * 2 + 1) % MAX_NODE_DEGREE;

            toff++;

            if (i < v1)
            {
                if (new_node(conn, i, v1))
                {
                    conn[i][j] = v1;
                    update_cost (conn, cost_matrix, pq, i, j, r, z);

                    insert_hash (nvertices, conn[i][j], i, i * MAX_NODE_DEGREE + j);
                 }
                 else
                 {
                     k = j;

                     remove_hash (nvertices, conn[i][k], i);

                     while (k < MAX_NODE_DEGREE - 1 && conn[i][k] != -1)
                     {
                         conn[i][k] = conn[i][k + 1];

                         update_cost (conn, cost_matrix, pq, i, k, r, z);

                         if (conn[i][k] != -1)
                         {
                             replace_hash (nvertices, conn[i][k], i, i * MAX_NODE_DEGREE + k);
                             k++;
                         }
                     }

                 }
            }
            else if (i > v1)
            {
                k = 0;
                while (k < MAX_NODE_DEGREE && conn[v1][k] != -1)
                {
                    k++;
                }

                if (k < MAX_NODE_DEGREE)
                {
                    if (new_node(conn, v1, i))
                    {
                        conn[v1][k] = i;

                        update_cost (conn, cost_matrix, pq, v1, k, r, z);
                        insert_hash (nvertices, conn[v1][k], v1, v1 * MAX_NODE_DEGREE + k);
                    }

                    k = j;

                    remove_hash (nvertices, conn[i][k], i);
                    while (k < MAX_NODE_DEGREE - 1 && conn[i][k] != -1)
                    {
                        conn[i][k] = conn[i][k + 1];

                        update_cost (conn, cost_matrix, pq, i, k, r, z);

                        if (conn[i][k] != -1)
                        {
                            replace_hash (nvertices, conn[i][k], i, i * MAX_NODE_DEGREE + k);

                            k++;
                        }
                    }

                    if (k == MAX_NODE_DEGREE - 1)
                    {
                        remove_hash (nvertices, conn[i][k], i);

                        conn[i][k] = -1;

                        update_cost (conn, cost_matrix, pq, i, k, r, z);
                    }
                }
            }
        }

        free (temp_table);

#endif
#if 0
        for (i = 0; i < v2; i++)
        {
            j = 0;

            if (i < v1)
            {
                while (j < MAX_NODE_DEGREE && conn[i][j] != -1)
                {
                    if (conn[i][j] == v1)
                    {
int * ptmp = g_hash_table_lookup (nodes_ght[v1], &i);
assert (* ptmp == (i * MAX_NODE_DEGREE + j));

                        update_cost (conn, cost_matrix, pq, i, j, r, z);
                        j++;
                    }
                    else if (conn[i][j] == v2)
                    {
int * ptmp = g_hash_table_lookup (nodes_ght[v2], &i);
assert (* ptmp == (i * MAX_NODE_DEGREE + j));

                        if (new_node(conn, i, v1))
                        {
                            conn[i][j] = v1;
                            update_cost (conn, cost_matrix, pq, i, j, r, z);

                            insert_hash (nvertices, conn[i][j], i, i * MAX_NODE_DEGREE + j);
                            j++;
                        }
                        else
                        {
                            k = j;

                            remove_hash (nvertices, conn[i][k], i);

                            while (k < MAX_NODE_DEGREE - 1 && conn[i][k] != -1)
                            {
                                conn[i][k] = conn[i][k + 1];

                                update_cost (conn, cost_matrix, pq, i, k, r, z);

                                if (conn[i][k] != -1)
                                {
                                    replace_hash (nvertices, conn[i][k], i, i * MAX_NODE_DEGREE + k);
                                    k++;
                                }
                            }
                        }
                    }
                    else
                    {
                        j++;
                    }
                }
            }
            else if (i > v1)
            {
                while (j < MAX_NODE_DEGREE && conn[i][j] != -1)
                {
                    if (conn[i][j] == v2)
                    {
int * ptmp = g_hash_table_lookup (nodes_ght[v2], &i);
assert (* ptmp == (i * MAX_NODE_DEGREE + j));

                        k = 0;
                        while (k < MAX_NODE_DEGREE && conn[v1][k] != -1)
                        {
                            k++;
                        }

                        if (k < MAX_NODE_DEGREE)
                        {
                            if (new_node(conn, v1, i))
                            {
                                conn[v1][k] = i;

                                update_cost (conn, cost_matrix, pq, v1, k, r, z);
                                insert_hash (nvertices, conn[v1][k], v1, v1 * MAX_NODE_DEGREE + k);
                            }

                            k = j;

                            remove_hash (nvertices, conn[i][k], i);

                            while (k < MAX_NODE_DEGREE - 1 && conn[i][k] != -1)
                            {
                                conn[i][k] = conn[i][k + 1];

                                update_cost (conn, cost_matrix, pq, i, k, r, z);

                                if (conn[i][k] != -1)
                                {
                                    replace_hash (nvertices, conn[i][k], i, i * MAX_NODE_DEGREE + k);

                                    k++;
                                }
                            }

                            if (k == MAX_NODE_DEGREE - 1)
                            {
                                remove_hash (nvertices, conn[i][k], i);

                                conn[i][k] = -1;

                                update_cost (conn, cost_matrix, pq, i, k, r, z);
                            }
                        }

                    }

                    j++;
                }
            }
        }
#endif

        vertices_cut++;
double t1 = MPI_Wtime();
#if 0
        min_idx = find_mincost (conn, cost_matrix, nvertices, r, z);
#endif
        node_t * pq_min = pqueue_pop(pq);
        min_idx = pq_min->val;

        int pq_v1 = min_idx / MAX_NODE_DEGREE;
        assert (pq_v1 >=0 && pq_v1 < nvertices);

        free (cost_matrix[pq_v1][min_idx % MAX_NODE_DEGREE].pq_node);
        cost_matrix[pq_v1][min_idx % MAX_NODE_DEGREE].pq_node = 0;
    }

double t2 = MPI_Wtime();
//printf ("time = %f\n", t2 - t0);
#if 0
    for (int i = 0; i < nvertices; i++)
    {
        printf ("(%d): ", i);
        for (int j = 0; j < 20; j++)
        {
            printf ("%d ", conn[i][j]);
        }
        printf ("\n");
    }
#endif

    * nvertices_new = nvertices - vertices_cut;
    //printf ("nvertices_old = %d, nvertices_new = %d\n", nvertices, * nvertices_new);

    prep_mesh (conn, nvertices, * nvertices_new);

    r_new = (double *) malloc ((* nvertices_new) * 8);
    z_new = (double *) malloc ((* nvertices_new) * 8);
    field_new = (double *) malloc ((*nvertices_new) * 8);
    assert (r_new && z_new && field_new);
   

    int * nodes_cut = build_nodes_cut_list (conn, nvertices, * nvertices_new);

    build_field (conn, nvertices, * nvertices_new, nodes_cut,
                 r, z, field, 
                 r_new, z_new, field_new);

    * nmesh_new = build_mesh (conn, nvertices, * nvertices_new, 
                              nmesh, nodes_cut, &mesh_new);

    free_nodes_cut_list (nodes_cut);

    * r_reduced = r_new;
    * z_reduced = z_new;
    * field_reduced = field_new;
    * mesh_reduced = mesh_new;

    free_conn (conn, nvertices);
    free_cost_matrix (cost_matrix, nvertices);
    pqueue_free (pq);

    for (int i = 0; i < nvertices; i++)
    {
        g_hash_table_destroy (nodes_ght[i]);
    }

    free (nodes_ght);

    free (r);
    free (z);
    free (field);
}

void thresholding (double * r, double * z, double * field,
                   int nelems, int * mesh, int nmesh,
                   double ** pnewr, double ** pnewz, 
                   double ** pnewfield, int ** pnewmesh,
                   int * pnewsize, int * pntaggedCells)
{
    int ntaggedCells = 0, cell_cnt = 0, newsize = 0;
    double surface_size = 0.0;

    for (int m = 0; m < nmesh; m++)
    {
        int n1 = * (mesh + m * 3);
        int n2 = * (mesh + m * 3 + 1);
        int n3 = * (mesh + m * 3 + 2);

        /* Gradient formular from Mark 
        grad u = u1 [y2-y3, x3-x2] + u2 [y3-y1, x1-x3] + u3 [y1-y2,x2-x1]
        */

        double avg_mag = (field[n1] + field[n2] + field[n3]) / 3;

        //grad[n1] = grad[n2] = grad[n3] = grad_mag;

        if (avg_mag > threshold)
        {
            ntaggedCells++;

            double a = sqrt (pow (r[n1] - r[n2], 2) + pow (z[n1] - z[n2], 2));
            double b = sqrt (pow (r[n1] - r[n3], 2) + pow (z[n1] - z[n3], 2));
            double c = sqrt (pow (r[n2] - r[n3], 2) + pow (z[n2] - z[n3], 2));
            double p = (a + b + c) / 2;

            surface_size += sqrt (p * (p - a) * (p - b) * (p - c));
        }
    }  // loop through the node connectivity array

    printf ("ntaggedCells = %d, orginal surface = %f\n", ntaggedCells, surface_size);

    double * newz = (double *) malloc (ntaggedCells * 3 * 8);
    double * newr = (double *) malloc (ntaggedCells * 3 * 8);
    double * newfield = (double *) malloc (ntaggedCells * 3 * 8);
    int * newmesh = (int *) malloc (ntaggedCells * 3 * 4);
    assert (newz && newr && newfield && newmesh);

    for (int m = 0; m < nmesh; m++)
    {
        int n1 = * (mesh + m * 3);
        int n2 = * (mesh + m * 3 + 1);
        int n3 = * (mesh + m * 3 + 2);

        double avg_mag = (field[n1] + field[n2] + field[n3]) / 3;

        if (avg_mag > threshold)
        {
            int tri1 = insert_node (newz, newr, newfield, &newsize,
                                    z[n1], r[n1], field[n1]);
            int tri2 = insert_node (newz, newr, newfield, &newsize,
                                    z[n2], r[n2], field[n2]);
            int tri3 = insert_node (newz, newr, newfield, &newsize,
                                    z[n3], r[n3], field[n3]);
            newmesh[cell_cnt++] = tri1;
            newmesh[cell_cnt++] = tri2;
            newmesh[cell_cnt++] = tri3;
        }
    }  // loop through the node connectivity

    * pnewz = newz;
    * pnewr = newr;
    * pnewfield = newfield;
    * pnewmesh = newmesh;
    * pntaggedCells = ntaggedCells;
    * pnewsize = newsize;
}

void extract_high_gradient (double * r, double * z, double * field,
                            int nelems, int * mesh, int nmesh,
                            double ** pnewr, double ** pnewz, 
                            double ** pnewfield, int ** pnewmesh,
                            int * pnewsize, int * pntaggedCells)
{
    int ntaggedCells = 0, cell_cnt = 0, newsize = 0;
    double * grad = (double *) malloc (nelems * 8);

    for (int m = 0; m < nmesh; m++)
    {
        int n1 = * (mesh + m * 3);
        int n2 = * (mesh + m * 3 + 1);
        int n3 = * (mesh + m * 3 + 2);

        /* Gradient formular from Mark 
        grad u = u1 [y2-y3, x3-x2] + u2 [y3-y1, x1-x3] + u3 [y1-y2,x2-x1]
        */

        double grad_z = field[n1] * (z[n2] - z[n3]) + field[n2] * (z[n3] - z[n1]) + field[n3]* (z[n1] - z[n2]);
        double grad_r = field[n1] * (r[n3] - r[n2]) + field[n2] * (r[n1] - r[n3]) + field[n3]* (r[n2] - r[n1]);
        double grad_mag = sqrt (pow (grad_z, 2) + pow (grad_r, 2));

        grad[n1] = grad[n2] = grad[n3] = grad_mag;

        if (grad_mag > threshold)
        {
            ntaggedCells++;
        }
    }  // loop through the node connectivity array

//  printf ("level = %d, ntaggedCells = %d\n", l, ntaggedCells);

    double * newz = (double *) malloc (ntaggedCells * 3 * 8);
    double * newr = (double *) malloc (ntaggedCells * 3 * 8);
    double * newfield = (double *) malloc (ntaggedCells * 3 * 8);
    int * newmesh = (int *) malloc (ntaggedCells * 3 * 4);
    assert (newz && newr && newfield && newmesh);

    for (int m = 0; m < nmesh; m++)
    {
        int n1 = * (mesh + m * 3);
        int n2 = * (mesh + m * 3 + 1);
        int n3 = * (mesh + m * 3 + 2);

        double grad_z = field[n1] * (z[n2] - z[n3]) + field[n2] * (z[n3] - z[n1]) + field[n3]* (z[n1] - z[n2]);
        double grad_r = field[n1] * (r[n3] - r[n2]) + field[n2] * (r[n1] - r[n3]) + field[n3]* (r[n2] - r[n1]);
        double grad_mag = sqrt (pow (grad_z, 2) + pow (grad_r, 2));

        grad[n1] = grad[n2] = grad[n3] = grad_mag;

        if (grad_mag > threshold)
        {
            int tri1 = insert_node (newz, newr, newfield, &newsize,
                                    z[n1], r[n1], field[n1]);
            int tri2 = insert_node (newz, newr, newfield, &newsize,
                                    z[n2], r[n2], field[n2]);
            int tri3 = insert_node (newz, newr, newfield, &newsize,
                                    z[n3], r[n3], field[n3]);
            newmesh[cell_cnt++] = tri1;
            newmesh[cell_cnt++] = tri2;
            newmesh[cell_cnt++] = tri3;
        }
    }  // loop through the node connectivity

    * pnewz = newz;
    * pnewr = newr;
    * pnewfield = newfield;
    * pnewmesh = newmesh;
    * pntaggedCells = ntaggedCells;
    * pnewsize = newsize;

    free (grad);
}

void extract_features (double * r, double * z, double * field,
                       int nelems, int * mesh, int nmesh,
                       double ** pnewr, double ** pnewz,
                       double ** pnewfield, int ** pnewmesh,
                       int * pnewsize, int * pntaggedCells)
{
    if (thresh_type == ABSOLUTE)
    {
        thresholding (r, z, field,
                      nelems, mesh, nmesh,
                      pnewr, pnewz,
                      pnewfield, pnewmesh,
                      pnewsize, pntaggedCells);
    }
    else
    {
        extract_high_gradient (r, z, field,
                      nelems, mesh, nmesh,
                      pnewr, pnewz,
                      pnewfield, pnewmesh,
                      pnewsize, pntaggedCells);
    }
}


#define DEFINE_VAR_LEVEL(varname, l, type)         \
    adios_common_define_var (md->level[l].grp      \
                            ,varname               \
                            ,var->path             \
                            ,type                  \
                            ,new_local_dimensions  \
                            ,new_global_dimensions \
                            ,new_local_offsets     \
                            )


void adios_sirius_adaptive_write (struct adios_file_struct * fd
                                 ,struct adios_var_struct * v
                                 ,const void * data
                                 ,struct adios_method_struct * method
                                 )

{
    struct adios_sa_data_struct * md = (struct adios_sa_data_struct *)
            method->method_data;
    struct var_struct * var;
    int i, l, ndims = count_dimensions (v->dimensions);
    int mesh_ndims;
    int type_size = adios_get_type_size (v->type,data);
    uint64_t varsize;
    uint64_t ldims[16], offsets[16], gdims[16];
    uint64_t mesh_ldims[16], mesh_offsets[16], mesh_gdims[16];
    uint64_t coord[16], lcoord[16], rcoord[16];
    uint64_t element, nelems, mesh_nelems;
    uint32_t ntaggedCells = 0;
    double * newz, * newr, * newfield;
    int * newmesh;
    int newsize = 0, cell_cnt = 0;
    double * r_reduced = 0, * z_reduced = 0;
    double * data_reduced = 0, * delta = 0;
    double * delta_compr = 0;
#ifdef TEST_REDUCTION
    double * test_field = 0;
#endif
    int nvertices_new;
    int * mesh_reduced = 0, nmesh_reduced;
    int compr_size = 0;

    for (l = 0; l < nlevels; l++)
    {
        if (alloc_var_struct (md, l) != err_no_error)
        {
            return;
        }

        md->level[l].varcnt++;
        var = md->level[l].vars;

        if (ndims == 0)
        {
            var->multidim = adios_flag_no;
            var->data = malloc(type_size);
            memcpy (var->data, data, type_size);

            // name the variable just like the original
            var->name = strdup (v->name);
            var->path = strdup (v->path);
            var->type = v->type;
            var->size = type_size;

            adios_common_define_var (md->level[l].grp
                                    ,var->name
                                    ,var->path
                                    ,var->type
                                    ,""
                                    ,""
                                    ,""
                                    );
        }
        else
        {
            //get the number of elements
            nelems = get_var_dimensions (v, ndims, gdims, ldims, offsets);
            varsize = nelems * type_size;

            if (v->type == adios_double)
            {
                // name the variable
                int len = 5 + strlen (v->name);

                var->name = (char *) malloc (len);
                sprintf (var->name, "%s/L%d", v->name, l);
                var->path = strdup (v->path);
                var->type = adios_double;

                if (l == 0)
                {
                    var->data = (void *) data;
                    var->global_dimensions = print_dimensions (1, gdims);
                    var->local_dimensions = print_dimensions (1, ldims);
                    var->local_offsets = print_dimensions (1, offsets);
                    var->size = varsize;
                }
                else if (l == 1)
                {
                    var->data = (void *) data;

                    if (!strcmp (v->name, "R")
                     || !strcmp (v->name, "Z"))
                    {
                        ldims[0] = gdims[0] = newsize;
                        offsets[0] = 0;

                        var->global_dimensions = print_dimensions (1, gdims);
                        var->local_dimensions = print_dimensions (1, ldims);
                        var->local_offsets = print_dimensions (1, offsets);

                    }
                    else
                    {
                        var->global_dimensions = print_dimensions (1, gdims);
                        var->local_dimensions = print_dimensions (1, ldims);
                        var->local_offsets = print_dimensions (1, offsets);
                    }
                } 
                else if (l == 2)
                {
                    var->data = (void *) data;
                     
                    if (!strcmp (v->name, "R")
                     || !strcmp (v->name, "Z"))
                    {
                        ldims[0] = gdims[0] = nvertices_new;
                        offsets[0] = 0;

                        var->global_dimensions = print_dimensions (1, gdims);
                        var->local_dimensions = print_dimensions (1, ldims);
                        var->local_offsets = print_dimensions (1, offsets);

                    }
                    else
                    {
                        var->global_dimensions = print_dimensions (1, gdims);
                        var->local_dimensions = print_dimensions (1, ldims);
                        var->local_offsets = print_dimensions (1, offsets);
                    }
                }

                if (!strcmp (v->name, "dpot"))
                {
                    if (l == 0)
                    {
                        struct adios_var_struct 
                            * mesh = adios_find_var_by_name (fd->group, "mesh");

                        assert (mesh);

                        mesh_ndims = count_dimensions (mesh->dimensions);
                        mesh_nelems = get_var_dimensions (mesh, 
                                                          mesh_ndims, 
                                                          mesh_gdims, 
                                                          mesh_ldims,  
                                                          mesh_offsets
                                                         );

                        assert (mesh_ldims[1] == 3);

                        struct adios_var_struct * R = adios_find_var_by_name (fd->group, "R");
                        assert (R);


                        struct adios_var_struct * Z = adios_find_var_by_name (fd->group, "Z");
                        assert (Z);
#if 0
                        // Decimation for level 0
                        decimate (R->data, Z->data, data, 
                                  nelems, mesh->data, mesh_ldims[0],
                                  &r_reduced, &z_reduced, &data_reduced, 
                                  &nvertices_new, &mesh_reduced, 
                                  &nmesh_reduced
                                 );
#endif

#ifdef DUMP_FEATURE 
                        extract_features ((double *) R->data, 
                                          (double *) Z->data, 
                                          (double *) data, 
                                          nelems, (int *) mesh->data,
                                          mesh_ldims[0],
                                          &newr, &newz, &newfield, &newmesh,
                                          &newsize, &ntaggedCells);
#endif
#if 0
                        for (int m = 0; m < mesh_ldims[0]; m++)
                        {
                            int n1 = * ((int *) mesh->data + m * 3);
                            int n2 = * ((int *) mesh->data + m * 3 + 1);
                            int n3 = * ((int *) mesh->data + m * 3 + 2);

                            double * field = (double *) data;
                            double * r = (double *) R->data;
                            double * z = (double *) Z->data;
                            /* Gradient formular from Mark 
                               grad u = u1 [y2-y3, x3-x2] + u2 [y3-y1, x1-x3] + u3 [y1-y2,x2-x1]
                             */

                            double grad_z = field[n1] * (z[n2] - z[n3]) + field[n2] * (z[n3] - z[n1]) + field[n3]* (z[n1] - z[n2]);
                            double grad_r = field[n1] * (r[n3] - r[n2]) + field[n2] * (r[n1] - r[n3]) + field[n3]* (r[n2] - r[n1]);
                            double grad_mag = sqrt (pow (grad_z, 2) + pow (grad_r, 2));

                            grad[n1] = grad[n2] = grad[n3] = grad_mag;

                            //TODO: To add threshold stuff
                            if (grad_mag > threshold)
                            {
                                ntaggedCells++;
                            }
                        }  // loop through the node connectivity array

//                        printf ("level = %d, ntaggedCells = %d\n", l, ntaggedCells);

                        newz = (double *) malloc (ntaggedCells * 3 * 8);
                        newr = (double *) malloc (ntaggedCells * 3 * 8);
                        newfield = (double *) malloc (ntaggedCells * 3 * 8);
                        newmesh = (int *) malloc (ntaggedCells * 3 * 4);
                        assert (newz && newr && newfield && newmesh);

                        for (int m = 0; m < mesh_ldims[0]; m++)
                        {
                            int n1 = * ((int *) mesh->data + m * 3);
                            int n2 = * ((int *) mesh->data + m * 3 + 1);
                            int n3 = * ((int *) mesh->data + m * 3 + 2);

                            double * field = (double *) data;
                            double * r = (double *) R->data;
                            double * z = (double *) Z->data;

                            double grad_z = field[n1] * (z[n2] - z[n3]) + field[n2] * (z[n3] - z[n1]) + field[n3]* (z[n1] - z[n2]);
                            double grad_r = field[n1] * (r[n3] - r[n2]) + field[n2] * (r[n1] - r[n3]) + field[n3]* (r[n2] - r[n1]);
                            double grad_mag = sqrt (pow (grad_z, 2) + pow (grad_r, 2));

                            grad[n1] = grad[n2] = grad[n3] = grad_mag;

                            //TODO: To add threshold stuff
                            if (grad_mag > threshold)
                            {
                                int tri1 = insert_node (newz, newr, newfield, &newsize,
                                     z[n1], r[n1], field[n1]);
                                int tri2 = insert_node (newz, newr, newfield, &newsize,
                                     z[n2], r[n2], field[n2]);
                                int tri3 = insert_node (newz, newr, newfield, &newsize,
                                     z[n3], r[n3], field[n3]);
                                newmesh[cell_cnt++] = tri1;
                                newmesh[cell_cnt++] = tri2;
                                newmesh[cell_cnt++] = tri3;
                            }
                        }  // loop through the node connectivity
#endif

                        // Decimation for level 0
                        decimate ((double *) R->data, (double *) Z->data, (double *) data, 
                                  nelems, (int *) mesh->data, mesh_ldims[0],
                                  &r_reduced, &z_reduced, &data_reduced, 
                                  &nvertices_new, &mesh_reduced, 
                                  &nmesh_reduced
                                 );

                        calc_area (r_reduced, z_reduced, data_reduced,
                                   nvertices_new, mesh_reduced,
                                   nmesh_reduced
                                  );

                        if (save_delta)
                        {
                            get_delta ((double *) R->data, 
                                       (double *) Z->data, 
                                       (double *) data,
                                       nelems, (int *) mesh->data, 
                                       mesh_ldims[0],
                                       r_reduced, z_reduced, data_reduced,
                                       nvertices_new, mesh_reduced, 
                                       nmesh_reduced, &delta
                                      );

                            if (compress_delta)
                            {
                                compr_size = compress (delta, 
                                                       nelems, 
                                                       compr_tolerance, 
                                                       &delta_compr
                                                      );
                            }
                        }
#ifdef TEST_REDUCTION
                        if (save_delta)
                        {
                            if (compress_delta)
                            {
                                decompress (delta, nelems, 1e-3,
                                            delta_compr, compr_size);
                            }

                            test_delta ((double *) R->data,
                                       (double *) Z->data,
                                       (double *) data,
                                       nelems, (int *) mesh->data, mesh_ldims[0],
                                       r_reduced, z_reduced, data_reduced,
                                       nvertices_new, mesh_reduced, nmesh_reduced,
                                       delta, &test_field
                                      );
                        }
#endif
                    }
                    else if (l == 1)
                    {
                    }
                    else if (l == 2)
                    {
                    }
                }  // if dpot
            } // if double
            else
            {
                /* Not a double */
                if (l == 0)
                {
                    //only in level 0 do we need to store this variable
                    int len = 5 + strlen (v->name);

                    var->name = (char *) malloc (len);
                    sprintf (var->name, "%s/L%d", v->name, l);

                    var->path = strdup (v->path);
                    var->type = v->type;
                    var->size = varsize;
                    //FIXME
                    var->data = (void*)v->data;
                    var->global_dimensions = print_dimensions (ndims, gdims);
                    var->local_dimensions = print_dimensions (ndims, ldims);
                    var->local_offsets = print_dimensions (ndims, offsets);
                }
                else
                {
                    var->size = 0;
                }
            }

            if (var->size > 0)
            {
                adios_common_define_var (md->level[l].grp
                                        ,var->name
                                        ,var->path
                                        ,var->type
                                        ,var->local_dimensions
                                        ,var->global_dimensions
                                        ,var->local_offsets
                                        );

                if (!strcmp (v->name, "dpot") && l == 0)
                {
                    char * new_global_dimensions;
                    char * new_local_dimensions;
                    char * new_local_offsets;
                    uint64_t new_gdims[16];
                    uint64_t new_ldims[16];
                    uint64_t new_offsets[16];
#ifdef DUMP_FEATURE
                    new_gdims[0] = newsize;
                    new_ldims[0] = newsize;
                    new_offsets[0] = 0;

                    new_global_dimensions = print_dimensions (1, new_gdims);
                    new_local_dimensions = print_dimensions (1, new_ldims);
                    new_local_offsets = print_dimensions (1, new_offsets);
                    
                    DEFINE_VAR_LEVEL("R/L2",2,adios_double);
                    DEFINE_VAR_LEVEL("Z/L2",2,adios_double);
                    DEFINE_VAR_LEVEL("dpot/L2",2,adios_double);

                    new_ldims[0] = ntaggedCells;
                    new_ldims[1] = 3;
                    new_local_dimensions = print_dimensions (2, new_ldims);
                    new_offsets[0] = 0;
                    new_offsets[1] = 0;
                    new_local_offsets = print_dimensions (2, new_offsets);
                    new_global_dimensions = "";

                    DEFINE_VAR_LEVEL("mesh/L2",2,adios_integer);
#endif
///////////
                    new_gdims[0] = nvertices_new;
                    new_ldims[0] = nvertices_new;
                    new_offsets[0] = 0;

                    new_global_dimensions = print_dimensions (1, new_gdims);
                    new_local_dimensions = print_dimensions (1, new_ldims);
                    new_local_offsets = print_dimensions (1, new_offsets);

                    DEFINE_VAR_LEVEL("R/L1",1,adios_double);
                    DEFINE_VAR_LEVEL("Z/L1",1,adios_double);
                    DEFINE_VAR_LEVEL("dpot/L1",1,adios_double);

                    new_ldims[0] = nmesh_reduced;
                    new_ldims[1] = 3;
                    new_local_dimensions = print_dimensions (2, new_ldims);
                    new_offsets[0] = 0;
                    new_offsets[1] = 0;
                    new_local_offsets = print_dimensions (2, new_offsets);
                    new_global_dimensions = "";

                    DEFINE_VAR_LEVEL("mesh/L1",1,adios_integer);

                    if (compress_delta)
                    {
                        new_gdims[0] = compr_size;
                        new_ldims[0] = compr_size;
                        new_offsets[0] = 0;
                    }
                    else
                    {
                        new_gdims[0] = nelems;
                        new_ldims[0] = nelems;
                        new_offsets[0] = 0;
                    }

                    new_global_dimensions = print_dimensions (1, new_gdims);
                    new_local_dimensions = print_dimensions (1, new_ldims);
                    new_local_offsets = print_dimensions (1, new_offsets);

                    if (compress_delta)
                    {
                        DEFINE_VAR_LEVEL("delta/L0",0,adios_byte);
                    }
                    else
                    {
                        DEFINE_VAR_LEVEL("delta/L0",0,adios_double);
                    }
#ifdef TEST_REDUCTION
                    DEFINE_VAR_LEVEL("full_field/L1",1,adios_double);
#endif
                }
            }

        }

        md->level[l].totalsize += var->size;

        if ( (!strcmp (v->name, "R") || !strcmp (v->name, "Z") 
             || !strcmp (v->name, "mesh") || !strcmp (v->name, "dpot")) 
           && (l == 1 || l == 2))
        {
            // do not write R and Z is level 1
        }
        else
        {
            // write it out
            if (md->level[l].vars->size > 0)
            {
                if (strcmp (v->name, "dpot"))
                {
                    do_write (md->level[l].fd, var->name, var->data);
                }

                if (!strcmp (v->name, "dpot") && l == 0)
                {
#ifdef DUMP_FEATURE
                    do_write (md->level[2].fd, "R/L2", newr);
                    free (newr);

                    do_write (md->level[2].fd, "Z/L2", newz);
                    free (newz);

                    do_write (md->level[2].fd, "mesh/L2", newmesh);
                    free (newmesh);

                    do_write (md->level[2].fd, "dpot/L2", newfield);
                    free (newfield);
#endif
                    do_write (md->level[1].fd, "R/L1", r_reduced);
                    free (r_reduced);

                    do_write (md->level[1].fd, "Z/L1", z_reduced);
                    free (z_reduced);

                    do_write (md->level[1].fd, "mesh/L1", mesh_reduced);
                    free (mesh_reduced);

                    do_write (md->level[1].fd, "dpot/L1", data_reduced);
                    free (data_reduced);

                    if (save_delta)
                    {
                        if (!compress_delta)
                        {
                            do_write (md->level[0].fd, "delta/L0", delta);
                        }
                        else
                        {
                            do_write (md->level[0].fd, "delta/L0", 
                                      delta_compr);
                            free (delta_compr);
                        }

                        free (delta);
                    }
                    else
                    {
                        // write the full dpot data.
                        do_write (md->level[l].fd, var->name, var->data);
                    }
#ifdef TEST_REDUCTION
                    do_write (md->level[1].fd, "full_field/L1", test_field);
                    free (test_field);
#endif
                }
            }
        } // if

    } // for levels

}

void adios_sirius_adaptive_read (struct adios_file_struct * fd
                        ,struct adios_var_struct * v, void * buffer
                        ,uint64_t buffer_size
                        ,struct adios_method_struct * method
    )

{
}

void adios_sirius_adaptive_buffer_overflow (struct adios_file_struct * fd,
                                   struct adios_method_struct * method)
{
    struct adios_sa_data_struct * md = (struct adios_sa_data_struct *)
        method->method_data;
    log_error ("rank %d: SIRIUS method only works with complete buffering of data between adios_open() "
               "and adios_close(). Variables that do not fit into the buffer will not be "
               "written by this method to file %s\n", md->rank, fd->name);
}

#define FREE(v) if (v!=NULL) {free(v); v=NULL;}

void release_resource_at_close (struct adios_sa_data_struct * md)
{
    int l;
    for (l=0; l < nlevels; l++)
    {
        FREE (md->level[l].filename);
        FREE (md->level[l].grp_name);

        struct var_struct *next;
        struct var_struct *vars = md->level[l].vars_head;
        while (vars)
        {
            next=vars->next;
            FREE(vars->data);
            FREE(vars->local_dimensions);
            FREE(vars->global_dimensions);
            FREE(vars->local_offsets);
            FREE(vars);
            vars=next;
        }
        md->level[l].varcnt = 0;
    }
}

void * threaded_call_common_close(void *lp)
{
    struct level_struct *level = (struct level_struct *)lp;
    if(level == NULL)
        return NULL;
    common_adios_close (level->fd);
    return NULL;
}

void adios_sirius_adaptive_close (struct adios_file_struct * fd
                         ,struct adios_method_struct * method
    )
{
    struct adios_sa_data_struct * md = (struct adios_sa_data_struct *)
        method->method_data;

    switch (fd->mode)
    {
        case adios_mode_read:
        {
            adios_error (err_invalid_file_mode, "SIRIUS method: Read mode is not supported.\n");
            break;
        }
        case adios_mode_append:
        case adios_mode_update:
        case adios_mode_write:
        {
            int l;
            for (l=0; l < nlevels; l++)
            {
                common_adios_close (md->level[l].fd);
            }

            release_resource_at_close (md);
            break;
        }
        default:
        {
            adios_error (err_invalid_file_mode, "SIRIUS method: Unknown file mode requested: %d\n", fd->mode);
            break;
        }
    }

    return;
}

void adios_sirius_adaptive_get_write_buffer (struct adios_file_struct * fd
                                    ,struct adios_var_struct * v
                                    ,uint64_t * size
                                    ,void ** buffer
                                    ,struct adios_method_struct * method
    )
{
}

void adios_sirius_adaptive_finalize (int mype, struct adios_method_struct * method)
{
    int l;
    for (l=0; l < nlevels; l++)
    {
        if (io_method[l])
            FREE (io_method[l]);
        if (io_parameters[l])
            FREE (io_parameters[l]);
        if (io_paths[l])
            FREE (io_paths[l]);
    }
}

void adios_sirius_adaptive_end_iteration (struct adios_method_struct * method)
{
}

void adios_sirius_adaptive_start_calculation (struct adios_method_struct * method)
{
}

void adios_sirius_adaptive_stop_calculation (struct adios_method_struct * method)
{
}
