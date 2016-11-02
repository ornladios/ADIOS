//passed the test for 1024 cores and level-3 spatial aggregation

#include <unistd.h>
#include <fcntl.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <errno.h>
#include <assert.h>
#include <sys/stat.h>
// xml parser
#include <mxml.h>

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

#include "mpi.h"

#define MAXLEVEL 10
#define Y_POS 0
#define X_POS 1

static char *io_method[MAXLEVEL]; //the IO methods for data output for each level
static char *io_parameters[MAXLEVEL]; //the IO method parameters
static char *io_paths[MAXLEVEL]; //the IO method output paths (prefix to filename)
static int nlevels=2; // Number of levels
static char *splitter_type;

double threshold = 0.0001;
double fillratio = 0.28;

typedef struct Box
{
    uint64_t lower[16];
    uint64_t upper[16];
} Box;

void print_2d_box (Box * b)
{
    printf ("Box: (%ld, %ld), (%ld, %ld)\n", b->lower[0], b->lower[1], b->upper[0], b->upper[1]);
}

typedef struct BoxList
{
    Box * box;
    void * data;
    struct BoxList * next;
} BoxList;

BoxList * bl_head = 0;

void print_list (BoxList * h)
{
    printf ("Box List: \n");
    while (h)
    {
        printf ("(%lu, %lu):(%lu, %lu)   ", h->box->lower[0], h->box->lower[1],
                h->box->upper[0], h->box->upper[1]);
        h = h->next;
    }
    printf ("\n");
}

void box_append (BoxList ** head, BoxList * bl)
{
    assert (head);
    assert (bl);

    if (* head == 0)
    {
        * head = bl;
    }
    else
    {
        BoxList * h = * head;
        while (h->next)
        {
            h = h->next;
        }
        
        h->next = bl;
    }
}

// merge boxes if two are adjacent to each other
void box_merge (BoxList ** head)
{
    BoxList * h = * head, * n, * t;
    int should_merge;

    while (h)
    {
        n = h->next;
        while (n)
        {
            should_merge = 0;

            // check left/right
            if (n->box->lower[Y_POS] == h->box->lower[Y_POS])
            {
                if ((n->box->upper[X_POS] + 1) == h->box->lower[X_POS]
                        && n->box->upper[Y_POS] == h->box->upper[Y_POS])
                {
                    // yes
                    should_merge = 1;
                    h->box->lower[Y_POS] = n->box->lower[Y_POS];
                    h->box->lower[X_POS] = n->box->lower[X_POS];
                }

                if (n->box->lower[X_POS] == (h->box->upper[X_POS] + 1)
                        && n->box->upper[Y_POS] == h->box->upper[Y_POS])
                {
                    // yes
                    should_merge = 1;
                    h->box->upper[Y_POS] = n->box->upper[Y_POS];
                    h->box->upper[X_POS] = n->box->upper[X_POS];
                }
            } 
            // check upper/lower
            else if (n->box->lower[X_POS] == h->box->lower[X_POS])
            {
                if ((n->box->upper[Y_POS] + 1) == h->box->lower[Y_POS]
                        && n->box->upper[X_POS] == h->box->upper[X_POS])
                {
                    // yes
                    should_merge = 1;
                    h->box->lower[Y_POS] = n->box->lower[Y_POS];
                    h->box->lower[X_POS] = n->box->lower[X_POS];
                }
                if (n->box->lower[Y_POS] == (h->box->upper[Y_POS] + 1)
                         && n->box->upper[X_POS] == h->box->upper[X_POS])
                {
                    // yes
                    should_merge = 1;
                    h->box->upper[Y_POS] = n->box->upper[Y_POS];
                    h->box->upper[X_POS] = n->box->upper[X_POS];
                }
            }
       
            if (should_merge)
            {
                h->next = n->next;
                free (n->box);
                t = n->next;
                free (n);
                n = t;
            }
            else
            {
                n = n->next;
            }
        }

        h = h->next;
    }
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

    splitter_type = strdup("float");
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
        }
        else if (!strcasecmp(p->name, "type"))
        {
            errno = 0;
            free(splitter_type);
            splitter_type = strdup(p->value);            
            fprintf(stderr, "set param type = %s\n", splitter_type);
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

#if 1
BoxList * create_partitions (uint8_t ** tags, uint64_t ndims, uint64_t * ldims, Box * domain)
{
    Box * minbox = (Box *) malloc (sizeof(Box));
    assert (minbox);

    if (ndims == 2)
    {
        int i, j;
        int64_t minx = domain->upper[1] + 1, miny = domain->upper[0] + 1; 
        int64_t maxx = domain->lower[1] - 1, maxy = domain->lower[0] - 1;

        //printf ("domain = (%lu, %lu), (%lu, %lu)\n", domain->lower[0], domain->lower[1], domain->upper[0], domain->upper[1]);
        // first, create a min box
        uint32_t n_tagged_cells = 0, n_cells = 0;

        // check fill ratio
        for (j = domain->lower[0]; j <= domain->upper[0]; j++)
        {
            for (i = domain->lower[1]; i <= domain->upper[1]; i++)
            {
                if (tags[j][i])
                {
                    //printf ("(%d,%d) is tagged.\n", j, i);
                    n_tagged_cells++;
                    if (i < minx)
                        minx = i;
                    if (j < miny)
                        miny = j;
                    if (i > maxx)
                        maxx = i;
                    if (j > maxy)
                        maxy = j;
                }
            }
        }
        if (n_tagged_cells == 0)
        {
            //printf ("no tagged cells.\n");
            return 0;
        }

        minbox->lower[0] = miny;
        minbox->lower[1] = minx;
        minbox->upper[0] = maxy;
        minbox->upper[1] = maxx;

        n_cells = (minbox->upper[0] - minbox->lower[0] + 1) * (minbox->upper[1] - minbox->lower[1] + 1);
        //printf ("minbox = (%d, %d), (%d, %d), fillratio = %f\n", 
        //        minbox->lower[0], minbox->lower[1], minbox->upper[0], minbox->upper[1], (double)n_tagged_cells / (double)n_cells);

        if (((double)n_tagged_cells / (double)n_cells >= fillratio) 
                || (minbox->upper[0] - minbox->lower[0] + 1) <= 2
                || (minbox->upper[1] - minbox->lower[1] + 1) <= 2) 
        {
            // good
            BoxList * bl = (BoxList *) malloc (sizeof (BoxList));
            assert (bl);
            bl->box = minbox;
            bl->data = 0;
            bl->next = 0;

            return bl;
        }
        else
        {
            printf ("This box is under filled.\n");
            uint64_t hy = (minbox->upper[Y_POS] - minbox->lower[Y_POS] + 1) / 2;
            uint64_t hx = (minbox->upper[X_POS] - minbox->lower[X_POS] + 1) / 2;
            //printf ("hy = %lu, hx = %lu\n", hy, hx);

            Box ul;
            ul.lower[Y_POS] = minbox->lower[Y_POS] + hy;
            ul.lower[X_POS] = minbox->lower[X_POS];
            ul.upper[Y_POS] = minbox->upper[Y_POS];
            ul.upper[X_POS] = minbox->lower[X_POS] + hx - 1;

            Box ll;
            ll.lower[Y_POS] = minbox->lower[Y_POS];
            ll.lower[X_POS] = minbox->lower[X_POS];
            ll.upper[Y_POS] = minbox->lower[Y_POS] + hy - 1;
            ll.upper[X_POS] = minbox->lower[X_POS] + hx - 1;
            
            Box ur;
            ur.lower[Y_POS] = minbox->lower[Y_POS] + hy;
            ur.lower[X_POS] = minbox->lower[X_POS] + hx;
            ur.upper[Y_POS] = minbox->upper[Y_POS];
            ur.upper[X_POS] = minbox->upper[X_POS];

            Box lr;
            lr.lower[Y_POS] = minbox->lower[Y_POS];
            lr.lower[X_POS] = minbox->lower[X_POS] + hx;
            lr.upper[Y_POS] = minbox->lower[Y_POS] + hy - 1;
            lr.upper[X_POS] = minbox->upper[X_POS];

            // need to further split the domain
            BoxList * new_bl;
            if (new_bl = create_partitions (tags, ndims, ldims, &ul))
                box_append (&bl_head, new_bl);

            if (new_bl = create_partitions (tags, ndims, ldims, &ur))
                box_append (&bl_head, new_bl);

            if (new_bl = create_partitions (tags, ndims, ldims, &ll))
                box_append (&bl_head, new_bl);

            if (new_bl = create_partitions (tags, ndims, ldims, &lr))
                box_append (&bl_head, new_bl);

            return 0;
        }
    }
    else if (ndims == 3)
    {
        printf ("3D is NOT supported yet.\n");
    }
}

#endif

uint8_t ** alloc_2d_uint8_array (uint64_t * ldims)
{
    int i, j;
    uint8_t ** p = (uint8_t *) malloc (ldims[0] * 8);
    assert (p);

    for (i = 0; i < ldims[0]; i++)
    {
        p[i] = (uint8_t *) malloc (ldims[1] * 1);
        assert (p[i]);
    }

    return p;
}


void free_2d_uint8_array (uint8_t ** p, uint64_t * ldims)
{
    int i, j;

    for (i = 0; i < ldims[0]; i++)
    {
        free (p[i]);
    }

    free (p);
}

void copy_box_data (void * data, uint64_t * ldims, BoxList * bl_head)
{
    BoxList * h = bl_head;
    uint64_t i;

    while (h)
    {
        uint64_t y = h->box->upper[Y_POS] - h->box->lower[Y_POS] + 1;
        uint64_t x = h->box->upper[X_POS] - h->box->lower[X_POS] + 1;

        h->data = malloc (x * y * 8);
        assert (h->data);

        for (i = h->box->lower[Y_POS]; i <= h->box->upper[Y_POS]; i++)
        {
            memcpy ((double *) h->data + (i - h->box->lower[Y_POS]) * x, 
                    (double *) data 
                        + (h->box->lower[Y_POS] + i - h->box->lower[Y_POS]) 
                          * ldims[X_POS] + h->box->lower[X_POS], 
                    x * 8); 
        }

        h = h->next;
    }
}

int insert_node (double * newz, double * newr, double * newfield, int * size,
                  double z, double r, double field)
{
    int found;
printf ("(%e, %e)\n", z, r);
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
                    var->data = data;
                    var->global_dimensions = print_dimensions (1, gdims);
                    var->local_dimensions = print_dimensions (1, ldims);
                    var->local_offsets = print_dimensions (1, offsets);
                    var->size = varsize;
                }
                else if (l == 1)
                {
                    var->data = data;

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

                if (!strcmp (v->name, "dpot"))
                {
                    if (l == 0)
                    {
                        double * grad = malloc (nelems * 8);
                        assert (grad);

                        struct adios_var_struct 
                            * mesh = adios_find_var_by_name (fd->group, "mesh");

                        if (!mesh)
                        {
                            adios_error (err_invalid_varname, 
                                 "Bad var name (ignored) in SIRIUS_ADAPTIVE"
                                 " adios_write(): %s\n", mesh->name);
                            return 1;
                        }

                        mesh_ndims = count_dimensions (mesh->dimensions);
                        mesh_nelems = get_var_dimensions (mesh, 
                                                          mesh_ndims, 
                                                          mesh_gdims, 
                                                          mesh_ldims,  
                                                          mesh_offsets
                                                         );

                        if (mesh_ldims[1] != 3)
                        {
                            printf ("The mesh is incorrect!\n");
                            return 1;
                        }

                        struct adios_var_struct * R = adios_find_var_by_name (fd->group, "R");

                        if (!R)
                        {
                            adios_error (err_invalid_varname,
                                 "Bad var name (ignored) in SIRIUS_ADAPTIVE"
                                 " adios_write(): %s\n", R->name);
                            return 1;
                        }

                        struct adios_var_struct * Z = adios_find_var_by_name (fd->group, "Z");

                        if (!Z)
                        {
                            adios_error (err_invalid_varname,
                                 "Bad var name (ignored) in SIRIUS_ADAPTIVE"
                                 " adios_write(): %s\n", Z->name);
                            return 1;
                        }

                        for (int m = 0; m < mesh_ldims[0]; m++)
                        {
                            int n1 = * ((int *) mesh->data + m * 3);
                            int n2 = * ((int *) mesh->data + m * 3 + 1);
                            int n3 = * ((int *) mesh->data + m * 3 + 2);

                            double * field = data;
                            double * r = R->data;
                            double * z = Z->data;
                            /* Gradient formular from Mark 
                               grad u = u1 [y2-y3, x3-x2] + u2 [y3-y1, x1-x3] + u3 [y1-y2,x2-x1]
                             */

                            double grad_z = field[n1] * (z[n2] - z[n3]) + field[n2] * (z[n3] - z[n1]) + field[n3]* (z[n1] - z[n2]);
                            double grad_r = field[n1] * (r[n3] - r[n2]) + field[n2] * (r[n1] - r[n3]) + field[n3]* (r[n2] - r[n1]);
                            double grad_mag = sqrt (grad_z * grad_z + grad_r * grad_r);

                            grad[n1] = grad[n2] = grad[n3] = grad_mag;

                            //TODO: To add threshold stuff
                            if (grad_mag > 0.2)
                            {
                                ntaggedCells++;
                            }
                        }  // loop through the node connectivity array

                        printf ("level = %d, ntaggedCells = %d\n", l, ntaggedCells);
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

                            double * field = data;
                            double * r = R->data;
                            double * z = Z->data;

                            double grad_z = field[n1] * (z[n2] - z[n3]) + field[n2] * (z[n3] - z[n1]) + field[n3]* (z[n1] - z[n2]);
                            double grad_r = field[n1] * (r[n3] - r[n2]) + field[n2] * (r[n1] - r[n3]) + field[n3]* (r[n2] - r[n1]);
                            double grad_mag = sqrt (grad_z * grad_z + grad_r * grad_r);

                            grad[n1] = grad[n2] = grad[n3] = grad_mag;

                            //TODO: To add threshold stuff
                            if (grad_mag > 0.2)
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
                    }
                    else if (l == 1)
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
                    var->name = strdup (v->name);
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

                    new_gdims[0] = newsize;
                    new_ldims[0] = newsize;
                    new_offsets[0] = 0;

                    new_global_dimensions = print_dimensions (1, &new_gdims);
                    new_local_dimensions = print_dimensions (1, &new_ldims);
                    new_local_offsets = print_dimensions (1, &new_offsets);

                    adios_common_define_var (md->level[1].grp
                                        ,"R/L1"
                                        ,var->path
                                        ,var->type
                                        ,new_local_dimensions
                                        ,new_global_dimensions
                                        ,new_local_offsets
                                        );

                    adios_common_define_var (md->level[1].grp
                                        , "Z/L1"
                                        ,var->path
                                        ,var->type
                                        ,new_local_dimensions
                                        ,new_global_dimensions
                                        ,new_local_offsets
                                        );

                    adios_common_define_var (md->level[1].grp
                                        , "dpot/L1"
                                        ,var->path
                                        ,var->type
                                        ,new_local_dimensions
                                        ,new_global_dimensions
                                        ,new_local_offsets
                                        );

                    new_ldims[0] = ntaggedCells;
                    new_ldims[1] = 3;
                    new_local_dimensions = print_dimensions (2, &new_ldims);
                    new_offsets[0] = 0;
                    new_offsets[1] = 0;
                    new_local_offsets = print_dimensions (2, &new_offsets);

                    adios_common_define_var (md->level[1].grp
                                        , "mesh/L1"
                                        ,var->path
                                        ,adios_integer
                                        ,new_local_dimensions
                                        ,""
                                        ,new_local_offsets
                                        );

                }
            }

        }

        md->level[l].totalsize += var->size;

        if ( (!strcmp (v->name, "R") || !strcmp (v->name, "Z") 
             || !strcmp (v->name, "mesh") || !strcmp (v->name, "dpot")) 
           && l == 1)
        {
            // do not write R and Z is level 1
        }
        else
        {
            // write it out
            if (md->level[l].vars->size > 0)
            {
                do_write (md->level[l].fd, var->name, var->data);

                if (!strcmp (v->name, "dpot") && l == 0)
                {
                    do_write (md->level[1].fd, "R/L1", newr);
                    do_write (md->level[1].fd, "Z/L1", newz);
                    do_write (md->level[1].fd, "mesh/L1", newmesh);
                    do_write (md->level[1].fd, "dpot/L1", newfield);
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
