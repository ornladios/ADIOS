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
#include "split-interface.h"
//threaded support
#include <pthread.h>

#define MAXLEVEL 10

static char *io_method[MAXLEVEL]; //the IO methods for data output for each level
static char *io_parameters[MAXLEVEL]; //the IO method parameters
static char *io_paths[MAXLEVEL]; //the IO method output paths (prefix to filename)
static int nlevels=1; // Number of storage levels
static int no_thread = 0;
static char *splitter_type;

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

struct adios_sirius_data_struct
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
    int ret;
    ret = adios_common_declare_group (id, name, adios_flag_no
                                      ,""
                                      ,""
                                      ,time_index
                                      ,adios_flag_no
        );
    if (ret == 1) {
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
    return adios_common_select_method_by_group_id (0, method, parameters, group ,base_path, 0);
}

static void define_iogroups(struct adios_method_struct * method)
{
    int len;
    int l;
    struct adios_sirius_data_struct * md = (struct adios_sirius_data_struct *) method->method_data;
    
    for (l=0; l < nlevels; l++)
    {
        len=5+strlen(method->group->name); //new groupname= tg_groupname
        md->level[l].grp_name=(char *)malloc(len);
        memset(md->level[l].grp_name, 0x00, len);
        sprintf(md->level[l].grp_name, "%s_L%d",method->group->name, l);
        declare_group (&(md->level[l].grp), md->level[l].grp_name, "", adios_flag_yes);
        select_method (md->level[l].grp, io_method[l], io_parameters[l],"");
    }
}

static int convert_file_mode(enum ADIOS_METHOD_MODE mode, char * file_mode)
{
    switch (mode)
    {
        case adios_mode_read:
            strcpy(file_mode,"r");
            break;

        case adios_mode_write:
            strcpy(file_mode,"w");
            break;

        case adios_mode_append:
            strcpy(file_mode,"a");
            break;

        case adios_mode_update:
            strcpy(file_mode,"u");
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


void adios_sirius_init(const PairStruct * parameters,
                       struct adios_method_struct * method)
{
    struct adios_sirius_data_struct * md = (struct adios_sirius_data_struct *)
        method->method_data;

    method->method_data = malloc (sizeof (struct adios_sirius_data_struct));
    md = (struct adios_sirius_data_struct *) method->method_data;

    char *threaded_writes = getenv("DOTHREAD");
    if(threaded_writes == NULL)
        no_thread = 1;
    else
        no_thread = 0;

    initialize_splitter();
    init_output_parameters(parameters);
}


static void init_method_parameters(struct adios_sirius_data_struct * md)
{
    int l;
    for(l=0; l < nlevels; l++) {
        md->level[l].varcnt=0;
        md->level[l].vars=NULL;
        md->level[l].vars_head=NULL;
        md->level[l].fd = 0;
        md->level[l].filename = NULL;
        md->level[l].grp_name = NULL;
        md->level[l].totalsize = 0;
    }
}


int adios_sirius_open (struct adios_file_struct * fd
                       ,struct adios_method_struct * method, MPI_Comm comm)
{

    struct adios_sirius_data_struct * md = (struct adios_sirius_data_struct *)
        method->method_data;
    char mode[2];
    int l;

    switch (fd->mode)
    {
        case adios_mode_read:
        {
            adios_error (err_invalid_file_mode, "SIRIUS method: Read mode is not supported.\n");
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

            //need to get the parameters form XML
            //init_output_parameters(method->parameters);
            init_method_parameters(md);
            define_iogroups(method);

            for (l=0; l < nlevels; l++)
            {
                //check if the directory exists and create it if it doesn't
                struct stat sb;
                if((stat(io_paths[l], &sb) != 0) || !S_ISDIR(sb.st_mode))
                {
                    //directory doesn't exist
                    //FIXME: there is a case where something already exists but
                    //isn't a directory. Hard to imagine though so I am ignoring
                    //it for the time beingmdki
                    mkdir(io_paths[l], 0700);
                }
                md->level[l].filename = malloc (strlen(io_paths[l]) + strlen(fd->name) + 1);
                sprintf (md->level[l].filename, "%s/%s", io_paths[l], fd->name);
                convert_file_mode(fd->mode, mode);
                common_adios_open( &(md->level[l].fd), md->level[l].grp_name, md->level[l].filename, mode, comm);
            }
            break;
        }
        default:
        {
            adios_error (err_invalid_file_mode, "SIRIUS method: Unknown file mode requested: %d\n", fd->mode);
            return adios_flag_no;
        }
    }

    return 1;
}

enum BUFFERING_STRATEGY adios_sirius_should_buffer (struct adios_file_struct * fd
                                                    ,struct adios_method_struct * method)
{
    //this method handles its own buffering
    return no_buffering;
}


//initial variable structure
static void init_var (struct var_struct *var)
{
    var->name = NULL;
    var->path = NULL;
    var->type = adios_unknown;
    var->next=NULL;
    var->global_dimensions= (char *) calloc (128, sizeof(char));
    var->local_dimensions = (char *) calloc (128, sizeof(char));
    var->local_offsets= (char *) calloc (128, sizeof(char));
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


static enum ADIOS_ERRCODES alloc_var (struct adios_sirius_data_struct * md, int level)
{
    struct var_struct *var = (struct var_struct *) malloc (sizeof(struct var_struct));
    if (!var) {
        adios_error (err_no_memory, "No memory to allocate yet another var in SIRIUS method\n");
        return err_no_memory;
    }

    var->prev=md->level[level].vars;
    var->next=NULL;
    if (md->level[level].varcnt == 0)
    {
        //assign the header of the variable list
        md->level[level].vars_head = var;
    }
    md->level[level].vars = var;

    // initialize the variable structure
    init_var (md->level[level].vars);
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

void adios_sirius_write (struct adios_file_struct * fd
                         ,struct adios_var_struct * v
                         ,const void * data
                         ,struct adios_method_struct * method
    )

{
    struct adios_sirius_data_struct * md = (struct adios_sirius_data_struct *)
        method->method_data;
    struct var_struct *var;

    int l;
    //for each variable
    //we need to create a splitdouble
    //then decide the levels
    int ndims=count_dimensions(v->dimensions);
    int type_size=adios_get_type_size(v->type,data);
    uint64_t varsize;
    uint64_t ldims[16], offsets[16], gdims[16], alldims[48];
    uint64_t nelems;
    void *splithandle = NULL;
    
    for (l=0; l < nlevels; l++)
    {
        if (alloc_var (md, l) != err_no_error)
            return;

        /*
          if(md->level[l].varcnt==0) {
          md->level[l].vars = (struct var_struct *) malloc (sizeof(struct var_struct));
          md->level[l].vars->prev=NULL;
          md->level[l].vars_head=md->level[l].vars; //assign the header of the variable list
          }
          else {
          tmp=md->level[l].vars;
          md->level[l].vars->next = (struct var_struct*) malloc (sizeof(struct var_struct));
          md->level[l].vars = md->level[l].vars->next;
          md->level[l].vars->prev=tmp;
          }

          //initial the variable structure
          init_var(var, v, ndims);
        */

        md->level[l].varcnt++;
        var = md->level[l].vars;

        //number of the dimensions of this variable

        if (!ndims)
        {
            // copy any scalar to each level?
            var->multidim=adios_flag_no;
            var->data=malloc(type_size);
            memcpy(var->data, data, type_size);

            // name the variable just like the original
            var->name = strdup (v->name);
            var->path = strdup (v->path);
            var->type = v->type;
            var->size = type_size;

            adios_common_define_var (md->level[l].grp, var->name, var->path, var->type, "", "", "");
        }
        else
        {
            //get the number of elements
            nelems = get_var_dimensions (v, ndims, gdims, ldims, offsets);
            varsize = nelems * type_size;
        
            
            //check if the type is a double
            if(v->type == adios_double)
            {
                //we can split this double up
                //check if we have already allocated the splitter
                if(splithandle == NULL)                    
                {
                    //we can split this
                    //initialize the splitter
                    splithandle = initialize_splitter_from_name(v->name, splitter_type, nelems);
                    //now feed the data into the splitter
                    split_doubles(data, splithandle);
                    //now data has been split
                }

                // name the variable
                int len = 5 + strlen (v->name);
                int splitter_type_int = get_splitter_type_from_name(splitter_type);
                var->name = (char *) malloc (len);
                sprintf(var->name, "%s/L%d",v->name, l);
                var->path = strdup (v->path);
                var->type = adios_byte;

                int i; // I guess the default settings on Titan is C89
                for(i = 0; i < ndims && i < 16; i ++)
                {
                    alldims[i] = ldims[i];                                        
                }
                for(i = 0; i < ndims && i < 16; i ++)
                {
                    alldims[i+ndims] = gdims[i];                                        
                }
                for(i = 0; i < ndims && i < 16; i ++)
                {
                    alldims[i+ndims*2] = offsets[i];                                        
                }
                
                fprintf(stderr, "splitter type = %d\n", splitter_type_int);
                adios_common_define_attribute_byvalue (md->level[l].grp, "dimensions", var->name,
                                                       adios_integer, 1, &ndims);
                adios_common_define_attribute_byvalue (md->level[l].grp, "splittype", var->name,
                                                       adios_integer, 1,
                                                       &splitter_type_int);
                adios_common_define_attribute_byvalue (md->level[l].grp, "dims", var->name, adios_long, (ndims*3),
                                                       alldims);
                                                       
                                                       
                                                       

                if(l == 0)
                {
                    //level 0 so we write out the top
                    //now get the top byte array
                    var->data = (void*)get_split_top(splithandle, nelems, &var->size);
            
                    /* FIXME: decompose here the data into multiple levels */
                    // right now just write everything at every level
                    gdims[0] = ldims[0] = var->size;
                    offsets[0] = 0;
                    var->global_dimensions = print_dimensions (1, gdims);
                    var->local_dimensions  = print_dimensions (1, ldims);
                    var->local_offsets     = print_dimensions (1, offsets);

                    // fprintf(stderr, "%s\n%s\n%s\n", var->global_dimensions,
                    //         var->local_dimensions, var->local_offsets);
                    

                    
                    /* End of decomposition code */
        
                }
                else if(l == 1)
                {
                    //now get the top byte array
                    var->data = get_split_bot(splithandle, nelems, &var->size);

                    //we write out the bottom
                    gdims[0] = ldims[0] = var->size;
                    offsets[0] = 0;
                    var->global_dimensions = print_dimensions (1, gdims);
                    var->local_dimensions  = print_dimensions (1, ldims);
                    var->local_offsets     = print_dimensions (1, offsets);
                    
                    
                }
                else
                {
                    //too many levels
                    //print some sort of error but no real need to break out here
                }
                
            }
            else
            {
                //we can't really do anything with non-doubles at the
                //moment
                //we just redefine the variable exactly as it came in
                //we store all non-split vars in the
                //level 0 file
                if(l == 0)
                {
                    //only in level 0 do we need to store this variable
                    var->name = strdup (v->name);
                    var->path = strdup (v->path);
                    var->type = v->type;
                    var->size = varsize;
                    var->data = (void*)v->data;
                    var->global_dimensions = print_dimensions (ndims, gdims);
                    var->local_dimensions  = print_dimensions (ndims, ldims);
                    var->local_offsets     = print_dimensions (ndims, offsets);
                }
                else
                {
                    //do not output to the lower layers
                    var->size = 0;
                }
            }


            
            
            if (var->size > 0)
            {
                // define the variable for this level
                adios_common_define_var (md->level[l].grp, var->name, var->path, var->type,
                                         var->local_dimensions, var->global_dimensions, var->local_offsets);
                signed char t = (char) v->type;
                adios_common_define_attribute_byvalue (md->level[l].grp, "type", var->name, adios_byte, 1, &t);
                adios_common_define_attribute_byvalue (md->level[l].grp, "level", var->name, adios_integer, 1, &l);
            }
        }

        md->level[l].totalsize += var->size;

        // write it out
        if (md->level[l].vars->size > 0)
        {
            do_write (md->level[l].fd, var->name, var->data);
        }

        /* FIXME: remove pointer to original user data since we just pointed to it to write it
         * at all levels
         */
        if (ndims > 0)
            var->data = NULL;

        
    }

    
    if(splithandle == NULL)
        free_splitter(splithandle);
        
}

void adios_sirius_read (struct adios_file_struct * fd
                        ,struct adios_var_struct * v, void * buffer
                        ,uint64_t buffer_size
                        ,struct adios_method_struct * method
    )

{
}

void adios_sirius_buffer_overflow (struct adios_file_struct * fd,
                                   struct adios_method_struct * method)
{
    struct adios_sirius_data_struct * md = (struct adios_sirius_data_struct *)
        method->method_data;
    log_error ("rank %d: SIRIUS method only works with complete buffering of data between adios_open() "
               "and adios_close(). Variables that do not fit into the buffer will not be "
               "written by this method to file %s\n", md->rank, fd->name);
}

#define FREE(v) if (v!=NULL) {free(v); v=NULL;}

void release_resource_at_close (struct adios_sirius_data_struct * md)
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

void adios_sirius_close (struct adios_file_struct * fd
                         ,struct adios_method_struct * method
    )
{
    struct adios_sirius_data_struct * md = (struct adios_sirius_data_struct *)
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
                if(no_thread == 1)
                    common_adios_close (md->level[l].fd);
                else
                    pthread_create(&md->level[l].thread,
                                   NULL, threaded_call_common_close, (void*)&md->level[l]);
            }
            if(no_thread == 0)
                for (l=0; l < nlevels; l++)
                {
                
                    pthread_join(md->level[l].thread, NULL);
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

void adios_sirius_get_write_buffer (struct adios_file_struct * fd
                                    ,struct adios_var_struct * v
                                    ,uint64_t * size
                                    ,void ** buffer
                                    ,struct adios_method_struct * method
    )
{
}

void adios_sirius_finalize (int mype, struct adios_method_struct * method)
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

void adios_sirius_end_iteration (struct adios_method_struct * method)
{
}

void adios_sirius_start_calculation (struct adios_method_struct * method)
{
}

void adios_sirius_stop_calculation (struct adios_method_struct * method)
{
}
