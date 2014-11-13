//passed the test for 1024 cores and level-3 spatial aggregation

#include <unistd.h>
#include <fcntl.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <assert.h>

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

//transport to store arrays as internal and boundary blocks
//boundary width is a parameter passed to init
//when an array is passed to write we will split it into two parts
//1. internal array. Same number of dimensions but all the dimensions are reduced by
//2*w. The offset is moved by w in each dimension
//2. 

#define SFLAGS S_IRUSR | S_IWUSR | S_IRGRP | S_IWGRP | S_IROTH | S_IWOTH
#define WFLAGS O_WRONLY | O_CREAT | O_TRUNC
#define RFLAGS O_RDONLY
#define AFLAGS O_WRONLY | O_CREAT

typedef struct adios_regen_data_struct
{
	struct adios_bp_buffer_struct_v1 b;

	int bwidth; //number of boundary cells
	int lasttime; //last timstep with a full checkpoint

	int btime;
	int bcount;
	int btotal;

	uint64_t vars_start, vars_header_size;
	uint64_t b_start, b_header_size;

	struct adios_group_struct *group;
	
}regendata;

typedef struct adios_regen_boundary_struct
{
	struct adios_var_struct *parent;
	struct adios_var_struct *self;
	uint32_t parent_time;
	uint32_t current_time;
}bvar;

extern struct adios_transport_struct * adios_transports;
static short adios_regen_initialized = 0;

//helper functions
static int internal_regen_write_open(struct adios_file_struct * fd,
                                     regendata *r,
                                     struct adios_method_struct * method,
                                     MPI_Comm comm);
static int internal_regen_read_open(struct adios_file_struct * fd,
                                     regendata *r,
                                     struct adios_method_struct * method,
                                     MPI_Comm comm);
	
                                     
void adios_regen_init (const PairStruct * parameters, 
                       struct adios_method_struct * method)
{
	regendata *r = NULL;
	const PairStruct *p = parameters;
	int i  = 0;
	
	if(adios_regen_initialized != 0)
	{
		return;
	}
	adios_regen_initialized = 1;


	while(p)
	{
		if(!strcasecmp(p->name, "btime"))
		{
			r->btime = strtol(p->value, NULL, 10);
		}
		p=p->next;
	}

	r = (regendata*)malloc(sizeof(regendata));
	method->method_data = r;

	adios_buffer_struct_init(&r->b);
	
	// r->boundary = (struct adios_bp_buffer_struct_v1*)
	// 	malloc(sizeof(struct adios_bp_buffer_struct_v1)*r->btime);
	
	// for(i = 0; i < r->btime; i ++)
	// 	adios_buffer_struct_init(&r->boundary[i]);
	r->btime = 100; //by default we keep 100 boundaries in memory

	r->vars_start = r->vars_header_size = r->bcount = r->btotal = 0;
	r->b_start = r->b_header_size = 0;

	log_debug("Regen initialized with btime = %d\n", r->btime);
	                       
}
int adios_regen_open (struct adios_file_struct * fd,
                      struct adios_method_struct * method,
                      MPI_Comm comm)
{
	regendata *r = (regendata*) method->method_data;
	char *name = NULL;
    struct stat s;

    //allocate the name for the file
	name = malloc(strlen (method->base_path) + strlen (fd->name) + 1);
	sprintf (name, "%s%s", method->base_path, fd->name);

	//check if exist and what the size of the file is
	if (stat (name, &s) == 0)
	{
		//same for both data structures?
        r->b.file_size = s.st_size;
	}        

	
	if(fd->mode == adios_mode_write)
	{
		//open file for writing
		r->b.f = open(name, WFLAGS, SFLAGS);
		if(r->b.f == -1)
		{
			fprintf (stderr, "adios_posix1_open failed for "
			         "base_path %s, name %s\n",
			         method->base_path, fd->name);

			free (name);
			return 0;			
		}

		fd->base_offset = 0;
		fd->pg_start_in_file = 0;
		
		internal_regen_write_open(fd, r, method, comm);
	}
	else if(fd->mode == adios_mode_read)
	{
		internal_regen_read_open(fd, r, method, comm);
	}
	else if(fd->mode == adios_mode_append)
	{
		internal_regen_write_open(fd, r, method, comm);
	}
	
}
enum ADIOS_FLAG adios_regen_should_buffer (struct adios_file_struct * fd,
                                           struct adios_method_struct * method)
{
	return adios_flag_no; //no buffering in common
}

void adios_regen_write (struct adios_file_struct * fd, 
                        struct adios_var_struct * v, 
                        void * data, 
                        struct adios_method_struct * method)
{
	
}
void adios_regen_get_write_buffer (struct adios_file_struct * fd, 
                                   struct adios_var_struct * v, 
                                   uint64_t * size, 
                                   void ** buffer, 
                                   struct adios_method_struct * method)
{
	
}

void adios_regen_read (struct adios_file_struct * fd, 
                       struct adios_var_struct * v, 
                       void * buffer, 
                       uint64_t buffer_size, 
                      struct adios_method_struct * method)
{
	
}

void adios_regen_close (struct adios_file_struct * fd,
                        struct adios_method_struct * method)
{
	
}

void adios_regen_finalize (int mype, struct adios_method_struct * method)
{
	
}

void adios_regen_end_iteration (struct adios_method_struct * method)
{
}

void adios_regen_start_calculation (struct adios_method_struct * method)
{
}
void adios_regen_stop_calculation (struct adios_method_struct * method)
{
}

static char * dup_path (const char *path)
{
    char * p = NULL;
    int len;
    if (!path)
        return strdup("");
    len = strlen (path);
    /* remove trailing / characters */
    while (len > 1 && path[len-1] == '/') {
        /* ends with '/' and it is not a single '/' */
        len--;
    }
    p = malloc (len+1);
    if (!p)
        return NULL;
    strncpy (p, path, len);
    p[len] = '\0';
    return p;
}


static int internal_regen_write_open(struct adios_file_struct * fd,
                                     regendata *r,
                                     struct adios_method_struct * method,
                                     MPI_Comm comm)
{
	struct adios_group_struct *group = method->group;
	struct adios_var_struct *v;
	char *tempath;
	int len;
	char *name = "boundary";

	v = group->vars;
	while(v)
	{
		if(v->dimensions)//variable has dimensions
		{
			if(v->name && v->path)
			{
				fprintf(stderr, "v->name = %s\t", v->name);
				fprintf(stderr, "v->path = %s\n", v->path);

				//now create a new path from v->path/v->name/
				len = strlen(v->name) + strlen(v->path);
				tempath=(char*)malloc(sizeof(char)*(len+3));
				snprintf(tempath, len+3, "/%s/%s", v->name, v->path);
				free(tempath);
			}
		}

		v=v->next;
	}
	
	return 0;
}

static int internal_regen_read_open(struct adios_file_struct * fd,
                                     regendata *r,
                                     struct adios_method_struct * method,
                                     MPI_Comm comm)
{
	return -1;
}
