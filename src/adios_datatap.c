#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>
#include "adios.h"
#include "adios_transport_hooks.h"
#include "adios_bp_v1.h"
#include "adios_internals.h"

#include <sys/stat.h>
#include <sys/queue.h>

#ifdef NO_DATATAP
#if NO_DATATAP == 0
#include <ffs.h>
#if HAVE_PTL == 1
#include <thin_portal.h>
#elif HAVE_INFINIBAND == 1
#include <thin_ib.h>
#endif


#define STARTINGSIZE 16

//static LIST_HEAD(listhead, fm_structure) globallist;

#define OPLEN 4
static	char OP[OPLEN] = {'+', '-', '*', '/'};
static 	char *OP_REP[OPLEN] = {"_plus_", "_minus_", "_mult_", "_div_"};

int rank = -1;

typedef struct nametable_
{
    char *originalname;
    char *mangledname;
    LIST_ENTRY(nametable_) entries;
}nametable;

/*
static inline double getlocaltime()
{
    struct timeval t;
    gettimeofday(&t, NULL);
    double dt = tv_time_to_secs(&t);
    return dt;
}
*/

struct fm_structure
{
    FMFormatRec *format;
    int size; //in bytes - no padding
    unsigned char * buffer;
    int snd_count;
    IOhandle *s;
    FMFormat ioformat;
    LIST_HEAD(tablehead, nametable_) namelist;
};

IOhandle *s = NULL;


typedef struct datatap_method_data
{
    int opencount;
    int initialized;	
    int cycle_id;
    char *pfile;
    struct fm_structure *fm;
}dmd;

extern MPI_Comm adios_mpi_comm_world;
int initialized = 0;


static char * getFixedName(char *name)
{
    char *tempname = (char*)malloc(sizeof(char)* 255);
    tempname = strdup(name);
		
		
    char *oldname = strdup(name);
    char *loc = NULL;
		
    do
    {
	for(int i = 0; i < OPLEN; i ++)
	{
	    //checking operator OP[i]
	    loc = strchr(oldname, OP[i]);
	    if(loc == NULL)
		continue;
	    *loc = 0;
	    snprintf(tempname, 255, "%s%s%s", oldname, OP_REP[i], &loc[1]);
	    free(oldname);
	    oldname = strdup(tempname);
	}
    }while(loc != NULL);
		
	
	
    return tempname;
}


static char *findFixedName(struct fm_structure *fm, char *name)
{
    nametable *node;
	
    for(node = fm->namelist.lh_first; node != NULL; node = node->entries.le_next)
    {
	if(!strcmp(node->originalname, name))
	{
	    //matched
	    return node->mangledname;
	}
		   
    }
	
    return name;
}


extern void adios_datatap_init (const char *params, struct adios_method_struct *method)
{

    if(method->method_data != NULL)
    {
	dmd *mdata = (dmd*)method->method_data;
	if(mdata->initialized == 1)
	    return;
    }
    
    method->method_data = (void*) malloc(sizeof(struct datatap_method_data));
    dmd *mdata = (dmd*)method->method_data;
    memset(mdata, 0, sizeof(dmd));
    
    mdata->opencount = 0;
    mdata->initialized = 1;
    if(params != NULL && strlen(params) > 1)
    {
	//contains the file name of the file to read?
	mdata->pfile = strdup(params);
    }
    else
	mdata->pfile = strdup("params");

    MPI_Comm_rank (MPI_COMM_WORLD, &rank);

}



extern int adios_datatap_open (struct adios_file_struct * fd, 
				struct adios_method_struct * method)
{
    if(fd == NULL || method == NULL)
    {
	fprintf(stderr, "Bad input parameters\n");
	return 0;
    }
	
    dmd *mdata = (dmd*)method->method_data;
	
    if(mdata != NULL)
    {
	if(mdata->initialized == 0)
	{
	    fprintf(stderr, "method not initialized properly\n");
	    return 0;
	}
    }
    else
    {
	fprintf(stderr, "method not initialized\n");
	return 0;
    }

    if(mdata->fm != NULL)
    {
	return 0;
    }
    


    struct adios_group_struct *t = method->group;
    if(t == NULL)
    {
	fprintf(stderr, "group is not initialized properly\n");
	return 0;
    }
	
    struct adios_var_struct *fields = t->vars;
    if(fields == NULL)
    {
	fprintf(stderr, "adios vars not initalized properly in the group\n");
	return 0;
    }
	
    //iterate through all the types
    //create a format rec
    FMFormatRec *format = (FMFormatRec*)malloc(sizeof(FMFormatRec) * 2);
    if(format == NULL)
    {
	perror("memory allocation failed");
	return 0;
    }
		
    memset(format, 0, sizeof(FMFormatRec) *2);


    struct fm_structure *current_fm = (struct fm_structure *)malloc(sizeof(struct fm_structure));
    if(current_fm == NULL)
    {
	perror("memory allocation failed");
	return 0;
    }
		
    memset(current_fm, 0, sizeof(struct fm_structure));
	
    LIST_INIT(&current_fm->namelist);
	
	
    //associate the FMFormat rec with the fm_structure
    current_fm->format = format;
    format->format_name = strdup(t->name);

    //allocate field list
    if(t->var_count == 0)
    {
	fprintf(stderr, "no variables in this group - possibly an error\n");
	return  0;
		
    }
		
    FMFieldList field_list = (FMFieldList)malloc(sizeof(FMField) * (t->var_count+ 1));
    if(field_list == NULL)
    {
	perror("memory allocation failed");
	return 0;
    }
		
    //keep count of the total number of fields
    int fieldno = 0;

    //associate field list
    format->field_list = field_list;

    //for each type look through all the fields
    struct adios_var_struct *f;
    for(f = t->vars; f != NULL; f = f->next, fieldno++)
    {
	//make the field list
	//check name for + - * / (operators) and replace them
	char *tempname = getFixedName(f->name);
		
	if(strcmp(tempname, f->name))
	{
	    //strings don't match
	    //add to name list
	    nametable *namenode = (nametable*)malloc(sizeof(nametable));
	    namenode->originalname = strdup(f->name);
	    namenode->mangledname = strdup(tempname);
			
	    LIST_INSERT_HEAD(&current_fm->namelist, namenode, entries);
	}


	//
	field_list[fieldno].field_name = strdup(tempname);
	free(tempname);
		
	if(!f->dimensions)
	{
	    switch(f->type)
	    {
	    case adios_unknown:
		fprintf(stderr, "bad type error\n");
		fieldno --;				
		break;
				
	    case adios_integer:
		field_list[fieldno].field_type = strdup("integer");
		field_list[fieldno].field_size = sizeof(int);
		field_list[fieldno].field_offset = current_fm->size;
		current_fm->size += sizeof(int);
		break;
				
	    case adios_long:
		field_list[fieldno].field_type = strdup("integer");
		field_list[fieldno].field_size = sizeof(long long);
		field_list[fieldno].field_offset = current_fm->size;
		current_fm->size += sizeof(long long);
		break;

	    case adios_real:
		field_list[fieldno].field_type = strdup("float");
		field_list[fieldno].field_size = sizeof(float);
		field_list[fieldno].field_offset = current_fm->size;
		current_fm->size += sizeof(float);
		break;
				
	    case adios_string:
		field_list[fieldno].field_type = strdup("string");
		field_list[fieldno].field_size = sizeof(char);
		field_list[fieldno].field_offset = current_fm->size;
		current_fm->size += sizeof(unsigned char*);
		break;
				
	    case adios_double:
		field_list[fieldno].field_type = strdup("float");
		field_list[fieldno].field_size = sizeof(double);
		field_list[fieldno].field_offset = current_fm->size;
		current_fm->size += sizeof(double);
		break;
				

	    case adios_byte:
		field_list[fieldno].field_type = strdup("char");
		field_list[fieldno].field_size = sizeof(char);
		field_list[fieldno].field_offset = current_fm->size;
		current_fm->size += sizeof(char);
		break;
				
	    default:
		fprintf(stderr, "unknown type error %d\n", f->type);
		fieldno--;
		break;
	    }
	}
	else
	{
	    //its a vector!
	    //find out the dimensions by walking the dimension list
	    struct adios_dimension_struct *d = f->dimensions;
#define DIMSIZE 10240
	    char dims[DIMSIZE] = {0};
#define ELSIZE 256
	    char el[ELSIZE] ={0};
		
	    struct adios_var_struct *var = NULL;
	    
				
	    //create the dimension thingy
	    for(; d != NULL;d = d->next)
	    {
		//for each dimension just take the upper_bound
		var = adios_find_var_by_id(fields, d->dimension.id);
		
		if(var)
		{
		    //findFixedName returns the mangled name from the original name
		    char *name = findFixedName(current_fm, var->name);
		    snprintf(el,ELSIZE,"[%s]", name);
		}
		else //its a number
		    snprintf(el,ELSIZE, "[%d]", d->dimension.rank);

		strncat(dims,el,DIMSIZE);
	    }
						
	    switch(f->type)
	    {
	    case adios_unknown:
		fprintf(stderr, "bad type error\n");
		fieldno --;				
		break;
				
	    case adios_integer:
		field_list[fieldno].field_type = (char*)malloc(sizeof(char)*255);
		snprintf((char*)field_list[fieldno].field_type, 255,
			 "integer%s", dims);
		field_list[fieldno].field_size = sizeof(int);
		field_list[fieldno].field_offset = current_fm->size;
		current_fm->size += sizeof(int*);
		break;

	    case adios_long:
		field_list[fieldno].field_type = (char*)malloc(sizeof(char)*255);
		snprintf((char*)field_list[fieldno].field_type, 255,
			 "integer%s", dims);
		field_list[fieldno].field_size = sizeof(long long);
		field_list[fieldno].field_offset = current_fm->size;
		current_fm->size += sizeof(int*);
		break;

				
	    case adios_real:
		field_list[fieldno].field_type = (char*)malloc(sizeof(char)*255);
		snprintf((char*)field_list[fieldno].field_type, 255,
			 "float%s", dims);
		field_list[fieldno].field_size = sizeof(float);
		field_list[fieldno].field_offset = current_fm->size;
		current_fm->size += sizeof(int*);
		break;
				
	    case adios_string:
		field_list[fieldno].field_type = strdup("string");
		field_list[fieldno].field_size = sizeof(char);
		field_list[fieldno].field_offset = current_fm->size;
		current_fm->size += sizeof(unsigned char*);
		break;
				
	    case adios_double:
		field_list[fieldno].field_type = (char*)malloc(sizeof(char)*255);
		snprintf((char*)field_list[fieldno].field_type, 255,
			 "float%s", dims);
		field_list[fieldno].field_size = sizeof(double);
		field_list[fieldno].field_offset = current_fm->size;
		current_fm->size += sizeof(double*);
		break;

	    case adios_byte:
		field_list[fieldno].field_type = (char*)malloc(sizeof(char)*255);
		snprintf((char*)field_list[fieldno].field_type, 255,
			 "char%s", dims);
		field_list[fieldno].field_size = sizeof(char);
		field_list[fieldno].field_offset = current_fm->size;
		current_fm->size += sizeof(char*);
		break;
				
	    default:
		fprintf(stderr, "unknown type error\n", f->type);
		fieldno--;
		break;
	    }
				
	}
    }
    //terminate the the fieldlist
    for(; fieldno < (t->var_count + 1); fieldno++)
    {
	field_list[fieldno].field_type = NULL;
	field_list[fieldno].field_name = NULL;
	field_list[fieldno].field_offset = 0;
	field_list[fieldno].field_size = 0;
    }

    current_fm->format->struct_size = current_fm->size;
    current_fm->buffer = (unsigned char*)malloc(current_fm->size);
    memset(current_fm->buffer, 0, current_fm->size);


//    current_fm->s = s;

#if HAVE_PTL == 1
//defined(__CRAYXT_COMPUTE_LINUX_TARGET)
//fprintf(stderr, "im here rank = %d %s %d\n", rank, __FILE__, __LINE__);
    current_fm->s = InitIOFromFile("param", rank);
    current_fm->s->rank = rank;
    current_fm->ioformat = register_data(current_fm->s, current_fm->format);
#elif HAVE_INFINIBAND == 1
    current_fm->s = EVthin_ib_InitIOFile("param", 1);
    current_fm->ioformat = EVthin_ib_registerData(current_fm->s, current_fm->format);
#endif

    current_fm->snd_count = 0;

    mdata->fm = current_fm;

    return 1;
    
}

static FMField* internal_find_field(char *name, FMFieldList flist)
{
    FMField *f = flist;
    while(f->field_name != NULL && strcmp(f->field_name, name))
    {
	f++;
    }
	
    return f;
}


extern void adios_datatap_write(struct adios_file_struct *fd,
				struct adios_var_struct *f,
				void *data,
				struct adios_method_struct *method)
{
    struct fm_structure *fm;
    dmd *mdata = (dmd*)method->method_data;
    
		
    struct adios_group_struct *group = method->group;
    
    fm = mdata->fm;
    
    if(group == NULL || fm == NULL)
    {
	fprintf(stderr, "group or fm is null - improperly initialized\n");
	return;
	
    }
    

    FMFieldList flist = fm->format->field_list;
    FMField *field = NULL;

    char *fixedname = findFixedName(fm, f->name);
	
    field = internal_find_field(fixedname, flist);
    if(field != NULL)
    {
	if(!f->dimensions)
	{
	    //scalar quantity
	    if(data)
	    {
		//why wouldn't it have data?
		memcpy(&fm->buffer[field->field_offset], data, 
		       field->field_size);

	    }
	    else
	    {
		fprintf(stderr, "no data for  scalar %s\n", f->name);
						
	    }
					
					
	}
	else
	{
	    //vector quantity
	    if(data)
	    {
		//we just need to copy the pointer stored in f->data
		memcpy(&fm->buffer[field->field_offset], &data,
		       sizeof(void*));

	    }
	    else
	    {
		fprintf(stderr, "no data for vector %s\n", f->name);	
	    }
	}				
    }
}

static void internal_adios_datatap_write (struct adios_file_struct * fd, 
					  struct adios_method_struct *method);

extern void adios_datatap_close(struct adios_file_struct *fd,
	struct adios_method_struct *method)
{

    dmd *mdata = method->method_data;
    
    
    if(!mdata->initialized)
    {
	return;
    }
	

    if(fd->mode & adios_mode_write)
    {
	internal_adios_datatap_write(fd, method);
    }
    
}


void internal_adios_datatap_write (struct adios_file_struct * fd, struct adios_method_struct *method)
{

	
    if(fd == NULL)
    {
	fprintf(stderr, "fd is null\n");
		
	return;
    }

    dmd *mdata = method->method_data;
    
	
    struct adios_group_struct *t = method->group;
	
    

    //initialize the globallist
	
    //iterate through all the types
    //find the correct format by name
    struct fm_structure *fm = mdata->fm;
		
    if(t == NULL || fm == NULL)
    {
	fprintf(stderr, "improperly initialized for write\n");
	
	return;
	
    }
    
//#if defined(__CRAYXT_COMPUTE_LINUX_TARGET)
#if HAVE_PTL == 1
    if(fm->snd_count > 0)
    {
	if(rank == 0)
	    fprintf(stderr, "Ending send\n");
	send_end(fm->s);
	fm->snd_count --;
	
    }
#elif HAVE_INFINIBAND == 1
    if(fm->snd_count > 0)
    {
	EVthin_ib_endSend(fm->s);
	fm->snd_count --;

    }
#endif
    

    //now that we have all the info lets call write on each of these data types
//#if defined(__CRAYXT_COMPUTE_LINUX_TARGET)
#if HAVE_PTL == 1
    start_send(fm->s, fm->buffer, fm->size, fm->ioformat, NULL);

#elif HAVE_INFINIBAND == 1
    EVthin_ib_startSend(fm->s, fm->buffer, fm->size, fm->ioformat);
#endif

    fm->snd_count ++;

//  		FMContext src_context = create_local_FMcontext(NULL);
// 		FMFormat ioformat = register_data_format(src_context, fm->format);
// 		int size =0;
// 		FFSBuffer encode_buffer = create_FFSBuffer();		
// 		char *xfer = FFSencode(encode_buffer, ioformat, fm->buffer, &size);
//  		FFSFile file = open_FFSfile("temp", "w");
//  		write_FFSfile(file, ioformat, fm->buffer);
//  		close_FFSfile(file);

}

extern void adios_datatap_finalize (int mype, struct adios_method_struct *method)
{

    dmd *mdata = method->method_data;
    
	
    struct adios_group_struct *t = method->group;
	
    

    //initialize the globallist
	
    //iterate through all the types
    //find the correct format by name
    struct fm_structure *fm = mdata->fm;
		
    if(t == NULL || fm == NULL)
    {
	fprintf(stderr, "improperly initialized for finalize %p %p\n", t, fm);
	
	return;
	
    }
    
    
    
#if HAVE_PTL == 1
    if(fm->snd_count > 0)
    {
	send_end(fm->s);
    }
#elif HAVE_INFINIBAND == 1
    if(fm->snd_count > 0)
    {
	EVthin_ib_endSend(fm->s);
    }
#endif

    char buffer[255];
    mkdir("client", S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    snprintf(buffer, 255, "client/client-%d", fm->s->rank);

    outputTimingInfo(buffer);
}


extern void adios_datatap_end_iteration (struct adios_method_struct *method)
{
    struct fm_structure *fm;
    dmd *mdata = (dmd*)method->method_data;
    fm = mdata->fm;
    
    if(fm == NULL)
	return;
    
	
#if HAVE_PTL == 1
    startIter(fm->s);
#elif HAVE_INFINIBAND == 1
    EVthin_ib_startIter(fm->s);
#endif
	
    mdata->cycle_id = 0;
    
}

extern void adios_datatap_start_calculation (struct adios_method_struct *method)
{
    struct fm_structure *fm;
    dmd *mdata = (dmd*)method->method_data;
    fm = mdata->fm;
	
    if(fm == NULL)
	return;
 

#if HAVE_PTL == 1
    startCompute(fm->s, mdata->cycle_id);
#elif HAVE_INFINIBAND == 1
    EVthin_ib_startCompute(fm->s, mdata->cycle_id);	
#endif
	
    
}

extern void adios_datatap_stop_calculation (struct adios_method_struct * method)
{
    struct fm_structure *fm;
    dmd *mdata = (dmd*)method->method_data;
    fm = mdata->fm;

    if(fm == NULL)
	return;
 
	
#if HAVE_PTL == 1
    endCompute(fm->s, mdata->cycle_id);
#elif HAVE_INFINIBAND == 1
    EVthin_ib_endCompute(fm->s, mdata->cycle_id);
#endif

    mdata->cycle_id++;
    
}

extern void adios_datatap_get_write_buffer (struct adios_file_struct * fd,
					    struct adios_var_struct * f,
					    uint64_t * size,
					    void ** buffer,
					    struct adios_method_struct *method)
{
    fprintf (stderr, "adios_datatap_write_get_buffer: datatap disabled, "
	     "no portals support\n");
}

void adios_datatap_read (struct adios_file_struct * fd,
			 struct adios_var_struct * f,
			 void * buffer,
			 uint64_t buffer_size,
			 struct adios_method_struct *method)
{

}

extern enum ADIOS_FLAG adios_datatap_should_buffer(struct adios_file_struct * fd, 
					struct adios_method_struct * method, void * comm) 
{
    return adios_flag_no;
}					\

#else


void adios_datatap_init (const char * parameters  
                      ,struct adios_method_struct * method  
                      ) {}  
int adios_datatap_open (struct adios_file_struct * fd  
                     ,struct adios_method_struct * method  
                     ) {return 0;}  
enum ADIOS_FLAG adios_datatap_should_buffer (struct adios_file_struct * fd  
                                          ,struct adios_method_struct * method  
                                          ,void * comm  
                                          ) {return 0;}  
void adios_datatap_write (struct adios_file_struct * fd  
                       ,struct adios_var_struct * v  
                       ,void * data  
                       ,struct adios_method_struct * method  
                       ) {}  
void adios_datatap_get_write_buffer (struct adios_file_struct * fd  
                                  ,struct adios_var_struct * v  
                                  ,uint64_t * size  
                                  ,void ** buffer  
                                  ,struct adios_method_struct * method  
                                  ) {}  
void adios_datatap_read (struct adios_file_struct * fd  
                      ,struct adios_var_struct * v  
                      ,void * buffer  
                      ,uint64_t buffer_size  
                      ,struct adios_method_struct * method  
                      ) {}  
void adios_datatap_close (struct adios_file_struct * fd  
                       ,struct adios_method_struct * method  
                       ) {}  
void adios_datatap_finalize (int mype, struct adios_method_struct * method) {}  
void adios_datatap_end_iteration (struct adios_method_struct * method) {}  
void adios_datatap_start_calculation (struct adios_method_struct * method) {}  
void adios_datatap_stop_calculation (struct adios_method_struct * method) {}

#endif
#else


void adios_datatap_init (const char * parameters  
                      ,struct adios_method_struct * method  
                      ) {}  
int adios_datatap_open (struct adios_file_struct * fd  
                     ,struct adios_method_struct * method  
                     ) {return 0;}  
enum ADIOS_FLAG adios_datatap_should_buffer (struct adios_file_struct * fd  
                                          ,struct adios_method_struct * method  
                                          ,void * comm  
                                          ) {return adios_flag_unknown;}  
void adios_datatap_write (struct adios_file_struct * fd  
                       ,struct adios_var_struct * v  
                       ,void * data  
                       ,struct adios_method_struct * method  
                       ) {}  
void adios_datatap_get_write_buffer (struct adios_file_struct * fd  
                                  ,struct adios_var_struct * v  
                                  ,uint64_t * size  
                                  ,void ** buffer  
                                  ,struct adios_method_struct * method  
                                  ) {}  
void adios_datatap_read (struct adios_file_struct * fd  
                      ,struct adios_var_struct * v  
                      ,void * buffer  
                      ,uint64_t buffer_size  
                      ,struct adios_method_struct * method  
                      ) {}  
void adios_datatap_close (struct adios_file_struct * fd  
                       ,struct adios_method_struct * method  
                       ) {}  
void adios_datatap_finalize (int mype, struct adios_method_struct * method) {}  
void adios_datatap_end_iteration (struct adios_method_struct * method) {}  
void adios_datatap_start_calculation (struct adios_method_struct * method) {}  
void adios_datatap_stop_calculation (struct adios_method_struct * method) {}

#endif 

extern void adios_empty_init (const char * parameters  
                      ,struct adios_method_struct * method  
                      ) {}  
extern void adios_empty_open (struct adios_file_struct * fd  
                      ,struct adios_method_struct * method  
                      ) {}  
extern void adios_empty_write (struct adios_file_struct * fd  
                       ,struct adios_var_struct * v  
                       ,void * data  
                       ,struct adios_method_struct * method  
                       ) {}  
extern void adios_empty_get_write_buffer (struct adios_file_struct * fd  
                                  ,struct adios_var_struct * v  
                                  ,unsigned long long * size  
                                  ,void ** buffer  
                                  ,struct adios_method_struct * method  
                                  ) {}  
extern void adios_empty_read (struct adios_file_struct * fd  
                      ,struct adios_var_struct * v  
                      ,void * buffer  
                      ,struct adios_method_struct * method  
                      ) {}  
extern void adios_empty_close (struct adios_file_struct * fd  
                       ,struct adios_method_struct * method  
                       ) {}  
extern void adios_empty_finalize (int mype, struct adios_method_struct * method) {}  
extern void adios_empty_end_iteration (struct adios_method_struct * method) {}  
extern void adios_empty_start_calculation (struct adios_method_struct * method) {}  
extern void adios_empty_stop_calculation (struct adios_method_struct * method) {}




FILE *prof_fp;
struct timeval prof_tm;

/* This should be called only after adios_init */
void prof_init_()
{
    int my_rank = 0, mysize= 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &mysize);

    #include <stdlib.h>
    
    srand(my_rank);
    srand48(my_rank);

    int tenpercent = (0.01 *(double)mysize);

    int rand1 = lrand48() % tenpercent;
    
    int rand2 = lrand48()% mysize;
    
    
    if (rand2 > rand1 && my_rank != 0)
    {
	prof_fp = NULL;
	return;
    }

    char prof_file[128];
    sprintf(prof_file,"profile/prof.%d", my_rank);
    prof_fp = fopen(prof_file, "w");
    if(!prof_fp){
	fprintf(stderr, "Unable to open profiling file: %s\n", prof_file);
	prof_fp = stdout;
    }
    else {
	fseek(prof_fp, 0, SEEK_END);
    }
	
}

void prof_finish_()
{
    if(prof_fp == NULL)
	return;
    
    if(fclose(prof_fp)){
	fprintf(stderr, "Unable to close profiling file\n");		
    }
	
}

const char *prof_type[3] = {"COMP", "COMM", "IO"};

/* Assuming null terminated strings as of now! */
void prof_evt_start_(const char *event, int *type)
{
    if(prof_fp == NULL)
	return;
    gettimeofday(&prof_tm, 0);
    fprintf(prof_fp, "%s:S \t%.0f \t%s\n", event,
	    prof_tm.tv_sec * 1e6 + prof_tm.tv_usec,
	    prof_type[*type]);
}

void prof_evt_end_(const char *event)
{
    if(prof_fp == NULL)
	return;
    gettimeofday(&prof_tm, 0);
    fprintf(prof_fp, "%s:E \t%.0f\n", event, 
	    prof_tm.tv_sec * 1e6 + prof_tm.tv_usec);
}
