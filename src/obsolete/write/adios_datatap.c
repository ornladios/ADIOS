#include "config.h"

#if NO_DATATAP == 0

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ffs.h>
#include <atl.h>
#include "public/adios.h"
#include "core/adios_internals.h"
#include "core/adios_transport_hooks.h"
#include "core/util.h"

#include <sys/queue.h>
#if HAVE_PORTALS == 1
#include <thin_portal.h>
#elif HAVE_INFINIBAND == 1
#include <thin_ib.h>
#endif


#define STARTINGSIZE 16

//static LIST_HEAD(listhead, fm_structure) globallist;

#define OPLEN 4
static char OP[OPLEN] = { '+', '-', '*', '/' };
static char *OP_REP[OPLEN] = { "_plus_", "_minus_", "_mult_", "_div_" };

int rank = -1;

typedef struct nametable_
{
    char *originalname;
    char *mangledname;
    LIST_ENTRY(nametable_) entries;
} nametable;



typedef struct altname_
{
    char *name;
    FMField *field;                                               //contains the field list
    LIST_ENTRY(altname_) entries;
} altname;

typedef struct dimnames_
{
    char *name;
    LIST_HEAD(alts, altname_) altlist;
    LIST_ENTRY(dimnames_) entries;
} dimnames;


int mpisize  = 1;

struct fm_structure
{
    FMFormatRec *format;
    int size;                                                   //in bytes - no padding
    unsigned char *buffer;
    int snd_count;
    IOhandle *s;
    FMFormat ioformat;
    attr_list alist;    

    LIST_HEAD(tablehead, nametable_) namelist;
    LIST_HEAD(dims, dimnames_) dimlist;
};


typedef struct datatap_method_data
{
    int opencount;
    int initialized;
    int cycle_id;
    char *pfile;
    struct fm_structure *fm;
} dmd;

MPI_Comm adios_mpi_comm_world = MPI_COMM_WORLD;
int initialized = 0;


static char *
getFixedName(char *name)
{
    char *tempname = (char *) malloc(sizeof(char) * 255);
    tempname = strdup(name);


    char *oldname = strdup(name);
    char *loc = NULL;
    int i;

    do
    {
        for (i = 0; i < OPLEN; i++)
        {
    //checking operator OP[i]
            loc = strchr(oldname, OP[i]);
            if (loc == NULL)
                continue;
            *loc = 0;
            snprintf(tempname, 255, "%s%s%s", oldname, OP_REP[i], &loc[1]);
            free(oldname);
            oldname = strdup(tempname);
        }
    }
    while (loc != NULL);



    return tempname;
}



//return a list of all the names associated with the variable
static char **
getAltName(char *variable, int *count)
{
    if (count == NULL)
        return NULL;

    return NULL;

}


static char *
findFixedName(struct fm_structure *fm, char *name)
{
    nametable *node;

    for (node = fm->namelist.lh_first; node != NULL;
         node = node->entries.le_next)
    {
        if (!strcmp(node->originalname, name))
        {
    //matched
            return node->mangledname;
        }

    }

    return name;
}


extern void
adios_datatap_init(const PairStruct *params, struct adios_method_struct *method)
{
    if (method->method_data != NULL)
    {
        dmd *mdata = (dmd *) method->method_data;
        if (mdata->initialized == 1)
            return;
    }

    method->method_data = (void *) malloc(sizeof(struct datatap_method_data));
    dmd *mdata = (dmd *) method->method_data;
    memset(mdata, 0, sizeof(dmd));

    mdata->opencount = 0;
    mdata->initialized = 1;
    if (md->parameters != NULL && strlen(md->parameters) > 1)
    {
//contains the file name of the file to read?
        mdata->pfile = strdup(md->parameters);
    }
    else
        mdata->pfile = strdup("params");

    MPI_Comm_rank(adios_mpi_comm_world, &rank);
    MPI_Comm_size(adios_mpi_comm_world, &mpisize);

    fflush(stderr);

}

static altname *
findAltName(struct fm_structure *current_fm, char *dimname, char *varname)
{
    int len = strlen(dimname) + strlen(varname) + 2;
    char *aname = (char *) malloc(sizeof(char) * len);
    strcpy(aname, dimname);
    strcat(aname, "_");
    strcat(aname, varname);
    dimnames *d = NULL;

    for (d = current_fm->dimlist.lh_first; d != NULL; d = d->entries.le_next)
    {
        if (!strcmp(d->name, dimname))
        {
    //matched
            break;
        }

    }

    if (d == NULL)
    {
        d = (dimnames *) malloc(sizeof(dimnames));
        d->name = dimname;
        LIST_INIT(&d->altlist);
        LIST_INSERT_HEAD(&current_fm->dimlist, d, entries);
    }

    FMField *field = (FMField *) malloc(sizeof(FMField));

    altname *a = (altname *) malloc(sizeof(altname));
    a->name = aname;
    a->field = field;
    field->field_name = strdup(aname);
    field->field_type = strdup("integer");
    field->field_size = sizeof(int);
    field->field_offset = -1;


    LIST_INSERT_HEAD(&d->altlist, a, entries);

    return a;
}



extern int
adios_datatap_open(struct adios_file_struct *fd,
           struct adios_method_struct *method, void*comm)
{
    if (fd == NULL || method == NULL)
    {
        fprintf(stderr, "Bad input parameters\n");
        return -1;
    }

    dmd *mdata = (dmd *) method->method_data;

    if (mdata != NULL)
    {
        if (mdata->initialized == 0)
        {
            fprintf(stderr, "method not initialized properly\n");
            return;
        }
    }
    else
    {
        fprintf(stderr, "method not initialized\n");
        return -1;
    }

    if (mdata->fm != NULL)
    {
        return -1;
    }



    struct adios_group_struct *t = method->group;
    if (t == NULL)
    {
        fprintf(stderr, "group is not initialized properly\n");
        return -1;
    }

    struct adios_var_struct *fields = t->vars;
    if (fields == NULL)
    {
        fprintf(stderr, "adios vars not initalized properly in the group\n");
        return -1;
    }

//iterate through all the types
//create a format rec
    FMFormatRec *format = (FMFormatRec *) malloc(sizeof(FMFormatRec) * 2);
    if (format == NULL)
    {
        perror("memory allocation failed");
        return -1;
    }

    memset(format, 0, sizeof(FMFormatRec) * 2);


    struct fm_structure *current_fm =
        (struct fm_structure *) malloc(sizeof(struct fm_structure));
    if (current_fm == NULL)
    {
        perror("memory allocation failed");
        return -1;
    }

    memset(current_fm, 0, sizeof(struct fm_structure));
    current_fm->alist = create_attr_list();
    set_int_attr(current_fm->alist, attr_atom_from_string("mpisize"), mpisize);
    

    LIST_INIT(&current_fm->namelist);
    LIST_INIT(&current_fm->dimlist);



//associate the FMFormat rec with the fm_structure
    current_fm->format = format;
    format->format_name = strdup(t->name);

//allocate field list
    if (t->var_count == 0)
    {
        fprintf(stderr, "no variables in this group - possibly an error\n");
        return -1;

    }

    int altvarcount = 0;

    FMFieldList field_list =
        (FMFieldList) malloc(sizeof(FMField) * (t->var_count + 1));
    if (field_list == NULL)
    {
        perror("memory allocation failed");
        return -1;
    }

//keep count of the total number of fields
    int fieldno = 0;


//for each type look through all the fields
    struct adios_var_struct *f;
    for (f = t->vars; f != NULL; f = f->next, fieldno++)
    {
//make the field list
//check name for + - * / (operators) and replace them
        char *tempname = getFixedName(f->name);


        if (strcmp(tempname, f->name))
        {
    //strings don't match
    //add to name list
            nametable *namenode = (nametable *) malloc(sizeof(nametable));
            namenode->originalname = strdup(f->name);
            namenode->mangledname = strdup(tempname);

            LIST_INSERT_HEAD(&current_fm->namelist, namenode, entries);
        }


//
        field_list[fieldno].field_name = strdup(tempname);
        free(tempname);

        if (!f->dimensions)
        {
            switch (f->type)
            {
            case adios_unknown:
                fprintf(stderr, "bad type error\n");
                fieldno--;
                break;

            case adios_integer:
                field_list[fieldno].field_type = strdup("integer");
                field_list[fieldno].field_size = sizeof(int);
                field_list[fieldno].field_offset = current_fm->size;
                current_fm->size += sizeof(int);
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
                current_fm->size += sizeof(unsigned char *);
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
            char dims[DIMSIZE] = { 0 };
#define ELSIZE 256
            char el[ELSIZE] = { 0 };


    //create the dimension thingy
            for (; d != NULL; d = d->next)
            {
                //for each dimension just take the upper_bound
                if (d->dimension.id)
                {
                    //findFixedName returns the mangled name from the original name
                    struct adios_var_struct *tmp_var = adios_find_var_by_id(t->vars, d->dimension.id);
                    char *name =
                    findFixedName(current_fm, 
                              tmp_var->name);
                    //create the alternate name for this variable and the array its defining
                    altname *a = findAltName(current_fm, name,
                                 (char*)field_list[fieldno].field_name);
                    //altname is a new variable that we need to add to the field list
                    altvarcount++;

                    snprintf(el, ELSIZE, "[%s]", a->name);
//                    fprintf(stderr, "%s\t", el);
                    //offset_increment is just the size of the pointer

                }
                else            //its a number
                {
                    //if its a number the offset_increment will be the size of the variable*rank
                    snprintf(el, ELSIZE, "[%d]", d->dimension.rank);
                }
                strncat(dims, el, DIMSIZE);
            }
//             fprintf(stderr, "%s\n", dims);

            while(current_fm->size % 8 != 0)
            {
                current_fm->size ++;                    
            }

            switch (f->type)
            {
            case adios_unknown:
                fprintf(stderr, "bad type error\n");
                fieldno--;
                break;

            case adios_integer:
                field_list[fieldno].field_type =
                    (char *) malloc(sizeof(char) * 255);
                snprintf((char *) field_list[fieldno].field_type, 255,
                         "integer%s", dims);
                field_list[fieldno].field_size = sizeof(int);
                
                field_list[fieldno].field_offset = current_fm->size;
                current_fm->size += sizeof(void *);
                break;

            case adios_real:
                field_list[fieldno].field_type =
                    (char *) malloc(sizeof(char) * 255);
                snprintf((char *) field_list[fieldno].field_type, 255,
                         "float%s", dims);
                field_list[fieldno].field_size = sizeof(float);
                field_list[fieldno].field_offset = current_fm->size;
                current_fm->size += sizeof(void *);
                break;

            case adios_string:
                field_list[fieldno].field_type = strdup("string");
                field_list[fieldno].field_size = sizeof(char);
                field_list[fieldno].field_offset = current_fm->size;
                current_fm->size += sizeof(void *);
                break;

            case adios_double:
                field_list[fieldno].field_type =
                    (char *) malloc(sizeof(char) * 255);
                snprintf((char *) field_list[fieldno].field_type, 255,
                         "float%s", dims);
                field_list[fieldno].field_size = sizeof(double);
                field_list[fieldno].field_offset = current_fm->size;
                current_fm->size += sizeof(void *);
                break;

            case adios_byte:
                field_list[fieldno].field_type =
                    (char *) malloc(sizeof(char) * 255);
                snprintf((char *) field_list[fieldno].field_type, 255, "char%s",
                         dims);
                field_list[fieldno].field_size = sizeof(char);
                field_list[fieldno].field_offset = current_fm->size;
                current_fm->size += sizeof(void *);
                break;

            default:
                fprintf(stderr, "unknown type error %d\n", f->type);
                fieldno--;
                break;
            }

        }

//      fprintf(formatfile, "%s, %s, %d, %d\n", field_list[fieldno].field_name, field_list[fieldno].field_type,field_list[fieldno].field_size,field_list[fieldno].field_offset); 


    }

    dimnames *d = NULL;
//     fprintf(stderr,
//             "============================\n\tdims and alts\n============================\n");
//     fprintf(stderr, "altvarcount = %d\n", altvarcount);
    field_list =
        (FMFieldList) realloc(field_list,
                              sizeof(FMField) * (altvarcount + t->var_count +
                                                 1));

    for (d = current_fm->dimlist.lh_first; d != NULL; d = d->entries.le_next)
    {
//         fprintf(stderr, "%s\t", d->name);
        altname *a = NULL;
        for (a = d->altlist.lh_first; a != NULL; a = a->entries.le_next)
        {
//             fprintf(stderr, "%s\t", a->name);
            a->field->field_offset = current_fm->size;
            current_fm->size += sizeof(int);
            memcpy(&field_list[fieldno], a->field, sizeof(FMField));
            fieldno++;

        }
//         fprintf(stderr, "\n");

    }

//     fprintf(stderr,
//             "============================\n\tdims and alts\n============================\n");



//terminate the the fieldlist
    for (; fieldno < (t->var_count + 1+altvarcount); fieldno++)
    {
        field_list[fieldno].field_type = NULL;
        field_list[fieldno].field_name = NULL;
        field_list[fieldno].field_offset = 0;
        field_list[fieldno].field_size = 0;
    }

//     fprintf(stderr, "=======================\n");


//dump field list to check
    {

        FMField *f = field_list;
        int x = 0;

        for (x = 0; x < fieldno; x++)
        {
            f = &field_list[x];
            if (f == NULL || f->field_name == NULL || f->field_size == 0)
                break;

//              fprintf(stderr, "%d: %s %s %d %d\n",
//                      x, f->field_name, f->field_type, f->field_size,
//                      f->field_offset);

        }

    }

//associate field list
    format->field_list = field_list;

    current_fm->format->struct_size = current_fm->size;
    current_fm->buffer = (unsigned char *) malloc(current_fm->size);
    memset(current_fm->buffer, 0, current_fm->size);

#if HAVE_PORTALS == 1
//defined(__CRAYXT_COMPUTE_LINUX_TARGET)
    current_fm->s = InitIOFromFile("param", rank);
    current_fm->s->rank = rank;
    current_fm->ioformat = register_data(current_fm->s, current_fm->format);
#elif HAVE_INFINIBAND == 1
    current_fm->s = EVthin_ib_InitIOFile("param", 1, rank);
    current_fm->ioformat =
        EVthin_ib_registerData(current_fm->s, current_fm->format);
    {
        //we can read the code from here
         char codebuffer[1024*1024];
        char readbuffer[1024*1024];
        
        char *filename = getenv("FILTER");
        
        if(filename != NULL)
        {
            FILE *codefile  = fopen(filename, "r");
            if(codefile != NULL) 
            {
                fread(readbuffer, sizeof(char), 1024*1024, codefile);
                fclose(codefile);
            }
        
            if(!strcmp(filename, "warp_stat.c"))
            {    
//                sprintf(codebuffer, readbuffer, 0);                        
                set_code(current_fm->s, readbuffer);
            }
            else if(!strcmp(filename, "warp_bb.c"))
            {
//                sprintf(codebuffer, readbuffer);                        
                set_code(current_fm->s, readbuffer);                
            }
            else if(!strcmp(filename, "warp_order.c"))
            {
                set_code(current_fm->s, readbuffer);                                            
            }
            else if(!strcmp(filename, "warp_tag.c"))
            {
                set_code(current_fm->s, readbuffer);                                            
            }
            else if(!strcmp(filename, "warp_bbs.c"))
            {
                set_code(current_fm->s, readbuffer);                                            
            }

        }
    }
    
#endif

    current_fm->snd_count = 0;

    mdata->fm = current_fm;

    return 0;
    
}

static FMField *
internal_find_field(char *name, FMFieldList flist)
{
    FMField *f = flist;
    while (f->field_name != NULL && strcmp(f->field_name, name))
    {
        f++;
    }

    return f;
}
 
extern enum ADIOS_FLAG adios_datatap_should_buffer (struct adios_file_struct * fd
                                                   ,struct adios_method_struct * method)
{
}

extern void
adios_datatap_write(struct adios_file_struct *fd,
                    struct adios_var_struct *f,
                    const void *data, struct adios_method_struct *method)
{
    struct fm_structure *fm;
    dmd *mdata = (dmd *) method->method_data;


    struct adios_group_struct *group = method->group;

    fm = mdata->fm;

    if (group == NULL || fm == NULL)
    {
        fprintf(stderr, "group or fm is null - improperly initialized\n");
        return;

    }


    FMFieldList flist = fm->format->field_list;
    FMField *field = NULL;

    char *fixedname = findFixedName(fm, f->name);

    field = internal_find_field(fixedname, flist);
    if (field != NULL)
    {
        if (!f->dimensions)
        {
            //scalar quantity
            if (data)
            {
                //why wouldn't it have data?
                memcpy(&fm->buffer[field->field_offset], data,
                       field->field_size);

                //scalar quantities can have altnames also so assign those
                if(field->field_name != NULL)
                {
                    
                    dimnames *d = NULL;
                    for (d = fm->dimlist.lh_first; d != NULL;
                         d = d->entries.le_next)
                    {
                        if (!strcmp(d->name, field->field_name))
                        {
                            //matches
                            //check if there are altnames
                            altname *a = NULL;
                            for (a = d->altlist.lh_first; a != NULL;
                                 a = a->entries.le_next)
                            {
                                //use the altname field to get the data into the buffer
                                memcpy(&fm->buffer[a->field->field_offset], data,
                                       a->field->field_size);
                                //debug
//                                  int *testingint = (int*)&fm->buffer[a->field->field_offset];
                                
//                                   fprintf(stderr, "writing %s to %s at %d %d\n",
//                                           f->name, a->name, a->field->field_offset,
//                                          (int)*testingint);

                            }
                        }
                    }
                }
            }

            else
            {
                fprintf(stderr, "no data for  scalar %s\n", f->name);

            }


        }
        else
        {
            //vector quantity
            if (data)
            {
                //we just need to copy the pointer stored in f->data
                memcpy(&fm->buffer[field->field_offset], &data, sizeof(void *));

            }
            else
            {
                fprintf(stderr, "no data for vector %s\n", f->name);
            }
        }
    }
}

static void internal_adios_datatap_write(struct adios_file_struct *fd,
                                         struct adios_method_struct *method);

extern void
adios_datatap_close(struct adios_file_struct *fd,
                    struct adios_method_struct *method)
{

    dmd *mdata = method->method_data;


    if (!mdata->initialized)
    {
        return;
    }


    if (fd->mode & adios_mode_write)
    {
        internal_adios_datatap_write(fd, method);
    }

}


void
internal_adios_datatap_write(struct adios_file_struct *fd,
                             struct adios_method_struct *method)
{


    if (fd == NULL)
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

    if (t == NULL || fm == NULL)
    {
        fprintf(stderr, "improperly initialized for write\n");

        return;

    }

//     fprintf(stderr, "\t---\n");
//     fflush(stderr);


//          FMContext src_context = create_local_FMcontext(NULL);
//           FMFormat ioformat = register_data_format(src_context, fm->format);
//          int size =0;
//           FFSBuffer encode_buffer = create_FFSBuffer();       
//          char *xfer = FFSencode(encode_buffer, ioformat, fm->buffer, &size);
//           FFSFile file = open_FFSfile("/tmp/temp", "w");
//           write_FFSfile(file, ioformat, fm->buffer);
//           close_FFSfile(file);

//#if defined(__CRAYXT_COMPUTE_LINUX_TARGET)
#if HAVE_PORTALS == 1
    if (fm->snd_count > 0)
    {
        send_end(fm->s);
    }
#elif HAVE_INFINIBAND == 1
    if (fm->snd_count > 0)
    {
        EVthin_ib_endSend(fm->s);
    }
#endif


//now that we have all the info lets call write on each of these data types
//#if defined(__CRAYXT_COMPUTE_LINUX_TARGET)
#if HAVE_PORTALS == 1
    start_send(fm->s, fm->buffer, fm->size, fm->ioformat, NULL);

#elif HAVE_INFINIBAND == 1
//     fprintf(stderr, "size = %d\n", fm->size);
    
    EVthin_ib_startSend(fm->s, fm->buffer, fm->size, fm->ioformat, fm->alist);
#endif

    fm->snd_count++;


}

extern void
adios_datatap_finalize(int mype, struct adios_method_struct *method)
{

    dmd *mdata = method->method_data;


    struct adios_group_struct *t = method->group;



//initialize the globallist

//iterate through all the types
//find the correct format by name
    struct fm_structure *fm = mdata->fm;

    if (t == NULL || fm == NULL)
    {
        fprintf(stderr, "improperly initialized for finalize %p %p\n", t, fm);

        return;

    }

#if HAVE_PORTALS == 1
    if (fm->snd_count > 0)
    {
        char buffer[128];
        sprintf(buffer, "client-%d", rank);

        outputTimingInfo(buffer);

        send_end(fm->s);

    }
#elif HAVE_INFINIBAND == 1
    if (fm->snd_count > 0)
    {

        char buffer[128];
        sprintf(buffer, "client-%d", rank);

        EVthin_ib_outputTimingInfo(buffer);

        EVthin_ib_endSend(fm->s);
    }
#endif

}


extern void
adios_datatap_end_iteration(struct adios_method_struct *method)
{
    struct fm_structure *fm;
    dmd *mdata = (dmd *) method->method_data;
    fm = mdata->fm;

    if (fm == NULL)
        return;


#if HAVE_PORTALS == 1
    startIter(fm->s);
#elif HAVE_INFINIBAND == 1
    EVthin_ib_startIter(fm->s);
#endif

    mdata->cycle_id = 0;

}

extern void
adios_datatap_start_calculation(struct adios_method_struct *method)
{
    struct fm_structure *fm;
    dmd *mdata = (dmd *) method->method_data;
    fm = mdata->fm;

    if (fm == NULL)
        return;


#if HAVE_PORTALS == 1
    startCompute(fm->s, mdata->cycle_id);
#elif HAVE_INFINIBAND == 1
    EVthin_ib_startCompute(fm->s, mdata->cycle_id);
#endif


}

extern void
adios_datatap_stop_calculation(struct adios_method_struct *method)
{
    struct fm_structure *fm;
    dmd *mdata = (dmd *) method->method_data;
    fm = mdata->fm;

    if (fm == NULL)
        return;


#if HAVE_PORTALS == 1
    endCompute(fm->s, mdata->cycle_id);
#elif HAVE_INFINIBAND == 1
    EVthin_ib_endCompute(fm->s, mdata->cycle_id);
#endif

    mdata->cycle_id++;

}

extern void
adios_datatap_get_write_buffer(struct adios_file_struct *fd,
                   struct adios_var_struct *f,
                   uint64_t *size,
                   void **buffer,
                   struct adios_method_struct *method)
{
    fprintf(stderr, "adios_datatap_write_get_buffer: datatap disabled, "
            "no portals support\n");
}


void
adios_datatap_read(struct adios_file_struct *fd,
           struct adios_var_struct *f,
           void *buffer, uint64_t buffer_size,
           struct adios_method_struct *method)
{

}

#else


void
adios_datatap_read(struct adios_file_struct *fd,
                   struct adios_var_struct *f,
                   void *buffer, struct adios_method_struct *method)
{

}

extern void
adios_datatap_get_write_buffer(struct adios_file_struct *fd,
                               struct adios_var_struct *f,
                               unsigned long long *size,
                               void **buffer,
                               struct adios_method_struct *method)
{
}

extern void
adios_datatap_stop_calculation(struct adios_method_struct *method)
{
}

extern void
adios_datatap_start_calculation(struct adios_method_struct *method)
{
}

extern void
adios_datatap_end_iteration(struct adios_method_struct *method)
{
}

extern void
adios_datatap_finalize(int mype, struct adios_method_struct *method)
{
}

extern void
adios_datatap_close(struct adios_file_struct *fd,
                    struct adios_method_struct *method)
{
}

extern void
adios_datatap_write(struct adios_file_struct *fd,
                    struct adios_var_struct *f,
                    void *data, struct adios_method_struct *method)
{
}

extern void
adios_datatap_open(struct adios_file_struct *fd,
                   struct adios_method_struct *method)
{
}

extern void
adios_datatap_init(const PairStruct *params, struct adios_method_struct *method)
{
}

enum ADIOS_FLAG adios_datatap_should_buffer (struct adios_file_struct * fd
                                            ,struct adios_method_struct * method)
{
}
#endif

/*    FILE *formatfile = fopen("/tmp/formatfile","a+");
    if(formatfile == NULL) {
        fprintf(stderr, "can not open format file\n");
        exit(0);
    }
    fprintf(formatfile, "\n\n\n");
    fclose(formatfile);

*/
