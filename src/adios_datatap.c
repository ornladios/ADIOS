#include "config.h"

#if NO_DATATAP == 0

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ffs.h>
#include <atl.h>
#include <mpi.h>
#include "adios.h"
#include "adios_internals.h"
#include "adios_transport_hooks.h"
#include "bp_utils.h"
// added by zf2 (moved getFixedname and get_full_path_name to globals.c)
#include "globals.h"

#include <sys/queue.h>
#if HAVE_PORTALS == 1
#define QUEUE_H 1
#include <thin_portal.h>
#elif HAVE_INFINIBAND == 1
#include <thin_ib.h>
#endif

//#include "memwatch.h"



#define STARTINGSIZE 16

//static LIST_HEAD(listhead, fm_structure) globallist;

// added by zf2: moved to globals.c
//#define OPLEN 4
//static char OP[OPLEN] = { '+', '-', '*', '/' };
//static char *OP_REP[OPLEN] = { "_plus_", "_minus_", "_mult_", "_div_" };

// added by zf2: rank and mpisize will be obsolete. We will use a per-group mpi comm/rank/size
int rank = -1;
int mpisize  = 1;

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

typedef struct _var_to_write
{
    struct adios_var_struct *var;
    struct _var_to_write *next;
} var_to_write;

typedef struct _vars_list
{
    var_to_write *head;
    var_to_write *tail;
} vars_list;


struct fm_structure
{
    FMFormatRec *format;
    int size;                                                   //in bytes - no padding
    unsigned char *buffer;
    int snd_count;  // TODO: this count for adios_close
    IOhandle *s;
    FMFormat ioformat;
    attr_list alist;    

    // added by zf2
    vars_list vars_to_write;
    char *var_info_buffer;
    int var_info_buffer_size;
    int num_vars;
        
    LIST_HEAD(tablehead, nametable_) namelist;
    LIST_HEAD(dims, dimnames_) dimlist;
};


typedef struct datatap_method_data
{
    int opencount;  // TODO: this count for adios_open
    int initialized;
    int cycle_id;
    char *pfile;
    struct fm_structure *fm;
} dmd;

MPI_Comm adios_mpi_comm_world = MPI_COMM_WORLD;
int initialized = 0;


static int add_var_to_list(vars_list *var_list, struct adios_var_struct *var)
{
    var_to_write *t = (var_to_write *) malloc(sizeof(var_to_write));
    if(!t) {
        fprintf(stderr, "cannot allocate memory.\n");
        return -1;
    }
    t->var = var;
    t->next= NULL;

    if(var_list->tail) {
        var_list->tail->next = t;
    }
    var_list->tail = t;
    if(!var_list->head) {
        var_list->head = t;
    }
    return 0;
}

static void reset_var_list(vars_list *var_list)
{
    var_list->head = var_list->tail = NULL;
}

static struct adios_var_struct *next_var_from_list(vars_list *var_list)
{
    struct adios_var_struct *v;
    if(!var_list->head) {
        return NULL;
    }
    else {
        var_to_write *t;
        t = var_list->head;
        v = var_list->head->var;
        var_list->head = var_list->head->next;
        free(t);
        return v;
    }
}

// added by zf2: moved to globals.c

//static char *
//getFixedName(char *name)
//{
//    char *tempname = (char *) malloc(sizeof(char) * 255);
//    tempname = strdup(name);
//
//
//    char *oldname = strdup(name);
//    char *loc = NULL;
//    int i;
//
//    do
//    {
//        for (i = 0; i < OPLEN; i++)
//        {
//    //checking operator OP[i]
//            loc = strchr(oldname, OP[i]);
//            if (loc == NULL)
//                continue;
//            *loc = 0;
//            snprintf(tempname, 255, "%s%s%s", oldname, OP_REP[i], &loc[1]);
//            free(oldname);
//            oldname = strdup(tempname);
//        }
//    }
//    while (loc != NULL);
//
//
//
//    return tempname;
//}



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
fprintf(stderr, "im here orig %s mangle %s name %s %s:%d\n", node->originalname,node->mangledname,name,__FILE__,__LINE__);
            return node->mangledname;
        }

    }

fprintf(stderr, "im here orig %s mangle %s name %s %s:%d\n", node->originalname,node->mangledname,name,__FILE__,__LINE__);
    return name;
}


extern void
adios_datatap_init(const char *params, struct adios_method_struct *method)
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

fprintf(stderr, "im here param file %s %s:%d\n", params, __FILE__, __LINE__);
    mdata->opencount = 0;
    mdata->initialized = 1;
    if (params != NULL && strlen(params) >= 1)
    {
//contains the file name of the file to read?
        // added by zf2: here we use config.xml to pass reader application id
        mdata->pfile = strdup(params);
    }
    else {
        mdata->pfile = NULL;
    }
    
    MPI_Comm_rank(adios_mpi_comm_world, &rank);
    MPI_Comm_size(adios_mpi_comm_world, &mpisize);

    fflush(stderr);

}

static altname *
findAltName(struct fm_structure *current_fm, char *dimname, char *varname)
{
fprintf(stderr, "im here dim %s var %s %s:%d\n", dimname, varname, __FILE__, __LINE__);
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

        // added by zf2
        //d->name = dimname;
        d->name = strdup(dimname);
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


fprintf(stderr, "im here dim %s mangle %s %s:%d\n", dimname, aname, __FILE__, __LINE__);
    LIST_INSERT_HEAD(&d->altlist, a, entries);

    return a;
}



extern int
adios_datatap_open(struct adios_file_struct *fd,
           struct adios_method_struct *method, void*comm)
{
double t1 = 0, t2, t3;
double t11, t12, t13, t14;
double t_start = getlocaltime();

t11=getlocaltime();

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
            return -1;
        }
    }
    else
    {
        fprintf(stderr, "method not initialized\n");
        return -1;
    }

    if (mdata->fm != NULL)
    {
        // added by zf2:
        if(strlen(fd->name) >= MAX_FILE_NAME_LEN) {
            fprintf(stderr, "file name exceed %d\n", MAX_FILE_NAME_LEN);
            return -1;
        }

        memcpy(mdata->fm->s->fname, fd->name, strlen(fd->name));
        mdata->fm->s->fname[strlen(fd->name)] = '\0';
        mdata->fm->s->timestep = mdata->fm->snd_count; // TODO: starts from 0
        //mdata->fm->s->timestep = fd->group->time_index; // starts from 1
        mdata->opencount ++;
        return 0;
        //return -1;
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
    
    // added by zf2
    reset_var_list(&(current_fm->vars_to_write));
    current_fm->var_info_buffer_size = 8;
    current_fm->num_vars = 0;
    current_fm->snd_count = 0; // first time 


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

t12=getlocaltime();

//for each type look through all the fields
    struct adios_var_struct *f;
    for (f = t->vars; f != NULL; f = f->next, fieldno++)
    {
        //make the field list
        //check name for + - * / (operators) and replace them
        //char *tempname = getFixedName(f->name);

        // added by zf2: we use full path name to distinghish identical var names
        char *tempname;
        char *full_path_name = NULL;
        int name_changed = 0;
//        if(t->all_unique_var_names == adios_flag_no) {
t2=getlocaltime();
fprintf(stderr, "im here %s %s:%d\n", f->name,  __FILE__,__LINE__);
        full_path_name = get_full_path_name(f->name, f->path);
fprintf(stderr, "im here %s %s:%d\n", f->name,  __FILE__,__LINE__);
t1+=getlocaltime()-t2;
        //tempname = full_path_name;
        tempname = getFixedName(full_path_name);
  
//            // NOTE: "/" will be replaced with "_div_"
//            tempname = getFixedName(full_path_name);
//            name_changed = 1;
//            //free(full_path_name);
//        }
//        else {
//            tempname = getFixedName(f->name);
//            name_changed = strcmp(tempname, f->name);
//        }

//        if (name_changed)
//        {
            //strings don't match
            //add to name list
        nametable *namenode = (nametable *) malloc(sizeof(nametable));
//            if(full_path_name) {
        namenode->originalname = full_path_name;
//            } 
//            else {
//                namenode->originalname = strdup(f->name);  
//            }
        namenode->mangledname = strdup(tempname);

fprintf(stderr, "im here orig %s mangle %s %s:%d\n", namenode->originalname, namenode->mangledname, __FILE__, __LINE__);
        LIST_INSERT_HEAD(&current_fm->namelist, namenode, entries);
//        }


//
        field_list[fieldno].field_name = tempname;

        // added by zf2 
        int num_dims = 0;
        struct adios_dimension_struct *dim = f->dimensions;
        for(; dim != NULL; dim = dim->next) {
            if(dim->dimension.time_index != adios_flag_yes) {
                num_dims ++;
            }
        }
        
        if(!num_dims)
        {              
            while(current_fm->size % 8 != 0)
            {
                current_fm->size ++;
            }



            // scalar 
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
                    field_list[fieldno].field_size = sizeof(char *); 
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

            // added by zf2: here we need to handle a special case where 
            // all dimensions are number, such as a[32][32][32], in this
            // case the offset should be calculated to include the whole
            // array "in place" otherwise FFSencode will fail
            int dim_are_numbers = 1;

            //create the dimension thingy
            for (; d != NULL; d = d->next)
            {
                // added by zf2: time index 
                if(d->dimension.time_index == adios_flag_yes) {
                    snprintf(el, ELSIZE, "[1]");
                    continue;
                }
  
                //for each dimension just take the upper_bound
                if (d->dimension.id)
                {
                    //findFixedName returns the mangled name from the original name
                    struct adios_var_struct *tmp_var = adios_find_var_by_id(t->vars, d->dimension.id);

                    // added by zf2: we use full path name here
                    char *name;
 //                   if(t->all_unique_var_names == adios_flag_no) {
t2=getlocaltime();
fprintf(stderr, "im here %s %s %s %s %s:%d\n", f->name, field_list[fieldno].field_name, tmp_var->name, tmp_var->path, __FILE__,__LINE__);
                    char *full_pathname = get_full_path_name(tmp_var->name, tmp_var->path);
t1+=getlocaltime()-t2;
fprintf(stderr, "im here %s %s fixedname %s %s:%d\n", f->name, field_list[fieldno].field_name, full_pathname, __FILE__,__LINE__);
                    name = findFixedName(current_fm, full_pathname);
fprintf(stderr, "im here %s %s fixedname %s %s %s:%d\n", f->name, field_list[fieldno].field_name, full_pathname, name, __FILE__,__LINE__);
                    free(full_pathname);
                    //name = full_pathname;

//                    }
//                    else {
//                        name = findFixedName(current_fm, tmp_var->name);
//                    } 

                    //create the alternate name for this variable and the array its defining
                    altname *a = findAltName(current_fm, name,
                                 (char*)field_list[fieldno].field_name);
                    //altname is a new variable that we need to add to the field list
                    altvarcount++;

                    snprintf(el, ELSIZE, "[%s]", a->name);
fprintf(stderr, "im here %s %s:%d\n",a->name, __FILE__,__LINE__);
//                    fprintf(stderr, "%s\t", el);
                    //offset_increment is just the size of the pointer

                    // added by zf2:
                    dim_are_numbers = 0;
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

            // added by zf2: we need to put the whole array in place
            uint64_t array_size = 1;
            if(dim_are_numbers) {
                for (d = f->dimensions; d != NULL; d = d->next) {
fprintf(stderr, "im here %s %lu %s:%d\n", f->name, d->dimension.rank, __FILE__,__LINE__);
                    if(d->dimension.time_index == adios_flag_yes) {
                        continue;
                    }
                    array_size *= d->dimension.rank;
                }
fprintf(stderr, "im here %s %lu %s:%d\n", f->name, array_size, __FILE__,__LINE__);
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
                    
                    // added by zf2
                    if(!dim_are_numbers) {
                        current_fm->size += sizeof(void *);
                    }
                    else {
                        current_fm->size += array_size * field_list[fieldno].field_size; 
                    }
                    break;

                case adios_real:
                    field_list[fieldno].field_type =
                        (char *) malloc(sizeof(char) * 255);
                    snprintf((char *) field_list[fieldno].field_type, 255,
                         "float%s", dims);
                    field_list[fieldno].field_size = sizeof(float);
                    field_list[fieldno].field_offset = current_fm->size;

                    if(!dim_are_numbers) {
                        current_fm->size += sizeof(void *);
                    }
                    else {
                        current_fm->size += array_size * field_list[fieldno].field_size;
                    }
                    break;

                case adios_string:
                    field_list[fieldno].field_type = strdup("string");
                    field_list[fieldno].field_size = sizeof(char *);
                    field_list[fieldno].field_offset = current_fm->size;

                    if(!dim_are_numbers) {
                        current_fm->size += sizeof(void *);
                    }
                    else {
                        current_fm->size += array_size * field_list[fieldno].field_size;
                    }
                    break;

                case adios_double:
                    field_list[fieldno].field_type =
                        (char *) malloc(sizeof(char) * 255);
                    snprintf((char *) field_list[fieldno].field_type, 255,
                         "float%s", dims);
                    field_list[fieldno].field_size = sizeof(double);
                    field_list[fieldno].field_offset = current_fm->size;

                    if(!dim_are_numbers) {
                        current_fm->size += sizeof(void *);
                    }
                    else {
                        current_fm->size += array_size * field_list[fieldno].field_size;
                    }
                    break;

                case adios_byte:
                    field_list[fieldno].field_type =
                        (char *) malloc(sizeof(char) * 255);
                    snprintf((char *) field_list[fieldno].field_type, 255, "char%s",
                         dims);
                    field_list[fieldno].field_size = sizeof(char);
                    field_list[fieldno].field_offset = current_fm->size;

                    if(!dim_are_numbers) {
                        current_fm->size += sizeof(void *);
                    }
                    else {
                        current_fm->size += array_size * field_list[fieldno].field_size;
                    }
                    break;

                default:
                    fprintf(stderr, "unknown type error %d\n", f->type);
                    fieldno--;
                    break;
            }

        }

      fprintf(stderr, "rank=%d %s, %s, %d, %d\n", rank,field_list[fieldno].field_name, field_list[fieldno].field_type,field_list[fieldno].field_size,field_list[fieldno].field_offset); 


    }
t13=getlocaltime();

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

t14=getlocaltime();


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

              fprintf(stderr, "%d: %s %s %d %d\n",
                      x, f->field_name, f->field_type, f->field_size,
                      f->field_offset);

        }

    }

double t_end = getlocaltime();
fprintf(stderr, "im here rank=%d %f %f t1 %f %s:%d\n", rank, t_start,t_end,t1,__FILE__,__LINE__);
fprintf(stderr, "im here rank=%d %f %f %f %f %s:%d\n", rank, t11,t12,t13,t14,__FILE__,__LINE__);



//associate field list
    format->field_list = field_list;

    current_fm->format->struct_size = current_fm->size;
    current_fm->buffer = (unsigned char *) malloc(current_fm->size);
    memset(current_fm->buffer, 0, current_fm->size);

#if HAVE_PORTALS == 1
//defined(__CRAYXT_COMPUTE_LINUX_TARGET)
    // added by zf2: for now just assume we have only one reader application
    char param_file[50];
    if(mdata->pfile) { // TODO: this should be per-group based not per method
        sprintf(param_file, "datatap_param%s\0", mdata->pfile);                    
    }
    else {
        strcpy(param_file, "datap_param");         
    }
fprintf(stderr, "im here param file %s %s:%d\n", param_file, __FILE__, __LINE__);
    current_fm->s = InitIOFromFile(param_file, rank,mpisize);
    current_fm->s->rank = rank;
    current_fm->ioformat = register_data(current_fm->s, current_fm->format);
    
    // added by zf2
    if(strlen(fd->name) >= MAX_FILE_NAME_LEN) {
        fprintf(stderr, "file name exceed %d\n", MAX_FILE_NAME_LEN);                    
        return -1;
    }
    memcpy(current_fm->s->fname, fd->name, strlen(fd->name)); 
    current_fm->s->fname[strlen(fd->name)] = '\0';
    current_fm->s->timestep = current_fm->snd_count; // TODO: starts from 0
    //current_fm->s->timestep = fd->group->time_index; // starts from 1

    
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

    // added by zf2
    mdata->opencount ++;

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
                    void *data, struct adios_method_struct *method)
{
fprintf(stderr, "im here %s:%d\n",__FILE__,__LINE__);
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

    // added by zf2
    char *fixedname;
//    if(group->all_unique_var_names == adios_flag_no) {
    char *full_path_name = get_full_path_name(f->name, f->path);
//    fixedname = full_path_name;

    fixedname = findFixedName(fm, full_path_name);
    free(full_path_name);
//    }
//    else {
//        fixedname = findFixedName(fm, f->name);
//    }

    //char *fixedname = findFixedName(fm, f->name);

    field = internal_find_field(fixedname, flist);

//    free(fixedname);

    if (field != NULL)
    {

        // added by zf2
        int num_dims = 0;
        struct adios_dimension_struct *dim = f->dimensions;
        for(; dim != NULL; dim = dim->next) {
            if(dim->dimension.time_index != adios_flag_yes) {
                num_dims ++;
            }
        }

        // string should be treat differently
        if(!strcmp(field->field_type, "string")) {
            if(data) {
                int array_size = bp_get_type_size(f->type, data);
fprintf(stderr, "im here var %s var info size %d %s:%d\n",f->name,array_size,__FILE__,__LINE__);

#if 0
                //we just need to copy the pointer stored in f->data

                memcpy(&fm->buffer[field->field_offset], &data, sizeof(void *));
#else

                // allocate a temp buffer
                void *temp_buffer = malloc(array_size);
                if(!temp_buffer) {
                    fprintf(stderr, "cannot allocate memory. %s:%d\n", __FILE__, __LINE__);
                    return;
                }

                // copy data
                memcpy(temp_buffer, data, array_size);

                // record the buffer address
                memcpy(&fm->buffer[field->field_offset], &temp_buffer, sizeof(void *));
#endif

                // added by zf2
fprintf(stderr, "im here var %s var info size %d %s:%d\n",f->name,fm->var_info_buffer_size,__FILE__,__LINE__);
                add_var_to_list(&(fm->vars_to_write), f);
                fm->num_vars ++;
                fm->var_info_buffer_size += 4; // size of var info
                fm->var_info_buffer_size += strlen(f->name) +
                    strlen(f->path) + 2 + 8;  // var name and path
                fm->var_info_buffer_size += 4; // id
                fm->var_info_buffer_size += 4; // type
                fm->var_info_buffer_size += 4; // whether have time dim
                fm->var_info_buffer_size += 4; // dim = 0
                fm->var_info_buffer_size += array_size; // data value
                //fm->var_info_buffer_size += 24; // TODO: let's just support plain string now
fprintf(stderr, "im here var %s var info size %d %s:%d\n",f->name,fm->var_info_buffer_size,__FILE__,__LINE__);
            }
            else {
                fprintf(stderr, "no data for string %s\n", f->name);
            }
        }
        else if (!num_dims)
        {
            //scalar quantity
            if (data)
            {
                //why wouldn't it have data?

                memcpy(&fm->buffer[field->field_offset], data,
                       field->field_size);
                // added by zf2
fprintf(stderr, "im here var %s var info size %d %s:%d\n",f->name,fm->var_info_buffer_size,__FILE__,__LINE__);
                add_var_to_list(&(fm->vars_to_write), f); 
                fm->num_vars ++;
                fm->var_info_buffer_size += 4; // size of var info
                fm->var_info_buffer_size += strlen(f->name) +
                    strlen(f->path) + 2 + 8;  // var name and path
                fm->var_info_buffer_size += 4; // id
                fm->var_info_buffer_size += 4; // type
                fm->var_info_buffer_size += 4; // whether hav time dim
                fm->var_info_buffer_size += 4; // dim = 0
                fm->var_info_buffer_size += bp_get_type_size(f->type, data); // data value
fprintf(stderr, "im here var %s var info size %d %s:%d\n",f->name,fm->var_info_buffer_size,__FILE__,__LINE__);

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
                                //int *testingint = (int*)&fm->buffer[a->field->field_offset];
                                
                                //   fprintf(stderr, "writing %s to %s at %d %d\n",
                                //           f->name, a->name, a->field->field_offset,
                                //          (int)*testingint);

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

// added by zf2: we make an extra copy to deal with temporary Fortran array

#if 0          
                //we just need to copy the pointer stored in f->data

                memcpy(&fm->buffer[field->field_offset], &data, sizeof(void *));
#else 
                // first figure out the dimension values 
                uint64_t array_size = bp_get_type_size(f->type, data);
                struct adios_dimension_struct *d = f->dimensions;
                
                // also we need to determine if this array is static or dynamic 
                int dim_are_numbers = 1;
fprintf(stderr, "im here var %s size %d %s:%d\n", f->name, array_size, __FILE__,__LINE__);

                for(; d != NULL; d = d->next) {
                    if(d->dimension.time_index == adios_flag_yes) {
                        continue;
                    }

                    //for each dimension just take the upper_bound
                    if (d->dimension.id) {
                        // it's dynamic array
                        dim_are_numbers = 0;                        

                        struct adios_var_struct *dim_var = adios_find_var_by_id(group->vars, d->dimension.id);
                        char *full_pathname = get_full_path_name(dim_var->name, dim_var->path);
                        char *dimname = getFixedName(full_pathname); 
                        free(full_pathname);
                        FMField *dim_field = internal_find_field(dimname, flist);
                         
                        // get the value
                        switch(dim_var->type) {
                            case adios_short:
                                {
                                short int temp;
                                memcpy(&temp, &fm->buffer[dim_field->field_offset], dim_field->field_size);
                                array_size *= temp;
                                break;
                                }
                            case adios_unsigned_short:
                                {
                                unsigned short int temp;
                                memcpy(&temp, &fm->buffer[dim_field->field_offset], dim_field->field_size);
                                array_size *= temp;
                                break;
                                }
                            case adios_integer:
                                {
                                int temp;
                                memcpy(&temp, &fm->buffer[dim_field->field_offset], dim_field->field_size);
                                array_size *= temp;
                                break;
                                }
                            case adios_unsigned_integer:
                                {
                                unsigned int temp;
                                memcpy(&temp, &fm->buffer[dim_field->field_offset], dim_field->field_size);
                                array_size *= temp;
                                break;
                                }
                            case adios_long:
                                {
                                long int temp;
                                memcpy(&temp, &fm->buffer[dim_field->field_offset], dim_field->field_size);
                                array_size *= temp;
                                break;
                                }
                            case adios_unsigned_long:
                                {
                                unsigned long int temp;
                                memcpy(&temp, &fm->buffer[dim_field->field_offset], dim_field->field_size);
                                array_size *= temp;
                                break;
                                }
                            default:
                                fprintf(stderr, "cannot support data type %d as dimension.", dim_var->type);
                        }
 
                        free(dimname);
                    }
                    else {   //its a number
                        array_size *= d->dimension.rank;
                    }
                }

fprintf(stderr, "im here var %s size %d %s:%d\n", f->name, array_size, __FILE__,__LINE__);
                if(dim_are_numbers) { // copy statis array in place
fprintf(stderr, "im here var %s size %d %s:%d\n", f->name, array_size, __FILE__,__LINE__);
                    memcpy(&fm->buffer[field->field_offset], data, array_size);
                } else {
                    // allocate a temp buffer   
fprintf(stderr, "im here var %s size %d %s:%d\n", f->name, array_size, __FILE__,__LINE__);
                    void *temp_buffer = malloc(array_size);      
                    if(!temp_buffer) {
                        fprintf(stderr, "cannot allocate memory. %s:%d\n", __FILE__, __LINE__);  
                        return; 
                    }

                    // copy data
                    memcpy(temp_buffer, data, array_size);

fprintf(stderr, "im here var %s offset %d addr %p %p %s:%d\n", f->name, field->field_offset, &fm->buffer[field->field_offset], temp_buffer, __FILE__,__LINE__);
                    // record the buffer address
                    memcpy(&fm->buffer[field->field_offset], &temp_buffer, sizeof(void *));
                } 
#endif
                // added by zf2
fprintf(stderr, "im here var %s var info size %d %s:%d\n",f->name,fm->var_info_buffer_size,__FILE__,__LINE__);
                add_var_to_list(&(fm->vars_to_write), f); 
                fm->num_vars ++;
                fm->var_info_buffer_size += 4; // size of var info
                fm->var_info_buffer_size += strlen(f->name) + 
                    strlen(f->path) + 2 + 8;  // var name and path
                fm->var_info_buffer_size += 4; // id
                fm->var_info_buffer_size += 4; // type
                fm->var_info_buffer_size += 4; // whether hav time dim
                fm->var_info_buffer_size += 4; // dim = N
                for(d = f->dimensions; d != NULL; d = d->next) {
                    fm->var_info_buffer_size += 24;
                }                                                    
fprintf(stderr, "im here var %s var info size %d %s:%d\n",f->name,fm->var_info_buffer_size,__FILE__,__LINE__);
            }
            else
            {
                fprintf(stderr, "no data for vector %s\n", f->name);
            }
        }
    }
fprintf(stderr, "im here %s:%d\n",__FILE__,__LINE__);
}

static void internal_adios_datatap_write(struct adios_file_struct *fd,
                                         struct adios_method_struct *method);

extern void
adios_datatap_close(struct adios_file_struct *fd,
                    struct adios_method_struct *method)
{
fprintf(stderr, "im here %s:%d\n",__FILE__,__LINE__);
    dmd *mdata = method->method_data;


    if (!mdata->initialized)
    {
        return;
    }

    // added by zf2: we need to support both write (a new file name)
    // and append (same file name, new timestep)
    if ((fd->mode & adios_mode_write) || (fd->mode & adios_mode_append))
    {
fprintf(stderr, "im here %s:%d\n",__FILE__,__LINE__);
        internal_adios_datatap_write(fd, method);
fprintf(stderr, "im here %s:%d\n",__FILE__,__LINE__);
    }

fprintf(stderr, "im here %s:%d\n",__FILE__,__LINE__);
}

char * create_var_info_buffer(struct adios_group_struct *g, struct fm_structure *fm)
{
fprintf(stderr, "im here %s:%d\n",__FILE__,__LINE__);
double t_start = getlocaltime();


    // format of binpacked meta-data
    // 4: total size
    // 4: language
    // 4: size of group name
    // K: group name
    // 4: num of vars
    // 4: size of var info
    // 2: var id
    // 4: size of var name
    // N: var name 
    // 4: size of var path
    // M: var path
    // 4: var type
    // 4: -1: no time dimension; otherwise specify which dimension is time dimension
    // 4: number of dimensions
    // 24: local bounds/global bounds/global offset for 1st dim
    // 24: local bounds/global bounds/global offset for 2nd dim
    // ...
    // 24: local bounds/global bounds/global offset for last dim
fprintf(stderr, "im here var info size %d %s:%d\n",fm->var_info_buffer_size,__FILE__,__LINE__);
    
    //struct adios_var_struct *vars = g->vars;    
    fm->var_info_buffer_size += sizeof(enum ADIOS_FLAG); // host language
    fm->var_info_buffer_size += sizeof(int); // size of group name
    fm->var_info_buffer_size += strlen(g->name) + 1;

fprintf(stderr, "im here var info size %d %s:%d\n",fm->var_info_buffer_size,__FILE__,__LINE__);
    // TODO: alignment
    char *buf = (char *) malloc(fm->var_info_buffer_size);
    if(!buf) {
        fprintf(stderr, "Cannot allocate memory. %s:%d\n", __FILE__, __LINE__);      
        exit(-1); 
    }
    
    char *current_pos = buf;    
    memcpy(current_pos, &(fm->var_info_buffer_size), sizeof(int)); // total size
    current_pos += sizeof(int);

    memcpy(current_pos, &(g->adios_host_language_fortran), sizeof(enum ADIOS_FLAG)); // host language
    current_pos += sizeof(enum ADIOS_FLAG);
fprintf(stderr, "");

    int group_name_len = strlen(g->name) + 1; 
    memcpy(current_pos, &group_name_len, sizeof(int)); // size of group name   
    current_pos += sizeof(int);
    memcpy(current_pos, g->name, group_name_len);
    current_pos += group_name_len;

    memcpy(current_pos, &(fm->num_vars), sizeof(int)); // total num of vars
    current_pos += sizeof(int);

fprintf(stderr, "im here rank=%d %d %s:%d\n", rank, current_pos-buf,__FILE__,__LINE__);

    FMFieldList flist = fm->format->field_list;
    struct adios_var_struct *f;
    FMField *field = NULL;
    int var_info_size;
    int var_name_len, var_path_len;
    int var_type;
    //for (f = vars; f != NULL; f = f->next) {
    while(f=next_var_from_list(&(fm->vars_to_write))) {
        char *fixedname;
//        if(g->all_unique_var_names == adios_flag_no) {
        char *full_path_name = get_full_path_name(f->name, f->path);
        fixedname = findFixedName(fm, full_path_name);
        free(full_path_name);
fprintf(stderr, "im here varname %s %s:%d\n", f->name, __FILE__,__LINE__);
//            free(full_path_name);  
//        }
//        else {
//            fixedname = findFixedName(fm, f->name);
//        }
//        fixedname = full_path_name;

        field = internal_find_field(fixedname, flist);

//        free(fixedname);
fprintf(stderr, "im here varname %s %s:%d\n", f->name, __FILE__,__LINE__);

        var_name_len = strlen(f->name) + 1;
        var_path_len = strlen(f->path) + 1;
fprintf(stderr, "im here varname %s %s:%d\n", f->name, __FILE__,__LINE__);

        int ndims = 0;
        struct adios_dimension_struct *dim = f->dimensions;
        for(; dim != NULL; dim = dim->next) {
            if(dim->dimension.time_index != adios_flag_yes) {
                ndims ++;
            }
        }

        if(!ndims) { // scalars and strings           
            var_info_size = sizeof(int) + sizeof(int) + sizeof(int) * 2 + var_name_len + var_path_len
                + sizeof(enum ADIOS_DATATYPES) + sizeof(int) *2 + field->field_size;    
            memcpy(current_pos, &var_info_size, sizeof(int)); // size of var info
            current_pos += sizeof(int);
            int id = (int) f->id;
            memcpy(current_pos, &(id), sizeof(int)); // var id
            current_pos += sizeof(int);
            memcpy(current_pos, &var_name_len, sizeof(int)); // length of var name
            current_pos += sizeof(int);
fprintf(stderr, "im here id %d varname %d %s %s:%d\n",id, var_name_len,f->name, __FILE__,__LINE__);
            memcpy(current_pos, f->name, var_name_len); // var name
            current_pos += var_name_len;
            memcpy(current_pos, &var_path_len, sizeof(int)); // length of var path
            current_pos += sizeof(int);
            memcpy(current_pos, f->path, var_path_len); // var path
            current_pos += var_path_len;
            memcpy(current_pos, &(f->type), sizeof(enum ADIOS_DATATYPES)); // type
            current_pos += sizeof(enum ADIOS_DATATYPES);

            if(f->dimensions && f->dimensions->dimension.time_index == adios_flag_yes) {
                 *((int *)current_pos) = 0; // have time dimension
            }
            else {
                 *((int *)current_pos) = -1; // do not have time dimension
            }
            current_pos += sizeof(int);

            *((int *)current_pos) = 0; // no of dimensions
            current_pos += sizeof(int);

            if(f->type == adios_string) {
                char *str_data = *(char **)(&fm->buffer[field->field_offset]);
                memcpy(current_pos, str_data, strlen(str_data)+1); // data value
                current_pos += strlen(str_data)+1;
            } 
            else {
                memcpy(current_pos, &fm->buffer[field->field_offset], field->field_size); // data value
                current_pos += field->field_size; 
            }
fprintf(stderr, "im here rank=%d %d %s:%d\n", rank, current_pos-buf,__FILE__,__LINE__);
        }
        else { // arrays
            var_info_size = sizeof(int) + sizeof(int) + sizeof(int) * 2 + var_name_len + var_path_len
                + sizeof(enum ADIOS_DATATYPES) + sizeof(int) * 2; // we don't know ndims yet
            char *start_pos = current_pos;    
            //memcpy(current_pos, &var_info_size, sizeof(int)); // size of var info
            current_pos += sizeof(int);
            int id = (int) f->id;
            memcpy(current_pos, &(id), sizeof(int)); // var id
            current_pos += sizeof(int);
            memcpy(current_pos, &var_name_len, sizeof(int)); // length of var name
            current_pos += sizeof(int);
            memcpy(current_pos, f->name, var_name_len); // var name
            current_pos += var_name_len;
            memcpy(current_pos, &var_path_len, sizeof(int)); // length of var path
            current_pos += sizeof(int);
            memcpy(current_pos, f->path, var_path_len); // var path
            current_pos += var_path_len;
            memcpy(current_pos, &(f->type), sizeof(enum ADIOS_DATATYPES)); // type
            current_pos += sizeof(enum ADIOS_DATATYPES);   
fprintf(stderr, "im here id %d varname %d %s %s:%d\n",id, var_name_len,f->name, __FILE__,__LINE__);
        
            int *time_index = (int *) current_pos; // hold this to fill in later
            *time_index = -1; // -1 means no time dimension 
            current_pos += sizeof(int);

            int *no_dim_pos = (int *) current_pos; // hold this to fill in later
            current_pos += sizeof(int);

            uint64_t local_bound, global_bound, global_offset;
            int num_dims = 0;
            struct adios_dimension_struct *d = f->dimensions;
            for(; d != NULL; d = d->next) {
                num_dims ++;
                
                int is_timeindex = 0; 
                // local bounds
                if (d->dimension.time_index == adios_flag_yes) { // time index dimension
                    local_bound = 1; 
                    global_bound = 1; // TODO: we don't know how many timesteps in total so set to 1
                    global_offset = g->time_index - 1;
                    is_timeindex = 1;

                    // TODO: for fortran index starts from 1; for C index starts from 0 
                    *time_index = num_dims;
                }
                else if (d->dimension.id) {
                
                    // findFixedName returns the mangled name from the original name
                    struct adios_var_struct *dim_var = adios_find_var_by_id(g->vars, d->dimension.id);
                    char *dim_name;
//                    if(g->all_unique_var_names == adios_flag_no) {
                    char *full_path_name = get_full_path_name(dim_var->name, dim_var->path);
                    dim_name = findFixedName(fm, full_path_name);
//                        dim_name = full_path_name;

                    free(full_path_name);
//                    }
//                    else {
//                        dim_name = findFixedName(fm, dim_var->name);
//                    }
                     

                    int dim_size = bp_get_type_size(dim_var->type, NULL);
                    
                    // get the value of it                              
                    FMField *dim_field = internal_find_field(dim_name, flist);

//                    free(dim_name);

                    // TODO: cast data type
                    switch(dim_var->type) {
                        case adios_short:
                            {
                            short int temp;
                            memcpy(&temp, &fm->buffer[dim_field->field_offset], dim_size);
                            local_bound = (uint64_t) temp; 
                            break;
                            }
                        case adios_unsigned_short:
                            {
                            unsigned short int temp;
                            memcpy(&temp, &fm->buffer[dim_field->field_offset], dim_size);
                            local_bound = (uint64_t) temp; 
                            break;
                            }
                        case adios_integer:
                            {
                            int temp;
                            memcpy(&temp, &fm->buffer[dim_field->field_offset], dim_size);
                            local_bound = (uint64_t) temp; 
                            break;
                            }  
                        case adios_unsigned_integer:
                            {
                            unsigned int temp;
                            memcpy(&temp, &fm->buffer[dim_field->field_offset], dim_size);
                            local_bound = (uint64_t) temp; 
                            break;
                            }
                        case adios_long:
                            {
                            long int temp;
                            memcpy(&temp, &fm->buffer[dim_field->field_offset], dim_size);
                            local_bound = (uint64_t) temp; 
                            break;
                            } 
                        case adios_unsigned_long:
                            {
                            unsigned long int temp;
                            memcpy(&temp, &fm->buffer[dim_field->field_offset], dim_size);
                            local_bound = (uint64_t) temp; 
                            break;
                            }
                         default:
                            fprintf(stderr, "cannot support data type %d as dimension.", dim_var->type);    
                    }

                }
                else { //its a number
                    local_bound = d->dimension.rank;
                }                                  

                // global bounds
                if (is_timeindex) { // time index dimension
                    global_bound = 1; // TODO: we don't know how many timesteps in total so set to 1
                }
                else if (d->global_dimension.id) {
                
                    // findFixedName returns the mangled name from the original name
                    struct adios_var_struct *dim_var = adios_find_var_by_id(g->vars, d->global_dimension.id);
                    char *dim_name;
//                    if(g->all_unique_var_names == adios_flag_no) {
                    char *full_path_name = get_full_path_name(dim_var->name, dim_var->path);
//                    dim_name = full_path_name;

                    dim_name = findFixedName(fm, full_path_name);
//                        fprintf(stderr, "im here full %s fixed %s %s:%d\n",full_path_name,dim_name,__FILE__,__LINE__);
                    free(full_path_name);
//                    }
//                    else {
//                        dim_name = findFixedName(fm, dim_var->name);
//                    }

                    int dim_size = bp_get_type_size(dim_var->type, NULL);
                    
                    // get the value of it                              
                    FMField *dim_field = internal_find_field(dim_name, flist);

//                    free(dim_name);

                    // TODO: cast data type
                    switch(dim_var->type) {
                        case adios_short:
                            {
                            short int temp;
                            memcpy(&temp, &fm->buffer[dim_field->field_offset], dim_size);
                            global_bound = (uint64_t) temp; 
                            break;
                            } 
                        case adios_unsigned_short:
                            {
                            unsigned short temp;
                            memcpy(&temp, &fm->buffer[dim_field->field_offset], dim_size);
                            global_bound = (uint64_t) temp; 
                            break;
                            } 
                        case adios_integer:
                            {
                            int temp;
                            memcpy(&temp, &fm->buffer[dim_field->field_offset], dim_size);
                            global_bound = (uint64_t) temp; 
                            break;
                            } 
                        case adios_unsigned_integer:
                            {
                            unsigned int temp;
                            memcpy(&temp, &fm->buffer[dim_field->field_offset], dim_size);
                            global_bound = (uint64_t) temp; 
                            break;
                            } 
                        case adios_long:
                            {
                            long int temp;
                            memcpy(&temp, &fm->buffer[dim_field->field_offset], dim_size);
                            global_bound = (uint64_t) temp; 
                            break;
                            } 
                        case adios_unsigned_long:
                            {
                            unsigned long int temp;
                            memcpy(&temp, &fm->buffer[dim_field->field_offset], dim_size);
                            global_bound = (uint64_t) temp; 
                            break;
                            } 
                         default:
                            fprintf(stderr, "cannot support data type %d as dimension.", dim_var->type);    
                    }

                }
                else { //its a number
                    global_bound = d->global_dimension.rank;
                }                                  

                // global offset
                if (is_timeindex) { // time index dimension
                    global_offset = g->time_index - 1; 
                }
                else if (d->local_offset.id) {

                    // findFixedName returns the mangled name from the original name
                    struct adios_var_struct *dim_var = adios_find_var_by_id(g->vars, d->local_offset.id);
                    char *dim_name;
//                    if(g->all_unique_var_names == adios_flag_no) {
                    char *full_path_name = get_full_path_name(dim_var->name, dim_var->path);
//                    dim_name = full_path_name;

                    dim_name = findFixedName(fm, full_path_name);
//                        fprintf(stderr, "im here full %s fixed %s %s:%d\n",full_path_name,dim_name,__FILE__,__LINE__);
                    free(full_path_name);
//                    }
//                    else {
//                        dim_name = findFixedName(fm, dim_var->name);
//                    }

                    int dim_size = bp_get_type_size(dim_var->type, NULL);
                    
                    // get the value of it                              
                    FMField *dim_field = internal_find_field(dim_name, flist);

//                    free(dim_name); 

                    // TODO: cast data type
                    switch(dim_var->type) {
                        case adios_short:
                            {
                            short int temp;
                            memcpy(&temp, &fm->buffer[dim_field->field_offset], dim_size);
                            global_offset = (uint64_t) temp;
                            break;
                            }
                        case adios_unsigned_short:
                            {
                            unsigned short temp;
                            memcpy(&temp, &fm->buffer[dim_field->field_offset], dim_size);
                            global_offset = (uint64_t) temp;
                            break;
                            }
                        case adios_integer:
                            {
                            int temp;
                            memcpy(&temp, &fm->buffer[dim_field->field_offset], dim_size);
                            global_offset = (uint64_t) temp;
                            break;
                            }
                        case adios_unsigned_integer:
                            {
                            unsigned int temp;
                            memcpy(&temp, &fm->buffer[dim_field->field_offset], dim_size);
                            global_offset = (uint64_t) temp;
                            break;
                            }
                        case adios_long:
                            {
                            long int temp;
                            memcpy(&temp, &fm->buffer[dim_field->field_offset], dim_size);
                            global_offset = (uint64_t) temp;
                            break;
                            }
                        case adios_unsigned_long:
                            {
                            unsigned long int temp;
                            memcpy(&temp, &fm->buffer[dim_field->field_offset], dim_size);
                            global_offset = (uint64_t) temp;
                            break;
                            }
                         default:
                            fprintf(stderr, "cannot support data type %d as dimension.", dim_var->type);
                    }

                }
                else { //its a number
                    global_offset = d->local_offset.rank;
                }                                  
                
                // copy dimension info
                memcpy(current_pos, &local_bound, sizeof(uint64_t));
                current_pos += sizeof(uint64_t); 
                memcpy(current_pos, &global_bound, sizeof(uint64_t));
                current_pos += sizeof(uint64_t); 
                memcpy(current_pos, &global_offset, sizeof(uint64_t));
                current_pos += sizeof(uint64_t);               
            }                                                    
            
            *no_dim_pos = num_dims; // no of dimensions
            var_info_size += num_dims * sizeof(uint64_t) * 3;
            memcpy(start_pos, &var_info_size, sizeof(int)); // size of var info
fprintf(stderr, "im here rank=%d %d %s:%d\n", rank, current_pos-buf,__FILE__,__LINE__);
        } 
    }

double t_end = getlocaltime();
fprintf(stderr, "im here rank=%d %f %f %s:%d\n", rank, t_start,t_end,__FILE__,__LINE__);
fprintf(stderr, "im here rank=%d %d %s:%d\n", rank, current_pos-buf,__FILE__,__LINE__);

    return buf;
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
// changed by zf2
//    start_send(fm->s, fm->buffer, fm->size, fm->ioformat, NULL);

fprintf(stderr, "im here %s:%d\n",__FILE__,__LINE__);
    // pack the variable info into a contiguous buffer and attach to request
//    int set;
//    if(globals_adios_get_application_id(&set)!=2){
//        fm->var_info_buffer = create_var_info_buffer(t, fm); 
//    } 
//    else {
//        fm->var_info_buffer = (char *) malloc();
//    }

    fm->var_info_buffer = create_var_info_buffer(t, fm); 

fprintf(stderr, "im here %s:%d\n",__FILE__,__LINE__);
    //start_send(fm->s, fm->buffer, fm->size, fm->ioformat, NULL,fm->var_info_buffer); 
    start_send(fm->s, fm->buffer, fm->size, fm->ioformat, fm->alist,fm->var_info_buffer); 
    free(fm->var_info_buffer);

fprintf(stderr, "im here %s:%d\n",__FILE__,__LINE__);
    //recycle all temp arrays
    struct adios_var_struct *f;
    FMField *field = NULL;
    FMFieldList flist = fm->format->field_list;
    //for (f = vars; f != NULL; f = f->next) {
    while(f=next_var_from_list(&(fm->vars_to_write))) {
        char *full_path_name = get_full_path_name(f->name, f->path);
        field = internal_find_field(full_path_name, flist);
        free(full_path_name);

        int ndims = 0;
        struct adios_dimension_struct *dim = f->dimensions;
        int dim_are_numbers = 1;
        for(; dim != NULL; dim = dim->next) {
            if(dim->dimension.time_index != adios_flag_yes) {
                ndims ++;
            }
            if(dim->dimension.id) {
                dim_are_numbers = 0;
            } 
        }

        if(ndims) { // arrays
            if(dim_are_numbers) {
                void *temp_array;
                memcpy(&temp_array, &fm->buffer[field->field_offset], sizeof(void *));
                free(temp_array);   
            }  
        }        
        if(f->type == adios_string) {
            void *temp_array;
            memcpy(&temp_array, &fm->buffer[field->field_offset], sizeof(void *));
            free(temp_array);
        }
    } 

    reset_var_list(&(fm->vars_to_write));
    fm->var_info_buffer_size = 8;
    fm->num_vars = 0;
         
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

fprintf(stderr, "im here end simulation %d %s:%d\n", rank, __FILE__,__LINE__);
//        outputTimingInfo(buffer);

fprintf(stderr, "im here end simulation %d %s:%d\n", rank, __FILE__,__LINE__);
        send_end(fm->s);
fprintf(stderr, "im here end simulation %d %s:%d\n", rank, __FILE__,__LINE__);
        // TODO: we need a way to terminate writers:
        // we set some state flag and push to dt server side so it
        // knows it's the last piece and delivers "EOF" to reader app
        end_simulation(fm->s);
fprintf(stderr, "im here end simulation %d %s:%d\n", rank, __FILE__,__LINE__);

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
adios_datatap_init(const char *params, struct adios_method_struct *method)
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
