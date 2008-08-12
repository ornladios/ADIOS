#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include "netcdf.h"
#include "adios_types.h"
#include "adios_transport_hooks.h"
#include "adios_bp_v1.h"
#include "adios_internals.h"
#define ERR(e){if(e){printf("Error:%s\n",nc_strerror(e));return 2;}}
#define DIVIDER "========================================================\n"

struct dump_struct
{
    int do_dump;
    char * dump_var;
    enum ADIOS_FLAG host_language_fortran;
};

struct var_dim
{
    uint16_t id;
    uint64_t rank;
    int      nc_dimid;
};

void ncd_addtimedim( int ncid, int timestep_id)
{
    static int time_dimid; 
    nc_def_dim(ncid,"timestep",NC_UNLIMITED,&time_dimid);
    nc_enddef(ncid);
}

void copy_buffer(struct adios_bp_buffer_struct_v1 *dest
                ,struct adios_bp_buffer_struct_v1 *src) {

    memcpy (dest, src, sizeof(struct adios_bp_buffer_struct_v1));
}

int ncd_attr_str_ds (int ncid 
                    ,struct adios_attribute_struct_v1 * attribute
                    ,struct adios_bp_buffer_struct_v1 * ptr_buffer
                    ,int count
                    ,struct var_dim * var_dims
                    ,int var_dims_count) {
    int i;
    char fullname[255];
    char *path = attribute->path;
    char *name = attribute->name;
    char *new_path;
    int  valid,retval; 
    new_path = strdup (path);
    if ( path[0] == '/')
         new_path=new_path+1;
    if ( path[strlen(path)-1] = '/')
         new_path[strlen(new_path)-1]='\0';
    for ( i = 0; i < strlen (new_path); i++) {
        if ( new_path[i] == '[' || new_path[i] == ']' || new_path[i] == '/' || new_path[i] == '\\')
            new_path[i] = '_';
    }
    if (*new_path != '\0')
        sprintf (fullname, "%s_%s", new_path, name);
    else
        strcpy (fullname, name);
    valid = -1;
    nc_redef(ncid);
    retval=nc_inq_varid(ncid,new_path,&valid);
    //printf ("\t\tpath: %s %d\n", new_path, valid);
    if ( valid < 0)
        valid = NC_GLOBAL; 
    //nc_enddef(ncid);
    void *value = attribute->value;
    size_t len = 1; 
    enum ADIOS_DATATYES type =  attribute->type;
    struct adios_var_header_struct_v1 var_header;
    struct adios_var_payload_struct_v1 var_payload;  
    var_payload.payload = 0;
    if (attribute->is_var == adios_flag_yes) {
         for ( i = 0; i < count; i++) {
             adios_parse_var_data_header_v1 (ptr_buffer, &var_header);
             if ( var_header.id == attribute->var_id) {
                 printf("is_var: %s\n", var_header.name); 
                  struct  adios_dimension_struct_v1 * dims = var_header.dims; 
                  while (dims) {
                      if ( dims->dimension.var_id != 0 ) {
                           for (i = 0; i < var_dims_count; i++) {
                                if (var_dims [i].id == dims->dimension.var_id ){
                                    len *= var_dims [i]. rank;
			            break;
                                } 
                           }
                      }
                      else
                         len *= dims->dimension.rank;
                      dims = dims->next;
                 }
                 type = var_header.type;
                 var_payload.payload = malloc (var_header.payload_size);
                 adios_parse_var_data_payload_v1 (ptr_buffer, &var_header, &var_payload);
                 value = var_payload.payload;
              }
              else
                  adios_parse_var_data_payload_v1 (ptr_buffer, &var_header, NULL);

         }
    }
    switch (type) {
         case adios_unsigned_byte:
            retval=nc_put_att_uchar(ncid,valid,fullname,NC_BYTE,len,value);
            break;
         case adios_byte:
            retval=nc_put_att_schar(ncid,valid,fullname,NC_BYTE,len,value);
            break;
         case adios_string:
            retval=nc_put_att_text(ncid,valid,fullname, strlen(value),value);
            break;
         case adios_short:
            retval=nc_put_att_short(ncid,valid,fullname,NC_SHORT,len,value);
            break;
         case adios_integer:
            retval=nc_put_att_int(ncid,valid,fullname,NC_INT,len,value);
            break;
         case adios_long:
            retval=nc_put_att_long(ncid,valid,fullname,NC_LONG,len,value);
            break;
         case adios_real:
            retval=nc_put_att_float(ncid,valid,fullname,NC_FLOAT,len,value);
            break;
         case adios_double:
            retval=nc_put_att_double(ncid,valid,fullname,NC_DOUBLE,len,value);
            break;
         default:
            break;
     }
     ERR(retval);

    if ( var_payload.payload)
        free (var_payload.payload);
}

int ncd_dataset (int ncid
                ,struct adios_var_header_struct_v1 *ptr_var_header
                ,struct adios_var_payload_struct_v1 *ptr_var_payload
                ,struct var_dim *var_dims
                ,int var_dims_count) {

    char *name = ptr_var_header->name;
    char *path = ptr_var_header->path;
    char *new_path, dimname[255],fullname[255];
    enum ADIOS_DATATYPES type = ptr_var_header->type;
    enum ADIOS_FLAG is_dim = ptr_var_header->is_dim;
    void *val = ptr_var_payload->payload; 
    uint64_t payload_size = ptr_var_header->payload_size;
    struct adios_dimension_struct_v1 *dims = ptr_var_header->dims; 
    int maxrank = 0, i,j, valid=-1, nc_dimid=-1, retval;
    size_t rank = 0, start_dims[10],count_dims[10];
    int dimids[10];
    uint64_t dimvalue;
    new_path = strdup (path);
    if ( path[0] == '/')
         new_path=new_path+1;
    for ( i = 0; i < strlen (new_path); i++) {
        if ( new_path[i] == '[' || new_path[i] == ']' || new_path[i] == '/' || new_path[i] == '\\')
            new_path[i] = '_';
    }
    if (*new_path != '\0')
        sprintf (fullname, "%s_%s", new_path, name);
    else
        strcpy (fullname, name);
    if(dims)printf("write float %s %d %d XXXXXXX\n",fullname, dims->global_dimension.var_id, dims->global_dimension.rank);
    const char one_name[] = "one";
    static int onename_dimid = -1;
    val = ptr_var_payload->payload;
    if (dims) {
        while (dims) {
            if ( dims->global_dimension.var_id != 0 ) {
                    for (i = 0; i < var_dims_count; i++) {
                        if (var_dims [i].id == dims->global_dimension.var_id)
                            dimids[ rank]=var_dims [i]. nc_dimid;
                    }
            
                if (dims->dimension.var_id!=0 ) {
                    for (i = 0; i < var_dims_count; i++){
                        if (var_dims [i].id == dims->dimension.var_id)
                            count_dims [ rank]=var_dims [i].rank;
                    }
                }
                else
                    count_dims [ rank] = dims->dimension.rank;

                if ( dims->local_offset.var_id != 0 ) {
                    for (i = 0; i < var_dims_count; i++){
                         if (var_dims [i].id == dims->local_offset.var_id){
                             start_dims[rank]=var_dims [i]. rank; 
                             printf("\t-----global var: %d, %d in %d\n",start_dims[rank], i, var_dims_count);
                             break;
                         }
                    }
                }
                else{
                         start_dims[ rank]=dims->local_offset.rank;
                         printf("\t-----global constant: %d, %d in %d\n",start_dims[rank], i, var_dims_count);
                }
                printf(" \t global_domain: %s rank %d, start %d, count %d\n", fullname , rank, start_dims[rank],count_dims[rank]);
            }
            else if (dims->global_dimension.rank !=0 ) {
                printf(" \tconstant global_info: %s rank: %d\n",fullname, dimids[rank]);
                dimids[rank] = dims->global_dimension.rank;
                if (dims->dimension.var_id!=0 ) {
                    for (i = 0; i < var_dims_count; i++){
                         if (var_dims [i].id == dims->dimension.var_id)
                            count_dims[rank]=var_dims [i]. rank;
                    }
                }
                else
                    dimids[rank] = dims->dimension.rank;
                if (dims->local_offset.var_id!=0 ) {
                    for (i = 0; i < var_dims_count; i++){
                        if (var_dims [i].id == dims->local_offset.var_id)
                                 start_dims[rank]=var_dims [i]. rank;
                    }
                }
                else
                    start_dims[ rank]=dims->local_offset.rank;
            }
            else {
                printf(" \tno global_info: %s rank: %d\n",fullname, dimids[rank]);
                if (dims->dimension.var_id!=0 ) {
                    for (i = 0; i < var_dims_count; i++){
                        if (var_dims [i].id == dims->dimension.var_id){
                            start_dims[rank]=0;
                            count_dims[rank]=var_dims[i].rank;
                            dimids[rank]=var_dims[i].nc_dimid;
                        }
                    }
                }
                else
                {   
                    //printf ("Error, every dimension in netcdf need to have name!\n");
                    char dimname[255];
                    sprintf(dimname,"%s_%d", fullname,rank);
                    printf ("%s\n",dimname);
                    retval = nc_def_dim ( ncid, dimname, dims->dimension.rank, &nc_dimid);
                    dimids[rank]=nc_dimid;
                    count_dims[rank] = dims->dimension.rank;
                    start_dims[0] =0; 
                    ERR(retval);
                } 
            }
            ++rank;
            dims = dims->next;
        }
        val = ptr_var_payload->payload;    
        nc_redef(ncid);
        nc_inq_varid(ncid,fullname,&valid);
        printf(" \t--XXXXX----%s valid %d, start %d, count %d\n", fullname , valid, start_dims[1],count_dims[1]);
        switch(type) { 
        case adios_real:
            printf("write float %s %d\n",fullname, valid);
            if ( valid<0) {
                retval=nc_def_var(ncid,fullname,NC_FLOAT,rank,dimids,&valid);
                ERR(retval);
                printf(" \t---XXXXX---%s rank %d, start %d, count %d\n", fullname , rank, start_dims[1],count_dims[1]);
            }
            retval=nc_enddef(ncid);
            retval=nc_put_vara_float(ncid,valid,start_dims,count_dims,val);
            ERR(retval);
            break;
        case adios_double:
            if ( valid<0) 
                retval=nc_def_var(ncid,fullname,NC_DOUBLE,rank,dimids,&valid);
            retval=nc_enddef(ncid);
            retval=nc_put_vara_double(ncid,valid,start_dims,count_dims,val);
            break;
        case adios_long:
            if ( valid<0) 
                retval=nc_def_var(ncid,fullname,NC_LONG,rank,dimids,&valid);
            retval=nc_enddef(ncid);
            retval=nc_put_vara_long(ncid,valid,start_dims,count_dims,val);
            break;
        case adios_integer:
            if(valid<0)
            {
               retval=nc_def_var(ncid,fullname,NC_INT,rank,dimids,&valid);
               ERR(retval);
            }
            retval=nc_enddef(ncid);
            ERR(retval);
            printf("write integers\n");
            retval=nc_put_vara_int(ncid,valid,start_dims,count_dims,val);
            ERR(retval);
            break;
        default:
            retval=nc_enddef(ncid);
            break;
        }
    }

    else if (ptr_var_header->is_dim == adios_flag_yes) {
        for(j=0;j<var_dims_count;j++){
           if (var_dims [j].id==ptr_var_header->id) { 
               break; 
           }
        }
          
        nc_redef(ncid);
        retval=nc_inq_dimid ( ncid, fullname, &nc_dimid);
        if ( var_dims[j].rank == 0)
            return;
            printf("ncd: %s %d %d\n",fullname, var_dims[j].id, var_dims[j].rank);
        if ( nc_dimid < 0) {
           retval = nc_def_dim ( ncid, fullname, var_dims[j].rank, &nc_dimid);
           ERR(retval);
           var_dims [j].nc_dimid = nc_dimid;
        }
    }
    else {
        rank = 1;
        nc_redef(ncid);
        if (onename_dimid==-1)
        {
            retval=nc_def_dim (ncid, "one", 1, &onename_dimid);
            ERR(retval);
        }
        else {
            nc_redef(ncid);
            nc_inq_varid (ncid, fullname, &valid);
        }
        dimids[0]=onename_dimid;
        rank = 1;
    switch(type){ 
        case adios_real:
            if (valid < 0 ) {
               printf("\t ncd-scalar-real: %d %d %s\n",dimids[0],valid, fullname);
               retval=nc_def_var(ncid,fullname,NC_FLOAT,rank,dimids,&valid);
               ERR(retval);
            }
            retval=nc_enddef(ncid);
            ERR(retval);
            retval=nc_put_var_float(ncid,valid,val);
            break;
        case adios_double:
            if (valid < 0 ) {
               printf("\t ncd-scalar: %d %d %s\n",dimids[0],valid, fullname);
               retval=nc_def_var(ncid,fullname,NC_DOUBLE,rank,dimids,&valid);
               ERR(retval);
               retval=nc_enddef(ncid);
               ERR(retval);
            }
            retval=nc_enddef(ncid);
            retval=nc_put_var_double(ncid,valid,val);
            break;
        case adios_long:
            if (valid < 0 ) {
               printf("\t ncd-scalar: %d %d %s\n",dimids[0],valid, fullname);
               retval=nc_def_var(ncid,fullname,NC_LONG,rank,dimids,&valid);
               ERR(retval);
            }
            retval=nc_enddef(ncid);
            ERR(retval);
            retval=nc_def_var(ncid,fullname,NC_LONG,rank,dimids,&valid);
            retval=nc_enddef(ncid);
            retval=nc_put_var_long(ncid,valid,val);
            break;
        case adios_integer:
            if (valid < 0 ) {
               printf("\t ncd-scalar-int:  %s\n",fullname);
               //printf("\t ncd-scalar-int: %d %s\n",dimids[0],*((int *)ptr_var_payload->payload), fullname);
               retval=nc_def_var(ncid,fullname,NC_INT,rank,dimids,&valid);
               ERR(retval);
            }
            retval=nc_enddef(ncid);
            ERR(retval);
            retval=nc_put_var_int(ncid,valid,val);
            ERR(retval);
            break;
        default:
            retval=nc_enddef(ncid);
            break;
        }
    }
    return 0;
}
const char * value_to_string (enum ADIOS_DATATYPES type, void * data, uint64_t element);

int main (int argc, char ** argv)
{
    char * out_fname;
    char * var;
    int i = 0;
    int rc = 0;
    uint64_t element_size = 0;
    struct adios_bp_element_struct * element = NULL;
    struct dump_struct dump;
    if (argc < 2)
    {
        fprintf (stderr, "usage: %s <argv[1]_in> [argv[1]_out]\n"
                ,argv [0]
                );

        return -1;
    }

    if (argc > 2)
        out_fname = strdup (argv[2]);
    else 
    {
        out_fname = strdup (argv[1]);
        int size = strlen(argv[1]);
        out_fname [size-2] = 'n'; 
        out_fname [size-1] = 'c';
    }
        dump.do_dump = 0;
        dump.dump_var = 0;
    int ncid, retval;
    nc_create ( out_fname, NC_CLOBBER | NC_64BIT_OFFSET, &ncid);

    struct adios_bp_buffer_struct_v1 * b = 0;
    struct adios_bp_buffer_struct_v1 * b_0 = 0;
    uint32_t version = 0;

    b = malloc (sizeof (struct adios_bp_buffer_struct_v1));
    b_0 = malloc (sizeof (struct adios_bp_buffer_struct_v1));
    adios_buffer_struct_init (b);

    rc = adios_posix_open_read_internal (argv[1], "", b);
    if (!rc)
    {
        fprintf (stderr, "bpdump: file not found: %s\n", argv[1]);

        return -1;
    }

    adios_posix_read_version (b);
    adios_parse_version (b, &version);

    struct adios_index_process_group_struct_v1 * pg_root = 0;
    struct adios_index_process_group_struct_v1 * pg = 0;
    struct adios_index_var_struct_v1 * vars_root = 0;

    adios_posix_read_index_offsets (b);
    adios_parse_index_offsets_v1 (b);

    adios_posix_read_process_group_index (b);
    adios_parse_process_group_index_v1 (b, &pg_root);

    adios_posix_read_vars_index (b);
    adios_parse_vars_index_v1 (b, &vars_root);

    uint64_t element_num = 1;
    pg = pg_root;
    while (pg)
    {
        int var_dims_count = 0;
        struct var_dim * var_dims = 0;

        printf (DIVIDER);

        struct adios_process_group_header_struct_v1 pg_header;
        struct adios_vars_header_struct_v1 vars_header;
        struct adios_attributes_header_struct_v1 attrs_header;

        struct adios_var_header_struct_v1 var_header;
        struct adios_var_payload_struct_v1 var_payload;
        struct adios_attribute_struct_v1 attribute;

        // setup here to read the process group from (and size)
        b->read_pg_offset = pg->offset_in_file;
        if (pg->next)
        {
            b->read_pg_size =   pg->next->offset_in_file
                              - pg->offset_in_file;
        }
        else
        {
            b->read_pg_size =   b->pg_index_offset
                              - pg->offset_in_file;
        }

        adios_posix_read_process_group (b);
        adios_parse_process_group_header_v1 (b, &pg_header);
        printf ("\tTime Index Name: %d\n", pg_header.time_index_name);
#if 0
        if ( pg_header.time_index >= 0)
             ncd_addtimedim(ncid,pg_header.time_index);
#endif
        adios_parse_vars_header_v1 (b, &vars_header);

        dump.host_language_fortran = pg_header.host_language_fortran;
        int i,j;
        for (i = 0; i < vars_header.count; i++)
        {
            if(i==0)
               copy_buffer(b_0, b);
            var_payload.payload = 0;
            adios_parse_var_data_header_v1 (b, &var_header);

            if (   var_header.is_dim == adios_flag_yes
                ||    dump.do_dump
                   && !strcasecmp (dump.dump_var, var_header.name)
               )
            {
                var_payload.payload = malloc (var_header.payload_size);
                adios_parse_var_data_payload_v1 (b, &var_header, &var_payload);
            }
            else
            {
                var_payload.payload = malloc (var_header.payload_size);
                adios_parse_var_data_payload_v1 (b, &var_header, &var_payload);
                ncd_dataset(ncid,&var_header, &var_payload,var_dims,var_dims_count);
            }

            if (var_header.is_dim == adios_flag_yes)
            {
                int flag=0;
                var_dims = realloc (var_dims,   (var_dims_count + 1)
                                              * sizeof (struct var_dim)
                                   );
                
                var_dims [var_dims_count].id = -1;
                for(j=0;j<var_dims_count;j++){
                    if (var_dims [j].id==var_header.id) { 
                         var_dims [j]. rank == *(unsigned int *) var_payload.payload; 
                         flag = 1;
                         break; 
                    }
                }
                if(flag == 0){
                    var_dims [var_dims_count].id = var_header.id;
                    var_dims [var_dims_count].rank = *(unsigned int *)
                                                            var_payload.payload;
                    //printf("%s %d %d\n",var_header.name, var_dims[0].id,var_dims[0].rank);
                    var_dims_count++;
                }
                ncd_dataset(ncid,&var_header, &var_payload,var_dims,var_dims_count);
            }

            if (var_payload.payload)
            {
                free (var_payload.payload);
                var_payload.payload;
            }
            printf ("\n");
        }

        adios_parse_attributes_header_v1 (b, &attrs_header);

        for (i = 0; i < attrs_header.count; i++)
        {

            adios_parse_attribute_v1 (b, &attribute);
            ncd_attr_str_ds (ncid, &attribute, b_0, vars_header.count, var_dims, var_dims_count);
        }

        var_dims_count = 0;
        if (var_dims)
            free (var_dims);
        pg = pg->next;
    }
    printf (DIVIDER);
    printf ("End of %s\n", argv[1]);

    adios_posix_close_internal (b);
    free (b);
    free (b_0);
    nc_close (ncid);
    return 0;
}

const char * value_to_string (enum ADIOS_DATATYPES type, void * data, uint64_t element)
{
    static char s [100];
    s [0] = 0;

    switch (type)
    {
        case adios_unsigned_byte:
            sprintf (s, "%u", *(((uint8_t *) &data) + element));
            break;

        case adios_byte:
            sprintf (s, "%d", *(((int8_t *) &data) + element));
            break;

        case adios_short:
            sprintf (s, "%hd", *(((int8_t *) &data) + element));
            break;

        case adios_unsigned_short:
            sprintf (s, "%uh", *(((int8_t *) &data) + element));
            break;

        case adios_integer:
            sprintf (s, "%d", *(((int32_t *) &data) + element));
            break;

        case adios_unsigned_integer:
            sprintf (s, "%u", *(((uint32_t *) &data) + element));
            break;

        case adios_long:
            sprintf (s, "%lld", *(((int64_t *) &data) + element));
            break;

        case adios_unsigned_long:
            sprintf (s, "%llu", *(((uint64_t *) &data) + element));
            break;

        case adios_real:
            sprintf (s, "%e", *(((float *) &data) + element));
            break;

        case adios_double:
            sprintf (s, "%le", *(((double *) &data) + element));
            break;

        case adios_long_double:
            sprintf (s, "%Le", *(((long double *) &data) + element));
            break;

        case adios_string:
            sprintf (s, "%s", ((char *) data) + element);
            break;

        case adios_complex:
            sprintf (s, "(%f %f)", *((float *) data) + (element * 2 + 0)
                                 , *((float *) data) + (element * 2 + 1)
                    );
            break;

        case adios_double_complex:
            sprintf (s, "(%lf %lf)", *((double *) data) + (element * 2 + 0)
                                   , *((double *) data) + (element * 2 + 1)
                    );
            break;
    }

    return s;
}
