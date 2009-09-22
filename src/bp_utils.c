#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <sys/types.h>
#include <string.h>
#include "adios.h"
#include "bp_utils.h"
#include "adios_bp_v1.h"
#include "adios_errcodes.h"
#define BYTE_ALIGN 8
#define MINIFOOTER_SIZE 28
#define PG_MINIHEADER_SIZE 16


#ifdef DMALLOC
#include "dmalloc.h"
#endif

int adios_errno;
char aerr[ERRMSG_MAXLEN];
char *adios_errmsg() { return aerr; }

void error(int errno, char *fmt, ...) 
{
    va_list ap;
    adios_errno = errno;
    va_start(ap, fmt);
    (void) vsnprintf(aerr, ERRMSG_MAXLEN, fmt, ap);
    va_end(ap);
}

void bp_alloc_aligned (struct adios_bp_buffer_struct_v1 * b, uint64_t size)
{
        
    b->allocated_buff_ptr =  malloc (size + BYTE_ALIGN - 1);
    if (!b->allocated_buff_ptr)
    {
        error( err_no_memory, "Cannot allocate %llu bytes\n", size);

        b->buff = NULL;
        b->length = 0;

        return;
    }
    uint64_t p = (uint64_t) b->allocated_buff_ptr;
    b->buff = (char *) ((p + BYTE_ALIGN - 1) & ~(BYTE_ALIGN - 1));
    b->length = size;
}

void bp_realloc_aligned (struct adios_bp_buffer_struct_v1 * b
                            ,uint64_t size
                            )
{
    b->allocated_buff_ptr = realloc (b->allocated_buff_ptr
                                    ,size + BYTE_ALIGN - 1
                                    );
    if (!b->allocated_buff_ptr)
    {
        error( err_no_memory, "Cannot allocate %llu bytes\n", size);

        b->buff = NULL;
        b->length = 0;

        return;
    }
    uint64_t p = (uint64_t) b->allocated_buff_ptr;
    b->buff = (char *) ((p + BYTE_ALIGN - 1) & ~(BYTE_ALIGN - 1));
    b->length = size;
}

int bp_read_open (const char * filename,
          MPI_Comm comm,
          struct BP_FILE * fh)
{
    int  err;
    int  rank;

    MPI_Comm_rank (comm, &rank);

    // variable definition 
    MPI_Offset  file_size;

    // open a file by the multiple processors within the same
    // communicator
    err = MPI_File_open (comm, (char *) filename, MPI_MODE_RDONLY, 
            MPI_INFO_NULL, &(fh->mpi_fh));
    if (err != MPI_SUCCESS) {
        char e [MPI_MAX_ERROR_STRING];
        int len = 0;
        memset (e, 0, MPI_MAX_ERROR_STRING);
        MPI_Error_string (err, e, &len);
        error(err_MPI_open_error, "MPI open failed for %s: '%s'\n", filename, e);
        return adios_flag_no;
    }

    MPI_File_get_size (fh->mpi_fh, &file_size);
    fh->b->file_size = file_size;
    fh->mfooter.file_size = file_size;
     
    return 0;
}

int bp_read_minifooter (struct BP_FILE * bp_struct)
{
        struct adios_bp_buffer_struct_v1 * b = bp_struct->b;
        struct bp_minifooter * mh = &bp_struct->mfooter;

    uint64_t attrs_end = b->file_size - MINIFOOTER_SIZE;
    uint32_t test = 1;
        int r;

    MPI_Status status;

    if (!b->buff) {
        bp_alloc_aligned (b, MINIFOOTER_SIZE);
        if (!b->buff) {
            error (err_no_memory, "could not allocate %d bytes\n", MINIFOOTER_SIZE);
            return 1;
        }
        memset (b->buff, 0, MINIFOOTER_SIZE);
        b->offset = 0;
    }
    MPI_File_seek (bp_struct->mpi_fh, 
            (MPI_Offset) attrs_end, 
            MPI_SEEK_SET);

    MPI_File_read (bp_struct->mpi_fh, b->buff, MINIFOOTER_SIZE, 
            MPI_BYTE, &status);
    
    memset (&mh->pgs_index_offset, 0, MINIFOOTER_SIZE);
    memcpy (&mh->pgs_index_offset, b->buff, MINIFOOTER_SIZE);
    
    b->pg_index_offset = *(uint64_t *) (b->buff + b->offset);
    mh->pgs_index_offset = b->pg_index_offset;
    b->offset += 8;

    b->vars_index_offset = *(uint64_t *) (b->buff + b->offset);
    b->offset += 8;

    b->attrs_index_offset = *(uint64_t *) (b->buff + b->offset);
    b->offset += 8;

    /* get version id. Needs the bp->offset be pointing to the last 4 bytes of the buffer,
       so do not move this call upward in this source */
    adios_parse_version (b, &mh->version);
    mh->change_endianness = (b->change_endianness == adios_flag_yes);

    b->end_of_pgs = b->pg_index_offset;
    b->pg_size = b->vars_index_offset - b->pg_index_offset;
    b->vars_size = b->attrs_index_offset - b->vars_index_offset;
    b->attrs_size = attrs_end - b->attrs_index_offset;

    uint64_t footer_size = mh->file_size - mh->pgs_index_offset;
    bp_realloc_aligned (b, footer_size);
    MPI_File_seek (bp_struct->mpi_fh,
                        (MPI_Offset)  mh->pgs_index_offset,
                        MPI_SEEK_SET);
    MPI_File_read (bp_struct->mpi_fh, b->buff, footer_size,
            MPI_BYTE, &status);

    MPI_Get_count (&status, MPI_BYTE, &r);

    // reset the pointer to the beginning of buffer
    b->offset = 0;    
    return 0;
}

int bp_parse_pgs (struct BP_FILE * fh)
{
    struct bp_index_pg_struct_v1 ** root = &(fh->pgs_root); // need the pointer to it to malloc below
    struct adios_bp_buffer_struct_v1 * b = fh->b;
    struct bp_minifooter * mh = &(fh->mfooter);
    uint64_t i;

    memcpy (&mh->pgs_count, b->buff, PG_MINIHEADER_SIZE);

    b->offset += PG_MINIHEADER_SIZE;

    int j;

    uint64_t group_count = 0;

    char ** namelist;

    namelist = (char **) malloc(sizeof(char *)*mh->pgs_count);
    uint16_t * grpidlist = (uint16_t *) malloc(sizeof(uint16_t)*mh->pgs_count);

    for (i = 0; i < mh->pgs_count; i++) {
        uint16_t length_of_group;
        namelist[i] = 0;    
        // validate remaining length
        length_of_group = *(uint16_t *) (b->buff + b->offset);
        b->offset += 2;

        if (!*root)
        {
            *root = (struct bp_index_pg_struct_v1 *)
                malloc (sizeof(struct bp_index_pg_struct_v1));
            memset (*root, 0, sizeof(struct bp_index_pg_struct_v1));
            (*root)->next = 0;
        }
        uint16_t length_of_name;

        length_of_name = *(uint16_t *) (b->buff + b->offset);
        b->offset += 2;
        (*root)->group_name = (char *) malloc (length_of_name + 1);
        (*root)->group_name [length_of_name] = '\0';
        memcpy ((*root)->group_name, b->buff + b->offset, length_of_name);
        b->offset += length_of_name;
        
        
        if ( group_count == 0 ) { 
            namelist[group_count] = (char *) malloc (length_of_name + 1);
            strcpy (namelist[group_count], (*root)->group_name);
            ++group_count;
            grpidlist[i] = group_count-1;
        }
        else {
            for (j=0; j<group_count; j++) {
                if (!strcmp(namelist[j], (*root)->group_name)) {
                    break;
                }
            }
            if (j==group_count) {
                namelist[group_count] = (char *) malloc (length_of_name + 1);
                strcpy (namelist[group_count], (*root)->group_name);
                ++group_count;
                grpidlist[i] = group_count - 1;
            }
            else 
                grpidlist[i] = j;
                    
        }
            
        (*root)->adios_host_language_fortran =
            (*(b->buff + b->offset) == 'y' ? adios_flag_yes
             : adios_flag_no
            );
        b->offset += 1;

        (*root)->process_id = *(uint32_t *) (b->buff + b->offset);
        b->offset += 4;

        length_of_name = *(uint16_t *) (b->buff + b->offset);
        b->offset += 2;
        (*root)->time_index_name = (char *) malloc (length_of_name + 1);
        (*root)->time_index_name [length_of_name] = '\0';
        memcpy ((*root)->time_index_name, b->buff + b->offset, length_of_name);
        b->offset += length_of_name;

        (*root)->time_index = *(uint32_t *) (b->buff + b->offset);
        b->offset += 4;

        (*root)->offset_in_file = *(uint64_t *) (b->buff + b->offset);
        b->offset += 8;
        if (i == mh->pgs_count-1)
            mh->time_steps = (*root)->time_index;

        root = &(*root)->next;
    }

    /*
    root = &(fh->pgs_root);
    for (i = 0; i < mh->pgs_count; i++) {
        printf("%d\tpg pid=%d addr=%x next=%x\n",i, (*root)->process_id, *root, (*root)->next);
        root = &(*root)->next;
    }
    */

    uint64_t * pg_offsets = 0; 
    uint32_t * pg_pids = 0; 
    uint32_t *** time_index = 0;
    pg_offsets = (uint64_t *) 
        malloc (sizeof(uint64_t)*mh->pgs_count);
    pg_pids = (uint32_t *)
        malloc (sizeof(uint32_t)*mh->pgs_count);
    // time_index[0]: record which pg to start from per timestep per group
    // time_index[1]: record the # of pgs per timesteps per group 
    time_index = (uint32_t ***) malloc (sizeof(uint32_t **)*2);

    for (j=0;j<2;j++) {
        time_index[j] = (uint32_t **) 
            malloc (sizeof(uint32_t*)*group_count);
        //printf ("### time_index[%d]=%x  group_count=%d  #pgs=%d #ts=%d\n", j, time_index[j], group_count, mh->pgs_count,  mh->time_steps);
        for (i=0;i<group_count;i++) {
            if (mh->pgs_count < mh->time_steps) {
                /* FIXME: when can this happen? 
                   pgs = time_steps * number of writers, if there is 1 group only
                */
                time_index[j][i] = (uint32_t *) 
                    malloc (sizeof(uint32_t)*mh->pgs_count);
            } else {
                time_index[j][i] = (uint32_t *) 
                    malloc (sizeof(uint32_t)*mh->time_steps);
            }
        }
    }

    root = &(fh->pgs_root);
    uint32_t time_id, tidx_start, tidx_stop;
    uint64_t grpid = grpidlist[0];
    uint32_t pg_time_count = 0, first_pg;
    tidx_start = (*root)->time_index;
    time_id = tidx_start;
    first_pg = 0; /* The first pg for a given timestep and group */
    for (i = 0; i < mh->pgs_count; i++) {
        pg_pids [i] = (*root)->process_id;
        pg_offsets [i] = (*root)->offset_in_file;
        //printf ("### root->time_index=%d,  time_id=%d\n", (*root)->time_index, time_id);
        if ((*root)->time_index == time_id) {
            /* processing still the same timestep */
            if (grpid == grpidlist[i]) {
                /* processing still the same group */
                /* FIXME: is this the order in the file? time..groups or group..times? */
                pg_time_count += 1;
            } else {
                /* changing group: pg_time_count is for the current group the number of pgs of the same time */
                time_index [0][grpid][time_id-tidx_start] = first_pg;
                time_index [1][grpid][time_id-tidx_start] = pg_time_count;
                //printf ("#-- time_index[0][%d][%d]=%d\n", grpid, time_id-tidx_start, first_pg);
                //printf ("#   time_index[1][%d][%d]=%d\n", grpid, time_id-tidx_start, pg_time_count);
                grpid = grpidlist [i];    
                pg_time_count = 1;
                first_pg = i; // new group starts from this pg
            }
        }    
        else {
            /* change in timestep */
            if (group_count == 1) {
                /* single group in file (the most frequent case) */
                time_index [0][grpid][time_id-tidx_start] = first_pg;
                time_index [1][grpid][time_id-tidx_start] = pg_time_count;
                //printf ("### time_index[0][%d][%d]=%d\n", grpid, time_id-tidx_start, first_pg);
                //printf ("    time_index[1][%d][%d]=%d\n", grpid, time_id-tidx_start, pg_time_count);
                first_pg = i;
            }
            else {    
                if (grpid == grpidlist[i]) {
                    pg_time_count += 1;
                } else {
                    time_index [0][grpid][time_id-tidx_start] = first_pg;
                    time_index [1][grpid][time_id-tidx_start] = pg_time_count;
                    //printf ("#.. time_index[0][%d][%d]=%d\n", grpid, time_id-tidx_start, first_pg);
                    //printf ("    time_index[1][%d][%d]=%d\n", grpid, time_id-tidx_start, pg_time_count);
                    grpid = grpidlist [i];    
                    first_pg = i;
                }
            }
            time_id = (*root)->time_index;
            pg_time_count = 1;
        }
        root = &(*root)->next;
    }
    /* set last grp/time count to complete the above procedure */
    time_index [0][grpid][time_id-tidx_start] = first_pg;
    time_index [1][grpid][time_id-tidx_start] = pg_time_count;
    //printf ("#   time_index[0][%d][%d]=%d\n", grpid, time_id-tidx_start, first_pg);
    //printf ("    time_index[1][%d][%d]=%d\n", grpid, time_id-tidx_start, pg_time_count);
    tidx_stop = time_id;


    /* Copy group_count strings from namelist and then free up namelist */
    char ** grp_namelist;
 
    grp_namelist = (char **) malloc (sizeof(char*) * group_count);
    for (i=0;i<group_count;i++) {
        //grp_namelist[i] = (char *) malloc (strlen(namelist[i])+1);
        //strcpy(grp_namelist[i],namelist[i]);
        grp_namelist[i] = namelist[i];
    }
    free(namelist);
    
    // here we need:
    //        grp_namelist [ngroup]
    //    time_index   [2][ngroup][nprocess]
    //    pg_offsets   [npgs]

    free (pg_pids);

    fh->gvar_h = (struct BP_GROUP_VAR *) malloc (sizeof(struct BP_GROUP_VAR));
    fh->gvar_h->group_count = group_count;
    fh->gvar_h->pg_offsets = pg_offsets;
    fh->gvar_h->namelist = grp_namelist; 
    fh->gvar_h->time_index = time_index; 
    fh->gvar_h->group_id = 0;
    fh->gvar_h->var_offsets = 0;
    fh->gvar_h->var_namelist = 0;
    fh->gvar_h->var_counts_per_group = 0;

    fh->gattr_h = (struct BP_GROUP_ATTR *) malloc (sizeof(struct BP_GROUP_ATTR));
    fh->gattr_h->group_count = group_count;
    fh->gattr_h->namelist = grp_namelist;
    fh->gattr_h->group_id = 0;
    fh->gattr_h->attr_offsets = 0;
    fh->gattr_h->attr_namelist = 0;
    fh->gattr_h->attr_counts_per_group = 0;

    fh->tidx_start = tidx_start; 
    fh->tidx_stop= tidx_stop; 

    free(grpidlist);
    return;
}

int bp_parse_attrs (struct BP_FILE * fh)
{
    struct adios_bp_buffer_struct_v1 * b = fh->b;
    struct adios_index_attribute_struct_v1 ** attrs_root = &(fh->attrs_root);
    struct bp_minifooter * mh = &(fh->mfooter);
    struct adios_index_attribute_struct_v1 ** root;
    int i;

    if (b->length - b->offset < VARS_MINIHEADER_SIZE) {
        fprintf (stderr, "adios_parse_attrs_index_v1 requires a buffer "
                "of at least %d bytes.  Only %llu were provided\n"
                ,VARS_MINIHEADER_SIZE
                ,b->length - b->offset
            );

        return 1;
    }

    root = attrs_root;

    memcpy (&mh->attrs_count, b->buff + b->offset, 2);
    b->offset += 2;

    memcpy (&mh->attrs_length, b->buff + b->offset, 8);
    b->offset += 8;

    for (i = 0; i < mh->attrs_count; i++) {
        if (!*root)
        {
            *root = (struct adios_index_attribute_struct_v1 *)
                      malloc (sizeof (struct adios_index_attribute_struct_v1));
            (*root)->next = 0;
        }
        uint8_t flag;
        uint32_t attr_entry_length;
        uint16_t len;
        uint64_t characteristics_sets_count;
        uint64_t type_size;

        attr_entry_length = *(uint32_t *) (b->buff + b->offset);
        b->offset += 4;

        (*root)->id = *(uint16_t *) (b->buff + b->offset);
        b->offset += 2;

        len = *(uint16_t *) (b->buff + b->offset);
        b->offset += 2;

        (*root)->group_name = (char *) malloc (len + 1);
        (*root)->group_name [len] = '\0';
        strncpy ((*root)->group_name, b->buff + b->offset, len);
        b->offset += len;

        len = *(uint16_t *) (b->buff + b->offset);
        b->offset += 2;
        (*root)->attr_name = (char *) malloc (len + 1);
        (*root)->attr_name [len] = '\0';
        strncpy ((*root)->attr_name, b->buff + b->offset, len);
        b->offset += len;

        len = *(uint16_t *) (b->buff + b->offset);
        b->offset += 2;
        (*root)->attr_path = (char *) malloc (len + 1);
        (*root)->attr_path [len] = '\0';
        strncpy ((*root)->attr_path, b->buff + b->offset, len);
        b->offset += len;

        flag = *(b->buff + b->offset);
        (*root)->type = (enum ADIOS_DATATYPES) flag;

        b->offset += 1;
        type_size = adios_get_type_size ((*root)->type, "");
        if (type_size == -1)
        {
            type_size = 4;
            (*root)->type = adios_integer;
        }

        characteristics_sets_count = *(uint64_t *) (b->buff + b->offset);
        (*root)->characteristics_count = characteristics_sets_count;
        (*root)->characteristics_allocated = characteristics_sets_count;
        b->offset += 8;

        // validate remaining length: offsets_count * (8 + 2 * (size of type))
        uint64_t j;
        (*root)->characteristics = malloc (characteristics_sets_count
                       * sizeof (struct adios_index_characteristic_struct_v1)
                      );
        memset ((*root)->characteristics, 0
               ,  characteristics_sets_count
                * sizeof (struct adios_index_characteristic_struct_v1)
               );

        for (j = 0; j < characteristics_sets_count; j++)
        {
            uint8_t characteristic_set_count;
            uint32_t characteristic_set_length;
            uint8_t item = 0;

            characteristic_set_count = (uint8_t) *(b->buff + b->offset);
            b->offset += 1;

            characteristic_set_length = *(uint32_t *) (b->buff + b->offset);
            b->offset += 4;

            while (item < characteristic_set_count)
            {
                uint8_t flag;
                enum ADIOS_CHARACTERISTICS c;
                flag = *(b->buff + b->offset);
                c = (enum ADIOS_CHARACTERISTICS) flag;
                b->offset += 1;

                switch (c)
                {
                    case adios_characteristic_value:
                    {
                        uint16_t data_size;
                        void * data = 0;

                        if ((*root)->type == adios_string)
                        {
                            data_size = *(uint16_t *) (b->buff + b->offset);
                            b->offset += 2;
                        }
                        else
                        {
                            data_size = adios_get_type_size ((*root)->type, "");
                        }

                        data = malloc (data_size + 1);
                        ((char *) data) [data_size] = '\0';

                        if (!data)
                        {
                            fprintf (stderr, "cannot allocate %d bytes to "
                                             "copy scalar %s\n"
                                    ,data_size
                                    ,(*root)->attr_name
                                    );

                            return 1;
                        }

                        switch ((*root)->type)
                        {
                            case adios_byte:
                            case adios_short:
                            case adios_integer:
                            case adios_long:
                            case adios_unsigned_byte:
                            case adios_unsigned_short:
                            case adios_unsigned_integer:
                            case adios_unsigned_long:
                            case adios_real:
                            case adios_double:
                            case adios_long_double:
                            case adios_complex:
                            case adios_double_complex:
                                memcpy (data, (b->buff + b->offset), data_size);
                                b->offset += data_size;
                                break;

                            case adios_string:
                                memcpy (data, (b->buff + b->offset), data_size);
                                b->offset += data_size;
                                break;

                            default:
                                free (data);
                                data = 0;
                                break;
                        }

                        (*root)->characteristics [j].value = data;

                        break;
                    }

                    case adios_characteristic_offset:
                    {
                        uint64_t size = adios_get_type_size ((*root)->type, "");
                        (*root)->characteristics [j].offset =
                                            *(uint64_t *) (b->buff + b->offset);
                        b->offset += 8;

                        break;
                    }

                    case adios_characteristic_payload_offset:
                    {
                        uint64_t size = adios_get_type_size ((*root)->type, "");
                        (*root)->characteristics [j].payload_offset =
                                            *(uint64_t *) (b->buff + b->offset);
                        b->offset += 8;

                        break;
                    }

                    case adios_characteristic_var_id:
                    {
                        (*root)->characteristics [j].var_id =
                                            *(uint16_t *) (b->buff + b->offset);
                        b->offset += 2;

                        break;
                    }

                }
                item++;
            }
        }

        root = &(*root)->next;
    }

    root = attrs_root;
    uint16_t * attr_counts_per_group;
    uint16_t *  attr_gids;
    uint64_t ** attr_offsets;
    char ** attr_namelist;
    int grpid, j,cnt;

    attr_counts_per_group = (uint16_t *)
        malloc (sizeof(uint16_t) * fh->gattr_h->group_count);
    memset (attr_counts_per_group, 0, fh->gattr_h->group_count * sizeof(uint16_t));
    attr_gids = (uint16_t *) malloc (sizeof(uint16_t ) * mh->attrs_count);
    attr_namelist = (char **)malloc (sizeof(char*) * mh->attrs_count);

    attr_offsets = (uint64_t **) malloc (sizeof(uint64_t *) * mh->attrs_count);
    memset (attr_offsets, 0, mh->attrs_count * sizeof(uint64_t));

    for (i = 0; i < mh->attrs_count; i++) {
        struct adios_index_characteristic_dims_struct_v1 * pdims;
        for (grpid = 0; grpid < fh->gattr_h->group_count; grpid++) {
            if (!strcmp((*root)->group_name, fh->gattr_h->namelist[grpid])) {
                attr_counts_per_group [grpid]++;
                attr_gids [i] = grpid;
                break;
            }
        }
        // Full name of attributes: concatenate attr_path and attr_name
        // Always have / in the beginning of the full name
        if (strcmp ((*root)->attr_path,"/")) {
            attr_namelist [i] = (char *) malloc ( strlen((*root)->attr_name)
                    +strlen((*root)->attr_path) + 1 + 1
                    );
            strcpy(attr_namelist[i], (*root)->attr_path);
        }
        else {
            attr_namelist [i] = (char *) malloc ( strlen((*root)->attr_name)+1+1);
                        attr_namelist[i][0] = '\0';
        }
        strcat(attr_namelist[i], "/");
        strcat(attr_namelist[i], (*root)->attr_name);

        attr_offsets[i] = (uint64_t *) malloc (
                sizeof(uint64_t)*(*root)->characteristics_count);
        for (j=0;j < (*root)->characteristics_count;j++) {
            attr_offsets[i][j] = (*root)->characteristics [j].offset;
        }

        pdims = &(*root)->characteristics [0].dims;
        cnt = pdims->count;
        root = &(*root)->next;
    }
    //here is the asssumption that attr_gids is linearly increased
    free(attr_gids);

    fh->gattr_h->attr_namelist = attr_namelist;
    fh->gattr_h->attr_counts_per_group = attr_counts_per_group;
    fh->gattr_h->attr_offsets = attr_offsets;

}
int bp_parse_vars (struct BP_FILE * fh)
{
    struct adios_bp_buffer_struct_v1 * b = fh->b;
    struct adios_index_var_struct_v1 ** vars_root = &(fh->vars_root);
    struct bp_minifooter * mh = &(fh->mfooter);

    struct adios_index_var_struct_v1 ** root;

    if (b->length - b->offset < VARS_MINIHEADER_SIZE) {
        fprintf (stderr, "adios_parse_vars_index_v1 requires a buffer "
                "of at least %d bytes.  Only %llu were provided\n"
                ,VARS_MINIHEADER_SIZE
                ,b->length - b->offset
            );

        return 1;
    }

    root = vars_root;

    memcpy (&mh->vars_count, b->buff + b->offset, 2);
    b->offset += VARS_MINIHEADER_SIZE;
    // validate remaining length    
    int i;
    for (i = 0; i < mh->vars_count; i++) {
        if (!*root) {
            *root = (struct adios_index_var_struct_v1 *)
                malloc (sizeof (struct adios_index_var_struct_v1));
            (*root)->next = 0;
        }
        uint8_t flag;
        uint32_t var_entry_length;
        uint16_t len;
        uint64_t characteristics_sets_count;
        uint64_t type_size;

        var_entry_length = *(uint32_t *) (b->buff + b->offset);
        b->offset += 4;

        (*root)->id = *(uint16_t *) (b->buff + b->offset);
        b->offset += 2;

        len = *(uint16_t *) (b->buff + b->offset);
        b->offset += 2;
        (*root)->group_name = (char *) malloc (len + 1);
        (*root)->group_name [len] = '\0';
        strncpy ((*root)->group_name, b->buff + b->offset, len);
        b->offset += len;

        len = *(uint16_t *) (b->buff + b->offset);
        b->offset += 2;
        (*root)->var_name = (char *) malloc (len + 1);
        (*root)->var_name [len] = '\0';
        strncpy ((*root)->var_name, b->buff + b->offset, len);
        b->offset += len;

        len = *(uint16_t *) (b->buff + b->offset);
        b->offset += 2;
        (*root)->var_path = (char *) malloc (len + 1);
        (*root)->var_path [len] = '\0';
        strncpy ((*root)->var_path, b->buff + b->offset, len);
        b->offset += len;
        flag = *(b->buff + b->offset);
        (*root)->type = (enum ADIOS_DATATYPES) flag;
        b->offset += 1;
        type_size = bp_get_type_size ((*root)->type, "");

        characteristics_sets_count = *(uint64_t *) (b->buff + b->offset);
        (*root)->characteristics_count = characteristics_sets_count;
        (*root)->characteristics_allocated = characteristics_sets_count;
        b->offset += 8;

        // validate remaining length: offsets_count * 
        // (8 + 2 * (size of type))
        (*root)->characteristics = malloc (characteristics_sets_count
            * sizeof (struct adios_index_characteristic_struct_v1)
            );
        memset ((*root)->characteristics, 0
            ,  characteristics_sets_count
            * sizeof (struct adios_index_characteristic_struct_v1)
               );

        uint64_t j;
        for (j = 0; j < characteristics_sets_count; j++)
        {
            uint8_t characteristic_set_count;
            uint32_t characteristic_set_length;
            uint8_t item = 0;

            characteristic_set_count = (uint8_t) 
                           *(b->buff + b->offset);
            b->offset += 1;

            characteristic_set_length = *(uint32_t *) 
                            (b->buff + b->offset);
            b->offset += 4;
                
            while (item < characteristic_set_count) {
                bp_parse_characteristics (b, root, j);
                item++;
            }

        }

        root = &(*root)->next;
    }
    
    root = vars_root;
    uint16_t * var_counts_per_group;
     uint16_t *  var_gids;
    uint64_t ** var_offsets;
    char ** var_namelist;
    int grpid, j,cnt;

    var_counts_per_group = (uint16_t *) 
        malloc (sizeof(uint16_t)*fh->gvar_h->group_count);
    memset ( var_counts_per_group, 0, fh->gvar_h->group_count*sizeof(uint16_t));
    var_gids = (uint16_t *) malloc (sizeof(uint16_t )*mh->vars_count);
    var_namelist = (char **)malloc(sizeof(char*)*mh->vars_count);

    var_offsets = (uint64_t **) malloc (sizeof(uint64_t *)*mh->vars_count);
    memset ( var_offsets, 0, mh->vars_count*sizeof(uint64_t));

    for (i = 0; i < mh->vars_count; i++) {
        struct adios_index_characteristic_dims_struct_v1 * pdims;
        for (grpid=0;grpid<fh->gvar_h->group_count;grpid++) {
            if (!strcmp((*root)->group_name,fh->gvar_h->namelist[grpid])) {
                var_counts_per_group [grpid]++;
                var_gids [i] = grpid;
                break;
            }
        }
                // Full name of variable: concatenate var_path and var_name
                // Always have / in the beginning of the full name
        if (strcmp ((*root)->var_path,"/")) {
            var_namelist [i] = (char *) malloc ( strlen((*root)->var_name)
                    +strlen((*root)->var_path) + 1 + 1   // extra / and ending \0
                    );
            strcpy(var_namelist[i], (*root)->var_path);
        }
        else {
            var_namelist [i] = (char *) malloc ( strlen((*root)->var_name)+1+1);
            var_namelist[i][0] = '\0';
        }
        strcat(var_namelist[i], "/");
        strcat(var_namelist[i], (*root)->var_name);

        var_offsets[i] = (uint64_t *) malloc (
                sizeof(uint64_t)*(*root)->characteristics_count);
        for (j=0;j < (*root)->characteristics_count;j++) {
            var_offsets[i][j] = (*root)->characteristics [j].offset;
        }

        pdims = &(*root)->characteristics [0].dims;
        cnt = pdims->count;
        root = &(*root)->next;
    }
    //here is the asssumption that var_gids is linearly increased
     free( var_gids);
    fh->gvar_h->var_namelist = var_namelist;
    fh->gvar_h->var_counts_per_group=var_counts_per_group;
    fh->gvar_h->var_offsets = var_offsets;
    return 0;
    
    // here we need
    // number of group
    // number of vars in each group
    // the offsets for each vars
}

int bp_parse_characteristics (struct adios_bp_buffer_struct_v1 * b,
                    struct adios_index_var_struct_v1 ** root,
                  uint64_t j)
{
    uint8_t flag;
    enum ADIOS_CHARACTERISTICS c;

    flag = *(b->buff + b->offset);
    c = (enum ADIOS_CHARACTERISTICS) flag;
    b->offset += 1;

    switch (c) {

        case adios_characteristic_value:
        case adios_characteristic_min:
        case adios_characteristic_max:
        {
            void * data = 0;

            int data_size = 0;
            if ((*root)->type == adios_string) {
                data_size = *(uint16_t *) (b->buff + b->offset);
                b->offset += 2;
            }
            else {
                data_size = bp_get_type_size ((*root)->type, "");
            }

            bp_get_characteristics_data (&data, b->buff+b->offset,
                             data_size, (*root)->type);

            switch (c) {
                case adios_characteristic_value:
                    (*root)->characteristics [j].value = data;
                    break;

                case adios_characteristic_min:
                    (*root)->characteristics [j].min = data;
                    break;

                case adios_characteristic_max:
                    (*root)->characteristics [j].max = data;
                    break;
            }

            b->offset += data_size;

            break;
        }
        case adios_characteristic_offset: 
        {
            uint64_t size = bp_get_type_size ((*root)->type, "");
            (*root)->characteristics [j].offset =
                *(uint64_t *) (b->buff + b->offset);
            b->offset += 8;

            break;
        }
        case adios_characteristic_payload_offset: 
        {
            uint64_t size = bp_get_type_size ((*root)->type, "");
            (*root)->characteristics [j].payload_offset =
                *(uint64_t *) (b->buff + b->offset);
            b->offset += 8;

            break;
        }
        case adios_characteristic_dimensions:
        {
            uint16_t dims_length;

            (*root)->characteristics [j].dims.count =
                *(uint8_t *) (b->buff + b->offset);
            b->offset += 1;

            dims_length = *(uint16_t *) (b->buff + b->offset);
            b->offset += 2;

            (*root)->characteristics [j].dims.dims = (uint64_t *)
                malloc (dims_length);
            memcpy ((*root)->characteristics [j].dims.dims
                    ,(b->buff + b->offset)
                    ,dims_length
                   );
            b->offset += dims_length;

                        break;
        }
            default:
                fprintf (stderr, "Unknown characteristic:%d. skipped.\n", c);
                break;
    }
}

int bp_get_characteristics_data (void ** ptr_data,
                 void * buffer,
                 int data_size,
                 enum ADIOS_DATATYPES type)
{
    void * data = 0;
    switch (type) {
        case adios_byte:
        case adios_short:
        case adios_integer:
        case adios_long:
        case adios_unsigned_byte:
        case adios_unsigned_short:
        case adios_unsigned_integer:
        case adios_unsigned_long:
        case adios_real:
        case adios_double:
        case adios_long_double:
        case adios_complex:
        case adios_double_complex:
        {
            data = malloc (data_size);

            if (!data)
            {
                fprintf (stderr, "cannot allocate %d bytes "
                        ,data_size
                    );

                return 1;
            }

            memcpy (data, buffer, data_size);
            break;

        }
        case adios_string:
        {
            data = malloc (data_size + 1);

            if (!data)
            {
                fprintf (stderr, "cannot allocate %d bytes "
                        ,data_size
                    );

                return 1;
            }

            ((char *) data) [data_size] = '\0';
            memcpy (data, buffer, data_size);
            break;
        }
        default:
        {
            data = 0;
            break;
        }
    }
    *ptr_data = data;
    return 0;
}

/*
void bp_grouping ( struct BP_FILE * fh_p,
           uint64_t * gh_p)
{
    struct BP_FILE * fh = (struct BP_FILE *) fh_p;
    struct bp_index_pg_struct_v1 * pg_root = fh->pgs_root;
    struct bp_minifooter * mh = &fh->mfooter;    
    int i, j; 
    uint32_t time_id;
    uint64_t pg_time_count = 0;
    uint64_t * pg_offsets = (uint64_t *) 
        malloc (sizeof(uint64_t)*mh->pgs_count);
    uint32_t * pg_pids = (uint32_t *) 
        malloc (sizeof(uint32_t)*mh->pgs_count);
    uint64_t * time_index = (uint64_t *) 
        malloc (sizeof(uint64_t)*mh->time_steps);
    time_id = pg_root->time_index;

    uint16_t group_count = 0;
     
    for (i = 0; i < mh->pgs_count; i++) {
        pg_pids [i] = pg_root->process_id;
        pg_offsets [i] = pg_root->offset_in_file;
        if (pg_root->time_index == time_id) {
            pg_time_count += 1;
        }    
        else {
            time_index [time_id-1] = pg_time_count;
            time_id = pg_root->time_index;
            pg_time_count = 1;
        }
    
        pg_root = pg_root->next;
    }

    pg_root = fh->pgs_root;

    time_index [time_id-1] = pg_time_count;
    time_id = 0;
    for (i = 0; i < mh->time_steps; i++) {
        if (i > 0) 
            time_id += time_index[i-1];

    }
    struct adios_index_var_struct_v1 * vars = fh->vars_root;
    int vars_cnt = 0;
    while (vars) {
        if (!strcmp(vars->group_name, pg_root->group_name)) {
             printf ("%s %s %d %d %d %d\n",vars->var_name,
                vars->group_name,
                vars->characteristics_count,
                vars->characteristics->dims.count,
                vars->characteristics->var_id,
                vars->id
                );
            ++vars_cnt;    
        }
        vars = vars->next;    
    }
     printf("cnt=%d \n",vars_cnt);

    return;
}
*/

int bp_read_pgs (struct BP_FILE * bp_struct)
{
        struct adios_bp_buffer_struct_v1 * b = bp_struct->b;
    int r = 0;
    MPI_Status status;
    // init buffer for pg reading
    bp_realloc_aligned (b, b->pg_size);
    b->offset = 0;

    if (sizeof (char *) == 4) { 
        MPI_File_seek (bp_struct->mpi_fh, 
                (MPI_Offset) b->pg_index_offset, 
                MPI_SEEK_SET);

        MPI_File_read (bp_struct->mpi_fh, b->buff, 
                b->pg_size, MPI_BYTE, &status);
        MPI_Get_count (&status, MPI_BYTE, &r);
    }
    else { 
        MPI_File_seek (bp_struct->mpi_fh, 
                (MPI_Offset) b->pg_index_offset, 
                MPI_SEEK_SET);

        MPI_File_read (bp_struct->mpi_fh, b->buff, 
                b->pg_size, MPI_BYTE, &status);
        MPI_Get_count (&status, MPI_BYTE, &r);
    }
    if (r != b->pg_size)
        fprintf (stderr, "could not read %llu bytes. read only: %llu\n",
                b->pg_size, r);

    return 0;
}

int bp_read_vars (struct BP_FILE * bp_struct)
{
        struct adios_bp_buffer_struct_v1 * b = bp_struct->b;
    int r = 0;
    MPI_Status status;

    //init buffer for vars reading
    bp_realloc_aligned (b, b->vars_size);
    b->offset = 0;

    if (sizeof (char *) == 4) { 
        MPI_File_seek (bp_struct->mpi_fh, 
                (MPI_Offset) b->vars_index_offset, 
                MPI_SEEK_SET);

        MPI_File_read (bp_struct->mpi_fh, b->buff, 
                b->vars_size, MPI_BYTE, &status);
        MPI_Get_count (&status, MPI_BYTE, &r);
    }
    else { 
        MPI_File_seek (bp_struct->mpi_fh,
                (MPI_Offset) b->vars_index_offset, 
                MPI_SEEK_SET);

        MPI_File_read (bp_struct->mpi_fh, b->buff, 
                b->vars_size, MPI_BYTE, &status);
        MPI_Get_count (&status, MPI_BYTE, &r);
    }
    if (r != b->vars_size)
        fprintf (stderr, "could not read %llu bytes. read only: %llu\n",
                b->vars_size, r);

    return 0;
}

void print_pg_index (struct bp_index_pg_struct_v1 * pg_root,
        struct bp_minifooter * mh)
{
    int i, j; 
    uint32_t time_id;
    uint64_t pg_time_count = 0;
    uint64_t * pg_offsets = (uint64_t *) 
        malloc (sizeof(uint64_t)*mh->pgs_count);
    uint32_t * pg_pids = (uint32_t *) 
        malloc (sizeof(uint32_t)*mh->pgs_count);
    uint64_t * time_index = (uint64_t *) 
        malloc (sizeof(uint64_t)*mh->time_steps);
    time_id = pg_root->time_index;
    for (i = 0; i < mh->pgs_count; i++)
    {
        pg_pids [i] = pg_root->process_id;
        pg_offsets [i] = pg_root->offset_in_file;
        if (pg_root->time_index == time_id)
        {
            pg_time_count += 1;
        }    
        else {
            time_index [time_id-1] = pg_time_count;
            time_id = pg_root->time_index;
            pg_time_count = 1;
        }    
        pg_root = pg_root->next;
    }

    time_index [time_id-1] = pg_time_count;
    time_id = 0;
    for (i = 0; i < mh->time_steps; i++) {
        if (i > 0) 
            time_id += time_index[i-1];
    }
}

void print_vars_index_top (struct adios_index_var_struct_v1 * vars_root)
{
    printf("Variables (group) :\n");    
    while (vars_root) {
        if (!strcmp (vars_root->var_path, "/")) {
            printf ("\t %s", 
                vars_root->var_name
                );
        }
        else {
            printf ("\t %s/%s",
                vars_root->var_path,
                vars_root->var_name
                );
        }
        
        int j, cnt;
        struct adios_index_characteristic_dims_struct_v1 * pdims;
        pdims = &vars_root->characteristics [0].dims;
        cnt = pdims->count;
        if (cnt != 0) {
            printf (" (");
            for (j = 0; j < cnt; j++) { 
                if (j>0)
                    printf (", ");
                if (pdims->dims [j*3 + 1] != 0) {
                    printf ("%llu", pdims->dims [j*3 + 1]);
                }
                else {
                    printf ("%llu", pdims->dims [j*3 + 0]);
                }
            }
            printf (")");
        }
        printf("\n");
        vars_root = vars_root->next;
    }
}

void print_vars_index (struct adios_index_var_struct_v1 * vars_root)
{
    while (vars_root) {
        if (!strcmp (vars_root->var_path, "/")) {
            printf ("Var (Group) [ID]: /%s (%s) [%d]\n", vars_root->var_name
                   ,vars_root->group_name, vars_root->id
                   );
        }
        else {
            printf ("Var (Group) [ID]: %s/%s (%s) [%d]\n", vars_root->var_path
                   ,vars_root->var_name, vars_root->group_name, vars_root->id
                   );
        }
        printf ("\tVars Characteristics: %llu\n"
               ,vars_root->characteristics_count
               );
        uint64_t i;
        for (i = 0; i < vars_root->characteristics_count; i++) {
            printf ("\t\tOffset(%llu)", vars_root->characteristics [i].offset);
            if (vars_root->characteristics [i].min)
            {
                printf ("\t\tMin(%s)", value_to_string (vars_root->type
                                           ,vars_root->characteristics [i].min
                                           )
                       );
            }
            if (vars_root->characteristics [i].max)
            {
                printf ("\t\tMax(%s)", value_to_string (vars_root->type
                                           ,vars_root->characteristics [i].max
                                           )
                       );
            }
            if (vars_root->characteristics [i].value)
            {
                printf ("\t\tValue(%s)", value_to_string (vars_root->type
                                         ,vars_root->characteristics [i].value
                                         )
                       );
            }
            if (vars_root->characteristics [i].dims.count != 0) {
                int j;

                printf ("\t\tDims (l:g:o): (");
                for (j = 0; j < vars_root->characteristics [i].dims.count; j++)
                {
                    if (j != 0)
                        printf (",");
                    if (  vars_root->characteristics [i].dims.dims [j * 3 + 1]
                        != 0
                       )
                    {
                        printf ("%llu:%llu:%llu"
                         ,vars_root->characteristics [i].dims.dims [j * 3 + 0]
                         ,vars_root->characteristics [i].dims.dims [j * 3 + 1]
                         ,vars_root->characteristics [i].dims.dims [j * 3 + 2]
                               );
                    }
                    else
                    {
                        printf ("%llu"
                         ,vars_root->characteristics [i].dims.dims [j * 3 + 0]
                               );
                    }
                }
                printf (")");
            }
            printf ("\n");
        }

        vars_root = vars_root->next;
    }
}

/* Return 1 if a < b wrt. the given type, otherwise 0 */
int adios_lt(int type, void *a, void *b)
{
    double ar, ai, br, bi;
    long double ar2, ai2, br2, bi2;
    switch (type)
    {
        case adios_unsigned_byte:
            return *((uint8_t *) a) < *((uint8_t *) b);
        case adios_byte:
            return *((int8_t *) a) < *((int8_t *) b);
        case adios_short:
            return *((int16_t *) a) < *((int16_t *) b);
        case adios_unsigned_short:
            return *((uint16_t *) a) < *((uint16_t *) b);
        case adios_integer:
            return *((int32_t *) a) < *((int32_t *) b);
        case adios_unsigned_integer:
            return *((uint32_t *) a) < *((uint32_t *) b);
        case adios_long:
            return *((int64_t *) a) < *((int64_t *) b);
        case adios_unsigned_long:
            return *((uint64_t *) a) < *((uint64_t *) b);
        case adios_real:
            return *((float *) a) < *((float *) b);
        case adios_double:
            return *((double *) a) < *((double *) b);
        case adios_long_double:
            return *((long double *) a) < *((long double *) b);
        case adios_string:
            return ( strcmp( (char *) a, (char *) b) < 0);
        case adios_complex:
            ar = (double) ((float *) a)[0];
            ai = (double) ((float *) a)[1];
            br = (double) ((float *) b)[0];
            bi = (double) ((float *) b)[1];
            return ( ar*ar+ai*ai < br*br+bi*bi );
        case adios_double_complex:
            ar2 = (long double) ((double *) a)[0];
            ai2 = (long double) ((double *) a)[1];
            br2 = (long double) ((double *) b)[0];
            bi2 = (long double) ((double *) b)[1];
            return ( ar2*ar2+ai2*ai2 < br2*br2+bi2*bi2 );
        default:
            return 1;
    }

    return 0;
}

const char * value_to_string (enum ADIOS_DATATYPES type, void * data)
{
    static char s [100];
    s [0] = 0;


    switch (type)
    {
        case adios_unsigned_byte:
            sprintf (s, "%u", *(((uint8_t *) data)));
            break;

        case adios_byte:
            sprintf (s, "%d", *(((int8_t *) data)));
            break;

        case adios_short:
            sprintf (s, "%hd", *(((int16_t *) data)));
            break;

        case adios_unsigned_short:
            sprintf (s, "%uh", *(((uint16_t *) data)));
            break;

        case adios_integer:
            sprintf (s, "%d", *(((int32_t *) data)));
            break;

        case adios_unsigned_integer:
            sprintf (s, "%u", *(((uint32_t *) data)));
            break;

        case adios_long:
            sprintf (s, "%lld", *(((int64_t *) data)));
            break;

        case adios_unsigned_long:
            sprintf (s, "%llu", *(((uint64_t *) data)));
            break;

        case adios_real:
            sprintf (s, "%f", *(((float *) data)));
            break;

        case adios_double:
            sprintf (s, "%le", *(((double *) data)));
            break;

        case adios_long_double:
            sprintf (s, "%Le", *(((long double *) data)));
            break;

        case adios_string:
            sprintf (s, "%s", ((char *) data));
            break;

        case adios_complex:
            sprintf (s, "(%f %f)", *(((float *) data) + 0)
                                 , *(((float *) data) + 1)
                    );
            break;

        case adios_double_complex:
            sprintf (s, "(%lf %lf)", *(((double *) data) + 0)
                                   , *(((double *) data) + 1)
                    );
            break;
    }

    return s;
}
void print_var_header (struct adios_var_header_struct_v1 * var_header)
{
    int i = 0;
    printf ("\t\tVar Name (ID): %s (%d)\n", var_header->name, var_header->id);
    printf ("\t\tVar Path: %s\n", var_header->path);
    printf ("\t\tIs Dimension: %c\n"
           ,(var_header->is_dim == adios_flag_yes ? 'Y' : 'N')
           );
}

void copy_data (void *dst, void *src,
        int idim,
        int ndim,
        uint64_t* size_in_dset, 
        uint64_t* ldims, 
        const uint64_t * readsize, 
        uint64_t dst_stride, 
        uint64_t src_stride,
        uint64_t dst_offset, 
        uint64_t src_offset,
        uint64_t ele_num,
        int      size_of_type
        )
{
    unsigned int i, j;
    uint64_t dst_offset_new=0; 
    uint64_t src_offset_new=0;
    uint64_t src_step, dst_step;
    if (ndim-1==idim) {
        for (i=0;i<size_in_dset[idim];i++) {
            memcpy (dst + (i*dst_stride+dst_offset)*size_of_type,
                    src + (i*src_stride+src_offset)*size_of_type,
                    ele_num*size_of_type);
        }
        return;
    }

    for (i = 0; i<size_in_dset[idim];i++) {
        // get the different step granularity 
        // for each different reading pattern broke
        src_step = 1;
        dst_step = 1;
        for (j = idim+1; j <= ndim-1;j++) {
            src_step *= ldims[j];
            dst_step *= readsize[j];
        }
        src_offset_new =src_offset + i * src_stride * src_step;
        dst_offset_new = dst_offset + i * dst_stride * dst_step;
        copy_data ( dst, src, idim+1, ndim, size_in_dset,
                ldims,readsize, 
                dst_stride, src_stride,
                dst_offset_new, src_offset_new,
                ele_num, size_of_type);
    }
}

uint64_t bp_get_type_size (enum ADIOS_DATATYPES type, void * var)
{
    switch (type)
    {
        case adios_byte:
        case adios_unsigned_byte:
            return 1;

        case adios_string:
            if (!var)
                return 0;
            else
                return strlen ((char *) var);

        case adios_short:
        case adios_unsigned_short:
            return 2;

        case adios_integer:
        case adios_unsigned_integer:
            return 4;

        case adios_long:
        case adios_unsigned_long:
            return 8;

        case adios_real:
            return 4;

        case adios_double:
            return 8;

        case adios_long_double:
            return 16;

        case adios_complex:
            return 2 * 4;

        case adios_double_complex:
            return 2 * 8;

        default:
            return -1;
    }
}
void alloc_namelist (char ***namelist, int length)
{
        int j;

        *namelist = (char **) malloc(length*sizeof(char*));
        for (j=0;j<length;j++)
                (*namelist)[j] = (char *) malloc(255);

        return;
}

void free_namelist (char **namelist, int length)
{
        int i;
        if (namelist) {
                for (i=0;i<length;i++) {
                        if(namelist[i])
                                free(namelist[i]);
                }
                free(namelist);
        }
        return;
}
