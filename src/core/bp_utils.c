/* 
 * ADIOS is freely available under the terms of the BSD license described
 * in the COPYING file in the top level directory of this source distribution.
 *
 * Copyright (c) 2008 - 2009.  UT-BATTELLE, LLC. All rights reserved.
 */

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <stdarg.h>
#include <sys/types.h>
#include <string.h>
#include "public/adios.h"
#include "public/adios_read.h"
#include "public/adios_error.h"
#include "core/bp_utils.h"
#include "core/adios_bp_v1.h"
#include "core/adios_endianness.h"
#define BYTE_ALIGN 8
#define MINIFOOTER_SIZE 28


#ifdef DMALLOC
#include "dmalloc.h"
#endif

#define BUFREAD8(b,var)  var = (uint8_t) *(b->buff + b->offset); \
                         b->offset += 1;

#define BUFREAD16(b,var) var = *(uint16_t *) (b->buff + b->offset); \
                         if (b->change_endianness == adios_flag_yes) \
                             swap_16(var); \
                         b->offset += 2;

#define BUFREAD32(b,var) var = *(uint32_t *) (b->buff + b->offset); \
                         if (b->change_endianness == adios_flag_yes) \
                             swap_32(var); \
                         b->offset += 4;

#define BUFREAD64(b,var) var = *(uint64_t *) (b->buff + b->offset); \
                         if (b->change_endianness == adios_flag_yes) \
                             swap_64(var); \
                         b->offset += 8;

/* prototypes */
void * bp_read_data_from_buffer(struct adios_bp_buffer_struct_v1 *b, enum ADIOS_DATATYPES type);
int bp_parse_characteristics (struct adios_bp_buffer_struct_v1 * b, struct adios_index_var_struct_v1 ** root, uint64_t j);



void bp_alloc_aligned (struct adios_bp_buffer_struct_v1 * b, uint64_t size)
{
        
    b->allocated_buff_ptr =  malloc (size + BYTE_ALIGN - 1);
    if (!b->allocated_buff_ptr)
    {
        adios_error ( err_no_memory, "Cannot allocate %llu bytes\n", size);

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
        adios_error ( err_no_memory, "Cannot allocate %llu bytes\n", size);

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
        adios_error (err_file_open_error, "MPI open failed for %s: '%s'\n", filename, e);
        return adios_flag_no;
    }

    MPI_File_get_size (fh->mpi_fh, &file_size);
    fh->b->file_size = file_size;
    fh->mfooter.file_size = file_size;
     
    return 0;
}

/* This routine does the parallel bp file open and index parsing.
 */
int bp_open (const char * fname,
             MPI_Comm comm,
             struct BP_FILE * fh)
{
    int rank;
    uint64_t header_size;

    MPI_Comm_rank (comm, &rank);

    adios_buffer_struct_init (fh->b);

    if (bp_read_open (fname, comm, fh))
    {
        return -1;
    }

    /* Only rank 0 reads the footer and it broadcasts to all other processes */
    if (rank == 0)
    {
        if (bp_read_minifooter (fh))
        {
            return -1;
        }
    }

    /* Broadcast to all other processors */
    MPI_Bcast (&fh->mfooter, sizeof (struct bp_minifooter), MPI_BYTE, 0, comm);

    header_size = fh->mfooter.file_size-fh->mfooter.pgs_index_offset;

    if (rank != 0)
    {
        if (!fh->b->buff)
        {
            bp_alloc_aligned (fh->b, header_size);
            assert (fh->b->buff);

            memset (fh->b->buff, 0, header_size);
            fh->b->offset = 0;
        }
    }

    MPI_Barrier (comm);
    MPI_Bcast (fh->b->buff, fh->mfooter.file_size-fh->mfooter.pgs_index_offset, MPI_BYTE, 0, comm);

    /* Everyone parses the index on its own */
    bp_parse_pgs (fh);
    bp_parse_vars (fh);
    bp_parse_attrs (fh);

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
            adios_error (err_no_memory, "could not allocate %d bytes\n", MINIFOOTER_SIZE);
            return 1;
        }
        memset (b->buff, 0, MINIFOOTER_SIZE);
        b->offset = 0;
    }
    MPI_File_seek (bp_struct->mpi_fh, (MPI_Offset) attrs_end, MPI_SEEK_SET);
    MPI_File_read (bp_struct->mpi_fh, b->buff, MINIFOOTER_SIZE, MPI_BYTE, &status);
    
    /*memset (&mh->pgs_index_offset, 0, MINIFOOTER_SIZE);
    memcpy (&mh->pgs_index_offset, b->buff, MINIFOOTER_SIZE);*/
    
    /* get version id. Needs the bp->offset be pointing to the last 4 bytes of the buffer,
       It also sets b->change_endianness */
    /* Note that b is not sent over to processes, only the minifooter and then b->buff (the footer) */
    b->offset = MINIFOOTER_SIZE - 4;
    adios_parse_version (b, &mh->version);
    mh->change_endianness = b->change_endianness;
    b->offset = 0; // reset offset to beginning 

    BUFREAD64(b, b->pg_index_offset)
    mh->pgs_index_offset = b->pg_index_offset;

    BUFREAD64(b, b->vars_index_offset)
    mh->vars_index_offset = b->vars_index_offset;

    BUFREAD64(b, b->attrs_index_offset)
    mh->attrs_index_offset = b->attrs_index_offset;

    b->end_of_pgs = b->pg_index_offset;
    b->pg_size = b->vars_index_offset - b->pg_index_offset;
    b->vars_size = b->attrs_index_offset - b->vars_index_offset;
    b->attrs_size = attrs_end - b->attrs_index_offset;

    /* Read the whole footer */
    /* FIXME: including the last 28 bytes read already above and it seems that is not processed anymore */
    /* It will be sent to all processes */
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

/****************/
/* Parse GROUPS */
/****************/
int bp_parse_pgs (struct BP_FILE * fh)
{
    struct bp_index_pg_struct_v1 ** root = &(fh->pgs_root); // need the pointer to it to malloc below
    struct adios_bp_buffer_struct_v1 * b = fh->b;
    struct bp_minifooter * mh = &(fh->mfooter);
    uint64_t i;

    /* Note that at this point, many variables of b->* is unset (init'ed to 0).
       It's the minifooter which holds accurate information.
       b holds the footer data from the file */

    b->offset = 0;
    b->change_endianness = (enum ADIOS_FLAG) mh->change_endianness;

    BUFREAD64(b, mh->pgs_count)
    BUFREAD64(b, mh->pgs_length)

    int j;
    uint64_t group_count = 0;
    char ** namelist;
    char fortran_flag;

    namelist = (char **) malloc(sizeof(char *)*mh->pgs_count);
    uint16_t * grpidlist = (uint16_t *) malloc(sizeof(uint16_t)*mh->pgs_count);

    uint32_t tidx_start, tidx_stop; /* Determine first and last timestep in file */

    for (i = 0; i < mh->pgs_count; i++) {
        uint16_t length_of_group;
        namelist[i] = 0;    
        // validate remaining length
        BUFREAD16(b, length_of_group)

        if (!*root)
        {
            *root = (struct bp_index_pg_struct_v1 *)
                malloc (sizeof(struct bp_index_pg_struct_v1));
            memset (*root, 0, sizeof(struct bp_index_pg_struct_v1));
            (*root)->next = 0;
        }
        uint16_t length_of_name;

        BUFREAD16(b, length_of_name)
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
            
        BUFREAD8(b, fortran_flag)
        (*root)->adios_host_language_fortran =
            (fortran_flag == 'y' ? adios_flag_yes : adios_flag_no);

        BUFREAD32(b, (*root)->process_id)

        BUFREAD16(b, length_of_name)
        (*root)->time_index_name = (char *) malloc (length_of_name + 1);
        (*root)->time_index_name [length_of_name] = '\0';
        memcpy ((*root)->time_index_name, b->buff + b->offset, length_of_name);
        b->offset += length_of_name;

        BUFREAD32(b, (*root)->time_index)

        BUFREAD64(b, (*root)->offset_in_file)

        if (i == 0)
            tidx_start = (*root)->time_index;
        if (i == mh->pgs_count-1) {
            tidx_stop = (*root)->time_index;
            mh->time_steps = tidx_stop - tidx_start + 1;
        }

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
    uint64_t grpid = grpidlist[0];
    uint32_t pg_time_count = 0, first_pg;
    uint32_t time_id = tidx_start;
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
    return 0;
}


/********************/
/* Parse ATTRIBUTES */
/********************/
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

    BUFREAD16(b, mh->attrs_count)
    BUFREAD64(b, mh->attrs_length)

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
        int type_size;

        BUFREAD32(b, attr_entry_length)
        BUFREAD16(b, (*root)->id)

        BUFREAD16(b, len)
        (*root)->group_name = (char *) malloc (len + 1);
        (*root)->group_name [len] = '\0';
        strncpy ((*root)->group_name, b->buff + b->offset, len);
        b->offset += len;

        BUFREAD16(b, len)
        (*root)->attr_name = (char *) malloc (len + 1);
        (*root)->attr_name [len] = '\0';
        strncpy ((*root)->attr_name, b->buff + b->offset, len);
        b->offset += len;

        BUFREAD16(b, len)
        (*root)->attr_path = (char *) malloc (len + 1);
        (*root)->attr_path [len] = '\0';
        strncpy ((*root)->attr_path, b->buff + b->offset, len);
        b->offset += len;

        BUFREAD8(b, flag)
        (*root)->type = (enum ADIOS_DATATYPES) flag;

        type_size = bp_get_type_size ((*root)->type, "");
        
        if (type_size == -1)
        {
            //type_size = 4;
            (*root)->type = adios_unknown;
        }
        

        BUFREAD64(b, characteristics_sets_count)
        (*root)->characteristics_count = characteristics_sets_count;
        (*root)->characteristics_allocated = characteristics_sets_count;

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

            BUFREAD8(b, characteristic_set_count)
            BUFREAD32(b, characteristic_set_length)

            while (item < characteristic_set_count)
            {
                uint8_t flag;
                enum ADIOS_CHARACTERISTICS c;

                BUFREAD8(b, flag)
                c = (enum ADIOS_CHARACTERISTICS) flag;

                switch (c)
                {
                    case adios_characteristic_value:
                        (*root)->characteristics [j].value = bp_read_data_from_buffer(b, (*root)->type);
                        break;

                    case adios_characteristic_offset:
                        BUFREAD64(b, (*root)->characteristics [j].offset)
                        break;

                    case adios_characteristic_payload_offset:
                        BUFREAD64(b, (*root)->characteristics [j].payload_offset)
                        break;

                    case adios_characteristic_file_index:
                        BUFREAD32(b, (*root)->characteristics [j].file_index);
                        break;

                    case adios_characteristic_time_index:
                        BUFREAD32(b, (*root)->characteristics [j].time_index)
                        break;
                    case adios_characteristic_var_id:
                        BUFREAD16(b, (*root)->characteristics [j].var_id)
                        break;
                }
                item++;
            }
            /* Old BP files do not have time_index characteristics, so we
               set it here automatically: j div # of pgs per timestep
               Assumed that in old BP files, all pgs write each variable in each timestep.*/
            if ((*root)->characteristics [j].time_index == 0) {
                (*root)->characteristics [j].time_index = 
                     j / (mh->pgs_count / (fh->tidx_stop - fh->tidx_start + 1)) + 1;
                /*printf("OldBP: attr %s time_index set to %d\n", 
                        (*root)->attr_name,
                        (*root)->characteristics [j].time_index);*/
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
    memset (attr_offsets, 0, mh->attrs_count * sizeof(uint64_t *));

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

    return 0;
}

/*******************/
/* Parse VARIABLES */
/*******************/
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

    BUFREAD16(b, mh->vars_count)
    BUFREAD64(b, mh->vars_length)

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
        int type_size;

        BUFREAD32(b, var_entry_length)
        BUFREAD16(b, (*root)->id)

        BUFREAD16(b, len)
        (*root)->group_name = (char *) malloc (len + 1);
        (*root)->group_name [len] = '\0';
        strncpy ((*root)->group_name, b->buff + b->offset, len);
        b->offset += len;

        BUFREAD16(b, len)
        (*root)->var_name = (char *) malloc (len + 1);
        (*root)->var_name [len] = '\0';
        strncpy ((*root)->var_name, b->buff + b->offset, len);
        b->offset += len;

        BUFREAD16(b, len)
        (*root)->var_path = (char *) malloc (len + 1);
        (*root)->var_path [len] = '\0';
        strncpy ((*root)->var_path, b->buff + b->offset, len);
        b->offset += len;

        BUFREAD8(b, flag)
        (*root)->type = (enum ADIOS_DATATYPES) flag;
        type_size = bp_get_type_size ((*root)->type, "");

        BUFREAD64(b, characteristics_sets_count)
        (*root)->characteristics_count = characteristics_sets_count;
        (*root)->characteristics_allocated = characteristics_sets_count;

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

            BUFREAD8(b, characteristic_set_count)
            BUFREAD32(b, characteristic_set_length)
                
            while (item < characteristic_set_count) {
                bp_parse_characteristics (b, root, j);
                item++;
            }

            /* Old BP files do not have time_index characteristics, so we
               set it here automatically: j div # of pgs per timestep
               Assumed that in old BP files, all pgs write each variable in each timestep.*/
            if ((*root)->characteristics [j].time_index == 0) {
                (*root)->characteristics [j].time_index = 
                     j / (mh->pgs_count / (fh->tidx_stop - fh->tidx_start + 1)) + 1;
                /*printf("OldBP: var %s time_index set to %d\n",
                        (*root)->var_name, 
                        (*root)->characteristics [j].time_index);*/
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
    memset ( var_offsets, 0, mh->vars_count*sizeof(uint64_t *));

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
}

int bp_parse_characteristics (struct adios_bp_buffer_struct_v1 * b,
                    struct adios_index_var_struct_v1 ** root,
                  uint64_t j)
{
    uint8_t flag;
    enum ADIOS_CHARACTERISTICS c;

    BUFREAD8(b, flag)
    c = (enum ADIOS_CHARACTERISTICS) flag;

    switch (c) {

        case adios_characteristic_value:
            (*root)->characteristics [j].value = bp_read_data_from_buffer(b, (*root)->type);
            break;

		// NCSU - Adding in backward compatibility
		case adios_characteristic_max:
		{
			if (!((*root)->characteristics [j].stats))
			{
				(*root)->characteristics [j].stats = malloc (sizeof(struct adios_index_characteristics_stat_struct *));
				(*root)->characteristics [j].stats[0] = malloc (2 * sizeof(struct adios_index_characteristics_stat_struct));
				(*root)->characteristics [j].bitmap = 0;
			}
			(*root)->characteristics [j].bitmap |= (1 << adios_statistic_max);
			(*root)->characteristics [j].stats[0][adios_statistic_max].data = bp_read_data_from_buffer(b, (*root)->type);
			break;
		}

        // NCSU - Adding in backward compatibility
        case adios_characteristic_min:
        {
            if (!((*root)->characteristics [j].stats))
            {
                (*root)->characteristics [j].stats = malloc (sizeof(struct adios_index_characteristics_stat_struct *));
                (*root)->characteristics [j].stats[0] = malloc (2 * sizeof(struct adios_index_characteristics_stat_struct));
                (*root)->characteristics [j].bitmap = 0;
            }
            (*root)->characteristics [j].bitmap |= (1 << adios_statistic_min);
            (*root)->characteristics [j].stats[0][adios_statistic_min].data = bp_read_data_from_buffer(b, (*root)->type);
			break;
        }

		// NCSU - Parse the statistical information based in the bitmap
        case adios_characteristic_stat:
		{
			uint8_t i, c, idx; 
			uint8_t count = adios_get_stat_set_count ((*root)->type);
			uint16_t characteristic_size;

			(*root)->characteristics [j].stats = malloc (count * sizeof(struct adios_index_characteristics_stat_struct *));

			for (c = 0; c < count; c ++)
			{
				i = idx = 0;
				(*root)->characteristics [j].stats[c] = malloc (ADIOS_STAT_LENGTH * sizeof(struct adios_index_characteristics_stat_struct));

				while ((*root)->characteristics[j].bitmap >> i)
				{
					(*root)->characteristics [j].stats[c][i].data = 0;
					if (((*root)->characteristics[j].bitmap >> i) & 1)
					{
						if (i == adios_statistic_hist)
						{
							uint32_t bi;
							
							(*root)->characteristics [j].stats[c][idx].data = malloc (sizeof(struct adios_index_characteristics_hist_struct));
							struct adios_index_characteristics_hist_struct * hist = (*root)->characteristics [j].stats[c][idx].data; 

            				BUFREAD32(b, hist->num_breaks)
            				hist->min = * (double *) bp_read_data_from_buffer(b, adios_double);
            				hist->max = * (double *) bp_read_data_from_buffer(b, adios_double);

            				hist->frequencies = malloc((hist->num_breaks + 1) * adios_get_type_size(adios_unsigned_integer, ""));
            				for (bi = 0; bi <= hist->num_breaks; bi ++) {
            				    BUFREAD32(b, hist->frequencies[bi])
            				}

            				hist->breaks = malloc(hist->num_breaks * adios_get_type_size(adios_double, ""));
            				for (bi = 0; bi < hist->num_breaks; bi ++) {
            				    hist->breaks[bi] = * (double *) bp_read_data_from_buffer(b, adios_double);
            				}
						}
						else
						{
							characteristic_size = adios_get_stat_size((*root)->characteristics [j].stats[c][idx].data
                                                                                                 ,(*root)->type
                                                                                                 ,(enum ADIOS_STAT)i
                                                                                                 );
							(*root)->characteristics [j].stats[c][idx].data = malloc (characteristic_size);

							void * data = (*root)->characteristics [j].stats[c][idx].data;
			            	memcpy (data, (b->buff + b->offset), characteristic_size);
            				b->offset += characteristic_size;

            				if(b->change_endianness == adios_flag_yes) 
                        		swap_ptr(data, characteristic_size * 8);
						}
						idx ++;
					}
					i ++;
				}	
			}	
			break;
		}

		// NCSU - Statistics. Read the bitmap
		case adios_characteristic_bitmap:
			BUFREAD32(b, (*root)->characteristics [j].bitmap);
			break;

        case adios_characteristic_offset: 
            BUFREAD64(b, (*root)->characteristics [j].offset)
            break;

        case adios_characteristic_payload_offset: 
            BUFREAD64(b, (*root)->characteristics [j].payload_offset)
            break;

        case adios_characteristic_file_index:
            BUFREAD32(b, (*root)->characteristics [j].file_index)
            break;

        case adios_characteristic_time_index:
            BUFREAD32(b, (*root)->characteristics [j].time_index)
            break;

        case adios_characteristic_dimensions:
        {
            uint16_t dims_length, di, dims_num;
            BUFREAD8(b, (*root)->characteristics [j].dims.count)
            BUFREAD16(b, dims_length);

            (*root)->characteristics [j].dims.dims = (uint64_t *) malloc (dims_length);

            dims_num = dims_length / 8;
            for (di = 0; di < dims_num; di ++) {
                BUFREAD64(b, ((*root)->characteristics [j].dims.dims)[di]);
            }
            break;
        }

        default:
            fprintf (stderr, "Unknown characteristic:%d. skipped.\n", c);
            break;
    }
    return 0;
}

/* Seek to the specified step and prepare related fields in ADIOS_FILE structure */
int bp_seek_to_step (ADIOS_FILE * fp, int tostep)
{
    int i, j;
    struct BP_PROC * p = (struct BP_PROC *) fp->fh;
    struct BP_FILE * fh = p->fh;
    struct adios_index_var_struct_v1 * var_root;
    struct adios_index_attribute_struct_v1 * attr_root;

    if (tostep < fh->tidx_start || tostep > fh->tidx_stop)
    {
        adios_error (err_invalid_timestep, "Invalid step: %d\n", tostep);
        return -1;
    }

    /* Prepare vars */
    fp->nvars = 0;
    var_root = fh->vars_root;

    while (var_root)
    {
        for (i = 0; i < var_root->characteristics_count; i++)
        {
            if (var_root->characteristics[i].time_index == tostep)
            {
                fp->nvars++;
                break;    
            }
        }
      
        var_root = var_root->next;
    }

    alloc_namelist (&fp->var_namelist, fp->nvars);

    var_root = fh->vars_root;
    j = 0;
    while (var_root)
    {
        for (i = 0; i < var_root->characteristics_count; i++)
        {
            if (var_root->characteristics[i].time_index == tostep)
            {
                if (strcmp (var_root->var_path,"/"))
                {
                    fp->var_namelist[j] = (char *)malloc (strlen((var_root)->var_name)
                                          + strlen (var_root->var_path) + 1 + 1   // extra / and ending \0
                                          );
                    strcpy (fp->var_namelist[j], var_root->var_path);
                }
                else
                {
                    fp->var_namelist[j] = (char *) malloc (strlen (var_root->var_name) + 1 + 1);
                    fp->var_namelist[j][0] = '\0';
                }

                strcat (fp->var_namelist[j], "/");
                strcat (fp->var_namelist[j], var_root->var_name);

                j++;

                break;
            }
        }

        var_root = var_root->next;
    }

    /* Prepare attrs */
    fp->nattrs = 0;
    attr_root = fh->attrs_root;

    while (attr_root)
    {
        for (i = 0; i < attr_root->characteristics_count; i++)
        {
            if (attr_root->characteristics[i].time_index == tostep)
            {
                fp->nattrs++;
                break;
            }
        }

        attr_root = attr_root->next;
    }

    alloc_namelist (&fp->attr_namelist, fp->nattrs);

    attr_root = fh->attrs_root;
    j = 0;
    while (attr_root)
    {
        for (i = 0; i < attr_root->characteristics_count; i++)
        {
            if (attr_root->characteristics[i].time_index == tostep)
            {
                if (strcmp (attr_root->attr_path,"/"))
                {
                    fp->attr_namelist[j] = (char *)malloc (strlen((attr_root)->attr_name)
                                          + strlen (attr_root->attr_path) + 1 + 1   // extra / and ending \0
                                          );
                    strcpy (fp->attr_namelist[j], attr_root->attr_path);
                }
                else
                {
                    fp->attr_namelist[j] = (char *) malloc (strlen (attr_root->attr_name) + 1 + 1);
                    fp->attr_namelist[j][0] = '\0';
                }

                strcat (fp->attr_namelist[j], "/");
                strcat (fp->attr_namelist[j], attr_root->attr_name);

                j++;

                break;
            }
        }

        attr_root = attr_root->next;
    }

    fp->current_step = tostep;

    return 0;
}

/* get local and global dimensions and offsets from a variable characteristics 
   return: 1 = it is a global array, 0 = local array
*/
int bp_get_dimension_characteristics(struct adios_index_characteristic_struct_v1 *ch,
                                    uint64_t *ldims, uint64_t *gdims, uint64_t *offsets)
{
    int is_global = 0; // global array or just an array written by one process?
    int ndim = ch->dims.count;
    int k;

    for (k = 0; k < ndim; k++)
    {
        ldims[k] = ch->dims.dims[k * 3];
        gdims[k] = ch->dims.dims[k * 3 + 1];
        offsets[k] = ch->dims.dims[k * 3 + 2];
        is_global = is_global || gdims[k];
    }

    return is_global;
}

/* Fill out ndim and dims for the variable.
   ndim and dims shouldn't include 'time' dimension.
*/
void bp_get_dimensions (struct adios_index_var_struct_v1 * var_root, int file_is_fortran,
                        int * ndim, uint64_t ** dims)
{
    int i, j;
    int is_global; // global array or just an array written by one process?
    uint64_t ldims[32];
    uint64_t gdims[32];
    uint64_t offsets[32];

    /* Get dimension information */
    * ndim = var_root->characteristics [0].dims.count;
    * dims = 0;
    if (* ndim == 0)
    {
        /* done with this scalar variable */
        return ;
    }

    * dims = (uint64_t *) malloc (sizeof (uint64_t) * (* ndim));
    assert (* dims);

    memset (*dims, 0, sizeof (uint64_t) * (* ndim));

    is_global = bp_get_dimension_characteristics (&(var_root->characteristics[0]),
                                                  ldims, gdims, offsets);

    if (!is_global)
    {
        j = 0;
        /* local array */
        for (i = 0; i < * ndim; i++)
        {
            /* size of time dimension is always registered as 1 for an array */
            if (ldims[i] == 1 && var_root->characteristics_count > 1)
            {
            }
            else
            {
                (* dims)[j++] = ldims[i];
            }
        }
    }
    else
    {
        /* global array:
           time dimension: ldims=1, gdims=0
           in C array, it can only be the first dim
           in Fortran array, it can only be the last dim
           (always the slowest changing dim)
        */
        if (gdims[* ndim - 1] == 0) // with time
        {
            if (!file_is_fortran)
            {
                /* first dimension is the time (C array)
                 * ldims[0] = 1 but gdims does not contain time info and 
                 * gdims[0] is 1st data dimension and 
                 * gdims is shorter by one value than ldims in case of C.
                 * Therefore, gdims[*ndim-1] = 0 if there is a time dimension. 
                 */
                // error check
                if (* ndim > 1 && ldims[0] != 1)
                {
                    fprintf(stderr,"ADIOS Error: this is a BP file with C ordering but we didn't find"
                            "an array to have time dimension in the first dimension. l:g:o = (");
                    for (i = 0; i < * ndim; i++)
                    {
                        fprintf (stderr,"%llu:%llu:%llu%s", ldims[i], gdims[i], offsets[i], (i<*ndim-1 ? ", " : "") );
                    }

                    fprintf(stderr, ")\n");
                }
            }
            else
            {
                // last dimension is the time (Fortran array)
                if (* ndim > 1 && ldims[* ndim - 1] != 1)
                {
                    fprintf(stderr,"ADIOS Error: this is a BP file with Fortran array ordering but we didn't find"
                            "an array to have time dimension in the last dimension. l:g:o = (");
                    for (i = 0; i < * ndim; i++)
                    {
                        fprintf (stderr,"%llu:%llu:%llu%s", ldims[i], gdims[i], offsets[i], (i<*ndim-1 ? ", " : "") );
                    }

                    fprintf (stderr, ")\n");
                }
            }

            * ndim = * ndim - 1;
        }

        for (i = 0; i < * ndim; i++)
        {
            (* dims)[i] = gdims[i];
        }

    }
}

/* Get dimensions of a variable and flip them if swap_flag is set */
void bp_get_and_swap_dimensions (struct adios_index_var_struct_v1 *var_root, int file_is_fortran,
                                 int *ndim, uint64_t **dims, int swap_flag)
{
    int timedim = 0;

    bp_get_dimensions (var_root, file_is_fortran, ndim, dims);
    if (swap_flag)
    {
        /* dummy timedim */
        swap_order (*ndim, *dims, &timedim);
    }
}

/* Get # of steps of a variable */
int get_var_nsteps (struct adios_index_var_struct_v1 * var_root)
{
    uint64_t i;
    int nsteps = 0;
    int prev_step = -1;

    for (i = 0; i < var_root->characteristics_count; i++)
    {
        if (var_root->characteristics[i].time_index != prev_step)
        {
            prev_step = var_root->characteristics[i].time_index;
            nsteps ++;
        }
    }

    return nsteps;
}

/* Get # of blocks per step */
int * get_var_nblocks (struct adios_index_var_struct_v1 * var_root, int nsteps)
{
    int i, j, prev_step = -1;
    int * nblocks = (int *) malloc (sizeof (int) * nsteps);

    assert (nblocks); 

    memset (nblocks, 0, sizeof (int) * nsteps);

    j = -1;
    for (i = 0; i < var_root->characteristics_count; i++)
    {
        if (var_root->characteristics[i].time_index != prev_step)
        {
            j ++;
            prev_step = var_root->characteristics[i].time_index;
        }

        nblocks[j] ++;
    }

    return nblocks;
}

void * bp_read_data_from_buffer(struct adios_bp_buffer_struct_v1 *b, enum ADIOS_DATATYPES type)
{
    int16_t data_size;
    void * data = 0;

    if (type == adios_string) {
        BUFREAD16(b, data_size)
        data = malloc (data_size + 1);
    } else {
        data_size = bp_get_type_size (type, "");
        data = malloc (data_size);
    }

    if (!data) {
        fprintf (stderr, "bp_read_data_from_buffer: cannot allocate %d bytes\n",data_size);
        return 0;
    }

    switch (type)
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
            memcpy (data, (b->buff + b->offset), data_size);
            b->offset += data_size;
            if(b->change_endianness == adios_flag_yes) {
                switch (data_size) {
                    case 2:
                        swap_16_ptr(data);
                        break;
                    case 4:
                        swap_32_ptr(data);
                        break;
                    case 8:
                        swap_64_ptr(data);
                        break;
                    case 16:
                        swap_128_ptr(data);
                        break;
                }
            }
            break;

        case adios_complex:
            memcpy (data, (b->buff + b->offset), data_size);
            if(b->change_endianness == adios_flag_yes) {
                swap_32_ptr(data); // swap REAL part 4 bytes
                swap_32_ptr( ((char *)data) + 4); // swap IMG part 4 bytes
            }
            b->offset += data_size;
            break;

        case adios_double_complex:
            memcpy (data, (b->buff + b->offset), data_size);
            if(b->change_endianness == adios_flag_yes) {
                swap_64_ptr(data); // swap REAL part 8 bytes
                swap_64_ptr( ((char *)data) + 8); // swap IMG part 8 bytes
            }
            b->offset += data_size;
            break;

        case adios_string:
            memcpy (data, (b->buff + b->offset), data_size);
            b->offset += data_size;
            ((char *) data) [data_size] = '\0';
            break;

        default:
            free (data);
            data = 0;
            break;
    }
    return data;
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
/*
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
*/
/*
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
*/
/*
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
*/

/* Return 1 if a < b wrt the given type, otherwise 0 */
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

const char * bp_value_to_string (enum ADIOS_DATATYPES type, void * data)
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

/*
void print_var_header (struct adios_var_header_struct_v1 * var_header)
{
    int i = 0;
    printf ("\t\tVar Name (ID): %s (%d)\n", var_header->name, var_header->id);
    printf ("\t\tVar Path: %s\n", var_header->path);
    printf ("\t\tIs Dimension: %c\n"
           ,(var_header->is_dim == adios_flag_yes ? 'Y' : 'N')
           );
}
*/

int is_fortran_file (struct BP_FILE * fh)
{
    return (fh->pgs_root->adios_host_language_fortran == adios_flag_yes);
}

int has_subfiles (struct BP_FILE * fh)
{
    return (fh->mfooter.version & ADIOS_VERSION_HAVE_SUBFILE);
}

/****************************************************
  Find the var associated with the given variable id 
*****************************************************/
struct adios_index_var_struct_v1 * bp_find_var_byid (struct BP_FILE * fh, int varid)
{
    struct adios_index_var_struct_v1 * var_root = fh->vars_root;
    int i;

    for (i = 0; i < varid && var_root; i++)
    {
        var_root = var_root->next;
    }

    if (i != varid)
    {
        adios_error (err_corrupted_variable,
               "Variable id=%d is valid but was not found in internal data structures!",
               varid);
        return NULL;
    }

    return var_root;
}

/* Check whether an array is global */
int is_global_array (struct adios_index_characteristic_struct_v1 *ch)
{
    int is_global = 0; // global array or just an array written by one process?
    int ndim = ch->dims.count, k;

    for (k = 0; k < ndim; k ++)
    {
        is_global = is_global || ch->dims.dims[k*3 + 1];
    }

    return is_global;
}

/** Return the memory size of one data element of an adios type.
 *  If the type is adios_string, and the second argument is
 *  the string itself, it returns strlen(var)+1.
 *  For other types, it does not care about var and returns
 *  the size occupied by one element. 
 *
 *  Note that adios_internals:adios_get_type_size returns
 *  strlen(var) for strings. 
 */ 
int bp_get_type_size (enum ADIOS_DATATYPES type, void * var)
{
    switch (type)
    {
        case adios_byte:
        case adios_unsigned_byte:
            return 1;

        case adios_string:
            if (!var)
                return 1;
            else
                return strlen ((char *) var) + 1;

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

double bp_value_to_double (enum ADIOS_DATATYPES type, void * data)
{
	switch (type)
    {
        case adios_string:
            return 0;

        case adios_complex:
			return * ((float *) data);
		
        case adios_double_complex:
			return * ((double *) data);

        case adios_double:
			return * ((double *) data);

        case adios_long_double:
			return * ((long double *) data);

		case adios_unsigned_byte:
            return * ((uint8_t *) data);

        case adios_byte:
            return * ((int8_t *) data);

        case adios_short:
            return * ((int16_t *) data);

        case adios_unsigned_short:
            return * ((uint16_t *) data);

        case adios_integer:
            return * ((int32_t *) data);

        case adios_unsigned_integer:
            return * ((uint32_t *) data);

        case adios_long:
            return * ((int64_t *) data);

        case adios_unsigned_long:
            return * ((uint64_t *) data);

        case adios_real:
            return * ((float *) data);
    }
}
