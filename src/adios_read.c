#include <stdlib.h>
#include <string.h>
#include "adios.h"
#include "bp_utils.h"
#include "adios_read.h"
#define BYTE_ALIGN 8

static void alloc_aligned (struct adios_bp_buffer_struct_v1 * b, uint64_t size)
{
        
    b->allocated_buff_ptr =  malloc (size + BYTE_ALIGN - 1);
    if (!b->allocated_buff_ptr)
    {
        fprintf (stderr, "Cannot allocate: %llu\n", size);

        b->buff = 0;
        b->length = 0;

        return;
    }
    uint64_t p = (uint64_t) b->allocated_buff_ptr;
    b->buff = (char *) ((p + BYTE_ALIGN - 1) & ~(BYTE_ALIGN - 1));
    b->length = size;
}

static void realloc_aligned (struct adios_bp_buffer_struct_v1 * b
                            ,uint64_t size
                            )
{
    b->allocated_buff_ptr = realloc (b->allocated_buff_ptr
                                    ,size + BYTE_ALIGN - 1
                                    );
    if (!b->allocated_buff_ptr)
    {
        fprintf (stderr, "Cannot allocate: %llu\n", size);

        b->buff = 0;
        b->length = 0;

        return;
    }
    uint64_t p = (uint64_t) b->allocated_buff_ptr;
    b->buff = (char *) ((p + BYTE_ALIGN - 1) & ~(BYTE_ALIGN - 1));
    b->length = size;
}

int adios_fopen ( int64_t * fh_p,
                  const char * fname,
                  MPI_Comm comm
                )
{
    int rc, rank;    
    struct BP_FILE * fh = (struct BP_FILE *)
        malloc (sizeof (struct BP_FILE));
    fh->comm =  comm;
    fh->gh = 0;
    fh->pgs_root = 0;
    fh->vars_root = 0;
    fh->attrs_root = 0;
    fh->b = malloc (sizeof (struct adios_bp_buffer_struct_v1));

    adios_buffer_struct_init (fh->b);
    MPI_Comm_rank (comm, &rank);
    rc = bp_read_open (fname, comm, fh);
    *fh_p = (int64_t) fh;
    if ( rank == 0 ) {
        bp_read_minifooter (fh);
    }
    MPI_Bcast (&fh->mfooter, sizeof(struct bp_minifooter),MPI_BYTE, 0, comm);
    
    uint64_t header_size = fh->mfooter.file_size-fh->mfooter.pgs_index_offset;

    if ( rank != 0) {
        if (!fh->b->buff) {
            alloc_aligned (fh->b, header_size);
            memset (fh->b->buff, 0, header_size);
            if (!fh->b->buff)
                fprintf(stderr, "could not allocate %d bytes\n",
                        header_size);
            fh->b->offset = 0;
        }
    }
    MPI_Barrier(comm);
    MPI_Bcast (fh->b->buff, fh->mfooter.file_size-fh->mfooter.pgs_index_offset,
           MPI_BYTE, 0 , comm);
    bp_parse_pgs (fh);
    bp_parse_vars (fh);
    bp_parse_attrs (fh);
    return rc;
}

int adios_fclose ( int64_t fh_p)
{
    struct BP_FILE * fh = (struct BP_FILE *) fh_p;
    struct BP_GROUP_VAR * gh = fh->gh;
    int i,j;
    MPI_File mpi_fh = fh->mpi_fh;

    if (fh->mpi_fh) 
        MPI_File_close (&mpi_fh);
    if (fh->b)
        adios_posix_close_internal (fh->b);
    if (gh) {
        for (j=0;j<2;j++) { 
            for (i=0;i<gh->group_count;i++) {
                if (gh->time_index[j][i])
                    free(gh->time_index[j][i]);
            }
            if (gh->time_index[j])
                free(gh->time_index[j]);
        }
        free (gh->time_index);
    
        for (i=0;i<gh->group_count;i++) { 
            if (gh->namelist[i])
                free(gh->namelist[i]);
        }
        if (gh->namelist)
            free (gh->namelist);

        for (i=0;i<fh->mfooter.vars_count;i++) 
            if (gh->var_namelist[i])
                free(gh->var_namelist[i]);
        if (gh->var_namelist)
            free (gh->var_namelist);

        if (gh->var_counts_per_group)
            free(gh->var_counts_per_group);

        if (gh->pg_offsets)
            free (gh->pg_offsets);

        free (gh);
    }
    if (fh)
        free (fh);    
    return;
}

int adios_gopen ( int64_t fh_p,
                  int64_t * gh_p,
                  char * grpname)
{
    struct BP_FILE * fh = (struct BP_FILE *) fh_p;
    struct BP_GROUP * bp_gh = (struct BP_GROUP *)
                malloc(sizeof(struct BP_GROUP));
    bp_gh->var_current = 0; 
    int i, grpid, offset=0;

    for (grpid=0;grpid<(fh->gh->group_count);grpid++) {
        if (!strcmp(fh->gh->namelist[grpid], grpname))
            break; 
    }

    for(i=0;i<grpid;i++)
        offset += fh->gh->var_counts_per_group[i];

    bp_gh->group_id = grpid;
    bp_gh->offset = offset;
    bp_gh->count = fh->gh->var_counts_per_group[grpid];
    bp_gh->fh = fh;
    *gh_p = (uint64_t) bp_gh;
    return;
}
                   
int adios_gclose ( int64_t gh_p)
{
    struct BP_GROUP * gh = (struct BP_GROUP *) gh_p;

    if (!gh) {
        fprintf(stderr, "group handle is NULL!\n");
        return  -2;
    }
    else
        free (gh);

    return 0;
}

int adios_inq_file ( int64_t fh_p, BP_FILE_INFO *pfinfo)
{
    if (!fh_p) {
        fprintf(stderr, "file handle is NULL!\n");
        return -2;
    }
    struct BP_FILE * fh = (struct BP_FILE *) fh_p;
    if (!fh->gh) {
        fprintf(stderr, "file handle is NULL!\n");
        return -2;
    }
    int i;
    pfinfo->groups_count = fh->gh->group_count;
    pfinfo->vars_count = fh->mfooter.vars_count;
    pfinfo->attrs_count = fh->mfooter.attrs_count;
    pfinfo->tidx_start = fh->tidx_start;
    pfinfo->tidx_stop = fh->tidx_stop;
    if (!pfinfo->namelist_true)
        return 0;
    if (!pfinfo->group_namelist)
        alloc_namelist (&pfinfo->group_namelist,pfinfo->groups_count); 
    for (i=0;i<pfinfo->groups_count;i++) {
        if (!pfinfo->group_namelist[i]) {
            fprintf(stderr, 
                "buffer given is too small, only hold %d entries",
                i);
            return -1;
        }
        else 
            strcpy(pfinfo->group_namelist[i],fh->gh->namelist[i]);
    }
    return 0;
}
/* 
int adios_inq_file_t ( int64_t fh_p, int *ngroup, 
          int *nvar, int *nattr, int *nt, char **gnamelist) 
{
    if (!fh_p) {
        fprintf(stderr, "file handle is NULL!\n");
        return -2;
    }
    struct BP_FILE * fh = (struct BP_FILE *) fh_p;
    int i;
    *ngroup = fh->gh->group_count;
    *nvar = fh->mfooter.vars_count;
    *nattr = fh->mfooter.attrs_count;
    *nt = fh->mfooter.time_steps;
    if (!gnamelist)
        return 0;
    for (i=0;i<fh->gh->group_count;i++) {
        if (!gnamelist[i]) {
            fprintf(stderr, 
                "buffer given is too small, only hold %d entries",
                i);
            return -1;
        }
        else 
            strcpy(gnamelist[i],fh->gh->namelist[i]);
    }
    return 0;
}
*/

int adios_inq_group (int64_t gh_p, BP_GROUP_INFO * pginfo)
{
    struct BP_GROUP * gh = (struct BP_GROUP *) gh_p;
    if (!gh_p) {
        fprintf(stderr, "group handle is NULL!\n");
        return -3;
    }
    int i, offset;

    pginfo->vars_count = gh->count;
    
    if (!pginfo->namelist_true)
        return 0;

    offset = gh->offset;
    alloc_namelist (&(pginfo->var_namelist), pginfo->vars_count);
    for (i=0;i<pginfo->vars_count;i++) {
        if (!pginfo->var_namelist[i]) { 
            fprintf(stderr, 
                    "given buffer only can hold %d entries",
                    i);
            return -1;
        }
        else
            strcpy(pginfo->var_namelist[i], gh->fh->gh->var_namelist[i+offset]);
    }

    return 0;
}
#if 0
int adios_inq_group (int64_t gh_p, int *nvar, char ** vnamelist)
{
    struct BP_GROUP * gh = (struct BP_GROUP *) gh_p;
    if (!gh_p) {
        fprintf(stderr, "group handle is NULL!\n");
        return -3;
    }
    int i, offset;

    * nvar = gh->count;
    
    if (!vnamelist)
        return 0;

    offset = gh->offset;

    for (i=0;i<*nvar;i++) {
        if (!vnamelist[i]) { 
            fprintf(stderr, 
                    "given buffer only can hold %d entries",
                    i);
            return -1;
        }
        else
            strcpy(vnamelist[i], gh->fh->gh->var_namelist[i+offset]);
    }
    return 0;    
}
#endif    

/** Find a string (variable by full name) in a list of strings (variable name list)
  * from an 'offset' within the first 'count' elements.
  * Return the index of the variable in [offset .. offset+count-1] if found,
  *  or -1 if not found.
  */
static int find_var( char ** varnamelist, int offset, int count, char * varname) 
{
    // Find the variable: full path is stored with a starting / 
    // Like in HDF5, we need to match names given with or without the starting /
    // startpos is 0 or 1 to indicate if the argument has starting / or not
    int i, endp;
    int vstartpos = 0, fstartpos = 0; 
    if (varname[0] == '/')
        vstartpos = 1;
        
    endp = offset + count;
    for (i = offset; i < endp; i++) {
        if (varnamelist[i][0] == '/')
            fstartpos = 1;
        if (!strcmp( varnamelist[i]+fstartpos, varname+vstartpos )) {
            return i;
        }
    }
    return -1;
}

int adios_inq_var (int64_t gh_p, char * varname,
         int * type,
         int * ndim,
         int * is_timebased,
         int * dims)
{
    struct BP_GROUP * gh = (struct BP_GROUP *) gh_p;
    if (!gh_p) {
        fprintf(stderr, "group handle is NULL!\n");
        return -3;
    }
    struct BP_FILE * fh = gh->fh;
    if (!fh) {
        fprintf(stderr, "file handle is NULL!\n");
        return -2;
    }
    
    struct adios_index_var_struct_v1 * var_root;
    int i,k, var_id;
    gh->var_current = 0;
    var_root = fh->vars_root;
    *is_timebased = 0;
    if (!varname) {
        fprintf(stderr, "Error: Variable %s is NULL!\n");
        return -4;
    }

    // find variable in var list
    var_id = find_var(fh->gh->var_namelist, gh->offset, gh->count, varname);

    for (i=0;i<gh->group_id;i++)
        var_id -= fh->gh->var_counts_per_group[i];
    for(i=0;i<gh->offset;i++)
        var_root = var_root->next;
    if (var_id < 0) {
        fprintf(stderr, 
            "Error: Variable %s does not exist in the group %s!\n",
                    varname, fh->gh->namelist[gh->group_id]);
        return -4;
    }
    for (i=0;i<var_id;i++) {
        var_root = var_root->next;
    }

    gh->var_current = var_root;
    *type = var_root->type;
    if (!var_root->characteristics_count) {
        fprintf(stderr, 
            "Error: Variable %s does not exist in the file!\n",
                    varname);
        *ndim = -1;
        return -4;
    }
        
    *ndim = var_root->characteristics [0].dims.count;
    if (!dims || !(*ndim)) {
        return 0;
    }

    int time_flag = -1;
    uint64_t * gdims = (uint64_t *) malloc (sizeof(uint64_t) * (*ndim));
    uint64_t * ldims = (uint64_t *) malloc (sizeof(uint64_t) * (*ndim));
    int is_global=0;
    memset(dims,0,sizeof(int)*(*ndim));

    for (k=0;k<(*ndim);k++) {
        gdims[k]=var_root->characteristics[0].dims.dims[k*3+1];
        ldims[k]=var_root->characteristics[0].dims.dims[k*3];
    }
    for (i=0;i<(*ndim);i++) {
         is_global = is_global || gdims[i];
    }
    if (!is_global) {
        for (i=0;i<(*ndim);i++) {
            if (   ldims[i] == 1 
                  && var_root->characteristics_count > 1) {
                *is_timebased = 1;
                time_flag = i;
            }
            if (dims) {
                if (time_flag==0) {
                    if (i>0)
                        dims[i-1]=ldims[i];
                }
                else { 
                    dims[i]=ldims[i];
                }
            }
        }         
    }         
    else {
        if (ldims[0]==1 && ldims[*ndim - 1]!=1) {
            time_flag = 0;
            *is_timebased = 1;
        }
        else if (ldims[0]!=1 && ldims[*ndim - 1]==1) {
            time_flag = *ndim - 1;
            *is_timebased = 1;
        }
        if (*is_timebased) {
            for (i=0;i<(*ndim)-1;i++) {
                if (dims) 
                    dims[i]=gdims[i];
            }
        }
        else {
            for (i=0;i<(*ndim);i++) 
                dims[i]=gdims[i];
        }
    }
    free(gdims);
    free(ldims);
            
    if ((*is_timebased))
        *ndim = *ndim - 1;
    
    return 0;
}

int64_t adios_get_var (int64_t gh_p,
                       char * varname, 
                       void * var,
                       int  * start,
                       int  * readsize,
                       int timestep
                      )
{
    double  start_time, stop_time;
    int    i, j, var_id, idx;
    int    offset, count, start_idx=-1;
    struct BP_GROUP * gh = (struct BP_GROUP *) gh_p;
    struct BP_FILE * fh = (struct BP_FILE *) (gh->fh);
    struct adios_index_var_struct_v1 * var_root = fh->vars_root;
    struct adios_var_header_struct_v1 var_header;
    struct adios_var_payload_struct_v1 var_payload;
    enum ADIOS_FLAG adios_host_language_fortran = fh->pgs_root->adios_host_language_fortran;
    uint8_t  ndim;  
    uint64_t tmpreadsize;
    uint64_t size, * ldims, * offsets, * gdims;
    uint64_t datasize, nloop, dset_stride,var_stride, total_size=1;
    int readsize_tmp[10], start_tmp[10];
    MPI_Status status;

    var_id = find_var(fh->gh->var_namelist, gh->offset, gh->count, varname);
    for (i=0;i<gh->group_id;i++)
        var_id -= fh->gh->var_counts_per_group[i];

    for(i=0;i<gh->offset && var_root;i++)
        var_root = var_root->next;

    if (var_id<0) {
        fprintf(stderr, "Error: Variable %s does not exist in the group %s!\n",
                    varname, fh->gh->namelist[gh->group_id]);
        return -4;
    }

    for (i=0;i<var_id && var_root;i++) 
        var_root = var_root->next;
    
    if (i!=var_id) {
        fprintf(stderr, "Error: time step should start from 1 %d %d!\n",
                i, var_id);
        return -5; 
    }

    gh->var_current = var_root;

    if (timestep < 0) {
        fprintf(stderr, "Error: time step should start from 1!\n");
        return -5; 
    }
    if (timestep < fh->tidx_start) {
        fprintf(stderr, "Error: time step should start from 1:%d %d!\n",
                timestep, fh->tidx_start);
        return -5; 
    }
    // get the starting offset for the given time step
    offset = fh->gh->time_index[0][gh->group_id][timestep - fh->tidx_start];
    count = fh->gh->time_index[1][gh->group_id][timestep - fh->tidx_start];
    for (i=0;i<var_root->characteristics_count;i++) {
        if (   (  var_root->characteristics[i].offset 
                > fh->gh->pg_offsets[offset])
            && (  (i == var_root->characteristics_count-1) 
                ||(  var_root->characteristics[i].offset 
                   < fh->gh->pg_offsets[offset + 1]))
           ) {
            start_idx = i;
            break;
        }
    }

    if (start_idx<0) {
        fprintf(stderr,"Error: %s has no data at %d time step\n",
            varname, timestep);
        return -4;
    }
    ndim = var_root->characteristics[start_idx].dims.count;
    ldims = (uint64_t *) malloc (sizeof(uint64_t) * ndim);
    gdims = (uint64_t *) malloc (sizeof(uint64_t) * ndim);
    offsets = (uint64_t *) malloc (sizeof(uint64_t) * ndim);
    int * idx_table = (int *) malloc (sizeof(int) * count);
    
    // main for loop
    int time_flag = -1;
    int rank;
    int is_global=0, is_timebased = 0;
    MPI_Comm_rank(gh->fh->comm, &rank);

    for (i=0;i<ndim;i++) {
        gdims[i]=var_root->characteristics[0].dims.dims[i * 3 + 1];
        offsets[i]=var_root->characteristics[0].dims.dims[i * 3 + 2];
        ldims[i]=var_root->characteristics[0].dims.dims[i * 3];
        is_global = is_global || gdims[i];
    }

    if (!is_global) {
        for (i = 0; i < ndim; i++) {
            if (   ldims[i] == 1 
                  && var_root->characteristics_count > 1
                && ndim>1) {
                is_timebased = 1;
                time_flag = i;
            }
        }    
    }         
    else {
        if (ldims[0]==1 && ldims[ndim - 1]!=1) {
            time_flag = 0;
            is_timebased = 1;
        }
        else if (ldims[0]!=1 && ldims[ndim - 1]==1) {
            time_flag = ndim - 1;
            is_timebased = 1;
        }
    }
    if (is_timebased) {
        if ((ldims[0] != 1) && (ldims[ndim - 1] == 1)) 
            time_flag = ndim - 1;
        else if ((ldims[0] == 1) && (ldims[ndim - 1] != 1)) { 
            time_flag = 0;
        }
        ndim = ndim -1;
    }

    // if bp is written by fortran, flip read size and offset
    if (adios_host_language_fortran)
    {
        if (readsize)
            tmpreadsize = (uint64_t) readsize[ndim - 1];

        for (i = 0; i < ndim; i++)
        {
            readsize_tmp[i] = readsize[ndim - 1 - i];
            start_tmp[i] = start[ndim -1 - i];
        }
    }
    else
    {
        if (readsize)
            tmpreadsize = (uint64_t) readsize[0];

        for (i = 0; i < ndim; i++)
        {
            readsize_tmp[i] = readsize[i];
            start_tmp[i] = start[i];
        }
    }

    for (i = 0; i < ndim; i++)
        total_size *= readsize_tmp[i];

    // generate the list of pgs to be read from

    uint64_t read_offset = 0;
    int npg = 0;
    uint64_t tmpcount = 0;
    int size_of_type = bp_get_type_size (var_root->type, "");

    // actions
    if (count > var_root->characteristics_count)
        count = var_root->characteristics_count;
    for (idx = 0; idx < count; idx++) {
        int flag;
        datasize = 1;
        nloop = 1;
        var_stride = 1;
        dset_stride = 1;
        idx_table[idx] = 1;
        uint64_t payload_size = size_of_type;

        for (j = 0; j< ndim; j++) {
            if (adios_host_language_fortran)
            {
                offsets[j]=var_root->characteristics[start_idx + idx].dims.dims[(ndim -1 - j) * 3 + 2];
                if (!time_flag)
                    ldims[j]=var_root->characteristics[start_idx + idx].dims.dims[(ndim - 1 - j) * 3 + 3];
                else
                    ldims[j]=var_root->characteristics[start_idx + idx].dims.dims[(ndim - 1 - j) * 3];

                if (is_global)
                    gdims[j]=var_root->characteristics[start_idx + idx].dims.dims[(ndim - 1 - j) * 3 + 1];
                else
                    gdims[j]=ldims[j];
            }
            else
            {
                offsets[j]=var_root->characteristics[start_idx+idx].dims.dims[j*3+2];
                if (!time_flag)
                    ldims[j]=var_root->characteristics[start_idx+idx].dims.dims[j*3+3];
                else
                    ldims[j]=var_root->characteristics[start_idx+idx].dims.dims[j*3];

                if (is_global)
                    gdims[j]=var_root->characteristics[start_idx+idx].dims.dims[j*3+1];
                else
                    gdims[j]=ldims[j];
            }

            payload_size *= ldims [j];

            if ( (readsize_tmp[j] > gdims[j]) 
              || (start_tmp[j] > gdims[j]) 
              || (start_tmp[j] + readsize_tmp[j] > gdims[j])){
                fprintf(stderr, "Error: %s out of bound ("
                    "the data to read is [%llu,%llu]"
                    " but the actual data is [0,%llu])\n",
                    varname, start_tmp[j], readsize_tmp[j], gdims[j] - 1);
                return -3;
            }

            flag = (offsets[j] >= start_tmp[j] 
                    && offsets[j] < start_tmp[j] + readsize_tmp[j])
                || (offsets[j] < start_tmp[j]
                        && offsets[j] + ldims[j]>start_tmp[j] + readsize_tmp[j]) 
                || (offsets[j] + ldims[j] > start_tmp[j] 
                        && offsets[j] + ldims[j] <= start_tmp[j] + readsize_tmp[j]);
            idx_table [idx] = idx_table[idx] && flag;
        }
        
        if ( !idx_table[idx] ) {
            continue;
        }
        ++npg;

        int hole_break;
        for (i = ndim - 1; i > -1; i--) {
            if (offsets[i] == start_tmp[i] && ldims[i] == readsize_tmp[i]) {
                datasize *= ldims[i];
            }
            else
                break;
        }

        hole_break = i;
        uint64_t slice_offset = 0;
        uint64_t slice_size = 0;

#define MPI_FILE_READ_OPS                           \
        realloc_aligned(fh->b, slice_size);         \
        fh->b->offset = 0;                          \
                                                    \
        MPI_File_seek (fh->mpi_fh                   \
                      ,(MPI_Offset)slice_offset     \
                      ,MPI_SEEK_SET                 \
                      );                            \
                                                    \
        MPI_File_read (fh->mpi_fh                   \
                      ,fh->b->buff                  \
                      ,slice_size                   \
                      ,MPI_BYTE                     \
                      ,&status                      \
                      );                            \
        fh->b->offset = 0;                          \

//We also need to be able to read old .bp which doesn't have 'payload_offset'
#define MPI_FILE_READ_OPS1                                                                  \
        MPI_File_seek (fh->mpi_fh                                                           \
                      ,(MPI_Offset) var_root->characteristics[start_idx + idx].offset       \
                      ,MPI_SEEK_SET);                                                       \
        MPI_File_read (fh->mpi_fh, fh->b->buff, 8, MPI_BYTE, &status);                      \
        tmpcount= *((uint64_t*)fh->b->buff);                                                \
                                                                                            \
        realloc_aligned(fh->b, tmpcount + 8);                                               \
        fh->b->offset = 0;                                                                  \
                                                                                            \
        MPI_File_seek (fh->mpi_fh                                                           \
                      ,(MPI_Offset) (var_root->characteristics[start_idx + idx].offset)     \
                      ,MPI_SEEK_SET);                                                       \
        MPI_File_read (fh->mpi_fh, fh->b->buff, tmpcount + 8, MPI_BYTE, &status);           \
        fh->b->offset = 0;                                                                  \
        adios_parse_var_data_header_v1 (fh->b, &var_header);                                \

        if (hole_break == -1) {
            slice_size = payload_size;

            if (var_root->characteristics[start_idx + idx].payload_offset > 0) {
                slice_offset = var_root->characteristics[start_idx + idx].payload_offset;
                MPI_FILE_READ_OPS
            } else {
                slice_offset = 0;
                MPI_FILE_READ_OPS1
            }

            memcpy(var, fh->b->buff + fh->b->offset, slice_size);
        }
        else if (hole_break == 0) 
        {
            if (tmpreadsize > ldims[0])
                slice_size = payload_size;
            else
                slice_size = tmpreadsize * datasize * size_of_type;

            if (var_root->characteristics[start_idx + idx].payload_offset > 0) {
                slice_offset = var_root->characteristics[start_idx + idx].payload_offset;
                MPI_FILE_READ_OPS
            } else {
                slice_offset = 0;
                MPI_FILE_READ_OPS1
            }

            memcpy (var + read_offset, fh->b->buff + fh->b->offset, slice_size);

            read_offset +=  slice_size;
            if (tmpreadsize > ldims[0])
                tmpreadsize -= ldims[0];
        }
        else 
        {
            uint64_t stride_offset = 0;
            int isize;
            uint64_t size_in_dset[10];
            uint64_t offset_in_dset[10];
            uint64_t offset_in_var[10];
            memset(size_in_dset, 0 , 10 * 8);
            memset(offset_in_dset, 0 , 10 * 8);
            memset(offset_in_var, 0 , 10 * 8);
            int hit = 0;
            for ( i = 0; i < ndim ; i++) {
                isize = offsets[i] + ldims[i];
                if (start_tmp[i] >= offsets[i]) {
                    // head is in
                    if (start_tmp[i]<isize) {
                        if (start_tmp[i] + readsize_tmp[i] > isize)
                            size_in_dset[i] = isize - start_tmp[i];
                        else
                            size_in_dset[i] = readsize_tmp[i];
                        offset_in_dset[i] = start_tmp[i] - offsets[i];
                        offset_in_var[i] = 0;
                        hit = 1 + hit * 10;
                    }
                    else
                        hit = -1;
                }
                else {
                    // middle is in
                    if (isize < start_tmp[i] + readsize_tmp[i]) {
                        size_in_dset[i] = ldims[i];
                        hit = 2 + hit * 10;
                    }
                    else {
                        // tail is in
                        size_in_dset[i] = readsize_tmp[i] + start_tmp[i] - offsets[i];
                        hit = 3 + hit * 10;
                    }
                    offset_in_dset[i] = 0;
                    offset_in_var[i] = offsets[i] - start_tmp[i];
                }
            }

            datasize = 1;
            var_stride = 1;

            for ( i = ndim-1; i >= hole_break; i--) {
                datasize *= size_in_dset[i];
                dset_stride *= ldims[i];
                var_stride *= readsize_tmp[i];
            }

            uint64_t start_in_payload = 0, end_in_payload = 0, s = 1;
            for (i = ndim - 1; i > -1; i--) {
                start_in_payload += s * offset_in_dset[i] * size_of_type;
                end_in_payload += s * (offset_in_dset[i] + size_in_dset[i] - 1) * size_of_type;
                s *= ldims[i];
            }

            printf ("start in payload = %llu, end in payload = %llu\n", start_in_payload, end_in_payload); 
            slice_size = end_in_payload - start_in_payload + 1 * size_of_type;

            if (var_root->characteristics[start_idx + idx].payload_offset > 0) {
                slice_offset =  var_root->characteristics[start_idx + idx].payload_offset
                              + start_in_payload;
                MPI_FILE_READ_OPS

                for ( i = 0; i < ndim ; i++) {
                    offset_in_dset[i] = 0;
                }
            } else {
                slice_offset =  start_in_payload;
                MPI_FILE_READ_OPS1
            }

            uint64_t var_offset = 0;
            uint64_t dset_offset = 0;
            for ( i = 0; i < hole_break; i++) {
                nloop *= size_in_dset[i];
            }

            for ( i = 0; i < ndim ; i++) {
                var_offset = offset_in_var[i] + var_offset * readsize_tmp[i];
                dset_offset = offset_in_dset[i] + dset_offset * ldims[i];
            }

            copy_data (var
                      ,fh->b->buff + fh->b->offset
                      ,0
                      ,hole_break
                      ,size_in_dset
                      ,ldims
                      ,readsize
                      ,var_stride
                      ,dset_stride
                      ,var_offset
                      ,dset_offset
                      ,datasize
                      ,size_of_type 
                      );
        }
    }  // end of loop

    free (gdims);
    free (offsets);
    free (ldims);
    free (idx_table);

    return total_size * size_of_type;
}

const char * adios_type_to_string (int type)
{
    switch (type)
    {
        case adios_unsigned_byte:    return "unsigned byte";
        case adios_unsigned_short:   return "unsigned short";
        case adios_unsigned_integer: return "unsigned integer";
        case adios_unsigned_long:    return "unsigned long long";

        case adios_byte:             return "byte";
        case adios_short:            return "short";
        case adios_integer:          return "integer";
        case adios_long:             return "long long";

        case adios_real:             return "real";
        case adios_double:           return "double";
        case adios_long_double:      return "long double";

        case adios_string:           return "string";
        case adios_complex:          return "complex";
        case adios_double_complex:   return "double complex";

        default:
        {
            static char buf [50];
            sprintf (buf, "(unknown: %d)", type);
            return buf;
        }
    }
}
void adios_init_groupinfo(BP_GROUP_INFO * pginfo, int flag)
{
    if (pginfo) {
        memset (pginfo, 0, sizeof(BP_GROUP_INFO));
        pginfo->namelist_true = flag;
    }
    else 
        printf ("fileinfo is NULL\n");
    return;    

}

void adios_free_groupinfo (BP_GROUP_INFO * pginfo)
{
    free_namelist ((pginfo->var_namelist),pginfo->vars_count);
    return;
}

void adios_print_groupinfo (BP_GROUP_INFO *pginfo) 
{
    int i;
    printf ("---------------------------\n");
    printf ("     var information\n");
    printf ("---------------------------\n");
    printf ("    var id\tname\n");
    if (pginfo->var_namelist) {
        for (i=0; i<pginfo->vars_count; i++)
            printf("\t%d)\t%s\n", i, pginfo->var_namelist[i]);
    }
    return;
}

void adios_init_fileinfo(BP_FILE_INFO * pfinfo, int flag)
{
    if (pfinfo) {
        memset (pfinfo, 0, sizeof(BP_FILE_INFO));
        pfinfo->namelist_true = flag;
    }
    else 
        printf ("fileinfo is NULL\n");
    return;    
}

void adios_free_fileinfo (BP_FILE_INFO * pfinfo)
{
    free_namelist ((pfinfo->group_namelist),pfinfo->groups_count);
    return;
}

void adios_print_fileinfo (BP_FILE_INFO *pfinfo) 
{
    int i;
    printf ("---------------------------\n");
    printf ("     group information\n");
    printf ("---------------------------\n");
    printf ("\t# of groups:\t%d\n"
         "\t# of variables:\t%d\n"
        "\t# of attributes:%d\n"
        "\t# of timesteps:\t%d-->%d\n",
        pfinfo->groups_count,
        pfinfo->vars_count,
        pfinfo->attrs_count,
        pfinfo->tidx_start,
        pfinfo->tidx_stop);
    printf ("\t----------------\n");
    printf ("\t  group id\tname\n");
    if (pfinfo->group_namelist) {
        for (i=0; i<pfinfo->groups_count; i++)
            printf("\t\t%d)\t%s\n", i, pfinfo->group_namelist[i]);
    }
    return;
}

void adios_fopen_(int64_t * fh_p,
                  char * fname,
                  void * fcomm,
                  int * err,
                  int fname_len
                 )
{
    int64_t fh;
    MPI_Comm comm = MPI_Comm_f2c (*((int *) fcomm));
    *err = adios_fopen (&fh,fname,comm);
    * fh_p = fh;
}

void adios_fclose_( int64_t * fh, int * err)
{
    *err = adios_fclose (*fh);
}
void adios_inq_file_ ( int64_t * fh_p,
                       int * groups_count,
                       int * vars_count,
                       int * attrs_count,
                       int * tstart,
                       int * tstop,
                       void * gnamelist,
                       int * err,
                       int gnamelist_len)
{
    BP_FILE_INFO finfo;
    int i;
    adios_init_fileinfo ( &finfo, 1);
    adios_inq_file ( *fh_p, &finfo);
    * groups_count = finfo.groups_count;
    * vars_count = finfo.vars_count;
    * attrs_count = finfo.attrs_count;
    * tstart = finfo.tidx_start;
    * tstop = finfo.tidx_stop;
    * err = 0;
    for (i=0;i<*groups_count;i++)
        strcpy(gnamelist+i*gnamelist_len,finfo.group_namelist[i]); 
    adios_free_fileinfo(&finfo);
}

void adios_gopen_ ( int64_t * fh,
                    int64_t * gh_p,
                    char * grpname,
                    int * err,
                    int grpname_len)
{
    int64_t gh=0;
    *err=adios_gopen (*fh,&gh,grpname);
    *gh_p = gh;
}
void adios_gclose_( int64_t * gh, int * err)
{
    *err=adios_gclose(*gh);
}
void adios_inq_group_ (int64_t * gh, int *vcnt, void *vnamelist, 
        int *err, int vnamelist_len) 
{
    BP_GROUP_INFO ginfo;
    int i;
    adios_init_groupinfo ( &ginfo, 1);
    adios_inq_group ( *gh, &ginfo);
    * vcnt = ginfo.vars_count;
    * err = 0;
    for (i=0;i<*vcnt;i++) {
        strncpy(vnamelist+i*vnamelist_len,ginfo.var_namelist[i],
                strlen(ginfo.var_namelist[i]));
        *((char*)(vnamelist+i*vnamelist_len+strlen(ginfo.var_namelist[i])))='\0';
    } 
    adios_free_groupinfo(&ginfo);
}

void adios_get_var_ (int64_t * gh,
        char * varname,
        void * var,
        int  * start,
        int  * readsize,
        int * timestep,
        int64_t * err,
        int varname_len)
{
    *err=adios_get_var (*gh,varname, var,start,readsize,*timestep);
}

void adios_inq_var_ (int64_t * gh_p, char * varname,
        int * type,
        int * ndim,
        int * is_timebased,
        int * dims,
        int * err,
        int varname_len)
{
    int vtype, vndim, vtimed;
    int vdims[10], i;
    *err = adios_inq_var (*gh_p, varname, &vtype, &vndim, &vtimed, vdims);
    *type = vtype;
    *ndim = vndim;
    *is_timebased = vtimed;
    for (i=0;i<vndim;i++)
        dims[i] = vdims[i];
}
