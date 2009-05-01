#include <stdlib.h>
#include <string.h>
#include "adios.h"
#include "bp_utils.h"
#include "bp_read.h"
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

int bp_fopen ( int64_t * fh_p,
	        const char * fname,
		MPI_Comm comm
	      )
{
	int rc, rank;	
	struct BP_FILE * fh = (struct BP_FILE *)
                               malloc (sizeof (struct BP_FILE));
	fh->comm = comm;
	fh->gh = 0;
	fh->pgs_root = 0;
	fh->vars_root = 0;
	fh->attrs_root = 0;
	fh->b = malloc (sizeof (struct adios_bp_buffer_struct_v1));

	adios_buffer_struct_init (fh->b);

	rc = bp_read_open (fname, comm, fh);

	*fh_p = (int64_t) fh;
	
	MPI_Comm_rank (comm, &rank);
	if ( rank == 0 ) {
		//bp_read_minifooter (fh->b, &(fh->mfooter));
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

int bp_fclose ( int64_t fh_p)
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

int bp_gopen ( int64_t * gh_p,
		int64_t fh_p, 
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
	   			
int bp_gclose ( int64_t gh_p)
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

int bp_inq_file ( int64_t fh_p, BP_FILE_INFO *pfinfo)
{
	if (!fh_p) {
		fprintf(stderr, "file handle is NULL!\n");
		return -2;
	}
	struct BP_FILE * fh = (struct BP_FILE *) fh_p;
	int i;
	pfinfo->groups_count = fh->gh->group_count;
	pfinfo->vars_count = fh->mfooter.vars_count;
	pfinfo->attrs_count = fh->mfooter.attrs_count;
	pfinfo->tidx_start = fh->tidx_start;
	pfinfo->tidx_stop = fh->tidx_stop;
	if (!pfinfo->namelist_true)
		return 0;
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
int bp_inq_file_t ( int64_t fh_p, int *ngroup, 
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

int bp_inq_group (int64_t gh_p, BP_GROUP_INFO * pginfo)
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
int bp_inq_group (int64_t gh_p, int *nvar, char ** vnamelist)
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
        /* Find the variable: full path is stored with a starting / 
           Like in HDF5, we need to match names given with or without the starting /
           startpos is 0 or 1 to indicate if the argument has starting / or not
        */
        int i, endp;
        int startpos = 0; 
        if (varname[0] != '/')
        	startpos = 1;
        
        endp = offset + count;
	for (i = offset; i < endp; i++) {
                //printf("-- find_var: compare %s with %s, startpos=%d\n", varname, varnamelist[i], startpos);
		if (!strcmp( &(varnamelist[i][startpos]), varname)) {
			return i;
		}
	}
        return -1;
}

int bp_inq_var (int64_t gh_p, char * varname,
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
	for(i=0;i<gh->offset;i++)
		var_root = var_root->next;
	if (var_id<0) {
		fprintf(stderr, 
			"Error: Variable %s does not exist in the group %s!\n",
                	varname, fh->gh->namelist[gh->group_id]);
		return -4;
	}

	for (i=0;i<var_id;i++) 
		var_root = var_root->next;

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
#if 0
	printf("var: %s,id=%d %d %d %d\n", var_root->var_name,
        var_id,
		*ndim, 
		var_root->characteristics[0].value,
		var_root->characteristics_count);
#endif
	if (!dims || !(*ndim)) {
        //printf(stderr, "scalar or dims is NULL\n");
		return 0;
    }

    int time_flag = -1;
	uint64_t * gdims = (uint64_t *) malloc (sizeof(uint64_t) * (*ndim));
	uint64_t * ldims = (uint64_t *) malloc (sizeof(uint64_t) * (*ndim));
	memset(dims,0,sizeof(int)*(*ndim));

	for (k=0;k<(*ndim);k++) {
		gdims[k]=var_root->characteristics[0].dims.dims[k*3+1];
		ldims[k]=var_root->characteristics[0].dims.dims[k*3];
	}
	int is_global=0;
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
        if (ldims[0]==1 && ldims[*ndim]!=1) {
            time_flag = 0;
            *is_timebased = 1;
        }
        else if (ldims[0]!=1 && ldims[*ndim]==1) {
            time_flag = *ndim;
            *is_timebased = 1;
        }
        if (is_timebased) {
            for (i=0;i<(*ndim)-1;i++) {
                //printf("time_flag:%d\n",time_flag);
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

int bp_get_var (int64_t gh_p,
		 char * varname, 
		 void * var,
		 int  * start,
		 int  * readsize,
		 int timestep
		)
{
	
    double  start_time, stop_time;
	int    i, j, var_id, idx;
	int    flag, offset, count, start_idx=-1;
	struct BP_GROUP * gh = (struct BP_GROUP *) gh_p;
	struct BP_FILE * fh = (struct BP_FILE *) (gh->fh);
	struct adios_index_var_struct_v1 * var_root = fh->vars_root;
    struct adios_var_header_struct_v1 var_header;
    struct adios_var_payload_struct_v1 var_payload;
	uint8_t  ndim;  
  	uint64_t tmpreadsize = (uint64_t) readsize[0];
	uint64_t size, * ldims, * offsets, * gdims;
	uint64_t datasize, nloop, dset_stride,var_stride, total_size=1;
	MPI_Status status;

	var_id = find_var(fh->gh->var_namelist, gh->offset, gh->count, varname);
	if (var_id<0) {
		fprintf(stderr, "Error: Variable %s does not exist in the group %s!\n",
                	varname, fh->gh->namelist[gh->group_id]);
		return -4;
	}

	for (i=0;i<var_id && var_root;i++) 
		var_root = var_root->next;

    if (i!=var_id) {
		fprintf(stderr, "Error: time step should start from 1!");
		return -5; 
    }

	gh->var_current = var_root;

	if (timestep < 0) {
		fprintf(stderr, "Error: time step should start from 1!");
		return -5; 
	}
	if (timestep< fh->tidx_start) {
		fprintf(stderr, "Error: time step should start from 1!");
		return -5; 
	}
	// get the starting offset for the given time step
	offset = fh->gh->time_index[0][gh->group_id][timestep-fh->tidx_start];
    count = fh->gh->time_index[1][gh->group_id][timestep-fh->tidx_start];
    //printf("var_id:%d offset: %d count:%d varname: %s\n",
    //        var_id,offset,count,varname);

	for (i=0;i<var_root->characteristics_count;i++) {
        if (   (  var_root->characteristics[i].offset 
                > fh->gh->pg_offsets[offset])
            && (  (i == var_root->characteristics_count-1) 
                ||(  var_root->characteristics[i].offset 
                   < fh->gh->pg_offsets[offset+1]))
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
    int * idx_table = (int *) malloc (sizeof(int)*count);
	
	// main for loop
	int time_flag = -1;
	int rank;
	int is_global=0, is_timebased = 0;
	MPI_Comm_rank(gh->fh->comm, &rank);

	for (i=0;i<ndim;i++) {
		gdims[i]=var_root->characteristics[0].dims.dims[i*3+1];
		offsets[i]=var_root->characteristics[0].dims.dims[i*3+2];
		ldims[i]=var_root->characteristics[0].dims.dims[i*3];
	}

	for (i=0;i<ndim;i++) {
		is_global = is_global || gdims[i];
	}
	if (!is_global) {
		for (i=0;i<(ndim);i++) {
			if (   ldims[i] == 1 
		  	    && var_root->characteristics_count > 1) {
				is_timebased = 1;
				time_flag = i;
			}
		}		 
	}		 
	else {
        if (ldims[0]==1 && ldims[ndim]!=1) {
            time_flag = 0;
            is_timebased = 1;
        }
        else if (ldims[0]!=1 && ldims[ndim]==1) {
            time_flag = ndim;
            is_timebased = 1;
        }
	}
	if (is_timebased) {
        if ((ldims[0] != 1) && (ldims[ndim] == 1)) 
            time_flag == ndim;
        else if ((ldims[0] == 1) && (ldims[ndim] != 1)) { 
            time_flag = 0;
        }
		ndim = ndim -1;
    }

	// generate the list of pgs to be read from

	uint64_t read_offset = 0;
	int npg=0;
	int tmpcount = 0;
	int size_of_type = bp_get_type_size (var_root->type, "");

    /* actions */
	//printf("rank=%d offset count: %d %d \n",rank, offset, count);
	//printf("rank=%d ndim: %d \n",rank, ndim);

	for (idx = 0; idx < count; idx++) {
		datasize = 1;
		nloop = 1;
		var_stride = 1;
		dset_stride = 1;
		idx_table[idx] = 1;
		for (j=0;j<ndim;j++) {
			offsets[j]=var_root->characteristics[start_idx+idx].dims.dims[j*3+2];
	        if (!time_flag)
				ldims[j]=var_root->characteristics[start_idx+idx].dims.dims[j*3+3];
			else
				ldims[j]=var_root->characteristics[start_idx+idx].dims.dims[j*3];
			if (is_global)
				gdims[j]=var_root->characteristics[start_idx+idx].dims.dims[j*3+1];
			else
				gdims[j]=ldims[j];
            /*
		    printf("rank=%d var_name:%s global: %llu offset:%llu local:%llu\n",
                    rank, varname,gdims[j],offsets[j],ldims[j]);
            printf("timeflag=%d is_global=%d\n",time_flag,is_global);
            */
			if (readsize[j] > gdims[j]) {
				fprintf(stderr, "Error: %s out of bound ("
					"the size to read is %llu,"
					" but the actual size is %llu)\n",
					varname, readsize[j], gdims[j]);
				return -3;
			}

			flag = (offsets[j] >= start[j] 
					&& offsets[j] <start[j]+readsize[j])
				|| (offsets[j] < start[j]
						&& offsets[j]+ldims[j]>start[j]+readsize[j]) 
				|| (offsets[j]+ldims[j] > start[j] 
						&& offsets[j]+ldims[j] <= start[j]+readsize[j]);
			idx_table [idx] = idx_table[idx] && flag;
		}
		
		if ( !idx_table[idx] ) {
			continue;
		}
        ++npg;
        MPI_Barrier(MPI_COMM_WORLD);
        start_time = MPI_Wtime();
		MPI_File_seek (fh->mpi_fh, 
			       (MPI_Offset) var_root->characteristics[start_idx+idx].offset,
			       MPI_SEEK_SET);
		MPI_File_read (fh->mpi_fh, fh->b->buff, 4, MPI_BYTE, &status);
		tmpcount= *((uint64_t*)fh->b->buff);	
		realloc_aligned(fh->b, tmpcount+4);

		MPI_File_seek (fh->mpi_fh, 
			       (MPI_Offset) var_root->characteristics[start_idx+idx].offset,
			       MPI_SEEK_SET);
		MPI_File_read (fh->mpi_fh, fh->b->buff, tmpcount+4, MPI_BYTE, &status);
		//MPI_Get_count (&status, MPI_BYTE, &tmpcount);
	    //printf("read back: %d\n",tmpcount);	
		fh->b->offset = 0;
		adios_parse_var_data_header_v1 (fh->b, &var_header);

        MPI_Barrier(MPI_COMM_WORLD);
        stop_time = MPI_Wtime();
        printf("%f\t",stop_time-start_time);
		//--data filtering--//
		if (!read_offset)
			for (i = 0; i < ndim; i++) 
				total_size *= readsize[i];

		int hole_break;
		for (i=ndim-1;i>-1;i--) {
			if (ldims[i] == readsize[i])
				datasize *= ldims[i];
			else 
				break;
		}
		hole_break = i;
/*        
		printf("rank=%d hole break:%d ndim:%d\n",rank,hole_break,ndim);
        if (ndim == 1) {
		    printf("time_flag:%d ndim=%d\n",time_flag,ndim);
		    printf("local: %llu\n",ldims[0]);
		    printf("global: %llu\n",gdims[0]);
		    printf("offset: %llu\n",offsets[0]);
		    printf("readsize:%d\n",readsize[0]);
		    printf("start:%d\n",start[0]);
        }
        MPI_Barrier(gh->fh->comm);
*/        
        /*
        if (rank == 0 && ndim == 2) {
		    printf("time_flag:%d ndim=%d\n",time_flag,ndim);
		    printf("local: %d %d\n",ldims[0],ldims[1]);
		    printf("global: %d %d\n",gdims[0],gdims[1]);
		    printf("offset: %d %d\n",offsets[0],offsets[1]);
		    printf("readsize:%d %d\n",readsize[0],readsize[1]);
		    printf("start:%d %d\n",start[0],start[1]);
        }
        */
        /*
        if (rank == 0 && ndim = 3) {
            printf("time_flag:%d ndim=%d\n",time_flag,ndim);
            printf("local: %d %d %d\n",ldims[0],ldims[1],ldims[2]);
            printf("global: %d %d %d\n",gdims[0],gdims[1],gdims[2]);
            printf("offset: %d %d %d\n",offsets[0],offsets[1],offsets[2]);
            printf("start:%d %d %d\n",start[0],start[1],start[2]);
            printf("readsize:%d %d %d\n",readsize[0],readsize[1],readsize[2]);
        }
        */
        /* hole_break 
         * -1: the content in the pg is exactly the same as request size 
         * 0 : the content is the pg is the subset of the read request
         * >0: there is hole in content of the pg need to be read
         */
		if (hole_break==-1)
			memcpy(var, fh->b->buff+fh->b->offset,var_header.payload_size);
		else if (hole_break==0) {
        /*
          printf("rank:%d npg:%d tmpreadsize=%llu,"
                    "readsize=%d read_offset=%llu\n",
                    rank, npg, tmpreadsize,readsize[0],read_offset);
        */
			//if (readsize[0]>ldims[0]) {
			if (tmpreadsize > ldims[0]) {
				memcpy ( var+read_offset, fh->b->buff+fh->b->offset, 
                         var_header.payload_size);
				read_offset +=  var_header.payload_size;
                
                tmpreadsize -= ldims[0]; 
                //tmpreadsize = readsize[0]-read_offset/size_of_type; 
                //printf("read_offset=%lld,sizetype=%d, plsize=%lld\n",
                //        read_offset,size_of_type,var_header.payload_size);
			}
			else { 
				memcpy (var+read_offset, 
					fh->b->buff+fh->b->offset+tmpreadsize*datasize*size_of_type, 
					datasize*size_of_type);
			}
		}
		else {
			uint64_t stride_offset = 0;
			int isize;
			uint64_t size_in_dset[10];
			uint64_t offset_in_dset[10];
			uint64_t offset_in_var[10];
			memset(size_in_dset,0,10*8);
			memset(offset_in_dset,0,10*8);
			memset(offset_in_var,0,10*8);
			int hit=0;
			for ( i = 0; i < ndim ; i++) {
				isize = offsets[i] + ldims[i];
				if (start[i] >= offsets[i]) {
					// head is in
					if (start[i]<isize) {
						if (start[i]+readsize[i]>isize)
							size_in_dset[i] = isize - start[i];
						else 
							size_in_dset[i] = readsize[i];	
						offset_in_dset[i] = start[i]-offsets[i];
						offset_in_var[i] = 0;
						hit = 1+hit*10; 
					}
					else
						hit = -1;
				}
				else {
					// middle is in
					if (isize < start[i]+readsize[i]) { 
						size_in_dset[i] = ldims[i];
						hit = 2+hit*10;
					}
					else { 
						// tail is in
						size_in_dset[i] = readsize[i]+start[i]-offsets[i];
						hit = 3+hit*10;
					}
					offset_in_dset[i] = 0;
					offset_in_var[i] = offsets[i]-start[i]; 
				}
			}
			datasize = 1;
			var_stride = 1;	
			for ( i = ndim-1; i >= hole_break; i--) {
				datasize *= size_in_dset[i];
				dset_stride *= ldims[i];
				var_stride *= readsize[i];
			}

			uint64_t var_offset=0;
			uint64_t dset_offset=0;
			for ( i = 0; i < hole_break; i++) { 
				nloop *= size_in_dset[i];
			}
/*
			if(rank==-1) 
			{	
				printf("\n\nhits:%d hole_break=%d nloop=%d\n",hit,hole_break, nloop);	
				if (ndim==2) {
					printf("local: %llu %llu \n",ldims[0],ldims[1]);
					printf("offsets: %llu %llu \n",offsets[0],offsets[1]);
					printf("size  : %llu %llu \n",readsize[0],readsize[1]);
					printf("global: %llu %llu \n",gdims[0],gdims[1]);
				}
				if (ndim==3) {
					printf("local: %llu %llu %llu\n",ldims[0],ldims[1], ldims[2]);
					printf("global: %llu %llu %llu\n",gdims[0],gdims[1], gdims[2]);
					printf("offsets: %llu %llu %llu \n",offsets[0],offsets[1], offsets[2]);
				}
				printf("datasize=%d\n",datasize);
			}

*/
			for ( i = 0; i < ndim ; i++) {
				var_offset = offset_in_var[i] + var_offset * readsize[i];
				dset_offset = offset_in_dset[i] + dset_offset * ldims[i];
			}

			copy_data (var, fh->b->buff+fh->b->offset,
					0,
					hole_break,
					size_in_dset,
					ldims,
					readsize,
					var_stride,
					dset_stride,
					var_offset,
					dset_offset,
					datasize,
					size_of_type);
            //printf("finish read pg: %d\n",npg);
		}
	}  // end of loop
	free (gdims);
	free (offsets);
	free (ldims);
	free (idx_table);
	return 0;
}

void bp_fopen_( int64_t * fh,
		char * fname, 
		MPI_Comm comm,
		int * err,
		int fname_len)
{
	*err = bp_fopen (fh,fname,comm);
}

void bp_fclose_ ( int64_t * fh, int * err)
{
	*err = bp_fclose (*fh);
}

void bp_inq_file_ ( int64_t * fh_p, int * ngroup, 
		   int * nvar, int * nattr, 
		   int * ntime, char ** gnamelist, int * err)
{
	//bp_inq_file ( *fh_p, ngroup, nvar, nattr, ntime, gnamelist);
}

const char * bp_type_to_string (int type)
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
void bp_init_groupinfo(BP_GROUP_INFO * pginfo, int flag)
{
	if (pginfo) {
		memset (pginfo, 0, sizeof(BP_GROUP_INFO));
		pginfo->namelist_true = flag;
	}
	else 
		printf ("fileinfo is NULL\n");
	return;	

}

void bp_free_groupinfo (BP_GROUP_INFO * pginfo)
{
	free_namelist ((pginfo->var_namelist),pginfo->vars_count);
	return;
}

void bp_print_groupinfo (BP_GROUP_INFO *pginfo) 
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

void bp_init_fileinfo(BP_FILE_INFO * pfinfo, int flag)
{
	if (pfinfo) {
		memset (pfinfo, 0, sizeof(BP_FILE_INFO));
		pfinfo->namelist_true = flag;
	}
	else 
		printf ("fileinfo is NULL\n");
	return;	
}

void bp_free_fileinfo (BP_FILE_INFO * pfinfo)
{
	free_namelist ((pfinfo->group_namelist),pfinfo->groups_count);
	return;
}

void bp_print_fileinfo (BP_FILE_INFO *pfinfo) 
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

