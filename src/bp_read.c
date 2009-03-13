//#include "bp_types.h"
#include <stdlib.h>
#include "bp_utils.h"
#include "bp_read.h"
#include "adios_bp_v1.h"
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
	        char * fname,
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
		bp_read_minifooter (fh->b, &(fh->mfooter));
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
	bp_gh->fh = fh_p;
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

int bp_inq_file ( int64_t fh_p, int *ngroup, 
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
	*nt = fh->mfooter.time_steps;
	*nattr = 0;
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
	int i,k;
	gh->var_current = 0;
	var_root = fh->vars_root;

	for(i=0;i<gh->offset;i++)
		var_root = var_root->next;
	for (i=0;i<gh->count;i++) {
		if (!strcmp(varname, var_root->var_name)) {
	
			gh->var_current = var_root;
			*type = var_root->type;
			*is_timebased = 0;
			*ndim = var_root->characteristics [0].dims.count;
			if (!dims || !(*ndim))
				return;
			int time_flag = -1;
			uint64_t * gdims = (uint64_t *) malloc (sizeof(uint64_t) * (*ndim));
			uint64_t * ldims = (uint64_t *) malloc (sizeof(uint64_t) * (*ndim));
			memset(dims,0,sizeof(int)*(*ndim));
			//for (j=0;j<var_root->characteristics_count;j++)  
			for (k=0;k<(*ndim);k++) {
				gdims[k]=var_root->characteristics[0].dims.dims[k*3+1];
				ldims[k]=var_root->characteristics[0].dims.dims[k*3];
			}
			int is_global=0;
			for (k=0;i<(*ndim);i++) {
			 	is_global = is_global || gdims[k];
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
				for (i=0;i<(*ndim);i++) {
					if (gdims[i]==0) {
						time_flag = i;
						*is_timebased = 1;
					}
					if (dims) {
						if (time_flag==0) {
							if (i>0)
								dims[i-1]=gdims[i];
						}
						else
							dims[i]=gdims[i];
					}
				}
			}
					
			free(gdims);
			free(ldims);
			if ((*is_timebased))
				*ndim = *ndim - 1;
			
			break;
		}
		else
			var_root = var_root->next;
	}
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
	
	int i, j, var_id;
	struct BP_GROUP * gh = (struct BP_GROUP *) gh_p;
	struct BP_FILE * fh = (struct BP_FILE *) (gh->fh);
	struct adios_index_var_struct_v1 * var_root = fh->vars_root;

	for (i = 0; i < gh->count; i++) {
		if (!strcmp(fh->gh->var_namelist[i+gh->offset], varname)) {
			var_id = i+gh->offset;
			break;
		}
	}
	for (i=0;i<var_id;i++) 
		var_root = var_root->next;
	gh->var_current = var_root;

	int offset, count, start_idx;
	if (timestep < 1) {
		fprintf(stderr, "time index should start from 1!");
		return; 
	}
	offset = fh->gh->time_index[0][gh->group_id][timestep-1];
        count = fh->gh->time_index[1][gh->group_id][timestep-1];
	for (i=0;i<var_root->characteristics_count;i++) {
		if (var_root->characteristics[i].offset > fh->gh->pg_offsets[offset]) {
			start_idx = i;
			break;
		}
	}
	//printf("time step: %d %d %d\n", offset, count, start_idx);
        struct adios_var_header_struct_v1 var_header;
        struct adios_var_payload_struct_v1 var_payload;
	uint64_t size;	
	MPI_Status status;
 	int idx;
	uint8_t ndim;  
	uint64_t * ldims, * offsets, * gdims;
	int *idx_start, *idx_stop;
	ndim = var_root->characteristics[start_idx].dims.count;
	ldims = (uint64_t *) malloc (sizeof(uint64_t) * ndim);
	offsets = (uint64_t *) malloc (sizeof(uint64_t) * ndim);
	gdims = (uint64_t *) malloc (sizeof(uint64_t) * ndim);
        int * idx_table = (int *) malloc (sizeof(int)*count);
	
	int idx_count = 0;
  	int flag;
	uint64_t datasize, nloop, dset_stride,var_stride, total_size=1;

	// main for loop
	int rank;
	MPI_Comm_rank(gh->fh->comm, &rank);
	int time_flag = -1;
	for (i=0;i<ndim;i++) {
		gdims[i]=var_root->characteristics[0].dims.dims[i*3+1];
		offsets[i]=var_root->characteristics[0].dims.dims[i*3+2];
		ldims[i]=var_root->characteristics[0].dims.dims[i*3];
	}
	int is_global=0, is_timebased = 0;
	for (i=0;i<ndim;i++) {
		is_global = is_global || gdims[i];
	}
	if (!is_global) {
		for (i=0;i<ndim;i++) {
			if (   ldims[i] == 1
			    && var_root->characteristics_count > 1) {
				is_timebased = 1;
				time_flag = i;
			}
		}
	}
	else {
		for (i=0;i<ndim;i++) {
			if (gdims[i]==0) {
				time_flag = i;
				is_timebased = 1;
			}
		}
	}
	
	if (is_timebased)
		ndim = ndim -1;
		
	for (idx = 0; idx < count; idx++) {
		idx_table[idx] = 1;
/* 
		for (i=0;i<ndim;i++) { 
			gdims[i]=var_root->characteristics[start_idx+idx].dims.dims[i*3+1];
			if (gdims[i] == 0)
				break;
		}
		if (i != ndim) { 
			ndim = ndim - 1;
		}
*/	
		for (j=0;j<ndim;j++) {
			if (time_flag == 0) {
				ldims[j]=var_root->characteristics[start_idx+idx].dims.dims[j*3+3];
				if (is_global) 
					gdims[j]=var_root->characteristics[start_idx+idx].dims.dims[j*3+1];
				else
					gdims[j]=ldims[j];
				offsets[j]=var_root->characteristics[start_idx+idx].dims.dims[j*3+2];
			}
			else if (time_flag == ndim) {
				ldims[j]=var_root->characteristics[start_idx+idx].dims.dims[j*3];
				if (is_global) 
					gdims[j]=var_root->characteristics[start_idx+idx].dims.dims[j*3+1];
				else
					gdims[j]=ldims[j];
				offsets[j]=var_root->characteristics[start_idx+idx].dims.dims[j*3+2];
			}
				
			flag = (offsets[j] >= start[j] 
					&& offsets[j] <start[j]+readsize[j])
				|| (offsets[j] < start[j]
						&& offsets[j]+ldims[j]>start[j]+readsize[j]) 
				|| (offsets[j]+ldims[j] > start[j] 
						&& offsets[j]+ldims[j] <= start[j]+readsize[j]);
			idx_table [idx] = idx_table[idx] && flag;
		}
	}
//	for (i=0;i<count;i++)
//		printf("%d %d\n",i,idx_table[i]);
	uint64_t read_offset = 0;
	int npg=0;
	for (idx = 0; idx < count; idx++) {
		datasize = 1;
		nloop = 1;
		var_stride = 1;
		dset_stride = 1;
		if ( !idx_table[idx] ) {
			continue;
		}
		//printf("start to read npg=%d !\n", npg);
		++npg;
		MPI_File_seek (fh->mpi_fh, 
			       (MPI_Offset) var_root->characteristics[start_idx+idx].offset,
			       MPI_SEEK_SET);
		MPI_File_read (fh->mpi_fh, fh->b->buff, 4, MPI_BYTE, &status);
		int tmpcount = 0;
		MPI_Get_count (&status, MPI_BYTE, &tmpcount);
		tmpcount= *((uint64_t*)fh->b->buff);	
		realloc_aligned(fh->b, tmpcount+4);

		MPI_File_seek (fh->mpi_fh, 
			       (MPI_Offset) var_root->characteristics[start_idx+idx].offset,
			       MPI_SEEK_SET);
		MPI_File_read (fh->mpi_fh, fh->b->buff, tmpcount+4, MPI_BYTE, &status);
		MPI_Get_count (&status, MPI_BYTE, &tmpcount);
		fh->b->offset = 0;
		adios_parse_var_data_header_v1 (fh->b, &var_header);

		//print_var_header (&var_header);
		int size_of_type = bp_get_type_size (var_root->type, "");

		struct   adios_dimension_struct_v1 * d = var_header.dims;
		d = var_header.dims;
/*
		int time_flag = -1;
		while (d) {
			if (d->dimension.time_index==adios_flag_yes) {
				++time_flag;
				break;
			}
			else {
				d=d->next;
			}
		}

		if (time_flag == -1)
			ndim = var_root->characteristics[start_idx].dims.count;
		else
			ndim = var_root->characteristics[start_idx].dims.count-1;
		//d = var_header.dims;
*/

		if (time_flag == 0) { 
			for (j=0;j<ndim;j++) {
				offsets[j]=var_root->characteristics[start_idx+idx].dims.dims[j*3+2];
				ldims[j]=var_root->characteristics[start_idx+idx].dims.dims[j*3+3];
				if (is_global)
					gdims[j]=var_root->characteristics[start_idx+idx].dims.dims[j*3+1];
				else
					gdims[j]=ldims[j];
			}
		}
		else {
			for (j=0;j<ndim;j++) {
				offsets[j]=var_root->characteristics[start_idx+idx].dims.dims[j*3+2];
				ldims[j]=var_root->characteristics[start_idx+idx].dims.dims[j*3];
				if (is_global)
					gdims[j]=var_root->characteristics[start_idx+idx].dims.dims[j*3+1];
				else
					gdims[j]=ldims[j];
			}
		}

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
/*
		printf("time_flag:%d ndim=%d\n",time_flag,ndim);
		printf("local: %d %d\n",ldims[0],ldims[1]);
		printf("global: %d %d\n",gdims[0],gdims[1]);
		printf("offset: %d %d\n",offsets[0],offsets[1]);
		printf("readsize:%d %d\n",readsize[0],readsize[1]);
		printf("start:%d %d\n",start[0],start[1]);
*/
		hole_break = i;	

		if (hole_break==-1)
			memcpy(var, fh->b->buff+fh->b->offset,var_header.payload_size);
		else if (hole_break==0) {
			if (readsize[0]>ldims[0]) {
				memcpy(var+read_offset, fh->b->buff+fh->b->offset, var_header.payload_size);
				read_offset +=  var_header.payload_size;	
			}
			else { 
				printf("datasize:%d \n",datasize);
				memcpy (var+read_offset, 
					fh->b->buff+fh->b->offset+start[0]*datasize*size_of_type, 
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
						//printf("head!\n");
						if (start[i]+readsize[i]>isize)
							size_in_dset[i] = isize - start[i];
						else 
							size_in_dset[i] = readsize[i];	
						offset_in_dset[i] = start[i]-offsets[i];
						offset_in_var[i] = 0;
						hit = 1+hit*10; 
						//offset_in_var[i] = start[i]; 
					}
					else
						hit = -1;
				}
				else {
					// middle is in
					if (isize < start[i]+readsize[i]) { 
						size_in_dset[i] = ldims[i];
						hit = 2+hit*10;
						//printf("middle!\n");
					}
					else { 
						//printf("tail!\n");
						size_in_dset[i] = readsize[i]+start[i]-offsets[i];
						hit = 3+hit*10;
						// tail is in
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
/*
			if(rank==-1) {
				//for(i=0;i<ndim;i++){
				//	printf("%llu %llu %llu\n", ldims[i], gdims[i], offsets[i]);	
//					printf("rank %d: size_in_dset=%d, offset_in_dset=%d, offset_in_var=%d\n",
//							rank, size_in_dset[i], offset_in_dset[i], offset_in_var[i]);
					//printf("rank %d: start=%d,offsets=%d isize=%d ldims=%d\n",
					//		rank, start[i], offsets[i],isize, ldims[i]);
//
				//}
				printf("rank %d: var_offset=%llu, var_stride=%llu, " 
					"dset_offset=%llu dset_stride=%llu\n",
					rank,var_offset,var_stride,dset_offset,dset_stride);
			}
*/
		}
/*
		if(rank==-1) 
			printf("how many blocks:%d \n",npg);
*/
	}  // end of loop

	return 0;
}

void bp_fread_( int64_t * fh,
		const char * fname, 
		MPI_Comm comm,
		int * err,
		int fname_len)
{
	*err = bp_fread (fh,fname,comm);
}

void bp_fclose_ ( int64_t * fh, int * err)
{
	*err = bp_fclose (*fh);
}

void bp_inq_file_ ( int64_t * fh_p, int * ngroup, 
		   int * nvar, int * nattr, 
		   int * ntime, char ** gnamelist, int * err)
{
	bp_inq_file ( fh_p, ngroup, nvar, nattr, ntime, gnamelist);
}
