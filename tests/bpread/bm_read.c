#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include "mpi.h"
#include "bp_read.h"

int main (int argc, char ** argv)
{
	char * filename;
	int    rank, pe_size, rc, version;	
	struct bp_index_pg_struct_v1 * pg_root = 0;
	struct bp_index_pg_struct_v1 * pg = 0;

	struct  adios_index_var_struct_v1 * vars_root = 0;
	struct  adios_index_attribute_struct_v1 * attrs_root = 0;
	int64_t fh=0, gh=0;
	int     ngroup, nvar, nattr, ntstep;
	int     i,j,k;	

	int	NX = 10, G_NX=30;
	int	NY = 10, G_NY=30;
	int	NZ = 10, G_NZ=30;
	int  	dims [10];
	void    * var = NULL;
	MPI_Comm comm; 

	filename = argv[1];
	comm = MPI_COMM_WORLD;
	MPI_Init (&argc, &argv);
	MPI_Comm_rank (comm, &rank);
	MPI_Comm_size (comm, &pe_size);

	bp_fopen (&fh, filename, comm);
	// fh = file handle, filename = name of file, comm is the 
	// communicator for the processors which will read in the data.
	// in this call, the header is actually being read, and this is done by proc 0, and
	// then it is broadcast to all of the processors.
	BP_FILE_INFO finfo;
	/*
	   BP_FILE_INFO's memmbers:
		uint16_t namelist_true;
		uint16_t groups_count;
		uint16_t vars_count;
		uint16_t attrs_count;
		uint32_t tidx_start;
		uint32_t tidx_stop;
		uint32_t version;
		uint32_t file_size;
		char     ** group_namelist;
	 */
	bp_init_fileinfo ( &finfo, 1);
	// initialize the FILE_INFO structure
	//       IN: 
	// 	    struct BP_FILE_INFO * 
 	// 	    int flag: 1 means allocate the momory for group_namelist
        //		      0 means no allocation 
	// * bp_free_fileinfo (&finfo) need to be paired with this function call
	bp_inq_file (fh, &finfo);
	
	if (rank == 0) 
		bp_print_fileinfo (&finfo);

	BP_GROUP_INFO ginfo;
	bp_gopen (&gh, fh, finfo.group_namelist[0]);
	bp_init_groupinfo (&ginfo,1);
	bp_inq_group (gh, &ginfo);
	if (rank == 0) 
		bp_print_groupinfo (&ginfo);

	int type, ndim, time_flag;
	bp_inq_var (gh, ginfo.var_namelist[7], &type, &ndim, &time_flag, dims);
	char *type_str = bp_type_to_string (type);
        printf("%s:\n\t\tis_timebased: %d\n\t\ttype: %s\n\t\tdimensions:",
		ginfo.var_namelist[7], time_flag, type_str);
	for (i=0;i<ndim;i++) {
		printf(" [%d]", dims[i]); 
	}
	printf("\n");
	// nowe we are getting the type, the number of dimensions, and the dimensions
	// dims is stored as the ordering of how it was written.
	
	int start[3], size[3];

	start[0]=0;
	start[1]=0;
	size[0]=3;
	size[1]=3;
	var = malloc (sizeof(double) * size[0]*size[1]);

	// time step should be no less than zero
	// vnamelist[14] is not written at time step 0
	// so the function returns as error
	bp_get_var (gh, ginfo.var_namelist[7], var, start, size, 4);
	
	printf("data:\n");
	k=0;
	for (j=0;j<size[0];j++) {
		for (i=0;i<size[1];i++){
			printf("%f\t",*(double*)(var+8*k));
			k++;
		}
		printf("\n");
	}
	printf("\n");
 	if (var)
		free(var);

	bp_free_groupinfo (&ginfo);
	bp_free_fileinfo (&finfo);

	bp_gclose (gh);
	bp_fclose (fh);
	MPI_Finalize ();

	return 0;
}

