#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include "mpi.h"
#include "adios_read.h"

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
        if (!filename)
            filename = "testbp_c.bp";
	comm = MPI_COMM_WORLD;
	MPI_Init (&argc, &argv);
	MPI_Comm_rank (comm, &rank);
	MPI_Comm_size (comm, &pe_size);
	adios_fopen (&fh, filename, comm);
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
	adios_init_fileinfo ( &finfo, 1);
	// initialize the FILE_INFO structure
	//       IN: 
	// 	    struct BP_FILE_INFO * 
 	// 	    int flag: 1 means allocate the momory for group_namelist
        //		      0 means no allocation 
	// * adios_free_fileinfo (&finfo) need to be paired with this function call
	adios_inq_file (fh, &finfo);
	
	if (rank == 0) 
		adios_print_fileinfo (&finfo);

	BP_GROUP_INFO ginfo;
	adios_gopen (fh, &gh, finfo.group_namelist[0]);
	adios_init_groupinfo (&ginfo,1);
	adios_inq_group (gh, &ginfo);
	if (rank == 0) 
		adios_print_groupinfo (&ginfo);

	int type, ndim, time_flag;
	int start[3], size[3];
	char *type_str;
    double time;
    // get scalar data
	adios_inq_var (gh, "int_1D", &type, &ndim, &time_flag, dims);
	type_str = bp_type_to_string (type);
    printf("%s:\n\t\tis_timebased: %d\n\t\ttype: %s\n",
		    "int_1D", time_flag, type_str);
    if (ndim==0) 
        printf("\t\tscalar");
    else 
            printf("\t\tdimension: [%d]", dims[0]);
    printf ("\n");
    size[0] = 1;
    start[0] = 0;
	adios_get_var (gh, "int_1D", &nvar, start, size, 1);
    printf("\t\tvalue: %d\n",(nvar));
    // get 2D data
	adios_inq_var (gh, "int_2D", &type, &ndim, &time_flag, dims);
    printf("%s:\n\t\tis_timebased: %d\n\t\ttype: %s\n\t\tdimensions:",
    		"int_2D", time_flag, type_str);
    if (ndim==0) {
        printf("\t\tscalar\n");
    }
    else {
        for (i=0;i<ndim;i++) {
            printf(" [%d]", dims[i]); 
        }
	}
	printf("\n");
	start[0]=0;
	start[1]=0;
	size[0]=10;
	size[1]=2;
	var = malloc (sizeof(int) * size[0]*size[1]);

	// time step should be no less than zero
	// vnamelist[14] is not written at time step 0
	// so the function returns as error
	adios_get_var (gh, "int_2D", var, start, size, 1);
	k=0;
    printf("\t\t[%d:%d, %d:%d]", 
           start[0], size[0],
           start[1], size[1]
           );
    for (j=0;j<size[0];j++) {
        printf("\n\t\t");
        for (i=0;i<size[1];i++){
            printf("%d  ",*(int*)(var+4*k));
            k++;
        }
    }
	printf("\n");

	adios_inq_var (gh, "int_3D", &type, &ndim, &time_flag, dims);
    printf("%s:\n\t\tis_timebased: %d\n\t\ttype: %s\n\t\tdimensions:",
    		"int_3D", time_flag, type_str);
    if (ndim==0) {
        printf("\t\tscalar\n");
    }
    else {
        for (i=0;i<ndim;i++) {
            printf(" [%d]", dims[i]); 
        }
	    printf("\n");
	}
    
    if (var)
        free (var);
	var = malloc ( sizeof(int) * dims[0]*dims[1]*dims[2]);
	start[0]=0;
	start[1]=0;
	start[2]=0;
	size[0]=2;
	size[1]=5;
	size[2]=3;
	adios_get_var (gh, "int_3D", var, start, size, 1);
	k=0;
    printf("\t\t[%d:%d, %d:%d, %d:%d]", 
           start[0], size[0],
           start[1], size[1],
           start[2], size[2]
           );
    for (nvar=0;nvar<size[0];nvar++) {
        for (j=0;j<size[1];j++) {
            printf("\n\t\t");
            for (i=0;i<size[2];i++){
                printf("%d  ",*(int*)(var+4*k));
                k++;
            }
        }
        printf("\n");
    }

    printf("\ndone\n");

    if (var)
		free(var);
        
	adios_free_groupinfo (&ginfo);
	adios_free_fileinfo (&finfo);

	adios_gclose (gh);
	adios_fclose (fh);
	MPI_Finalize ();

	return 0;
}

