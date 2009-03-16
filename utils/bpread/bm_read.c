
#include <stdio.h>
#include <sys/types.h>
#include "mpi.h"
#include "bp_read.h"
#include "bp_types.h"
#define VARS_MINIHEADER_SIZE 10
#define DIVIDER "========================================================\n"
#define SUBDIVIDER "------------------------------------\n"
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

void print_namelist (char **namelist, int length)
{
	int i;
	for (i=0; i<length; i++)
		printf("\t%d: \t%s\n", i, namelist[i]);
	return;
}

int main (int argc, char ** argv)
{
	char * filename;
	int    rank, pe_size, rc, version;	
	struct bp_index_pg_struct_v1 * pg_root = 0;
	struct bp_index_pg_struct_v1 * pg = 0;

	struct  adios_index_var_struct_v1 * vars_root = 0;
	struct  adios_index_attribute_struct_v1 * attrs_root = 0;
	int64_t fh, gh;
	int     ngroup, nvar, nattr, ntstep;
	int     i,j,k;	

	int	NX = 10, G_NX=30;
	int	NY = 10, G_NY=30;
	int	NZ = 10, G_NZ=30;
	int  	dims [10];
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
	bp_inq_file (fh, &ngroup, &nvar, &nattr, &ntstep, 0);
	// this API gets the number of groups, variables, timesteps, and attributes from
	// the file
	// all parameters are integers, besides the last parameter, which is an array of 
	// strings
	if (rank == 0) {
                printf(DIVIDER);
                printf("bp_inq_file():\n"); 
                printf ("\t# of groups: %d\n", ngroup);
                printf ("\t# of variables: %d\n", nvar);
                printf ("\t# of attributes: %d\n", nattr);
                printf ("\t# of time steps: %d\n", ntstep);
        }

	char ** gnamelist = NULL;
	char ** vnamelist = NULL;

	alloc_namelist (&gnamelist, ngroup);

	bp_inq_file (fh, &ngroup, &nvar, &nattr, &ntstep, gnamelist);
	// if the last parameter is not zero, then this routine will return the list of 
	// group names into gnamelist (the last parameter)
	if (rank == 0) {
		printf(SUBDIVIDER);
		printf("%d groups in the group!\n", ngroup);
		print_namelist (gnamelist, ngroup);
	}

	bp_gopen (&gh, fh, gnamelist[0]);
	// this opens up the specific group, given by the name (last parameter).
	// returns the second parameter.
	bp_inq_group (gh, &nvar, 0);
	// this gets the number of variables inside the group
	// if we give the last parameter as an array of strings, then we can return this
	// which is the name of variables in the group.
	alloc_namelist (&vnamelist, nvar);
	bp_inq_group (gh, &nvar, vnamelist);
	// this rotuine is now returning the name of the variables.

	if (rank == 0) {
		printf(SUBDIVIDER);
		printf("%d vars in the group!\n", nvar);
		print_namelist (vnamelist, nvar);
	}

	int type, ndim, time_flag;
			
	bp_inq_var (gh, vnamelist[14], &type, &ndim, &time_flag, dims);
	
	char *type_str = bp_type_to_string (type);
        printf("%s:\n"
		"\t\tis_timebased: %d\n"
		"\t\ttype: %s\n\t\tdimensions:",
		vnamelist[14], time_flag, type_str);
	for (i=0;i<ndim;i++) {
		printf(" [%d]", dims[i]); 
	}
	printf("\n");

	// nowe we are getting the type, the number of dimensions, and the dimensions
	// dims is stored as the ordering of how it was written.
	
	int start[3], size[3];

	start[0]=0;
	size[0]=50;

	void * var = NULL;
	var = malloc (sizeof(double) * size[0]);

	// time step should be no less than zero
	// vnamelist[14] is not written at time step 0
	// so the function returns as error
	bp_get_var (gh, vnamelist[14], var, start, size, 0);
	
	bp_get_var (gh, vnamelist[14], var, start, size, 1);

	printf("data:\n");
	for (i=0;i<size[0];i++)
		printf("%lf\t",*(double*)(var+8*i));
	printf("\n");

	if (var)
		free(var);
	if (vnamelist)
		free_namelist (vnamelist, nvar);
	if (gnamelist)
		free_namelist (gnamelist, ngroup);

	bp_gclose (gh);


	bp_fclose (fh);

	MPI_Finalize ();


	return 0;
}

