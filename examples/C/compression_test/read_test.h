
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "mpi.h"
#include "adios_read.h"

#define VAR_NAME	"/var/temperature"
#define MAX_DIM	8

typedef int (*test_func)(const char*, const char*, const uint64_t*, const uint64_t*,  const uint64_t*, const uint64_t, const int, double**, uint64_t*);

int g_ndims = 0;
uint64_t g_dim[MAX_DIM] = {0};
uint64_t g_total = 1;

int get_dim_info(char** all_file_names, char* var_name)
{
	char* file_name = all_file_names[0];
	enum ADIOS_READ_METHOD  method = ADIOS_READ_METHOD_BP;
    MPI_Comm comm = MPI_COMM_WORLD;
	
	ADIOS_FILE *f = adios_read_open_file(file_name, method, comm);
    ADIOS_VARINFO *varinfo = adios_inq_var(f, var_name);
	
	g_ndims = varinfo->ndim;
	int d = 0;
	for(d = 0; d < g_ndims && d < MAX_DIM; d++)
	{
		g_dim[d] = varinfo->dims[d];
		g_total *= varinfo->dims[d];
	}
	
	adios_free_varinfo(varinfo);
    adios_read_close(f);
    adios_read_finalize_method(ADIOS_READ_METHOD_BP);
	return 0;
}

int generate_random_box(uint64_t* starts, uint64_t* counts)
{
	srand(time(NULL));
	int d = 0;
	for(d = 0; d < g_ndims && d < MAX_DIM; d++)
	{
		starts[d] = rand() % (g_dim[d] - 1);
		counts[d] = rand() % (g_dim[d] - 1 - starts[d]) + 1;
		printf("sel dim[%d]: %d %d\n", d, starts[d], counts[d]);
	}
	return 0;
}

int generate_random_points(uint64_t** points, uint64_t* npoints)
{
	srand(time(NULL));
	*npoints = rand() % g_total + 1;	
	*points  = (uint64_t *)malloc((*npoints) * g_ndims * sizeof (uint64_t));
	
	printf("num of points %d %d\n", *npoints, g_total);
	
	int i = 0;
	for (i = 0; i < *npoints; i ++)
	{
		int d = 0;
		for (d = 0; d < g_ndims; d++)
		{
			(*points)[i * g_ndims + d] = rand() % g_dim[d];
		}
	}
	
	return 0;
}

int generate_random_writeblock(int* writeblock_index)
{
	*writeblock_index = rand() % g_dim[0];
	return 0;
}

int test_all_files(char** all_file_names, int file_count, test_func func)
{
	uint64_t starts[MAX_DIM] = {0};
	uint64_t counts[MAX_DIM] = {0};
	generate_random_box(starts, counts);
		
	uint64_t npoints = 0;
	uint64_t* points = NULL;
	generate_random_points(&points, &npoints);
	
	int writeblock_index = 0;
	generate_random_writeblock(&writeblock_index);
	
	double** buffers = (double**)malloc(file_count * sizeof(double*));
	uint64_t* nresults_all = (uint64_t*)malloc(file_count * sizeof(uint64_t));
	
	int f = 0;
	for(f = 0; f < file_count; f++)
	{
		double* data;
		uint64_t nresults;
		int rtn = func(all_file_names[f], VAR_NAME, 
						starts, counts, 
						points, npoints,
						writeblock_index,						
						&data, &nresults);
		if(rtn != 0)
		{
			printf("test failed\n");
			return -1;
		}
		
		buffers[f] = data;
		nresults_all[f] = nresults;
		// printf("%s read in\n", all_file_names[f]);
	}
	
	uint64_t npoints_test = nresults_all[0];
	for(f = 1; f < file_count; f++)
	{
		if(nresults_all[f] != npoints_test)
		{
			printf("test failed\n");
			return -1;
		}
	}
	
	for(f = 1; f < file_count; f++)
	{
		// int i = 0;
		// for(i = 0; i < npoints_test; i++)
		// {
			// printf("%f ", buffers[0][i]);
		// }
		// printf("\n");
		// for(i = 0; i < npoints_test; i++)
		// {
			// printf("%f ", buffers[f][i]);
		// }
		// printf("\n");
		
		if(memcmp(buffers[0], buffers[f], npoints_test * sizeof(double)) != 0)
		{
			printf("test failed\n");
			return -1;
		}
	}
	
	printf("test succ\n");
	return 0;
}



