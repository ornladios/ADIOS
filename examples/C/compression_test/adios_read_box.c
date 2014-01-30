
#include "read_test.h"

int read_bounding_box(const char* file_name, const char* var_name, 
						const uint64_t* starts, const uint64_t* counts,
						const uint64_t* points, const uint64_t npoints,
						int writeblock_index,
						double** data, uint64_t* nreslts)
{
	enum ADIOS_READ_METHOD  method = ADIOS_READ_METHOD_BP;
    MPI_Comm comm = MPI_COMM_WORLD;
	
	ADIOS_FILE *f = adios_read_open_file(file_name, method, comm);
	// adios_group_view(f,-1);
    ADIOS_VARINFO *varinfo = adios_inq_var(f, var_name);
	
	if(varinfo->ndim != g_ndims)
	{
		printf("test failed\n");
		adios_free_varinfo(varinfo);
		adios_read_close(f);
		adios_read_finalize_method(ADIOS_READ_METHOD_BP);
		return -1;
	}

    assert (varinfo);
    adios_inq_var_blockinfo (f, varinfo);
	
	ADIOS_SELECTION* sel1 = adios_selection_boundingbox(varinfo->ndim, starts, counts);
	
	*nreslts = 1;
	int i = 0;
	for (i = 0; i < varinfo->ndim; i++) 
	{
        (*nreslts) *= counts[i];
    }
	
	*data = (double*)malloc((*nreslts) * sizeof (double));
	
	adios_schedule_read(f, sel1, var_name, 0, 1, *data);
	adios_perform_reads(f, 1);
	
	adios_selection_delete(sel1);
	adios_free_varinfo(varinfo);
    adios_read_close(f);
    adios_read_finalize_method(ADIOS_READ_METHOD_BP);
	
	return 0;
}

int main(int argc, char** argv)
{
	MPI_Init (&argc, &argv);
	
	// const char* all_file_names[2] = {"output/mpi_none.bp", "output/mpi_zlib.bp"};
	int file_count = argc - 1;
	char** all_file_names = argv + 1;
	
	get_dim_info(all_file_names, VAR_NAME);
	
	printf("[%d]: ", g_ndims);
	int i = 0;
	for(i = 0; i < g_ndims; i++)
	{
		printf("%lld ", g_dim[i]);
	}
	printf("\n");
	
	// for(i = 0; i < file_count; i++)
	// {
		// printf("%s\n", all_file_names[i]);
	// }
	
	int rtn = test_all_files(all_file_names, file_count, read_bounding_box);
	
    MPI_Finalize ();
	
	return rtn;
}
