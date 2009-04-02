#include <stdio.h>
#include <sys/types.h>
#include "mpi.h"
#include "bp_types.h"
#define VARS_MINIHEADER_SIZE 10
#define DIVIDER "========================================================\n"
#define SUBDIVIDER "------------------------------------\n"

//static void alloc_aligned (struct adios_bp_buffer_struct_v1 * b, uint64_t size);
//static void realloc_aligned (struct adios_bp_buffer_struct_v1 * b, uint64_t size);
int bp_parse_characteristics (struct adios_bp_buffer_struct_v1 * b,
		  	      struct adios_index_var_struct_v1 ** root,
			      uint64_t j);
int bp_get_characteristics_data (void ** ptr_data,
				 void * buffer,
				 int data_size,
				 enum ADIOS_DATATYPES type);
int bp_read_close (struct adios_bp_buffer_struct_v1 * b);
int bp_read_open (const char * filename,
	 	  MPI_Comm comm, 
		  struct BP_FILE * fh);

const char * value_to_string (enum ADIOS_DATATYPES type, void * data);
void bp_grouping ( struct BP_FILE * fh, uint64_t * gh);
uint64_t bp_get_type_size (enum ADIOS_DATATYPES type, void * var);

void print_process_group_index (
                         struct adios_index_process_group_struct_v1 * pg_root
                         );

