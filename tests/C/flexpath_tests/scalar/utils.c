/**
 * utils.c
 *
 *  Created on: Jul 5, 2013
 *  Author: Magda Slawinska aka Magic Magg magg dot gatech at gmail.com
 *
 *  Utility functions.
 */
#include <stdio.h>
#include "misc.h"

/**
 * @param program_name The name of the program
 * @param program_desc The description of the program
 * @return 0 Returns always 0
 */
int usage(char *program_name, char *program_desc){
	printf("Usage: mpirun -np nprocs %s\n", program_name);
	printf("  %s\n", program_desc);
	printf("  METHOD: ");
#ifdef FLEXPATH_METHOD
	printf("FLEXPATH/ADIOS_READ_METHOD_FLEXPATH" );
#else
	printf("MPI/ADIOS_READ_METHOD_BP");
#endif
	printf("  will be used.\n");
	printf("  Please check " XML_ADIOS_INIT_FILENAME "\n");
	printf("  nprocs Number of processors you want to use\n");
	printf("E.g. mpirun -np 5 %s\n", program_name);

	return 0;
}
