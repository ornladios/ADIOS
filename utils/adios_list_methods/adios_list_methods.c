/*
 * in the COPYING file in the top level directory of this source distribution.
 *
 * Copyright (c) 2008 - 2009.  UT-BATTELLE, LLC. All rights reserved.
 */

/* ADIOS list_methods utility
 *  list available 
 *    write methods
 *    read methods
 *    transform methods
 *
 * This is a sequential program but compiled with MPI to use libadios.a
 */

#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <inttypes.h>
#include <unistd.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <getopt.h>
#include <errno.h>

#include <math.h>     // NAN
#include <libgen.h>   // basename
#include <regex.h>    // regular expression matching
#include <fnmatch.h>  // shell pattern matching

#ifdef WRITE
#include "public/adios.h"
#endif
#include "public/adios_read.h"
#include "public/adios_transform_methods.h"
#include "public/adios_query.h"


int main (int argc, char ** argv) {
    int  rank = 0, i;
#ifndef _NOMPI  // added only to enable compiling with Score-P and other libs substituting MPI 
    MPI_Comm comm = MPI_COMM_WORLD;
    MPI_Init(&argc,&argv);
    MPI_Comm_rank(comm,&rank);
#endif


#ifdef WRITE
    adios_init_noxml(MPI_COMM_WORLD);
#endif

#ifndef _NOMPI
    adios_read_init_method(ADIOS_READ_METHOD_BP, MPI_COMM_WORLD, "");
#else
    adios_read_init_method(ADIOS_READ_METHOD_BP, 1, "");
#endif

    if(rank==0) {

#ifdef WRITE
    	printf ("Available write methods (in XML <method> element or in adios_select_method()):\n");
    	ADIOS_AVAILABLE_WRITE_METHODS * wm = adios_available_write_methods();
    	if (wm) {
    		for (i = 0; i < wm->nmethods; i++) {
    			printf("    \"%s\"\n", wm->name[i]);
    		}
    		adios_available_write_methods_free(wm);
    	}
#endif

        printf ("Available read methods (constants after #include \"adios_read.h\"):\n");
        ADIOS_AVAILABLE_READ_METHODS * rm = adios_available_read_methods();
        if (rm) {
        	for (i = 0; i < rm->nmethods; i++) {
        		printf("    \"%s\"\n", rm->name[i]);
        	}
        	adios_available_read_methods_free(rm);
        }

        printf ("Available data transformation methods (in XML transform tags in <var> elements):\n");
        ADIOS_AVAILABLE_TRANSFORM_METHODS * t = adios_available_transform_methods();
        if (t) {
            for (i=0; i<t->ntransforms; i++)
            {
                printf("    \"%s\"\t: %s\n",  t->name[i], t->description[i]);
            }
            adios_available_transform_methods_free(t);
        }


        printf ("Available query methods (in adios_query_set_method()):\n");
        ADIOS_AVAILABLE_QUERY_METHODS * qm = adios_available_query_methods();
        if (qm) {
        	for (i = 0; i < qm->nmethods; i++) {
        		printf("    \"%s\"\n", qm->name[i]);
        	}
        	adios_available_query_methods_free(qm);
        }
    }

#ifndef _NOMPI  // added only to enable compiling with Score-P and other libs substituting MPI 
    MPI_Barrier(comm);
    MPI_Finalize();
#endif
    return(0);
}
