/* 
 * ADIOS is freely available under the terms of the BSD license described
 * in the COPYING file in the top level directory of this source distribution.
 *
 * Copyright (c) 2008 - 2009.  UT-BATTELLE, LLC. All rights reserved.
 */

#include <stdlib.h>
#include <math.h>
#include <string.h>

// mpi
//#include "mpi.h"

// xml parser
#include <mxml.h>

// dart
#include <sys/uio.h>
#include "adios.h"
#include "adios_transport_hooks.h"
#include "adios_internals.h"
#include "adios_internals_mxml.h"

#define STR_LEN 1000

//*****************************************
// Check the groups as follows:
//*****************************************
void check_groups (struct adios_group_list_struct * groups)
{
// no extra checks needed
}

//*****************************************
// Check the groups as follows:
// 1. Make sure all method groups match real groups
//*****************************************
void check_methods (struct adios_method_list_struct * methods
                   ,struct adios_group_list_struct * groups
                   )
{
// no extra checks needed
}

int main (int argc, char ** argv)
{
    struct adios_method_list_struct * methods = 0;
    struct adios_group_list_struct * groups = 0;
    char * filename;
    MPI_Comm comm = 0; //dummy comm

    if (argc < 2)
        filename = "config.xml";
    else
        filename = argv [1];

    if (!adios_parse_config (filename, comm))
        return 1;

    methods = adios_get_methods ();
    groups = adios_get_groups ();

    check_groups (groups);
    check_methods (methods, groups);

    adios_cleanup ();

    return 0;
}
