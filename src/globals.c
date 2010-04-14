/* 
 * ADIOS is freely available under the terms of the BSD license described
 * in the COPYING file in the top level directory of this source distribution.
 *
 * Copyright (c) 2008 - 2009.  UT-BATTELLE, LLC. All rights reserved.
 */

/*
 *   Functions, constants globally for both the Write and Read API
 */

#include "globals.h"

static int globals_adios_appid = -1;
static int globals_adios_was_set = 0;
void globals_adios_set_application_id (int id)
{
    globals_adios_appid = id;
    globals_adios_was_set = 1;
}

int globals_adios_get_application_id (int *was_set) 
{
    *was_set = globals_adios_was_set;
    return globals_adios_appid;
}

#ifdef HAVE_DART
enum DART_CONNECTION { dart_disconnected = 0, 
                       dart_connected_from_reader = 1,
                       dart_connected_from_writer = 2,
                       dart_connected_from_both = 3
                     };
static enum DART_CONNECTION globals_adios_connected_to_dart = dart_disconnected;

void globals_adios_set_dart_connected_from_reader()
{ 
    if (globals_adios_connected_to_dart == dart_disconnected)
        globals_adios_connected_to_dart = dart_connected_from_reader;
    else if (globals_adios_connected_to_dart == dart_connected_from_writer)
        globals_adios_connected_to_dart = dart_connected_from_both;
}
void globals_adios_set_dart_disconnected_from_reader()
{ 
    if (globals_adios_connected_to_dart == dart_connected_from_reader)
        globals_adios_connected_to_dart = dart_disconnected;
    else if (globals_adios_connected_to_dart == dart_connected_from_both)
        globals_adios_connected_to_dart = dart_connected_from_writer;
}
void globals_adios_set_dart_connected_from_writer()
{
    if (globals_adios_connected_to_dart == dart_disconnected)
        globals_adios_connected_to_dart = dart_connected_from_writer;
    else if (globals_adios_connected_to_dart == dart_connected_from_reader)
        globals_adios_connected_to_dart = dart_connected_from_both;
}
void globals_adios_set_dart_disconnected_from_writer()
{ 
    if (globals_adios_connected_to_dart == dart_connected_from_writer)
        globals_adios_connected_to_dart = dart_disconnected;
    else if (globals_adios_connected_to_dart == dart_connected_from_both)
        globals_adios_connected_to_dart = dart_connected_from_reader;
}
int  globals_adios_is_dart_connected()
{ 
    return (globals_adios_connected_to_dart != dart_disconnected);
}
int  globals_adios_is_dart_connected_from_reader()
{ 
    return (globals_adios_connected_to_dart == dart_connected_from_reader || 
            globals_adios_connected_to_dart == dart_connected_from_both);
}
int  globals_adios_is_dart_connected_from_writer()
{ 
    return (globals_adios_connected_to_dart == dart_connected_from_writer || 
            globals_adios_connected_to_dart == dart_connected_from_both);
}
int  globals_adios_is_dart_connected_from_both()
{
    return (globals_adios_connected_to_dart == dart_connected_from_both);
}
#endif
