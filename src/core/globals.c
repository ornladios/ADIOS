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

#ifdef HAVE_DIMES
enum DIMES_CONNECTION { dimes_disconnected = 0,
                        dimes_connected_from_reader = 1,
                        dimes_connected_from_writer = 2,
                        dimes_connected_from_both = 3
                     };
static enum DIMES_CONNECTION globals_adios_connected_to_dimes = dimes_disconnected;

void globals_adios_set_dimes_connected_from_reader()
{
    if (globals_adios_connected_to_dimes == dimes_disconnected)
        globals_adios_connected_to_dimes = dimes_connected_from_reader;
    else if (globals_adios_connected_to_dimes == dimes_connected_from_writer)
        globals_adios_connected_to_dimes = dimes_connected_from_both;
}
void globals_adios_set_dimes_disconnected_from_reader()
{
    if (globals_adios_connected_to_dimes == dimes_connected_from_reader)
        globals_adios_connected_to_dimes = dimes_disconnected;
    else if (globals_adios_connected_to_dimes == dimes_connected_from_both)
        globals_adios_connected_to_dimes = dimes_connected_from_writer;
}
void globals_adios_set_dimes_connected_from_writer()
{
    if (globals_adios_connected_to_dimes == dimes_disconnected)
        globals_adios_connected_to_dimes = dimes_connected_from_writer;
    else if (globals_adios_connected_to_dimes == dimes_connected_from_reader)
        globals_adios_connected_to_dimes = dimes_connected_from_both;
}
void globals_adios_set_dimes_disconnected_from_writer()
{
    if (globals_adios_connected_to_dimes == dimes_connected_from_writer)
        globals_adios_connected_to_dimes = dimes_disconnected;
    else if (globals_adios_connected_to_dimes == dimes_connected_from_both)
        globals_adios_connected_to_dimes = dimes_connected_from_reader;
}
int  globals_adios_is_dimes_connected()
{
    return (globals_adios_connected_to_dimes != dimes_disconnected);
}
int  globals_adios_is_dimes_connected_from_reader()
{
    return (globals_adios_connected_to_dimes == dimes_connected_from_reader ||
            globals_adios_connected_to_dimes == dimes_connected_from_both);
}
int  globals_adios_is_dimes_connected_from_writer()
{
    return (globals_adios_connected_to_dimes == dimes_connected_from_writer ||
            globals_adios_connected_to_dimes == dimes_connected_from_both);
}
int  globals_adios_is_dimes_connected_from_both()
{
    return (globals_adios_connected_to_dimes == dimes_connected_from_both);
}
#endif

#if NO_DATATAP == 0

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#define OPLEN 2
static char OP[OPLEN] = { '(', ')' };
static char *OP_REP[OPLEN] = { "_PPLT_", "_PPRT_" };

char *getFixedName(char *name)
{
    char *tempname = (char *) malloc(sizeof(char) * 255);
    snprintf(tempname, 255, "%s\0", name);
    char *oldname = strdup(name);
    char *loc = NULL;
    int i;

    do
    {
        for (i = 0; i < OPLEN; i++)
        {
    //checking operator OP[i]
            loc = strchr(oldname, OP[i]);
            if (loc == NULL)
                continue;
            *loc = 0;
            snprintf(tempname, 255, "%s%s%s\0", oldname, OP_REP[i], &loc[1]);
            free(oldname);
            oldname = strdup(tempname);
        }
    }
    while (loc != NULL);
    free(oldname);

fprintf(stderr, "im here %s %s %s:%d\n", name, tempname, __FILE__,__LINE__);
    return tempname;
}

char *get_full_path_name(char *name, char *path)
{
    char *full_pathname = (char *) malloc(strlen(name)+strlen(path)+2);
    if(!full_pathname) {
        fprintf(stderr, "cannot allocate memory. %s:%d\n", __FILE__,__LINE__);
        return NULL;
    }
    if (!strcmp (path, "/")) {
        sprintf (full_pathname, "/%s\0", name);
    }
    else {
        sprintf (full_pathname, "%s/%s\0", path, name);
    }
fprintf(stderr, "im here %s %s %s %s:%d\n", name, path, full_pathname,__FILE__,__LINE__);
    return full_pathname;
}

#endif

