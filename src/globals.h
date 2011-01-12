/* 
 * ADIOS is freely available under the terms of the BSD license described
 * in the COPYING file in the top level directory of this source distribution.
 *
 * Copyright (c) 2008 - 2009.  UT-BATTELLE, LLC. All rights reserved.
 */

/*
 *   Functions, constants globally for both the Write and Read API
 */
#ifndef __GLOBALS_H__
#define __GLOBALS_H__

#include "../config.h"

/** Set an application ID for this program. 
 *  This function is necessary for methods who needs a unique ID from each participating applications.
 *  Currently, this is the DART method for code coupling of independent applications.
 *
 *  This function is called from the applicatin through adios_set_application_id()
 */
void globals_adios_set_application_id (int id);


/** Get the application ID set by the application itself.
  * Returns the ID set by the application, -1 if it was not set.
  * It also returns a boolean was_set, 0 if it was not set, 1 otherwise.
  */
int globals_adios_get_application_id (int *was_set);


/* Note: would be nice a <string, int> map for arbitrary globals */
#ifdef HAVE_DART
void globals_adios_set_dart_connected_from_reader();
void globals_adios_set_dart_disconnected_from_reader();
void globals_adios_set_dart_connected_from_writer();
void globals_adios_set_dart_disconnected_from_writer();
int  globals_adios_is_dart_connected(); // from any
int  globals_adios_is_dart_connected_from_reader();
int  globals_adios_is_dart_connected_from_writer();
int  globals_adios_is_dart_connected_from_both();
#endif /* HAVE_DART */

#ifdef HAVE_DIMES
void globals_adios_set_dimes_connected_from_reader();
void globals_adios_set_dimes_disconnected_from_reader();
void globals_adios_set_dimes_connected_from_writer();
void globals_adios_set_dimes_disconnected_from_writer();
int  globals_adios_is_dimes_connected(); // from any
int  globals_adios_is_dimes_connected_from_reader();
int  globals_adios_is_dimes_connected_from_writer();
int  globals_adios_is_dimes_connected_from_both();
#endif

#if NO_DATATAP == 0

char *getFixedName(char *name);
char *get_full_path_name(char *name, char *path);

#endif

#endif  /*__GLOBALS_H__*/
