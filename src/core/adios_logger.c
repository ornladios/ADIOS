/* 
 * ADIOS is freely available under the terms of the BSD license described
 * in the COPYING file in the top level directory of this source distribution.
 *
 * Copyright (c) 2008 - 2009.  UT-BATTELLE, LLC. All rights reserved.
 */

/* Logger functions */
#include <stdio.h>
#include <string.h>
#include <errno.h>
#include "core/adios_logger.h"

FILE *adios_logf = 0;
int adios_verbose_level = 1; // ERROR level (0 = no logs)
char *adios_log_names[4] = {"ERROR","WARN","INFO","DBG"};


void adios_logger_init (char *logpath, int verbose_level, int rank)
{ 
    if (!logpath || !strcmp(logpath, "stderr")) { 
        adios_logf = stderr; 
    } else if (!strcmp(logpath, "stdout")) { 
        adios_logf = stdout; 
    } else {  
        char path[256]; 
        if (rank >= 0) 
            snprintf (path, 256, "%s.%d", logpath, rank); 
        else 
            strncpy (path, logpath, 256); 
        adios_logf = fopen (path, "w"); 
        if (!adios_logf) {  
            fprintf (stderr, "Logger file %s cannot be opened. Use stderr for logging.\n" 
                             "       errno=%d: %s\n", path, errno, strerror(errno)); 
            adios_logf = stderr; 
        } 
    } 
    adios_verbose_level = verbose_level; 
}   
    
void adios_logger_finalize() 
{ 
    if (!adios_logf && adios_logf != stdout && adios_logf != stderr) 
        fclose(adios_logf); 
}
    
/*
void adios_logger(int verbose_level, ...) 
{
    if (adios_verbose_level >= verbose_level) { 
        fprintf (adios_logf, __VA_ARGS__); 
        fflush(adios_logf); 
    }
}
*/
