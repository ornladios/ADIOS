#ifndef __ARDMA_LOGGER_H__
#define __ARDMA_LOGGER_H__

/* Logger functions */
#include <stdio.h>
#include <string.h>
#include <errno.h>

// These variables must be declared in every implementation of 
// ardma_server.h/ardma_client.h
extern FILE *ardma_logf;  
extern int ardma_verbose_level; 

inline void ardma_logger_init (char *logpath, int verbose_level, int rank)
{   
    if (!logpath || !strcmp(logpath, "stderr")) {
        ardma_logf = stderr;
    } else if (!strcmp(logpath, "stdout")) {
        ardma_logf = stdout;
    } else { 
        char path[256];
        if (rank >= 0)
            snprintf (path, 256, "%s.%d", logpath, rank);
        else
            strncpy (path, logpath, 256);
        ardma_logf = fopen (path, "w");
        if (!ardma_logf) {           
            fprintf (stderr, "Logger file %s cannot be opened. Use stderr for logging.\n"
                             "       errno=%d: %s\n", path, errno, strerror(errno));
            ardma_logf = stderr;     
        }
    }
    ardma_verbose_level = verbose_level;
}   
    
inline void ardma_logger_finalize()
{   
    if (!ardma_logf && ardma_logf != stdout && ardma_logf != stderr)
        fclose(ardma_logf);
}   
    
//#define ardma_logger(verbose_level, ...) if (verbose >= verbose_level) fprintf (stderr, __VA_ARGS__); 
  
#define ardma_logger(verbose_level, ...) {\
    if (ardma_verbose_level >= verbose_level) { \
        fprintf (ardma_logf, __VA_ARGS__); \
        fflush(ardma_logf); \
    }\
}

#define log_error(...) ardma_logger(0, __VA_ARGS__)
#define log_warn(...) ardma_logger(0, __VA_ARGS__)
#define log_info(...) ardma_logger(1, __VA_ARGS__)
#define log_debug(...) ardma_logger(2, __VA_ARGS__)

#endif
