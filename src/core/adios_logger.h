#ifndef __ADIOS_LOGGER_H__
#define __ADIOS_LOGGER_H__

/* Logger functions */
#include <stdio.h>
#include <string.h>

extern FILE *adios_logf;
extern int adios_verbose_level;
extern char *adios_log_names[4];

void adios_logger_init (char *logpath, int verbose_level, int rank);
void adios_logger_finalize();

#define  adios_logger(verbose_level, ...) { \
    if (adios_verbose_level >= verbose_level) { \
        if (!adios_logf) adios_logf=stderr; \
        fprintf (adios_logf, "ADIOS %s: ",adios_log_names[verbose_level]); \
        fprintf (adios_logf, __VA_ARGS__); \
        fflush(adios_logf);\
    }\
}


#define log_error(...) adios_logger(0, __VA_ARGS__)
#define log_warn(...) adios_logger(1, __VA_ARGS__)
#define log_info(...) adios_logger(2, __VA_ARGS__)
#define log_debug(...) adios_logger(3, __VA_ARGS__)

#endif
