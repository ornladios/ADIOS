#ifndef ADIOS_INTERNALS_MXML_H
#define ADIOS_INTERNALS_MXML_H

int adios_parse_config (const char * config);
int adios_local_config ();
int adios_common_select_method (int priority, const char * method
                               ,const char * parameters, const char * group 
                               ,const char * base_path, int iters
                               );
void adios_cleanup ();

#endif
