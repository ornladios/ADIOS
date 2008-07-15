#ifndef BR_UTILS_H
#define BR_UTILS_H

#include "binpack-general.h"
#include "binpack-utils.h"

typedef void (* ADIOS_BR_PRE_FETCH) (struct adios_bp_element_struct * element, void ** buffer, uint64_t * buffer_size, void * private_data);
typedef void (* ADIOS_BR_POST_FETCH) (struct adios_bp_element_struct * element, void * buffer, uint64_t buffer_size, void * private_data);

long long br_fopen (char * filename);
void br_fclose (long long handle);

long long br_bopen (void * buffer, uint64_t buf_len);
void br_bclose (long long handle);

void br_free_element (struct adios_bp_element_struct * element);

int br_get_next_element_specific (long long handle, ADIOS_BR_PRE_FETCH pre, ADIOS_BR_POST_FETCH post, void * private_data , struct adios_bp_element_struct ** element);

int br_get_next_element_general (long long handle, void * buffer, uint64_t buffer_size, struct adios_bp_element_struct ** element);

#endif
