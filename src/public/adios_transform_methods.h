#ifndef ADIOS_TRANSFORM_METHODS_H
#define ADIOS_TRANSFORM_METHODS_H

#ifdef __cplusplus
extern "C" {
#endif

typedef struct
{
    int ntransforms;            // number of transformations
    char ** name;               // array of transform names
    char ** description;        // array of descriptions
} ADIOS_AVAILABLE_TRANSFORM_METHODS;

/* Provide the names of transformations available in the running application
 */
ADIOS_AVAILABLE_TRANSFORM_METHODS * adios_available_transform_methods();

void adios_available_transform_methods_free (ADIOS_AVAILABLE_TRANSFORM_METHODS *);

#ifdef __cplusplus
}
#endif

#endif
