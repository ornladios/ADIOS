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
} ADIOS_AVAILABLE_TRANSFORMS;

/* Provide the names of transformations available in the running application
 */
ADIOS_AVAILABLE_TRANSFORMS * adios_available_transforms();

void adios_available_transforms_free (ADIOS_AVAILABLE_TRANSFORMS *);

#ifdef __cplusplus
}
#endif

#endif
