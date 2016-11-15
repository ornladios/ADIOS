#ifndef ADIOS_TRANSFORM_METHODS_H
#define ADIOS_TRANSFORM_METHODS_H

#ifdef __cplusplus
extern "C" {
#endif

typedef struct
{
    int ntransforms; 		// number of transformations
    char ** name;     		// array of transform names
    char ** description;	// array of descriptions
} ADIOS_IMPLEMENTED_TRANSFORMS;

/* Return an array of transformations available in the running application
 * The return value is the size of the array.
 */
ADIOS_IMPLEMENTED_TRANSFORMS * adios_implemented_transforms();

void adios_implemented_transforms_free (ADIOS_IMPLEMENTED_TRANSFORMS *t);

#ifdef __cplusplus
}
#endif

#endif
