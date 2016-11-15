#include <stdlib.h>
#include <string.h>
#include "public/adios_transform_methods.h"
#include "core/transforms/adios_transforms_hooks.h" 
#include "core/transforms/adios_transforms_hooks_read.h"
#include "core/transforms/adios_transforms_read.h"

ADIOS_IMPLEMENTED_TRANSFORMS * adios_implemented_transforms()
{
    int i, n;
    n = 0;
    for (i = (int)adios_transform_none; i < num_adios_transform_types; i++) {
        if (adios_transform_is_implemented((enum ADIOS_TRANSFORM_TYPE)i)) {
            n++;
        }
    }

    if (n == 0)
        return NULL;

    ADIOS_IMPLEMENTED_TRANSFORMS * t = (ADIOS_IMPLEMENTED_TRANSFORMS *) malloc (sizeof(ADIOS_IMPLEMENTED_TRANSFORMS));
    if (!t)
        return NULL;

    t->name    = (char**) malloc (n*sizeof(char*));
    t->description = (char**) malloc (n*sizeof(char*));
    t->ntransforms = n;

    n = 0;
    for (i = (int)adios_transform_none; i < num_adios_transform_types; i++) {    
        if (adios_transform_is_implemented((enum ADIOS_TRANSFORM_TYPE)i)) {
            t->name[n] = strdup (adios_transform_plugin_primary_xml_alias((enum ADIOS_TRANSFORM_TYPE)i));
            t->description[n] = strdup (adios_transform_plugin_desc((enum ADIOS_TRANSFORM_TYPE)i));
            n++;
        }
    }
    return t;
}

void adios_implemented_transforms_free (ADIOS_IMPLEMENTED_TRANSFORMS *t)
{
	int i;
	if (t)
	{
	    for (i=0; i < t->ntransforms; i++)
	    {
	        if (t->name[i]) {
	            free (t->name[i]);
	            t->name[i] = NULL;
	        }
	        if (t->description[i]) {
	            free (t->description[i]);
	            t->description[i] = NULL;
	        }
	    }
	free (t);
	}
}

