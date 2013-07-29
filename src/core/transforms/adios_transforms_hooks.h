/*
 * adios_transforms_hooks.h
 *
 *  Created on: Feb 14, 2013
 *      Author: David A. Boyuka II
 */

#ifndef ADIOS_TRANSFORMS_HOOKS_H_
#define ADIOS_TRANSFORMS_HOOKS_H_

// Build the ADIOS_TRANSFORM_TYPE enum
#include "core/transforms/plugindetect/detect_plugin_types.h"

// Transform info functions

const char * adios_transform_plugin_uuid(enum ADIOS_TRANSFORM_TYPE transform_type);
const char * adios_transform_plugin_desc(enum ADIOS_TRANSFORM_TYPE transform_type);
int adios_transform_plugin_num_xml_aliases(enum ADIOS_TRANSFORM_TYPE transform_type);
const char ** adios_transform_plugin_xml_aliases(enum ADIOS_TRANSFORM_TYPE transform_type);
const char * adios_transform_plugin_primary_xml_alias(enum ADIOS_TRANSFORM_TYPE transform_type);

// Transform XML alias <-> ID conversion

/*
 * @param xml_alias the name of a transform type as specified in the ADIOS XML
 * @return the ADIOS_TRANSFORM_TYPE corresponding to that alias, or
 *         adios_transform_unknown if it does not match any registered
 *         transform type
 */
enum ADIOS_TRANSFORM_TYPE adios_transform_find_type_by_xml_alias(const char *xml_alias);

/*
 * @param transform_type a transform type
 * @return a name that can be used to refer to that transform type within
 *         the ADIOS XML file, or "" if transform_type is not valid.
 */
const char * adios_transform_xml_alias_by_type(enum ADIOS_TRANSFORM_TYPE transform_type);

// Other inspection functions

int is_transform_type_valid(enum ADIOS_TRANSFORM_TYPE transform_type);

#endif /* ADIOS_TRANSFORMS_HOOKS_H_ */
