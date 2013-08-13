/*
 * plugin_info_types.h
 *
 *  Created on: Feb 15, 2013
 *      Author: David A. Boyuka II
 */

#ifndef PLUGIN_INFO_TYPES_H_
#define PLUGIN_INFO_TYPES_H_

typedef struct {
    enum ADIOS_TRANSFORM_TYPE type;
    const char *uid;
    const char *description;
} adios_transform_plugin_info_t;

typedef struct {
    enum ADIOS_TRANSFORM_TYPE type;

    int xmlAliasCount;
    const char **xmlAliases;
} adios_transform_method_xml_aliases_t;

#endif /* PLUGIN_INFO_TYPES_H_ */
