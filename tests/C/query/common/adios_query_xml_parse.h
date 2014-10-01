/*
 * adios_query_xml_parse.h
 *
 *  Created on: Sep 30, 2014
 *      Author: Houjun Tang
 */
#ifndef ADIOS_QUERY_XML_PARSE_H_
#define ADIOS_QUERY_XML_PARSE_H_

#include "adios_query.h"

/*
 * Parses the given XML file to extract the ADIOS_QUERY object described therein.
 */
ADIOS_QUERY * parseXml(const char *inputxml, ADIOS_FILE* f);

#endif /* ADIOS_QUERY_XML_PARSE_H_ */
