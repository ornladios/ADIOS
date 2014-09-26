/*
 * generate_indexed_dataset.c
 *
 *  Created on: Sep 26, 2014
 *      Author: David A. Boyuka II
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct {
	const char *group_name;
	unsigned int buffer_size_mb;
	const char *write_transport_method;
	int ndim;
	const char **varnames;
	const char **types;
} dataset_xml_spec_t;

dataset_xml_spec_t DATASET_1 = {
	"S3D",
	128,
	"MPI",
	2,
	(const char*[]){ "temp", NULL },
	(const char*[]){ "float", NULL },
};

static void build_dimension_var_list(int ndim, const char *dimvar_base, char *outbuf) {
	int i;
	for (i = 0; i < ndim; i++)
		outbuf += sprintf(outbuf, i == 0 ? "%s%d" : ",%s%d", dimvar_base, i);
}

static void produce_xml(
		FILE *outfile,
		const dataset_xml_spec_t *xml_spec,
		const char *transform_name)
{
	static const char *HEADER_XML =
	"<?xml version=\"1.0\"?>\n"
	"<adios-config host-language=\"C\">\n"
	"	<adios-group name=\"%s\" coordination-communicator=\"comm\">\n";

	static const char *DIMVAR_XML =
	"		<var name=\"N%d\" type=\"integer\"/>\n"
	"		<var name=\"D%d\" type=\"integer\"/>\n"
	"		<var name=\"O%d\" type=\"integer\"/>\n";

	static const char *GLOBALBOUNDS_HEADER_XML =
	"		<global-bounds dimensions=\"%s\" offsets=\"%s\">\n";

	static const char *VAR_XML =
	"			<var name=\"%s\" type=\"%s\" dimensions=\"%s\" transform=\"%s\"/>\n";

	static const char *GLOBALBOUNDS_FOOTER_XML =
	"		</global-bounds>\n";

	static const char *FOOTER_XML =
	"	</adios-group>\n"
	"	<method group=\"%s\" method=\"%s\"/>\n"
	"	<buffer size-MB=\"%u\" allocate-time=\"now\"/>\n"
	"</adios-config>\n";

	// Actual function begins
	int i;
	char dimvar_list_buf1[256]; // for storing D0,D1,D2,...,Dn
	char dimvar_list_buf2[256]; // for storing D0,D1,D2,...,Dn

	fprintf(outfile, HEADER_XML, xml_spec->group_name);
	for (i = 0; i < xml_spec->ndim; i++)
		fprintf(outfile, DIMVAR_XML, i, i, i);

	build_dimension_var_list(xml_spec->ndim, "N", dimvar_list_buf1);
	build_dimension_var_list(xml_spec->ndim, "O", dimvar_list_buf2);
	fprintf(outfile, GLOBALBOUNDS_HEADER_XML, dimvar_list_buf1, dimvar_list_buf2);

	const char **varnames = xml_spec->varnames;
	const char **types = xml_spec->types;
	while (*varnames) {
		build_dimension_var_list(xml_spec->ndim, "D", dimvar_list_buf1);
		fprintf(outfile, VAR_XML, *varnames, *types, dimvar_list_buf1, transform_name);

		++varnames;
		++types;
	}

	fprintf(outfile, GLOBALBOUNDS_FOOTER_XML);
	fprintf(outfile, FOOTER_XML, xml_spec->group_name, xml_spec->write_transport_method, xml_spec->buffer_size_mb);
}

void usage_and_exit() {
	fprintf(stderr, "Usage: build_indexed_dataset <dataset-id> [<transform-type>]\n");
	fprintf(stderr, "\n");
	fprintf(stderr, "  dataset-id: the ID of a standard dataset packaged in this executable:\n");
	fprintf(stderr, "    - 'DS1': 1 variable, 2D array of floats\n");
	fprintf(stderr, "\n");
	fprintf(stderr, "  transform-type: the data transform to apply (default: none)\n");
	exit(1);
}

int main(int argc, char **argv) {
	if (argc < 2 || argc > 3)
		usage_and_exit();

	const char *dataset_id = argv[1];
	const char *transform_name = (argc >= 3) ? argv[2] : "none";

	// Select the dataset by dataset ID
	dataset_xml_spec_t *xml_spec = NULL;
	if (strcasecmp(dataset_id, "DS1") == 0) {
		xml_spec = &DATASET_1;
	} else {
		fprintf(stderr, "Error: '%s' does not name a dataset packaged in this executable\n");
		usage_and_exit();
	}

	produce_xml(stdout, xml_spec, transform_name);
}
