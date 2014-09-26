/*
 * generate_indexed_dataset.c
 *
 *  Created on: Sep 26, 2014
 *      Author: David A. Boyuka II
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <assert.h>

#include <mpi.h>
#include <adios.h>

typedef enum { DATASET_1 } DATASET_ID;

typedef struct {
	const char *group_name;
	unsigned int buffer_size_mb;
	const char *write_transport_method;
	int ndim;
	const char **varnames;
	const char **types;
} dataset_xml_spec_t;

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

extern void adios_pin_timestep(uint32_t ts); // Not in the standard header, but accessible
void build_dataset(const char *filename_prefix, const dataset_xml_spec_t *xml_spec, const char *transform_name, int num_timesteps, int num_pgs_per_timestep, const uint64_t *groupsizes, const void **vardatas) {
	int timestep, pg_in_timestep;
	char xml_filename[strlen(filename_prefix) + strlen(".xml") + 1];
	char bp_filename[strlen(filename_prefix) + strlen(".bp") + 1];

	// Construct the XML and BP filenames
	sprintf(xml_filename, "%s.xml", filename_prefix);
	sprintf(bp_filename, "%s.bp", filename_prefix);

	// Write out the XML file
	FILE *xml_out = fopen(xml_filename, "w");
	assert(xml_out);
	produce_xml(xml_out, xml_spec, transform_name);
	fclose(xml_out);

	// Write out the BP file
	adios_init(xml_filename, MPI_COMM_WORLD);

	// Compute the groupsize contribution of the dimension scalars
	const base_groupsize = xml_spec->ndim * 3 * 4; // *3 for 3 scalars (N, D, O) *4 for sizeof(adios_integer) (not sure how what function in the User API to call to get this programatically

	// For each timestep, for each PG in that timestep, write out all variables using the provided vardata buffers
	int64_t adios_file;
	for (timestep = 0; timestep < num_timesteps; ++timestep) {
		for (pg_in_timestep = 0; pg_in_timestep < num_pgs_per_timestep; ++pg_in_timestep) {
			// (Re-)open the file in write or append mode, depending on whether or not this is the first PG written
			const int is_first_pg = (timestep == 0 && pg_in_timestep == 0);
			adios_open(&adios_file, xml_spec->group_name, bp_filename, is_first_pg ? "w" : "a", MPI_COMM_WORLD);

			// Pin the timestep to allow multiple adios_open/adios_close cycles to write
			// to the same timestep (this simulates a parallel file write with fewer core)
			adios_pin_timestep(timestep);

			// Compute the group size
			uint64_t out_groupsize;
			adios_group_size(adios_file, base_groupsize + *groupsizes, &out_groupsize);

			// Write each variable
			const char **varnames = xml_spec->varnames;
			while (*varnames != NULL) {
				adios_write(adios_file, *varnames, (void*)*vardatas); // (void*) to get rid of compiler complaining about constness
				++varnames;
				++vardatas;
			}

			// Close the file to commit it
			adios_close(adios_file);

			++groupsizes;
		}
	}
}

void build_dataset_1(const char *filename_prefix, const char *transform_name) {
	static const char *VARNAMES[] = { "temp", NULL };
	static const char *VARTYPES[] = { "float", NULL };

	static const dataset_xml_spec_t XML_SPEC = {
		.group_name = "S3D",
		.buffer_size_mb = 128,
		.write_transport_method = "MPI",
		.ndim = 2,
		.varnames = VARNAMES,
		.types = VARTYPES,
	};

	static const float TEMP_DATA[] = {
		//  1.00000000     1.00003052     2.00000000     2.00006104
			0x1.000000p+0, 0x1.000200p+0, 0x1.000000p+1, 0x1.000200p+1,
		//  2.00012207     30.00000000    30.00048828    30.00097656
			0x1.000400p+1, 0x1.e00000p+4, 0x1.e00200p+4, 0x1.e00400p+4,
		//  30.00146484    50.00000000    50.00097656    50.00195312
			0x1.e00600p+4, 0x1.900000p+5, 0x1.900200p+5, 0x1.900400p+5,
		//  50.00292969    50.00390625    50.00488281    50.00585938
			0x1.900600p+5, 0x1.900800p+5, 0x1.900a00p+5, 0x1.900c00p+5,
	};

	// DS1 is simple: it has just 1 PG
	static const int NUM_TIMESTEPS = 1;
	static const int NUM_PGS_PER_TIMESTEP = 1;
	static const uint64_t PG_GROUPSIZES[] = {
		sizeof(TEMP_DATA),
	};
	static const void *PG_DATAS[] = {
		TEMP_DATA,
	};

	build_dataset(filename_prefix, &XML_SPEC, transform_name, NUM_TIMESTEPS, NUM_PGS_PER_TIMESTEP, PG_GROUPSIZES, PG_DATAS);
}

void usage_and_exit() {
	fprintf(stderr, "Usage: build_indexed_dataset <dataset-id> <filename-prefix> [<transform-type>]\n");
	fprintf(stderr, "\n");
	fprintf(stderr, "  dataset-id: the ID of a standard dataset packaged in this executable:\n");
	fprintf(stderr, "    - 'DS1': 1 variable, 2D array of floats\n");
	fprintf(stderr, "\n");
	fprintf(stderr, "  filename-prefix: the filename prefix for the XML and BP files to be generated.\n");
	fprintf(stderr, "    - Example: filename-prefix = 'some/path/myfile':\n");
	fprintf(stderr, "            -> produces some/path/myfile.xml, some/path/myfile.bp:\n");
	fprintf(stderr, "\n");
	fprintf(stderr, "  transform-type: the data transform to apply (default: none)\n");
	exit(1);
}

int main(int argc, char **argv) {
	if (argc < 3 || argc > 4)
		usage_and_exit();

	const char *dataset_id = argv[1];
	const char *path = argv[2];
	const char *transform_name = (argc >= 4) ? argv[3] : "none";

	DATASET_ID dataset;
	// Select the dataset by dataset ID
	if (strcasecmp(dataset_id, "DS1") == 0) {
		dataset = DATASET_1;
	} else {
		fprintf(stderr, "Error: '%s' does not name a dataset packaged in this executable\n");
		usage_and_exit();
	}

	MPI_Init(&argc, &argv);

	switch (dataset) {
	case DATASET_1:
		build_dataset_1(path, transform_name);
		break;
	}

	MPI_Finalize();
	return 0;
}
