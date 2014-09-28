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
#include <adios_types.h>

typedef enum { DATASET_1, DATASET_2 } DATASET_ID;

typedef struct {
	const char *group_name;
	unsigned int buffer_size_mb;
	const char *write_transport_method;
	int ndim;
	int nvar;
	const char **varnames;
	const enum ADIOS_DATATYPES *vartypes;
} dataset_xml_spec_t;

// Imported from adios_internal.c
extern const char * adios_type_to_string_int(int type);                    // converts enum ADIOS_DATATYPES to string
extern uint64_t adios_get_type_size(enum ADIOS_DATATYPES type, void *var); // returns the size in bytes of a given enum ADIOS_DATATYPES

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

	for (i = 0; i < xml_spec->nvar; ++i) {
		build_dimension_var_list(xml_spec->ndim, "D", dimvar_list_buf1);
		fprintf(outfile, VAR_XML, xml_spec->varnames[i], adios_type_to_string_int(xml_spec->vartypes[i]), dimvar_list_buf1, transform_name);
	}

	fprintf(outfile, GLOBALBOUNDS_FOOTER_XML);
	fprintf(outfile, FOOTER_XML, xml_spec->group_name, xml_spec->write_transport_method, xml_spec->buffer_size_mb);
}

static void write_adios_dimension_scalars(int64_t fd, const char *dimvar_basename, int ndim, const uint64_t *dims) {
	int i;
	char dimvar_name[32];
	for (i = 0; i < ndim; ++i) {
		sprintf(dimvar_name, "%s%d", dimvar_basename, i);
		adios_write(fd, dimvar_name, (void*)dims);
		++dims;
	}
}

typedef struct {
	int num_ts;
	int num_pgs_per_ts;
	const uint64_t *global_dims;
} dataset_global_spec_t;

typedef struct {
	const uint64_t *pg_dim;
	const uint64_t *pg_offset;
	const void **vardata;
} dataset_pg_spec_t;

static uint64_t dim_prod(int ndim, const uint64_t *dims) {
	uint64_t pg_gridsize = 1;
	while (ndim--)
		pg_gridsize *= *dims++;
	return pg_gridsize;
}

static uint64_t compute_groupsize(uint64_t base_groupsize, const dataset_xml_spec_t *xml_spec, const dataset_pg_spec_t *pg) {
	int var, dim;

	// Compute the number of points contained in this PG
	uint64_t pg_gridsize = dim_prod(xml_spec->ndim, pg->pg_dim);

	// Compute the sum of the datatype sizes across all variables defined in this PG
	uint64_t total_var_datatypes_size = 0;
	for (var = 0; var < xml_spec->nvar; ++var)
		total_var_datatypes_size += adios_get_type_size(xml_spec->vartypes[var], NULL);

	// The final group size is the product of the number of points and the number of bytes per point, plus the base groupsize
	return base_groupsize + pg_gridsize * total_var_datatypes_size;
}

extern void adios_pin_timestep(uint32_t ts); // Not in the standard header, but accessible
void build_dataset_from_specs(
		const char *filename_prefix,
		const char *transform_name,
		const dataset_xml_spec_t *xml_spec,
		const dataset_global_spec_t *global_spec,
		/*const*/ dataset_pg_spec_t pg_specs[global_spec->num_ts][global_spec->num_pgs_per_ts]) // Not const because C has an corner case here (http://c-faq.com/ansi/constmismatch.html)
{
	int var;
	char xml_filename[strlen(filename_prefix) + strlen(".xml") + 1];
	char bp_filename[strlen(filename_prefix) + strlen(".bp") + 1];
	int timestep, pg_in_timestep;
	char dimvar[32];

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
	for (timestep = 0; timestep < global_spec->num_ts; ++timestep) {
		for (pg_in_timestep = 0; pg_in_timestep < global_spec->num_pgs_per_ts; ++pg_in_timestep) {
			// (Re-)open the file in write or append mode, depending on whether or not this is the first PG written
			const int is_first_pg = (timestep == 0 && pg_in_timestep == 0);
			adios_open(&adios_file, xml_spec->group_name, bp_filename, is_first_pg ? "w" : "a", MPI_COMM_WORLD);

			// Pin the timestep to allow multiple adios_open/adios_close cycles to write
			// to the same timestep (this simulates a parallel file write with fewer core)
			adios_pin_timestep(timestep + 1); // +1 because we want the timesteps to be 1-based

			const dataset_pg_spec_t *pg_spec = &pg_specs[timestep][pg_in_timestep];

			// Compute the group size
			uint64_t groupsize = compute_groupsize(base_groupsize, xml_spec, pg_spec);
			uint64_t out_groupsize;
			adios_group_size(adios_file, groupsize, &out_groupsize);

			write_adios_dimension_scalars(adios_file, "N", xml_spec->ndim, global_spec->global_dims);
			write_adios_dimension_scalars(adios_file, "D", xml_spec->ndim, pg_spec->pg_dim);
			write_adios_dimension_scalars(adios_file, "O", xml_spec->ndim, pg_spec->pg_offset);

			// Write each variable
			for (var = 0; var < xml_spec->nvar; ++var) {
				adios_write(adios_file, xml_spec->varnames[var], (void*)pg_spec->vardata[var]); // (void*) to get rid of compiler complaining about constness
			}

			// Close the file to commit it
			adios_close(adios_file);
		}
	}
}

// NOTE: varblocks_by_var is actually a 1D array varblocks_by_var[var] of 3D "arrays" varblockdata[ts][pg][point_in_pg] of type xml_spec->vartypes[var] (however, points per PG varies)
void collect_varblocks_by_pg(
		const dataset_xml_spec_t *xml_spec,
		const dataset_global_spec_t *global_spec,
		const uint64_t pg_dims[global_spec->num_ts][global_spec->num_pgs_per_ts][xml_spec->ndim],
		const void **varblocks_by_var,
		const void *out_varblocks_by_pg[global_spec->num_ts][global_spec->num_pgs_per_ts][xml_spec->nvar])
{
	const int num_ts = global_spec->num_ts;
	const int num_pgs_per_ts = global_spec->num_pgs_per_ts;
	const int num_vars = xml_spec->nvar;
	const int ndim = xml_spec->ndim;
	int ts, pg, var;

	// Some maths
	const uint64_t varblocks_per_ts = num_vars * num_pgs_per_ts;

	// Cache the datatype size for each variable, and create a data pointer that we can advance
	int var_typesizes[num_vars];
	const char *varblock_datas[num_vars];
	for (var = 0; var < num_vars; ++var) {
		var_typesizes[var] = adios_get_type_size(xml_spec->vartypes[var], NULL);
		varblock_datas[var] = (const char *)varblocks_by_var[var];
	}

	// Iterate over all varblocks (var in pg in timestep), get a pointer to the data buffer for that
	// varblock, and assign it to the proper place in the output array of PG data buffers
	for (ts = 0; ts < num_ts; ++ts) {
		for (pg = 0; pg < num_pgs_per_ts; ++pg) {
			// Compute the points per varblock in this PG
			const uint64_t points_per_varblock = dim_prod(ndim, pg_dims[ts][pg]);

			for (var = 0; var < num_vars; ++var) {
				const int var_typesize = var_typesizes[var];
				const uint64_t varblock_size = points_per_varblock * var_typesize;

				// Get the data for this varblock, and advance the pointer by one varblock
				const char *data_for_varblock = varblock_datas[var];
				varblock_datas[var] += varblock_size;

				// Finally, assign the pointer
				out_varblocks_by_pg[ts][pg][var] = data_for_varblock;
			}
		}
	}
}

// NOTE: pg_dims and pg_offsets are really a 2D arrays pd_dims[pg][dim] of dimension lengths
// NOTE: pg_datas is really a 2D array pg_datas[pg][var] of varblock buffers
// NOTE: pg_specs is really a 1D array of pg_specs[pg] of dataset_pg_spec_t's
void collect_pg_specs(
		int num_ts, int num_pgs_per_ts, int ndim, int nvar,
		const uint64_t pg_dims[num_ts][num_pgs_per_ts][ndim],
		const uint64_t pg_offsets[num_ts][num_pgs_per_ts][ndim],
		const void *pg_datas[num_ts][num_pgs_per_ts][nvar],
		dataset_pg_spec_t pg_specs[num_ts][num_pgs_per_ts]) {
	int ts, pg;
	for (ts = 0; ts < num_ts; ++ts) {
		for (pg = 0; pg < num_pgs_per_ts; ++pg) {
			pg_specs[ts][pg].pg_dim = pg_dims[ts][pg];
			pg_specs[ts][pg].pg_offset = pg_offsets[ts][pg];
			pg_specs[ts][pg].vardata = pg_datas[ts][pg];
		}
	}
}

void build_dataset_from_varblocks_by_var(
		const char *filename_prefix,
		const char *transform_name,
		const dataset_xml_spec_t *xml_spec,
		const dataset_global_spec_t *global_spec,
		const uint64_t pg_dims[global_spec->num_ts][global_spec->num_pgs_per_ts][xml_spec->ndim],
		const uint64_t pg_offsets[global_spec->num_ts][global_spec->num_pgs_per_ts][xml_spec->ndim],
		const void *varblocks_by_var[xml_spec->nvar])
{
	const uint64_t num_pgs = global_spec->num_ts * global_spec->num_pgs_per_ts;

	// Array for repackaging global per-variable data in per-PG data
	const void *varblocks_by_pg[global_spec->num_ts][global_spec->num_pgs_per_ts][xml_spec->nvar];
	// Array for repackaging pieces of per-PG inforamtion into an array of per-PG specification structs
	dataset_pg_spec_t pg_specs[global_spec->num_ts][global_spec->num_pgs_per_ts];

	collect_varblocks_by_pg(xml_spec, global_spec, pg_dims, varblocks_by_var, varblocks_by_pg);
	collect_pg_specs(global_spec->num_ts, global_spec->num_pgs_per_ts, xml_spec->ndim, xml_spec->nvar, pg_dims, pg_offsets, varblocks_by_pg, pg_specs);

	build_dataset_from_specs(filename_prefix, transform_name, xml_spec, global_spec, pg_specs);
}

void build_dataset_1(const char *filename_prefix, const char *transform_name) {
	// Basic dataset information
	// NOTE: we have to use an anonymous enum here to define these constants, since
	// C is picky and doesn't consider a static const int "const enough" to use
	// as an array length (e.g., if these were static const ints, it would not compile)
	enum {
		NUM_DIMS = 2,
		NUM_TS = 1,
		NUM_PGS_PER_TS = 1,
		NUM_VARS = 1,
		NUM_PGS = NUM_TS * NUM_PGS_PER_TS,
	};

	// Variable names/types
	static const char *VARNAMES[NUM_VARS]					= { "temp"     };
	static const enum ADIOS_DATATYPES VARTYPES[NUM_VARS]	= { adios_real };

	// Global and PG dimensions/offsets
	static const uint64_t GLOBAL_DIMS                         [NUM_DIMS] = { 4, 4 };
	static const uint64_t PG_DIMS	  [NUM_TS][NUM_PGS_PER_TS][NUM_DIMS] = { { { 4, 4 } } };
	static const uint64_t PG_OFFSETS  [NUM_TS][NUM_PGS_PER_TS][NUM_DIMS] = { { { 0, 0 } } };

	// Variable data (we can use [TS][PG][16] here because every PG is the same size, 16)
	static const float TEMP_DATA[NUM_TS][NUM_PGS_PER_TS][16] = {
		// Timestep 1
		// PG 0 in timestep 1
		//  1.00000000     1.00003052     2.00000000     2.00006104
			0x1.000000p+0, 0x1.000200p+0, 0x1.000000p+1, 0x1.000200p+1,
		//  2.00012207     30.00000000    30.00048828    30.00097656
			0x1.000400p+1, 0x1.e00000p+4, 0x1.e00200p+4, 0x1.e00400p+4,
		//  30.00146484    50.00000000    50.00097656    50.00195312
			0x1.e00600p+4, 0x1.900000p+5, 0x1.900200p+5, 0x1.900400p+5,
		//  50.00292969    50.00390625    50.00488281    50.00585938
			0x1.900600p+5, 0x1.900800p+5, 0x1.900a00p+5, 0x1.900c00p+5,
	};

	static const void *VARBLOCKS_BY_VAR[NUM_VARS] = {
		TEMP_DATA,
	};

	// Now, collect all this information into specification structs
	// File specification
	static const dataset_xml_spec_t XML_SPEC = {
		.group_name = "S3D",
		.buffer_size_mb = 128,
		.write_transport_method = "MPI",
		.ndim = NUM_DIMS,
		.nvar = NUM_VARS,
		.varnames = VARNAMES,
		.vartypes = VARTYPES,
	};

	// Global space specification
	static const dataset_global_spec_t GLOBAL_SPEC = {
		.num_ts = NUM_TS,
		.num_pgs_per_ts = NUM_PGS_PER_TS,
		.global_dims = GLOBAL_DIMS,
	};

	// Finally, invoke the dataset builder with this information
	build_dataset_from_varblocks_by_var(filename_prefix, transform_name, &XML_SPEC, &GLOBAL_SPEC, PG_DIMS, PG_OFFSETS, (const void **)VARBLOCKS_BY_VAR);
}

void build_dataset_2(const char *filename_prefix, const char *transform_name) {
	// Basic dataset information
	// NOTE: we have to use an anonymous enum here to define these constants, since
	// C is picky and doesn't consider a static const int "const enough" to use
	// as an array length (e.g., if these were static const ints, it would not compile)
	enum {
		NUM_DIMS = 3,
		NUM_TS = 2,
		NUM_PGS_PER_TS = 8,
		NUM_VARS = 1,
		NUM_PGS = NUM_TS * NUM_PGS_PER_TS,
	};

	// Variable names/types
	static const char *VARNAMES[NUM_VARS]					= { "temp"     };
	static const enum ADIOS_DATATYPES VARTYPES[NUM_VARS]	= { adios_real };

	// Global and PG dimensions/offsets
	static const uint64_t GLOBAL_DIMS                         [NUM_DIMS] = { 4, 4, 4 };
	static const uint64_t PG_DIMS	  [NUM_TS][NUM_PGS_PER_TS][NUM_DIMS] = {
		{ { 2, 2, 2 }, { 2, 2, 2 }, { 2, 2, 2 }, { 2, 2, 2 }, { 2, 2, 2 }, { 2, 2, 2 }, { 2, 2, 2 }, { 2, 2, 2 } }, // Timestep 1
		{ { 2, 2, 2 }, { 2, 2, 2 }, { 2, 2, 2 }, { 2, 2, 2 }, { 2, 2, 2 }, { 2, 2, 2 }, { 2, 2, 2 }, { 2, 2, 2 } }, // Timestep 2
	};
	static const uint64_t PG_OFFSETS  [NUM_TS][NUM_PGS_PER_TS][NUM_DIMS] = {
		{ { 0, 0, 0 }, { 0, 0, 2 }, { 0, 2, 0 }, { 0, 2, 2 }, { 2, 0, 0 }, { 2, 0, 2 }, { 2, 2, 0 }, { 2, 2, 2 } }, // Timestep 1
		{ { 0, 0, 0 }, { 0, 0, 2 }, { 0, 2, 0 }, { 0, 2, 2 }, { 2, 0, 0 }, { 2, 0, 2 }, { 2, 2, 0 }, { 2, 2, 2 } }, // Timestep 2
	};

	// Variable data (we can use [TS][PG][8] here because every PG is the same size, 8)
	static const float TEMP_DATA[NUM_TS][NUM_PGS_PER_TS][8] = {
		{ // Timestep 1
			// PG 0 in timestep 1
			{ 0x1.122f92p-5, 0x1.4e8798p-4, 0x1.d3cddep-1, 0x1.e16ae6p-1, 0x1.554fdep-1, 0x1.21c5a2p-1, 0x1.a8c114p-1, 0x1.8576bep-1, },
			// PG 1 in timestep 1
			{ 0x1.51e224p-2, 0x1.1a16eep-1, 0x1.ee4930p-1, 0x1.1d0e76p-1, 0x1.e3a098p-1, 0x1.fbf468p-1, 0x1.195c60p-2, 0x1.b53a84p-4, },
			// PG 2 in timestep 1
			{ 0x1.619b00p-1, 0x1.211a0ep-1, 0x1.022282p-2, 0x1.3ca72ap-2, 0x1.032826p-2, 0x1.8ea7eep-2, 0x1.2adcccp-1, 0x1.15f8c4p-3, },
			// PG 3 in timestep 1
			{ 0x1.b0a05ap-2, 0x1.7c8df4p-1, 0x1.eb4ba0p-4, 0x1.d8142ap-1, 0x1.818200p-1, 0x1.efe95cp-2, 0x1.979e1cp-1, 0x1.457d2ap-2, },
			// PG 4 in timestep 1
			{ 0x1.a66e56p-3, 0x1.f6a944p-1, 0x1.b9681ap-3, 0x1.8d1b8ap-1, 0x1.d6fc1cp-4, 0x1.d9f6a6p-3, 0x1.41fbc4p-1, 0x1.87e004p-1, },
			// PG 5 in timestep 1
			{ 0x1.0021acp-2, 0x1.c03452p-3, 0x1.c6fc40p-1, 0x1.86fa2ap-1, 0x1.f4e5e4p-1, 0x1.21594cp-2, 0x1.750404p-1, 0x1.75c0a6p-2, },
			// PG 6 in timestep 1
			{ 0x1.45eb00p-1, 0x1.d14d2ap-2, 0x1.f7962ap-1, 0x1.c3425ep-2, 0x1.5a94cap-3, 0x1.64ecfap-2, 0x1.1c473ep-1, 0x1.985a90p-6, },
			// PG 7 in timestep 1
			{ 0x1.ba2cb6p-1, 0x1.098d06p-1, 0x1.08cca2p-1, 0x1.65d43ep-2, 0x1.4682d0p-2, 0x1.d3bf16p-1, 0x1.25ce96p-1, 0x1.34e536p-2, },
		},
		{ // Timestep 2
			// PG 0 in timestep 1
			{ 0x1.34e536p-2, 0x1.25ce96p-1, 0x1.d3bf16p-1, 0x1.4682d0p-2, 0x1.65d43ep-2, 0x1.08cca2p-1, 0x1.098d06p-1, 0x1.ba2cb6p-1, },
			// PG 1 in timestep 1
			{ 0x1.985a90p-6, 0x1.1c473ep-1, 0x1.64ecfap-2, 0x1.5a94cap-3, 0x1.c3425ep-2, 0x1.f7962ap-1, 0x1.d14d2ap-2, 0x1.45eb00p-1, },
			// PG 2 in timestep 1
			{ 0x1.75c0a6p-2, 0x1.750404p-1, 0x1.21594cp-2, 0x1.f4e5e4p-1, 0x1.86fa2ap-1, 0x1.c6fc40p-1, 0x1.c03452p-3, 0x1.0021acp-2, },
			// PG 3 in timestep 1
			{ 0x1.87e004p-1, 0x1.41fbc4p-1, 0x1.d9f6a6p-3, 0x1.d6fc1cp-4, 0x1.8d1b8ap-1, 0x1.b9681ap-3, 0x1.f6a944p-1, 0x1.a66e56p-3, },
			// PG 4 in timestep 1
			{ 0x1.457d2ap-2, 0x1.979e1cp-1, 0x1.efe95cp-2, 0x1.818200p-1, 0x1.d8142ap-1, 0x1.eb4ba0p-4, 0x1.7c8df4p-1, 0x1.b0a05ap-2, },
			// PG 5 in timestep 1
			{ 0x1.15f8c4p-3, 0x1.2adcccp-1, 0x1.8ea7eep-2, 0x1.032826p-2, 0x1.3ca72ap-2, 0x1.022282p-2, 0x1.211a0ep-1, 0x1.619b00p-1, },
			// PG 6 in timestep 1
			{ 0x1.b53a84p-4, 0x1.195c60p-2, 0x1.fbf468p-1, 0x1.e3a098p-1, 0x1.1d0e76p-1, 0x1.ee4930p-1, 0x1.1a16eep-1, 0x1.51e224p-2, },
			// PG 7 in timestep 1
			{ 0x1.8576bep-1, 0x1.a8c114p-1, 0x1.21c5a2p-1, 0x1.554fdep-1, 0x1.e16ae6p-1, 0x1.d3cddep-1, 0x1.4e8798p-4, 0x1.122f92p-5, },
		},
	};

	static const void *VARBLOCKS_BY_VAR[NUM_VARS] = {
		TEMP_DATA,
	};

	// Now, collect all this information into specification structs
	// File specification
	static const dataset_xml_spec_t XML_SPEC = {
		.group_name = "S3D",
		.buffer_size_mb = 128,
		.write_transport_method = "MPI",
		.ndim = NUM_DIMS,
		.nvar = NUM_VARS,
		.varnames = VARNAMES,
		.vartypes = VARTYPES,
	};

	// Global space specification
	static const dataset_global_spec_t GLOBAL_SPEC = {
		.num_ts = NUM_TS,
		.num_pgs_per_ts = NUM_PGS_PER_TS,
		.global_dims = GLOBAL_DIMS,
	};

	// Finally, invoke the dataset builder with this information
	build_dataset_from_varblocks_by_var(filename_prefix, transform_name, &XML_SPEC, &GLOBAL_SPEC, PG_DIMS, PG_OFFSETS, (const void **)VARBLOCKS_BY_VAR);
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
	} else if (strcasecmp(dataset_id, "DS2") == 0) {
		dataset = DATASET_2;
	} else {
		fprintf(stderr, "Error: '%s' does not name a dataset packaged in this executable\n");
		usage_and_exit();
	}

	MPI_Init(&argc, &argv);

	switch (dataset) {
	case DATASET_1:
		build_dataset_1(path, transform_name);
		break;
	case DATASET_2:
		build_dataset_2(path, transform_name);
		break;
	}

	MPI_Finalize();
	return 0;
}
