/*
 * adios_transform_zfp_common.h
 *
 * 	Author: Eric Suchyta
 * 	contact: eric.d.suchyta@gmail.com
 *
 */

#ifndef ADIOS_TRANSFORM_ZFP_COMMON_H_
#define ADIOS_TRANSFORM_ZFP_COMMON_H_


/* general C stuff */
#include <stdint.h> 	// uint64_t
#include <stdio.h>	// NULL, sprintf, sscanf
#include <stdlib.h>	// NULL, malloc, free
#include <string.h>	// memcpy, strcmp, strlen


/* Were in the template included from ADIOS. Not necessarily sure if they're all strictly needed. */
#include "core/transforms/adios_transforms_common.h"
#include "core/transforms/adios_transforms_write.h"
#include "core/transforms/adios_transforms_hooks_write.h"
#include "core/transforms/adios_transforms_util.h"
#include "core/transforms/adios_transforms_hooks_read.h"
#include "core/transforms/adios_transforms_reqgroup.h"


/* Extra ADIOS headers that weren't added in the template */
#include "adios_logger.h" 		// log_warn()
#include "public/adios_error.h"		// adios_error
#include "public/adios_types.h"		// adios datatypes


/* Sort of a hack to deal with ADIOS and ZFP both using TYPES_H */
#ifdef TYPES_H
#undef TYPES_H
#endif 


/* ZFP specific */
#include "zfp.h"
#define ZFP_STRSIZE 256		// size for string variables


/* Transform metadata */
struct zfp_metadata 
{
	uint64_t usize;			// uncompressed size;
	uint64_t csize;			// compressed size	
	uint cmode;			// compression mode
	char ctol[ZFP_STRSIZE];		// string of "tolerance"
	char name[ZFP_STRSIZE];		// variable name
};


/* Attributes related to ZFP's operation */
struct zfp_buffer
{
	bool error;					// true if error
	
	char name[ZFP_STRSIZE]; 			// Name of variable
	zfp_type type;					// data type
	uint mode;					// 0 = accuracy, 1 = precsion, 2 = rate
	char ctol[ZFP_STRSIZE];				// string for "tolerance"
	uint ndims; 					// number of dimensions
	uint* dims;					// array of dimension sizes

	zfp_field* field;				// Point to the array that zfp will compress
	zfp_stream* zstream;				// Connect field to this as input and to bitstream as output
	bitstream* bstream;				// Stores the "sorted" bits. Send this into the output ADIOS buffer
	size_t buffsize;				// Output bufffer size
};


/* Set a few defaults */
static void init_zfp_buffer(struct zfp_buffer* zbuff, char* name)
{
	strcpy(zbuff->name, name);
	zbuff->error = false;
	return;
}


/* adios to zfp datatype */
static int zfp_get_datatype(struct zfp_buffer* zbuff, enum ADIOS_DATATYPES type)
{
	if (type == adios_double) 
	{
		zbuff->type = zfp_type_double;
	}
	else if (type == adios_real) 
	{
		zbuff->type = zfp_type_float;
	}
	else 
	{
		adios_error(err_unspecified, "ZFP does not handle the type of variable %s. "
				"Supported types are adios_double, adios_real.\n",
				zbuff->name);
		zbuff->error = true;
		return 0;
	}
	return 1;
}


/* Configure ZFP according to the user's specified mode and tolerance */
static void zfp_initialize(void* array, struct zfp_buffer* zbuff)
{
	zbuff->zstream = zfp_stream_open(NULL);

	/* set up the field dimensionality */
	if (zbuff->ndims == 1)
	{
		zbuff->field = zfp_field_1d(array, zbuff->type, zbuff->dims[0]);
	}
	else if (zbuff->ndims == 2)
	{
		zbuff->field = zfp_field_2d(array, zbuff->type, zbuff->dims[0], zbuff->dims[1]);
	}
	else if (zbuff->ndims == 3)
	{
		zbuff->field = zfp_field_3d(array, zbuff->type, zbuff->dims[0], zbuff->dims[1], zbuff->dims[2]);
	}
	else 
	{
		adios_error(err_invalid_dimension, "ZFP does not handle the %u dimensional variable %s. "
			"Only 1, 2, and 3 dimensions are handled.\n",
			zbuff->ndims, zbuff->name);
		zbuff->error = true;
		return;
	}

	
	/* Which mode to use. Alrady checked that the mode is okay upstream. */
	if (zbuff->mode == 0) 	// accuracy
	{
		double tol;
		int success = sscanf(zbuff->ctol, "%lf", &tol);
		if (success != 1) 
		{
			adios_error(err_invalid_argument, "Error in accuracy specification for variable %s: %s. "
					"Provide a double value.\n",
					zbuff->name, zbuff->ctol);
			zbuff->error = true;
			return;
		}
	       	zfp_stream_set_accuracy(zbuff->zstream, tol, zbuff->type);
	}
	else if (zbuff->mode == 1) 	// precision
	{
		uint tol;
		long int ct;
		char* end;
		ct = strtol(zbuff->ctol, &end, 10);
		if (ct == 0)
		{
			adios_error(err_invalid_argument, "Error in precision specification for variable %s: %s. "
					"Provide an integer value.\n",
					zbuff->name, zbuff->ctol);
			zbuff->error = true;
			return;
		}

		if (*end != '\0')
		{
			log_warn("A float was given for ZFP precision for variable %s: %s "
					"-- the value was cast to an integer. ZFP accepts integer precisions.",
					zbuff->name, zbuff->ctol);
		}
		tol = (uint) ct;


		zfp_stream_set_precision(zbuff->zstream, tol, zbuff->type);
	}
	else if (zbuff->mode == 2) 	// rate
	{
		double tol;
		int success = sscanf(zbuff->ctol, "%lf", &tol);
		if (success != 1)
		{
		       	adios_error(err_invalid_argument, "Error in rate specification for variable %s: %s. "
					"Provide a double value.\n",
					zbuff->name, zbuff->ctol);
			zbuff->error = true;
			return;
		}
		zfp_stream_set_rate(zbuff->zstream, tol, zbuff->type, zbuff->ndims, 0);  // I don't know what the 0 is.
	}
	
	zbuff->buffsize = zfp_stream_maximum_size(zbuff->zstream, zbuff->field);
}


/* Do the bit streaming */
static void zfp_streaming(struct zfp_buffer* zbuff, void* abuff, bool decompress, uint64_t* finalsize)
{

	/* associate bit stream with allocated buffer */
	zbuff->bstream = stream_open(abuff, zbuff->buffsize);
	zfp_stream_set_bit_stream(zbuff->zstream, zbuff->bstream);
	zfp_stream_rewind(zbuff->zstream);

	/* (de)compress array */
	if (decompress)
	{
		int success = zfp_decompress(zbuff->zstream, zbuff->field);
		if (!success)
		{
			adios_error(err_transform_failure, "ZFP decompression failed for variable %s\n", zbuff->name);
			zbuff->error = true;
			return;
		}
	}
	else 
	{
		*finalsize = (uint64_t) zfp_compress(zbuff->zstream, zbuff->field);
		if (! *finalsize)
		{
			adios_error(err_transform_failure, "ZFP compression failed for variable %s\n", zbuff->name);
			zbuff->error = true;
			return;
		}
	}

	/* clean up */
	zfp_field_free(zbuff->field);
	zfp_stream_close(zbuff->zstream);
	stream_close(zbuff->bstream);
	free(zbuff->dims);
}


#endif /* ADIOS_TRANSFORM_ZFP_COMMON_H_ */
