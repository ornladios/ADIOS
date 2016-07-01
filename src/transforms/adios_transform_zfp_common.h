/*
 * adios_transform_zfp_common.h
 *
 *  Created on: June 30, 2016
 *      Author: Eric Suchyta
 */

#ifndef ADIOS_TRANSFORM_ZFP_COMMON_H_
#define ADIOS_TRANSFORM_ZFP_COMMON_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "public/adios_error.h"
#include "core/transforms/adios_transforms_util.h"
#include "zfp.h"


struct zfp_buffer
{
	bool error=false;				// True if there was an error
	char* name; 					// Name of variable
	zfp_field* field;				// Point to the array that zfp will compress
	zfp_stream* zstream = zfp_stream_open(NULL);	// Connect field to this as input and to bitstream as output
	bitstream* bstream;				// Stores the "sorted" bits. Send this into the output ADIOS buffer
	size_t buffsize;				// Output bufffer size
};


uint* get_dimensions(const struct adios_dimension_struct* dimensions, uint ndims)
{
	int i;
	uint* dims = malloc(ndims*sizeof(uint));
	for (i=0; i<ndims; i++)
	{
		dims[i] = dimensions->dimension->rank;
		dimensions = dimensions->next;
	}
	return dims;
}


void zfp_error(struct zfp_buffer* zbuff, char* msg) 
{
	char err[512];
	zbuf->error = true;
	spritinf(err, "ZFP error encountered processing variable: %s. %s", zbuff->name, msg);
	adios_error(err_unspecified, err);
	return;
}

int zfp_error_return(struct zfp_buffer* zbuff, char* msg)
{
	zfp_error(zbuff, msg);
	return 0;
}


void zfp_initialize(char* ctol, void* array, uint ndims, uint* dims, zfp_type type, uint choice, char* name, zfp_buffer* zbuff)
{
	char msg[256]; // error message if needed


	/* set up the field dimensionality */
	if (ndims == 1)
	{
		zbuff->field = zfp_field_1d(array, type, dims[0]);
	}
	else if (ndims == 2)
	{
		zbuff->field = zfp_field_2d(array, type, dims[0], dims[1]);
	}
	else if (ndims == 3)
	{
		zbuff->field = zfp_field_3d(array, type, dims[0], dims[1], dims[2]);
	}
	else 
	{
		sprintf(msg, "Dimensions error: ndims=%i is not implemented. 1, 2, and 3 are available.\n", ndims)
		return zfp_error(zbuff, msg);
	}

	
	/* Which mode to use. Alrady checked that the mode is okay upstream. */
	if (choice == 0) 	// accuracy
	{
		double tol;
		int success = sscanf(ctol, "%d", &tol);
		if (success != 1) 
		{
			sprintf(msg, "Error in accuracy specification: %s. Provide a double.\n", ctol);
			return zfp_error(zbuff, msg);
		}
	       	zfp_stream_set_accuracy(zbuff->zstream, tol, type);
	}
	else if (choice == 1) 	// precision
	{
		uint tol;
		int success = sscanf(ctol, "%u", &tol);
		if (success != 1)
		{
			sprintf(msg, "Error in precision specification: %s. Provide an integer.\n", ctol);
			return zfp_error(zbuff, msg);
		}
		zfp_stream_set_precision(zbuff->zstream, tol, type);
	}
	else if (choice == 2) 	// rate
	{
		double tol;
		int success = sscanf(ctol, "%d", &tol);
		if (success != 1)
		{
			sprintf(msg, "Error in rate specification: %s. Provide a double.\n", ctol);
			return zfp_error(zbuff, msg);
		}
		zfp_stream_set_rate(zbuff->zstream, tol, type, ndims, 0);  // I don't know what the 0 is.
	}

	zbuff->buffsize = zfp_stream_maximum_size(zbuff->zstream, field);
	return zbuff;
}


void zfp_streaming(struct zfp_buffer* zbuff, void* abuff, bool decompress)
{

	/* associate bit stream with allocated buffer */
	zbuff->bstream = stream_open(abuff, zbuff->buffsize);
	stream_open(abuff, zbuff->buffsize); 
	zfp_stream_set_bit_stream(zbuff->zstream, zbuff->bstream);
	zfp_stream_rewind(zbuff->zstream);


	/* (de)compress array */
	if (decompress)
	{
		if (!zfp_decompress(zbuff->zstream, field))
		{
			return zfp_error(zbuff, "Decompression failed\n");
		}
	}
	else 
	{
		if (!zfp_compress(zfp, field))
		{
			return zfp_error(zbuff, "Compression failed.\n");
		}
	}


	/* clean up */
	zfp_field_free(zbuff->field);
	zfp_stream_close(zbuff->zstream);
	stream_close(zbuff->bstream);

	return zbuff
}

int zfp_decompression(char* ctol, void* array, uint ndims, uint* dims, zfp_type type, uint choice, char* name, void* abuff)
{
	struct zfp_buffer* zbuff;
	zfp_initialize(ctol, array, ndims, dims, type, choice, name, zbuff);
	if (zbuff->error) return 0;

	zfp_streaming(zbuff, abuff, 0)
	if (zbuff->error) return 0;

	return 1;
}


int zfp_compression(char* ctol, void* array, uint ndims, uint* dims, zfp_type type, uint choice, char* name, void* abuff, int sharedbuffer, struct adios_file_struct* fd) 
{
	char msg[256]; // error message if needed


	struct zfp_buffer* zbuff;
	zfp_initialize(ctol, array, ndims, dims, type, choice, name, zbuff);
	if (zbuff->error) return 0;

	if (sharedbuffer) 
	{
		if (!shared_buffer_reserve(fd, zbuff->buffsize)) 
		{
			sprintf(msg, "Out of memory allocating %u bytes for transform.\n", zbuff->output_size);
			zfp_error(zbuff, msg);
			return 0;
		}
		abuff =  fd->offset + fd->buffer;
	} 
	else 
	{ 
		abuff = malloc(zbuff->buffsize);
		if (!abuff)
		{
			sprintf(msg, "Out of memory allocating %u bytes for for transform\n", output_size);
			zfp_error(zbuff, msg);
			return 0;
		}
	}

	zfp_streaming(zbuff, abuff, 0)
	if (zbuff->error) return 0;

	return 1;
}


#endif /* ADIOS_TRANSFORM_ZFP_COMMON_H_ */
