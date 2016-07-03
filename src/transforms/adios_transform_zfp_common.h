/*
 * adios_transform_zfp_common.h
 *
 * 	Author: Eric Suchyta
 * 	contact: eric.d.suchyta@gmail.com
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


/* Extra ADIOS headers that weren't added in the template */
#include "adios_logger.h" 		// log_warn()
#include "public/adios_error.h"		// adios_error
#include "public/adios_types.h"		// adios datatypes


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


/* Basically giving an address and a size */
void* zfp_read_metadata_var(void* pos, size_t size, size_t &offset)
{
	offset += size;
	return pos + offset - size;
}


/* Basically giving an address and a size */
void zfp_write_metadata_var(char* pos, void* towrite, size_t size, size_t &offset)
{
	memcpy(pos + offset, towrite, size);
	offset += size;
	return;
}


/* Read each memory location and cast to the correct type */
struct zfp_metadata* zfp_read_metadata(adios_transform_pg_read_request *completed_pg_reqgroup)
{
	struct zfp_metadata* metadata;
	void* pos = completed_pg_reqgroup->transform_metadata;
	
	metadata->usize = *((uint64_t*)zfp_read_metadata_var(pos, sizeof(uint64_t), offset));
	metadata->csize = *((uint64_t*)zfp_read_metadata_var(pos, sizeof(uint64_t), offset));
	metadata->cmode = *((uint*)zfp_read_metadata_var(pos, sizeof(uint), offset));
	metadata->ctol = *((char*)zfp_read_metadata_var(pos, ZFP_STRSIZE, offset));
	metadata->name = *((char*)zfp_read_metadata_var(pos, ZFP_STRSIZE, offset));
	return metadata;
}


/* Attributes related to ZFP's operation */
struct zfp_buffer
{
	bool error;					// true if error
	char msg[ZFP_STRSIZE];				// error/warning message, if needed
	
	char name[ZFP_STRSIZE]; 			// Name of variable
	zfp_type type;					// data type
	uint mode;					// 0 = accuracy, 1 = precsion, 2 = rate
	char ctol[ZFP_STRSIZE];				// string for "tolerance"
	uint ndims; 					// number of dimensions
	uint* dims;					// array of dimension sizes

	zfp_field* field;				// Point to the array that zfp will compress
	zfp_stream* zstream = zfp_stream_open(NULL);	// Connect field to this as input and to bitstream as output
	bitstream* bstream;				// Stores the "sorted" bits. Send this into the output ADIOS buffer
	size_t buffsize;				// Output bufffer size
};


/* Get the dimensionality of the input data */
void get_dimensions(const struct adios_dimension_struct* dimensions, struct zfp_buffer* zbuff)
{
	int i;
	zbuff->dims = malloc(zbuff->ndims*sizeof(uint));
	for (i=0; i<zbuff->ndims; i++)
	{
		zbuff->dims[i] = dimensions->dimension->rank;
		dimensions = dimensions->next;
	}
}


/* Function for common way to log errors */
void zfp_error(struct zfp_buffer* zbuff) 
{
	zbuff->error = True;
	sprintf(zbuff->msg + strlen(zbuff->msg), "ERROR in zfp processing variable: %s. %s", zbuff->name, zbuff->msg);
	adios_error(err_unspecified, zbuff->msg);
	return;
}


/* Function for common way to log warning */
void zfp_warn(struct zfp_buffer* zbuff) 
{
	sprintf(zbuff->msg + strlen(zbuff->msg), "WARNING in zfp processing variable: %s. %s", zbuff->name, zbuff->msg);
	log_warn(err_unspecified, zbuff->msg);
	return;
}


/* adios to zfp datatype */
void zfp_get_datatype(struct zfp_buffer* zbuff, enum ADIOS_DATATYPES type)
{
	if (type == adios_double) 
	{
		zbuff->type = zfp_type_double;
	}
	else if (type == adios_float) 
	{
		zbuff->type = zfp_type_float;
	}
	else 
	{
		sprintf(zbuff->msg, "A datatype ZFP does not understand was given for compression. Understood types are adios_double, adios_float.");
		zfp_error(zbuff);
		return 0;
	}
	return 1;
}


/* Configure ZFP according to the user's specified mode and tolerance */
void zfp_initialize(void* array, zfp_buffer* zbuff)
{

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
		sprintf(zbuff->msg, "Dimensions error: ndims=%i is not implemented. 1, 2, and 3 are available.\n", zbuff->ndims);
		return zfp_error(zbuff);
	}

	
	/* Which mode to use. Alrady checked that the mode is okay upstream. */
	if (zbuff->mode == 0) 	// accuracy
	{
		double tol;
		int success = sscanf(zbuff->ctol, "%d", &tol);
		if (success != 1) 
		{
			sprintf(zbuff->msg, "Error in accuracy specification: %s. Provide a double.\n", zbuff->ctol);
			return zfp_error(zbuff);
		}
	       	zfp_stream_set_accuracy(zbuff->zstream, tol, zbuff->type);
	}
	else if (zbuff->mode == 1) 	// precision
	{
		uint tol;
		int success = sscanf(zbuff->ctol, "%u", &tol);
		if (success != 1)
		{
			sprintf(zbuff->msg, "Error in precision specification: %s. Provide an integer.\n", zbuff->ctol);
			return zfp_error(zbuff);
		}
		zfp_stream_set_precision(zbuff->zstream, tol, zbuff->type);
	}
	else if (zbuf->mode == 2) 	// rate
	{
		double tol;
		int success = sscanf(zbuff->ctol, "%d", &tol);
		if (success != 1)
		{
			sprintf(zbuff->msg, "Error in rate specification: %s. Provide a double.\n", zbuff->ctol);
			return zfp_error(zbuff);
		}
		zfp_stream_set_rate(zbuff->zstream, tol, zbuff->type, zbuff->ndims, 0);  // I don't know what the 0 is.
	}

	zbuff->buffsize = zfp_stream_maximum_size(zbuff->zstream, field);
	return zbuff;
}


/* Do the bit streaming */
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
		/*
		if (!zfp_decompress(zbuff->zstream, zbuff->field))
		{
			sprintf(zbuff->msg, "Decompression failed\n");
			return zfp_error(zbuff);
		}
		*/
	}
	else 
	{
		if (!zfp_compress(zbuff->zstream, zbuff->field))
		{
			sprintf(zbuff->msg, "Compression failed.\n")
			return zfp_error(zbuff);
		}
	}


	/* clean up */
	zfp_field_free(zbuff->field);
	zfp_stream_close(zbuff->zstream);
	stream_close(zbuff->bstream);

	return zbuff
}



/* This is called in the main transform-level function. 
 * In a nutshell: compress array, using the settings in the other args to configure the compression. Connect to the ADIOS output buffer.
 */
int zfp_compression(zfp_buffer* zbuff, void* array, void* abuff, int &asize, int sharedbuffer, struct adios_file_struct* fd) 
{
	zfp_initialize(array, zbuff);
	if (zbuff->error) 
	{
		return 0;
	}


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


	zfp_streaming(zbuff, abuff, 0);
	if (zbuff->error)
	{
		return 0;
	}

	
	asize = (uint64_t) zbuff->buffsize;
	return 1;
}


/* This is called in the main transform-level function. 
 * In a nutshell: decompress array, using (undoing) the settings in the other args. Connect to the ADIOS buffer.
 */
int zfp_decompression(zfp_buffer* zbuff, void* uarray, void* carray, uint64_t buffsize)
{
	zfp_initialize(uarray, zbuff);
	if (zbuff->error)
	{
		return 0;
	}
	if (zbuff->buffsize != buffsize)
	{
		sprintf(zbuff->msg, "ZFP thinks compressed size should be %u" \
				"bytes. ADIOS thinks compressed size should be %" PRIu64 \
				"bytes. Likely corruption.\n", zbuff->buffsize, buffsize);
		zfp_warn(zbuff);
	}


	zfp_streaming(zbuff, carray, 1)
	if (zbuff->error)
	{
		return 0;
	}

	return 1;
}


#endif /* ADIOS_TRANSFORM_ZFP_COMMON_H_ */
