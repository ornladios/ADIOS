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

#include <public/adios_error.h>
#include "zfp.h"


void get_dimensions(const struct adios_dimension_struct* dimensions, uint ndims)
{
	int i;
	int* dims = malloc(ndims*sizeof(int));
	for (i=0; i<ndims; i++)
	{
		dims[i] = dimensions->dimension->rank;
		dimensions = dimensions->next;
	}
	return 
}

void zfp_error(int &status, char* msg) 
{
	adios_error(err_unspecified, msg);
	status = 1;
}


int _zfp_comp(char* ctol, void* array, uint ndims, uint* dims, zfp_type type, uint choice, uint decompress, void* abuffer) 
{

	int status = 0;					// return value: 0 = success, 1 = error
	zfp_field* field;				// Roughly speaking this "stores" the array that zfp will compress
	zfp_stream* zfp = zfp_stream_open(NULL);	// Connect field to this as input, and bitstream as output
	bitstream* stream;				// Roughly speaking this stores the sorted bits. Send this into the output ADIOS abuffer


	/* set up the field dimensionality */
	if      (ndims==1) field = zfp_field_1d(array, type, dims[0]);
	else if (ndims==2) field = zfp_field_2d(array, type, dims[0], dims[1]);
	else if (ndims==3) field = zfp_field_3d(array, type, dims[0], dims[1], dims[2]); 	
	else zfp_error( sprintf("ZFP dimensions error: ndims=%i is not implemented. 1, 2, and 3 are available\n", ndims) );

	
	/* Which mode to use */
	if (choice==0) 		// accuracy
	{
		double tol;
		int success = sscanf(ctol, "%d", &tol);
		if (success != 1) zfp_error("ZFP error in accuracy specification: %s. Provide a double.", ctol);
	       	zfp_stream_set_accuracy(zfp, tol, type);
	}
	else if (choice==1) 	// precision
	{
		uint tol;
		int success = sscanf(ctol, "%u", &tol);
		if (success != 1) zfp_error("ZFP error in precision specification: %s. Provide an integer.", ctol);
		zfp_stream_set_precision(zfp, tol, type);
	}
	else if (choice==2) 	// rate
	{
		double tol;
		int success = sscanf(ctol, "%d", &tol);
		if (success != 1) zfp_error("ZFP error in rate specification: %s. Provide a double.", ctol);
		zfp_stream_set_rate(zfp, tol, type, ndims, 0);  // The fourth arg is number of dimensions. I don't know what the 0 is.
	}


	/* allocate buffer for compressed data */
	size_t bufsize = zfp_stream_maximum_size(zfp, field);
	abuffer = malloc(bufsize);


	/* associate bit stream with allocated buffer */
	stream = stream_open(abuffer, bufsize);
	stream_open(abuffer, bufsize); 
	zfp_stream_set_bit_stream(zfp, stream);
	zfp_stream_rewind(zfp);


	/* (de)compress array */
	if (decompress) if (!zfp_decompress(zfp, field)) zfp_error("ZFP decompression failed\n");
	else if (!zfp_compress(zfp, field)) zfp_error("ZFP compression failed.\n");


	/* clean up */
	zfp_field_free(field);
	zfp_stream_close(zfp);
	stream_close(stream);

	return status;
}


#endif /* ADIOS_TRANSFORM_ZFP_COMMON_H_ */
