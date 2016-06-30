/*
 * adios_transform_zfp_common.h
 *
 *  Created on: June 30, 2016
 *      Author: Eric Suchyta
 */

#ifndef ADIOS_TRANSFORM_ZFP_COMMON_H_
#define ADIOS_TRANSFORM_ZFP_COMMON_H_

#include <stdint.h>
#include <assert.h>
#include <public/adios_error.h>
//#include <assert.h>

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "zfp.h"


void zfp_error(int &status, char* msg) {
	adios_error(err_unspecified, msg);
	status = 1;
}


int _zfp_allocate(void* tol, void* array, int ndims, int* dims, zfp_type type, uint choice, uint decompress, void* outbuffer) {
	int status = 0;		// return value: 0 = success, 1 = error
	zfp_field* field;	// Roughly speaking this "stores" the array that zfp will compress
	zfp_stream* zfp;	// Connect field to this as input, and bitstream as output
	bitstream* stream;	// Roughly speaking this stores the sorted bits. Send this into the output ADIOS outbuffer

	/* set up the field you're compressing */
	if      (ndims==1) field = zfp_field_1d(array, type, dims[0]);
	else if (ndims==2) field = zfp_field_2d(array, type, dims[0], dims[1]);
	else if (ndims==3) field = zfp_field_3d(array, type, dims[0], dims[1], dims[2]); 	
	else zfp_error( sprintf("ZFP dimensions error: ndims=%i is not implemented. 1, 2, and 3 are available\n", ndims) );

	/* set compression mode and parameters via one of three functions */
	zfp = zfp_stream_open(NULL);
	if      (choice==0) zfp_stream_set_accuracy(zfp, *((double*)tol), type);
	else if (choice==1) zfp_stream_set_precision(zfp, *((uint*)tol), type);
	else if (choice==2) zfp_stream_set_rate(zfp, *((double*)tol), type, ndims, 0);  // The fourth arg is number of dimensions. I don't know what the 0 is.

	/* allocate buffer for compressed data */
	size_t bufsize = zfp_stream_maximum_size(zfp, field);	// byte size of compressed buffer
	outbuffer = malloc(bufsize);

	/* associate bit stream with allocated buffer */
	stream = stream_open(outbuffer, bufsize);
	stream_open(outbuffer, bufsize); zfp_stream_set_bit_stream(zfp, stream);
	zfp_stream_rewind(zfp);


	if (decompress) {
		// read compressed stream and decompress array
		FILE* file = fopen("test_zfp.out", "r");
		int zfpsize = fread(buffer, 1, bufsize, file);
		if (!zfp_decompress(zfp, field)) {
			fprintf(stderr, "decompression failed\n");
			status = 1;
		}
	}

	else {
		/* compress array into stream */
		int zfpsize = zfp_compress(zfp, field);
		if (!zfpsize) zfp_error("ZFP compression failed at the zfp_compress step.\n");
	}


	/* clean up */
	zfp_field_free(field);
	zfp_stream_close(zfp);
	stream_close(stream);
	free(buffer);
	if (!decompress) free(array);

	return status;
}

int zfp_compress_accuracy(double tolerance, void* array, uint ndims, uint* dims, zfp_type type) { 
	return _zfp_compress_choose(&tolerance, array, ndims, dims, type, 0, 0);
}

int zfp_decompress_accuracy(double tolerance, void* array, uint ndims, uint* dims, zfp_type type) { 
	return _zfp_compress_choose(&tolerance, array, ndims, dims, type, 0, 1);
}


int zfp_compress_precision(uint precision, void* array, uint ndims, uint* dims, zfp_type type) { 
	return _zfp_compress_choose(&precision, array, ndims, dims, type, 1, 0);
}

int zfp_decompress_precision(uint precision, void* array, uint ndims, uint* dims, zfp_type type) { 
	return _zfp_compress_choose(&precision, array, ndims, dims, type, 1, 1);
}


int zfp_compress_rate(double rate, void* array, uint ndims, uint* dims, zfp_type type) {
	return _zfp_compress_choose(&rate, array, ndims, dims, type, 2, 0);
}

int zfp_decompress_rate(double rate, void* array, uint ndims, uint* dims, zfp_type type) {
	return _zfp_compress_choose(&rate, array, ndims, dims, type, 2, 1);
}


#endif /* ADIOS_TRANSFORM_ZFP_COMMON_H_ */
