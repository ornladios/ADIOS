/*
 * adios_transform_read.h
 *
 *  Created on: Jul 23, 2012
 *      Author: David A. Boyuka II
 */

#ifndef ADIOS_READ_TRANSFORMED_H_
#define ADIOS_READ_TRANSFORMED_H_

#include "adios_transforms_common.h"
#include "public/adios_read.h"

// NCSU ALACRITY-ADIOS - struct for original metadata
typedef struct {
    enum ADIOS_TRANSFORM_TYPE transform_type;

    uint16_t transform_metadata_len;
    void *transform_metadata;
    int should_free_transform_metadata; // Used internally by read method and free

    enum ADIOS_DATATYPES orig_type;

    int orig_ndim;
    uint64_t *orig_dims;

    int orig_global;

    ADIOS_VARBLOCK *orig_blockinfo;
} ADIOS_TRANSINFO;

void adios_free_transinfo (ADIOS_TRANSINFO *ti);

#endif /* ADIOS_READ_TRANSFORMED_H_ */
