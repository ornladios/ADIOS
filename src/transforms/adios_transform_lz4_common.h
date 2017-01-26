/*
 * adios_transform_lz4_common.h
 *
 * 	Author: Rene Widera
 * 	contact: r.widera@hzdr.de
 */
#ifndef ADIOS_TRANSFORM_LZ4_COMMON_H
#define ADIOS_TRANSFORM_LZ4_COMMON_H

#include "lz4.h"

/** LZ4 type to define sizes */
typedef int adiosLz4Size_t;

/** largest allowed input size (in byte) for the native LZ4 compression call */
#define ADIOS_LZ4_MAX_INPUT_SIZE LZ4_MAX_INPUT_SIZE

#endif /* ADIOS_TRANSFORM_LZ4_COMMON_H */
