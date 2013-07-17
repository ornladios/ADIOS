/**
 * misc.h
 *
 *  Created on: Jul 5, 2013
 *  Author: Magda Slawinska aka Magic Magg magg dot gatech at gmail.com
 */

#ifndef MISC_H_
#define MISC_H_

//! if defined (not commented) the flexpath method will
//! be used; otherwise the  ADIOS_READ_METHOD_BP and MPI
//! if you are switching methods, make sure that arrays.xml
//! contains the correct method
//! you can use the -DFLEXPATH_METHOD in make CFLAGS="-DFLEXPATH_METHOD"
//! to turn this on as well
//#define FLEXPATH_METHOD 1

//! the name of the file to be written test values
#define FILE_NAME "test.bp"

//! size of the X dimension
#define NX_DIM 10

//! the xml containing configuration of ADIOS
#ifdef FLEXPATH_METHOD
#define XML_ADIOS_INIT_FILENAME "test_config_flex.xml"
#define METHOD ADIOS_READ_METHOD_FLEXPATH
#else
#define XML_ADIOS_INIT_FILENAME "test_config_mpi.xml"
#define METHOD ADIOS_READ_METHOD_BP
#endif

//! options how verbose is ADIOS (see  adios_read_init_method)
//! 0=quiet, ..., 4=debug
#define ADIOS_OPTIONS "verbose=4; show hidden_attrs"

//! defines if the test passed
#define TEST_PASSED 0
//! defines if the test failed
#define TEST_FAILED -1

//! indicates that the program
#define PROGRAM_ERROR -2

#endif /* MISC_H_ */
