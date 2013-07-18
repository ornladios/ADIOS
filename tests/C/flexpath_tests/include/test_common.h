/**
 * @file: test_common.h
 * @author: Magda Slawinska, aka Magic Magg, magg dot gatech at gmail dot com
 * @date: Dec 19, 2012
 * Modified: Jan 03, 2013
 * The utility macros
 */

#ifndef TEST_COMMON_H_
#define TEST_COMMON_H_

#include <stdio.h>
#include <stdlib.h>


#define DBG_TEST_FAILED_STR		"TEST_FAILED"
#define DBG_TEST_PASSED_STR		"TEST_PASSED"

#define DEBUG

//! Debug printing verbosity
#define DBG_LEVEL   DBG_DEBUG

// New debug messaging state. There is no sense of a "level" for debugging. Each of these define the
// purpose of the messages and is enabled/disabled per file

//! system cannot continue, e.g. malloc
#define DBG_ERROR   0
//!
#define DBG_CRITICAL 1
//! some serious problems
#define DBG_WARNING 2
#define DBG_MESSAGE 3
//! messages about state or configuration; high-level flow
#define DBG_INFO    4
//!  func args, variable values, etc; full flow, may slow system down
#define DBG_DEBUG   5

#define DBG_ERROR_STR 		"ERROR\t"
#define DBG_CRITICAL_STR 	"CRITICAL\t"
#define DBG_WARNING_STR 	"WARNING\t"
#define DBG_MESSAGE_STR 	"MESSAGE\t"
#define DBG_INFO_STR		"INFO\t"
#define DBG_DEBUG_STR		"DEBUG\t"
#define DBG_TODO_STR		"TODO\t"

#define p_test_failed(fmt, args...)                             \
    do {                                                  \
         printf("%s %s:%s:%d: " fmt, DBG_TEST_FAILED_STR,  __FILE__, __FUNCTION__, __LINE__, ##args);  \
         fflush(stdout);  											\
    } while(0)

#define p_test_passed(fmt, args...)                             \
    do {                                                  \
         printf("%s %s:%s:%d: " fmt, DBG_TEST_PASSED_STR,  __FILE__, __FUNCTION__, __LINE__, ##args);  \
         fflush(stdout);  											\
    } while(0)



//! @todo do something like that but smarter without unnecessary copying
#define p_error(fmt, args...)                             				\
    do {                                                                \
        if((DBG_ERROR) <= DBG_LEVEL) {                                  \
            printf("%s(%d) %s:%s:%d: " fmt, DBG_ERROR_STR, (DBG_ERROR), __FILE__, __FUNCTION__, __LINE__, ##args);   \
            fflush(stdout);  											\
        }                                                               \
    } while(0)


#define p_warn(fmt, args...) \
	do {                                                                \
        if((DBG_WARNING) <= DBG_LEVEL) {                                      \
        	printf("%s(%d) %s:%s:%d: " fmt, DBG_ERROR_STR, (DBG_ERROR), __FILE__, __FUNCTION__, __LINE__, ##args);   \
            fflush(stdout);												\
        }                                                               \
    } while(0)

#define p_info(fmt, args...) \
	do {                                                                \
        if((DBG_INFO) <= DBG_LEVEL) {                                      \
        	printf("%s(%d) %s:%s:%d: " fmt, DBG_ERROR_STR, (DBG_ERROR), __FILE__, __FUNCTION__, __LINE__, ##args);   \
            fflush(stdout);												\
        }                                                               \
    } while(0)

#define p_debug(fmt, args...) \
	do {                                                                \
        if((DBG_DEBUG) <= DBG_LEVEL) {                                      \
        	printf("%s(%d) %s:%s:%d: " fmt, DBG_ERROR_STR, (DBG_ERROR), __FILE__, __FUNCTION__, __LINE__, ##args);   \
            fflush(stdout);												\
        }                                                               \
    } while(0)


// ------------------------------------
// other debugging util macro

/**
 * warns if the pointer is null
 * @param ptr the pointer to be checked
 * @param mesg The message to be displayed to the user
 */
#define ptr_null_warn(ptr, mesg)	\
	if( NULL == ptr )				\
		p_warn( "The pointer is NULL. %s\n", mesg);

// ADIOS UTILS
#define CLOSE_ADIOS \
	do { 														\
		adios_read_close(adios_handle);  						\
		adios_read_finalize_method(method);						\
		MPI_Finalize(); \
	} while (0)

#define JUST_CLEAN \
	do {								\
		adios_selection_delete(sel);	\
		sel = NULL;						\
		free(t);						\
		t = NULL;						\
	} while (0)

#define CLEAN_ON_ERROR_AND_CLOSE_ADIOS 	\
	do {								\
		JUST_CLEAN;						\
		CLOSE_ADIOS;					\
	} while (0)


#endif /* TEST_COMMON_H_ */
