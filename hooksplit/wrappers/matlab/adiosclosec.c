/* 
 * ADIOS is freely available under the terms of the BSD license described
 * in the COPYING file in the top level directory of this source distribution.
 *
 * Copyright (c) 2008 - 2009.  UT-BATTELLE, LLC. All rights reserved.
 */

/*=================================================================
 * adiosclosec.c - Close an ADIOS file
 *
 * Input: fh, ghs,  Verbose
 *    fh:       int64 adios file handler
 *    ghs:      array of int64 group handlers
 *    Verbose:  numeric (double)
 *
 * $Revision: 1.0 $  $Date: 2009/08/05 12:53:41 $
 * Author: Norbert Podhorszki <pnorbert@ornl.gov>
 *=================================================================*/

#include "mex.h"
#include "adios_read.h"
#include "adios_types.h"

static int verbose=0;

mxClassID adiostypeToMatlabClass(int type, mxComplexity *complexity );
void errorCheck(int nlhs, int nrhs, const mxArray *prhs[]);
char* getString(const mxArray *mxstr);


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    int status;
    char msg[512];                               /* error messages from function calls */
    mwSize ngroups;
    int64_t *ghs, *int64p;               /* file and group handlers and temp pointers */
    int i;

    errorCheck(nlhs, nrhs, prhs);

    /*mexPrintf("nrhs=%d  nlhs=%d\n", nrhs, nlhs);*/
    
    /***********************/
    /* 0. get verbose parameter first */
    verbose = (int)mxGetScalar(prhs[2]);
    if (verbose) mexPrintf("Verbosity level: %d\n", verbose);

    /* 1. get file handler */
    int64p = (int64_t *) mxGetData(prhs[0]);
    if (verbose) mexPrintf("File handler: \"%lld\"\n", *int64p);

    /* 2. get group handlers */
    ghs = (int64_t *) mxGetData(prhs[1]);
    ngroups=mxGetNumberOfElements(prhs[1]);
    if (verbose) mexPrintf("Number of group handlers: \"%d\"\n", ngroups);

    for (i=0; i<ngroups; i++) {
        if (verbose) mexPrintf("Close group handler: \"%lld\"\n", ghs[i]);
        adios_gclose((ADIOS_GROUP *)ghs[i]);
    }
    if (verbose) mexPrintf("Close file handler: \"%lld\"\n", *int64p);
    adios_fclose((ADIOS_FILE *) *int64p);

    if (verbose) mexPrintf("return from adiosclosec\n");
}


void errorCheck(int nlhs, int nrhs, const mxArray *prhs[]){
    /* Assume that we are called from adiosread.m which checks the arguments already */
    /* Check for proper number of arguments. */
    
    if ( nrhs != 3 ) {
        mexErrMsgIdAndTxt("MATLAB:adiosclosec:rhs","This function needs exactly 3 arguments: FileHandler, GroupHandlers, Verbose");
    }
    
    if ( !mxIsInt64(prhs[0]) ) {
        mexErrMsgIdAndTxt("MATLAB:adiosclosec:rhs","First arg must be an int64 handler to an ADIOS file .");
    } 
    
    if ( !mxIsInt64(prhs[1]) ) {
        mexErrMsgIdAndTxt("MATLAB:adiosclosec:rhs","Second arg must be an int64 adios group handler.");
    } 

    if ( !mxIsNumeric(prhs[2]) ) {
        mexErrMsgIdAndTxt("MATLAB:adiosclosec:rhs","Third arg must be a number.");
    } 
    
    if ( nlhs > 0 ) {
        mexErrMsgIdAndTxt("MATLAB:adiosclosec:lhs","This function does not have output arguments.");
    }
    
#if !defined(MX_COMPAT_32)
    /* Make sure that it is safe to cast dim to mwSize when using largeArrayDims.*/
    if ( dim > MWSIZE_MAX ) {
        mexErrMsgIdAndTxt("MATLAB:adiosclosec:dimensionTooLarge",
                          "The input dimension, %.0f, is larger than the maximum value of mwSize, %u, when built with largeArrayDims.", dim, MWSIZE_MAX);
    }
#endif
 }


/** Make a C char* string from a Matlab string */
char* getString(const mxArray *mxstr) 
{
    mwSize buflen;
    char   *str;
    /* Allocate enough memory to hold the converted string. */
    buflen = mxGetNumberOfElements(mxstr) + 1;
    str = mxCalloc(buflen, sizeof(char));
    /* Copy the string data from string_array_ptr and place it into buf. */
    if (mxGetString(mxstr, str, buflen) != 0)
        mexErrMsgTxt("Could not convert string data from the file name.");
    return str;
}

