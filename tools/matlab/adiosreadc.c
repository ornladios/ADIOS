/*=================================================================
 * adiosreadc.c - read a variable/attribute from an ADIOS file
 *
 * Input: File, Group, Path, Time, Offset, Count, Verbose
 *    File:     string (name) or int64 (handler)
 *    Group:    string (name) or int64 (handler) or int32 (index)
 *    Path:     string (variable/attribute name with full path)
 *    Time:     int32  (timestep, 1..)
 *    Offset:   int32 array 
 *    Count:    int32 array
 *    Verbose:  numeric (double)
 * Output: The variable/attribute as mxArray
 *
 *
 * Copyright 2009 Oak Ridge National Laboratory
 * $Revision: 1.0 $  $Date: 2009/08/05 12:53:41 $
 * Author: Norbert Podhorszki <pnorbert@ornl.gov>
 *=================================================================*/

#include "mex.h"
#include "adios_readutil.h"
#include "adios_types.h"

static int verbose=0;

mxArray* readdata( int64_t fh, int64_t gh, const char *path, int in_noffsets, 
                   const int32_t *in_offsets, const int32_t *in_counts,
                   int32_t in_timestep, int *isvar);
void errorCheck(int nlhs, int nrhs, const mxArray *prhs[]);
char* getString(const mxArray *mxstr);
mxArray* createMatlabArray( int adiostype, int ndims, int *dims); 
void recalc_offsets( int ndims, int *dims, mwSize in_noffsets, 
                     const int32_t *in_offsets, const int32_t *in_counts,
                     int *offsets, int *counts);
void swap_order(int n, int *array);



void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    char *fname;        /* file name */
    char *gname;        /* group name */
    char *path;         /* name of variable or attribute */
    mwSize in_noffsets;    /* number of offsets provided, = number of counts */
    int32_t *in_offsets;    /* offset array provided */
    int32_t *in_counts;     /* count array provided */
    int status;
    char msg[512];        /* error messages from function calls */

    int64_t fh, gh, *int64p; /* file and group handlers and temp pointers */
    int32_t timestep, *int32p;
    int  haveToCloseFile, haveToCloseGroup, readAttributesToo;

    int  isvar;             /* 0: attribute, 1: variable */
    adios_readutil_namelist *attrlist;  /* linked list of attrs directly under a path */
    adios_readutil_namelist *attr;      /* one attr from the list */
    int  nattrs;            /* number of attrs directly under a path */

    int i;
    mxArray *out;        /* output data array */
    const char *field_names[] = {"name", "value"}; /* 2nd output arg's struct fields */
    mwSize attr_struct_dims[2];   /* dimensions for 2nd arg struct array: 1-by-sth */

    errorCheck(nlhs, nrhs, prhs);

    /*mexPrintf("nrhs=%d  nlhs=%d\n", nrhs, nlhs);*/
    
    /***********************/
    /* 1. get verbose parameter first */
    verbose = (int)mxGetScalar(prhs[6]);
    if (verbose) mexPrintf("Verbosity level: %d\n", verbose);
    adios_readutil_setverbosity(verbose);

    /***********************/
    /* 1. get file handler */
    if (mxIsChar(prhs[0])) {
        fname = getString( (mxArray *)prhs[0]);
        if (verbose) mexPrintf("File name: \"%s\"\n", fname);

        /* Open ADIOS file now */
        status = adios_readutil_fopen (fname, &fh, msg); 
        if (status != 0) {
            mexErrMsgIdAndTxt("MATLAB:adiosreadc:open",msg);
        }
        haveToCloseFile = 1;

    } else  { /* int64 handler provided */
        int64p = (int64_t *) mxGetData(prhs[0]); 
        fh = *int64p;
        haveToCloseFile = 0;
    } 


    /************************/
    /* 2. get group handler */
    if (mxIsChar(prhs[1])) {
        gname = getString( (mxArray *)prhs[1]);
        if (verbose) mexPrintf("Group name: \"%s\"\n", gname );

        /* Open group by name */
        status = adios_readutil_gopen_byname (gname, fh, &gh, msg); 
        if (status != 0) {
            mexErrMsgIdAndTxt("MATLAB:adiosreadc:open",msg);
        }
        haveToCloseGroup = 1;

    } else if (mxIsInt32(prhs[1])) { /* group index provided */
        int32p  = (int32_t *) mxGetData(prhs[1]); 
        if (verbose) mexPrintf("Group index: %d\n", *int32p);

        /* Open group by index */
        status = adios_readutil_gopen_byindex (*int32p, fh, &gh, msg); 
        if (status != 0) {
            mexErrMsgIdAndTxt("MATLAB:adiosreadc:open",msg);
        }
        haveToCloseGroup = 1;

    } else  { /* int64 handler provided */
        int64p = (int64_t *) mxGetData(prhs[1]); 
        gh = *int64p;
        haveToCloseGroup = 0;
    } 


    /*****************************************************************************************/
    /* 3. get other arguments: char variable name, int32 time, int32 in_offsets[], int32 in_counts[] */
    path     = getString(prhs[2]);
    if (verbose) mexPrintf("Variable name: \"%s\"\n", path );
    int32p   = (int32_t *) mxGetData(prhs[3]); 
    timestep = *int32p;
    if (verbose) mexPrintf("Timestep: %d\n", timestep );
    in_noffsets = mxGetNumberOfElements(prhs[4]);
    in_offsets  = (int32_t *) mxGetData(prhs[4]); 
    in_counts   = (int32_t *) mxGetData(prhs[5]); 


    /*********************************************************************/
    /* 4. read in variable/attribute or variable with all its attributes */
    out = readdata( fh, gh, path, in_noffsets, in_offsets, in_counts, timestep, &isvar);
    if ( nlhs >= 1 ) {
        plhs[0] = out;
    }

    /********************************************************/
    /* 5. read in all attributes of a variable if requested */
    if ( nlhs == 2 ) {
        if (!isvar) {
            mexErrMsgIdAndTxt("MATLAB:adiosreadc:read",
               "Second output argument can be provided only for variables. Path %s refers to an attribute.",
               path);
        } else {
            /* Read all attributes directly under the given path (= path/name) */
            if (verbose) mexPrintf("Read in attributes under the path %s\n", path);
            nattrs = adios_readutil_getdirectattributes(gh, path, &attrlist); 
            if (verbose) mexPrintf("Found %d matching attributes\n", nattrs);
            attr_struct_dims[0] = 1;
            attr_struct_dims[1] = nattrs;
            /* Create a 1-by-n array of structs with fields {name,value}. */ 
            plhs[1] = mxCreateStructArray(2, attr_struct_dims, 2, field_names);
            if (verbose) mexPrintf("Created struct array, 1-by-%d\n", attr_struct_dims[1]);
            i=0;
            attr=attrlist;
            while (attr != NULL) {
                mxSetFieldByNumber(plhs[1],i,0,mxCreateString(attr->name));
                out = readdata( fh, gh, attr->name, 0, NULL, NULL, timestep, &isvar);
                mxSetFieldByNumber(plhs[1],i,1,out);
                if (verbose) mexPrintf("Added attr: %s\n", attr->name);
                i++;
                attr=attr->next;
            }
            adios_readutil_freenamelist(attrlist);
        }
    }


    /**************************************************/
    /* 6. close group and file if opened in this call */
    if (haveToCloseGroup) {
        if (verbose) mexPrintf("Close group\n");
        adios_readutil_gclose(gh);
        mxFree(gname);
    }

    if (haveToCloseFile) {
        if (verbose) mexPrintf("Close file\n");
        adios_readutil_fclose(fh);
        mxFree(fname);
    }


    mxFree(path);

}

mxArray* readdata( int64_t fh, int64_t gh, const char *path, mwSize in_noffsets,
                   const int32_t *in_offsets, const int32_t *in_counts,
                   int32_t in_timestep, int *isvar) 
{
    void *data;             /* content of variable/attr read in */
    int  ndims, dims[16];   /* size info for variable/attr read in */
    int  type;              /* type of variable/attr read in */
    int32_t timestep;
    int offsets[16], counts[16]; /* extended offsets/counts */
    int  status, i;
    mxArray * out;
    char msg[512];

    /* read in variable/attribute or variable with all its attributes */
    if (verbose) mexPrintf("Get info on var/attr: %s\n", path);
    /* get type/size info on variable */
    status = adios_readutil_getdatainfo( gh, path, &ndims, dims, &type, msg, isvar);
    if (verbose) mexPrintf("Got info on var/attr %s: ndims=%d type=%s\n", path, ndims, adios_readutil_type_to_string(type));
    /* get type/size info on variable */
    if (status != 0) {
        mexErrMsgIdAndTxt("MATLAB:adiosreadc:read",msg);
    }
    if (verbose) { 
        mexPrintf("%s %s C dimensions %d [", (*isvar ? "variable" : "attribute"), path, ndims);
        for (i=0; i<ndims; i++)
            mexPrintf("%d ", dims[i]);
        mexPrintf("]\n");
    }

    /* Flip dimensions from ADIOS-read-api/C/row-major order to Matlab/Fortran/column-major order */
    swap_order(ndims, dims);

    /* extend offsets/counts if needed, change from 1..N to 0..N-1 indexing and
       interpret negative indices */
    recalc_offsets( ndims, dims, in_noffsets, in_offsets, in_counts, offsets, counts);

    /* interpret negative timestep */
    timestep = adios_readutil_calctimestep(fh, in_timestep);
    if (verbose) mexPrintf("Adjusted timestep: %d\n", timestep );

    /* create Matlab array with the appropriate type and size */
    if (type == adios_string) {
        /* Matlab string != char array of C, so handle separately */
        data = (void *) mxCalloc(dims[0], sizeof(char));
    } else {
        if (verbose) { 
            mexPrintf("Create %d-D Matlab array [", ndims);
            for (i=0; i<ndims; i++)
                mexPrintf("%d ", counts[i]);
            mexPrintf("]\n");
        }
        out = createMatlabArray( type, ndims, counts); 
        data = (void *) mxGetData(out);
    }

    /* Flip offsets/counts from Matlab/Fortran/column-major order to ADIOS-read-api/C/row-major order */
    swap_order(ndims, offsets);
    swap_order(ndims, counts);

    /* read in data */
    if (verbose) mexPrintf("Read in data\n");
    
    status = adios_readutil_readdata( gh, path, *isvar, offsets, counts, timestep, data, msg);
    if (status != 0) {
        mexErrMsgIdAndTxt("MATLAB:adiosreadc:read",msg);
    }

    if (type == adios_string) {
        out = mxCreateString( (char *)data);
        mxFree(data);
    }

    return out;
}
    

void errorCheck(int nlhs, int nrhs, const mxArray *prhs[]){
    /* Assume that we are called from adiosread.m which checks the arguments already */
    /* Check for proper number of arguments. */
    
    if ( nrhs != 7 ) {
        mexErrMsgIdAndTxt("MATLAB:adiosreadc:rhs","This function needs exactly 7 arguments: File, Group, Varpath, Time, Offsets, Counts, Verbose");
    }
    
    if ( !mxIsChar(prhs[0]) && !mxIsInt64(prhs[0]) ) {
        mexErrMsgIdAndTxt("MATLAB:adiosreadc:rhs","First arg must be either a string or an int64 handler.");
    } 
    
    if ( nlhs > 2 ) {
        mexErrMsgIdAndTxt("MATLAB:adiosreadc:lhs","Too many output arguments.");
    }
    
#if !defined(MX_COMPAT_32)
    /* Make sure that it is safe to cast dim to mwSize when using largeArrayDims.*/
    if ( dim > MWSIZE_MAX ) {
        mexErrMsgIdAndTxt("MATLAB:adiosreadc:dimensionTooLarge",
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

/** return the appropriate class for an adios type */
mxClassID adiostypeToMatlabClass(int type, int *isComplex) 
{
    *isComplex = 0;
    switch( type ) {
        case adios_unsigned_byte:
            return mxUINT8_CLASS;
        case adios_byte:
            return mxINT8_CLASS;
        case adios_string:
            return mxCHAR_CLASS;
               
        case adios_unsigned_short:
            return mxUINT16_CLASS;
        case adios_short:
            return mxINT16_CLASS;

        case adios_unsigned_integer:
            return mxUINT32_CLASS;
        case adios_integer:
            return mxINT32_CLASS;

        case adios_unsigned_long:
            return mxUINT64_CLASS;
        case adios_long:
            return mxINT64_CLASS;

        case adios_real: 
            return mxSINGLE_CLASS;
        case adios_double:
            return mxDOUBLE_CLASS;
             
        case adios_complex:     /* 8 bytes */
            *isComplex = 1;
            return mxSINGLE_CLASS;
        case adios_double_complex: /*  16 bytes */
            *isComplex = 2;
            return mxDOUBLE_CLASS;

        case adios_long_double: /* 16 bytes */
        default:
            mexErrMsgIdAndTxt("MATLAB:adiosreadc.c:dimensionTooLarge",
                 "Adios type %d (%s) not supported in visit.\n",
                 type, adios_readutil_type_to_string(type));
            break;
    }
    return 0; /* just to avoid warnings. never executed */
}

/* create an N-dim array with given type and dimensions */
mxArray* createMatlabArray( int adiostype, int ndims, int *dims)
{
    mxClassID mxtype;
    mxArray *arr;
    mwSize  mNdims;
    mwSize  mDims[16];
    mxComplexity ComplexFlag;
    int     i;
    int     isTypeComplex;

    /* convert ints to mwSizes */
    if (ndims > 0) {
        mNdims = (mwSize) ndims;
        for (i=0; i<ndims; i++) 
            mDims[i] = (mwSize)dims[i];
    } else {
        /* 0 dim: scalar value -> 1-by-1 Matlab array */
        mNdims = 2;
        mDims[0] = 1;
        mDims[1] = 1;
    }

    /* get type */
    mxtype = adiostypeToMatlabClass(adiostype, &isTypeComplex);
    if (isTypeComplex) ComplexFlag = mxCOMPLEX;
    else               ComplexFlag = mxREAL;

    /* create array */
    arr = mxCreateNumericArray( mNdims, mDims, mxtype, ComplexFlag);

    if (verbose) mexPrintf("Array for adios type %s is created\n", adios_readutil_type_to_string(adiostype));

    return arr;
}


/** - extend offset/count arrays to ndims if needed 
    - recalculate "from 1" Matlab indices to "from 0" C indices
    - recalculate the negative indices 
    !!! Provide the output arrays in the caller !!!
*/
void recalc_offsets( int ndims, int *dims, mwSize in_noffsets, 
                     const int32_t *in_offsets, const int32_t *in_counts,
                     int *offsets, int *counts)
{
    int i;
    for (i=0; i<ndims; i++) {
        if ((mwSize)i < in_noffsets) {
            if (in_offsets[i] < 0) /* negative offset means last-|offset| */
                offsets[i] = dims[i] + (int) in_offsets[i];
            else 
                offsets[i] = (int) in_offsets[i] - 1; /* C index start from 0 */

            if (in_counts[i] < 0) /* negative count means last-|count|+1-start */
                counts[i]  = dims[i] + (int) in_counts[i] - offsets[i] + 1;
            else
                counts[i]  = (int) in_counts[i]; 
        } else {
            /* extend offset/count array to match variable's dimensions */
            if (verbose) mexPrintf("Extend offset/counts for dim %d: offset=%d count=%d\n", i, 0, dims[i] );
            offsets[i] = 0;
            counts[i]  = dims[i];
        }

    }
}

/* Reverse the order in an array in place.
   use swapping from Matlab/Fortran/column-major order to ADIOS-read-api/C/row-major order and back
*/
void swap_order(int n, int *array)
{
    int i, tmp;
    for (i=0; i<n/2; i++) {
        tmp = array[i];
        array[i] = array[n-1-i];
        array[n-1-i] = tmp;
    }
}




