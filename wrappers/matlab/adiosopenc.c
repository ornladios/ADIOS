/* 
 * ADIOS is freely available under the terms of the BSD license described
 * in the COPYING file in the top level directory of this source distribution.
 *
 * Copyright (c) 2008 - 2009.  UT-BATTELLE, LLC. All rights reserved.
 */

/*=================================================================
 * adiosopenc.c - Open an ADIOS file
 *
 * Input: File, Verbose
 *    File:     string (name) or int64 (handler)
 *    Verbose:  numeric (double)
 *
 * Output: Information structure
 *     Name        file path
 *     FileHandler int64 file handler
 *     TimeStart   First timestep in file (usually = 1)
 *     TimeSteps   Number of timesteps in file, always at least 1
 *     Groups      Adios groups in the file. Usually 1 group is in a file.
 *                 This is a structure array of 
 *
 *        Name          group name
 *        GroupHandler  int64 group handler
 *        Variables     structure array of variables
 *           Name          path of variable
 *           Type          Matlab type class of data
 *           Dims          Array of dimensions
 *           Timedim       The time dimension, 0 if there is no time varying 
 *                         part of the variable
 *           GlobalMin     global minimum  of the variable (1-by-1 mxArray)
 *           GlobalMax     global maximum of the variable
 *           
 *        Attributes  structure array of attributes
 *           Name          path of attribute
 *           Type          Matlab type class of data
 *           Value         attribute value (mxArray)
 *
 *
 *
 * $Revision: 1.0 $  $Date: 2009/08/05 12:53:41 $
 * Author: Norbert Podhorszki <pnorbert@ornl.gov>
 *=================================================================*/

#include <string.h>              /* memcpy */
#include "mex.h"
#include "adios_types.h"
#include "adios_read.h"
#include "adios_types.h"

static int verbose=0;

mxClassID adiostypeToMatlabClass(int type, mxComplexity *complexity );
mxArray* valueToMatlabValue( void * data, mxClassID mxtype, mxComplexity complexFlag);
void errorCheck(int nlhs, int nrhs, const mxArray *prhs[]);
char* getString(const mxArray *mxstr);
static void swap_order(int n, uint64_t *array);


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    char *fname;                                 /* file name */
    int status;
    char msg[512];                               /* error messages from function calls */
    int32_t *int32p;

    ADIOS_FILE *fp;                             /* File handler structure */
    ADIOS_GROUP *gp;                            /* Group handler structure */
    ADIOS_VARINFO *vinfo;                       /* Variable info structure */
    int mpi_comm_dummy;                         /* ADIOS read API needs an MPI communicator */
    int32_t timedim;                            /* time dimension is translated C 0.. -> Matlab 1.. */
    int asize;                                  /* Attribute size info */
    enum ADIOS_DATATYPES adiostype;             /* Attribute type info */
    void *data;                                 /* Attributes return their values */
    mwSize mDims[] = {1};                       /* dimensions to create a scalar mxArray */  

    int gi,vi,ai, i;                            /* loop variables for groups, vars and attrs */
    mxArray *arr;                               /* temp array for constructions */
    mxArray *groups, *vars, *attrs;             /* struct arrays under top struct */
    mxClassID mxtype;                           /* matlab type (of vars and attrs) */
    mxComplexity complexFlag; 

    /* Output structure definition */
    const char *top_field_names[] = {"Name", "FileHandler", "TimeStart", "TimeSteps", "Groups"}; /* top level struct fields */
    mwSize ntopfields = 5;
    mwSize top_struct_dims[] = {1,1};   /* dimensions for top level struct array: 1-by-1 */
    int top_field_Name;
    int top_field_FileHandler;
    int top_field_TimeStart;
    int top_field_TimeSteps;
    int top_field_Groups;

    const char *group_field_names[] = {"Name", "GroupHandler", "Variables", "Attributes"}; /* group level struct fields */
    mwSize ngroupfields = 4;
    mwSize group_struct_dims[2];        /* dimensions for group level struct array: 1-by-sth */
    int group_field_Name;
    int group_field_GroupHandler;
    int group_field_Variables;
    int group_field_Attributes;

    const char *var_field_names[] = {"Name","Type","Dims", "Timedim", "GlobalMin", "GlobalMax"}; /* variable level struct fields */
    mwSize nvarfields = 6;
    mwSize var_struct_dims[2];        /* dimensions for variable level struct array: 1-by-sth */
    int var_field_Name;
    int var_field_Type;
    int var_field_Dims;
    int var_field_Timedim;
    int var_field_GlobalMin;
    int var_field_GlobalMax;

    const char *attr_field_names[] = {"Name","Type","Value"}; /* attribute level struct fields */
    mwSize nattrfields = 3;
    mwSize attr_struct_dims[2];        /* dimensions for attribute level struct array: 1-by-sth */
    int attr_field_Name;
    int attr_field_Type;
    int attr_field_Value;

    errorCheck(nlhs, nrhs, prhs);

    /*mexPrintf("nrhs=%d  nlhs=%d\n", nrhs, nlhs);*/
    
    /***********************/
    /* 0. get verbose parameter first */
    verbose = (int)mxGetScalar(prhs[1]);
    if (verbose) mexPrintf("Verbosity level: %d\n", verbose);

    /* 1. get file handler */
    if (mxIsChar(prhs[0])) {
        fname = getString( (mxArray *)prhs[0]);
        if (verbose) mexPrintf("File name: \"%s\"\n", fname);
    } 

    /**************************************/
    /* Open ADIOS file now and get groups */
    fp = adios_fopen (fname, mpi_comm_dummy);
    if (fp == NULL) {
       mexErrMsgIdAndTxt("MATLAB:adiosopenc:open",adios_errmsg());
    }
    if (verbose) mexPrintf("Number of adios groups: %d, fp=%lld\n", fp->groups_count, (int64_t) fp);

    /******************************/
    /* Create top level structure */
    if (verbose) mexPrintf("Create top struct array, 1-by-1\n");
    plhs[0] = mxCreateStructArray(2, top_struct_dims, ntopfields, top_field_names);
    top_field_Name = mxGetFieldNumber(plhs[0],"Name");
    top_field_FileHandler = mxGetFieldNumber(plhs[0],"FileHandler");
    top_field_TimeStart = mxGetFieldNumber(plhs[0],"TimeStart");
    top_field_TimeSteps = mxGetFieldNumber(plhs[0],"TimeSteps");
    top_field_Groups = mxGetFieldNumber(plhs[0],"Groups");
    mxSetFieldByNumber(plhs[0],0,top_field_Name,mxCreateString(fname));
    arr = valueToMatlabValue(&fp, mxINT64_CLASS, mxREAL);
    mxSetFieldByNumber(plhs[0],0,top_field_FileHandler,arr);
    arr = mxCreateNumericMatrix( 1, 1, mxINT32_CLASS, mxREAL);
    *(int32_t *)mxGetData(arr) = fp->tidx_start;
    mxSetFieldByNumber(plhs[0],0,top_field_TimeStart,arr);
    arr = mxCreateNumericMatrix( 1, 1, mxINT32_CLASS, mxREAL);
    *(int32_t *)mxGetData(arr) = fp->ntimesteps;
    mxSetFieldByNumber(plhs[0],0,top_field_TimeSteps,arr);
    /* Create top.Groups structure array */
    if (verbose) mexPrintf("Create top.Groups struct array, 1-by-%d\n",fp->groups_count);
    group_struct_dims[0] = 1;
    group_struct_dims[1] = fp->groups_count;
    groups = mxCreateStructArray(2, group_struct_dims, ngroupfields, group_field_names);
    mxSetFieldByNumber(plhs[0],0,top_field_Groups,groups);


    /****************************/
    /* Fill in Groups structure */
    group_field_Name = mxGetFieldNumber(groups,"Name");
    group_field_GroupHandler = mxGetFieldNumber(groups,"GroupHandler");
    group_field_Variables = mxGetFieldNumber(groups,"Variables");
    group_field_Attributes = mxGetFieldNumber(groups,"Attributes");

    for (gi = 0; gi < fp->groups_count; gi++) {
        /* Get info of one group: handler, list of vars, list of attrs */
        if (verbose) mexPrintf("Group %s: get info\n", fp->group_namelist[gi]);
        gp = adios_gopen_byid(fp, gi);
        if (gp == NULL) {
           mexErrMsgIdAndTxt("MATLAB:adiosopenc:groupinfo",adios_errmsg());
        }
        if (verbose) mexPrintf("    %d variables and %d attributes\n", gp->vars_count, gp->attrs_count);
        /* Group fields for gi-th group */
        mxSetFieldByNumber(groups,gi,group_field_Name,mxCreateString(fp->group_namelist[gi]));
        mexPrintf("Group gp=%lld id=%d vcnt=%d\n", (int64_t)gp, gp->grpid, gp->vars_count);
        arr = valueToMatlabValue(&gp, mxINT64_CLASS, mxREAL);
        mxSetFieldByNumber(groups,gi,group_field_GroupHandler,arr);
        /* Create top.Groups.Variables structure array */
        var_struct_dims[0] = 1;
        var_struct_dims[1] = gp->vars_count;
        vars = mxCreateStructArray(2, var_struct_dims, nvarfields, var_field_names);
        mxSetFieldByNumber(groups,gi,group_field_Variables,vars);
        var_field_Name = mxGetFieldNumber(vars,"Name");
        var_field_Type = mxGetFieldNumber(vars,"Type");
        var_field_Dims = mxGetFieldNumber(vars,"Dims");
        var_field_Timedim = mxGetFieldNumber(vars,"Timedim");
        var_field_GlobalMin = mxGetFieldNumber(vars,"GlobalMin");
        var_field_GlobalMax = mxGetFieldNumber(vars,"GlobalMax");
        /* Create top.Groups.Attributes structure array */
        attr_struct_dims[0] = 1;
        attr_struct_dims[1] = gp->attrs_count;
        attrs = mxCreateStructArray(2, attr_struct_dims, nattrfields, attr_field_names);
        mxSetFieldByNumber(groups,gi,group_field_Attributes,attrs);
        attr_field_Name = mxGetFieldNumber(attrs,"Name");
        attr_field_Type = mxGetFieldNumber(attrs,"Type");
        attr_field_Value = mxGetFieldNumber(attrs,"Value");

        /******************************/
        /* Add variables to the group */

        if (verbose>1) mexPrintf("    Variables\n");
        for (vi=0; vi < gp->vars_count; vi++) {
            vinfo = adios_inq_var_byid( gp, vi);
            if (vinfo == NULL) {
               mexErrMsgIdAndTxt("MATLAB:adiosopenc:varinfo",adios_errmsg());
            }
            /* Flip dimensions from ADIOS-read-api/C/row-major order to Matlab/Fortran/column-major order */
            swap_order(vinfo->ndim, vinfo->dims);

            if (verbose>1) {
                mexPrintf("      %s: ndims=%d, adios type=%s, timedim=%d dimensions [", 
                gp->var_namelist[vi], vinfo->ndim, adios_type_to_string(vinfo->type), vinfo->timedim);
                for (i=0; i<vinfo->ndim; i++)
                    mexPrintf("%lld ", vinfo->dims[i]);
                mexPrintf("]\n");
            }
            /* field NAME */
            mxSetFieldByNumber(vars,vi,var_field_Name,mxCreateString(gp->var_namelist[vi]));
            /* field TYPE */
            mxtype = adiostypeToMatlabClass(vinfo->type, &complexFlag);
            arr = mxCreateNumericMatrix( 1, 1, mxtype, complexFlag);
            mxSetFieldByNumber(vars,vi,var_field_Type,mxCreateString(mxGetClassName(arr)));
            mxDestroyArray(arr);
            /* field DIMS */
            if (vinfo->ndim > 0) {
                arr = mxCreateNumericMatrix( 1, vinfo->ndim, mxINT32_CLASS, mxREAL);
                int32p = (int32_t *)mxGetData(arr);
                for (i=0; i<vinfo->ndim; i++) 
                    int32p[i] = (int32_t) vinfo->dims[i];
            } else {
                arr = mxCreateNumericMatrix( 0, 0, mxINT32_CLASS, mxREAL);
            }
            mxSetFieldByNumber(vars,vi,var_field_Dims,arr);

            /* field TIMEDIM */
            /* Timedim is -1,0...ndim-1 in C, 0,1..ndim in Matlab */
            timedim = vinfo->timedim + 1;
            arr = valueToMatlabValue((void *)(&timedim), mxINT32_CLASS, mxREAL);
            mxSetFieldByNumber(vars,vi,var_field_Timedim,arr);

            /* field GLOBALMIN */
            arr = valueToMatlabValue(vinfo->gmin, mxtype, complexFlag);
            mxSetFieldByNumber(vars,vi,var_field_GlobalMin,arr);

            /* field GLOBALMAX */
            arr = valueToMatlabValue(vinfo->gmax, mxtype, complexFlag);
            mxSetFieldByNumber(vars,vi,var_field_GlobalMax,arr);

            adios_free_varinfo(vinfo);
        }



        /******************************/
        /* Add attributes to the group */

        if (verbose>1) mexPrintf("    Attributes\n");
        for (ai=0; ai < gp->attrs_count; ai++) {
            status = adios_get_attr_byid( gp, ai, &adiostype, &asize, &data);
            if (status != 0) {
               mexErrMsgIdAndTxt("MATLAB:adiosopenc:varinfo",adios_errmsg());
            }
            if (verbose>1) 
                mexPrintf("      %s: adios type=%s, size=%d\n", 
                gp->attr_namelist[ai], adios_type_to_string(adiostype), asize);
            /* field NAME */
            mxSetFieldByNumber(attrs,ai,attr_field_Name,mxCreateString(gp->attr_namelist[ai]));
            /* field TYPE */
            mxtype = adiostypeToMatlabClass(adiostype, &complexFlag);
            arr = mxCreateNumericMatrix( 1, 1, mxtype, complexFlag);
            mxSetFieldByNumber(attrs,ai,attr_field_Type,mxCreateString(mxGetClassName(arr)));
            mxDestroyArray(arr);
            /* field VALUE */
            arr = valueToMatlabValue(data, mxtype, complexFlag);
            mxSetFieldByNumber(attrs,ai,attr_field_Value,arr);

            free(data); /* we do not store attribute values yet */

        }

        if (verbose>1) mexPrintf("    finished defining group\n");
    }


    if (verbose) mexPrintf("return from adiosopenc\n");
}

mxArray * valueToMatlabValue( void * data, mxClassID mxtype, mxComplexity complexFlag)
{
    /* copies values in all cases, so one can free(data) later */
    mxArray *arr;
    if (data == NULL) {
        arr = mxCreateString("undefined");
    } else if (mxtype == mxCHAR_CLASS) {
        arr = mxCreateString((char *)data);
    } else if (complexFlag == mxCOMPLEX) {
        arr = mxCreateDoubleMatrix( 1, 1, mxCOMPLEX);
        if (mxtype == mxSINGLE_CLASS) {
            *(double *)mxGetPr(arr) = ((float *)data)[0];
            *(double *)mxGetPi(arr) = ((float *)data)[1];
        } else {
            *(double *)mxGetPr(arr) = ((double *)data)[0];
            *(double *)mxGetPi(arr) = ((double *)data)[1];
        }
    } else {
        arr = mxCreateNumericMatrix( 1, 1, mxtype, mxREAL);
        memcpy( mxGetData(arr), data, mxGetElementSize(arr));
    }
    return arr;
}

void errorCheck(int nlhs, int nrhs, const mxArray *prhs[]){
    /* Assume that we are called from adiosread.m which checks the arguments already */
    /* Check for proper number of arguments. */
    
    if ( nrhs != 2 ) {
        mexErrMsgIdAndTxt("MATLAB:adiosopenc:rhs","This function needs exactly 2 arguments: File, Verbose");
    }
    
    if ( !mxIsChar(prhs[0]) ) {
        mexErrMsgIdAndTxt("MATLAB:adiosopenc:rhs","First arg must be a string.");
    } 
    
    if ( !mxIsNumeric(prhs[1]) ) {
        mexErrMsgIdAndTxt("MATLAB:adiosopenc:rhs","Second arg must be a number.");
    } 
    
    if ( nlhs > 1 ) {
        mexErrMsgIdAndTxt("MATLAB:adiosopenc:lhs","Too many output arguments.");
    }
    
#if !defined(MX_COMPAT_32)
    /* Make sure that it is safe to cast dim to mwSize when using largeArrayDims.*/
    if ( dim > MWSIZE_MAX ) {
        mexErrMsgIdAndTxt("MATLAB:adiosopenc:dimensionTooLarge",
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

/** return the appropriate class for an adios type (and complexity too) */
mxClassID adiostypeToMatlabClass(enum ADIOS_DATATYPES type, mxComplexity *complexity ) 
{
    *complexity = mxREAL;
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
            *complexity = mxCOMPLEX;
            return mxSINGLE_CLASS;
        case adios_double_complex: /*  16 bytes */
            *complexity = mxCOMPLEX;
            return mxDOUBLE_CLASS;

        case adios_long_double: /* 16 bytes */
        default:
            mexErrMsgIdAndTxt("MATLAB:adiosopenc.c:dimensionTooLarge",
                 "Adios type %d (%s) not supported in visit.\n",
                 type, adios_type_to_string(type));
            break;
    }
    return 0; /* just to avoid warnings. never executed */
}


/* Reverse the order in an array in place.
   use swapping from Matlab/Fortran/column-major order to ADIOS-read-api/C/row-major order and back
*/
static void swap_order(int n, uint64_t *array)
{
    int i, tmp;
    for (i=0; i<n/2; i++) {
        tmp = array[i];
        array[i] = array[n-1-i];
        array[n-1-i] = tmp;
    }
}


