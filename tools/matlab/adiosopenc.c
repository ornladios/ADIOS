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
 *     TimeStart   First timestep in file
 *     TimeEnd     Last timestep in file
 *     Groups      Adios groups in the file. Usually 1 group is in a file.
 *                 This is a structure array of 
 *
 *        Name          group name
 *        GroupHandler  int64 group handler
 *        Variables     structure array of variables
 *           Name          path of variable
 *           Type          Matlab type class of data
 *           Dims          Array of dimensions
 *           Timed         BOOL, true: several timesteps are stored in file. 
 *           
 *        Attributes  structure array of attributes
 *           Name          path of attribute
 *           Type          Matlab type class of data
 *           Size          1 or length of string in case of string attributes.
 *
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

mxClassID adiostypeToMatlabClass(int type, mxComplexity *complexity );
void errorCheck(int nlhs, int nrhs, const mxArray *prhs[]);
char* getString(const mxArray *mxstr);


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    char *fname;                                 /* file name */
    int status;
    char msg[512];                               /* error messages from function calls */
    int64_t fh, gh, *int64p;                     /* file and group handlers and temp pointers */
    int32_t *int32p;

    /* results from adios_readutil_open */
    int  ngroups;                                /* number of adios groups */
    adios_readutil_namelist *groupnames, *gname; /* linked list of group names */
    int32_t tstart, tend;                        /* first and last timesteps */

    /* results from adios_readutil_groupinfo */
    int  nvars;                                  /* number of adios groups */
    adios_readutil_namelist *varnames, *vname;   /* linked list of variable names */
    int  nattrs;                                 /* number of adios groups */
    adios_readutil_namelist *attrnames, *aname;  /* linked list of attribute names */

    /* results from adios_readutil_varinfo and _attrinfo */
    int  ndims, dims[16];                        /* size info for variable/attr read in */
    int  adiostype, timed;                       /* type of variable/attr read in */
    int  asize;                                  /* size of attribute (1 or string length) */

    int gi,vi,ai, i;                            /* loop variables for groups, vars and attrs */
    mxArray *arr;                               /* temp array for constructions */
    mxArray *groups, *vars, *attrs;             /* struct arrays under top struct */
    mxClassID mxtype;                           /* matlab type (of vars and attrs) */
    mxComplexity complexFlag; 

    /* Output structure definition */
    const char *top_field_names[] = {"Name", "FileHandler", "TimeStart", "TimeEnd", "Groups"}; /* top level struct fields */
    mwSize ntopfields = 5;
    mwSize top_struct_dims[] = {1,1};   /* dimensions for top level struct array: 1-by-1 */
    int top_field_Name;
    int top_field_FileHandler;
    int top_field_TimeStart;
    int top_field_TimeEnd;
    int top_field_Groups;

    const char *group_field_names[] = {"Name", "GroupHandler", "Variables", "Attributes"}; /* group level struct fields */
    mwSize ngroupfields = 4;
    mwSize group_struct_dims[2];        /* dimensions for group level struct array: 1-by-sth */
    int group_field_Name;
    int group_field_GroupHandler;
    int group_field_Variables;
    int group_field_Attributes;

    const char *var_field_names[] = {"Name","Type","Dims", "Timed"}; /* variable level struct fields */
    mwSize nvarfields = 4;
    mwSize var_struct_dims[2];        /* dimensions for variable level struct array: 1-by-sth */
    int var_field_Name;
    int var_field_Type;
    int var_field_Dims;
    int var_field_Timed;

    const char *attr_field_names[] = {"Name","Type","Size"}; /* attribute level struct fields */
    mwSize nattrfields = 3;
    mwSize attr_struct_dims[2];        /* dimensions for attribute level struct array: 1-by-sth */
    int attr_field_Name;
    int attr_field_Type;
    int attr_field_Size;

    errorCheck(nlhs, nrhs, prhs);

    /*mexPrintf("nrhs=%d  nlhs=%d\n", nrhs, nlhs);*/
    
    /***********************/
    /* 0. get verbose parameter first */
    verbose = (int)mxGetScalar(prhs[1]);
    if (verbose) mexPrintf("Verbosity level: %d\n", verbose);
    adios_readutil_setverbosity(verbose);

    /* 1. get file handler */
    if (mxIsChar(prhs[0])) {
        fname = getString( (mxArray *)prhs[0]);
        if (verbose) mexPrintf("File name: \"%s\"\n", fname);
    } 

    /**************************************/
    /* Open ADIOS file now and get groups */
    status = adios_readutil_open (fname, &fh, &ngroups, &groupnames, &tstart, &tend, msg); 
    if (status != 0) {
       mexErrMsgIdAndTxt("MATLAB:adiosopenc:open",msg);
    }
    if (verbose) mexPrintf("Number of adios groups: %d, fh=%lld\n", ngroups, fh);

    /******************************/
    /* Create top level structure */
    if (verbose) mexPrintf("Create top struct array, 1-by-1\n");
    plhs[0] = mxCreateStructArray(2, top_struct_dims, ntopfields, top_field_names);
    top_field_Name = mxGetFieldNumber(plhs[0],"Name");
    top_field_FileHandler = mxGetFieldNumber(plhs[0],"FileHandler");
    top_field_TimeStart = mxGetFieldNumber(plhs[0],"TimeStart");
    top_field_TimeEnd = mxGetFieldNumber(plhs[0],"TimeEnd");
    top_field_Groups = mxGetFieldNumber(plhs[0],"Groups");
    mxSetFieldByNumber(plhs[0],0,top_field_Name,mxCreateString(fname));
    arr = mxCreateNumericMatrix( 1, 1, mxINT64_CLASS, mxREAL);
    *(int64_t *)mxGetData(arr) = fh;
    mxSetFieldByNumber(plhs[0],0,top_field_FileHandler,arr);
    arr = mxCreateNumericMatrix( 1, 1, mxINT32_CLASS, mxREAL);
    *(int32_t *)mxGetData(arr) = tstart;
    mxSetFieldByNumber(plhs[0],0,top_field_TimeStart,arr);
    arr = mxCreateNumericMatrix( 1, 1, mxINT32_CLASS, mxREAL);
    *(int32_t *)mxGetData(arr) = tend;
    mxSetFieldByNumber(plhs[0],0,top_field_TimeEnd,arr);
    /* Create top.Groups structure array */
    if (verbose) mexPrintf("Create top.Groups struct array, 1-by-%d\n",ngroups);
    group_struct_dims[0] = 1;
    group_struct_dims[1] = ngroups;
    groups = mxCreateStructArray(2, group_struct_dims, ngroupfields, group_field_names);
    mxSetFieldByNumber(plhs[0],0,top_field_Groups,groups);


    /****************************/
    /* Fill in Groups structure */
    group_field_Name = mxGetFieldNumber(groups,"Name");
    group_field_GroupHandler = mxGetFieldNumber(groups,"GroupHandler");
    group_field_Variables = mxGetFieldNumber(groups,"Variables");
    group_field_Attributes = mxGetFieldNumber(groups,"Attributes");

    gname=groupnames;
    gi = 0;
    while (gname != NULL) {
        /* Get info of one group: handler, list of vars, list of attrs */
        if (verbose) mexPrintf("Group %s: get info\n", gname->name);
        status = adios_readutil_groupinfo (fh, gname->name, &gh, &nvars, &varnames, &nattrs, &attrnames, msg); 
        if (status != 0) {
           mexErrMsgIdAndTxt("MATLAB:adiosopenc:groupinfo",msg);
        }
        if (verbose) mexPrintf("    %d variables and %d attributes\n", nvars, nattrs);
        /* Group fields for gi-th group */
        mxSetFieldByNumber(groups,gi,group_field_Name,mxCreateString(gname->name));
        arr = mxCreateNumericMatrix( 1, 1, mxINT64_CLASS, mxREAL);
        *(int64_t *)mxGetData(arr) = gh;
        mxSetFieldByNumber(groups,gi,group_field_GroupHandler,arr);
        /* Create top.Groups.Variables structure array */
        var_struct_dims[0] = 1;
        var_struct_dims[1] = nvars;
        vars = mxCreateStructArray(2, var_struct_dims, nvarfields, var_field_names);
        mxSetFieldByNumber(groups,gi,group_field_Variables,vars);
        var_field_Name = mxGetFieldNumber(vars,"Name");
        var_field_Type = mxGetFieldNumber(vars,"Type");
        var_field_Dims = mxGetFieldNumber(vars,"Dims");
        var_field_Timed = mxGetFieldNumber(vars,"Timed");
        /* Create top.Groups.Attributes structure array */
        attr_struct_dims[0] = 1;
        attr_struct_dims[1] = nattrs;
        attrs = mxCreateStructArray(2, attr_struct_dims, nattrfields, attr_field_names);
        mxSetFieldByNumber(groups,gi,group_field_Attributes,attrs);
        attr_field_Name = mxGetFieldNumber(attrs,"Name");
        attr_field_Type = mxGetFieldNumber(attrs,"Type");
        attr_field_Size = mxGetFieldNumber(attrs,"Size");

        /******************************/
        /* Add variables to the group */

        vname = varnames;
        vi = 0;
        if (verbose>1) mexPrintf("    Variables\n");
        while (vname != NULL) {
            status = adios_readutil_getvarinfo( gh, vname->name, &ndims, dims, &adiostype, &timed, msg);
            if (status != 0) {
               mexErrMsgIdAndTxt("MATLAB:adiosopenc:varinfo",msg);
            }
            if (verbose>1) mexPrintf("      %s: ndims=%d, adios type=%s, timed=%s\n", vname->name, ndims, adios_readutil_type_to_string(adiostype), (timed ? "true" : "false" ));
            /* field NAME */
            mxSetFieldByNumber(vars,vi,var_field_Name,mxCreateString(vname->name));
            /* field TYPE */
            mxtype = adiostypeToMatlabClass(adiostype, &complexFlag);
            arr = mxCreateNumericMatrix( 1, 1, mxtype, complexFlag);
            mxSetFieldByNumber(vars,vi,var_field_Type,mxCreateString(mxGetClassName(arr)));
            mxDestroyArray(arr);
            /* field DIMS */
            if (ndims > 0) {
                arr = mxCreateNumericMatrix( 1, ndims, mxINT32_CLASS, mxREAL);
                int32p = (int32_t *)mxGetData(arr);
                for (i=0; i<ndims; i++) 
                    int32p[i] = (int32_t) dims[i];
            } else {
                arr = mxCreateNumericMatrix( 0, 0, mxINT32_CLASS, mxREAL);
            }

            mxSetFieldByNumber(vars,vi,var_field_Dims,arr);
            /* field TIMED */
            mxSetFieldByNumber(vars,vi,var_field_Timed,mxCreateLogicalScalar(timed));

            vi++;
            vname = vname->next;
        }



        /******************************/
        /* Add attributes to the group */

        aname = attrnames;
        ai = 0;
        if (verbose>1) mexPrintf("    Attributes\n");
        while (aname != NULL) {
            status = adios_readutil_getattrinfo( gh, aname->name, &adiostype, &asize, msg);
            if (status != 0) {
               mexErrMsgIdAndTxt("MATLAB:adiosopenc:varinfo",msg);
            }
            if (verbose>1) mexPrintf("      %s: adios type=%s, size=%d\n", aname->name, adios_readutil_type_to_string(adiostype), asize);
            /* field NAME */
            mxSetFieldByNumber(attrs,ai,attr_field_Name,mxCreateString(aname->name));
            /* field TYPE */
            mxtype = adiostypeToMatlabClass(adiostype, &complexFlag);
            arr = mxCreateNumericMatrix( 1, 1, mxtype, complexFlag);
            mxSetFieldByNumber(attrs,ai,attr_field_Type,mxCreateString(mxGetClassName(arr)));
            mxDestroyArray(arr);
            /* field SIZE */
            mxSetFieldByNumber(attrs,ai,attr_field_Size,mxCreateDoubleScalar((double)asize));

            ai++;
            aname = aname->next;
        }




        if (verbose>1) mexPrintf("    finished defining group\n");
        gi++;
        gname=gname->next;
    }


    adios_readutil_freenamelist(groupnames);
    adios_readutil_freenamelist(varnames);
    adios_readutil_freenamelist(attrnames);
    if (verbose) mexPrintf("return from adiosopenc\n");
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
mxClassID adiostypeToMatlabClass(int type, mxComplexity *complexity ) 
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
                 type, adios_readutil_type_to_string(type));
            break;
    }
    return 0; /* just to avoid warnings. never executed */
}




