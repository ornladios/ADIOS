/* List the content of a BP file.
 *
 * Copyright Oak Ridge National Laboratory 2009
 * Author: Norbert Podhorszki, pnorbert@ornl.gov
 *
 * TODO: -S handle int8[] as string
**/

#ifndef _GNU_SOURCE
#   define _GNU_SOURCE
#endif

#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <unistd.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <getopt.h>
#include <errno.h>
#include <limits.h>   // LONG_MAX
#include <math.h>     // NAN
#include <libgen.h>   // basename
#include <regex.h>    // regular expression matching
#include <fnmatch.h>  // shell pattern matching

#include "bpls.h"
#include "adios_read.h"
#include "adios_types.h"

// global variables 
// Values from the arguments or defaults
char *outpath;              // output files' starting path (can be extended with subdirs, names, indexes)
char *varmask[MAX_MASKS];   // can have many -var masks (either shell patterns or extended regular expressions)
char *grpmask;              // list only groups matching the mask
int  nmasks;                // number of masks specified
char *vfile;                // file containing variable (values) to plot
char *start;                // dimension spec starting points 
char *count;                // dimension spec counts
char format[32];            // format string for one data element (e.g. %6.2f)
bool formatgiven;           // true if format string is provided as argument

// Flags from arguments or defaults
bool dump;                 // dump data not just list info
bool output_xml;
bool use_regexp;           // use varmasks as regular expressions
bool sortnames;            // sort names before listing
bool listattrs;            // do list attributes too
bool readattrs;            // also read all attributes and print
bool readscalars;          // read all scalar variables and print
bool noindex;              // do no print array indices with data

// other global variables
char *prgname; /* argv[0] */
//long timefrom, timeto;
int  istart[MAX_DIMS], icount[MAX_DIMS], ndimsspecified=0;
regex_t varregex[MAX_MASKS]; // compiled regular expressions of varmask
regex_t grpregex;            // compiled regular expressions of grpmask
int  ncols = 6; // how many values to print in one row (only for -p)
int  verbose = 0;
FILE *outf;   // file to print to or stdout

struct option options[] = {
    {"help",                 no_argument,          NULL,    'h'},
    {"verbose",              no_argument,          NULL,    'v'},
    {"dump",                 no_argument,          NULL,    'd'},
    {"group",                no_argument,          NULL,    'g'},
    {"regexp",               no_argument,          NULL,    'e'},
    {"output",               required_argument,    NULL,    'o'},
    {"xml",                  no_argument,          NULL,    'x'},
    {"start",                required_argument,    NULL,    's'}, 
    {"count",                required_argument,    NULL,    'c'}, 
    {"noindex",              no_argument,          NULL,    'y'},
    {"sort",                 no_argument,          NULL,    't'},
    {"attrs",                no_argument,          NULL,    'a'},
    {"scalars",              no_argument,          NULL,    'l'},
    {"columns",              required_argument,    NULL,    'n'}, 
    {"format",               required_argument,    NULL,    'f'}, 
//    {"time",                 required_argument,    NULL,    't'}, 
    {NULL,                   0,                    NULL,    0}
};


static const char *optstring = "hveytaldg:o:x:s:c:n:f:";

// help function
void display_help() {
   //printf( "Usage: %s  \n", prgname);
   printf("usage: bpls [OPTIONS] file [mask1 mask2 ...]\n"
        "\nList/dump content of a BP file. \n"
        "A mask can be a shell pattern like with 'ls' e.g. \"*/x?\".\n"
        "Variables with multiple timesteps are reported with an extra dimensions.\n"
        "The time dimension is the first dimension then.\n"
        "\n"
        "  --scalars | -l           Print values of all scalars\n"
        "  --attrs   | -a           List/match attributes too\n"
        "  --sort    | -t           Sort names before listing\n"
        "  --group   | -g <mask>    List/dump groups matching the mask only\n"
        "  --dump    | -d           Dump matched variables/attributes\n"
        "                             To match attributes too, add option -a\n"
        "  --regexp  | -e           Treat masks as extended regular expressions\n"
        "  --output  | -o <path>    Print to a file instead of stdout\n"
/*
        "  --xml    | -x            # print as xml instead of ascii text\n"
*/
        "  --start   | -s \"spec\"    Offset indices in each dimension \n"
        "                             (default is 0 for all dimensions) \n"
        "                             <0 is handled as in python (-1 is last)\n"
        "  --count   | -c \"spec\"    Number of elements in each dimension\n"
        "                             -1 denotes 'until end' of dimension\n"
        "                             (default is -1 for all dimensions)\n"
        "  --noindex | -y           Print data without array indices\n"
        "  --columns | -n \"cols\"    Number of data elements per row to print\n"
        "  --format  | -f \"str\"     Format string to use for one data item in print\n"
        "                             instead of the default. E.g. \"%%6.3f\"\n"
/*
        "  --time    | -t N [M]      # print data for timesteps N..M only (or only N)\n"
        "                              default is to print all available timesteps\n"
*/
        "\n"
        "  Examples for slicing:\n"
        "  -s \"0,0,0\"    -c \"1,99,1\":  print 100 elements (of the 2nd dimension).\n"
        "  -s \"0,0\"      -c \"1,-1\":    print the whole 2nd dimension however large it is.\n"
        "  -s \"-1,-1\"    -c \"1,1\":     print the very last element (of a 2D array)\n"
        "\n"
        "Help options\n"
        "  --help    | -h               Print this help.\n"
        "  --verbose | -v               Print log about what this program is doing.\n"
        "                                 Use multiple -v to increase logging level.\n"
        "Typical use: bpls -lat <file>\n"
        );
}

/** Main */
int main( int argc, char *argv[] ) {
    int retval = 0;
    int i, timearg=false; 
    long int tmp;

    init_globals();

    ////prgname = strdup(argv[0]);

    /* other variables */
    char c, last_c='_';
    int last_opt = -1;
    /* Process the arguments */
    while ((c = getopt_long(argc, argv, optstring, options, NULL)) != -1) {
        switch (c) {
                case 'a':
                    listattrs=true;;
                    break;
                case 'c':
                    count = strndup(optarg,256);
                    break;
                case 'd':
                    dump = true;
                    break;
                case 'e':
                    use_regexp = true;
                    break;
                case 'f':
                    snprintf(format, sizeof(format), "%s", optarg);
                    formatgiven = true;
                    break;
                case 'g':
                    grpmask = strndup(optarg,256);
                    break;
                case 'h':
                    display_help();
                    return 0;
                case 'l':
                    readscalars = true;
                    readattrs = true;
                    break;
                case 'n':
                    tmp = strtol(optarg, (char **)NULL, 0);
                    if (errno) {
                        fprintf(stderr, "Error: could not convert --columns value: %s\n", optarg);
                        return 1;
                    }
                    ncols=tmp;
                    break;
                case 'o':
                    outpath = strndup(optarg,256);
                    break;
                case 's':
                    start = strndup(optarg,256);
                    break;
                case 't':
                    sortnames = true;
                    break;
                case 'x':
                    output_xml = true;
                    break;
                case 'y':
                    noindex = true;
                    break;
                case 'v':
                    verbose++;
                    break;
                /*
                case 't':
                    tmp = strtol(optarg, (char **)NULL, 0);
                    if (errno) {
                        fprintf(stderr, "Error: could not convert --time value: %s\n", optarg);
                        return 1;
                    }
                    timefrom = tmp; // 1st arg to --time
                    lastopt = 't';  // maybe there is a a 2nd arg too
                    break;
                */

                case 1:
                    /* This means a field is unknown, or could be multiple arg or bad arg*/
                    /*
                    if (last_c=='t') {  // --time 2nd arg (or not if not a number)
                        errno = 0;
                        tmp = strtol(optarg, (char **)NULL, 0);
                        if (!errno) {
                            timeto = tmp;
                            printf("Time set to %d - %d\n", timefrom, timeto);
                            timearg=true;
                        }
                    } 
                    */
                    if (!timearg) {
                        fprintf(stderr, "Unrecognized argument: %s\n", optarg);
                        return 1;
                    }
                    break;

                default:
                    printf("Processing default: %c\n", c);
                    break;
        } /* end switch */
        last_c = c;
    } /* end while */

    /* Check if we have a file defined */
    if (optind >= argc) {
        fprintf(stderr,"Missing file name\n");
        return 1;
    }
    vfile = strdup(argv[optind++]);
    while (optind < argc) {
        varmask[nmasks] = strndup(argv[optind++],256);
        nmasks++;
    }
 
    /* Process dimension specifications */
    if (start != NULL) parseDimSpec(start, istart);
    if (count != NULL) parseDimSpec(count, icount);

    // process the regular expressions
    if (use_regexp) {
        retval = compile_regexp_masks();
        if (retval)
            return retval;
    }

    if (verbose) 
        printSettings();
   
    retval = print_start(outpath);
    if (retval)
        return retval;

    /* Start working */
    retval = doList(vfile);

    print_stop();

   /* Free allocated memories */
   //myfree(prgname);
   myfree(outpath);
   for (i=0; i<nmasks; i++) {
        myfree(varmask[i]);
        regfree(&(varregex[i]));
   }
   myfree(vfile);

   return retval; 
}

void init_globals(void) {
    int i;
    // variables for arguments
    outpath              = NULL;
    for (i=0; i<MAX_MASKS; i++)
        varmask[i]       = NULL;
    nmasks               = 0;
    vfile                = NULL;
    start                = NULL;
    count                = NULL;
    verbose              = 0;
    ncols                = 6;    // by default when printing ascii, print "X Y", not X: Y1 Y2...
    dump                 = false;
    output_xml           = false;
    noindex              = false;
    sortnames            = false;
    listattrs            = false;
    readattrs            = false;
    readscalars          = false;
    //timefrom             = 1;
    //timeto               = -1;
    use_regexp           = false;
    formatgiven         = false;
    for (i=0; i<MAX_DIMS; i++) {
        istart[i]  = 0;
        icount[i]  = -1;  // read full var by default
    }
    ndimsspecified = 0;
}


#define PRINT_DIMS(str, v, n, loopvar) printf("%s = { ", str); \
    for (loopvar=0; loopvar<n;loopvar++) printf("%d ", v[loopvar]);    \
    printf("}")

void printSettings(void) {
    int i;
    printf("Settings:\n");
    printf("  masks : %d ", nmasks);
    for (i=0; i<nmasks; i++)
        printf("%s ", varmask[i]);
    printf("\n");
    printf("  file  : %s\n", vfile);
    printf("  output: %s\n", outpath);

    if (start != NULL) {
        PRINT_DIMS("  start", istart, ndimsspecified,i); printf("\n");
    }
    if (count != NULL) {
        PRINT_DIMS("  count", icount, ndimsspecified,i); printf("\n");
    }

    if (output_xml)
        printf("  output data in XML format\n");
}

void bpexit(int code, int64_t fh, int64_t gh) {
    if (gh > 0)
        adios_gclose (gh);
    if (fh > 0)
        adios_fclose (fh);
    exit(code);
}

int doList(const char *path) {
    int64_t fh, gh;
    int     ntsteps;

    BP_FILE_INFO finfo;
    BP_GROUP_INFO ginfo; // reused for each group
    int     gr, i, j, n;             // loop vars
    int     status;
    int     vartype, ndims, dims[MAX_DIMS]; // info about one variable 
    int     attrsize;                       // info about one attribute
    int     hastimesteps;             // variable is spread among timesteps in the file
    int     mpi_comm_dummy;
    bool    matches;
    int     len, maxlen;
    int     nVarsMatched=0;
    int     nGroupsMatched=0;
    int     retval;
    char    commentchar;
    char  **names;  // vars and attrs together, sorted or unsorted
    bool   *isVar;  // true for each var, false for each attr
    int     nNames; // number of vars + attrs
    
    if (verbose>1) printf("\nADIOS BP open: read header info from %s\n", path);

    // open the BP file
    status = adios_fopen (&fh, path, mpi_comm_dummy); 
    if (status != 0) {
	fprintf(stderr, "readadiosbp: error opening bp file %s\n", path);
	bpexit(7, 0, 0);
    }

    // get number of groups, variables, timesteps, and attributes 
    // all parameters are integers, 
    // besides the last parameter, which is an array of strings for holding the list of group names
    adios_init_fileinfo ( &finfo, 1);
    adios_inq_file (fh, &finfo);
    ntsteps = finfo.tidx_stop - finfo.tidx_start + 1;
    if (verbose) {
         printf ("File info:\n");
         printf ("  of groups: %d\n", finfo.groups_count);  //ngroups);
         printf ("  of variables: %d\n", finfo.vars_count); // nvars);
         printf ("  of attributes: %d\n", finfo.attrs_count); // nattrs);
         printf ("  time steps: %d - %d\n", finfo.tidx_start, finfo.tidx_stop);
         printf ("\n");
    }

    if (noindex) commentchar = ';';
    else         commentchar = ' ';

    // each group has to be handled separately
    for (gr=0; gr<finfo.groups_count; gr++) {
        if (!grpMatchesMask(finfo.group_namelist[gr]))
            continue;
        nGroupsMatched++;
        if (!dump) fprintf(outf, "Group %s:\n", finfo.group_namelist[gr]);
        adios_gopen (fh, &gh, finfo.group_namelist[gr]);
        adios_init_groupinfo( &ginfo, 1);
        // get variable info from group
        adios_inq_group (gh, &ginfo);

        //printf("Attrs:\n");
        //for (i=0; i<ginfo.attrs_count; i++) {
        //   printf("%s\n", ginfo.attr_namelist[i]);
        //}

        if (sortnames) {
            qsort(ginfo.var_namelist, ginfo.vars_count, sizeof(char *), cmpstringp);
            if (listattrs)
                qsort(ginfo.attr_namelist, ginfo.attrs_count, sizeof(char *), cmpstringp);
        }

        if (listattrs)
            nNames = ginfo.vars_count + ginfo.attrs_count;
        else 
            nNames = ginfo.vars_count;
        names = (char **) malloc (nNames * sizeof (char*)); // store only pointers
        isVar = (bool *) malloc (nNames * sizeof (bool));
        if (names == NULL || isVar == NULL) {
            fprintf(stderr, "Error: could not allocate char* and bool arrays of %d elements\n", nNames);
            return 5;
        }
        mergeLists(ginfo.vars_count, ginfo.var_namelist, ginfo.attrs_count, ginfo.attr_namelist,
                  names, isVar);

        // calculate max length of variable names in the first round
        maxlen = 4;
        for (n=0; n<nNames; n++) {
            len = strlen(names[n]);
            if (len > maxlen) maxlen = len;
        }
        
        /* VARIABLES */
        for (n=0; n<nNames; n++) {
            matches = false;
            if (isVar[n])  {
                adios_inq_var (gh, names[n], &vartype, &ndims, &hastimesteps, dims);
            } else {
                adios_inq_attr (gh, names[n], &vartype, &attrsize);
                ndims = 0;
                hastimesteps = 0;
            }

            matches = matchesAMask(names[n]);

            if (matches) {
                nVarsMatched++;
                // print definition of variable
                fprintf(outf,"%c %-*s  %-*s", commentchar, 8, bp_type_to_string(vartype), maxlen, names[n]); 
                if (!isVar[n]) {
                    // list (and print) attribute
                    if (readattrs || dump) {
                        fprintf(outf,"  attr   = ");
                        retval = readAttr(gh, names[n], vartype, attrsize);
                        if (retval) return retval;
                        matches = false; // already printed
                    } else {
                        fprintf(outf,"  attr\n");
                    }
                } else if (ndims > 0 || hastimesteps) {
                    if (hastimesteps) {
                        fprintf(outf,"  {%d", ntsteps);
                        j=0;
                    } else {
                        fprintf(outf,"  {%d", dims[0]);
                        j=1;
                    }
                    for (; j < ndims; j++)
                       fprintf(outf,", %d", dims[j]);
                    fprintf(outf,"}\n");
                } else {
                    if (readscalars || dump) {
                        fprintf(outf,"  scalar = ");
                        retval = readVar(gh, names[n], vartype, ndims, dims, hastimesteps, 
                                          finfo.tidx_start, finfo.tidx_stop); 
                        if (retval) return retval;
                        matches = false; // already printed
                    } else {
                        fprintf(outf,"  scalar\n");
                    }
                }

                if (dump) {
                    // print variable content 
                    if (isVar[n])
                        retval = readVar(gh, names[n], vartype, ndims, dims, hastimesteps, 
                                          finfo.tidx_start, finfo.tidx_stop); 
                    else 
                        retval = readAttr(gh, names[n], vartype, attrsize);
                    if (retval)
                        return retval;
                    fprintf(outf,"\n");
                }
            }
        }
        adios_free_groupinfo(&ginfo);
        adios_gclose (gh);
        free(names);
        free(isVar);
    }

    if (grpmask != NULL && nGroupsMatched == 0) {
        fprintf(stderr, "\nError: None of the groups matched the group mask you provided: %s\n", grpmask);
        return 4;
    }
    if (nmasks > 0 && nVarsMatched == 0) {
        fprintf(stderr, "\nError: None of the variables matched any name/regexp you provided\n");
        return 4;
    }
    adios_free_fileinfo(&finfo);
    adios_fclose (fh);
    return 0;
}                


int cmpstringp(const void *p1, const void *p2)
{
    /* The actual arguments to this function are "pointers to
       pointers to char", but strcmp() arguments are "pointers
       to char", hence the following cast plus dereference */
    return strcmp(* (char * const *) p1, * (char * const *) p2);
}
/** Merge listV with listA if listattrs=true, otherwise just return listV.
    If sortnames=true, quicksort the result list.
*/
void mergeLists(int nV, char **listV, int nA, char **listA, char **mlist, bool *isVar) 
{
    int v, a;
    if (sortnames && listattrs) {
        // merge sort the two lists
        v = 0;
        a = 0;
        while (v < nV || a < nA) {
            if (a < nA && (v >= nV || strcmp(listV[v], listA[a]) > 0)) {
                // fully consumed var list or 
                // next item in attr list is less than next item in var list
                mlist[v+a] = listA[a];
                isVar[v+a] = false;
                a++;
            } else {
                mlist[v+a] = listV[v];
                isVar[v+a] = true;
                v++;
            }
        }
    } else {
        // first add vars then attrs (if asked)
        for (v=0; v<nV; v++) {
                mlist[v] = listV[v];
                isVar[v] = true;
        }
        if (listattrs) {
            for (a=0; a<nA; a++) {
                    mlist[a+nV] = listA[a];
                    isVar[a+nV] = false;
            }
        }
    }
}

int getTypeInfo( int adiosvartype, int* elemsize)
{
    switch(adiosvartype) {
        case adios_unsigned_byte:
            *elemsize = 1;
            break;
        case adios_byte:
            *elemsize = 1;
            break;
        case adios_string:
            *elemsize = 1;
            break;
            
        case adios_unsigned_short:  
            *elemsize = 2;
            break;
        case adios_short:
            *elemsize = 2;
            break;
            
        case adios_unsigned_integer:
            *elemsize = 4;
            break;
        case adios_integer:    
            *elemsize = 4;
            break;

        case adios_unsigned_long:
            *elemsize = 8;
            break;
        case adios_long:        
            *elemsize = 8;
            break;

        case adios_real:
            *elemsize = 4;
            break;

        case adios_double:
            *elemsize = 8;
            break;

        case adios_complex:  
            *elemsize = 8;
            break;

        case adios_double_complex:
            *elemsize = 16;
            break;

        case adios_long_double: // do not know how to print
            // *elemsize = 16;
        default:
	    return 1;
    }
    return 0;
}

/** Read data of a variable and print 
  * Return: 0: ok, != 0 on error
  */
int readVar(int64_t gh, char *name, int vartype, int ndims, int *dims, int hastimesteps, int tidx_start, int tidx_stop) 
{
    int i,j;
    int start_t[MAX_DIMS], count_t[MAX_DIMS]; // processed <0 values in start/count
    int s[MAX_DIMS], c[MAX_DIMS]; // for block reading of smaller chunks
    uint64_t nelems;         // number of elements to read
    int elemsize;            // size in bytes of one element
    int st, ct;
    int timeFrom, timeSteps; // to handle time dimension if present
    void *data;
    uint64_t sum;           // working var to sum up things
    int  maxreadn;          // max number of elements to read once up to a limit (10MB of data)
    int  actualreadn;       // our decision how much to read at once
    int  readn[MAX_DIMS];   // how big chunk to read in in each dimension?
    int64_t bytes_read;     // retval from adios_get_var()
    bool incdim;            // used in incremental reading in

    if (getTypeInfo(vartype, &elemsize)) {
        fprintf(stderr, "Adios type %d (%s) not supported in bpls. var=%s\n", 
                vartype, bp_type_to_string(vartype), name);
        return 10;
    }

    // create the counter arrays with the appropriate lengths
    // transfer start and count arrays to format dependent arrays

    // handle timesteps first
    if (hastimesteps) {
        if (istart[0] < 0)  // negative index means lasttime-|index|+1
            timeFrom = tidx_stop+istart[i]+1;
        else
            timeFrom = tidx_start+istart[0];
        if (icount[0] < 0)     // negative index means timelast-|index|+1-start
            timeSteps = tidx_stop+icount[0]+1-timeFrom+tidx_start;
        else
            timeSteps = icount[0];
        i=1;  // process istart[], icount[] from second index
    } else {
        timeFrom = tidx_start;  // default: no timesteps in file, so read in 1 piece (from 0)
        timeSteps = 1; 
        i=0;
    }
    nelems  = 1;
    // index j for start_t/count_t (lags behind i with 1 if 'hastimesteps')
    for (j=0; j<ndims; j++) {
        if (istart[i] < 0)  // negative index means last-|index|
            st = dims[j]+istart[i];
        else
            st = istart[i];
        if (icount[i] < 0)  // negative index means last-|index|+1-start
            ct = dims[j]+icount[i]+1-st;
        else
            ct = icount[i];

        if (verbose>2) 
            printf("    j=%d, st=%d ct=%d\n", j, st, ct);

        start_t[j] = st;
        count_t[j] = ct;
        nelems *= ct;
        if (verbose>1) 
            printf("    s[%d]=%d, c[%d]=%d, n=%lld\n", j, start_t[j], j, count_t[j], nelems);
        i++;
    }

    if (verbose>1) {
        printf(" total size of data to read %d times = %lld\n", timeSteps, nelems*elemsize);
    }
    
    maxreadn = MAX_BUFFERSIZE/elemsize;
    if (nelems < maxreadn)
        maxreadn = nelems;
   
    // allocate data array
    data = (void *) malloc (maxreadn*elemsize+8); // +8 for just to be sure

    // determine strategy how to read in:
    //  - at once
    //  - loop over 1st dimension
    //  - loop over 1st & 2nd dimension
    //  - etc
    if (verbose>1) printf("Read size strategy:\n");
    sum = (uint64_t) 1;
    actualreadn = (uint64_t) 1;
    for (i=ndims-1; i>=0; i--) {
        if (sum >= (uint64_t) maxreadn) {
            readn[i] = 1;
        } else {
            readn[i] = maxreadn / (int)sum; // sum is small for 4 bytes here
            // this may be over the max count for this dimension
            if (readn[i] > count_t[i]) 
                readn[i] = count_t[i];
        }
        if (verbose>1) printf("    dim %d: read %d elements\n", i, readn[i]);
        sum = sum * (uint64_t) count_t[i];
        actualreadn = actualreadn * readn[i];
    }
    if (verbose>1) printf("    read %d elements at once, %lld in total (nelems=%lld)\n", actualreadn, sum, nelems);

    // Start reading. Outer loop is by timesteps.
    for (i=0; i<timeSteps; i++) {
        // init s and c
        for (j=0; j<ndims; j++) {
            s[j]=start_t[j];
            c[j]=readn[j];
        }

        sum = 0;
        while (sum < nelems) {
            
            // how many elements do we read in next?
            actualreadn = 1;
            for (j=0; j<ndims; j++) 
                actualreadn *= c[j];

            if (verbose>2) {
                printf("adios_get_var name=%s timestep=%d", name, timeFrom+i);
                PRINT_DIMS("  start", s, ndims, j); 
                PRINT_DIMS("  count", c, ndims, j); 
                printf("  read %d elems\n", actualreadn);
            }

            // read a slice finally
            bytes_read = adios_get_var (gh, name, data, s, c, timeFrom+i); 

            if (bytes_read < 0) {
                fprintf(stderr, "Error when reading variable %s\n", name);
                free(data);
                return 11;
            }

            if (verbose>2) printf("  read %lld bytes\n", bytes_read);

            // print slice
            print_data(data, vartype, elemsize, s, c, ndims, dims, timeFrom+i, hastimesteps); 

            // prepare for next read
            sum += actualreadn;
            incdim=true; // largest dim should be increased 
            for (j=ndims-1; j>=0; j--) {
                if (incdim) {
                    if (s[j]+c[j] == start_t[j]+count_t[j]) {
                        // reached the end of this dimension
                        s[j] = start_t[j];
                        c[j] = readn[j];
                        incdim = true; // next smaller dim can increase too
                    } else {
                        // move up in this dimension up to total count
                        s[j] += readn[j];
                        if (s[j]+c[j] > start_t[j]+count_t[j]) {
                            // do not reach over the limit
                            c[j] = start_t[j]+count_t[j]-s[j];
                        }
                        incdim = false;
                    }
                }
            }
        } // end while sum < nelems
    } // end for timesteps
    print_endline();
    

    free(data);
    return 0;
}

int readAttr( int64_t gh, char *name, int attrtype, int size) 
{
    void *data;
    int elemsize, readsize;

    // size already has the size, getTypeInfo is used for type checking only
    if (getTypeInfo(attrtype, &elemsize)) {
        fprintf(stderr, "Adios type %d (%s) not supported in bpls. attr=%s\n", 
                attrtype, bp_type_to_string(attrtype), name);
        return 10;
    }

    data = (void *) malloc (size);
    if (data == NULL) {
        fprintf(stderr, "Error: could not allocate %d bytes to read an attribute\n");
        return 11;
    }

    readsize = adios_get_attr(gh, name, data);
    if (readsize != size) {
        fprintf(stderr, "Warning: attribute reading reported different size (%d) "
                        "than we allocated (%d) based on the inqury results\n", readsize, size);
    }

    print_data(data, attrtype, size, NULL, NULL, 0, NULL, 0, false); 
    print_endline();

    free(data);
    return 0;
}

/*
bool conformsToDimSpec(VarInfo vi) {
    if (timeDimidx > vi.ndims) {
            fprintf(stderr,"The time dimension (%ld) is out of the dimensions of this variable %s (%d). Skip it.\n", 
                    timeDimidx, vi.name, vi.ndims);
            return false;
    }
    if (timeDimidx == 0 && vi.ndims == 1) { 
            fprintf(stderr,"Variable %s has only 1 dimension and time is specified. Skip it.\n", vi.name);
            return false;
    }
    return true;
}
*/

bool matchesAMask(char *name) {
    int excode;
    int i;
    int startpos=0; // to match with starting / or without
    regmatch_t pmatch[1] = {{ (regoff_t) -1, (regoff_t) -1}};

    if (nmasks == 0) 
        return true;

    for (i=0; i<nmasks; i++) {
        if (use_regexp) {
            excode = regexec (&(varregex[i]), name, 1, pmatch, 0);
            if (name[0] == '/') // have to check if it matches from the second character too
                startpos = 1;
            if (excode == 0 &&                  // matches
                (pmatch[0].rm_so == 0 || pmatch[0].rm_so == startpos) &&         // from the beginning
                pmatch[0].rm_eo == strlen(name) // to the very end of the name
               ) {
                if (verbose>1)
                    printf("Name %s matches regexp %i %s\n", name, i, varmask[i]);
                //printf("Match from %d to %d\n", (int) pmatch[0].rm_so, (int) pmatch[0].rm_eo);
                return true;
            }
        } else {
            // use shell pattern matching
            if (varmask[i][0] != '/' && name[0] == '/')
                startpos = 1;
            if ( fnmatch( varmask[i], name+startpos, FNM_FILE_NAME) == 0) {
                if (verbose>1)
                    printf("Name %s matches varmask %i %s\n", name, i, varmask[i]);
                return true; 
            }
        }
    }
    return false;
}

/* return true if mask is null */
bool grpMatchesMask(char *name) {
    int startpos=0;
    int excode;
    regmatch_t pmatch[1] = {{ (regoff_t) -1, (regoff_t) -1}};

    if (grpmask == NULL)
        return true;

    if (use_regexp) {
        excode = regexec (&(grpregex), name, 1, pmatch, 0);
        if (name[0] == '/') // have to check if it matches from the second character too
            startpos = 1;
        if (excode == 0 &&                  // matches
            (pmatch[0].rm_so == 0 || pmatch[0].rm_so == startpos) &&         // from the beginning
            pmatch[0].rm_eo == strlen(name) // to the very end of the name
            ) {
            if (verbose>1)
                printf("Name %s matches regexp %s\n", name, grpmask);
            //printf("Match from %d to %d\n", (int) pmatch[0].rm_so, (int) pmatch[0].rm_eo);
            return true;
        }
    } else {
        // use shell pattern matching
        if (grpmask[0] != '/' && name[0] == '/')
            startpos = 1;
        if ( fnmatch( grpmask, name+startpos, FNM_FILE_NAME) == 0) {
            if (verbose>1)
                printf("Name %s matches groupmask %s\n", name, grpmask);
            return true; 
        }
    }
    return false;
}


int  print_start(const char *fname) {
    if ( fname == NULL) {
        outf = stdout;
    } else {
        if ((outf = fopen(fname,"w")) == NULL) {
            fprintf(stderr, "Error at opening for writing file %s: %s\n",
                    fname, strerror(errno));
            return 30;
        }
    }
    return 0;
}

void print_stop() {
    fclose(outf);
}

static int nextcol=0;  // column index to start with (can have lines split in two calls)

int print_data(void *data, int adiosvartype, int elemsize,
                int *s, int* c, int ndims, int *dims, int time, bool hastimesteps)
{
    int i,item, steps;
    char idxstr[128], vstr[128], buf[16];
    int ids[MAX_DIMS], idsdims; // current indices
    bool roll;

    // init current indices
    if (hastimesteps) {
        ids[0] = time-1; // time in adios starts with 1, we show from 0 as dimension
        idsdims = ndims+1;
    } else {
        idsdims = ndims;
    }
    steps = 1;
    for (i=1; i<=ndims; i++) {
        ids[idsdims-i] = s[ndims-i];
        steps *= c[ndims-i];
    }
    
    item = 0; // index to *data 
    // loop through each data item and print value
    while (item < steps) {

        // print indices if needed into idxstr;
        idxstr[0] = '\0'; // empty idx string
        if (nextcol == 0) {
            if (!noindex && idsdims > 0) {
                sprintf(idxstr,"  (%d",ids[0]);
                for (i=1; i<idsdims; i++) {
                    sprintf(buf,",%d",ids[i]);
                    strcat(idxstr, buf);
                }
                strcat(idxstr,")    ");
            }
        }

        // print next data item into vstr
        switch(adiosvartype) {
            case adios_unsigned_byte:
                snprintf(vstr,127,(formatgiven ? format : "%hhu "), ((unsigned char *) data)[item]);
                break;
            case adios_byte:
                snprintf(vstr,127,(formatgiven ? format : "%hhd "), ((char *) data)[item]);
                break;
            case adios_string:
                snprintf(vstr,127,(formatgiven ? format : "\"%s\""), ((char *) data)+item);
                break;
        
            case adios_unsigned_short:  
                snprintf(vstr,127,(formatgiven ? format : "%hu "), ((unsigned short *) data)[item]);
                break;
            case adios_short:
                snprintf(vstr,127,(formatgiven ? format : "%hd "), ((short *) data)[item]);
                break;
        
            case adios_unsigned_integer:
                snprintf(vstr,127,(formatgiven ? format : "%u "), ((unsigned int *) data)[item]);
                break;
            case adios_integer:    
                snprintf(vstr,127,(formatgiven ? format : "%d "), ((int *) data)[item]);
                break;

            case adios_unsigned_long:
                snprintf(vstr,127,(formatgiven ? format : "%llu "), ((unsigned long long *) data)[item]);
                break;
            case adios_long:        
                snprintf(vstr,127,(formatgiven ? format : "%lld "), ((long long *) data)[item]);
                break;

            case adios_real:
                snprintf(vstr,127,(formatgiven ? format : "%g "), ((float *) data)[item]);
                break;

            case adios_double:
                snprintf(vstr,127,(formatgiven ? format : "%g "), ((double *) data)[item]);
                break;

            /*
            case adios_long_double:
                snprintf(vstr,127,(formatgiven ? format : "%g "), ((double *) data)[item]);
                break;
            */

            case adios_complex:  
                snprintf(vstr,127,(formatgiven ? format : "(%g,i%g) "), ((float *) data)[2*item], ((float *) data)[2*item+1]);
                break;

            case adios_double_complex:
                snprintf(vstr,127,(formatgiven ? format : "(%g,i%g)" ), ((double *) data)[2*item], ((double *) data)[2*item+1]);
                break;
        } // end switch

        // print item
        fprintf(outf, "%s%s", idxstr, vstr);

        // increment/reset column index
        nextcol++;
        if (nextcol == ncols) {
            fprintf(outf,"\n");
            nextcol = 0;
        }
            
        // increment indices
        item++;
        roll = true;
        for (i=1; i<=ndims; i++) {
            if (roll) {
                if (ids[idsdims-i] == s[ndims-i]+c[ndims-i]-1 ) {
                    // last index in this dimension, roll upward
                    ids[idsdims-i] = s[ndims-i];
                } else {
                    ids[idsdims-i]++;
                    roll = false;
                }
            }
        }
    }
}

void print_endline(void) 
{
    if (nextcol != 0)
        fprintf(outf,"\n");
    nextcol = 0;
}



// parse a string "0, 3; 027" into an integer array
// of [0,3,27] 
// exits if parsing failes
void parseDimSpec(char *str, int *dims)
{
    char *token, *saveptr;
    char *s;  // copy of s; strtok modifies the string
    int  i=0;

    s = strndup(str, 256);
    token = strtok_r(s, " ,;x\t\n", &saveptr);
    while (token != NULL && i < MAX_DIMS) {
        //printf("\t|%s|", token);
        errno = 0;
        dims[i] = strtol(token, (char **)NULL, 0);
        if (errno) {
            fprintf(stderr, "Error: could not convert field into a value: %s from \"%s\"\n", token, str);
            exit(200);
        }

        // get next item
        token = strtok_r(NULL, " ,;x\t\n", &saveptr);
        i++;
    }
    //if (i>0) printf("\n");

    if (i > ndimsspecified) ndimsspecified = i;

    // check if number of dims specified is larger than we can handle
    if (token != NULL) {
        fprintf(stderr, "Error: More dimensions specified in \"%s\" than we can handle (%d)\n", str, MAX_DIMS);
        exit(200);
    }
}

int compile_regexp_masks(void)
{
    int i, errcode;
    char buf[256];
    for (i=0; i<nmasks; i++) {
        errcode = regcomp( &(varregex[i]), varmask[i], REG_EXTENDED);
        if (errcode) {
            regerror(errcode, &(varregex[i]), buf, sizeof(buf));
            fprintf(stderr, "Error: var %s is an invalid extended regular expression: %s\n", varmask[i], buf);
            return 2;
        }
    }
    if (grpmask != NULL) {
        errcode = regcomp( &(grpregex), grpmask, REG_EXTENDED);
        if (errcode) {
            regerror(errcode, &(grpregex), buf, sizeof(buf));
            fprintf(stderr, "Error: var %s is an invalid extended regular expression: %s\n", grpmask, buf);
            return 2;
        }
    }
    return 0;
}
