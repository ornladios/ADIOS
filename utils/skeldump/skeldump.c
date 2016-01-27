/*
 * ADIOS is freely available under the terms of the BSD license described
 * in the COPYING file in the top level directory of this source distribution.
 *
 * Copyright (c) 2008 - 2009.  UT-BATTELLE, LLC. All rights reserved.
 */

/* 
 * Extract metadata needed by skel to recreate I/O pattern of a BP file.
 *   This is a prototype version that contains much unnecessary stuff that
 *   is leftover from bpls but not needed here. It should be cleaned up 
 *   eventually.
 *
 * Adapted from bpls.c
 *
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


#include "skeldump.h"
#include "bp_utils.h"
#include "adios_read.h"
#include "adios_types.h"
 

#ifdef DMALLOC
#include "dmalloc.h"
#endif

// global variables 
// Values from the arguments or defaults
char *outpath;              // output files' starting path (can be extended with subdirs, names, indexes)
char *varmask[MAX_MASKS];   // can have many -var masks (either shell patterns or extended regular expressions)
char *grpmask;              // list only groups matching the mask
int  nmasks;                // number of masks specified
char *vfile;                // file name to bpls
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
bool attrsonly;            // do list attributes only
bool readattrs;            // also read all attributes and print
bool longopt;              // -l is turned on
bool timestep;
bool noindex;              // do no print array indices with data
bool printByteAsChar;      // print 8 bit integer arrays as string
bool plot;                 // dump histogram related information
bool hidden_attrs;         // show hidden attrs in BP file
bool show_decomp;          // show decomposition of arrays

// other global variables
char *prgname; /* argv[0] */
//long timefrom, timeto;
int  istart[MAX_DIMS], icount[MAX_DIMS], ndimsspecified=0;
regex_t varregex[MAX_MASKS]; // compiled regular expressions of varmask
regex_t grpregex;            // compiled regular expressions of grpmask
int  ncols = 6; // how many values to print in one row (only for -p)
int  verbose = 0;
FILE *outf;   // file to print to or stdout
char commentchar;

struct option options[] = {
    {"help",                 no_argument,          NULL,    'h'},
    {"verbose",              no_argument,          NULL,    'v'},
    {"dump",                 no_argument,          NULL,    'd'},
    {"group",                no_argument,          NULL,    'g'},
//    {"regexp",               no_argument,          NULL,    'e'},
//    {"plot",                 no_argument,          NULL,    'p'},
    {"output",               required_argument,    NULL,    'o'},
//    {"xml",                  no_argument,          NULL,    'x'},
//    {"start",                required_argument,    NULL,    's'}, 
//    {"count",                required_argument,    NULL,    'c'}, 
//    {"noindex",              no_argument,          NULL,    'y'},
//    {"sort",                 no_argument,          NULL,    'r'},
//    {"timestep",             no_argument,          NULL,    't'},
//    {"attrs",                no_argument,          NULL,    'a'},
//    {"attrsonly",            no_argument,          NULL,    'A'},
//   {"long",                 no_argument,          NULL,    'l'},
//    {"string",               no_argument,          NULL,    'S'},
//    {"columns",              required_argument,    NULL,    'n'}, 
//    {"format",               required_argument,    NULL,    'f'}, 
    {"hidden_attrs",         no_argument,          (int*)&hidden_attrs,    true},
//    {"decomp",               no_argument,          NULL,    'D'},
    //    {"time",                 required_argument,    NULL,    't'}, 
    {NULL,                   0,                    NULL,    0}
};


static const char *optstring = "hvepyrtaAldSDg:o:x:s:c:n:f:";

// help function
void display_help() {
    //printf( "Usage: %s  \n", prgname);
    printf("usage: skeldump [OPTIONS] file\n"
//    printf("usage: skeldump [OPTIONS] file [mask1 mask2 ...]\n"
            "\nList/dump content of a BP file. \n"
//            "A mask can be a shell pattern like with 'ls' e.g. \"*/x?\".\n"
/*
            "Variables with multiple timesteps are reported with an extra dimensions.\n"
            "The time dimension is the first dimension then.\n"
            "\n"
            "  --long      | -l           Print values of all scalars and attributes and\n"
            "                               min/max values of arrays (no overhead to get them!)\n"
            "  --attrs     | -a           List/match attributes too\n"
            "  --attrsonly | -A           List attributes only\n"
            "  --sort      | -r           Sort names before listing\n"
            "  --timestep  | -t           Print values of timestep elements\n"
            "  --group     | -g <mask>    List/dump groups matching the mask only\n"
            "  --dump      | -d           Dump matched variables/attributes\n"
            "                               To match attributes too, add option -a\n"
            "  --regexp    | -e           Treat masks as extended regular expressions\n"
            "  --plot      | -p           Dumps the histogram information that can be read by gnuplot\n"
            "  --output    | -o <path>    Print to a file instead of stdout\n"
               "  --xml    | -x            # print as xml instead of ascii text\n"
            "  --start     | -s \"spec\"    Offset indices in each dimension \n"
            "                               (default is 0 for all dimensions) \n"
            "                               <0 is handled as in python (-1 is last)\n"
            "  --count     | -c \"spec\"    Number of elements in each dimension\n"
            "                               -1 denotes 'until end' of dimension\n"
            "                               (default is -1 for all dimensions)\n"
            "  --noindex   | -y           Print data without array indices\n"
            "  --string    | -S           Print 8bit integer arrays as strings\n"
            "  --columns   | -n \"cols\"    Number of data elements per row to print\n"
            "  --format    | -f \"str\"     Format string to use for one data item in print\n"
            "                               instead of the default. E.g. \"%%6.3f\"\n"
            "  --hidden_attrs             Show hidden ADIOS attributes in the file\n"
            "  --decomp    | -D           Show decomposition of variables as layed out in file\n"
               "  --time    | -t N [M]      # print data for timesteps N..M only (or only N)\n"
               "                              default is to print all available timesteps\n"
            "\n"
            "  Examples for slicing:\n"
            "  -s \"0,0,0\"   -c \"1,99,1\":  Print 100 elements (of the 2nd dimension).\n"
            "  -s \"0,0\"     -c \"1,-1\":    Print the whole 2nd dimension however large it is.\n"
            "  -s \"-1,-1\"   -c \"1,1\":     Print the very last element (of a 2D array)\n"
*/
            "\n"
            "Help options\n"
            "  --help      | -h           Print this help.\n"
            "  --verbose   | -v           Print log about what this program is doing.\n"
            "                               Use multiple -v to increase logging level.\n"
            "Typical use: skeldump <file> > <outfile>\n"
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
    int c;
    /* Process the arguments */
    while ((c = getopt_long(argc, argv, optstring, options, NULL)) != -1) {
        switch (c) {
            case 'a':
                listattrs=true;
                break;
            case 'A':
                listattrs=true;
                attrsonly=true;
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
            case 'r':
                //sortnames = true;
                break;
            case 'l':
                longopt = true;
                readattrs = true;
                break;
            case 'n':
                errno = 0;
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
            case 'p':
                plot = true; 
                break;
            case 's':
                start = strndup(optarg,256);
                break;
            case 'S':
                printByteAsChar = true;
                break;
            case 't':
                timestep = true;
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
            case 'D':
                show_decomp = true;
                break;
                /*
                   case 't':
                   errno = 0;
                   tmp = strtol(optarg, (char **)NULL, 0);
                   if (errno) {
                   fprintf(stderr, "Error: could not convert --time value: %s\n", optarg);
                   return 1;
                   }
                   timefrom = tmp; // 1st arg to --time
                   lastopt = 't';  // maybe there is a a 2nd arg too
                   break;
                 */

            //case 200:
            //    hidden_attrs = true;
            //    break;

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

    if (noindex) commentchar = ';';
    else         commentchar = ' ';


    if (verbose>1) 
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
    timestep             = false;
    sortnames            = false;
    listattrs            = false;
    attrsonly            = false;
    readattrs            = false;
    longopt              = false;
    //timefrom             = 1;
    //timeto               = -1;
    use_regexp           = false;
    plot                 = false;
    hidden_attrs         = false;
    formatgiven          = false;
    printByteAsChar      = false;
    show_decomp          = false;
    for (i=0; i<MAX_DIMS; i++) {
        istart[i]  = 0;
        icount[i]  = -1;  // read full var by default
    }
    ndimsspecified = 0;
}


#define PRINT_DIMS(str, v, n, loopvar) printf("%s = { ", str); \
    for (loopvar=0; loopvar<n;loopvar++) printf("%d ", v[loopvar]);    \
printf("}")

#define PRINT_DIMS64(str, v64, n, loopvar) printf("%s = { ", str); \
    for (loopvar=0; loopvar<n;loopvar++) printf("%" PRId64, v64[loopvar]);    \
printf("}")

void printSettings(void) {
    int i;
    printf("Settings :\n");
    printf("  masks  : %d ", nmasks);
    for (i=0; i<nmasks; i++)
        printf("%s ", varmask[i]);
    printf("\n");
    printf("  file   : %s\n", vfile);
    if (grpmask)
        printf("  groups : %s\n", grpmask);
    printf("  output : %s\n", (outpath ? outpath : "stdout"));

    if (start != NULL) {
        PRINT_DIMS("  start", istart, ndimsspecified,i); printf("\n");
    }
    if (count != NULL) {
        PRINT_DIMS("  count", icount, ndimsspecified,i); printf("\n");
    }

    if (longopt)
        printf("      -l : show scalar values and min/max/avg of arrays\n");
    if (sortnames)
        printf("      -r : sort names before listing\n");
    if (attrsonly)
        printf("      -A : list attributes only\n");
    else if (listattrs)
        printf("      -a : list attributes too\n");
    if (dump)
        printf("      -d : dump matching variables and attributes\n");
    if (use_regexp)
        printf("      -e : handle masks as regular expressions\n");
    if (formatgiven)
        printf("      -f : dump using printf format \"%s\"\n", format);
    if (output_xml)
        printf("      -x : output data in XML format\n");
    if (show_decomp)
        printf("      -D : show decomposition of variables in the file\n");
    if (hidden_attrs)
        printf("         : show hidden attributes in the file\n");
}

    void bpexit(int code, ADIOS_FILE *fp) {
        if (fp > 0)
            adios_read_close (fp);
        exit(code);
    }

void print_file_size(uint64_t size)
{
    static const int  sn=7;
    static const char *sm[]={"bytes", "KB", "MB", "GB", "TB", "PB", "EB"};
    uint64_t s = size, r;
    int idx = 0;
    while ( s/1024 > 0 && idx < sn) {
        r = s%1024; 
        s = s/1024;
        idx++;
    }
    if (r > 511)
        s++;
    printf ("  file size:     %" PRId64 " %s\n", s, sm[idx]);

}


static inline int ndigits (int n) 
{
    static char digitstr[32];
    return snprintf (digitstr, 32, "%d", n);
}


int     nVarsMatched=0;

int doList_group (ADIOS_FILE *fp)
{
    ADIOS_VARINFO *vi; 
    ADIOS_VARINFO **vis; 
    enum ADIOS_DATATYPES vartype;
    int     i, j, n;             // loop vars
    int     attrsize;                       // info about one attribute
    bool    matches;
    int     len, maxlen, maxtypelen;
    int     retval;
    char  **names;  // vars and attrs together, sorted or unsorted
    bool   *isVar;  // true for each var, false for each attr
    int     nNames; // number of vars + attrs
    void   *value;  // scalar value is returned by get_attr
    bool    timed;  // variable has multiple timesteps


    if (attrsonly)
        nNames = fp->nattrs;
    else if (listattrs)
        nNames = fp->nvars + fp->nattrs;
    else 
        nNames = fp->nvars;

    /*
       if (sortnames) {
       if (!attrsonly)
       qsort(gp->var_namelist, gp->vars_count, sizeof(char *), cmpstringp);
       if (listattrs)
       qsort(gp->attr_namelist, gp->attrs_count, sizeof(char *), cmpstringp);
       }
     */

    names = (char **) malloc (nNames * sizeof (char*)); // store only pointers
    isVar = (bool *) malloc (nNames * sizeof (bool));
    vis   = (ADIOS_VARINFO **) malloc (nNames * sizeof (ADIOS_VARINFO*));
    if (names == NULL || isVar == NULL || vis == NULL) {
        fprintf(stderr, "Error: could not allocate char* and bool arrays of %d elements\n", nNames);
        return 5;
    }
    mergeLists(fp->nvars, fp->var_namelist, fp->nattrs, fp->attr_namelist, names, isVar);

    // calculate max length of variable names in the first round
    maxlen = 4;
    for (n=0; n<nNames; n++) {
        len = strlen(names[n]);
        if (len > maxlen) maxlen = len;
    }

    // Get VARINFO's and attr types and calculate max length of type names 
    maxtypelen = 7;
    for (n=0; n<nNames; n++) {
        if (isVar[n])  {
            vis[n] = adios_inq_var (fp, names[n]);
            if (!vis[n]) {
                fprintf(stderr, "Error: %s\n", adios_errmsg());
            }
            vartype = vis[n]->type;
        } else {
            retval = adios_get_attr (fp, names[n], &vartype, &attrsize, &value);
            if (retval) {
                fprintf(stderr, "Error: %s\n", adios_errmsg());
            }
        }
        len = strlen(adios_type_to_string(vartype));
        if (len > maxtypelen) maxtypelen = len;
    }

    /* LANGUAGE */
    // Let's peek at the first pg to see what language was used to create this file.
    // Assume the file is homogenous.
#if 1    
    BP_PROC * p = (BP_PROC *) fp->fh;
    BP_FILE * fh = (BP_FILE *) p->fh;
    fprintf (outf, "lang: %s\n", 
             fh->pgs_root->adios_host_language_fortran==adios_flag_yes?"Fortran":"C");
#endif

    fprintf (outf, "procs: %" PRIu64 "\n", fh->mfooter.pgs_count);
    fprintf (outf, "group: group1\n");


    /* VARIABLES */
    fprintf (outf,"variables: [\n");
    for (n=0; n<nNames; n++) {
        matches = false;
        if (isVar[n])  {
            vi = vis[n];
            vartype = vi->type;
            //timed = adios_read_bp_is_var_timed(fp, vi->varid);
            timed = (vi->nsteps > 1);
        } else {
            retval = adios_get_attr (fp, names[n], &vartype, &attrsize, &value);
            if (retval) {
                fprintf(stderr, "Error: %s\n", adios_errmsg());
            }
        }

        matches = matchesAMask(names[n]);

        if (matches) {
            nVarsMatched++;
            fprintf (outf, "  {\n");

            // print definition of variable
            fprintf (outf, "    name: \"%s\",\n", names[n] );
            fprintf (outf, "    type: %s,\n", adios_type_to_string(vartype) );
            if (vartype == adios_string)
            {
                fprintf (outf, "    len: %i,\n", 1024);
            }
            // fprintf(outf,"%c %-*s  %-*s", commentchar, maxtypelen, 
            //        adios_type_to_string(vartype), maxlen, names[n]); 
            if (!isVar[n]) {
                // list (and print) attribute
                if (readattrs || dump) {
                    fprintf(outf,"  attr   = ");
                    print_data(value, 0, vartype, false); 
                    fprintf(outf,"\n");
                    matches = false; // already printed
                } else {
                    fprintf(outf,"  attr\n");
                }
            } else if (!vi) { 
                // after error
                fprintf(outf, "\n");
            } else if (vi->ndim > 0 || timed) {
                // array
                //fprintf(outf,"  {%s%d", (vi->timedim==0 ? "T-": ""),vi->dims[0]);

                //fprintf(outf,"    ");
                if (timed) 
                    fprintf(outf, "%d*", vi->nsteps);


                adios_inq_var_blockinfo (fp, vi);
                fprintf(outf, "    dims: ");
                if (vi->ndim > 0) {
                    if (vi->blockinfo)
                        fprintf(outf,"[%" PRId64, vi->blockinfo[0].count[0]);
                    else
                        fprintf(outf,"[%" PRId64, vi->dims[0]);
                    for (j=1; j < vi->ndim; j++) {
                        if (vi->blockinfo)
                            fprintf(outf,", %" PRId64, vi->blockinfo[0].count[j]);
                        else
                            fprintf(outf,", %" PRId64, vi->dims[j]);
                    }
                    fprintf(outf,"],");
                } else {
                    fprintf(outf,"scalar");
                }

                if (longopt || plot) {
                    adios_inq_var_stat (fp, vi, timestep && timed, 0);
                }

                if (plot && vi->statistics && vi->statistics->histogram) {
                    print_data_hist(vi, &names[n][1]);
                }


                if (longopt && vi->statistics) {

                    if(timestep == false || timed == false ) {

                        fprintf(outf," = ");
                        print_data(vi->statistics->min, 0, vartype, false); 

                        fprintf(outf,"/ ");
                        print_data(vi->statistics->max, 0, vartype, false); 

                        if(vartype == adios_complex || vartype == adios_double_complex) {

                            fprintf(outf,"/ ");
                            print_data(vi->statistics->avg, 0, adios_double_complex, false);
                        } else {

                            fprintf(outf,"/ ");
                            print_data(vi->statistics->avg, 0, adios_double, false);
                        }

                        if(vartype == adios_complex || vartype == adios_double_complex) {

                            fprintf(outf,"/ ");
                            print_data(vi->statistics->std_dev, 0, adios_double_complex, false);
                        } else {

                            fprintf(outf,"/ ");
                            print_data(vi->statistics->std_dev, 0, adios_double, false);
                        }

                        //fprintf(outf," {MIN / MAX / AVG / STD_DEV} ");
                    } else {
                        int time_start = 0, time_end = vi->nsteps;

                        if (start != NULL) {
                            if (istart[0] >= 0)
                                time_start = istart[0];
                            else
                                time_start = vi->nsteps - 1 + istart[0];
                        }

                        if (count != NULL) {
                            if(icount[0] > 0)
                                time_end = time_start + icount[0];
                            else
                                time_end = vi->nsteps + icount[0] + 1;
                        }

                        if (time_start < 0 || time_start >= vi->nsteps) {
                            fprintf (stderr, "Error when reading variable %s. errno=%d : Variable (id=%d) has no data at %d time step\n", names[n], 15, vi->varid, time_start);
                            bpexit(15,fp);
                        }

                        if (time_end < 0 || time_end > vi->nsteps) {
                            fprintf (stderr, "Error when reading variable %s. errno=%d : Variable (id=%d) has no data at %d time step\n", names[n], 15, vi->varid, time_end);
                            bpexit(16,fp);
                        }

                        static char *indent_char = " ";
                        int indent_len=11;

                        /* Start - Print the headers of statistics first */
                        fprintf(outf, "\n%-*s", indent_len+7, indent_char);
                        fprintf(outf, "%10s  ", "MIN");
                        fprintf(outf, "%10s  ", "MAX");
                        fprintf(outf, "%10s  ", "AVG");
                        fprintf(outf, "%10s  ", "STD DEV");

                        /* End - Print the headers of statistics first */

                        void *min, *max, *avg, *std_dev;
                        enum ADIOS_DATATYPES vt = vartype;
                        struct ADIOS_STAT_STEP *s = vi->statistics->steps;
                        if (vi->type == adios_complex || vi->type == adios_double_complex)
                            vt = adios_double;
                        fprintf(outf, "\n%-*sglobal:", indent_len, indent_char);
                        print_data_characteristics (vi->statistics->min, 
                                vi->statistics->max, 
                                vi->statistics->avg, 
                                vi->statistics->std_dev, 
                                vt, false);

                        for(i = time_start; i < time_end; i++) {
                            min = max = avg = std_dev = 0;
                            if (s->maxs && s->maxs[i]) max = s->maxs[i];
                            if (s->mins && s->mins[i]) min = s->mins[i];
                            if (s->avgs && s->avgs[i]) avg = s->avgs[i];
                            if (s->std_devs && s->std_devs[i]) std_dev = s->std_devs[i];

                            // Align the output, previous lines has atleast (maxlen + strlen(names[n])) characters
                            // Better way to printf N spaces?
                            fprintf(outf, "\n%-*st%-5d:", indent_len, indent_char, i);
                            print_data_characteristics(min, max, avg, std_dev, vt, false);
                        }
                        fprintf(outf, "\n");
                    }
                } // longopt && vi->statistics 
                fprintf(outf,"\n");

                if (/*show_decomp*/1) { // Use the bpls -D mechanism to get domain decomposition info
                    adios_inq_var_blockinfo (fp, vi);
                    print_decomp(vi);
                }

            } else {
                // scalar
                fprintf(outf,"    dims: scalar");
                if (vi->value) {
                    fprintf (outf, ",\n    value: \"");
                    print_data (vi->value, 0, vartype, false);
                    fprintf (outf, "\"");
                }
                if (longopt && vi->value) {
                    fprintf(outf," = ");
                    print_data(vi->value, 0, vartype, false); 
                    matches = false; // already printed
                }
                fprintf(outf,"\n");

                if (show_decomp) {
                    adios_inq_var_blockinfo (fp, vi);
                    print_decomp(vi);
                }
            }
            fprintf (outf, "  },\n");
        }

        if (matches && dump) {
            // print variable content 
            if (isVar[n])
                retval = readVar(fp, vi, names[n], timed);
            if (retval && retval != 10) // do not return after unsupported type
                return retval;
            fprintf(outf,"\n");
        }

        //if (isVar[n])
        //    adios_free_varinfo(vi);
        //else
        if (!isVar[n])
            free(value);

    }
    /* Free ADIOS_VARINFOs */
    fprintf (outf, "]\n");
    for (n=0; n<nNames; n++) {
        if (isVar[n])  {
            adios_free_varinfo(vis[n]);
        }
    }
    free(names);
    free(isVar);
    return 0;
}                


int doList(const char *path) {
    ADIOS_FILE  *fp;
    int     grpid;     // loop vars
    int     status;
    int     mpi_comm_dummy=0;
    int     nGroupsMatched=0;
    int     nGroups; // number of groups
    char  **group_namelist;
    char    init_params[128];

    if (verbose>1) printf("\nADIOS BP open: read header info from %s\n", path);

    // initialize BP reader
    strcpy (init_params, "verbose=2");
    if (hidden_attrs)
        strcat (init_params, ";show_hidden_attrs");
    status = adios_read_init_method (ADIOS_READ_METHOD_BP, mpi_comm_dummy, init_params);
    if (status) {
        fprintf(stderr, "Error: %s\n", adios_errmsg());
        bpexit(6, 0);
    }

    // open the BP file
    fp = adios_read_open_file (path, ADIOS_READ_METHOD_BP, mpi_comm_dummy); 
    if (fp == NULL) {
        fprintf(stderr, "Error: %s\n", adios_errmsg());
        bpexit(7, 0);
    }

    // get number of groups
    nGroups = adios_get_grouplist (fp, &group_namelist);

    //, variables, timesteps, and attributes 
    // all parameters are integers, 
    // besides the last parameter, which is an array of strings for holding the list of group names
    //ntsteps = fp->tidx_stop - fp->tidx_start + 1;
    if (verbose) {
        printf ("File info:\n");
        printf ("  of groups:     %d\n", nGroups);
        printf ("  of variables:  %d\n", fp->nvars);
        printf ("  of attributes: %d\n", fp->nattrs);
        printf ("  time steps:    %d - %d\n", fp->current_step, fp->last_step);
        print_file_size(fp->file_size);
        printf ("  bp version:    %d\n", fp->version);
        printf ("  endianness:    %s\n", (fp->endianness ? "Big Endian" : "Little Endian"));
        if (longopt) 
            printf ("  statistics:    Min / Max / Avg / Std_dev\n");
        printf ("\n");
    }


    if (grpmask) {
        // each group has to be handled separately
        for (grpid=0; grpid < nGroups; grpid++) {
            if (!grpMatchesMask(group_namelist[grpid]))
                continue;
            nGroupsMatched++;
            if (!dump) fprintf(outf, "Group %s:\n", group_namelist[grpid]);
            status = adios_group_view (fp, grpid);
            if (status) {
                fprintf(stderr, "Error: %s\n", adios_errmsg());
                bpexit(8, fp);
            }

            doList_group (fp);

            adios_group_view (fp, -1); // reset full view (for next group view)
        }
    } else {
        doList_group (fp);
    }

    if (grpmask != NULL && nGroupsMatched == 0) {
        fprintf(stderr, "\nError: None of the groups matched the group mask you provided: %s\n", grpmask);
        return 4;
    }
    if (nmasks > 0 && nVarsMatched == 0) {
        fprintf(stderr, "\nError: None of the variables matched any name/regexp you provided\n");
        return 4;
    }
    adios_read_close (fp);
    return 0;
}



int print_data_hist(ADIOS_VARINFO * vi, char * varname)
{
    char hist_file[256], gnuplot_file[256];
    int i;
    char xtics[512], str[512];
    FILE * out_hist, * out_plot;
    struct ADIOS_HIST *h = vi->statistics->histogram;

    memcpy(hist_file,  varname, strlen(varname) + 1);    
    strcat(hist_file, ".hist");

    if ((out_hist = fopen(hist_file,"w")) == NULL) {
        fprintf(stderr, "Error at opening for writing file %s: %s\n",
                hist_file, strerror(errno));
        return 30;
    }

    memcpy(gnuplot_file,  varname, strlen(varname) + 1);    
    strcat(gnuplot_file, ".gpl");

    if ((out_plot = fopen(gnuplot_file,"w")) == NULL) {
        fprintf(stderr, "Error at opening for writing file %s: %s\n",
                gnuplot_file, strerror(errno));
        return 30;
    }

    xtics[0] = '\0';
    strcat(xtics, "set xtics offset start axis (");
    for (i = 0; i <= h->num_breaks; i++)
    {
        if (i == 0)
        {
            fprintf(out_hist, "-Inf %.2lf %u\n", h->breaks[i], h->gfrequencies[i]);
            sprintf(str, "\"-Inf\" pos(%d)", i); 
        }
        else if (i < h->num_breaks)
        {
            fprintf(out_hist, "%.2lf %.2lf %u\n", h->breaks[i - 1], h->breaks[i], h->gfrequencies[i]);
            sprintf(str, ", \"%.2lf\" pos(%d)", h->breaks[i - 1], i); 
        }
        else     
        {
            fprintf(out_hist, "%.2lf Inf %u\n", h->breaks[i], h->gfrequencies[i]);
            sprintf(str, ", \"Inf\" pos(%d)", i); 
        }
        strcat(xtics, str);
    }
    strcat(xtics, ")\n");

    fprintf(out_plot, "start = -0.5\npos(x) = start + x * 1\nset boxwidth 1\nset style fill solid border 5#5lt6#6\n");
    fprintf(out_plot, "%s", xtics);
    fprintf(out_plot, "plot '%s' using 3 smooth frequency w boxes\n", hist_file);
    fprintf(out_plot, "pause -1 'Press Enter to quit'\n");
    return 0;
}

int cmpstringp(const void *p1, const void *p2)
{
    /* The actual arguments to this function are "pointers to
       pointers to char", but strcmp() arguments are "pointers
       to char", hence the following cast plus dereference */
    return strcmp(* (char * const *) p1, * (char * const *) p2);
}
/** Merge listV with listA if listattrs=true, return listA if attrsonly=true, otherwise just return listV.
  If sortnames=true, quicksort the result list.
 */
void mergeLists(int nV, char **listV, int nA, char **listA, char **mlist, bool *isVar) 
{
    int v, a, idx;
    if (sortnames && listattrs && !attrsonly) {
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
        // first add vars then attrs (if ask ed)
        idx = 0;
        if (!attrsonly) {
            for (v=0; v<nV; v++) {
                mlist[idx] = listV[v];
                isVar[idx] = true;
                idx++;
            }
        }
        if (listattrs) {
            for (a=0; a<nA; a++) {
                mlist[idx] = listA[a];
                isVar[idx] = false;
                idx++;
            }
        }
    }
}

int getTypeInfo( enum ADIOS_DATATYPES adiosvartype, int* elemsize)
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
            //*elemsize = 16;
        default:
            return 1;
    }
    return 0;
}

/** Read data of a variable and print 
 * Return: 0: ok, != 0 on error
 */
int readVar(ADIOS_FILE *fp, ADIOS_VARINFO *vi, const char * name, bool timed)
{
    int i,j;
    uint64_t start_t[MAX_DIMS], count_t[MAX_DIMS]; // processed <0 values in start/count
    uint64_t s[MAX_DIMS], c[MAX_DIMS]; // for block reading of smaller chunks
    int tdims;               // number of dimensions including time
    int tidx;                // 0 or 1 to account for time dimension
    uint64_t nelems;         // number of elements to read
    int elemsize;            // size in bytes of one element
    uint64_t st, ct;
    void *data;
    uint64_t sum;           // working var to sum up things
    int  maxreadn;          // max number of elements to read once up to a limit (10MB of data)
    int  actualreadn;       // our decision how much to read at once
    int  readn[MAX_DIMS];   // how big chunk to read in in each dimension?
    int  status;            
    bool incdim;            // used in incremental reading in
    ADIOS_SELECTION * sel;  // boundnig box to read
    int ndigits_dims[32];        // # of digits (to print) of each dimension 

    if (getTypeInfo(vi->type, &elemsize)) {
        fprintf(stderr, "Adios type %d (%s) not supported in bpls. var=%s\n", 
                vi->type, adios_type_to_string(vi->type), name);
        return 10;
    }

    // create the counter arrays with the appropriate lengths
    // transfer start and count arrays to format dependent arrays

    nelems = 1;
    tidx = 0;

    if (timed) {
        if (istart[0] < 0)  // negative index means last-|index|
            st = vi->nsteps+istart[0];
        else
            st = istart[0];
        if (icount[0] < 0)  // negative index means last-|index|+1-start
            ct = vi->nsteps+icount[0]+1-st;
        else
            ct = icount[0];

        if (verbose>2) 
            printf("    j=0, st=%" PRIu64 " ct=%" PRIu64 "\n", st, ct);

        start_t[0] = st;
        count_t[0] = ct;
        nelems *= ct;
        if (verbose>1) 
            printf("    s[0]=%" PRIu64 ", c[0]=%" PRIu64 ", n=%" PRIu64 "\n",
                    start_t[0], count_t[0], nelems);
        
        tidx = 1;
    }
    tdims = vi->ndim + tidx;

    for (j=0; j<vi->ndim; j++) {
        if (istart[j+tidx] < 0)  // negative index means last-|index|
            st = vi->dims[j]+istart[j+tidx];
        else
            st = istart[j+tidx];
        if (icount[j+tidx] < 0)  // negative index means last-|index|+1-start
            ct = vi->dims[j]+icount[j+tidx]+1-st;
        else
            ct = icount[j+tidx];

        if (verbose>2) 
            printf("    j=%d, st=%" PRIu64 " ct=%" PRIu64 "\n", j+tidx, st, ct);

        start_t[j+tidx] = st;
        count_t[j+tidx] = ct;
        nelems *= ct;
        if (verbose>1) 
            printf("    s[%d]=%" PRIu64 ", c[%d]=%" PRIu64 ", n=%" PRIu64 "\n",
                    j+tidx, start_t[j+tidx], j+tidx, count_t[j+tidx], nelems);
    }

    if (verbose>1) {
        printf(" total size of data to read = %" PRIu64 "\n", nelems*elemsize);
    }

    print_slice_info(vi->ndim, vi->dims, timed, vi->nsteps, start_t, count_t);

    maxreadn = MAX_BUFFERSIZE/elemsize;
    if (nelems < maxreadn)
        maxreadn = nelems;

    // special case: string. Need to use different elemsize
    if (vi->type == adios_string) {
        if (vi->value)
            elemsize = strlen(vi->value)+1;
        maxreadn = elemsize;
    }

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
    for (i=tdims-1; i>=0; i--) {
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
    if (verbose>1) printf("    read %d elements at once, %" PRId64 " in total (nelems=%" PRId64 ")\n", actualreadn, sum, nelems);


    // init s and c
    // and calculate ndigits_dims
    for (j=0; j<tdims; j++) {
        s[j]=start_t[j];
        c[j]=readn[j];

        ndigits_dims[j] = ndigits (start_t[j]+count_t[j]-1); // -1: dim=100 results in 2 digits (0..99)
    }

    // read until read all 'nelems' elements
    sum = 0;
    while (sum < nelems) {

        // how many elements do we read in next?
        actualreadn = 1;
        for (j=0; j<tdims; j++) 
            actualreadn *= c[j];

        if (verbose>2) {
            printf("adios_read_var name=%s ", name);
            PRINT_DIMS64("  start", s, tdims, j);
            PRINT_DIMS64("  count", c, tdims, j);
            printf("  read %d elems\n", actualreadn);
        }

        // read a slice finally
        sel = adios_selection_boundingbox (vi->ndim, s+tidx, c+tidx);
        if (timed) {
            status = adios_schedule_read_byid (fp, sel, vi->varid, s[0], c[0], data); 
        } else {
            status = adios_schedule_read_byid (fp, sel, vi->varid, 0, 1, data); 
        }

        if (status < 0) {
            fprintf(stderr, "Error when scheduling variable %s for reading. errno=%d : %s \n", name, adios_errno, adios_errmsg());
            free(data);
            return 11;
        }

        status = adios_perform_reads (fp, 1); // blocking read performed here
        if (status < 0) {
            fprintf(stderr, "Error when reading variable %s. errno=%d : %s \n", name, adios_errno, adios_errmsg());
            free(data);
            return 11;
        }

        //if (verbose>2) printf("  read %" PRId64 " bytes\n", bytes_read);

        // print slice
        print_dataset(data, vi->type, s, c, tdims, ndigits_dims); 

        // prepare for next read
        sum += actualreadn;
        incdim=true; // largest dim should be increased 
        for (j=tdims-1; j>=0; j--) {
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
    print_endline();


    free(data);
    return 0;
}


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

void print_slice_info(int ndim, uint64_t *dims, int timed, int nsteps, uint64_t *s, uint64_t *c)
{
    // print the slice info in indexing is on and 
    // not the complete variable is read
    int i;
    bool isaslice = false;
    int tidx = (timed == true);
    int tdim = ndim + tidx;
    if (timed) {
        if (c[0] < nsteps) 
            isaslice = true;
    }
    for (i=0; i<ndim; i++) {
        if (c[i+tidx] < dims[i])
            isaslice = true;
    }
    if (isaslice) {
        fprintf(outf,"%c   slice (%" PRId64 ":%" PRId64, commentchar, s[0], s[0]+c[0]-1);
        for (i=1; i<tdim; i++) {
            fprintf(outf,", %" PRId64 ":%" PRId64, s[i], s[i]+c[i]-1);
        }
        fprintf(outf,")\n");
    }
}

int print_data_as_string(void * data, int maxlen, enum ADIOS_DATATYPES adiosvartype)
{
    char *str = (char *)data;
    int len = maxlen;
    switch(adiosvartype) {
        case adios_unsigned_byte:
        case adios_byte:
        case adios_string:
            while ( str[len-1] == 0) { len--; }   // go backwards on ascii 0s
            if (len < maxlen) {
                // it's a C string with terminating \0
                fprintf(outf,"\"%s\"", str);
            } else {
                // fortran VARCHAR, lets trim from right padded zeros
                while ( str[len-1] == ' ') {len--;}
                fprintf(outf,"\"%*.*s\"", len, len, (char *) data);
                if (len < maxlen)
                    fprintf(outf," + %d spaces", maxlen-len);
            }
            break;
        default:
            fprintf(stderr, "Error in bpls code: cannot use print_data_as_string() for type \"%s\"\n", 
                    adios_type_to_string(adiosvartype));
            return -1;
            break;
    }
    return 0;
}

int print_data_characteristics(void * min, void * max, double * avg, double * std_dev, enum ADIOS_DATATYPES adiosvartype, bool allowformat)
{
    bool f = formatgiven && allowformat;

    switch(adiosvartype) {
        case adios_unsigned_byte:
            if (min) fprintf(outf,(f ? format : "%10hhu  "), * ((unsigned char *) min));
            else fprintf(outf, "      null  ");
            if (max) fprintf(outf,(f ? format : "%10hhu  "), * ((unsigned char *) max));
            else fprintf(outf, "      null  ");
            if (avg) fprintf(outf, "%10.2f  ", * avg);
            else fprintf(outf, "      null  ");
            if (std_dev) fprintf(outf, "%10.2f  ", * std_dev);
            else fprintf(outf, "      null  ");
            break;
        case adios_byte:
            if (min) fprintf(outf,(f ? format : "%10hhd  "), * ((char *) min));
            else fprintf(outf, "      null  ");
            if (max) fprintf(outf,(f ? format : "%10hhd  "), * ((char *) max));
            else fprintf(outf, "      null  ");
            if (avg) fprintf(outf, "%10.2f  ", * avg);
            else fprintf(outf, "      null  ");
            if (std_dev) fprintf(outf, "%10.2f  ", * std_dev);
            else fprintf(outf, "      null  ");
            break;
        case adios_string:
        case adios_string_array:
            break;

        case adios_unsigned_short:
            if (min) fprintf(outf,(f ? format : "%10hu  "), (* (unsigned short *) min));
            else fprintf(outf, "      null  ");
            if (max) fprintf(outf,(f ? format : "%10hu  "), (* (unsigned short *) max));
            else fprintf(outf, "      null  ");
            if (avg) fprintf(outf, "%10.2f  ", * avg);
            else fprintf(outf, "      null  ");
            if (std_dev) fprintf(outf, "%10.2f  ", * std_dev);
            else fprintf(outf, "      null  ");
            break;
        case adios_short:
            if (min) fprintf(outf,(f ? format : "%10hd  "), (* (short *) min));
            else fprintf(outf, "      null  ");
            if (max) fprintf(outf,(f ? format : "%10hd  "), (* (short *) max));
            else fprintf(outf, "      null  ");
            if (avg) fprintf(outf, "%10.2f  ", * avg);
            else fprintf(outf, "      null  ");
            if (std_dev) fprintf(outf, "%10.2f  ", * std_dev);
            else fprintf(outf, "      null  ");
            break;

        case adios_unsigned_integer:
            if (min) fprintf(outf,(f ? format : "%10u  "), (* (unsigned int *) min));
            else fprintf(outf, "      null  ");
            if (max) fprintf(outf,(f ? format : "%10u  "), (* (unsigned int *) max));
            else fprintf(outf, "      null  ");
            if (avg) fprintf(outf, "%10.2f  ", * avg);
            else fprintf(outf, "      null  ");
            if (std_dev) fprintf(outf, "%10.2f  ", * std_dev);
            else fprintf(outf, "      null  ");
            break;
        case adios_integer:
            if (min) fprintf(outf,(f ? format : "%10d  "), (* (int *) min));
            else fprintf(outf, "      null  ");
            if (max) fprintf(outf,(f ? format : "%10d  "), (* (int *) max));
            else fprintf(outf, "      null  ");
            if (avg) fprintf(outf, "%10.2f  ", * avg);
            else fprintf(outf, "      null  ");
            if (std_dev) fprintf(outf, "%10.2f  ", * std_dev);
            else fprintf(outf, "      null  ");
            break;

        case adios_unsigned_long:
            if (min) fprintf(outf,(f ? format : "%10llu  "), (* (unsigned long long *) min));
            else fprintf(outf, "      null  ");
            if (max) fprintf(outf,(f ? format : "%10llu  "), (* (unsigned long long *) max));
            else fprintf(outf, "      null  ");
            if (avg) fprintf(outf, "%10.2f  ", * avg);
            else fprintf(outf, "      null  ");
            if (std_dev) fprintf(outf, "%10.2f  ", * std_dev);
            else fprintf(outf, "      null  ");
            break;
        case adios_long:
            if (min) fprintf(outf,(f ? format : "%10lld  "), (* (long long *) min));
            else fprintf(outf, "      null  ");
            if (max) fprintf(outf,(f ? format : "%10lld  "), (* (long long *) max));
            else fprintf(outf, "      null  ");
            if (avg) fprintf(outf, "%10.2f  ", * avg);
            else fprintf(outf, "      null  ");
            if (std_dev) fprintf(outf, "%10.2f  ", * std_dev);
            else fprintf(outf, "      null  ");
            break;

        case adios_real:
            if (min) fprintf(outf,(f ? format : "%10.2g  "), (* (float *) min));
            else fprintf(outf, "      null  ");
            if (max) fprintf(outf,(f ? format : "%10.2g  "), (* (float *) max));
            else fprintf(outf, "      null  ");
            if (avg) fprintf(outf, "%10.2g  ", * avg);
            else fprintf(outf, "      null  ");
            if (std_dev) fprintf(outf, "%10.2g  ", * std_dev);
            else fprintf(outf, "      null  ");
            break;
        case adios_double:
            if (min) fprintf(outf,(f ? format : "%10.2g  "), (* (double *) min));
            else fprintf(outf, "      null  ");
            if (max) fprintf(outf,(f ? format : "%10.2g  "), (* (double *) max));
            else fprintf(outf, "      null  ");
            if (avg) fprintf(outf, "%10.2g  ", * avg);
            else fprintf(outf, "      null  ");
            if (std_dev) fprintf(outf, "%10.2g  ", * std_dev);
            else fprintf(outf, "      null  ");
            break;

        case adios_long_double:
            //fprintf(outf,(f ? format : "%g "), ((double *) data)[item]);
            fprintf(outf, "%s", (f ? format : "????????"));
            break;

            // TO DO
            /*
               case adios_complex:
               fprintf(outf,(f ? format : "(%g,i%g) "), ((float *) data)[2*item], ((float *) data)[2*item+1]);
               break;

               case adios_double_complex:
               fprintf(outf,(f ? format : "(%g,i%g)" ), ((double *) data)[2*item], ((double *) data)[2*item+1]);
               break;
             */
        default:
            break;
    } // end switch
    return 0;
}

int print_data(void *data, int item, enum ADIOS_DATATYPES adiosvartype, bool allowformat)
{
    bool f = formatgiven && allowformat;
    if (data == NULL) {
        fprintf(outf, "null ");
        return 0;
    }
    // print next data item into vstr
    switch(adiosvartype) {
        case adios_unsigned_byte:
            fprintf(outf,(f ? format : "%hhu "), ((unsigned char *) data)[item]);
            break;
        case adios_byte:
            fprintf(outf,(f ? format : "%hhd "), ((signed char *) data)[item]);
            break;
        case adios_string:
            fprintf(outf,(f ? format : "%s"), ((char *) data)+item);
            break;
        case adios_string_array:
            // we expect one elemet of the array here
            fprintf(outf,(f ? format : "\"%s\""), *((char **)data+item));
            break;

        case adios_unsigned_short:  
            fprintf(outf,(f ? format : "%hu "), ((unsigned short *) data)[item]);
            break;
        case adios_short:
            fprintf(outf,(f ? format : "%hd "), ((signed short *) data)[item]);
            break;

        case adios_unsigned_integer:
            fprintf(outf,(f ? format : "%u "), ((unsigned int *) data)[item]);
            break;
        case adios_integer:    
            fprintf(outf,(f ? format : "%d "), ((signed int *) data)[item]);
            break;

        case adios_unsigned_long:
            fprintf(outf,(f ? format : "%llu "), ((unsigned long long *) data)[item]);
            break;
        case adios_long:        
            fprintf(outf,(f ? format : "%lld "), ((signed long long *) data)[item]);
            break;

        case adios_real:
            fprintf(outf,(f ? format : "%g "), ((float *) data)[item]);
            break;

        case adios_double:
            fprintf(outf,(f ? format : "%g "), ((double *) data)[item]);
            break;


        case adios_long_double:
            //fprintf(outf,(f ? format : "%g "), ((double *) data)[item]);
            fprintf(outf, "%s", (f ? format : "????????"));
            break;


        case adios_complex:  
            fprintf(outf,(f ? format : "(%g,i%g) "), ((float *) data)[2*item], ((float *) data)[2*item+1]);
            break;

        case adios_double_complex:
            fprintf(outf,(f ? format : "(%g,i%g)" ), ((double *) data)[2*item], ((double *) data)[2*item+1]);
            break;
        default:
            break;
    } // end switch
    return 0;
}

int print_dataset(void *data, enum ADIOS_DATATYPES adiosvartype, 
        uint64_t *s, uint64_t *c, int tdims, int *ndigits)
{
    int i,item, steps;
    char idxstr[128], buf[16];
    uint64_t ids[MAX_DIMS];  // current indices
    bool roll;

    // init current indices
    steps = 1;
    for (i=0; i<tdims; i++) {
        ids[i] = s[i];
        steps *= c[i];
    }

    item = 0; // index to *data 
    // loop through each data item and print value
    while (item < steps) {

        // print indices if needed into idxstr;
        idxstr[0] = '\0'; // empty idx string
        if (nextcol == 0) {
            if (!noindex && tdims > 0) {
                sprintf(idxstr,"    (%*" PRId64, ndigits[0], ids[0]);
                for (i=1; i<tdims; i++) {
                    sprintf(buf,",%*" PRId64, ndigits[i], ids[i]);
                    strcat(idxstr, buf);
                }
                strcat(idxstr,")    ");
            }
        }

        // print item
        fprintf(outf, "%s", idxstr);
        if (printByteAsChar && (adiosvartype == adios_byte || adiosvartype == adios_unsigned_byte)) {
            /* special case: k-D byte array printed as (k-1)D array of strings */
            if (tdims == 0) {
                print_data_as_string(data, steps, adiosvartype);
            } else {
                print_data_as_string(data+item, c[tdims-1], adiosvartype); // print data of last dim as string
                item += c[tdims-1]-1; // will be ++-ed once below
                ids[tdims-1] = s[tdims-1]+c[tdims-1]-1; // will be rolled below
            }
            nextcol = ncols-1; // force new line, will be ++-ed once below
        } else {
            print_data(data, item, adiosvartype, true);
        }

        // increment/reset column index
        nextcol++;
        if (nextcol == ncols) {
            fprintf(outf,"\n");
            nextcol = 0;
        }

        // increment indices
        item++;
        roll = true;
        for (i=tdims-1; i>=0; i--) {
            if (roll) {
                if (ids[i] == s[i]+c[i]-1 ) {
                    // last index in this dimension, roll upward
                    ids[i] = s[i];
                } else {
                    ids[i]++;
                    roll = false;
                }
            }
        }
    }
    return 0;
}

void print_endline(void) 
{
    if (nextcol != 0)
        fprintf(outf,"\n");
    nextcol = 0;
}


int print_decomp(ADIOS_VARINFO *vi)
{
    /* Print block info */
    int i,j,k;
    int ndigits_nsteps = ndigits (vi->nsteps-1);
    if (vi->ndim == 0) 
    {
        // scalars
        for (i=0; i < vi->nsteps; i++) {
            fprintf(outf, "        step %*d: ", ndigits_nsteps, i);
            fprintf(outf, "%d instances available\n", vi->nblocks[i]);
        }
    } 
    else 
    {
        // arrays
        for (i=0; i < /*vi->nsteps*/1; i++) { // For now, just look at the first step xx
            fprintf(outf, "    decomposition: [");
            fprintf(outf,"\n");
            for (j=0; j < vi->nblocks[i]; j++) {
                fprintf(outf,"        [");
                for (k=0; k < vi->ndim; k++) {
                    fprintf(outf, "[%" PRId64 ",%" PRId64 "]",
                            vi->blockinfo[j].start[k],
                            vi->blockinfo[j].start[k] + vi->blockinfo[j].count[k]-1);
                    if (k < vi->ndim-1)
                        fprintf(outf, ", ");
                }
                fprintf(outf, "]");
                if (j < vi->nblocks[i]-1)
                    fprintf(outf, ",");
                fprintf(outf, "\n");
            }
            fprintf(outf,"    ]\n");
        }
    }
    return 0;
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
