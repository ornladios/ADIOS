#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "hw-utils.h"

void print_usage() 
{
    printf("bp2h5:\n");
    printf("    Convert a bp file to h5 file.\n");
    printf("\nUSAGE: bp2h5 [OPTION] bp_file [h5_file]\n\n");
    printf("    bp2h5 converts a bp file named by bp_file to a h5 file named by h5_file.\n");
    printf("    If h5_file is not specified, the generated h5 file is named by the basename\n"
           "    of bp_file suffixed by \".h5\".\n");
    printf("\nOPTION:\n");
    printf("    --scalar_as_array\n");
    printf("        (Optional) Write a scalar variable or attribute in a single-element \n"
           "        array. If not specified, a scalar variable/attribute is written in a\n"
           "        scalar dataspace.\n");
    printf("    -V, --verbose\n");
    printf("        (Optional) Print detailed information during conversion.\n");
    printf("    -h, --help\n");
    printf("        (Optional) Print usage information and exit.\n");
    printf("\n");
}

/*
 * parse_cmdline() function parses command line arguments. 
 * It returns 0 if no error is encountered and -1 otherwise.
 */
int parse_cmdline(int argc, char **argv
                 ,char **bp_filename
				 ,char **h5_filename
				 ,enum scalar_convention *scalar_as_array
				 ,enum verbose_level *verb                   
                 ) 
{
    int i = 1;
    int found_bp_file = 0;
    int found_h5_file = 0;
    int found_config_file = 0;
    
    *scalar_as_array = USE_SCALAR;
    *verb = NO_INFO;

    if(argc < 2) {
        print_usage();
        return -1;
    }
    
    while (i < argc) {
        if(!strcmp(argv[i], "--scalar-as_array")) {
            *scalar_as_array = USE_SINGLE_ELE_ARRAY;
        }
        else if(!strcmp(argv[i], "--verbose") || !strcmp(argv[i], "-V")) {
            *verb = DEBUG_INFO;
        }
        else if(!strcmp(argv[i], "-h") || !strcmp(argv[i], "--help")) {
            print_usage();
            exit(0);
        }
        else if(!found_bp_file) {
            *bp_filename = argv[i];
            found_bp_file = 1;
        }
        else if(found_bp_file && !found_h5_file) {
            *h5_filename = argv[i];
            found_h5_file = 1;
        }
        else {
            // unknown arguments
            fprintf(stderr, "Error in parsing command line: unknown argument %s\n\n", argv[i]);
            print_usage();
            return -1;
        }
        i ++;
    }
    
    if(!found_bp_file) {
        fprintf(stderr, "Error in parsing command line: bp_file not provided\n");
        print_usage();
        return -1;
    }

    return 0;
}

int main (int argc, char ** argv)
{
    char *bp_filename = NULL;
    char *h5_filename = NULL;
    enum scalar_convention scalar_as_array;
    enum verbose_level verb;
    
    // parse cmdline options
    if(parse_cmdline(argc, argv, &bp_filename, &h5_filename, &scalar_as_array, &verb)) {
        return -1;
    }

    initialize_bp2h5(USE_FORTRAN, USE_FORTRAN, USE_FORTRAN, USE_FORTRAN, scalar_as_array, verb);

    // generate h5 file
    return hw_makeh5(bp_filename, h5_filename);
}
