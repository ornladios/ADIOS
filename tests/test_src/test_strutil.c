#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "core/strutil.h"

/* Test the the strutil functions in ADIOS.
 * trimL, trimR, trimLR, tokenize_dimensions, cleanup_dimensions
 */

struct TestRecord {
    char * str;
    char * trimL;
    char * trimR;
    char * trimLR;
};
typedef struct TestRecord TestRecord;

// Single letter
const TestRecord tr1 = {.str="a ", .trimL="a ", .trimR="a", .trimLR="a"};
const char * dims1[1] = {"a"};

// Two elements
const TestRecord tr2 = {.str    = " abc , efg  ",
                        .trimL  =  "abc , efg  ",
                        .trimR  = " abc , efg",
                        .trimLR =  "abc , efg"
                       };
const char * dims2[2] = {"abc","efg"};


// Single element in multiple lines
const TestRecord tr3 = {.str    = " \n \t a \t \n ",
                        .trimL  =        "a \t \n ",
                        .trimR  = " \n \t a",
                        .trimLR =        "a"
                       };
const char * dims3[1] = {"a"};

// Multiple elements in multiple lines
const TestRecord tr4 = {.str    = " \t a\tb \t, c\td ",
                        .trimL  =     "a\tb \t, c\td ",
                        .trimR  = " \t a\tb \t, c\td",
                        .trimLR =     "a\tb \t, c\td"
                       };
const char * dims4[2] = {"a\tb", "c\td"};

// Empty string
const TestRecord tr5 = {.str    = "",
                        .trimL  = "",
                        .trimR  = "",
                        .trimLR = ""
                       };

void print_elements (char ** dims, int ndims)
{
    int i;
    printf("  elements   = {");
    for (i=0; i < ndims; i++)
    {
        printf ("%s", dims[i]);
        if (i < ndims - 1)
            printf (",");
    }
    printf ("}\n");
}

int dotest (const TestRecord * tr, const char ** dims, int ndims)
{
    int nerrors = 0, i;
    char * res, *s;
    printf("text = [%s]\n", tr->str);

    s = strdup (tr->str);
    res = a2s_trimL(s);
    printf("  left  trim = [%s]\n", res);
    if (!res || strcmp (res, tr->trimL))
    {
        printf("   ERROR: left trim [%s] does not match expected name [%s]\n", res, tr->trimL);
        nerrors++;
    }
    free (s);

    s = strdup (tr->str);
    res = a2s_trimR(s);
    printf("  right trim = [%s]\n", res);
    if (!res || strcmp (res, tr->trimR))
    {
        printf("   ERROR: right trim [%s] does not match expected name [%s]\n", res, tr->trimR);
        nerrors++;
    }
    free (s);

    s = strdup (tr->str);
    res = a2s_trimLR(s);
    printf("  full  trim = [%s]\n", res);
    if (!res || strcmp (res, tr->trimLR))
    {
        printf("   ERROR: trim [%s] does not match expected name [%s]\n", res, tr->trimLR);
        nerrors++;
    }
    free (s);

    char ** d;
    int count;
    a2s_tokenize_dimensions (tr->str, &d, &count);
    print_elements (d, count);
    if (ndims != count)
    {
        printf("   ERROR: number of comma-separated elements [%d] does not match expected number [%d]\n", count, ndims);
        nerrors++;
    }
    for (i=0; i<ndims; i++)
    {
        if ( strcmp(d[i], dims[i]))
        {
            printf("   ERROR: comma-separated element No %d [%s] does not match expected element [%s]\n",
                    i, d[i], dims[i]);
            nerrors++;
        }
    }

    a2s_cleanup_dimensions (d, count);
    printf ("------------------------------------------------------\n");
    return nerrors;
}

int main (int argc, char ** argv)
{
    int nerrors = 0;
    printf("\n============= Test string utility functions =========\n");

    nerrors += dotest (&tr1, dims1, 1);
    nerrors += dotest (&tr2, dims2, 2);
    nerrors += dotest (&tr3, dims3, 1);
    nerrors += dotest (&tr4, dims4, 2);
    nerrors += dotest (&tr5, NULL, 0);

    printf("\nNumber of errors in this test: %d\n", nerrors);
    return nerrors;
}

