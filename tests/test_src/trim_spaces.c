#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "core/strutil.h"

/* Test the the trim_spaces function in ADIOS.
 */

const char text1[] = "abc";
const char result1[]  = "abc";

const char text2[] = " abc";
const char result2[]  = "abc";

const char text3[] = "abc ";
const char result3[]  = "abc";

const char text4[] = " abc ";
const char result4[]  = "abc";

const char text5[] = "a bc";
const char result5[]  = "abc";

const char text6[] = "ab c";
const char result6[]  = "abc";

const char text7[] = "  a  b  c  ";
const char result7[]  = "abc";

int dotest (const char * text, const char *expected)
{
    int nerrors = 0;
    printf("text=[%s]\n", text);
    char * str = a2s_trim_spaces (text);

    /* Check value */
    if (strcmp(expected,str))
    {
        printf("   ERROR: text [%s] does not match expected result [%s]\n", str, expected);
        nerrors++;
    }

    return nerrors;
}

int main (int argc, char ** argv)
{
    int nerrors = 0;
    printf("\n============= Trim spaces test =========\n");

    nerrors += dotest (text1, result1);
    nerrors += dotest (text2, result2);
    nerrors += dotest (text3, result3);
    nerrors += dotest (text4, result4);
    nerrors += dotest (text5, result5);
    nerrors += dotest (text6, result6);
    nerrors += dotest (text7, result7);

    printf("\nNumber of errors in this test: %d\n", nerrors);
    return nerrors;
}

