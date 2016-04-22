#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "core/strutil.h"

/* Test the the PairStruct functions in ADIOS.
 * Convert text to name-value pairs.
 */

// Single pair
const char text1[] = "a=5";
PairStruct p1 = {.name="a", .value="5", .next=NULL};

// Two pairs
const char text2[] = "b=4 ; othername=othervalue;";
PairStruct p2_2 = {.name="othername", .value="othervalue", .next=NULL};
PairStruct p2   = {.name="b", .value="4", .next=&p2_2};

// Single pair in multiple lines
const char text3[] = "\n  a  \n"
                     "\n  =  \n"
                     "\n  5  \n";
PairStruct p3 = {.name="a", .value="5", .next=NULL};

// Multiple pairs in multiple lines, and a TAB
const char text4[] = "a=1; \n"
                     "\tb=2; \n"
                     "c=3";
PairStruct p4_c = {.name="c", .value="3", .next=NULL};
PairStruct p4_b = {.name="b", .value="2", .next=&p4_c};
PairStruct p4   = {.name="a", .value="1", .next=&p4_b};

// Single pair, value in ""
const char text5[] = "stringvalue=\"hello\"";
PairStruct p5 = {.name="stringvalue", .value="\"hello\"", .next=NULL};

// Single pair, value is a long string in "" containing name-value pairs delimited by ;
const char text6[] = "stringvalue=\"a=1;b=2;c=3\"";
PairStruct p6 = {.name="stringvalue", .value="\"a=1;b=2;c=3\"", .next=NULL};

// Single pair, NULL value
const char text7[] = "novalue";
PairStruct p7 = {.name="novalue", .value=NULL, .next=NULL};

// Multiple pairs, some value with " and ;, and multiple lines
const char text8[] =
        "method=MPI_AGGREGATE;parameters=\"num_aggregators=3;num_ost=2;\";"
        "method=MPI;parameters=\"\";"
        "method=POSIX;parameters=\"debug\";verbose";
PairStruct p8_verb = {.name="verbose",      .value=NULL, .next=NULL};
PairStruct p8_p3   = {.name="parameters",   .value="\"debug\"",                         .next=&p8_verb};
PairStruct p8_m3   = {.name="method",       .value="POSIX",                             .next=&p8_p3};
PairStruct p8_p2   = {.name="parameters",   .value="\"\"",                              .next=&p8_m3};
PairStruct p8_m2   = {.name="method",       .value="MPI",                               .next=&p8_p2};
PairStruct p8_p1   = {.name="parameters",   .value="\"num_aggregators=3;num_ost=2;\"",  .next=&p8_m2};
PairStruct p8      = {.name="method",       .value="MPI_AGGREGATE",                     .next=&p8_p1};

// Single pair, value with opening " and without closing "
const char text9[] = "badexample=\"no closing quote    ";
PairStruct p9 = {.name="badexample", .value="\"no closing quote", .next=NULL};

void print_name_value_pairs (PairStruct * pairs)
{
    printf("pairs:\n");
    while (pairs) {
        printf ("   %s = %s\n", pairs->name, pairs->value);
        pairs=pairs->next;
    }
}

int dotest (const char * text, PairStruct *p)
{
    int nerrors = 0;
    printf("text=[%s]\n", text);
    PairStruct *q = a2s_text_to_name_value_pairs(text);

    /* Get each element and check value */
    PairStruct *p1 = p;
    PairStruct *q1 = q;
    while (p1 != NULL && q1 != NULL)
    {
        if (strcmp(p1->name,q1->name))
        {
            printf("   ERROR: name [%s] does not match expected name [%s]\n", q1->name, p1->name);
            nerrors++;
        }

        if (p1->value == NULL && q1->value == NULL) {
            // this is okay
            ;
        }
        else if (p1->value == NULL || q1->value == NULL) {
            printf("   ERROR: for name [%s], value [%s] does not match expected value [%s]\n", q1->name, q1->value, p1->value);
                        nerrors++;
        }
        else if (strcmp(p1->value,q1->value))
        {
            printf("   ERROR: for name [%s], value [%s] does not match expected value [%s]\n", q1->name, q1->value, p1->value);
            nerrors++;
        }
        p1 = p1->next;
        q1 = q1->next;
    }
    print_name_value_pairs (q);
    printf ("------------------------------------------------------\n");
    return nerrors;
}

int main (int argc, char ** argv)
{
    int nerrors = 0;
    printf("\n============= Text to name-value pairs conversion test =========\n");

    nerrors += dotest (text1, &p1);
    nerrors += dotest (text2, &p2);
    nerrors += dotest (text3, &p3);
    nerrors += dotest (text4, &p4);
    nerrors += dotest (text5, &p5);
    nerrors += dotest (text6, &p6);
    nerrors += dotest (text7, &p7);
    nerrors += dotest (text8, &p8);
    nerrors += dotest (text9, &p9);

    printf("\nNumber of errors in this test: %d\n", nerrors);
    return nerrors;
}

