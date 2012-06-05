#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <limits.h>
#include <errno.h>

#include "util.h"


static char * remove_whitespace (char *start, char *end) 
{
    char *s = start;
    char *e = end;
    int orig_len = (int) (e-s);
    int final_len;
    char *res;
    // remove front whitespace (but do not go far beyond the end)
    while (s <= e && 
           (*s==' ' || *s=='\t' || *s=='\n')
          ) s++;
    if (s <= e) { // there is some text 
        // remove tail whitespace
        while (s <= e && 
               (*e==' ' || *e=='\t' || *e=='\n')
              ) e--;
        // create result 
        final_len = e - s + 1; //  length of result 
        if (final_len > 0) {
            res = (char *) malloc (final_len + 1); // allocate space s..e and \0
            memcpy(res, s, final_len);
            res[final_len] = 0;
        } else {
            // "   = something" patterns end here
            res = NULL;
        }
    } else {
        // no non-whitespace character found
        res = NULL;
    }
    return res;
}


/* Split a line at = sign into name and value pair
   Remove " ", TAB and Newline from around names and values
   Return NULL for name and value if there is no = sign in line
   Return newly allocated strings otherwise
   Used by: esimmon_internal_text_to_name_value_pairs
 */
static void splitnamevalue (const char * line, int linelen,  char **name, char **value)
{
    char *equal; // position of first = sign in line

    equal = strchr (line, '=');
    if (equal && equal != line) {
        /* 1. name */
        // from first char to before =
        *name = remove_whitespace ((char*)line, equal-1);
        //printf ("      --name=[%s]\n", *name);
        /* 2. value */
        // from after = to the last character of line
        *value = remove_whitespace (equal+1, (char*)line+linelen-1);
        //printf ("      --value=[%s]\n", *value);

    } else if (equal != line) {
        /* check if it as name without = value statement */
        *name = remove_whitespace ((char*)line, (char*)line+linelen-1);
        //printf ("      --name only=[%s]\n", *name);
        *value = NULL;
    } else { 
        // funny text starting with =. E.g. "=value" 
        *name = NULL;
        *value = NULL;
    }
}

PairStruct * text_to_name_value_pairs (char * text)
{
    /* Process a multi-line and/or ;-separated text and create a list
       of name=value pairs from each line which has a 
           name = value 
       pattern. Whitespaces are removed. 
         "X = 1
          Y = 2"  
       is not valid because of missing ';', but
          "X=1; Y=5;
          Z=apple"  
       is valid
    */
    char *name, *value; 
    char *item, *delim;
    int len;
    char line[256];
    PairStruct *res = NULL, *last = NULL, *pair;

    if (!text) return res;

    item  = text; 
    while (item) {
        delim = strchr (item, ';');
        if (delim) 
            len = (int) (delim-item); 
        else 
            len = strlen (item);

        strncpy (line, item, len);
        line[len] = '\0';

        //printf ("    --Line=[%s]\n", line);
        splitnamevalue(line, len, &name, &value);
        if (name) {
            pair = (PairStruct *) malloc (sizeof(PairStruct));
            pair->name = name;
            pair->value = value;
            pair->next = NULL;
            if (last) {
                last->next = pair;
                last = pair;
            } else {
                res = pair; 
                last = pair;
            }
        }
        if (delim && delim+1 != 0)
            item = delim+1;
        else
            item = NULL;
    }
    return res;
}


void free_name_value_pairs (PairStruct * pairs)
{
    PairStruct *p;
    while (pairs) {
        free(pairs->name);
        free(pairs->value);
        p = pairs;
        pairs=pairs->next;
        free(p);
    }
}



