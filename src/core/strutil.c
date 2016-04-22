#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <limits.h>
#include <errno.h>
#include <assert.h>
#include <ctype.h>

#include "config.h"
#include "core/strutil.h"

// copy an array of strings with allocation, return pointer
// also return the sum of string lengths in 'total_length'
char ** a2s_dup_string_array (char ** v, int nelems, int * total_length) 
{
    *total_length = 0;

    if (!v || nelems < 1)
        return NULL;

    char ** p = malloc (nelems*sizeof(char*));
    if (!p) return NULL;

    int i, len;
    for (i=0; i<nelems; i++) {
        if (v[i]) {
            len = strlen (v[i]) + 1;
            p[i] = malloc (len*sizeof(char));
            if (p[i])
                memcpy (p[i], v[i], len);
            *total_length += len;
        } else {
            p[i] = NULL;
        }
    }
    return p;
}

void a2s_free_string_array (char ** v, int nelems) 
{
    int i;
    for (i=0; i<nelems; i++) {
        if (v[i]) free (v[i]);
        v[i] = 0;
    }
    free (v);
}

void a2s_alloc_namelist (char ***namelist, int length)
{
    int j;

    *namelist = (char **) malloc(length*sizeof(char*));
    for (j=0;j<length;j++)
        (*namelist)[j] = (char *) malloc(255);

    return;
}

void a2s_free_namelist (char **namelist, int length)
{
    int i;
    if (namelist) {
        for (i=0;i<length;i++) {
            if(namelist[i])
                free(namelist[i]);
            namelist[i] = NULL;
        }
        free(namelist);
    }
    return;
}

// remove leading white spaces
// it returns a pointer inside str, it does not allocate new memory
char * a2s_trimL (char * str)
{
    if (!str)
        return str;
    char * b = str;
    while (isspace(*b))
    {
        b++;
    }
    return b;
}

// remove trailing white spaces
// it returns str, it does not allocate new memory, just shortens the strings
char * a2s_trimR (char * str)
{
    if (!str)
        return str;
    int len = strlen(str);
    if (!len)
        return str;
    char * t = str+len-1;
    while (isspace(*t))
    {
        *t = '\0';
        t--;
    }
    return str;
}


char * a2s_trimLR (char * str)
{
    if (!str)
        return str;
    int len = strlen(str);
    if (!len)
        return str;

    // trim front
    char * b = str;
    while (isspace(*b))
    {
        b++;
    }

    // trim trailing whitespaces
    char * t = str+len-1;
    while (isspace(*t))
    {
        *t = '\0';
        t--;
    }
    return b;
}

void a2s_tokenize_dimensions (const char * str, char *** tokens, int * count)
{
    *count = 0;
    *tokens = 0;
    if (!str)
        return;

    char * dims[32];
    char * save_str = strdup (str);
    char * t = save_str;
    int i;

    for (t = strtok(save_str, ",");
         t;
         t = strtok(NULL, ","))
    {
        t = a2s_trimLR(t);
        dims[*count] = strdup (t);
        (*count)++;
    }

    if (*count)
    {
        *tokens = (char **) malloc (sizeof (char **) * *count);
        for (i = 0; i < *count; i++)
        {
            (*tokens) [i] = strdup (dims[i]);
        }
    }

    free (save_str);
}

void a2s_cleanup_dimensions (char ** tokens, int count)
{
    int i;
    for (i = 0; i < count; i++)
    {
        free (tokens[i]);
    }
    if (tokens)
        free (tokens);
}

void trim_spaces (char * str)
{
    char * t = str, * p = NULL;
    while (*t != '\0')
    {
        if (*t == ' ')
        {
            p = t + 1;
            strcpy (t, p);
        }
        else
            t++;
    }

}

/*******************************************************
   Processing parameter lists
**********************************************************/
static char * remove_whitespace (char *start, char *end) 
{
    char *s = start;
    char *e = end;
    //int orig_len = (int) (e-s);
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

PairStruct * a2s_text_to_name_value_pairs (const char * text)
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
    char *item, *delim, *quote1, *quote2;
    int len;
    char line[256];
    PairStruct *res = NULL, *last = NULL, *pair;

    if (!text) return res;

    item  = (char *)text; 
    while (item) {
        /* Before cutting the string up at ;, look for "
         * A value might be a long string containing ;
         * e.g.  methodparams = "aggr=5;ost=3;verbose=1"; method = MPI_AGGREGATE;
         */
        quote1 = strchr (item, '"');
        delim  = strchr (item, ';');
        if (quote1 && delim && delim > quote1) // ; is after opening "
        {
            if (quote1+1) // string continues
            {
                quote2 = strchr (quote1+1, '"'); // find closing "

                if (quote2)
                {
                    delim = strchr (quote2, ';');
                }
                else
                {
                    // This is illegal text, no closing "
                    // let's just pass the whole thing so result will look like e.g. name = "this is an unclosed string
                }
            }
        }

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

void a2s_free_name_value_pairs (PairStruct * pairs)
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
