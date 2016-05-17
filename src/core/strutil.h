#ifndef STRUTIL_H_
#define STRUTIL_H_

// copy an array of strings with allocation, return pointer
// also return the sum of string lengths in 'total_length'
char ** dup_string_array (char ** v, int nelems, int * total_length);
void free_string_array (char ** v, int nelems);

void alloc_namelist (char ***namelist, int length);
void free_namelist (char **namelist, int length);

char * trim_spaces (const char * str);
void tokenize_dimensions (const char * str, char *** tokens, int * count);

/*******************************************************
   Processing parameter lists
**********************************************************/
/*
   Process a ;-separated and possibly multi-line text and 
   create a list of name=value pairs from each 
   item which has a "name=value" pattern. Whitespaces are removed. 
   Input is not modified. Space is allocated;
   Also, simple "name" or "name=" patterns are processed and 
   returned with value=NULL. 
*/
struct PairStruct {
    char * name;
    char * value;
    struct PairStruct * next;
};
typedef struct PairStruct PairStruct;

PairStruct * text_to_name_value_pairs (const char * text);
void free_name_value_pairs (PairStruct * pairs);


#endif
