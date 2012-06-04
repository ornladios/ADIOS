#ifndef PARAMS_H_
#define PARAMS_H_


struct PairStruct {
    char * name;
    char * value;
    struct PairStruct * next;
};
typedef struct PairStruct PairStruct;

/* Process a ;-separated and possibly multi-line text and 
   create a list of name=value pairs from each 
   item which has a "name=value" pattern. Whitespaces are removed. 
   Input is not modified. Space is allocated;
   Also, simple "name" or "name=" patterns are processed and 
   returned with value=NULL. 
*/
PairStruct * text_to_name_value_pairs (char * text);
void free_name_value_pairs (PairStruct * pairs);

#endif
