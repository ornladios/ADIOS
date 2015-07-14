#ifndef ADIOS_LINK_H
#define ADIOS_LINK_H

#ifdef __cplusplus
extern "C" {
#endif

enum ADIOS_LINK_TYPE {
    LINK_VAR = 1,
    LINK_IMAGE = 2
    /* expand supported link types here */
};

typedef struct
{
    int id;
    char * name;
    int nrefs;
    enum ADIOS_LINK_TYPE * type;
    char ** ref_names;    /* the linked variable name referred from this var */
    char ** ref_files;    /* full path, 0 means link from the same file, otherwise link from externel file */
} ADIOS_LINK;

#ifdef __cplusplus
}
#endif

#endif
