/*
 * parse_config() function parse XML config file. 
 * The code is larged borrowed from adios_parse_config() in ADIOS.
 * Right now it only read the host-language parameter.
 *
 */
#include <string.h>
#include "adios_internals.h"
#include "hw-utils.h"

// xml parser
#include <mxml.h>

/*
 * parse_config() function parse XML config file specified by config_filename. 
 * The parsed information is stored in config_params.
 * It returns 0 if no error is encountered and -1 otherwise. 
 */
int parse_config(char *config_filename, adios_config_file *config_params)
{
    FILE * fp = 0;
    mxml_node_t * doc = NULL;
    mxml_node_t * node = NULL;
    mxml_node_t * root = NULL;
    int saw_datagroup = 0;
    int saw_method = 0;
    int saw_buffer = 0;

    fp = fopen (config_filename, "r");
    if (!fp)
    {
        fprintf (stderr, "Error in parsing config file: missing config file %s\n", config_filename);

        return -1;
    }
    doc = mxmlLoadFile (NULL, fp, MXML_TEXT_CALLBACK);
    fclose (fp);
    if (!doc)
    {
        fprintf (stderr, "Error in parsing config file: unknown error in parsing XML\n"
                         "Did you remember to start the file with\n"
                         "<?xml version=\"1.0\"?>");

        return -1;
    }

    root = mxmlWalkNext (doc, doc, MXML_DESCEND); // get rid of the xml version
    root = mxmlWalkNext (root, doc, MXML_DESCEND); // should be our root tag
    
    while (!strncmp (root->value.element.name, "!--", 3)) // skip comments
    {
        root = mxmlWalkNext (root, doc, MXML_NO_DESCEND);
    }

    if (strcmp (root->value.element.name, "adios-config"))
    {
        fprintf (stderr, "Error in parsing config file: invalid root xml element: %s\n"
                ,root->value.element.name
                );

        mxmlRelease (doc);

        return -1;
    }
    else
    {
        const char * host_language = NULL;
        host_language = mxmlElementGetAttr (root, "host-language");
        if (!host_language)
        {
            host_language = "Fortran";
        }

        if (!strcasecmp (host_language, "Fortran"))
        {
            config_params->adios_host_language_fortran = USE_FORTRAN; 
        }
        else
        {
            if (!strcasecmp (host_language, "C")) 
            {
                config_params->adios_host_language_fortran = USE_C;
            }
            else
            {
                fprintf (stderr, "Error in parsing config file: invalid host-language %s"
                        ,host_language
                        );

                mxmlRelease (doc);

                return -1;
            }
        }
    }

    mxmlRelease (doc);

    return 0;
}

