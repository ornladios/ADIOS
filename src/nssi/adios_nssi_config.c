/*
 * ADIOS is freely available under the terms of the BSD license described
 * in the COPYING file in the top level directory of this source distribution.
 *
 * Copyright (c) 2008 - 2009.  UT-BATTELLE, LLC. All rights reserved.
 */

#include <string.h>
#include <stdlib.h>
#include <stdint.h>
#include <sys/stat.h> /* struct stat */

// xml parser
#include <mxml.h>

#include "public/adios_mpi.h"
#include "nssi_client.h"
#include "adios_nssi_args.h"
#include "adios_nssi_config.h"

// this macro makes getting the attributes easier
// fix the bgp bugs
#define GET_ATTR(n,attr,var,en)                              \
if (!strcasecmp (n, attr->name)) {                           \
    if (!var)                                                \
    {                                                        \
        var = attr->value;                                   \
        continue;                                            \
    }                                                        \
    else                                                     \
    {                                                        \
        fprintf (stderr, "xml: duplicate attribute %s on %s (ignored)",n,en); \
        continue;                                            \
    }                                                        \
}


static int parseGroup(mxml_node_t * node, struct adios_nssi_config *config)
{
    mxml_node_t * n;
    int i;
    int service_count=0;
    int service_num=0;

//    printf("enter parseGroup\n");

    for (n = mxmlWalkNext (node, node, MXML_DESCEND)
        ;n
        ;n = mxmlWalkNext (n, node, MXML_NO_DESCEND)
        )
    {
        if (n->type != MXML_ELEMENT)
        {
            continue;
        }

        if (!strcasecmp (n->value.element.name, "staging-service"))
        {
            service_count++;
        }
    }

    config->num_servers=service_count;
    config->nssi_server_ids=(nssi_remote_pid *)calloc(service_count, sizeof(nssi_remote_pid));


    service_num=0;
    for (n = mxmlWalkNext (node, node, MXML_DESCEND)
        ;n
        ;n = mxmlWalkNext (n, node, MXML_NO_DESCEND)
        )
    {
        if (n->type != MXML_ELEMENT)
        {
            continue;
        }

        if (!strcasecmp (n->value.element.name, "staging-service"))
        {
            mxml_node_t * n1;   // used for global_bounds

            const char *nid = 0;
            const char *pid= 0;
            const char *hostname = 0;
            const char *port = 0;

            for (i = 0; i < n->value.element.num_attrs; i++)
            {
                mxml_attr_t * attr = &n->value.element.attrs [i];

                GET_ATTR("nid",attr,nid,"var")
                GET_ATTR("pid",attr,pid,"var")

                GET_ATTR("hostname",attr,hostname,"var")
                GET_ATTR("port",attr,port,"var")

                fprintf (stderr, "config.xml: unknown attribute '%s' on %s "
                                 "(ignored)\n"
                        ,attr->name
                        ,"staging-service"
                        );
            }

            config->nssi_server_ids[service_num].nid = (nssi_pid)atoll(nid);
            config->nssi_server_ids[service_num].pid = (nssi_pid)atoll(pid);
            config->nssi_server_ids[service_num].hostname[0]='\0';
            strncpy(config->nssi_server_ids[service_num].hostname, hostname, NSSI_HOSTNAME_LEN);
    #ifdef HAVE_INFINIBAND
            config->nssi_server_ids[service_num].addr = inet_addr(config->nssi_server_ids[service_num]->hostname);
    #endif
            config->nssi_server_ids[service_num].port = (nssi_pid)atoi(port);

//            log_debug(netcdf_config_debug_level, "service (nid=%lld, pid=%llu, hostname=%s, addr=%d, port=%d)",
//                    (long long)id->nid, (long long)id->pid, id->hostname, id->addr, id->port);

//            printf("staging-service: service_num(%d) nid(%lld) pid(%llu) hostname(%s) port(%d)\n",
//                    service_num,
//                    config->nssi_server_ids[service_num].nid,
//                    config->nssi_server_ids[service_num].pid,
//                    config->nssi_server_ids[service_num].hostname,
//                    config->nssi_server_ids[service_num].port);

            service_num++;

        } else {
            if (!strncmp (n->value.element.name, "!--", 3)) // a comment
            {
                continue;
            }
            else
            {
                fprintf (stderr, "config.xml: invalid xml element: '%s'\n"
                        ,n->value.element.name
                        );

                return 0;
            }
        }
    }

    return 1;
}


int parse_nssi_config(const char *config_file, struct adios_nssi_config *config)
{
    FILE * fp = 0;
    mxml_node_t * doc = NULL;
    mxml_node_t * node = NULL;
    mxml_node_t * root = NULL;
    int saw_staging_group = 0;

    int i;

    char * buffer = NULL;
//#if HAVE_MPI
    int buffer_size = 0;
    int rank;
    MPI_Comm_rank (MPI_COMM_WORLD, &rank);
    if (rank == 0)
    {
//#endif
        fp = fopen (config_file, "r");
        if (!fp)
        {
            fprintf (stderr, "missing config file %s\n", config_file);

            return 0;
        }
        struct stat s;
        if (stat (config_file, &s) == 0)
        {
            buffer = malloc (s.st_size + 1);
            buffer [s.st_size] = 0;
        }
        if (buffer)
        {
            size_t bytes_read = fread (buffer, 1, s.st_size, fp);
            if (bytes_read != s.st_size)
            {
                fprintf (stderr, "error reading config file: %s. Expected %d Got %d\n"
                        ,config_file, s.st_size, bytes_read );

                return 0;
            }
        }
        else
        {
            fprintf (stderr, "error allocating %d for reading config.\n"
                    ,s.st_size + 1
                    );

            return 0;
        }
        fclose (fp);
//#if HAVE_MPI
        buffer_size = s.st_size;
        MPI_Bcast (&buffer_size, 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast (buffer, buffer_size, MPI_BYTE, 0, MPI_COMM_WORLD);
    }
    else
    {
        MPI_Bcast (&buffer_size, 1, MPI_INT, 0, MPI_COMM_WORLD);
        buffer = malloc (buffer_size + 1);
        if (!buffer)
        {
            fprintf (stderr, "cannot allocate %d bytes to receive config file\n"
                    ,buffer_size + 1
                    );

            return 0;
        }
        MPI_Bcast (buffer, buffer_size, MPI_BYTE, 0, MPI_COMM_WORLD);
        buffer [buffer_size] = 0;
    }
//#endif

    doc = mxmlLoadString (NULL, buffer, MXML_TEXT_CALLBACK);
    free (buffer);
    buffer = NULL;

    if (!doc)
    {
        fprintf (stderr, "config.xml: unknown error parsing XML "
                         "(probably structural)\n"
                         "Did you remember to start the file with\n"
                         "<?xml version=\"1.0\"?>\n");

        return 0;
    }

    root = doc;

    while (root && root->type != MXML_ELEMENT)
    {
        root = mxmlWalkNext (root, doc, MXML_DESCEND);
    }

    while (!strncmp (root->value.element.name, "!--", 3))
    {
        root = mxmlWalkNext (root, doc, MXML_NO_DESCEND);
        root = mxmlWalkNext (root, doc, MXML_NO_DESCEND);
    }

    if (strcasecmp (root->value.element.name, "nssi-config"))
    {
        if (strncmp (root->value.element.name, "?xml", 4))
        {
            fprintf (stderr, "config.xml: invalid root xml element: %s\n"
                    ,root->value.element.name
                    );

            mxmlRelease (doc);

            return 0;
        }
        else
        {
            while (!strncmp (root->value.element.name, "!--", 3))
            {
                root = mxmlWalkNext (root, doc, MXML_NO_DESCEND);
            }

            root = mxmlWalkNext (root, doc, MXML_DESCEND);  // skip ver num
            root = mxmlWalkNext (root, doc, MXML_NO_DESCEND);  // get next
            while (!strncmp (root->value.element.name, "!--", 3))
            {
                root = mxmlWalkNext (root, doc, MXML_NO_DESCEND);
                root = mxmlWalkNext (root, doc, MXML_NO_DESCEND);
            }
        }
    }
    else
    {
        printf("it is nssi-config\n");
    }


    for (node = mxmlWalkNext (root, doc, MXML_DESCEND_FIRST)
        ;node
        ;node = mxmlWalkNext (node, root, MXML_NO_DESCEND)
        )
    {
        if (node->type != MXML_ELEMENT)
        {
            continue;
        }

        if (!strcasecmp (node->value.element.name, "staging-group"))
        {
            const char *write_type = 0;

            for (i = 0; i < node->value.element.num_attrs; i++)
            {
                mxml_attr_t *attr = &node->value.element.attrs[i];

                GET_ATTR("write-type",attr,write_type,"var")

                fprintf (stderr, "config.xml: unknown attribute '%s' on %s "
                                 "(ignored)\n"
                        ,attr->name
                        ,"staging-service"
                        );
            }

            if (!strcmp(write_type, "WRITE_DIRECT")) {
                config->write_type = WRITE_DIRECT;
//                log_debug(LOG_ALL, "using %s", write_type);
            } else if (!strcmp(write_type, "WRITE_AGGREGATE_INDEPENDENT")) {
                config->write_type = WRITE_AGGREGATE_INDEPENDENT;
//                log_debug(LOG_ALL, "using %s", write_type);
            } else if (!strcmp(write_type, "WRITE_AGGREGATE_COLLECTIVE")) {
                config->write_type = WRITE_AGGREGATE_COLLECTIVE;
//                log_debug(LOG_ALL, "using %s", write_type);
            } else if (!strcmp(write_type, "WRITE_CACHING_INDEPENDENT")) {
                config->write_type = WRITE_CACHING_INDEPENDENT;
//                log_debug(LOG_ALL, "using %s", write_type);
            } else if (!strcmp(write_type, "WRITE_CACHING_COLLECTIVE")) {
                config->write_type = WRITE_CACHING_COLLECTIVE;
//                log_debug(LOG_ALL, "using %s", write_type);
            }

            if (!parseGroup(node, config))
                break;
            saw_staging_group = 1;
        }
        else
        {
            if (!strncmp (node->value.element.name, "!--", 3))
            {
                continue;
            }
            else
            {
                fprintf (stderr, "config.xml: invalid element: %s\n"
                        ,node->value.element.name
                );

                break;
            }
        }
    }

    mxmlRelease (doc);

    if (!saw_staging_group)
    {
        fprintf (stderr, "config.xml: must define at least 1 staging-group in "
                         "config.xml\n"
                );

        return 0;
    }

//    for (i=0;i<config->num_servers;i++) {
//        printf("staging-service: service_num(%d) nid(%lld) pid(%llu) hostname(%s) port(%d)\n",
//                i,
//                config->nssi_server_ids[i].nid,
//                config->nssi_server_ids[i].pid,
//                config->nssi_server_ids[i].hostname,
//                config->nssi_server_ids[i].port);
//    }

    return 1;
}

void free_nssi_config(struct adios_nssi_config *config)
{
    if (config->nssi_server_ids) {
        free(config->nssi_server_ids);
        config->nssi_server_ids=NULL;
    }
}
