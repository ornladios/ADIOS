/*
 * ADIOS is freely available under the terms of the BSD license described
 * in the COPYING file in the top level directory of this source distribution.
 *
 * Copyright (c) 2008 - 2009.  UT-BATTELLE, LLC. All rights reserved.
 */

#ifndef ADIOS_NSSI_CONFIG_H
#define ADIOS_NSSI_CONFIG_H

/**
 * @brief A structure to represent the configuration of
 * NSSI core services.
 */
struct adios_nssi_config {

    /** @brief Number of available storage servers */
    int num_servers;

    /** @brief storage service IDs */
    nssi_remote_pid *nssi_server_ids;

    /** @brief The type of write operation the client wishes the server to perform */
    enum write_type write_type;
};


int parse_nssi_config(const char *config_file, struct adios_nssi_config *config);
void free_nssi_config(struct adios_nssi_config *config);

#endif
