adios_groupsize = 4 \
                + 4 \
                + 4 \
                + 8 * (2) * (NX) \
                + 8 * (1) * (NX);
adios_group_size (adios_handle, adios_groupsize, &adios_totalsize);
adios_write (adios_handle, "NX", &NX);
adios_write (adios_handle, "size", &size);
adios_write (adios_handle, "rank", &rank);
adios_write (adios_handle, "complex", c);
adios_write (adios_handle, "temperature", t);
