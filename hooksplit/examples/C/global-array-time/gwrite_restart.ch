adios_groupsize = 4 \
                + 4 \
                + 4 \
                + 4 \
                + 8 * (1) * (NX) \
                + 8 * (1) * (NY);
adios_group_size (adios_handle, adios_groupsize, &adios_totalsize);
adios_write (adios_handle, "NX", &NX);
adios_write (adios_handle, "NY", &NY);
adios_write (adios_handle, "size", &size);
adios_write (adios_handle, "rank", &rank);
adios_write (adios_handle, "temperature", t);
adios_write (adios_handle, "pressure", p);
