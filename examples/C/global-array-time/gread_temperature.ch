adios_group_size (adios_handle, adios_groupsize, &adios_totalsize, &comm);
adios_buf_size = 4;
adios_read (adios_handle, "NX", &NX, adios_buf_size);
adios_buf_size = 4;
adios_read (adios_handle, "size", &size, adios_buf_size);
adios_buf_size = 4;
adios_read (adios_handle, "rank", &rank, adios_buf_size);
adios_buf_size = 8;
adios_read (adios_handle, "temperature", temperature, adios_buf_size);
