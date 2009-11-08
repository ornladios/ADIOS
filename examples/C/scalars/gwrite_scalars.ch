adios_groupsize = 1 \
                + 2 \
                + 4 \
                + 8 \
                + 1 \
                + 2 \
                + 4 \
                + 8 \
                + 4 \
                + 8 \
                + strlen(v11) \
                + 8 \
                + 16;
adios_group_size (adios_handle, adios_groupsize, &adios_totalsize);
adios_write (adios_handle, "var_byte", &v1);
adios_write (adios_handle, "var_short", &v2);
adios_write (adios_handle, "var_int", &v3);
adios_write (adios_handle, "var_long", &v4);
adios_write (adios_handle, "var_ubyte", &v5);
adios_write (adios_handle, "var_ushort", &v6);
adios_write (adios_handle, "var_uint", &v7);
adios_write (adios_handle, "var_ulong", &v8);
adios_write (adios_handle, "var_real", &v9);
adios_write (adios_handle, "var_double", &v10);
adios_write (adios_handle, "var_string", v11);
adios_write (adios_handle, "var_complex", &v12);
adios_write (adios_handle, "var_double_complex", &v13);
