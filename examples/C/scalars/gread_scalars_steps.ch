s = adios_selection_writeblock (rank);
adios_perform_reads (fp, 1);
adios_selection_delete (s);
