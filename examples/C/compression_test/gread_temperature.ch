s = adios_selection_writeblock (rank);
adios_schedule_read (fp, s, "temperature", 1, 1, t);
adios_perform_reads (fp, 1);
adios_selection_delete (s);
