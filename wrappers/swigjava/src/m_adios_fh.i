%ignore _ADIOS_FILE::fh;
%extend _ADIOS_FILE {
  int64_t getFh() {
    return (int64_t) $self->fh;
  }
}
