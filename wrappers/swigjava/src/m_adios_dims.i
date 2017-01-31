%ignore _ADIOS_VARINFO::dims;
%extend _ADIOS_VARINFO {
  int64_t getDims(int i) {
    return (int64_t) $self->dims[i];
  }
}
