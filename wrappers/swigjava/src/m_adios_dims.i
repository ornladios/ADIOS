%ignore _ADIOS_VARINFO::dims;
%extend _ADIOS_VARINFO {
  uint64_t getDims(int i) {
    return $self->dims[i];
  }
}
