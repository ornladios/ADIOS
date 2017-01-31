
%define WRITEVAR_ARR(CTYPE)
%apply CTYPE[] {CTYPE *arr};
%inline %{
int adios_write_##CTYPE##_arr (int64_t fd_p, const char * name, CTYPE *arr) {
  return adios_write (fd_p, name, arr);
}
%}
%enddef

%define WRITEVAR(CTYPE)
%inline %{
int adios_write_##CTYPE (int64_t fd_p, const char * name, CTYPE value) {
  return adios_write (fd_p, name, &value);
}
%}
%enddef

WRITEVAR_ARR(double)
WRITEVAR_ARR(int)

WRITEVAR(double)
WRITEVAR(int)
