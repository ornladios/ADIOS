%{
#define JCALL0(func, jenv) (*jenv)->func(jenv)
#define JCALL1(func, jenv, ar1) (*jenv)->func(jenv, ar1)
#define JCALL2(func, jenv, ar1, ar2) (*jenv)->func(jenv, ar1, ar2)
#define JCALL3(func, jenv, ar1, ar2, ar3) (*jenv)->func(jenv, ar1, ar2, ar3)
#define JCALL4(func, jenv, ar1, ar2, ar3, ar4) (*jenv)->func(jenv, ar1, ar2, ar3, ar4)
%}

%define READVAR(CTYPE, CAP, JTYPE)
%typemap(jni) (void* CTYPE##_arr) "JTYPE"
%typemap(jstype) (void* CTYPE##_arr) "CTYPE[]"
%typemap(jtype) (void* CTYPE##_arr) "CTYPE[]"
%typemap(javain) (void* CTYPE##_arr) "$javainput"
%typemap(in) (void* CTYPE##_arr) {
  $1 = JCALL2(Get##CAP##ArrayElements, jenv, $input, 0);
}
%typemap(freearg) (void* CTYPE##_arr) {
  JCALL3(Release##CAP##ArrayElements, jenv, $input, $1, 0);
}

%inline %{
int readvar_##CTYPE(ADIOS_FILE *fh, int varid, const int64_t start[], const int64_t count[], void* CTYPE##_arr) {
  int result, i;
  ADIOS_VARINFO * v;
  ADIOS_SELECTION * sel;

  v = adios_inq_var_byid (fh, varid);

  printf("=== readvar_val ===\n");
  printf("fh: %p\n", fh);
  printf("varid: %d\n", varid);
  printf("start: ");
  for (i=0; i<v->ndim; i++)
    printf("%lld ", start[i]);
  printf("\n");
  printf("count: ");
  for (i=0; i<v->ndim; i++)
    printf("%lld ", count[i]);
  printf("\n");

  sel = adios_selection_boundingbox (v->ndim, (const uint64_t *)start, (const uint64_t *)count);
  result = adios_schedule_read_byid (fh, sel, varid, 0, 1, CTYPE##_arr);
  adios_perform_reads (fh, 1);

  adios_selection_delete (sel);
  adios_free_varinfo (v);
  return result;
}
%}
%enddef

READVAR(double, Double, jdoubleArray)
READVAR(int, Int, jintArray)
