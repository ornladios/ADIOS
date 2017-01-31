%module adioslib
%include "stdint.i"
%include "typemaps.i"
%include "arrays_java.i"
%include "m_adios_namelist.i"
%include "m_adios_dims.i"
%include "m_adios_fh.i"

%apply int64_t *OUTPUT { int64_t *fd };
%apply char **STRING_ARRAY { char **var_namelist }
%apply char **STRING_ARRAY { char **attr_namelist }

/*
%typemap(out) uint64_t *dims {
      *(struct Bar **)&$result = $1;
      printf("Here... %p\n", $1);
      *(uint64_t **)&jresult = result;
}
*/

%{
#include "adios.h"
#include "adios_read_v2.h"
%}

// pakcage-private
//%pragma(java) jniclassclassmodifiers="class"

// compiler constants
%include "enums.swg"
/*
%typemap(javain) enum SWIGTYPE "$javainput.ordinal()"
%typemap(javaout) enum SWIGTYPE {
    return $javaclassname.class.getEnumConstants()[$jnicall];
  }
%typemap(javabody) enum SWIGTYPE ""
*/

%typedef struct _ADIOS_FILE ADIOS_FILE;
%typedef struct _ADIOS_VARINFO ADIOS_VARINFO;
%rename(ADIOS_FILE) _ADIOS_FILE;
%rename(ADIOS_VARINFO) _ADIOS_VARINFO;

%ignore adios_reset_dimension_order;
%ignore adios_stat_cor;
%ignore adios_stat_cov;

%include "mpidummy.h"
%include "adios_types.h"
%include "adios_read_v2.h"
%include "adios.h"

%include "m_adios_read.i"
%include "m_adios_write.i"
