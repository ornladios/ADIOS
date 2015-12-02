#
# Configuration values from configure script
#

# Install directory
ADIOS_DIR="@prefix@"

# Flags to build code using ADIOS write API (and read API)
ADIOS_INC="-I${ADIOS_DIR}/include @ADIOSLIB_CPPFLAGS@ @ADIOSLIB_CFLAGS@"
ADIOS_CLIB="${ADIOS_DIR}/lib/libadios.a @ADIOSLIB_LDADD_M@ @LDFLAGS@ @LIBS@"
ADIOS_FLIB="${ADIOS_DIR}/lib/libadiosf.a @ADIOSLIB_LDADD_M@ @LDFLAGS@ @LIBS@"
ADIOS_V1_FLIB="${ADIOS_DIR}/lib/libadiosf_v1.a @ADIOSLIB_LDADD_M@ @LIBS@"

# Flags to build code using ADIOS read API only
ADIOSREAD_INC="-I${ADIOS_DIR}/include @ADIOSREADLIB_CPPFLAGS@ @ADIOSREADLIB_CFLAGS@"
ADIOSREAD_V1_INC="-I${ADIOS_DIR}/include @MACRODEFFLAG@ADIOS_USE_READ_API_1 @ADIOSREADLIB_CPPFLAGS@ @ADIOSREADLIB_CFLAGS@"
ADIOSREAD_CLIB="${ADIOS_DIR}/lib/libadiosread.a @ADIOSREADLIB_LDADD_M@"
ADIOSREAD_FLIB="${ADIOS_DIR}/lib/libadiosreadf.a @ADIOSREADLIB_LDADD_M@"
ADIOSREAD_V1_FLIB="${ADIOS_DIR}/lib/libadiosreadf_v1.a @ADIOSREADLIB_LDADD_M@"

# Flags to build code using ADIOS read API only in a sequential code (no MPI)
ADIOSREAD_SEQ_INC="-I${ADIOS_DIR}/include @ADIOSREADLIB_SEQ_CPPFLAGS@ @ADIOSREADLIB_SEQ_CFLAGS@"
ADIOSREAD_SEQ_V1_INC="-I${ADIOS_DIR}/include @MACRODEFFLAG@ADIOS_USE_READ_API_1 @ADIOSREADLIB_SEQ_CPPFLAGS@ @ADIOSREADLIB_SEQ_CFLAGS@"
ADIOSREAD_SEQ_CLIB="${ADIOS_DIR}/lib/libadiosread_nompi.a @ADIOSREADLIB_SEQ_LDADD_M@"
ADIOSREAD_SEQ_FLIB="${ADIOS_DIR}/lib/libadiosreadf_nompi.a @ADIOSREADLIB_SEQ_LDADD_M@"
ADIOSREAD_SEQ_V1_FLIB="${ADIOS_DIR}/lib/libadiosreadf_nompi_v1.a @ADIOSREADLIB_SEQ_LDADD_M@"

#Flags to build code using ADIOS write API in a sequential code (no MPI)
ADIOS_SEQ_INC="-I${ADIOS_DIR}/include @ADIOSLIB_SEQ_CPPFLAGS@ @ADIOSLIB_SEQ_CFLAGS@"
ADIOS_SEQ_CLIB="${ADIOS_DIR}/lib/libadios_nompi.a @ADIOSLIB_SEQ_LDADD_M@ @LIBS@"
ADIOS_SEQ_FLIB="${ADIOS_DIR}/lib/libadiosf_nompi.a @ADIOSLIB_SEQ_LDADD_M@ @LIBS@"
ADIOS_SEQ_V1_FLIB="${ADIOS_DIR}/lib/libadiosf_nompi_v1.a @ADIOSLIB_SEQ_LDADD_M@ @LIBS@"

#The following flags are not used. It is only for internal utilities of ADIOS, using libadios_internal_nompi.a
ADIOS_INT_INC="-I${ADIOS_DIR}/include @ADIOSLIB_INT_CPPFLAGS@ @ADIOSLIB_INT_CFLAGS@"
ADIOS_INT_CLIB="${ADIOS_DIR}/lib/libadios_internal_nompi.a @ADIOSLIB_INT_LDADD_M@ @LIBS@"

VERSIONSTRING="@VERSION@"

