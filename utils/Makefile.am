SUBDIRS= gpp bpdump bp2ascii bpsplit bpls skeldump adios_list_methods bpmeta bprecover


if BUILD_WRITE
SUBDIRS += adios_lint 
if BUILD_MPI
    SUBDIRS += bp2bp bpdiff
if HAVE_FASTBIT
        SUBDIRS +=fastbit
endif HAVE_FASTBIT
endif BUILD_MPI
endif BUILD_WRITE

if HAVE_HDF5
SUBDIRS += bp2h5
endif

if HAVE_NETCDF
SUBDIRS +=bp2ncd
endif

SUBDIRS +=skel
