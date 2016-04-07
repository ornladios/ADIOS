SUBDIRS= gpp bpdump bp2ascii bpsplit bpls bpdiff skeldump adios_list_methods bpmeta bprecover

if BUILD_WRITE
SUBDIRS += adios_lint bp2bp 
if HAVE_FASTBIT
    SUBDIRS +=fastbit
endif
endif

if HAVE_HDF5
SUBDIRS += bp2h5
endif

if HAVE_NETCDF
SUBDIRS +=bp2ncd
endif

SUBDIRS +=skel
