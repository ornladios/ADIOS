SUBDIRS= gpp bpdump bp2ascii bpsplit bpls bpdiff skeldump list_methods bpmeta
if BUILD_WRITE
SUBDIRS += adios_lint bp2bp 
endif
if HAVE_HDF5
SUBDIRS += bp2h5
endif
if HAVE_NETCDF
SUBDIRS +=bp2ncd
endif
if HAVE_FASTBIT
SUBDIRS +=query
endif
SUBDIRS +=skel
