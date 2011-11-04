INCLUDES = $(all_includes)

AUTOMAKE_OPTIONS = no-dependencies

lib_LIBRARIES = libardma_client.a libardma_server.a


if HAVE_PORTALS
libardma_client_a_SOURCES = portals_client.c portals_common.c
libardma_client_a_CPPFLAGS = $(PORTALS_CPPFLAGS)
libardma_client_a_CFLAGS = $(PORTALS_CFLAGS)
libardma_server_a_SOURCES = portals_server.c portals_common.c
libardma_server_a_CPPFLAGS = $(PORTALS_CPPFLAGS)
libardma_server_a_CFLAGS = $(PORTALS_CFLAGS)
endif

if HAVE_RDMACM
libardma_client_a_SOURCES = ib_client.c 
libardma_client_a_CPPFLAGS = $(RDMACM_CPPFLAGS) $(INFINIBAND_CPPFLAGS)
libardma_client_a_CFLAGS = $(RDMACM_CFLAGS) $(INFINIBAND_CFLAGS)
libardma_server_a_SOURCES = ib_server.c 
libardma_server_a_CPPFLAGS = $(RDMACM_CPPFLAGS) $(INFINIBAND_CPPFLAGS)
libardma_server_a_CFLAGS = $(RDMA_CFLAGS) $(INFINIBAND_CFLAGS)
endif

CC=${MPICC}


