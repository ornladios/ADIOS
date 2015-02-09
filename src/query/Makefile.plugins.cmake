if(HAVE_FASTBIT)
set(query_method_SOURCES ${query_method_SOURCES} query/query_fastbit.c)
set(query_method_SOURCES ${query_method_SOURCES} query/fastbit_adios.c)
endif() # HAVE_FASTBIT

if(HAVE_ALACRITY)
set(query_method_SOURCES ${query_method_SOURCES} query/query_alac.c)
endif() # HAVE_ALACRITY 
