
#include "adios_read.h"
#include "adios_read_transformed.h"

void adios_free_transinfo(ADIOS_VARINFO *vi, ADIOS_TRANSINFO *ti) {
    common_read_free_transinfo(ti);
}
