#include "config.h"
#ifdef BUILD_WITH_CMAKE
  #include "FC.h"
#endif

#include <stdio.h>
#include <stdlib.h>


/* defined in skel_xml_output.c */
int common_skel_write_coarse_xml_data (double open_time, double write_time, double close_time, double total_time);


int FC_FUNC_(skel_write_coarse_xml_data_f, SKEL_WRITE_COARSE_XML_DATA_F) (double *open_time, double *write_time, double *close_time, double *total_time)
{
    return common_skel_write_coarse_xml_data (*open_time, *write_time, *close_time, *total_time);
}


