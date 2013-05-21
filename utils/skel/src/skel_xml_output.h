#ifndef SKEL_XML_OUTPUT_H
#define SKEL_XML_OUTPUT_H

#include "config.h"

int FC_FUNC_(skel_write_coarse_xml_data_f, SKEL_WRITE_COARSE_XML_DATA_F) (double *open_time, double *write_time, double *close_time, double *total_time);

int skel_write_coarse_xml_data (double open_time, double write_time, double close_time, double total_time);

#endif
