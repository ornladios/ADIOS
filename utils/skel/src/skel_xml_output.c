#include "skel_xml_output.h"
#include <stdio.h>

int FC_FUNC_(skel_write_coarse_xml_data_f, SKEL_WRITE_COARSE_XML_DATA_F) ()
{
    return common_skel_write_coarse_xml_data ();
}


int skel_write_coarse_xml_data ()
{
    return common_skel_write_coarse_xml_data ();
}


int common_skel_write_coarse_xml_data ()
{
    fprintf (stdout, "In common_skel_write_coarse_xml_data\n");

    return 0;
}
