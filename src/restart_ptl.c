#include "restart_ptl.h"
#include <stdlib.h>
#include <stdio.h>

#define MAXSEND 100

static IOhandle *h = NULL;
static IOFormat ioformat;
static restart_data *rdata = NULL;

void sinit_ ()
{
	h = InitIOFromFile ("param");
	if (h == NULL)
	{
		printf ("something went wrong\n");
		return;
	}

	ioformat = register_data (h, restart_data_format_list);
	if (ioformat == NULL)
	{
		printf ("register data failed\n");
		free (h);
		h = NULL;
		return;
	}
	rdata = NULL;
}

void swait_ ()
{
	if (rdata == NULL)
	{
		return;
	}
	send_end (h);
	if (rdata)
		free (rdata);
}

void ssend_ (char * filename, int * length, char * buffer)
{
	if (h == NULL)
		return;
        if (!rdata)
		rdata = (restart_data *) malloc (sizeof(restart_data));

        rdata->filename = filename;
	rdata->length = *length;
	rdata->buffer = buffer;

	start_send (h, rdata, sizeof(restart_data), ioformat);
}

void  terminate_ ()
{
	outputTimingInfo ("benchmarkresults");
}
