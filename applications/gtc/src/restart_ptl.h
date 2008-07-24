#include <thin_portal.h>

typedef struct restart_data_
{
    char * filename;
    int length;
    char * buffer;
} restart_data, *restart_data_p;

static  IOField restart_data_field_list[] =
{
	{"filename", "string", sizeof (char *), IOOffset (restart_data_p, filename)},
	{"length", "integer", sizeof(int), IOOffset(restart_data_p, length)},
	{"buffer","char[length]",sizeof(char),IOOffset(restart_data_p,buffer)},
	{NULL, NULL, 0, 0}
}; 
static IOFormatRec restart_data_format_list[] = 
{
	{"restart_format", restart_data_field_list},
	{NULL, NULL}
};
