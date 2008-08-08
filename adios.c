#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <unistd.h>
#include <stdint.h>

// xml parser
#include <mxml.h>

// Chen's encoder
#include "bw-utils.h"
#include "br-utils.h"

#include "adios.h"
#include "adios_transport_hooks.h"
#include "adios_bp_v1.h"
#include "adios_internals.h"

#define STR_LEN 1000

extern struct adios_transport_struct * adios_transports;

///////////////////////////////////////////////////////////////////////////////
static int common_adios_init (const char * config)
{
    // parse the config file
    return adios_parse_config (config);
}

int adios_init (const char * config)
{
    return common_adios_init (config);
}

void adios_init_ (const char * config, int * err, int config_size)
{
    char buf1 [STR_LEN] = "";

    adios_extract_string (buf1, config, config_size);

    *err = common_adios_init (buf1);
}

///////////////////////////////////////////////////////////////////////////////
static int common_adios_finalize (int mype)
{
    struct adios_method_list_struct * m;

    for (m = adios_get_methods (); m; m = m->next)
    {
        if (   m->method->m != ADIOS_METHOD_UNKNOWN
            && m->method->m != ADIOS_METHOD_NULL
            && adios_transports [m->method->m].adios_finalize_fn
           )
        {
            adios_transports [m->method->m].adios_finalize_fn (mype, m->method);
        }
    }

    return 0;
}

int adios_finalize (int mype)
{
    return common_adios_finalize (mype);
}

void adios_finalize_ (int * mype, int * err)
{
    *err = common_adios_finalize (*mype);
}

///////////////////////////////////////////////////////////////////////////////
static int common_adios_allocate_buffer ()
{
    return adios_set_buffer_size ();
}

int adios_allocate_buffer ()
{
    return common_adios_allocate_buffer ();
}

void adios_allocate_buffer_ (int * err)
{
    *err = common_adios_allocate_buffer ();
}

///////////////////////////////////////////////////////////////////////////////
static int common_adios_open (long long * fd, const char * group_name
                             ,const char * name, const char * file_mode
                             )
{
    long long group_id = 0;
    struct adios_file_struct * fd_p = (struct adios_file_struct *)
                                  malloc (sizeof (struct adios_file_struct));
    struct adios_group_struct * g = 0;
    struct adios_method_list_struct * methods = 0;
    enum ADIOS_METHOD_MODE mode;

    adios_common_get_group (&group_id, group_name);
    g = (struct adios_group_struct *) group_id;
    methods = g->methods;

    if (!strcasecmp (file_mode, "r"))
        mode = adios_mode_read;
    else
        if (!strcasecmp (file_mode, "w"))
            mode = adios_mode_write;
        else
            if (!strcasecmp (file_mode, "a"))
                mode = adios_mode_append;
            else
                if (!strcasecmp (file_mode, "u"))
                    mode = adios_mode_update;
                else
                {
                    fprintf (stderr, "adios_open: unknown file mode: %s\n"
                            ,file_mode
                            );

                    *fd = 0;

                    return 1;
                }

    fd_p->name = strdup (name);
    fd_p->group = g;
    fd_p->mode = mode;
    fd_p->data_size = 0;
    fd_p->buffer = 0;
    fd_p->offset = 0;
    fd_p->bytes_written = 0;
    fd_p->buffer_size = 0;
    fd_p->vars_start = 0;
    fd_p->vars_written = 0;
    fd_p->write_size_bytes = 0;
    fd_p->base_offset = 0;
    fd_p->pg_start_in_file = 0;

    while (methods)
    {
        if (   methods->method->m != ADIOS_METHOD_UNKNOWN
            && methods->method->m != ADIOS_METHOD_NULL
            && adios_transports [methods->method->m].adios_open_fn
           )
        {
            adios_transports [methods->method->m].adios_open_fn
                                                 (fd_p, methods->method);
        }

        methods = methods->next;
    }

    *fd = (long long) fd_p;

    return 0;
}

int adios_open (long long * fd, const char * group_name, const char * name
               ,const char * mode
               )
{
    return common_adios_open (fd, group_name, name, mode);
}

void adios_open_ (long long * fd, const char * group_name, const char * name
                 ,const char * mode, int * err
                 ,int group_name_size, int name_size, int mode_size
                 )
{
    char buf1 [STR_LEN] = "";
    char buf2 [STR_LEN] = "";
    char buf3 [STR_LEN] = "";

    adios_extract_string (buf1, group_name, group_name_size);
    adios_extract_string (buf2, name, name_size);
    adios_extract_string (buf3, mode, mode_size);

    *err = common_adios_open (fd, buf1, buf2, buf3);
}

///////////////////////////////////////////////////////////////////////////////
static int common_adios_group_size (long long fd_p
                                   ,uint64_t data_size
                                   ,uint64_t * total_size
                                   ,void * comm
                                   )
{
    struct adios_file_struct * fd = (struct adios_file_struct *) fd_p;
    struct adios_method_list_struct * m = fd->group->methods;

    fd->write_size_bytes = data_size;

    uint64_t overhead = adios_calc_overhead_v1 (fd);

    *total_size = data_size + overhead;

    // try to reserve a buffer using the adios_method_buffer_alloc
    // if it does not give the correct amount, overflow.  Make sure
    // the amount given is big enough for the first part of the file.

    fd->write_size_bytes += overhead;

    uint64_t allocated = adios_method_buffer_alloc (fd->write_size_bytes);
    if (allocated != fd->write_size_bytes)
    {
        fd->shared_buffer = adios_flag_no;
    }
    else
    {
        fd->shared_buffer = adios_flag_yes;
    }

    // call each transport method to coordinate the write and handle
    // if an overflow is detected.
    // now tell each transport attached that it is being written
    while (m)
    {
        enum ADIOS_FLAG should_buffer = adios_flag_yes;
        if (   m->method->m != ADIOS_METHOD_UNKNOWN
            && m->method->m != ADIOS_METHOD_NULL
            && adios_transports [m->method->m].adios_should_buffer_fn
           )
        {
            should_buffer = adios_transports [m->method->m].
                                            adios_should_buffer_fn (fd
                                                                   ,m->method
                                                                   ,comm
                                                                   );
        }

        if (should_buffer == adios_flag_no)     // can't write directly since
            fd->shared_buffer = adios_flag_no;  // some might want to share

        m = m->next;
    }

    if (fd->shared_buffer == adios_flag_no)
    {
        adios_method_buffer_free (allocated);
        fd->buffer = 0;
        fd->offset = 0;
        fd->bytes_written = 0;
    }
    else
    {
        fd->buffer = malloc (fd->write_size_bytes);
        fd->offset = 0;
        fd->bytes_written = 0;
        if (!fd->buffer)
        {
            fprintf (stderr, "Cannot allocate %lld bytes for buffered "
                             "output.\n"
                    );

            return 1;
        }
        else
        {
            // write the process group header
            adios_write_process_group_header_v1 (fd, *total_size);

            // setup for writing vars
            adios_write_open_vars_v1 (fd);
        }
    }

    // each var will be added to the buffer by the adios_write calls
    // attributes will be added by adios_close

    return 0;
}

int adios_group_size (long long fd_p, uint64_t data_size
                     ,uint64_t * total_size, void * comm
                     )
{
    return common_adios_group_size (fd_p, data_size, total_size, comm);
}

void adios_group_size_ (long long * fd_p, int64_t * data_size
                       ,int64_t * total_size, void * comm, int * err
                       )
{
    *err = common_adios_group_size (*fd_p, (uint64_t) *data_size
                                   ,(uint64_t *) total_size, comm
                                   );
}

///////////////////////////////////////////////////////////////////////////////
static int common_adios_write (long long fd_p, const char * name, void * var)
{
    struct adios_file_struct * fd = (struct adios_file_struct *) fd_p;
    struct adios_var_struct * v = fd->group->vars;
    struct adios_method_list_struct * m = fd->group->methods;

    v = adios_find_var_by_name (v, name, fd->group->all_unique_var_names);

    if (!v)
    {
        fprintf (stderr, "Bad var name (ignored): '%s'\n", name);

        return 1;
    }

    if (fd->mode == adios_mode_read)
    {
        if (   strcasecmp (name, fd->group->group_comm)
            && v->is_dim != adios_flag_yes
           )
        {
            fprintf (stderr, "write attempted on %s in %s.  This was opened for"
                             " read\n"
                    ,name , fd->name
                    );

            return 1;
        }
    }

    if (!var)
    {
        fprintf (stderr, "Invalid data: %s\n", name);

        return 1;
    }

    v->data = var;
    if (fd->shared_buffer == adios_flag_yes)
    {
        // var payload sent for sizing information
        adios_write_var_header_v1 (fd, v);

        // generate characteristics (like min and max)
        adios_generate_var_characteristics_v1 (fd, v);

        // write these characteristics
        adios_write_var_characteristics_v1 (fd, v);

        // write payload
        adios_write_var_payload_v1 (fd, v);
    }

    // now tell each transport attached that it is being written
    while (m)
    {
        if (   m->method->m != ADIOS_METHOD_UNKNOWN
            && m->method->m != ADIOS_METHOD_NULL
            && adios_transports [m->method->m].adios_write_fn
           )
        {
            adios_transports [m->method->m].adios_write_fn
                                   (fd, v, var, m->method);
        }

        m = m->next;
    }

    if (v->is_dim == adios_flag_no)
    {
        v->data = 0;
    }
    else
    {
        int element_size = adios_get_type_size (v->type, var);
        v->data = malloc (element_size);

        if (!v->data)
        {
            fprintf (stderr, "cannot allocate %d bytes to copy scalar %s\n"
                    ,element_size
                    ,v->name
                    );

            return 0;
        }

        switch (v->type)
        {
            case adios_byte:
            case adios_short:
            case adios_integer:
            case adios_long:
            case adios_unsigned_byte:
            case adios_unsigned_short:
            case adios_unsigned_integer:
            case adios_unsigned_long:
            case adios_real:
            case adios_double:
            case adios_long_double:
            case adios_string:
            case adios_complex:
            case adios_double_complex:
                memcpy ((char *) v->data, (char *) var, element_size);
                break;

            default:
                v->data = 0;
                break;
        }
    }

    return 0;
}

int adios_write (long long fd_p, const char * name, void * var)
{
    return common_adios_write (fd_p, name, var);
}

void adios_write_ (long long * fd_p, const char * name, void * var, int * err
                  ,int name_size
                  )
{
    char buf1 [STR_LEN] = "";

    adios_extract_string (buf1, name, name_size);

    *err = common_adios_write (*fd_p, buf1, var);
}


///////////////////////////////////////////////////////////////////////////////
static int common_adios_get_write_buffer (long long fd_p, const char * name
                                         ,uint64_t * size
                                         ,void ** buffer
                                         )
{
    struct adios_file_struct * fd = (struct adios_file_struct *) fd_p;
    struct adios_var_struct * v = fd->group->vars;
    struct adios_method_list_struct * m = fd->group->methods;

    v = adios_find_var_by_name (v, name, fd->group->all_unique_var_names);

    if (!v)
    {
        fprintf (stderr
                ,"Bad var name (ignored): '%s' (%c%c%c)\n"
                ,name, name[0], name[1], name[2]
                );

        return 1;
    }

    if (fd->mode == adios_mode_read)
    {
        fprintf (stderr, "write attempted on %s in %s.  This was opened for"
                         " read\n"
                ,name , fd->name
                );

        return 1;
    }

    // since we are only getting one buffer, get it from the first
    // transport method that can provide it.
    while (m)
    {
        if (   m->method->m != ADIOS_METHOD_UNKNOWN
            && m->method->m != ADIOS_METHOD_NULL
            && adios_transports [m->method->m].adios_write_fn
           )
        {
            adios_transports [m->method->m].adios_get_write_buffer_fn
                                (fd, v, size, buffer, m->method);
            m = 0;
        }
        else
            m = m->next;
    }

    return 0;
}

int adios_get_write_buffer (long long fd_p, const char * name
                           ,uint64_t * size
                           ,void ** buffer
                           )
{
    return common_adios_get_write_buffer (fd_p, name, size, buffer);
}

void adios_get_write_buffer_ (long long * fd_p, const char * name
                             ,uint64_t * size
                             ,void ** buffer, int * err, int name_size
                             )
{
    char buf1 [STR_LEN] = "";

    adios_extract_string (buf1, name, name_size);

    *err = common_adios_get_write_buffer (*fd_p, buf1, size, buffer);
}

///////////////////////////////////////////////////////////////////////////////
static int common_adios_read (long long fd_p, const char * name, void * buffer
                             ,uint64_t buffer_size
                             )
{
    struct adios_file_struct * fd = (struct adios_file_struct *) fd_p;
    struct adios_var_struct * v;

    if (!(fd->mode == adios_mode_read))
    {
        fprintf (stderr, "read attempted on %s which was opened for write\n"
                ,fd->name
                );

        return 1;
    }

    v = adios_find_var_by_name (fd->group->vars, name
                               ,fd->group->all_unique_var_names
                               );
    if (v)
    {
        struct adios_method_list_struct * m = fd->group->methods;

        // since can only read from one place into the buffer,
        // read from the first transport method that can
        while (m)
        {
            if (   m->method->m != ADIOS_METHOD_UNKNOWN
                && m->method->m != ADIOS_METHOD_NULL
                && adios_transports [m->method->m].adios_read_fn
               )
            {
                adios_transports [m->method->m].adios_read_fn
                                     (fd, v, buffer, m->method);
                m = 0;
            }
            else
                m = m->next;
	}
    }
    else
    {
        fprintf (stderr, "var %s in file %s not found on read\n"
                ,name, fd->name
                );

        return 1;
    }

    return 0;
}

int adios_read (long long fd_p, const char * name, void * buffer
               ,uint64_t buffer_size
               )
{
    return common_adios_read (fd_p, name, buffer, buffer_size);
}

void adios_read_ (long long * fd_p, const char * name, void * buffer
                 ,long long * buffer_size, int * err, int name_size
                 )
{
    char buf1 [STR_LEN] = "";

    adios_extract_string (buf1, name, name_size);

    *err = common_adios_read (*fd_p, buf1, buffer, *buffer_size);
}

///////////////////////////////////////////////////////////////////////////////
static int common_adios_set_path (long long fd_p, const char * path)
{
    struct adios_file_struct * fd = (struct adios_file_struct *) fd_p;
    struct adios_group_struct * t = fd->group;
    struct adios_var_struct * v = t->vars;
    struct adios_attribute_struct * a = t->attributes;

    while (v)
    {
        if (v->path)
        {
            free (v->path);
        }

        v->path = strdup (path);

        v = v->next;
    }

    while (a)
    {
        if (a->path)
        {
            free (a->path);
        }

        a->path = strdup (path);

        a = a->next;
    }

    return 0;
}

int adios_set_path (long long fd_p, const char * path)
{
    return common_adios_set_path (fd_p, path);
}

void adios_set_path_ (long long * fd_p, const char * path, int * err
                     ,int path_size
                     )
{
    char buf1 [STR_LEN] = "";

    adios_extract_string (buf1, path, path_size);

    *err = common_adios_set_path (*fd_p, buf1);
}

///////////////////////////////////////////////////////////////////////////////
static int common_adios_set_path_var (long long fd_p, const char * path
                                     ,const char * name
                                     )
{
    struct adios_file_struct * fd = (struct adios_file_struct *) fd_p;
    struct adios_group_struct * t = fd->group;
    struct adios_var_struct * v = t->vars;

    // check for vars and then attributes
    v = adios_find_var_by_name (t->vars, name, fd->group->all_unique_var_names);

    if (v)
    {
        if (v->path)
        {
            free (v->path);
        }

        v->path = strdup (path);
    }
    else
    {
        fprintf (stderr, "adios_set_path_var (path=%s, var=%s): var not found\n"
                ,path, name
                );

        return 1;
    }

    return 0;
}

int adios_set_path_var (long long fd_p, const char * path, const char * name)
{
    return common_adios_set_path_var (fd_p, path, name);
}

void adios_set_path_var_ (long long * fd_p, const char * path
                         ,const char * name, int * err, int path_size
                         ,int name_size
                         )
{
    char buf1 [STR_LEN] = "";
    char buf2 [STR_LEN] = "";

    adios_extract_string (buf1, path, path_size);
    adios_extract_string (buf2, name, name_size);

    *err = common_adios_set_path_var (*fd_p, buf1, buf2);
}

///////////////////////////////////////////////////////////////////////////////
static int common_adios_end_iteration ()
{
    struct adios_method_list_struct * m;

    for (m = adios_get_methods (); m; m = m->next)
    {
        if (   m->method->m != ADIOS_METHOD_UNKNOWN
            && m->method->m != ADIOS_METHOD_NULL
            && adios_transports [m->method->m].adios_end_iteration_fn
           )
        {
            adios_transports [m->method->m].adios_end_iteration_fn
                                                (m->method);
        }
    }

    return 0;
}

// hint that we reached the end of an iteration (for asynchronous pacing)
int adios_end_iteration ()
{
    return common_adios_end_iteration ();
}

// hint that we reached the end of an iteration (for asynchronous pacing)
void adios_end_iteration_ (int * err)
{
    *err = common_adios_end_iteration ();
}

///////////////////////////////////////////////////////////////////////////////
static int common_adios_start_calculation ()
{
    struct adios_method_list_struct * m;

    for (m = adios_get_methods (); m; m = m->next)
    {
        if (   m->method->m != ADIOS_METHOD_UNKNOWN
            && m->method->m != ADIOS_METHOD_NULL
            && adios_transports [m->method->m].adios_start_calculation_fn
           )
        {
            adios_transports [m->method->m].adios_start_calculation_fn
                                                  (m->method);
        }
    }

    return 0;
}

// hint to start communicating
int adios_start_calculation ()
{
    return common_adios_start_calculation ();
}

// hint to start communicating
void adios_start_calculation_ (int * err)
{
    *err = common_adios_start_calculation ();
}

///////////////////////////////////////////////////////////////////////////////
static int common_adios_stop_calculation ()
{
    struct adios_method_list_struct * m;

    for (m = adios_get_methods (); m; m = m->next)
    {
        if (   m->method->m != ADIOS_METHOD_UNKNOWN
            && m->method->m != ADIOS_METHOD_NULL
            && adios_transports [m->method->m].adios_stop_calculation_fn
           )
        {
            adios_transports [m->method->m].adios_stop_calculation_fn
                                                   (m->method);
        }
    }

    return 0;
}

// hint to stop communicating
int adios_stop_calculation ()
{
    return common_adios_stop_calculation ();
}

// hint to stop communicating
void adios_stop_calculation_ (int * err)
{
    *err = common_adios_stop_calculation ();
}

///////////////////////////////////////////////////////////////////////////////
static int common_adios_close (long long fd_p)
{
    struct adios_file_struct * fd = (struct adios_file_struct *) fd_p;
    struct adios_method_list_struct * m = fd->group->methods;
    struct adios_attribute_struct * a = fd->group->attributes;
    struct adios_var_struct * v = fd->group->vars;

    if (fd->shared_buffer == adios_flag_yes)
    {
        adios_write_close_vars_v1 (fd);

        adios_write_open_attributes_v1 (fd);

        while (a)
        {
            adios_write_attribute_v1 (fd, a);

            a = a->next;
        }

        adios_write_close_attributes_v1 (fd);
    }

    // in order to get the index assembled, we need to do it in the
    // transport once we have collected all of the pieces

    // now tell all of the transports to write the buffer during close
    for (;m; m = m->next)
    {
        if (   m->method->m != ADIOS_METHOD_UNKNOWN
            && m->method->m != ADIOS_METHOD_NULL
            && adios_transports [m->method->m].adios_close_fn
           )
        {
            adios_transports [m->method->m].adios_close_fn
                                 (fd, m->method);
        }
    }

    while (v)
    {
        v->write_offset = 0;
        if (v->data)
            free (v->data);

        v->data = 0;

        v = v->next;
    }

    free ((void *) fd_p);

    return 0;
}

int adios_close (long long fd_p)
{
    return common_adios_close (fd_p);
}

void adios_close_ (long long * fd_p, int * err)
{
    *err = common_adios_close (*fd_p);
}

//////////////////////////////////////////////////////////////////////////////
// Methods normally only called by the XML parser
//////////////////////////////////////////////////////////////////////////////

// adios_common_declare_group is in adios_internals.c

///////////////////////////////////////////////////////////////////////////////
// group a list of vars into a composite group
int adios_declare_group (long long * id, const char * name
                        ,const char * coordination_comm
                        ,const char * coordination_var
                        ,const char * time_index
                        )
{
    return adios_common_declare_group (id, name, adios_flag_no
                                      ,coordination_comm
                                      ,coordination_var
                                      ,time_index
                                      );
}

void adios_declare_group_ (long long * id, const char * name
                          ,const char * coordination_comm
                          ,const char * coordination_var
                          ,const char * time_index, int * err
                          ,int name_size, int coordination_comm_size
                          ,int coordination_var_size, int time_index_size
                          )
{
    char buf1 [STR_LEN] = "";
    char buf2 [STR_LEN] = "";
    char buf3 [STR_LEN] = "";
    char buf4 [STR_LEN] = "";

    adios_extract_string (buf1, name, name_size);
    adios_extract_string (buf2, coordination_comm, coordination_comm_size);
    adios_extract_string (buf3, coordination_var, coordination_var_size);
    adios_extract_string (buf4, time_index, time_index_size);

    *err = adios_common_declare_group (id, buf1, adios_flag_yes, buf2
                                      ,buf3, buf4
                                      );
}

///////////////////////////////////////////////////////////////////////////////

// adios_common_define_var is in adios_internals.c

// declare a single var as an entry in a group
int adios_define_var (long long group_id, const char * name
                     ,const char * path, int type
                     ,const char * dimensions
                     ,const char * global_dimensions
                     ,const char * local_offsets
                     )
{
    return adios_common_define_var (group_id, name, path
                                   ,(enum ADIOS_DATATYPES) type
                                   ,dimensions
                                   ,global_dimensions, local_offsets
                                   );
}

// declare a single var as an entry in a group
void adios_define_var_ (long long * group_id, const char * name
                       ,const char * path, int * type
                       ,const char * dimensions
                       ,const char * global_dimensions
                       ,const char * local_offsets, int * err
                       ,int name_size, int path_size, int dimensions_size
                       ,int global_dimensions_size, int local_offsets_size
                       )
{
    char buf1 [STR_LEN] = "";
    char buf2 [STR_LEN] = "";
    char buf3 [STR_LEN] = "";
    char buf4 [STR_LEN] = "";
    char buf5 [STR_LEN] = "";

    adios_extract_string (buf1, name, name_size);
    adios_extract_string (buf2, path, path_size);
    adios_extract_string (buf3, dimensions, dimensions_size);
    adios_extract_string (buf4, global_dimensions, global_dimensions_size);
    adios_extract_string (buf5, local_offsets, local_offsets_size);

    *err = adios_common_define_var (*group_id, buf1, buf2
                                   ,(enum ADIOS_DATATYPES) *type
                                   ,buf3, buf4, buf5
                                   );
}

///////////////////////////////////////////////////////////////////////////////

// adios_common_define_attribute is in adios_internals.c

int adios_define_attribute (long long group, const char * name
                           ,const char * path, enum ADIOS_DATATYPES type
                           ,const char * value, const char * var
                           )
{
    return adios_common_define_attribute (group, name, path, type, value, var);
}

void adios_define_attribute_ (long long * group, const char * name
                             ,const char * path, int type, const char * value
                             ,const char * var, int * err
                             ,int name_size, int path_size, int value_size
                             ,int var_size
                             )
{
    char buf1 [STR_LEN] = "";
    char buf2 [STR_LEN] = "";
    char buf3 [STR_LEN] = "";
    char buf4 [STR_LEN] = "";

    adios_extract_string (buf1, name, name_size);
    adios_extract_string (buf2, path, path_size);
    adios_extract_string (buf3, value, value_size);
    adios_extract_string (buf4, var, var_size);

    *err = adios_common_define_attribute (*group, buf1, buf2
                                         ,(enum ADIOS_DATATYPES) type, buf3
                                         ,buf4
                                         );
}

///////////////////////////////////////////////////////////////////////////////

// adios_common_select_method is in adios_internals.c

int adios_select_method (int priority, const char * method
                        ,const char * parameters, const char * group
                        ,const char * base_path, int iters
                        )
{
    return adios_common_select_method (priority, method, parameters, group
                                      ,base_path, iters
                                      );
}

void adios_select_method_ (int * priority, const char * method
                          ,const char * parameters, const char * group
                          ,const char * base_path, int * iters, int * err
                          ,int method_size, int parameters_size
                          ,int group_size, int base_path_size
                          )
{
    char buf1 [STR_LEN] = "";
    char buf2 [STR_LEN] = "";
    char buf3 [STR_LEN] = "";
    char buf4 [STR_LEN] = "";

    adios_extract_string (buf1, method, method_size);
    adios_extract_string (buf2, parameters, parameters_size);
    adios_extract_string (buf3, group, group_size);
    adios_extract_string (buf4, base_path, base_path_size);

    *err = adios_common_select_method (*priority, buf1, buf2, buf3, buf4
                                      ,*iters
                                      );
}
