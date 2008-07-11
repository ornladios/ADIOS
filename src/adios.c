#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <unistd.h>

// xml parser
#include <mxml.h>

// Chen's encoder
#include "bw-utils.h"
#include "br-utils.h"

#include "adios.h"
#include "adios_transport_hooks.h"
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

void adios_init_ (const char * config, int err, int config_size)
{
    char buf1 [STR_LEN] = "";

    adios_extract_string (buf1, config, config_size);

    err = common_adios_init (buf1);
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

void adios_finalize_ (int * mype, int err)
{
    err = common_adios_finalize (*mype);
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

void adios_allocate_buffer_ (int err)
{
    err = common_adios_allocate_buffer ();
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

    if (!strcmp (file_mode, "r"))
        mode = adios_mode_read;
    else
        if (!strcmp (file_mode, "w"))
            mode = adios_mode_write;
        else
            if (!strcmp (file_mode, "a"))
                mode = adios_mode_append;
            else
                if (!strcmp (file_mode, "u"))
                    mode = adios_mode_update;
                else
                {
                    fprintf (stderr, "adios_open: unknown file mode: %s\n"
                            ,file_mode
                            );

                    *fd = 0;

                    return;
                }

    fd_p->name = strdup (name);
    fd_p->base_offset = 0;
    fd_p->offset = 0;
    fd_p->write_size = 0;
    fd_p->group = g;
    fd_p->mode = mode;

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
                 ,const char * mode, int err
                 ,int group_name_size, int name_size, int mode_size
                 )
{
    char buf1 [STR_LEN] = "";
    char buf2 [STR_LEN] = "";
    char buf3 [STR_LEN] = "";

    adios_extract_string (buf1, group_name, group_name_size);
    adios_extract_string (buf2, name, name_size);
    adios_extract_string (buf3, mode, mode_size);

    err = common_adios_open (fd, buf1, buf2, buf3);
}

///////////////////////////////////////////////////////////////////////////////
static int common_adios_group_size (long long fd_p, int nvars
                                   ,long long byte_size
                                   )
{
    struct adios_file_struct * fd = (struct adios_file_struct *) fd_p;

    fd->write_size = byte_size;

    return 0;
}

int adios_group_size (long long fd_p, int nvars, long long byte_size)
{
    return common_adios_group_size (fd_p, nvars, byte_size);
}

void adios_group_size_ (long long fd_p, int nvars, long long byte_size, int err)
{
    err = common_adios_group_size (fd_p, nvars, byte_size);
}

///////////////////////////////////////////////////////////////////////////////
static int common_adios_write (long long fd_p, const char * name, void * var)
{
    struct adios_file_struct * fd = (struct adios_file_struct *) fd_p;
    struct adios_var_struct * v = fd->group->vars;
    struct adios_method_list_struct * m = fd->group->methods;

    v = adios_find_var_by_name (v, name);

    if (!v)
    {
        fprintf (stderr, "Bad var name (ignored): '%s'\n", name);

        return 1;
    }

    if (fd->mode == adios_mode_read)
    {
        if (strcmp (name, fd->group->group_comm) && v->is_dim != adios_flag_yes)
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

    return 0;
}

int adios_write (long long fd_p, const char * name, void * var)
{
    return common_adios_write (fd_p, name, var);
}

void adios_write_ (long long * fd_p, const char * name, void * var, int err
                  ,int name_size
                  )
{
    char buf1 [STR_LEN] = "";

    adios_extract_string (buf1, name, name_size);

    err = common_adios_write (*fd_p, buf1, var);
}


///////////////////////////////////////////////////////////////////////////////
static int common_adios_get_write_buffer (long long fd_p, const char * name
                                         ,unsigned long long * size
                                         ,void ** buffer
                                         )
{
    struct adios_file_struct * fd = (struct adios_file_struct *) fd_p;
    struct adios_var_struct * v = fd->group->vars;
    struct adios_method_list_struct * m = fd->group->methods;

    v = adios_find_var_by_name (v, name);

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
                           ,unsigned long long * size
                           ,void ** buffer
                           )
{
    return common_adios_get_write_buffer (fd_p, name, size, buffer);
}

void adios_get_write_buffer_ (long long * fd_p, const char * name
                             ,unsigned long long * size
                             ,void ** buffer, int err, int name_size
                             )
{
    char buf1 [STR_LEN] = "";

    adios_extract_string (buf1, name, name_size);

    err = common_adios_get_write_buffer (*fd_p, buf1, size, buffer);
}

///////////////////////////////////////////////////////////////////////////////
static int common_adios_read (long long fd_p, const char * name, void * buffer)
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

    v = adios_find_var_by_name (fd->group->vars, name);
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

int adios_read (long long fd_p, const char * name, void * buffer)
{
    return common_adios_read (fd_p, name, buffer);
}

void adios_read_ (long long * fd_p, const char * name, void * buffer, int err
                 ,int name_size
                 )
{
    char buf1 [STR_LEN] = "";

    adios_extract_string (buf1, name, name_size);

    err = common_adios_read (*fd_p, buf1, buffer);
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

#if 0
    while (a)
    {
        if (a->path)
        {
            free (a->path);
        }

        a->path = strdup (path);

        a = a->next;
    }
#endif

    return 0;
}

int adios_set_path (long long fd_p, const char * path)
{
    return common_adios_set_path (fd_p, path);
}

void adios_set_path_ (long long * fd_p, const char * path, int err
                     ,int path_size
                     )
{
    char buf1 [STR_LEN] = "";

    adios_extract_string (buf1, path, path_size);

    err = common_adios_set_path (*fd_p, buf1);
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
    v = adios_find_var_by_name (t->vars, name);

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
                         ,const char * name, int err, int path_size
                         ,int name_size
                         )
{
    char buf1 [STR_LEN] = "";
    char buf2 [STR_LEN] = "";

    adios_extract_string (buf1, path, path_size);
    adios_extract_string (buf2, name, name_size);

    err = common_adios_set_path_var (*fd_p, buf1, buf2);
}

#if 0
///////////////////////////////////////////////////////////////////////////////
static int common_adios_get_data_size (long long fd_p
                                      ,unsigned long long * size
                                      )
{
    struct adios_file_struct * fd = (struct adios_file_struct *) fd_p;

    *size = adios_data_size (fd->group);

    return 0;
}

int adios_get_data_size (long long fd_p, unsigned long long * size)
{
    return common_adios_get_data_size (fd_p, size);
}

void adios_get_data_size_ (long long * fd_p, unsigned long long * size, int err)
{
    err = common_adios_get_data_size (*fd_p, size);
}
#endif

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
void adios_end_iteration_ (int err)
{
    err = common_adios_end_iteration ();
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
void adios_start_calculation_ (int err)
{
    err = common_adios_start_calculation ();
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
void adios_stop_calculation_ (int err)
{
    err = common_adios_stop_calculation ();
}

///////////////////////////////////////////////////////////////////////////////
static int common_adios_close (long long fd_p)
{
    struct adios_file_struct * fd = (struct adios_file_struct *) fd_p;
    struct adios_method_list_struct * m = fd->group->methods;

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

    free ((void *) fd_p);

    return 0;
}

int adios_close (long long fd_p)
{
    return common_adios_close (fd_p);
}

void adios_close_ (long long * fd_p, int err)
{
    err = common_adios_close (*fd_p);
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
                        )
{
    return adios_common_declare_group (id, name, adios_flag_no
                                      ,coordination_comm
                                      ,coordination_var
                                      );
}

void adios_declare_group_ (long long * id, const char * name
                          ,const char * coordination_comm
                          ,const char * coordination_var, int err
                          ,int name_size, int coordination_comm_size
                          ,int coordination_var_size
                          )
{
    char buf1 [STR_LEN] = "";
    char buf2 [STR_LEN] = "";
    char buf3 [STR_LEN] = "";

    adios_extract_string (buf1, name, name_size);
    adios_extract_string (buf2, coordination_comm, coordination_comm_size);
    adios_extract_string (buf3, coordination_var, coordination_var_size);

    err = adios_common_declare_group (id, buf1, adios_flag_yes, buf2, buf3);
}

///////////////////////////////////////////////////////////////////////////////

// adios_common_define_var is in adios_internals.c

// declare a single var as an entry in a group
int adios_define_var (long long group_id, const char * name
                     ,const char * path, int type
                     ,int copy_on_write
                     ,const char * dimensions
                     ,struct adios_global_bounds_struct * global_dimensions
                     )
{
    return adios_common_define_var (group_id, name, path, type
                                   ,copy_on_write, dimensions
                                   ,global_dimensions
                                   );
}

// declare a single var as an entry in a group
void adios_define_var_ (long long * group_id, const char * name
                       ,const char * path, int * type
                       ,int * copy_on_write
                       ,const char * dimensions
                       ,long long * global_dimensions, int err
                       ,int name_size, int path_size, int dimensions_size
                       )
{
    struct adios_global_bounds_struct ** b =
                   (struct adios_global_bounds_struct **) global_dimensions;

    char buf1 [STR_LEN] = "";
    char buf2 [STR_LEN] = "";
    char buf3 [STR_LEN] = "";

    adios_extract_string (buf1, name, name_size);
    adios_extract_string (buf2, path, path_size);
    adios_extract_string (buf3, dimensions, dimensions_size);

    err = adios_common_define_var (*group_id, buf1, buf2, *type
                                  ,*copy_on_write, buf3, *b
                                  );
}

///////////////////////////////////////////////////////////////////////////////

// adios_common_define_var is in adios_internals.c

int adios_define_global_bounds (long long group, const char * dimensions
                               ,const char * offsets
                               ,struct adios_global_bounds_struct ** b
                               )
{
    return adios_common_define_global_bounds (group, dimensions, offsets, b);
}

void adios_define_global_bounds_ (long long * group, const char * dimensions
                                 ,const char * offsets
                                 ,long long * global_bounds, int err
                                 ,int dimensions_size, int offsets_size
                                 )
{
    struct adios_global_bounds_struct ** b =
                  (struct adios_global_bounds_struct **) global_bounds;
    char buf1 [STR_LEN] = "";
    char buf2 [STR_LEN] = "";

    adios_extract_string (buf1, dimensions, dimensions_size);
    adios_extract_string (buf2, offsets, offsets_size);

    err = adios_common_define_global_bounds (*group, buf1, buf2, b);
}

///////////////////////////////////////////////////////////////////////////////

// adios_common_define_attribute is in adios_internals.c

int adios_define_attribute (long long group, const char * name
                           ,const char * path, const char * value
                           ,const char * var
                           )
{
    return adios_common_define_attribute (group, name, path, value, var);
}

void adios_define_attribute_ (long long * group, const char * name
                             ,const char * path, const char * value
                             ,const char * var, int err
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

    err = adios_common_define_attribute (*group, buf1, buf2, buf3, buf4);
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
                          ,const char * base_path, int * iters, int err
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

    err = adios_common_select_method (*priority, buf1, buf2, buf3, buf4
                                     ,*iters
                                     );
}
