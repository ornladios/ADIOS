#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <unistd.h>

// mpi
#include "mpi.h"

// xml parser
#include <mxml.h>

// Chen's encoder
#include "bw-utils.h"
#include "br-utils.h"

#include "adios.h"
#include "adios_transport_hooks.h"
#include "adios_internals.h"

#define STR_LEN 1000

extern MPI_Comm adios_mpi_comm_world;
extern MPI_Comm adios_mpi_comm_self;
extern MPI_Info adios_mpi_info;

extern struct adios_transport_struct * adios_transports;

///////////////////////////////////////////////////////////////////////////////
static void common_adios_init (const char * config, MPI_Comm world, MPI_Comm self
                              ,MPI_Info info_obj
                              )
{
    MPI_Comm_dup (world, &adios_mpi_comm_world);
    MPI_Comm_dup (self, &adios_mpi_comm_self);
    if (info_obj != MPI_INFO_NULL)
        MPI_Info_dup (info_obj, &adios_mpi_info);
    else
        adios_mpi_info = MPI_INFO_NULL;

    // parse the config file
    adios_parse_config (config);
}

void adios_init (const char * config, MPI_Comm world, MPI_Comm self
                 ,MPI_Info info_obj
                 )
{
    common_adios_init (config, world, self, info_obj);
}

void adios_init_ (const char * config, MPI_Fint * world, MPI_Fint * self
                 ,MPI_Fint * info_obj
                 ,int config_size
                 )
{
    char buf1 [STR_LEN] = "";

    adios_extract_string (buf1, config, config_size);

    common_adios_init (buf1, MPI_Comm_f2c (*world), MPI_Comm_f2c (*self)
                      ,MPI_Info_f2c (*info_obj)
                      );
}

///////////////////////////////////////////////////////////////////////////////
static void common_adios_finalize (int mype)
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

    if (adios_mpi_info != MPI_INFO_NULL)
    {
        MPI_Info_free (&adios_mpi_info);
    }
    MPI_Comm_free (&adios_mpi_comm_self);
    MPI_Comm_free (&adios_mpi_comm_world);
}

void adios_finalize (int mype)
{
    common_adios_finalize (mype);
}

void adios_finalize_ (int * mype)
{
    common_adios_finalize (*mype);
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

int adios_allocate_buffer_ ()
{
    return common_adios_allocate_buffer ();
}

///////////////////////////////////////////////////////////////////////////////

// adios_common_get_group is in adios_internals.c

void adios_get_group (long long * group_id, const char * name)
{
    adios_common_get_group (group_id, name);
}

void adios_get_group_ (long long * group_id, const char * name, int name_size)
{
    char buf1 [STR_LEN] = "";

    adios_extract_string (buf1, name, name_size);

    adios_common_get_group (group_id, buf1);
}

///////////////////////////////////////////////////////////////////////////////
static void common_adios_open (long long * fd, long long group
                              ,const char * name, enum ADIOS_METHOD_MODE mode
                              )
{
    struct adios_file_struct * fd_p = (struct adios_file_struct *)
                                  malloc (sizeof (struct adios_file_struct));
    struct adios_group_struct * g = (struct adios_group_struct *) group;

    struct adios_method_list_struct * methods = g->methods;

    fd_p->name = strdup (name);
    fd_p->base_offset = 0;
    fd_p->offset = 0;
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
}

void adios_open (long long * fd, long long group, const char * name)
{
    common_adios_open (fd, group, name, adios_mode_write);
}

void adios_open_ (long long * fd, long long * group, const char * name
                 ,int name_size
                 )
{
    char buf1 [STR_LEN] = "";

    adios_extract_string (buf1, name, name_size);

    common_adios_open (fd, *group, buf1, adios_mode_write);
}

///////////////////////////////////////////////////////////////////////////////
void adios_open_append (long long * fd, long long group, const char * name)
{
    common_adios_open (fd, group, name, adios_mode_write | adios_mode_append);
}

void adios_open_append_ (long long * fd, long long * group, const char * name
                        ,int name_size
                        )
{
    common_adios_open (fd, *group, name, adios_mode_write | adios_mode_append);
}

///////////////////////////////////////////////////////////////////////////////
void adios_open_read (long long * fd, long long group, const char * name)
{
    common_adios_open (fd, group, name, adios_mode_read);
}

void adios_open_read_ (long long * fd, long long * group, const char * name
                      ,int name_size
                      )
{
    char buf1 [STR_LEN] = "";

    adios_extract_string (buf1, name, name_size);

    common_adios_open (fd, *group, buf1, adios_mode_read);
}

///////////////////////////////////////////////////////////////////////////////
static void common_adios_set_offset (long long fd_p, long long offset)
{
    struct adios_file_struct * fd = (struct adios_file_struct *) fd_p;

    if (fd)
    {
        if (!(fd->mode & adios_mode_write))
        {
            fprintf (stderr, "write attempted on %s which was opened for read\n"
                    ,fd->name
                    );

            return;
        }

        fd->offset = offset;
    }
    else
    {
        fprintf (stderr, "adios_write_offset (), invalid handle");
    }
}

void adios_set_offset (long long fd_p, long long offset)
{
    common_adios_set_offset (fd_p, offset);
}

void adios_set_offset_ (long long * fd_p, long long * offset)
{
    common_adios_set_offset (*fd_p, *offset);
}

///////////////////////////////////////////////////////////////////////////////
static void common_adios_write (long long fd_p, const char * name, void * var)
{
    struct adios_file_struct * fd = (struct adios_file_struct *) fd_p;
    struct adios_var_struct * v = fd->group->vars;
    struct adios_method_list_struct * m = fd->group->methods;

    v = adios_find_var_by_name (v, name);

    if (!v)
    {
        fprintf (stderr, "Bad var name (ignored): '%s'\n", name);

        return;
    }

    if (fd->mode & adios_mode_read && v->write_or_not == adios_flag_yes)
    {
        fprintf (stderr, "write attempted on %s in %s.  This was opened for"
                         " read and this is not a informational item %d\n"
                ,name , fd->name, v->write_or_not
                );

        return;
    }

    if (!var)
    {
        fprintf (stderr, "Invalid data: %s\n", name);

        return;
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
}

void adios_write (long long fd_p, const char * name, void * var)
{
    common_adios_write (fd_p, name, var);
}

void adios_write_ (long long * fd_p, const char * name, void * var
                  ,int name_size
                  )
{
    char buf1 [STR_LEN] = "";

    adios_extract_string (buf1, name, name_size);

    common_adios_write (*fd_p, buf1, var);
}


///////////////////////////////////////////////////////////////////////////////
static void common_adios_get_write_buffer (long long fd_p, const char * name
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

        return;
    }

    if (fd->mode & adios_mode_read && v->write_or_not == adios_flag_yes)
    {
        fprintf (stderr, "write attempted on %s in %s.  This was opened for"
                         " read and this is not a informational item %d\n"
                ,name , fd->name, v->write_or_not
                );

        return;
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
}

void adios_get_write_buffer (long long fd_p, const char * name
                            ,unsigned long long * size
                            ,void ** buffer
                            )
{
    common_adios_get_write_buffer (fd_p, name, size, buffer);
}

void adios_get_write_buffer_ (long long * fd_p, const char * name
                             ,unsigned long long * size
                             ,void ** buffer, int name_size
                             )
{
    char buf1 [STR_LEN] = "";

    adios_extract_string (buf1, name, name_size);

    common_adios_get_write_buffer (*fd_p, buf1, size, buffer);
}

///////////////////////////////////////////////////////////////////////////////
static void common_adios_read (long long fd_p, const char * name, void * buffer)
{
    struct adios_file_struct * fd = (struct adios_file_struct *) fd_p;
    struct adios_var_struct * v;

    if (!(fd->mode & adios_mode_read))
    {
        fprintf (stderr, "read attempted on %s which was opened for write\n"
                ,fd->name
                );

        return;
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
    }
}

void adios_read (long long fd_p, const char * name, void * buffer)
{
    common_adios_read (fd_p, name, buffer);
}

void adios_read_ (long long * fd_p, const char * name, void * buffer
                 ,int name_size
                 )
{
    char buf1 [STR_LEN] = "";

    adios_extract_string (buf1, name, name_size);

    common_adios_read (*fd_p, buf1, buffer);
}

///////////////////////////////////////////////////////////////////////////////
static void common_adios_set_path (long long group, const char * path)
{
    struct adios_group_struct * t = (struct adios_group_struct *) group;
    struct adios_var_struct * v = t->vars;

    while (v)
    {
        if (v->path)
        {
            free (v->path);
        }

        v->path = strdup (path);

        v = v->next;
    }
}

void adios_set_path (long long group, const char * path)
{
    common_adios_set_path (group, path);
}

void adios_set_path_ (long long * group, const char * path, int path_size)
{
    char buf1 [STR_LEN] = "";

    adios_extract_string (buf1, path, path_size);

    common_adios_set_path (*group, buf1);
}

///////////////////////////////////////////////////////////////////////////////
static void common_adios_set_path_var (long long group, const char * path
                                      ,const char * name
                                      )
{
    struct adios_group_struct * t = (struct adios_group_struct *) group;
    struct adios_var_struct * v = t->vars;

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
    }
}

void adios_set_path_var (long long group, const char * path, const char * name)
{
    common_adios_set_path_var (group, path, name);
}

void adios_set_path_var_ (long long * group, const char * path
                         ,const char * name, int path_size, int name_size
                         )
{
    char buf1 [STR_LEN] = "";
    char buf2 [STR_LEN] = "";

    adios_extract_string (buf1, path, path_size);
    adios_extract_string (buf2, name, name_size);

    common_adios_set_path_var (*group, buf1, buf2);
}

///////////////////////////////////////////////////////////////////////////////
static void common_adios_get_data_size (long long fd_p
                                       ,unsigned long long * size
                                       )
{
    struct adios_file_struct * fd = (struct adios_file_struct *) fd_p;

    *size = adios_data_size (fd->group);
}

void adios_get_data_size (long long fd_p, unsigned long long * size)
{
    common_adios_get_data_size (fd_p, size);
}

void adios_get_data_size_ (long long * fd_p, unsigned long long * size)
{
    common_adios_get_data_size (*fd_p, size);
}

///////////////////////////////////////////////////////////////////////////////
static void common_adios_end_iteration ()
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
}

// hint that we reached the end of an iteration (for asynchronous pacing)
void adios_end_iteration ()
{
    common_adios_end_iteration ();
}

// hint that we reached the end of an iteration (for asynchronous pacing)
void adios_end_iteration_ ()
{
    common_adios_end_iteration ();
}

///////////////////////////////////////////////////////////////////////////////
static void common_adios_start_calculation ()
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
}

// hint to start communicating
void adios_start_calculation ()
{
    common_adios_start_calculation ();
}

// hint to start communicating
void adios_start_calculation_ ()
{
    common_adios_start_calculation ();
}

///////////////////////////////////////////////////////////////////////////////
static void common_adios_stop_calculation ()
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
}

// hint to stop communicating
void adios_stop_calculation ()
{
    common_adios_stop_calculation ();
}

// hint to stop communicating
void adios_stop_calculation_ ()
{
    common_adios_stop_calculation ();
}

///////////////////////////////////////////////////////////////////////////////
static void common_adios_close (long long fd_p)
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
}

void adios_close (long long fd_p)
{
    common_adios_close (fd_p);
}

void adios_close_ (long long * fd_p)
{
    common_adios_close (*fd_p);
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
    return adios_common_declare_group (id, name, coordination_comm
                                      ,coordination_var
                                      );
}

int adios_declare_group_ (long long * id, const char * name
                         ,const char * coordination_comm
                         ,const char * coordination_var
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

    return adios_common_declare_group (id, buf1, buf2, buf3);
}

///////////////////////////////////////////////////////////////////////////////

// adios_common_define_var is in adios_internals.c

// declare a single var as an entry in a group
int adios_define_var (long long group_id, const char * name
                     ,const char * path, int type
                     ,int write_or_not, int copy_on_write
                     ,const char * dimensions
                     ,struct adios_global_bounds_struct * global_dimensions
                     )
{
    return adios_common_define_var (group_id, name, path, type, write_or_not
                                   ,copy_on_write, dimensions
                                   ,global_dimensions
                                   );
}

// declare a single var as an entry in a group
int adios_define_var_ (long long * group_id, const char * name
                      ,const char * path, int * type
                      ,int * write_or_not, int * copy_on_write
                      ,const char * dimensions
                      ,long long * global_dimensions
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

    return adios_common_define_var (*group_id, buf1, buf2, *type
                                   ,*write_or_not
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

int adios_define_global_bounds_ (long long * group, const char * dimensions
                                ,const char * offsets
                                ,long long * global_bounds
                                ,int dimensions_size, int offsets_size
                                )
{
    struct adios_global_bounds_struct ** b =
                  (struct adios_global_bounds_struct **) global_bounds;
    char buf1 [STR_LEN] = "";
    char buf2 [STR_LEN] = "";

    adios_extract_string (buf1, dimensions, dimensions_size);
    adios_extract_string (buf2, offsets, offsets_size);

    return adios_common_define_global_bounds (*group, buf1, buf2, b);
}

///////////////////////////////////////////////////////////////////////////////

// adios_common_define_attribute is in adios_internals.c

//CJ
int adios_define_attribute (long long group, const char * name
                           ,const char * path, const char * value, int type
                           )
{
//CJ
    return adios_common_define_attribute (group, name, path, value,type);
}

//CJ
int adios_define_attribute_ (long long * group, const char * name
                            ,const char * path, const char * value, int * type
                            ,int name_size, int path_size, int value_size
                            )
{
    char buf1 [STR_LEN] = "";
    char buf2 [STR_LEN] = "";
    char buf3 [STR_LEN] = "";
//CJ
    adios_extract_string (buf1, name, name_size);
    adios_extract_string (buf2, path, path_size);
    adios_extract_string (buf3, value, value_size);
//CJ
    return adios_common_define_attribute (*group, buf1, buf2, buf3, *type);
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

int adios_select_method_ (int * priority, const char * method
                         ,const char * parameters, const char * group
                         ,const char * base_path, int * iters, int method_size
                         ,int parameters_size, int group_size, int base_path_size
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

    return adios_common_select_method (*priority, buf1, buf2, buf3, buf4, *iters);
}
