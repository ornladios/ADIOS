#include <unistd.h>
#include <fcntl.h>
#include "br-utils.h"

struct adios_bp_read_struct
{
    char type;  // 1 == file, 0 == buffer
    void * buf;
    uint64_t long buf_len;
    int fd;
    int offset;
};

long long br_fopen (char * filename)
{
    struct adios_bp_read_struct * b = (struct adios_bp_read_struct *)
                              malloc (sizeof (struct adios_bp_read_struct));

    if (!b)
    {
        fprintf (stderr, "Could not allocate memory in br_fopen ()\n");

        return 0;
    }

    b->type = 1;
    b->buf = 0;
    b->buf_len = 0;
    b->fd = open (filename, O_RDONLY);
    b->offset = 0;
    if (b->fd == -1)
    {
        fprintf (stderr, "Could not open %s\n", filename);

        free (b);

        return 0;
    }

    return (long long) b;
}

void br_fclose (long long handle)
{
    struct adios_bp_read_struct * b = (struct adios_bp_read_struct *) handle;

    if (b->fd != -1)
    {
        close (b->fd);
    }

    free (b);
}

long long br_bopen (void * buffer, uint64_t len)
{
    struct adios_bp_read_struct * b = (struct adios_bp_read_struct *)
                              malloc (sizeof (struct adios_bp_read_struct));

    if (!b)
    {
        fprintf (stderr, "Could not allocate memory in br_bopen ()\n");

        return 0;
    }

    b->type = 0;
    b->buf = buffer;
    b->buf_len = len;
    b->fd = -1;
    b->offset = 0;
    if (!b->buf)
    {
        fprintf (stderr, "invalid buffer provided to br_bopen\n");

        free (b);

        return 0;
    }

    return (long long) b;
}

void br_bclose (long long handle)
{
    struct adios_bp_read_struct * b = (struct adios_bp_read_struct *) handle;

    free (b);
}

void br_free_element (struct adios_bp_element_struct * element)
{
    if (element->dims)
        free (element->dims);
    if (element->name)
        free (element->name);
    if (element->path)
        free (element->path);

    free (element);
}

uint64_t br_read_buffer (struct adios_bp_read_struct * b
                        ,uint64_t read_size
                        ,void * buffer
                        ,uint64_t buffer_size
                        )
{
    int ret;

    switch (b->type)
    {
        case 0:
            if (buffer)
            {
                if (b->offset + read_size >= b->buf_len)
                {
                    read_size = b->buf_len - b->offset;
                }
                memcpy (buffer, b->buf + b->offset, read_size);
            }
            else
            {
                ; // do nothing
            }
            ret = read_size;
            break;

        case 1:
            if (buffer)
            {
                ret = read (b->fd, buffer, 1 * read_size);
            }
            else
            {
                lseek (b->fd, read_size, SEEK_CUR);
                ret = read_size;
            }
            break;
    }

    b->offset += read_size;

    return ret;
}

static int br_get_next_element (struct adios_bp_read_struct * b, ADIOS_BR_PRE_FETCH pre, ADIOS_BR_POST_FETCH post, void * private_data, void * buffer, uint64_t buffer_size, struct adios_bp_element_struct ** element)
{
    unsigned int var_step = 0;
    uint64_t total_size = 0;
    void * buf = 0;         // used in pre/post fetch calls
    uint64_t buf_size = 0;  // used in pre/post fetch calls
    int tag = NULL_TAG;
    unsigned long int size = 0;

    *element = (struct adios_bp_element_struct *)
                         calloc (1, sizeof (struct adios_bp_element_struct));

    if (!*element)
    {
        fprintf (stderr, "Failed to allocate memory for new element\n");

        return 0;
    }
    (*element)->data = buffer;
    (*element)->tag = NULL_TAG;

    total_size = br_read_buffer (b, sizeof (int), &var_step, sizeof (int));
    while (total_size < var_step)
    {
        br_read_buffer (b, sizeof (int), &tag, sizeof (int));
        total_size += sizeof (int);

        if ((*element)->tag == NULL_TAG)
        {
            (*element)->tag = tag;
        }

        br_read_buffer (b, sizeof (int), &size, sizeof (int));
        total_size += sizeof (int);

        if (buffer_size && size > buffer_size)
        {
            fprintf (stderr, "attempting to read element of %lu bytes, but "
                             "buffer only has %llu bytes of space\n"
                    ,size
                    ,buffer_size
                    );

            return 0;
        }

        switch (tag)
        {
            case SCR_TAG:
            case DST_TAG:
            case GRP_TAG:
            case DSTATRS_TAG:
            case DSTATRN_TAG:
            case GRPATRS_TAG:
            case GRPATRN_TAG:
                (*element)->name = malloc (size + 1);
                if (!(*element)->name)
                {
                    fprintf (stderr, "Cannot allocate %lu bytes during read\n"
                            ,size + 1
                            );

                    return 0;
                }
                br_read_buffer (b, size, (*element)->name, size + 1);
                (*element)->name [size] = '\0';
                break;

            case DIR_TAG:
                (*element)->path = malloc (size + 1);
                if (!(*element)->path)
                {
                    fprintf (stderr, "Cannot allocate %lu bytes during read\n"
                            ,size + 1
                            );

                    return 0;
                }
                br_read_buffer (b, size, (*element)->path, size + 1);
                (*element)->path [size] = '\0';
                break;

            case VAL_TAG:
                br_read_buffer (b, sizeof (int), &(*element)->type, sizeof (int));
                total_size += sizeof (int);

                (*element)->size = size;
                if ((*element)->type == bp_string)
                {
                    (*element)->size++; // account for the null
                }
                if (pre)
                {
                   pre (*element, &buf, &buf_size, private_data);
                }
                else
                {
                    buf = (*element)->data;
                    buf_size = buffer_size;
                }
                br_read_buffer (b, size, buf, buf_size);
                if ((*element)->type == bp_string)
                {
                    (*element)->size--; // account for the null adjustment above
                    if (buf)
                    {
                        ((char *) buf) [size] = '\0';
                    }
                }
                if (post)
                {
                    post (*element, buf, buf_size, private_data);
                }
                else
                {
                    // nothing to do
                }
                break;

            case DSTVAL_TAG:
                br_read_buffer (b, sizeof (int), &(*element)->ranks, sizeof (int));
                total_size += sizeof (int);
                (*element)->dims = (struct adios_bp_dimension_struct *)
                                 calloc ((*element)->ranks
                                        ,sizeof (struct adios_bp_dimension_struct)
                                        );
                if (!(*element)->dims)
                {
                    fprintf (stderr, "Cannot allocate %lu bytes during read\n"
                            ,  sizeof (struct adios_bp_dimension_struct)
                             * (*element)->ranks
                            );

                    return 0;
                }
                br_read_buffer (b
                               ,  sizeof (struct adios_bp_dimension_struct)
                                * (*element)->ranks
                               ,(*element)->dims
                               ,  sizeof (struct adios_bp_dimension_struct)
                                * (*element)->ranks
                               );
                total_size +=   sizeof (struct adios_bp_dimension_struct)
                              * (*element)->ranks;
                br_read_buffer (b, sizeof (int), &(*element)->type, sizeof (int));
                total_size += sizeof (int);

                (*element)->size = size;
                if (pre)
                {
                    pre (*element, &buf, &buf_size, private_data);
                }
                else
                {
                    buf = (*element)->data;
                    buf_size = buffer_size;
                }
                br_read_buffer (b, size, buf, buf_size);
                if (post)
                {
                    post (*element, buf, buf_size, private_data);
                }
                else
                {
                    // nothing to do
                }
                break;

            default:
                fprintf (stderr, "invalid tag during reading: %s(%d)\n"
                        ,adios_tag_to_string (tag)
                        ,tag
                        );

                return 0;
        }
        total_size += size;
    }

    return var_step;
}

int br_get_next_element_specific (long long handle, ADIOS_BR_PRE_FETCH pre, ADIOS_BR_POST_FETCH post, void * private_data, struct adios_bp_element_struct ** element)
{
    struct adios_bp_read_struct * b = (struct adios_bp_read_struct *) handle;

    return br_get_next_element (b, pre, post, private_data, 0, 0, element);
}

int br_get_next_element_general (long long handle, void * buffer, uint64_t buffer_size, struct adios_bp_element_struct ** element)
{
    struct adios_bp_read_struct * b = (struct adios_bp_read_struct *) handle;

    return br_get_next_element (b, 0, 0, 0, buffer, buffer_size, element);
}
