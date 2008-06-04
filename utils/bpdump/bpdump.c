#include "binpack-general.h"
#include "br-utils.h"

int print_dataset (int type, int ranks, struct adios_bp_dimension_struct * dims
                  ,void * data
                  );

struct dump_struct
{
    int dump;
    char * dump_var;
    unsigned long DATALEN;
    void * val;
};

static void pre_element_fetch (struct adios_bp_element_struct * element
                              ,void ** buffer, long long * buffer_size
                              ,void * private_data
                              )
{
    struct dump_struct * data = (struct dump_struct *) private_data;

    if (element->tag == DST_TAG)
    {
        if (data->dump == 0)
        {
            element->data = 0;
            *buffer = 0;
            *buffer_size = 0;

            return;
        }
    }

    element->data = data->val;
    *buffer = data->val;
    *buffer_size = data->DATALEN;
}

int main (int argc, char ** argv)
{
    char * filename;
    char * var;
    int i = 0;
    long long handle = 0;
    long long element_size = 0;
    struct adios_bp_element_struct * element = NULL;
    struct dump_struct data;
    data.DATALEN = 100 * 1024 * 1024;

    if (argc < 2)
    {
        fprintf (stderr, "usage: %s [-d [var]|--dump [var]] <filename>\n", argv [0]);

        return -1;
    }

    if (argv [1][0] && argv [1][0] == '-')
    {
        if (   !strcmp (argv [1], "-d")
            || !strcmp (argv [1], "--dump")
           )
        {
            data.dump = 1;
            if (argc > 2)
            {
                data.dump_var = argv [2];
                filename = argv [3];
                printf("%s %s\n",data.dump_var,filename);
            }
            else
            {
                data.dump_var = 0;
                filename = argv [2];
                printf("%s %s\n",data.dump_var,filename);
            }
        }
        else
        {
            fprintf (stderr, "usage: %s [-d [var]|--dump [var]] <filename>\n", argv [0]);

            return -1;
        }
    }
    else
    {
        filename = argv [1];
        data.dump = 0;
        data.dump_var = 0;
    }

    handle = br_fopen (filename);
    if (!handle)
    {
        fprintf (stderr, "file not found: %s\n", filename);

        return -1;
    }

    data.val = (char *) malloc (data.DATALEN);
    if (!data.val)
    {
        fprintf (stderr, "cannot allocate %ul for data buffer\n", data.DATALEN);

        return -1;
    }

    while (element_size = br_get_next_element_specific (handle
                                                       ,pre_element_fetch
                                                       ,0
                                                       ,&data
                                                       ,&element
                                                       )
          )
    {
        printf ("element size: %d\n", element_size);
        printf ("%s %s\n", adios_tag_to_string (element->tag), element->name);
        printf ("\tPath: %s\n", element->path);

        switch (element->tag)
        {
            case SCR_TAG:
            case DSTATRS_TAG:
            case DSTATRN_TAG:
            case GRPATRS_TAG:
            case GRPATRN_TAG:
                printf ("\tType: %s (%d)\n"
                       ,adios_type_to_string (element->type)
                       ,element->type
                       );
                switch (element->type)
                {
                    case bp_int: //adios_integer:
                        printf ("\t(int): %d\n", *((int *) element->data));
                        break;
                    case bp_float: //adios_real:
                        printf ("\t(float): %f\n", *((float *) element->data));
                        break;
                    case bp_string: //adios_string:
                        printf ("\t(char): %s\n", ((char *) element->data));
                        break;
                    case bp_double: //adios_double:
                        printf ("\t(double): %lf\n"
                               ,*((double *) element->data)
                               );
                        break;
                    case bp_uchar: //adios_byte:
                        printf ("\t(byte): %d\n"
                               ,*((unsigned char *) element->data)
                               );
                        break;
                    case bp_longlong: // adios_long
                        printf ("\t(long): %lld\n"
                               ,*((long long *) element->data)
                               );
                        break;
                    case bp_complex: // adios_complex
                        printf ("\t(complex): %lf %lf\n"
                               ,((double *) element->data) [0]
                               ,((double *) element->data) [1]
                               );
                        break;
                    default:
                        break;
                }
                break;

            case DST_TAG:
                printf ("\tType: %s (%d)\n"
                       ,adios_type_to_string (element->type)
                       ,element->type
                       );
                printf ("\tRanks: %u\n", element->ranks);
                if (element->dims)
                {
                    for (i = 0; i < element->ranks; i++)   
                    {                                      
                        printf ("\t\tGlobal Dim(%d): %d(%d):%d(%d):%d\n"
                               ,i
                               ,element->dims [i].local_bound
                               ,element->dims [i].global_bound
                               ,element->dims [i].global_offset
                           );
                    }
                }
                for (i = 0; i < element->ranks; i++)
                {
                    printf ("\t\tDim(%d): %d(%d):%d(%d):%d\n"
                           ,i
                           ,element->dims [i].local_bound
                           ,element->dims [i].global_bound
                           ,element->dims [i].global_offset
                           );
                }
                if (data.dump)
                {
                    printf ("\t%s: %s\n", data.dump_var, element->name);
                    if (!data.dump_var || !strcmp (data.dump_var, element->name))
                    {
                        printf ("\tdata.dump: %u\n", data.dump);
                        print_dataset (element->type
                                      ,element->ranks
                                      ,element->dims
                                      ,element->data
                                      );
                    }
                }
                break;

            case GRP_TAG:
            default:
                fprintf (stderr, "invalid tag\n");

                break;
        }
        br_free_element (element);
        printf ("-----------------------\n");
    }
    br_fclose (handle);


    return 0;
}

int print_dataset (int type, int ranks, struct adios_bp_dimension_struct * dims
                  ,void * data
                  )
{
    printf("My Value: \n");
    int use_element_count = 1;
    int total_element_count = 1;
    int e = 0; // which element we are printing
    int i,j;
    int position[ranks];
    for (i = 0; i < ranks; i++)
    {
        use_element_count *= (  dims [i].local_bound
                             );
        position[i]=dims[i].local_bound; 
        total_element_count *= dims [i].local_bound;
    } 
    printf ("DIMS: dims[%d][%d]\n",position[0],position[1]);
    printf ("ranks=%d\n",ranks);
    for (j = 0; j < total_element_count; j++)
    {
            switch (type)
            {
                case bp_uchar:
                    printf ("%c ", (((unsigned char *) data) [e]));
                    break;

                case bp_char:
                    printf ("%c ", (((char *) data) [e]));
                    break;

                case bp_int:
                    printf ("%d ", (((int *) data) [e]));
                    break;

                case bp_float:
                    printf ("(%d: %e) ", j,(((float *) data) [e]));
                    break;

                case bp_double:
                    printf ("(%d: %le) ",j, (((double *) data) [e]));
                    break;

                case bp_string:
                    break;

                case bp_longlong: // adios_long
                    printf ("%lld ", (((long long *) data) [e]));
                    break;

                case bp_complex: // adios_complex
                    printf ("(%lf %lf) ", ((double *) data) [e * 2 + 0]
                                        , ((double *) data) [e * 2 + 1]
                           );
                    break;
            }
            e++;
        }
}
