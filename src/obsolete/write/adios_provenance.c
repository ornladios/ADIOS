/* 
 * ADIOS is freely available under the terms of the BSD license described
 * in the COPYING file in the top level directory of this source distribution.
 *
 * Copyright (c) 2008 - 2009.  UT-BATTELLE, LLC. All rights reserved.
 */

#include <unistd.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

// xml parser
#include <mxml.h>

#include "public/adios_mpi.h" 
#include "core/adios_transport_hooks.h"
#include "core/adios_bp_v1.h"
#include "core/adios_internals.h"
#include "core/buffer.h"
#include "core/util.h"

/** A method to allow apps sending provenance info 
    - through a socket or
    - into a file.

    Information sent/written:
    - timestep
    - list of full paths of files (created by other methods)
    - index in bp binary format

    It is the world master (communicator in adios_init) process which communicates to the
    external world. If there are several groups in the application 
    communicating with each other, their masters will communicate to 
    the world master. 

 */

/** Implementation as intended
    --------------------------

    This method supports the NxM groups. 
    Information flows from N process --> M group writers --> 1 world master.
    Only the world master writes/sends data. 

    MPI is used for the communication, similarly as in the MPI method.
    
    Functions:
      init()
         Just initialize the provenance struct and decide what output method to use
      open()
         Open file/socket by world master
      close()
         Communicate indices N-->M-->1
         Write/send information 

 */

/** Current problems and todos
    --------------------------

    1. TODO: socket method
    2. TODO: collect all file names being written and send them out too.
    3. If the socket cannot be opened, we have to
       - store the info in a file and
       - send that "old" info to the socket in a later step if
         the socket will be available

 */

/* ======== PRIVATE FUNCTION DECLARATIONS ======== */
static void _check_other_methods(struct adios_file_struct * fd
                                ,struct adios_method_struct * method
                                );

static void _look_for_listener(char *env_var_name
        ,char *filename
        ,char *path
        ,char *alternate_path
        ,int   default_port
        ,char **host
        ,int  *port
        );

static int _connect_to_listener (char * host, int port);

static void _fix_offsets(struct adios_index_process_group_struct_v1 * pg_root 
                        ,struct adios_index_var_struct_v1 * vars_root
                         ,uint64_t offset);

static void _debug_print_index(char *ext
                              ,struct adios_index_process_group_struct_v1 * pg_root
                              ,struct adios_index_var_struct_v1 * vars_root
                              );


/* ======== PRIVATE VARIABLES AND TYPES ========== */
static int adios_provenance_initialized = 0;


enum ADIOS_PROVENANCE_OUTPUTMETHOD {adios_provenance_outputmethod_file = 1
                                   ,adios_provenance_outputmethod_socket = 2
};



/* If output method is file, the name of output file is this*/
static char ADIOS_PROVENANCE_OUTPUTFILE[] = "provenance.out";
/* If output method is socket, the file containing the socket info is this */
static char ADIOS_PROVENANCE_SOCKETINFOFILE[] = "provenance.conf";
/* If output method is socket, the environment variable containing the path of the socket info file is this */
static char ADIOS_PROVENANCE_SOCKETINFOENV[] = "ADIOS_PROVENANCE_SOCKETINFOFILE";

static char LOGHDR[]="provenance method:";

struct adios_PROVENANCE_data_struct
{
    MPI_Status status;   // status for MPI queries here to allocate memory 
                         //   in the beginning and not during its use 
    MPI_Comm group_comm; // communication for the group in an open
    int rank;            // rank of this process within this group
    int size;            // number of processes in this group

    MPI_Comm world_comm; // all processes calling adios_open..adios_close

    enum   ADIOS_PROVENANCE_OUTPUTMETHOD output_method;

    enum   ADIOS_IO_METHOD another_method;  // the other method in case of multiple methods

    struct adios_bp_buffer_struct_v1 b;

    struct adios_index_process_group_struct_v1 * old_pg_root; // process group index
    struct adios_index_var_struct_v1 * old_vars_root;         // vars index
    struct adios_index_attribute_struct_v1 * old_attrs_root;  // attributes index
};


/* ======== PUBLIC FUNCTIONS  ========== */


void adios_provenance_init (const PairStruct * parameters
                           ,struct adios_method_struct * method
                           )
{
    struct adios_PROVENANCE_data_struct * md = 0;
    if (!adios_provenance_initialized)
    {
        adios_provenance_initialized = 1;
    }
    method->method_data = malloc (sizeof (struct adios_PROVENANCE_data_struct));
    md = (struct adios_PROVENANCE_data_struct *) method->method_data;

    /* for testing we start with file output method */
    md->output_method = adios_provenance_outputmethod_file;
    /* FIXME: parse parameters and set output method info here */

    
    memset (&md->status, 0, sizeof (MPI_Status));
    md->rank = 0;  // will be set in ..open()
    md->size = 0;  // will be set in ..open()
    md->group_comm = MPI_COMM_NULL; // unused, adios_open sets current comm
    md->world_comm = method->init_comm; 
    md->old_pg_root = 0;
    md->old_vars_root = 0;
    md->old_attrs_root = 0;
    md->another_method = ADIOS_METHOD_NULL; // no other method, may be changed in open()

    adios_buffer_struct_init (&md->b);


}

    
    
int adios_provenance_open (struct adios_file_struct * fd
                          ,struct adios_method_struct * method, MPI_Comm comm
                          )
{
    struct adios_PROVENANCE_data_struct * md = (struct adios_PROVENANCE_data_struct *)
                                                    method->method_data;

    int wrank;

    // TODO: record the filename fd->filename, which is different for each group
    //char *name;
    //name = malloc (strlen (method->base_path) + strlen (fd->name) + 1);
    //sprintf (name, "%s%s", method->base_path, fd->name);
        

    _check_other_methods( fd, method); // set md->another_method
    
    // Get the MPI related info
    md->group_comm = comm;
    if (md->group_comm != MPI_COMM_NULL)
    {
        MPI_Comm_rank (md->group_comm, &md->rank);
        MPI_Comm_size (md->group_comm, &md->size);
    }
    else {
        md->size = 1;  // each process is a group 
    }

    MPI_Comm_rank (md->world_comm, &wrank);
    // The worldmaster can open the file/socket/dbconnection now. 
    if (wrank == 0) 
    {
        switch (fd->mode)
        {
            case adios_mode_write:
            case adios_mode_append:
            {
                switch (md->output_method)
                {
                    // Sending through a socket
                    case (adios_provenance_outputmethod_socket):
                    {
                        /* Look for info about a listener. (host:port formatted string in a file)
                           1. check file named in the env variable named in ADIOS_PROVENANCE_SOCKETINFOENV[] 
                           2. check file named in ADIOS_PROVENANCE_SOCKETINFOFILE[] in . 
                           3. check file named in ADIOS_PROVENANCE_SOCKETINFOFILE[] in $PBS_O_WORKDIR
                         */        
                        char *listener_host;
                        int   listener_port;
                        _look_for_listener(ADIOS_PROVENANCE_SOCKETINFOENV 
                                ,ADIOS_PROVENANCE_SOCKETINFOFILE
                                ,"."
                                ,getenv("PBS_O_WORKDIR")
                                ,-1
                                ,&listener_host
                                ,&listener_port
                                );
                        printf ("open(): socket info host:port=%s:%d.\n", listener_host, listener_port);
                        break;
                    } // case: to socket

                    // Writing into a file
                    case (adios_provenance_outputmethod_file):
                    {
                        md->b.f = open (ADIOS_PROVENANCE_OUTPUTFILE, O_RDWR);
                        if (md->b.f == -1)
                        {
                            md->b.f = open (ADIOS_PROVENANCE_OUTPUTFILE
                                    ,O_WRONLY | O_CREAT | O_TRUNC
#ifndef __APPLE__
| O_LARGEFILE
#endif
                                    ,S_IRUSR | S_IWUSR 
                                    | S_IRGRP | S_IWGRP 
                                    | S_IROTH | S_IWOTH
                                    );
                            if (md->b.f == -1)
                            {
                                fprintf (stderr, "%s adios_posix_open failed for %s\n", LOGHDR
                                        ,ADIOS_PROVENANCE_OUTPUTFILE
                                        );
                                return 0;
                            }
#ifndef __APPLE__
                            off_t offset = lseek64(md->b.f, 0, SEEK_END);
#else
                            off_t offset = lseek(md->b.f, 0, SEEK_END);
#endif
                            fd->base_offset = offset;
                            fd->pg_start_in_file = offset;
                        }
                        break; 
                    } // case: to file


                    default:
                    {
                        fprintf (stderr, "%s Unknown  mode: %d\n", LOGHDR, md->output_method);
                    } 
                } // switch md->output_method
                break;
            }
            case adios_mode_read:
            {
                fprintf (stderr
                        ,"%s adios_provenance_open: File mode READ is not supported for provenance method.\n"
                        , LOGHDR);
                return 0;
            }
            default:
            {
                fprintf (stderr, "%s Unknown file mode: %d\n", LOGHDR, fd->mode);
                return 0;
            }
        } // switch fd->mode

        // some debugging
        struct adios_method_list_struct * m = fd->group->methods;
        while (m)
        {
            printf("open(): Method %s is used\n", m->method->method);
            m = m->next;
        }

        
    } // (md->rank==0)
    return 1;
}


/** Method called from within adios_group_size(). This gives us the communicator. 
  The group of processes should coordinate here their offsets to the output file in
  case they are writing to a single file (as with the MPI method). In POSIX method, they
  write separate files. 
 */
enum ADIOS_FLAG adios_provenance_should_buffer (struct adios_file_struct * fd
                                               ,struct adios_method_struct * method
                                               )
{
    struct adios_PROVENANCE_data_struct * md = (struct adios_PROVENANCE_data_struct *)
                                                                method->method_data;

    /* If the other (real output) method is POSIX, we must not use the communicator
       provided here. People usually provide the world communicator in the config file
       but the POSIX method does not use it. Each process writes its own file
    */
    if (md->another_method == ADIOS_METHOD_POSIX) 
    {
        md->group_comm = MPI_COMM_NULL;
    }

    fd->group->process_id = md->rank;


    printf("should_buffer(): rank %d, fd->write_size_bytes=%ld group size=%d\n"
                ,md->rank, fd->write_size_bytes, md->size);

    return adios_flag_yes;
}

void adios_provenance_write (struct adios_file_struct * fd
                            ,struct adios_var_struct * v
                            ,void * data
                            ,struct adios_method_struct * method
                            )
{
    // FIXME: What is this code doing and why is it here?
    if (v->got_buffer == adios_flag_yes)
    {
        if (data != v->data)  // if the user didn't give back the same thing
        {
            if (v->free_data == adios_flag_yes)
            {
                free (v->data);
                adios_method_buffer_free (v->data_size);
            }
        }
        else
        {
            // we already saved all of the info, so we're ok.
            return;
        }
    }

}

void adios_provenance_get_write_buffer (struct adios_file_struct * fd
                                       ,struct adios_var_struct * v
                                       ,uint64_t * size
                                       ,void ** buffer
                                       ,struct adios_method_struct * method
                                       )
{
    // FIXME: What is this code doing and why is it here?
    uint64_t mem_allowed;

    if (*size == 0)
    {
        *buffer = 0;

        return;
    }

    if (v->data && v->free_data)
    {
        adios_method_buffer_free (v->data_size);
        free (v->data);
    }

    mem_allowed = adios_method_buffer_alloc (*size);
    if (mem_allowed == *size)
    {
        *buffer = malloc (*size);
        if (!*buffer)
        {
            adios_method_buffer_free (mem_allowed);
            fprintf (stderr, "%s Out of memory allocating %llu bytes for %s\n"
                    ,LOGHDR, *size, v->name
                    );
            v->got_buffer = adios_flag_no;
            v->free_data = adios_flag_no;
            v->data_size = 0;
            v->data = 0;
            *size = 0;
            *buffer = 0;
        }
        else
        {
            v->got_buffer = adios_flag_yes;
            v->free_data = adios_flag_yes;
            v->data_size = mem_allowed;
            v->data = *buffer;
        }
    }
    else
    {
        adios_method_buffer_free (mem_allowed);
        fprintf (stderr, "%s OVERFLOW: Cannot allocate requested buffer of %llu "
                         "bytes for %s\n"
                ,LOGHDR
                ,*size
                ,v->name
                );
        *size = 0;
        *buffer = 0;
    }
}

void adios_provenance_read (struct adios_file_struct * fd
                           ,struct adios_var_struct * v, void * buffer
                           ,uint64_t buffer_size
                           ,struct adios_method_struct * method
                           )
{
    fprintf (stderr 
            ,"%s adios_provenance_read: File mode READ is not supported for provenance method.\n"
            , LOGHDR);
    v->data = buffer;
    v->data_size = buffer_size;
}



/** In the close() function
  - group master collects all indices from the group
  - world master collects all group indices from group masters
  - convert the global index into text format
  - write/send text index
  */
void adios_provenance_close (struct adios_file_struct * fd
                            ,struct adios_method_struct * method
                            )
{
    struct adios_PROVENANCE_data_struct * md = (struct adios_PROVENANCE_data_struct *)
                                                 method->method_data;

    struct adios_index_process_group_struct_v1 * new_pg_root = 0;
    struct adios_index_var_struct_v1 * new_vars_root = 0;
    struct adios_index_attribute_struct_v1 * new_attrs_root = 0;


    switch (fd->mode)
    {
        // In writing mode we gather the indices to one process which then
        // writes/sends it.
        case adios_mode_write:
        case adios_mode_append:
        {
            char * buffer = 0;
            uint64_t buffer_size = 0;
            uint64_t buffer_offset = 0;
            uint64_t index_start = md->b.pg_index_offset;

            // build index
            adios_build_index_v1 (fd, &md->old_pg_root, &md->old_vars_root, &md->old_attrs_root);
            
            /* N->M step: gather the indices from the rest of the group.
               We end up with group masters having the group-wide index.
               If the other method is POSIX, we skip this step.
            */
            if (md->group_comm != MPI_COMM_NULL)
            {

                if (md->rank == 0)
                {
                    int ranks_sent = 1; // assume that we have sent to ourselves
                    buffer_size = 100 * 1024; // try 100k to start
                    buffer = malloc (buffer_size); // group master needs to allocate
                    buffer_offset = 0;
                    int count = 0;

                    // to calculate real offset 
                    //   end of master's is known already: fd->write_size_bytes
                    uint64_t offset_global = fd->write_size_bytes; 
                    uint64_t offset_current;

                    while (ranks_sent < md->size)
                    {
                        // receive from 'ranks_sent' rather than from any 
                        //   to maintain the order of index
                        MPI_Recv (buffer, buffer_size, MPI_BYTE, ranks_sent
                                 ,0, md->group_comm, &md->status
                                 );
                        MPI_Get_count (&md->status, MPI_BYTE, &count);
                        if (buffer_size <= count)
                        {
                            fprintf (stderr, "%s 100k buffer size was too small to receive "
                                             "index of one process.\n", LOGHDR);
                        }
                        else
                        {
                            // get the write_size_bytes value first
                            memcpy( &offset_current, buffer, sizeof(uint64_t));

                            printf("offset received from rank %d: %ld\n", ranks_sent, offset_current);


                            // Parse single index and merge into the global one
                            /* We need to use here the md->b bp_buffer, because
                               its many variables have been initialized.
                               However we need to save and restore its pointers.
                            */
                            char * buffer_save = md->b.buff;
                            uint64_t buffer_size_save = md->b.length;
                            uint64_t offset_save = md->b.offset;

                            md->b.buff = buffer + sizeof(uint64_t);
                            md->b.length = count - sizeof(uint64_t);
                            md->b.offset = 0;


                            adios_parse_process_group_index_v1 (&md->b, &new_pg_root);
                            adios_parse_vars_index_v1 (&md->b, &new_vars_root);
                            adios_parse_attributes_index_v1 (&md->b, &new_attrs_root);

                            /* 
                               We need to fix the offsets here. In MPI method,
                               processes write into a single file. The offsets 
                               are communicated in the should_buffer() function of
                               the MPI method that we did not do in this method.
                               So we need to do it here if this method is used
                               alone. 
                               FIXME: with MPI together, the offsets are already 
                               okay. Why? Who sets them to be correct in this method?
                            */
                            if (md->another_method == ADIOS_METHOD_NULL)
                            {
                                _fix_offsets( new_pg_root, new_vars_root, offset_global);
                            }
                            offset_global += offset_current;


                            // Merge n'th index into group-wide index 
                            adios_merge_index_v1 (&md->old_pg_root
                                                 ,&md->old_vars_root
                                                 ,&md->old_attrs_root
                                                 ,new_pg_root, new_vars_root, new_attrs_root
                                                 );
                            // new_pg_root, new_vars_root, new_attrs_root data struct has been merged
                            // into md->old_*_root, so they must not be freed.
                            new_pg_root = 0;
                            new_vars_root = 0;
                            new_attrs_root = 0;
                            md->b.buff = buffer_save;
                            md->b.length = buffer_size_save;
                            md->b.offset = offset_save;

                        }
                        ranks_sent++;
                    }
                }
                else 
                {
                    // we need to pass the fd->write_size_bytes besides the indices
                    buffer = malloc (sizeof(uint64_t));
                    memcpy( buffer, &(fd->write_size_bytes), sizeof(uint64_t));
                    buffer_offset += sizeof(fd->write_size_bytes);
                    // buffer will be reallocated in the next call to accomodate the index
                    adios_write_index_v1 (&buffer, &buffer_size, &buffer_offset, 0
                                         ,md->old_pg_root
                                         ,md->old_vars_root
                                         ,md->old_attrs_root
                                         );

                    printf("rank %d, offset sent: %ld, fd->base_offset=%ld\n"
                              ,md->rank, fd->write_size_bytes, fd->base_offset);

                    MPI_Send (buffer, buffer_offset, MPI_BYTE, 0, 0
                             ,md->group_comm
                             );
                }
            }
            else  // no communicator: each process wrote its own file
            {
                if (md->rank == 0) 
                    printf ("close(): N->N writing in other method\n");
            }


            // Here we have the groupwide indices at group masters
            // They need to send those indices to the world master
            // along with the file names that can be different for each group

            // get the world master
            int wrank;
            int wsize;
            MPI_Comm_rank (md->world_comm, &wrank);
            MPI_Comm_size (md->world_comm, &wsize);


            if (md->group_comm != md->world_comm ) 
            {
                /* FIXME: Here the group masters should send info to world master */
                /* World master collects the file names and the indices from the
                   group masters. Group sizes have to be communicated as well to 
                   be able to stop the loop.
                 */
                // buffer format:  int group size, int filename length, char[] filename and the indices
                if (wrank == 0) 
                {
                    int ranks_sent = md->size; // assume that we have sent to ourselves
                    buffer_size = 100 * 1024; // try 100k to start
                    if (!buffer) buffer = malloc (buffer_size); // allocate if not allocated earlier
                    buffer_offset = 0;
                    int count = 0;
                    // data to be received additionally to the indices:
                    int gsize; // group size received
                    int len;   // file name length
                    char *fname; // group file name

                    printf ("close(): wrank=%d, starts gathering. buffer size=%llu, buffer=%x\n"
                            ,wrank, buffer_size, buffer);

                    _debug_print_index("master",md->old_pg_root, md->old_vars_root);
                    
                    while (ranks_sent < wsize)
                    {
                        // receive from any because we do not know the group masters 
                        MPI_Recv (buffer, buffer_size, MPI_BYTE, MPI_ANY_SOURCE
                                 ,0, md->world_comm, &md->status
                                 );
                        MPI_Get_count (&md->status, MPI_BYTE, &count);
                        if (buffer_size <= count)
                        {
                            fprintf (stderr, "%s 100k buffer size was too small to receive "
                                             "index from one group.\n", LOGHDR);
                            if (count >= sizeof(int)) 
                            {
                                // we can get the group size at least
                                memcpy( &gsize, buffer, sizeof(int));
                            }
                            else
                            {
                                gsize = wsize; // FIXME: we do not know here how much to skip
                                // this setting will hang the sender processes
                                fprintf (stderr, "%s processes will hang because we do not know "
                                             "at the world master how many is sending data.\n", LOGHDR);
                            }
                        }
                        else
                        {
                            printf ("close(): wrank=%d, received msg, bytes=%d\n", wrank, count);
                            // get the additional values first
                            buffer_offset = 0;
                            memcpy( &gsize, buffer, sizeof(int));
                            buffer_offset += sizeof(int);
                            //printf("close(): wrank=%d, groupd size %d, offset=%ld\n",wrank, gsize, buffer_offset);
                            memcpy( &len, buffer+buffer_offset, sizeof(int));
                            buffer_offset += sizeof(int);
                            //printf("close(): wrank=%d, str len %d, offset=%ld\n",wrank, len, buffer_offset);
                            fname = malloc(len+1);
                            //printf("close(): wrank=%d, fname=%x\n",wrank, fname);
                            memcpy( fname, buffer+buffer_offset, len);
                            buffer_offset += len;
                            fname[len] = '\0'; // terminate string
                            //printf("close(): wrank=%d, filename=%s, offset=%ld\n",wrank, fname, buffer_offset);

                            printf("close(): wrank=%d, info received: group size %d, len %d, filename %s\n"
                                  ,wrank, gsize, len, fname);


                            // Parse single index and merge into the global one
                            /* We need to use here the md->b bp_buffer, because
                               its many variables have been initialized.
                               However we need to save and restore its pointers.
                            */
                            char * buffer_save = md->b.buff;
                            uint64_t buffer_size_save = md->b.length;
                            uint64_t offset_save = md->b.offset;

                            md->b.buff = buffer + buffer_offset;
                            md->b.length = count - buffer_offset;
                            md->b.offset = 0;


                            adios_parse_process_group_index_v1 (&md->b, &new_pg_root);
                            adios_parse_vars_index_v1 (&md->b, &new_vars_root);
                            adios_parse_attributes_index_v1 (&md->b, &new_attrs_root);

                            // debug
                            char name[128];
                            sprintf(name,"group%d",ranks_sent+gsize); 
                            _debug_print_index(name, new_pg_root, new_vars_root);

                            // Merge n'th index into group-wide index 
                            adios_merge_index_v1 (&md->old_pg_root
                                                 ,&md->old_vars_root
                                                 ,&md->old_attrs_root
                                                 ,new_pg_root, new_vars_root, new_attrs_root
                                                 );
                            // new_pg_root, new_vars_root, new_attrs_root data struct has been merged
                            // into md->old_*_root, so they must not be freed.
                            new_pg_root = 0;
                            new_vars_root = 0;
                            new_attrs_root = 0;
                            md->b.buff = buffer_save;
                            md->b.length = buffer_size_save;
                            md->b.offset = offset_save;
                        }

                        ranks_sent += gsize;
                    } //while
                }
                else if (md->rank == 0) // group master
                {
                    // we need to pass group size and filename
                    printf ("close(): wrank=%d, prepare msg\n", wrank);
                    int len = strlen(fd->name);
                    buffer_offset = 0;
                    buffer = malloc (2*sizeof(int)+len);  // 2 ints and a string
                    
                    memcpy( buffer, &(md->size), sizeof(int));
                    buffer_offset += sizeof(int);
                    memcpy( buffer+buffer_offset, &(len), sizeof(int));
                    buffer_offset += sizeof(int);
                    memcpy( buffer+buffer_offset, fd->name, len);
                    buffer_offset += len;

                    // buffer will be reallocated in the next call to accomodate the index
                    adios_write_index_v1 (&buffer, &buffer_size, &buffer_offset, 0
                                         ,md->old_pg_root
                                         ,md->old_vars_root
                                         ,md->old_attrs_root
                                         );

                    //printf("rank %d, offset sent: %ld, fd->base_offset=%ld\n"
                    //          ,md->rank, fd->write_size_bytes, fd->base_offset);

                    printf ("close(): wrank=%d, group size=%d, len=%d, file=%s, msg.len=%ld\n"
                            , wrank, md->size, len, fd->name, buffer_offset);
                    MPI_Send (buffer, buffer_offset, MPI_BYTE, 0, 0, md->world_comm);
                    printf ("close(): wrank=%d sent msg\n", wrank);
                }


                
            }
            else {
                //there is 1 group master only who is the world master
                if (wrank == 0)
                   printf("%s There is only one group and the master is the world master\n", LOGHDR);
            }

            // The world master finally can write/send the index
            if (wrank == 0)
            {
                adios_write_index_v1 (&buffer, &buffer_size, &buffer_offset, index_start
                                     ,md->old_pg_root, md->old_vars_root, md->old_attrs_root
                                     );
                adios_write_version_v1 (&buffer, &buffer_size, &buffer_offset);

                _debug_print_index("global", md->old_pg_root, md->old_vars_root);

                switch(md->output_method)
                {
                    // Writing into a file
                    case (adios_provenance_outputmethod_file):
                    {
#ifndef __APPLE__
                        lseek64 (md->b.f, 0, SEEK_END);
#else
                        lseek (md->b.f, 0, SEEK_END);
#endif
                        write (md->b.f, buffer, buffer_offset);
                        adios_posix_close_internal (&md->b);
                        
                        // some debugging logs
                        // FIXME: do we screw up buffer offset here?
                        /*
                        uint32_t version;
                        adios_parse_version(&md->b, &version);

                        struct adios_index_process_group_struct_v1 * pg_root = 0;
                        struct adios_index_process_group_struct_v1 * pg = 0;
                        struct adios_index_var_struct_v1 * vars_root = 0;

                        adios_parse_index_offsets_v1 (b);
                        

                        printf("%s Process wrank=%d, rank=%d writes provenance to file. \n"
                               "  version=%d, bigendian=%d\n"
                               ,LOGHDR, wrank, rank
                               ,version, md->b.change_endianness
                              );
                        */
                        
                        break;
                    }

                    // Sending through a socket
                    case (adios_provenance_outputmethod_socket):
                    {
                        fprintf (stderr, "%s provenance output method SOCKET is not supported yet.\n", LOGHDR);
                        break;
                    }

                    default:
                    {
                        fprintf (stderr, "%s Unknown provenance output mode: %d\n", LOGHDR, md->output_method);
                    } 
                } // switch (md->output_method)
            } // wrank == 0

            free (buffer);

            adios_clear_index_v1 (new_pg_root, new_vars_root, new_attrs_root);
            new_pg_root = 0;
            new_vars_root = 0;

            break;
        } // cases adios_mode_write and adios_mode_append

        case adios_mode_read:
        {
            fprintf (stderr, "%s  adios_provenance_close(): " 
                    "File mode READ is not supported for provenance method.\n", LOGHDR);
            break;
        }

        default:
        {
            fprintf (stderr, "%s Unknown file mode: %d\n", LOGHDR, fd->mode);
        }
    } // switch (fd->mode)


    // clean up provenance struct's MPI variables
    if (   md->group_comm != MPI_COMM_WORLD
        && md->group_comm != MPI_COMM_SELF
        && md->group_comm != MPI_COMM_NULL
       )
    {
        md->group_comm = MPI_COMM_NULL;
    }

    memset (&md->status, 0, sizeof (MPI_Status));
    md->group_comm = MPI_COMM_NULL;

    // clean up index on
    adios_clear_index_v1 (md->old_pg_root, md->old_vars_root, md->old_attrs_root);
    md->old_pg_root = 0;
    md->old_vars_root = 0;
    md->old_attrs_root = 0;
}

void adios_provenance_finalize (int mype, struct adios_method_struct * method)
{
// nothing to do here
    if (adios_provenance_initialized)
        adios_provenance_initialized = 0;
}

void adios_provenance_end_iteration (struct adios_method_struct * method)
{
}

void adios_provenance_start_calculation (struct adios_method_struct * method)
{
}

void adios_provenance_stop_calculation (struct adios_method_struct * method)
{
}




/* =========== PRIVATE FUNCTIONS ============== */

/** Check what other methods are in use. We need to do different things if
    this method is the only one, or there is MPI or there is POSIX besides.
    Other methods are not yet supported.
 */
static void _check_other_methods(struct adios_file_struct * fd
                                ,struct adios_method_struct * method
                                )
{
    struct adios_PROVENANCE_data_struct * md = (struct adios_PROVENANCE_data_struct *)
                                                            method->method_data;
    struct adios_method_list_struct * m = fd->group->methods;
    while (m)
    {
        printf("should_buffer(): Method %s is used\n", m->method->method);
        switch (m->method->m)
        {
            case (ADIOS_METHOD_POSIX):
            case (ADIOS_METHOD_MPI):
            {
                if (md->another_method != ADIOS_METHOD_NULL)
                    fprintf( stderr, "%s More than one additional methods in use. " 
                            "We consider the last one.\n", LOGHDR);
                md->another_method = m->method->m;
                break;
            }
            case (ADIOS_METHOD_PROVENANCE):
                break;
            default:
            {
                fprintf( stderr, "%s Method %s is not supported with PROVENANCE method.\n" 
                        ,LOGHDR, m->method->method);
            }
        }
        m = m->next;
    }
}


/** Parse string "host:port" and return host as string, port as int.
    If :port part is missing, set port to  default_port  
 */
static void _parse_hoststr (char *str, int default_port, char **host, int *port)
{
    char * colon;
    if (str == NULL) return;
    colon = strchr(str, ':');
    if (colon)
    {
        colon[0] = '\0';       // terminate str at :
        *host = strdup(str);   // host is what is before (in str)
        *port = atoi(colon+1); // port is in colon[1...]
    }
    else
    {
        *host = strdup(str);   // host is the whole string
        *port = default_port;
    }
}



/** Get host:port info from some file.
    Input:
      env_var_name    environment variable name that contains the path to a file to be checked
      filename        name of file to be looked for if not named in the env var.
      path,           first path to look for 'filename'
      alternate_path  alternate path to look for 'filename'
      default_port    value to be returned in 'port' if port info is not found

    Output:
      host            host name, string will be allocated
      port            port value

    Host will be NULL if information is not found
    port will be default_port if that part of info is not found
 */

static void _look_for_listener(char *env_var_name
                              ,char *filename
                              ,char *path
                              ,char *alternate_path
                              ,int   default_port
                              ,char **host
                              ,int  *port
                              )
{
    // first check the environment variable
    char *filepath = NULL;
    if (env_var_name)
    {
        filepath = getenv(env_var_name);
        printf("_look_for_listener(): env %s value is %s\n", env_var_name, filepath);
    }
    
    // check for the paths
    FILE *f = NULL;
    char line[256];
    if (filepath)              // try to open file pointed by the env var
    {
        f = fopen( filepath, "r");
    }
    filepath = malloc(256);
    if (!f && path)           // try to open file in first path option
    {
        snprintf(filepath, 256, "%s/%s", path, filename);
        f = fopen( filepath, "r");
    }
    if (!f && alternate_path) // try to open file in alternate path
    {
        snprintf(filepath, 256, "%s/%s", alternate_path, filename);
        f = fopen( filepath, "r");
    }
    if (f) printf("_look_for_listener(): opened file %s\n", filepath);
    free(filepath);

    // read file content (first line only) if a file found
    if (f)
    {
        if (fgets( line, 256, f))
        {
            _parse_hoststr( line, default_port, host, port);
        }
        fclose(f);
    }
    else // did not find anything
    {
        *host = NULL;
        *port = -1;
    }
} // _look_for_listener


#include "core/adios_socket.h"


/** Create a socket and connect to listener.
    Return socket (int), -1 on failure.
  */
static int _connect_to_listener (char * host, int port)
{
    struct sockaddr_in address;
    int socket;

    if (adios_create_socket(&socket) != 0)
    {
        fprintf(stderr, "%s Cannot create socket!\n", LOGHDR);
        return -1;
    }
    if (adios_set_socket_address( host, port, &address) != 0)
    {
        fprintf(stderr, "%s Cannot set socket address with host=%s port=%d!\n"
               ,LOGHDR, host, port);
        return -1;
    }
    if (adios_connect_socket(socket,&address) != 0)
    {
        fprintf(stderr, "%s Cannot connect to listener: %s!\n", LOGHDR, strerror(errno));
        return -1;
    }
    return socket;

}


/** Adjust the offsets in the group and variable index with the offset argument.
  */
static void _fix_offsets(struct adios_index_process_group_struct_v1 * pg_root 
                        ,struct adios_index_var_struct_v1 * vars_root
                        ,uint64_t offset)
{
    // Traverse the process group (list) and update offset in each item
    while (pg_root)
    {
        pg_root->offset_in_file += offset;
        pg_root = pg_root->next;
    }

    // Traverse the variable index (list) and update offsets in all entries in each item
    uint64_t i;
    while (vars_root)
    {
        for (i=0; i<vars_root->characteristics_count; i++) 
        {
            vars_root->characteristics[i].offset += offset;

        }
        vars_root = vars_root->next;
    }

}

/** used by _debug_print_index function
  */
static const char * value_to_string (enum ADIOS_DATATYPES type, void * data)
{
    static char s [100];
    s [0] = 0;

    switch (type)
    {
        case adios_unsigned_byte:
            sprintf (s, "%u", *(((uint8_t *) data)));
            break;

        case adios_byte:
            sprintf (s, "%d", *(((int8_t *) data)));
            break;

        case adios_short:
            sprintf (s, "%hd", *(((int8_t *) data)));
            break;

        case adios_unsigned_short:
            sprintf (s, "%uh", *(((int8_t *) data)));
            break;

        case adios_integer:
            sprintf (s, "%d", *(((int32_t *) data)));
            break;

        case adios_unsigned_integer:
            sprintf (s, "%u", *(((uint32_t *) data)));
            break;

        case adios_long:
            sprintf (s, "%lld", *(((int64_t *) data)));
            break;

        case adios_unsigned_long:
            sprintf (s, "%llu", *(((uint64_t *) data)));
            break;

        case adios_real:
            sprintf (s, "%f", *(((float *) data)));
            break;

        case adios_double:
            sprintf (s, "%le", *(((double *) data)));
            break;

        case adios_long_double:
            sprintf (s, "%Le", *(((long double *) data)));
            break;

        case adios_string:
            sprintf (s, "%s", ((char *) data));
            break;

        case adios_complex:
            sprintf (s, "(%f %f)", *(((float *) data) + 0)
                                 , *(((float *) data) + 1)
                    );
            break;

        case adios_double_complex:
            sprintf (s, "(%lf %lf)", *(((double *) data) + 0)
                                   , *(((double *) data) + 1)
                    );
            break;
    }

    return s;
}


/** Print the index into a file formatted as bpdump does 
    Output file name is: ./idx.'ext'
  */
static void _debug_print_index(char *ext
                              ,struct adios_index_process_group_struct_v1 * pg_root
                              ,struct adios_index_var_struct_v1 * vars_root
                              )
{

    FILE *f;
    char fname[256];

    snprintf(fname, 256, "idx.%s", ext);
    if ( (f = fopen(fname, "w")) == NULL )
    {
        fprintf(stderr,"%s _debug_print_index(): cannot create file %s\n", LOGHDR, fname);
        f = stderr;
    }
            
    fprintf(f, "============= Process Group Index %s =================\n", ext);
    while (pg_root)
    {
        fprintf (f, "Group: %s\n", pg_root->group_name);
        fprintf (f, "\tProcess ID: %d\n", pg_root->process_id);
        fprintf (f, "\tTime Name: %s\n", pg_root->time_index_name);
        fprintf (f, "\tTime: %d\n", pg_root->time_index);
        fprintf (f, "\tOffset in File: %llu\n", pg_root->offset_in_file);
        pg_root = pg_root->next;
    }

    fprintf(f, "============= Var Index %s =================\n", ext);
    while (vars_root)
    {
        if (!strcmp (vars_root->var_path, "/"))
        {
            fprintf (f, "Var (Group): /%s (%s)\n", vars_root->var_name ,vars_root->group_name);
        }
        else
        {
            fprintf (f, "Var (Group): %s/%s (%s)\n", vars_root->var_path ,vars_root->var_name, vars_root->group_name);
        }
        fprintf (f, "\tDatatype: %s\n", adios_type_to_string_int (vars_root->type));
        fprintf (f, "\tVars Entries: %llu\n", vars_root->characteristics_count);
        uint64_t i;
        for (i = 0; i < vars_root->characteristics_count; i++)
        {
            fprintf (f, "\t\tOffset(%llu)", vars_root->characteristics[i].offset);
        if ((vars_root->characteristics[i].stats) && (vars_root->characteristics[i].stats[0][adios_statistic_min].data))
        {
        fprintf (f, "\t\tMin(%s)", value_to_string (vars_root->type ,vars_root->characteristics[i].stats[0][adios_statistic_min].data));
        }
        if ((vars_root->characteristics[i].stats) && (vars_root->characteristics[i].stats[0][adios_statistic_max].data))
        {
        fprintf (f, "\t\tMin(%s)", value_to_string (vars_root->type ,vars_root->characteristics[i].stats[0][adios_statistic_max].data));
        }
            if (vars_root->characteristics[i].value)
            {
                fprintf (f, "\t\tValue(%s)", value_to_string (vars_root->type ,vars_root->characteristics[i].value));
            }
            if (vars_root->characteristics[i].dims.count != 0)
            {
                int j;

                fprintf (f, "\t\tDims (l:g:o): (");
                for (j = 0; j < vars_root->characteristics[i].dims.count; j++)
                {
                    if (j != 0)
                        fprintf (f, ",");
                    if (vars_root->characteristics[i].dims.dims [j * 3 + 1] != 0)
                    {
                        fprintf (f, "%llu:%llu:%llu"
                                ,vars_root->characteristics[i].dims.dims [j * 3 + 0]
                                ,vars_root->characteristics[i].dims.dims [j * 3 + 1]
                                ,vars_root->characteristics[i].dims.dims [j * 3 + 2]
                               );
                    }
                    else
                    {
                        fprintf (f, "%llu"
                                ,vars_root->characteristics[i].dims.dims [j * 3 + 0]
                               );
                    }
                }
                fprintf (f, ")");
            }
            fprintf (f, "\n");
        }

        vars_root = vars_root->next;
    }


    if (f != stderr) 
        fclose(f);

}
