#include <fcntl.h>
#include <stdio.h>
#include <sys/stat.h>
#include <unistd.h>
#include <errno.h>
#include <pthread.h>

#define MAX_BUF 1024*32


void create_pipe (char * pipe_name)
{
    // Create a named pipe for field data
    int rc = mkfifo(pipe_name, S_IRWXU | S_IRWXG | S_IRWXO);
    if (rc == -1)
    {
        printf ("mkfifo() error = %s\n", strerror(errno));
    }

}

void read_pipe (char * pipe_name)
{
    char buf[MAX_BUF];
printf ("reading %s\n", pipe_name);
    int rc, fd = open(pipe_name, O_RDONLY);
    if (fd < 0)
    {
        printf ("open error.\n");
        return -1;
    }

//    while (1)
//    {
        if ((rc = read (fd, buf, MAX_BUF)) > 0)
        {
            printf ("rc = %d\n", rc);
        }
        else
        {
 //           break;
        }
  //  }

    close(fd);
}

int main()
{
    int fd_field, fd_R, fd_Z, fd_mesh;
    int rc;
    char buf[MAX_BUF];
    pthread_t thread_field, thread_R, thread_Z, thread_mesh;

    remove ("/tmp/MdtmManPipes/field");
    remove ("/tmp/MdtmManPipes/R");
    remove ("/tmp/MdtmManPipes/Z");
    remove ("/tmp/MdtmManPipes/mesh");

    // Create a named pipe for field data
    create_pipe ("/tmp/MdtmManPipes/field");
    create_pipe ("/tmp/MdtmManPipes/R");
    create_pipe ("/tmp/MdtmManPipes/Z");
    create_pipe ("/tmp/MdtmManPipes/mesh");
#if 0
    read_pipe ("/tmp/MdtmManPipes/field");
    read_pipe ("/tmp/MdtmManPipes/R");
    read_pipe ("/tmp/MdtmManPipes/Z");
    read_pipe ("/tmp/MdtmManPipes/mesh");
#else
    char field_string[] = "/tmp/MdtmManPipes/field";
    char R_string[] = "/tmp/MdtmManPipes/R";
    char Z_string[] = "/tmp/MdtmManPipes/Z";
    char mesh_string[] = "/tmp/MdtmManPipes/mesh";

    pthread_create (&thread_field, NULL, (void *) &read_pipe, (void *) field_string);
    pthread_create (&thread_R, NULL, (void *) &read_pipe, (void *) R_string);
    pthread_create (&thread_Z, NULL, (void *) &read_pipe, (void *) Z_string);
    pthread_create (&thread_mesh, NULL, (void *) &read_pipe, (void *) mesh_string);

    pthread_join (thread_field, NULL);
    pthread_join (thread_R, NULL);
    pthread_join (thread_Z, NULL);
    pthread_join (thread_mesh, NULL);
#endif

    remove ("/tmp/MdtmManPipes/field");
    remove ("/tmp/MdtmManPipes/R");
    remove ("/tmp/MdtmManPipes/Z");
    remove ("/tmp/MdtmManPipes/mesh");

    return 0;
}

