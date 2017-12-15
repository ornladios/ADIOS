#include <fcntl.h>
#include <stdio.h>
#include <sys/stat.h>
#include <unistd.h>
#include <errno.h>

#define MAX_BUF 1024*1024*32


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
    int rc, fd = open(pipe_name, O_RDONLY);

    if (fd < 0)
    {
        printf ("open error.\n");
        return -1;
    }

    do
    {
        rc = read (fd, buf, MAX_BUF);
    } while (rc != 0);

    close(fd);
}

int main()
{
    int fd_field, fd_R, fd_Z, fd_mesh;
    int rc;
    char buf[MAX_BUF];

    remove ("/tmp/MdtmManPipes/field");
    remove ("/tmp/MdtmManPipes/R");
    remove ("/tmp/MdtmManPipes/Z");
    remove ("/tmp/MdtmManPipes/mesh");

    // Create a named pipe for field data
    create_pipe ("/tmp/MdtmManPipes/field");
    create_pipe ("/tmp/MdtmManPipes/R");
    create_pipe ("/tmp/MdtmManPipes/Z");
    create_pipe ("/tmp/MdtmManPipes/mesh");

    read_pipe ("/tmp/MdtmManPipes/field");
    read_pipe ("/tmp/MdtmManPipes/R");
    read_pipe ("/tmp/MdtmManPipes/Z");
    read_pipe ("/tmp/MdtmManPipes/mesh");

    return 0;
}

