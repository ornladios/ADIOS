#include <stdio.h>
#include <xpmem.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/mman.h>
#include <errno.h>
#include <string.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <mpi.h>

static int lerror;

#define PAGE_SIZE	sysconf(_SC_PAGE_SIZE)

#define SHARE_SIZE(X) (((X) < (PAGE_SIZE) ? (PAGE_SIZE) : (X)))

/* Used to specify size of /tmp/xpmem.share */
#define TMP_SHARE_SIZE	32
#define LOCK_INDEX	TMP_SHARE_SIZE - 1
#define COW_LOCK_INDEX	TMP_SHARE_SIZE - 2

typedef struct _data
{
	uint32_t version;
	uint64_t size;
	uint32_t readcount;
	uint32_t offset;
	uint32_t namelen;
	char name[1];
}data;


xpmem_segid_t make_share(char **data, size_t size)
{
	xpmem_segid_t segid;
	int ret;
	
	if(size < PAGE_SIZE)
		size = PAGE_SIZE;
	ret = posix_memalign((void**)data, PAGE_SIZE, size);
	if(ret != 0)
	{
		lerror = errno;
		fprintf(stderr, "error in posix_memalign %d %s\n",
		        lerror, strerror(lerror));		
		return -1;
	}

	memset((void*)*data, 0, size);
	
	segid = xpmem_make(*data, size, XPMEM_PERMIT_MODE, (void *)0666);
	if(segid == -1)
	{
		fprintf(stderr, "error in xpmem_make size = %u data = %p\n",
		        size, *data);
		return -1;
	}
	return segid;
}

int unmake_share(xpmem_segid_t segid, char *data)
{
	int ret;

	ret = xpmem_remove(segid);
	free(data);
	
	return ret;
}

char *attach_segid(xpmem_segid_t segid, int size, xpmem_apid_t *apid)
{
	struct xpmem_addr addr;
	char *buff;

	if(*apid <= 0)
	{
		*apid = xpmem_get(segid, XPMEM_RDWR, XPMEM_PERMIT_MODE, NULL);
		if(*apid == -1)
		{
			lerror = errno;
			return NULL;
		}
	}
	addr.apid = *apid;
	addr.offset = 0;
	buff = xpmem_attach(addr, size, NULL);
	if(buff == (void*)-1)
	{
		perror("xpmem_attach");
		return NULL;
	}

	return buff;
}

int write_segid(xpmem_segid_t segid)
{
	char *fname = "xpmem.share";
	int fd = open (fname, O_CREAT | O_RDWR | O_TRUNC, 0666);	
	if(fd <= 0)
	{
		perror("open");
		return -1;
	}
	lseek(fd, 0, SEEK_SET);
	write(fd, &segid, sizeof(xpmem_segid_t));
	close(fd);
	return 0;
}

int read_segid(xpmem_segid_t *segid)
{
	char *fname = "xpmem.share";
	int fd = open (fname, O_RDONLY, 0666);
	if(fd < 0)
	{
		perror("open");
		return -1;
	}
	lseek(fd, 0, SEEK_SET);
	read(fd, segid, sizeof(xpmem_segid_t));	
	printf("segid = %llx\n", *segid);
	close(fd);
	return 0;
}



int main(int argc, char **argv)
{
	int version;
	char *buffer = NULL;
	int size = 1024;   
	xpmem_segid_t xseg;
	int i = 0;
	xpmem_apid_t apid = 0;
	int         rank, rsize;
	data *d;
	
	// version = xpmem_version();
	// printf("version = %d page_size = %d\n", version, PAGE_SIZE);
	MPI_Init (&argc, &argv);
	MPI_Comm_rank (MPI_COMM_WORLD, &rank);
	MPI_Comm_size (MPI_COMM_WORLD, &rsize);
	
	if(rank == 0)
	{
		//this is a server
		xseg = make_share(&buffer, size);

		if(xseg == -1)
		{
			printf("error in creating segment\n");
			return -1;
		}
			
		write_segid(xseg);
		printf("segid = %llx\n", xseg);

		d = (data*)buffer;
		d->version = 0;
		d->size = 100;
		strcpy(d->name, "test_group");
		d->namelen = strlen(d->name);
		d->offset = (int)d->name - (int)&d->version + d->namelen;

		printf("version = %d, size = %d, namelen = %d offset = %d\n",
		       d->version, d->size, d->namelen, d->offset);

		sprintf(&buffer[d->offset], "This is a string I am putting in\n");

		d->version = 1;

		while(d->readcount != 1)
			usleep(10000);

		sprintf(&buffer[d->offset], "This is version 2\n");

		d->readcount = 0;
		d->version ++;


		while(d->readcount != 1)
			usleep(10000);
		
		unmake_share(xseg, buffer);
	}
	else
	{
		//this is a client
		do
		{
			read_segid(&xseg);
			printf("segid = %llx\n", xseg);
			d = (data*)attach_segid(xseg, sizeof(data), &apid);
		}while(d == NULL);

		printf("apid = %d, buffer = %p\n", apid, d);


		while(d->version != 1)
		{
			usleep(10000);
		}

		printf("version = %d, size = %d, namelen = %d offset = %d\n",
		       d->version, d->size, d->namelen, d->offset);

		size = d->size;
		
		xpmem_detach(d);
		
		buffer = attach_segid(xseg, size, &apid);
		
		if(buffer == NULL)
		{
			printf("unable to attach twice");
			return -1;
		}

		d= (data*)buffer;
		printf("apid = %d, buffer = %p\n", apid, buffer);

		
		printf("%s\n", &buffer[d->offset]);
		// for(i = 0; i < size; i ++)
		// {
		// 	if(buffer[i] == 0)
		// 	{
		// 		snprintf(&buffer[i], size - i, "\nNew string from rank %d\n", rank);
		// 		break;
		// 	}
		// 	putc(buffer[i], stdout);
			
		// }

		d->readcount ++;

		while(d->version != 2)
			usleep(10000);

		
		printf("%s\n", &buffer[d->offset]);

		d->readcount ++;
		
		xpmem_detach(buffer);
		xpmem_release(apid);
	}

	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize();
	return 0;
}
