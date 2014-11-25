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

typedef struct _memory
{
    int size;
    int req_size;
    int offset;
    short childcount;
    short ref;    
    unsigned char* buffer;
    struct _memory *parent;
    struct _memory *prev;
    struct _memory *next;
    struct _memory *child;
    LIST_ENTRY(_memory) entries;
}memory;


static memory* findMemory(IOhandle *h, int req_size, int check_size)
{

	memory *memtemp = NULL;
	memory *current = NULL;
	memory *prov = NULL;
	memory *tb = NULL;
	memory *next = NULL;
    

	int page_size = sysconf(_SC_PAGESIZE);

		
    
	for(memtemp = h->memlist.lh_first; memtemp != NULL;
	    memtemp = memtemp->entries.le_next)
	{

		if(memtemp && (memtemp->size - memtemp->offset) >= req_size)	
		{
			if(!prov || (prov->size - prov->offset) >= (memtemp->size - memtemp->offset))
			{
				prov = memtemp;		
			}
		
		}	
	}

	if(prov)
	{
		//we found a memory region
		//allign req size to 8 bits
		int aligned = req_size + (8-1)&~(8-1);
		//prov is hte parent
		tb = (memory*)malloc(sizeof(memory));
		memset(tb, 0, sizeof(memory));
	
		tb->parent = prov;
		if(prov->child == NULL)
		{
			prov->child = tb;	    
		}
		else
		{
			current = prov->child;
			while(current->next != NULL)
				current = current->next;
			current->next = tb;
			tb->prev = current;	    
			tb->next = NULL;	    
		}
		prov->childcount ++;
		tb->buffer = ptr_from_int64(prov->offset + int64_from_ptr(prov->buffer));
		tb->offset = aligned;
		tb->size = aligned;	
		prov->offset += aligned;

		LIST_INSERT_HEAD(&h->memlist, tb, entries);	

		tb->mr = prov->mr;	
		tb->ref ++;	
		return tb;

	}
	else
	{
		//no memory region suitable - allocate more memory if we can otherwise
		fprintf(stderr, "no suitable buffer found\n  - Allocating more\n");
		return NULL;	
	}
    
	fprintf(stderr, "\n");
    
	return NULL ;
	
}


static void free_func(void *cbd)
{

	memory *tb = (memory*)cbd;
	memory *temp;


	memory *memtemp = NULL;

	if(tb->childcount == 0)
	{
		//no children we can merge to the left
		if(tb->prev != NULL)
		{
			//we have a left!
			temp = tb->prev;
			temp->next = tb->next;
			temp->size += tb->size;
			LIST_REMOVE(tb, entries);
			if(tb->parent != NULL)
			{
				tb->parent->childcount --;	    
			}
			free(tb);	 
			return;	    
		}	
		else
		{
			tb->offset = 0;
			while(tb->next)
			{
				temp = tb->next;
				if(temp->childcount == 0 && temp->ref == 0)
				{
					//we can merge temp into tb also
					tb->next = temp->next;
					tb->size += temp->size;
					LIST_REMOVE(temp, entries);
					if(tb->parent)
					{
						tb->parent->childcount --;
					}
		    
					free(temp);		    
				}
				else 
					break;
		
			}	    
			if(tb->parent && tb->parent->childcount == 1)
			{
				//we have exhausted all the children
				//we merge upwards now
				memory *parent = tb->parent;
				parent->offset = 0;
				parent->child = NULL;
				parent->childcount = 0;
				LIST_REMOVE(tb, entries);
				free(tb);
			}
	    
		}
	
	}
    
    

}



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
