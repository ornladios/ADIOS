#ifndef _ADIOS_XPMEM_H_
#define _ADIOS_XPMEM_H_

#include <errno.h>
#include <stdio.h>
#include <xpmem.h>
#include <adios_logger.h>


static int lerror;

#define PAGE_SIZE	sysconf(_SC_PAGE_SIZE)

#define SHARE_SIZE(X) (((X) < (PAGE_SIZE) ? (PAGE_SIZE) : (X)))

/* Used to specify size of /tmp/xpmem.share */
#define TMP_SHARE_SIZE	32
#define LOCK_INDEX	TMP_SHARE_SIZE - 1
#define COW_LOCK_INDEX	TMP_SHARE_SIZE - 2

static uint64_t share_size = 10*1024*1024;
static uint64_t index_share_size = 1*1024*1024;

typedef struct _shared_data
{
	uint32_t version;
	uint64_t size;
	uint32_t finalized;
	uint32_t readcount;
	uint32_t offset;
	uint32_t reader;
	uint32_t fd;
	char buffer[1];
}shared_data;

#define ATOMIC_INCREMENT(X) X++;
#define ATOMIC_TEST_INCREMENT(X) if(X==0) X++;
#define ATOMIC_TEST_RESET(X) if(X!=0) X=0;

#define _USE_HOBBES_

#ifndef _USE_HOBBES_
static xpmem_segid_t make_share(char **data, size_t size)
{
	xpmem_segid_t segid;
	int ret;
	
	if(size < PAGE_SIZE)
		size = PAGE_SIZE;
	ret = posix_memalign((void**)data, PAGE_SIZE, size);
	if(ret != 0)
	{
		lerror = errno;
		log_debug(stderr, "error in posix_memalign %d %s\n",
		        lerror, strerror(lerror));		
		return -1;
	}

	memset((void*)*data, 0, size);
	
	segid = xpmem_make(*data, size, XPMEM_PERMIT_MODE, (void *)0600);
	if(segid == -1)
	{
	    lerror = errno;
	    log_debug(stderr, "error in posix_memalign %d %s\n",
		    lerror, strerror(lerror));		
		return -1;
	}
	return segid;
}

static int unmake_share(xpmem_segid_t segid, char *data)
{
	int ret;

	ret = xpmem_remove(segid);
	free(data);
	
	return ret;
}

static char *attach_segid(xpmem_segid_t segid, int size, xpmem_apid_t *apid)
{
	struct xpmem_addr addr;
	char *buff;

	// if(*apid <= 0)
	// {
		*apid = xpmem_get(segid, XPMEM_RDWR, XPMEM_PERMIT_MODE, NULL);
		if(*apid == -1)
		{
			lerror = errno;
			return NULL;
		}
	// }
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

static int write_segid(xpmem_segid_t segid, char *fname)
{
	int fd = open (fname, O_CREAT | O_RDWR | O_TRUNC, 0666);	
	if(fd <= 0)
	{
		perror("open");
		return -1;
	}
	lseek(fd, 0, SEEK_SET);
	log_debug(stderr, "segid = %lli\n", segid);
	write(fd, &segid, sizeof(xpmem_segid_t));
	close(fd);
	return 0;
}

static int read_segid(xpmem_segid_t *segid, char *fname)
{
	int fd = open (fname, O_RDONLY, 0666);
	if(fd < 0)
	{
		perror("open");
		return -1;
	}
	lseek(fd, 0, SEEK_SET);
	read(fd, segid, sizeof(xpmem_segid_t));
	log_debug(stderr, "segid = %lli\n", *segid);
	close(fd);
	return 0;
}

#else
//we use the hobbes specific xpmem functionality from Brian
#define WELL_KNOWN_ADIOS_ID 1

static xpmem_segid_t make_share(char **data, size_t size)
{
	xpmem_segid_t segid = WELL_KNOWN_ADIOS_ID;
	int ret;

	int flags = XPMEM_MEM_MODE | XPMEM_REQUEST_MODE | XPMEM_SIG_MODE; 

	//we will want to send this fd out of this function in some way
	int fd;
	
	
	if(size < PAGE_SIZE)
		size = PAGE_SIZE;
	ret = posix_memalign((void**)data, PAGE_SIZE, size);
	if(ret != 0)
	{
		lerror = errno;
		log_debug(stderr, "error in posix_memalign %d %s\n",
		        lerror, strerror(lerror));		
		return -1;
	}

	memset((void*)*data, 0, size);

	segid = xpmem_make_hobbes(*data, size, XPMEM_PERMIT_MODE, (void*)0600,
	                          flags, segid, &fd);
	
	if(segid == -1)
	{
	    lerror = errno;
	    log_debug(stderr, "error in xpmem make hobbes %d %s\n",
		    lerror, strerror(lerror));		
		return -1;
	}
	if(segid != WELL_KNOWN_ADIOS_ID)
	{
		log_debug(stderr, "unable to get well known segid\n");		
	}

	//store the fd in the segment
	((shared_data*)(*data))->fd = fd;

	log_debug("buffer = %p\tfd = %d\n", *data, ((shared_data*)(*data))->fd);
	
	return segid;
}

static int sleep_on_share(shared_data *d)
{
	fd_set fds;
	unsigned long irqs;
	int status;
	struct timeval tv;

	FD_ZERO(&fds);
	FD_SET(d->fd, &fds);

	tv.tv_sec = 1;
	tv.tv_usec = 0;

	status = select(1, &fds, NULL, NULL, &tv);

	return status;
		
}


static int unmake_share(xpmem_segid_t segid, char *data)
{
	int ret;

	ret = xpmem_remove(segid);
	free(data);
	
	return ret;
}

static char *attach_segid( xpmem_segid_t _segid, int size, xpmem_apid_t *apid)
{
	struct xpmem_addr addr;
	char *buff;
	//we ignore the _segid variable because we are using a
	//well known segid
	xpmem_segid_t segid = WELL_KNOWN_ADIOS_ID;

	// if(*apid <= 0)
	// {
		*apid = xpmem_get(segid, XPMEM_RDWR, XPMEM_PERMIT_MODE, NULL);
		if(*apid == -1)
		{
			lerror = errno;
			return NULL;
		}
	// }
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

//since we are using a well known segid it doesn't really matter for us to
//write the segid to the file
static int write_segid(xpmem_segid_t segid, char *fname)
{
	return 0;
}

static int read_segid(xpmem_segid_t *segid, char *fname)
{
	return 0;
}


#endif


#endif //adios_xpmem_h
