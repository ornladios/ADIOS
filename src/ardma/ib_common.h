#ifndef ARDMA_COMMON_IB_H
#define ARDMA_COMMON_IB_H


// client sends some data in the connection request to the server
struct ardma_connect_data {
    int rank; // local MPI rank of this client
    int size; // number of clients connecting to the server (all servers together)
    int nc;   // number of clients connecting to this particular server process
    int lrank; // the lowest rank which connects to the same server process 
               // lrank <= rank <= lrank+nc-1
};

#define MAXPATHLEN 255 
#define TEST_CALC_CHECKSUM 1    // calculate checksum of PG data buffers 

enum MSG_TYPE {
    ARDMA_MSG_REQUEST,        // client sends request to server to pull data
    ARDMA_MSG_ACK,            // server sends it after pull completed 
    ARDMA_MSG_UNKNOWN
};

// data sent when clients asks server to pull its data
struct ardma_request_message {
    uint32_t rank;             // client rank

    uint64_t pg_size;          // size of data block (process group buffer)
    uint64_t pg_addr;          // memory address 
    uint32_t pg_rkey;          // rdma key to this address
    uint64_t pg_cksum;         // hash of buffer, calc only if _TEST_CALC_CHECKSUM is defined otherwise=1

    uint64_t idx_size;         // size of data block (index buffer)
    uint64_t idx_addr;         // memory address 
    uint32_t idx_rkey;         // rdma key to this address

    char     path[MAXPATHLEN]; // name of output file
    uint32_t timestep;         // timestep of output
};

struct ardma_ack_message {
    int status; // 0 success, 1 failed
};

// union of all kinds of messages to have a fix size for all messages
struct ardma_message {
    enum MSG_TYPE type;
    union {
        struct ardma_request_message request;
        struct ardma_ack_message     ack;
    } msg;
};


/* PGI compiler has no bswap_64, we define here.
   Source taken from FFMPEG's libavutil/bswap.h, an LGPL co
   http://ffmpeg.org/doxygen/0.5/bswap_8h-source.html
*/
#include <byteswap.h>   // __bswap_64

#ifndef __GNUC__
#  warning __bswap_64 does not exist, define now in ib_common.h
static inline uint64_t __bswap_64(uint64_t x)
{
    union {
        uint64_t ll;
        uint32_t l[2];
    } w, r;
    w.ll = x;
    r.l[0] = bswap_32 (w.l[1]);
    r.l[1] = bswap_32 (w.l[0]);
    return r.ll;
}

#endif /* __GNUC__ */

#include <infiniband/arch.h> // htonll using __bswap_64



#if TEST_CALC_CHECKSUM
#if 0
void hashlittle2(
  const void *key,       /* the key to hash */
  size_t      length,    /* length of the key */
  uint32_t   *pc,        /* IN: primary initval, OUT: primary hash */
  uint32_t   *pb         /* IN: secondary initval, OUT: secondary hash */
  /* for a 64-bit value do something like "*pc + (((uint64_t)*pb)<<32) */
); /* implemented in lookup3.c, public domain software */

static inline uint64_t ardma_calc_checksum (char * buf, uint64_t size)
{
    uint32_t pc, pb;
    size_t length = size & 0xffffffffffffffff;
    fprintf(stderr, "ardma_calc_checksum, buf=%x, size=%llu, length=%lu\n", buf, size, length); 
    hashlittle2(buf, length, &pc, &pb);
    return (uint64_t)pc + (((uint64_t)pb)<<32);
}
#endif /* 0 */

static inline uint64_t ardma_calc_checksum (char * buf, uint64_t size)
{
    unsigned long hash = 5381;
    char* p = buf;
    uint64_t i = 0;

    for (i=0; i<size; i++) {
        hash = ((hash << 5) + hash) + *p; /* hash * 33 + c */
        p++;
    }

    i = (uint64_t) hash;
    return i;
}

#else
static inline uint64_t ardma_calc_checksum (char * buf, uint64_t size) { return 1; }
#endif  /* CALC_CHECKSUM */


#endif /* ARDMA_COMMON_IB_H */
