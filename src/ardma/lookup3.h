/* 
    lookup3.c, by Bob Jenkins, May 2006, Public Domain.

*/

#ifndef __ADIOS_RDMA_LOOKUP3_H__
#define __ADIOS_RDMA_LOOKUP3_H__

void hashlittle2(
  const void *key,       /* the key to hash */
  size_t      length,    /* length of the key */
  uint32_t   *pc,        /* IN: primary initval, OUT: primary hash */
  uint32_t   *pb         /* IN: secondary initval, OUT: secondary hash */
  /* for a 64-bit value do something like "*pc + (((uint64_t)*pb)<<32) */
);

#endif
