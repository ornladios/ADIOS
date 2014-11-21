/*
 * adios_infocache.c
 *
 *  Created on: Nov 21, 2014
 *      Author: David A. Boyuka II
 */

#include <stdlib.h>

// Utilities
static inline int min(int a, int b) { return a < b ? a : b; }
static inline int max(int a, int b) { return a > b ? a : b; }
#define MALLOC_ARRAY(arr,type,len) { (arr) = (type *)malloc((len) * sizeof(type)); }
#define CALLOC_ARRAY(arr,type,len) { (arr) = (type *)calloc((len), sizeof(type)); }
#define REALLOC_ARRAY(arr,type,len) { (arr) = (type *)realloc((arr), (len) * sizeof(type)); }

#define MALLOC(type, var) type var; MALLOC_ARRAY(var, type, 1);

#define FREE(p) {if (p){free(p); (p)=NULL;}}

#define INITIAL_INFOCACHE_SIZE 16

static void expand_infocache(adios_infocache *cache, int var_capacity) {
    int i;
    const int oldcap = cache->capacity;
    const int newcap = max(max(oldcap * 2, var_capacity), INITIAL_INFOCACHE_SIZE);

    if (oldcap == 0) {
        MALLOC_ARRAY(cache->varinfos, ADIOS_VARINFO, newcap);
        MALLOC_ARRAY(cache->transinfos, ADIOS_TRANSINFO, newcap);
    } else {
        REALLOC_ARRAY(cache->varinfos, ADIOS_VARINFO, newcap);
        REALLOC_ARRAY(cache->transinfos, ADIOS_TRANSINFO, newcap);
    }

    for (i = oldcap; i < newcap; i++) {
        cache->varinfos[i] = NULL;
        cache->transinfos[i] = NULL;
    }

    cache->capacity = newcap;
}

adios_infocache * adios_infocache_new() {
    MALLOC(adios_infocache *, cache);
    cache->capacity = 0;
    cache->varinfos = NULL;
    cache->transinfos = NULL;

    expand_infocache(cache, INITIAL_INFOCACHE_SIZE);
    return cache;
}

void adios_infocache_invalidate(adios_infocache *cache) {
    int i;
    for (i = 0; i < cache->capacity; i++) {
        if (cache->varinfos[i]) {
            if (cache->transinfos[i]) {
                common_read_free_transinfo(cache->varinfos[i], cache->transinfos[i]);
                cache->transinfos[i] = NULL;
            }
            common_read_free_varinfo(cache->varinfos[i]);
            cache->varinfos[i] = NULL;
        }
    }
}

void adios_infocache_free(adios_infocache **cache_ptr) {
    adios_infocache *cache = *cache_ptr;

    adios_infocache_invalidate(cache); // Frees all varinfos/transinfos
    FREE(cache->varinfos);
    FREE(cache->transinfos);
    cache->capacity = 0;
    FREE(*cache_ptr);
}

ADIOS_VARINFO * adios_infocache_inq_varinfo(const ADIOS_FILE *fp, adios_infocache *cache, int varid) {
    if (varid >= cache->capacity)
        expand_infocache(cache, varid);

    if (cache->varinfos[varid])
        return cache->varinfos[varid];
    else
        return cache->varinfos[varid] = common_read_inq_var_raw_byid(fp, varid);
}

ADIOS_TRANSINFO * adios_infocache_inq_transinfo(const ADIOS_FILE *fp, adios_infocache *cache, int varid) {
    if (varid >= cache->capacity)
        expand_infocache(cache, varid);

    if (cache->transinfos[varid]) {
        return cache->transinfos[varid];
    } else {
        ADIOS_VARINFO * vi = adios_infocache_inq_varinfo(fp, cache, varid);
        return cache->transinfos[varid] = common_read_inq_transinfo(fp, vi);
    }
}
