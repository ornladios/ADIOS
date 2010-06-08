/*
 * aggregation.cpp
 *
 *  Created on: Mar 9, 2009
 *      Author: thkorde
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <assert.h>

#include <nssi_server.h>

#include <algorithm>
#include <map>
#include <list>

using namespace std;

#include "aggregation.h"


struct ltstr
{
  bool operator()(const char* s1, const char* s2) const
  {
    return strcmp(s1, s2) < 0;
  }
};



typedef struct {
    int ahead_count;
    int behind_count;
    int same_count;
    int no_match_count;
} chunk_location_count_t;

typedef list<aggregation_chunk_details_t *> chunk_details_t;
typedef list<aggregation_chunk_details_t *>::iterator chunk_details_iterator_t;

typedef struct {
    aggregation_chunk_details_t *details;
    chunk_details_t              component_chunks;
} aggregation_chunk_t;

typedef list<aggregation_chunk_t *> chunks_t;
typedef list<aggregation_chunk_t *>::iterator chunks_iterator_t;


typedef struct {
    char      var_path[ADIOS_PATH_MAX];
    char      var_name[ADIOS_PATH_MAX];
    chunks_t *chunks;
    int       chunks_received;
} per_var_details_t;

typedef map<const char *, per_var_details_t *, ltstr> var_map_t;
typedef map<const char *, per_var_details_t *, ltstr>::iterator var_map_iterator_t;
typedef pair<const char *, per_var_details_t *> var_map_pair_t;

typedef struct {
    int        fd;
    var_map_t  vars;
    write_type type; /* direct, aggregate independent, aggregate collective */
} file_details_t;

static map<int, file_details_t *> open_file_map;
typedef map<int, file_details_t *>::iterator open_file_map_iterator_t;
typedef pair<int, file_details_t *> open_file_map_pair_t;



static int DEBUG=0;




bool compare_chunks_for_aggregation(const aggregation_chunk_t* c1, const aggregation_chunk_t* c2)
{
    aggregation_chunk_details_t *details1=c1->details;
    aggregation_chunk_details_t *details2=c2->details;

    for (int i=0;i<details1->ndims;i++) {
        if (details1->count[i] < details2->count[i]) {
            return true;
        } else if (details1->count[i] > details2->count[i]) {
            return false;
        }
    }
    for (int i=0;i<details1->ndims;i++) {
        if (details1->offset[i] < details2->offset[i]) {
            return true;
        } else if (details1->offset[i] > details2->offset[i]) {
            return false;
        }
    }
    return false;
}

bool compare_chunks_for_caching(const aggregation_chunk_t* c1, const aggregation_chunk_t* c2)
{
    aggregation_chunk_details_t *details1=c1->details;
    aggregation_chunk_details_t *details2=c2->details;

    for (int i=0;i<details1->ndims;i++) {
        if (details1->offset[i] < details2->offset[i]) {
            return true;
        } else if (details1->offset[i] > details2->offset[i]) {
            return false;
        }
    }
    return false;
}

file_details_t *new_open_file(const int fd)
{
    file_details_t *details=NULL;

    details = new file_details_t;
    details->fd = fd;

    return details;
}

int use_aggregation(const int fd)
{
    file_details_t *details=NULL;

    details = open_file_map[fd];
    if ((details->type == WRITE_AGGREGATE_INDEPENDENT) ||
        (details->type == WRITE_AGGREGATE_COLLECTIVE)) {
        return TRUE;
    }

    return FALSE;
}

int use_caching(const int fd)
{
    file_details_t *details=NULL;

    details = open_file_map[fd];
    if ((details->type == WRITE_CACHING_INDEPENDENT) ||
        (details->type == WRITE_CACHING_COLLECTIVE)) {
        return TRUE;
    }

    return FALSE;
}

int use_collective(const int fd)
{
    file_details_t *details=NULL;

    details = open_file_map[fd];
    if ((details->type == WRITE_AGGREGATE_COLLECTIVE) ||
        (details->type == WRITE_CACHING_COLLECTIVE)) {
        return TRUE;
    }

    return FALSE;
}

int use_independent(const int fd)
{
    file_details_t *details=NULL;

    details = open_file_map[fd];
    if ((details->type == WRITE_AGGREGATE_INDEPENDENT) ||
        (details->type == WRITE_CACHING_INDEPENDENT)) {
        return TRUE;
    }

    return FALSE;
}

int use_direct(const int fd)
{
    file_details_t *details=NULL;

    details = open_file_map[fd];
    if (details->type == WRITE_DIRECT) {
        return TRUE;
    }

    return FALSE;
}

int getTypeSize(
        enum ADIOS_DATATYPES type,
        void *val)
{
    switch (type)
    {
    case adios_byte:
    case adios_unsigned_byte:
        return 1;

    case adios_string:
        return strlen ((char *) val);

    case adios_short:
    case adios_unsigned_short:
        return 2;

    case adios_integer:
    case adios_unsigned_integer:
        return 4;

    case adios_real:
        return 4;

    case adios_long:
    case adios_unsigned_long:
        return 8;

    case adios_double:
        return 8;

    case adios_long_double:
        return 16;

    case adios_complex:
        return 2 * 4;

    case adios_double_complex:
        return 2 * 8;

    default:
        return -1;
    }
}

void add_file(const int fd,
            const write_type write_type)
{
    file_details_t *details=NULL;

    details = open_file_map[fd];
    if (details == NULL) {
        details=new_open_file(fd);
        open_file_map[fd]=details;
    }

    details->type = write_type;
}

void add_chunk(aggregation_chunk_details_t *chunk_details)
{
    file_details_t  *file_details=NULL;

    if (DEBUG > 3) printf("adding chunk: fd(%d) var_name(%s)\n", chunk_details->fd, chunk_details->var_name);

    chunk_details->atype_size = getTypeSize(chunk_details->atype, chunk_details->buf);
    if ((chunk_details->len > 0) && (chunk_details->num_elements > 0)) {
        if ((chunk_details->len/chunk_details->num_elements) != chunk_details->atype_size) {
            printf("datatype size conflict: (%lu/%d)==%lu is not equal to %d\n",
                    chunk_details->len, chunk_details->num_elements, chunk_details->len/chunk_details->num_elements, chunk_details->atype_size);
            print_chunk(chunk_details);
        }
    }

    file_details = open_file_map[chunk_details->fd];
    if (file_details == NULL) {
        printf("failed to add chunk.  cannot aggregate.\n");
        return;
    }
    per_var_details_t *var_details = file_details->vars[chunk_details->var_name];
    if (var_details == NULL) {
//        if (DEBUG > 3) printf("var_details don't exist for %s\n", chunk_details->var_name);
        var_details=new per_var_details_t;
        strcpy(var_details->var_path, chunk_details->var_path);
        strcpy(var_details->var_name, chunk_details->var_name);
        var_details->chunks = new chunks_t;
        var_details->chunks_received=0;
        file_details->vars[chunk_details->var_name]=var_details;
    } else {
//        if (DEBUG > 3) printf("var_details already exist for %s\n", chunk_details->var_name);
    }
    aggregation_chunk_t *chunk=new aggregation_chunk_t;
    chunk->details = chunk_details;
    var_details->chunks->push_back(chunk);
    var_details->chunks_received++;

    return;
}

void destroy_chunk(aggregation_chunk_details_t *details)
{
    free(details->offset);
    free(details->count);
    for (int i=0;i<details->ndims;i++) {
        free(details->offset_path[i]);
        free(details->offset_name[i]);
        free(details->count_path[i]);
        free(details->count_name[i]);
    }
    free(details->offset_path);
    free(details->offset_name);
    free(details->count_path);
    free(details->count_name);
//    if (DEBUG > 3) printf("freeing details->buf(%p)\n", details->buf);
    free(details->buf);
    delete details;
}

void cleanup_aggregation_chunks(const int fd)
{
    file_details_t  *details=NULL;
    var_map_iterator_t var_iter;
    per_var_details_t *var_details=NULL;
    aggregation_chunk_t *chunk=NULL;
    chunks_iterator_t chunks_iter;
    chunk_details_iterator_t component_iter;

//    if (DEBUG > 3) printf("entered cleanup_aggregation_chunks\n");
//    if (DEBUG > 3) printf("cleaning up - fd(%d)\n", fd);

    details = open_file_map[fd];
    if (details == NULL) {
        return;
    }
//    var_iter = details->vars.begin();
//    for (; var_iter != details->vars.end(); ++var_iter) {
//        var_details = var_iter->second;
//        if (var_details != NULL) {
//            if (DEBUG > 3) printf("var_details first(%p) second(%s)\n", var_iter->first, var_details->var_name);
//        } else {
//            if (DEBUG > 3) printf("var_details is NULL\n");
//        }
//    }
    var_iter = details->vars.begin();
    for (; var_iter != details->vars.end();) {
        var_details = var_iter->second;
        if (var_details != NULL) {
//            cleanup_aggregation_chunks(fd, var_details->var_name);
            chunks_iter = var_details->chunks->begin();
            for (;chunks_iter != var_details->chunks->end(); ++chunks_iter) {
                chunk = *chunks_iter;
                component_iter = chunk->component_chunks.begin();
                for (;component_iter != chunk->component_chunks.end(); ++component_iter) {
//                    if (DEBUG > 3) printf("cleanup - destroying component\n");
                    destroy_chunk(*component_iter);
                }
                chunk->component_chunks.clear();
//                if (DEBUG > 3) printf("cleanup - destroying details\n");
                destroy_chunk(chunk->details);
                delete chunk;
            }
            var_details->chunks->clear();
            var_details->chunks_received=0;
            delete var_details;
        } else {
//            if (DEBUG > 3) printf("cannot cleanup - var_details is NULL\n");
        }
        details->vars.erase(var_iter++);
    }
//    details->vars.clear();

    for(var_map_iterator_t vars_iter=details->vars.begin(); vars_iter!=details->vars.end(); ++vars_iter) {
        per_var_details_t *pvd=vars_iter->second;
        if (pvd != NULL) {
            if (DEBUG > 3) printf("var_details first(%p) second(%s)\n", vars_iter->first, vars_iter->second->var_name);
        } else {
            if (DEBUG > 3) printf("var_details is NULL\n");
        }
    }
}

void cleanup_aggregation_chunks(const int fd, const char *var_name)
{
    file_details_t  *details=NULL;
    aggregation_chunk_t *chunk=NULL;
    chunks_iterator_t chunks_iter;
    chunk_details_iterator_t component_iter;
    var_map_iterator_t vars_iter;

//    if (DEBUG > 3) printf("cleaning up - fd(%d) var_name(%s)\n", fd, var_name);

    // for each variable, iterate over the chunks and destroy them

    details = open_file_map[fd];

    per_var_details_t *var_details = details->vars[var_name];
    if (var_details != NULL) {
        chunks_iter = var_details->chunks->begin();
        for (;chunks_iter != var_details->chunks->end(); ++chunks_iter) {
            chunk = *chunks_iter;
            component_iter = chunk->component_chunks.begin();
            for (;component_iter != chunk->component_chunks.end(); ++component_iter) {
//                if (DEBUG > 3) printf("cleanup - destroying component\n");
                destroy_chunk(*component_iter);
            }
            chunk->component_chunks.clear();
//            if (DEBUG > 3) printf("cleanup - destroying details\n");
            destroy_chunk(chunk->details);
            delete chunk;
        }
        var_details->chunks->clear();
        var_details->chunks_received=0;
        delete var_details;
    } else {
//        if (DEBUG > 3) printf("cleanup failed - var_details is NULL (%s)\n", var_name);
    }

    var_map_iterator_t iter=details->vars.find(var_name);
    if (iter != details->vars.end()) {
//        if (DEBUG > 3) printf("erasing var_details with iter\n");
        details->vars.erase(iter);
    } else {
//        if (DEBUG > 3) printf("cannot erase var_details with iter.  var_details not found.\n");
    }

}

static void recursive_print_chunk(aggregation_chunk_details_t *details, int offset, int *index, int current_dim)
{
    int my_offset=0;
    char tmp_str[20];
    char out_str[1024];
    int remaining=1023;

    if (current_dim < details->ndims-1) {
        for (int i=0;i<details->count[current_dim];i++) {
            my_offset = index[current_dim];
            for (int i=current_dim+1;i<details->ndims;i++) {
                my_offset *= details->count[i];
            }

            index[current_dim+1]=0;
            recursive_print_chunk(details, offset+my_offset, index, current_dim+1);
            index[current_dim] += details->atype_size;
        }
        //if (DEBUG > 3) printf("-----------------------------\n");
    } else {
        if (details->buf == NULL) {
            if (DEBUG > 3) printf("details->buf == NULL\n");
        } else {
            out_str[0]='\0';
            for (int i=0;i<details->count[current_dim];i++) {
                my_offset = offset+index[current_dim];

//                if (i==0) if (DEBUG > 3) printf("[%d][%d][%d] (my_offset==%d)\n", index[0], index[1], index[2], my_offset);
                if ((details->atype == adios_byte) || (details->atype == adios_unsigned_byte)) {
                    sprintf(tmp_str, "%c, ", *(char *)(((char *)details->buf) + my_offset));
                }
                else if (details->atype == adios_short || details->atype == adios_unsigned_short) {
                    sprintf(tmp_str, "%hx, ", *(short *)(((char *)details->buf) + my_offset));
                }
                else if (details->atype == adios_integer || details->atype == adios_unsigned_integer) {
                    sprintf(tmp_str, "%x, ", *(int *)(((char *)details->buf) + my_offset));
                }
                else if (details->atype == adios_long || details->atype == adios_unsigned_long) {
                    sprintf(tmp_str, "%lx, ", *(int *)(((char *)details->buf) + my_offset));
                }
                else if (details->atype == adios_real) {
                    sprintf(tmp_str, "%f, ", *(float *)(((char *)details->buf) + my_offset));
                }
                else if (details->atype == adios_double) {
                    sprintf(tmp_str, "%f, ", *(double *)(((char *)details->buf) + my_offset));
                }
                strncat(out_str, tmp_str, remaining);
                remaining -= strlen(out_str);

                index[current_dim] += details->atype_size;
            }
//            if (DEBUG > 3) printf("[%d][%d][%d] (my_offset==%d)\n", index[0], index[1], index[2], my_offset);
            if (DEBUG > 3) printf("%s\n", out_str);
        }
    }
}

void print_chunk(aggregation_chunk_details_t *details)
{
    int *index=(int *)calloc(details->ndims, sizeof(int));
    char tmp_str[20];
    char out_str[1024];
    int remaining=1023;

    if (DEBUG > 3) printf("+++++++++++++++++++++++++++++\n");

    if (DEBUG > 3) printf("fd==%d\n", details->fd);
    if (DEBUG > 3) printf("var_path==%s\n", details->var_path);
    if (DEBUG > 3) printf("var_name==%s\n", details->var_name);
    if (DEBUG > 3) printf("ndims==%d\n", details->ndims);
    if (DEBUG > 3) printf("len==%ld\n", details->len);
    if (DEBUG > 3) printf("num_elements==%d\n", details->num_elements);
    out_str[0]='\0';
    remaining=1023;
    for (int i=0;(i<details->ndims) && (remaining>0);i++) {
        sprintf(tmp_str, "%ld,", details->offset[i]);
        strncat(out_str, tmp_str, remaining);
        remaining -= strlen(tmp_str);
    }
    if (DEBUG > 3) printf("offset[]==%s\n", out_str);
    out_str[0]='\0';
    remaining=1023;
    for (int i=0;(i<details->ndims) && (remaining>0);i++) {
        sprintf(tmp_str, "%ld,", details->count[i]);
        strncat(out_str, tmp_str, remaining);
        remaining -= strlen(tmp_str);
    }
    if (DEBUG > 3) printf("count[]==%s\n", out_str);
    if (DEBUG > 3) printf("buf==%p\n", details->buf);

//    int offset=0;
//    int current_dim=0;
//    recursive_print_chunk(details, offset, index, current_dim);
    if (DEBUG > 3) printf("+++++++++++++++++++++++++++++\n");

    free(index);
}

void print_chunk(aggregation_chunk_t *c)
{
    if (c->details == NULL) {
        if (DEBUG > 3) printf("chunk has no details.  perhaps it was aggregated into another chunk.\n");
        return;
    }
    print_chunk(c->details);
}

static void recursive_copy_chunk(aggregation_chunk_details_t *src,
                                 aggregation_chunk_details_t *dst,
                                 int src_offset,
                                 int dst_offset,
                                 int *src_index,
                                 int *dst_index,
                                 int current_dim)
{
    int my_src_offset=0;
    int my_dst_offset=0;

    if (current_dim < src->ndims-1) {
        for (int i=0;i<src->count[current_dim];i++) {
            my_src_offset = src_index[current_dim];
            my_dst_offset = dst_index[current_dim];
//            if (DEBUG > 3) printf("join_offset(%d) offset_diff[%d](%d)\n",
//                    join_offset, current_dim, src->offset[current_dim] - dst->offset[current_dim]);
            my_dst_offset += ((src->offset[current_dim] - dst->offset[current_dim]) * src->atype_size);
            for (int j=current_dim+1;j<src->ndims;j++) {
                my_src_offset *= src->count[j];
                my_dst_offset *= dst->count[j];
            }

            src_index[current_dim+1]=0;
            dst_index[current_dim+1]=0;
            recursive_copy_chunk(src, dst, src_offset+my_src_offset, dst_offset+my_dst_offset,
                                 src_index, dst_index, current_dim+1);
            src_index[current_dim] += src->atype_size;
            dst_index[current_dim] += dst->atype_size;
        }
    } else {
        dst_offset += ((src->offset[current_dim] - dst->offset[current_dim]) * src->atype_size);
        memcpy(((char *)dst->buf) + dst_offset,
               ((char *)src->buf) + src_offset,
               src->count[current_dim]*src->atype_size);
    }
}

static void recursive_aggregate_chunks(aggregation_chunk_t *src1,
                                       aggregation_chunk_t *src2,
                                       aggregation_chunk_t *dst)
{
    int *src_index=(int *)calloc(src1->details->ndims, sizeof(int));
    int *dst_index=(int *)calloc(dst->details->ndims, sizeof(int));
    int src_offset=0;
    int dst_offset=0;
    int current_dim=0;

    memset(src_index, 0, src1->details->ndims*sizeof(int));
    memset(dst_index, 0, dst->details->ndims*sizeof(int));
    recursive_copy_chunk(src1->details, dst->details, src_offset, dst_offset, src_index, dst_index, current_dim);
    memset(src_index, 0, src2->details->ndims*sizeof(int));
    memset(dst_index, 0, dst->details->ndims*sizeof(int));
    recursive_copy_chunk(src2->details, dst->details, src_offset, dst_offset, src_index, dst_index, current_dim);

    free(src_index);
    free(dst_index);
}

static void copy_chunk(aggregation_chunk_details_t *src,
                       aggregation_chunk_details_t *dst)
{
    int *src_index=(int *)calloc(src->ndims, sizeof(int));
    int *dst_index=(int *)calloc(dst->ndims, sizeof(int));
    int src_offset=0;
    int dst_offset=0;
    int current_dim=0;

    memset(src_index, 0, src->ndims*sizeof(int));
    memset(dst_index, 0, dst->ndims*sizeof(int));
    recursive_copy_chunk(src, dst, src_offset, dst_offset, src_index, dst_index, current_dim);

    free(src_index);
    free(dst_index);
}

aggregation_chunk_t *aggregate_chunks(aggregation_chunk_t *c1,
                                      aggregation_chunk_t *c2,
                                      int join_dim)
{
    aggregation_chunk_t *out=new aggregation_chunk_t;

    //if (DEBUG > 3) printf("entered aggregate_chunks\n");

    assert(c1->details->ndims == c2->details->ndims);
    assert(out != NULL);

    out->details = new aggregation_chunk_details_t;

    out->details->fd           = c1->details->fd;
    strcpy(out->details->var_path, c1->details->var_path);
    strcpy(out->details->var_name, c1->details->var_name);
    out->details->ndims        = c1->details->ndims;
//    out->details->buf          = calloc(c1->details->len+c2->details->len, c1->details->atype_size);
    out->details->buf          = NULL;
    out->details->atype        = c1->details->atype;
    out->details->len          = c1->details->len+c2->details->len;
    out->details->atype     = c1->details->atype;
    out->details->num_elements = c1->details->num_elements+c2->details->num_elements;
    out->details->atype_size   = c1->details->atype_size;
    out->details->offset_path  = (char **)calloc(c1->details->ndims, sizeof(char *));
    out->details->offset_name  = (char **)calloc(c1->details->ndims, sizeof(char *));
    out->details->offset       = (uint64_t *)calloc(c1->details->ndims, sizeof(uint64_t));
    out->details->count_path   = (char **)calloc(c1->details->ndims, sizeof(char *));
    out->details->count_name   = (char **)calloc(c1->details->ndims, sizeof(char *));
    out->details->count        = (uint64_t *)calloc(c1->details->ndims, sizeof(uint64_t));

    for (int i=0;i<c1->details->ndims;i++) {
        out->details->offset_path[i]  = strdup(c1->details->offset_path[i]);
        out->details->offset_name[i]  = strdup(c1->details->offset_name[i]);
        out->details->count_path[i]   = strdup(c1->details->count_path[i]);
        out->details->count_name[i]   = strdup(c1->details->count_name[i]);
    }
    memcpy(out->details->offset, c1->details->offset, c1->details->ndims*sizeof(uint64_t));
    memcpy(out->details->count, c1->details->count, c1->details->ndims*sizeof(uint64_t));
    out->details->count[join_dim] += c2->details->count[join_dim];

//    recursive_aggregate_chunks(c1, c2, out);

    if (c1->component_chunks.size() > 0) {
        out->component_chunks.merge(c1->component_chunks);
        c1->component_chunks.clear();
        destroy_chunk(c1->details);
    } else {
        out->component_chunks.push_back(c1->details);
    }
    c1->details = NULL;
    if (c2->component_chunks.size() > 0) {
        out->component_chunks.merge(c2->component_chunks);
        c2->component_chunks.clear();
        destroy_chunk(c2->details);
    } else {
        out->component_chunks.push_back(c2->details);
    }
    c2->details = NULL;

    assert(out != NULL);

    //if (DEBUG > 3) printf("finished\n");

    return(out);
}

/*
 * Aggregate a particular variable in the file.
 *
 * Aggregation rules:
 *  - dimension count must be equal
 *  - strides must be equal
 *  - counts on matching faces must be equal
 *  -
 *
 */
int try_aggregation(const int fd, const char *var_name)
{
    int aggregation_success=FALSE;

    file_details_t  *file_details=NULL;
    per_var_details_t *var_details=NULL;
    aggregation_chunk_t *base_chunk=NULL;
    aggregation_chunk_t *candidate_chunk=NULL;
    aggregation_chunk_t *new_chunk=NULL;
    chunks_iterator_t base_iter, candidate_iter;
    int *offset_diff;
    chunk_location_count_t chunk_location_count;
    int dim_with_movement=-1;

    chunks_t agg_chunks;

    int failed=FALSE;

    file_details = open_file_map[fd];
    if (file_details == NULL) {
//        if (DEBUG > 3) printf("agg failed for %s: file_details==NULL\n", var_name);
        return(aggregation_success);
    }
    var_details = file_details->vars[var_name];
    if (var_details == NULL) {
//        if (DEBUG > 3) printf("agg failed for %s: var_details==NULL\n", var_name);
        return(aggregation_success);
    }
    if (var_details->chunks->size() < 2) {
//        if (DEBUG > 3) printf("returning with chunk count(%d)\n", var_details->chunks->size());
        return(aggregation_success);
    }
    if (DEBUG > 3) printf("chunk count(%d)\n", var_details->chunks->size());


    if (DEBUG > 3) printf("trying aggregation - fd(%d) var_name(%s)\n", fd, var_name);

    var_details->chunks->sort(compare_chunks_for_aggregation);

    if (DEBUG > 4) {
        printf("*****************\n");
        printf("start aggregation (begin list)\n");
        printf("*****************\n");
        int chunk_count;
        aggregation_chunk_details_t **chunks = get_chunks(fd, var_name, &chunk_count);
        for (int i=0;i<chunk_count;i++) {
            print_chunk(chunks[i]);
        }
        free(chunks);
        printf("*****************\n");
        printf("start aggregation (end list)\n");
        printf("*****************\n");
    }

    int success_this_pass=TRUE;
    while (success_this_pass==TRUE) {
        success_this_pass=FALSE;

//        if (DEBUG > 3) printf("top: while loop\n");


        base_iter = var_details->chunks->begin();
        base_chunk = *base_iter;
        offset_diff=new int[base_chunk->details->ndims];
        for (;base_iter != var_details->chunks->end(); ++base_iter) {
//            if (DEBUG > 3) printf("top: base_iter loop\n");

            base_chunk = *base_iter;

            //if (base_chunk != NULL)      print_chunk(base_chunk);

            // look for a chunk that can be aggregated to the base chunk
            candidate_iter = base_iter;
            candidate_iter++;
            for (;candidate_iter != var_details->chunks->end(); ++candidate_iter) {
//                if (DEBUG > 3) printf("top: candidate_iter loop\n");

                candidate_chunk = *candidate_iter;

                //if (candidate_chunk != NULL) print_chunk(candidate_chunk);

                failed=FALSE;

                if (base_chunk->details->ndims != candidate_chunk->details->ndims) {
                    continue;
                }
//                if (candidate_chunk->details->offset[0] != base_chunk->details->offset[0]) {
//                    continue;
//                }
                for (int i=0; i<base_chunk->details->ndims; i++) {
                    offset_diff[i] = candidate_chunk->details->offset[i] - base_chunk->details->offset[i];
                }
                if (failed) continue;

                chunk_location_count.ahead_count=0;
                chunk_location_count.behind_count=0;
                chunk_location_count.same_count=0;
                chunk_location_count.no_match_count=0;
                int agg_dims=base_chunk->details->ndims; /* the number of dimensions to aggregate */
                for (int i=0; i<agg_dims; i++) {
                    if ((offset_diff[i] < 0) && (-offset_diff[i] == candidate_chunk->details->count[i])) {
                        // the candidate is "behind/below" and touching the base chunk in this dimension
                        chunk_location_count.behind_count++;
                        dim_with_movement=i;
                    } else if ((offset_diff[i] > 0) && (offset_diff[i] == base_chunk->details->count[i])) {
                        // the candidate is "ahead of/above" and touching the base chunk in this dimension
                        chunk_location_count.ahead_count++;
                        dim_with_movement=i;
                    } else if (offset_diff[i] == 0) {
                        // the candidate is "equal to" the base chunk in this dimension
                        chunk_location_count.same_count++;
                    } else {
                        // the candidate and the base chunk don't match in this dimension
                        chunk_location_count.no_match_count++;
                    }
                }

#ifdef DEBUG
                /*
                 * These tests can be interesting, but are not required to get the job done.
                 */
                if (chunk_location_count.no_match_count > 0) {
                    // no matching face found.  can't aggregate.
                    continue;
                }

                if (chunk_location_count.same_count == base_chunk->ndims) {
                    // base and candidate have same offset.  bad?  can't aggregate.
                    continue;
                }

                if (chunk_location_count.ahead_count > 1) {
                    // movement in more than one direction
                    continue;
                }
                if (chunk_location_count.behind_count > 1) {
                    // movement in more than one direction
                    continue;
                }
                if ((chunk_location_count.ahead_count > 0)  &&
                        (chunk_location_count.behind_count > 0)) {
                    // movement in more than one direction
                    continue;
                }

                if ((chunk_location_count.ahead_count == 0)  &&
                        (chunk_location_count.behind_count == 0)) {
                    // possible movement, but the chunks don't touch
                    continue;
                }
#endif

                // check that the matching faces have the same dimensions
                for (int i=0; i<base_chunk->details->ndims; i++) {
                    if ((i != dim_with_movement) &&
                        (base_chunk->details->count[i] != candidate_chunk->details->count[i])) {
                        failed=TRUE;
                        break;
                    }
                }
                if (failed) continue;

                /*
                 * Do NOT uncomment these print_chunk() lines in production code.
                 * They are *very* slow even if the debug level is set low and
                 * nothing is being logged.
                 */
//                netcdf_debug_level=LOG_ALL;
//                if (DEBUG > 3) printf("*****************\n");
//                if (DEBUG > 3) printf("base chunk\n");
//                if (DEBUG > 3) printf("*****************\n");
//                if (base_chunk != NULL)      print_chunk(base_chunk);
//                if (DEBUG > 3) printf("*****************\n");
//                if (DEBUG > 3) printf("candidate chunk\n");
//                if (DEBUG > 3) printf("*****************\n");
//                if (candidate_chunk != NULL) print_chunk(candidate_chunk);
//                netcdf_debug_level=old;

                if ((chunk_location_count.ahead_count == 1)  &&
                        (chunk_location_count.behind_count == 0) &&
                        (chunk_location_count.same_count == agg_dims-1)) {
                    // aggregation is base + candidate
                    new_chunk = aggregate_chunks(base_chunk, candidate_chunk, dim_with_movement);
                } else if ((chunk_location_count.ahead_count == 0)  &&
                        (chunk_location_count.behind_count == 1) &&
                        (chunk_location_count.same_count == agg_dims-1)) {
                    // aggregation is candidate + base
                    new_chunk = aggregate_chunks(candidate_chunk, base_chunk, dim_with_movement);
                } else {
                    // chunks aren't aligned
                    //if (DEBUG > 3) printf("**********\nchunks are not aligned\n**********\n");
                    continue;
                }

                assert(new_chunk != NULL);

                /*
                 * Do NOT uncomment these print_chunk() lines in production code.
                 * They are *very* slow even if the debug level is set low and
                 * nothing is being logged.
                 */
//                netcdf_debug_level=LOG_ALL;
//                if (DEBUG > 3) printf("*****************\n");
//                if (DEBUG > 3) printf("new chunk\n");
//                if (DEBUG > 3) printf("*****************\n");
//                if (new_chunk != NULL)       print_chunk(new_chunk);
//                netcdf_debug_level=old;

                var_details->chunks->remove(base_chunk);
                var_details->chunks->remove(candidate_chunk);
                delete base_chunk;
                delete candidate_chunk;

                agg_chunks.push_back(new_chunk);

                aggregation_success = TRUE;
                success_this_pass = TRUE;

                break;
            }
            if (success_this_pass == TRUE) break;
        }
        chunks_iterator_t agg_iter = agg_chunks.begin();
        for (;agg_iter != agg_chunks.end();agg_iter++) {
            var_details->chunks->push_back(*agg_iter);
        }
        agg_chunks.clear();

        delete[] offset_diff;
    }

    if (DEBUG > 4) {
        printf("*****************\n");
        printf("end aggregation (begin list)\n");
        printf("*****************\n");
        int chunk_count;
        aggregation_chunk_details_t **chunks = get_chunks(fd, var_name, &chunk_count);
        for (int i=0;i<chunk_count;i++) {
            print_chunk(chunks[i]);
        }
        free(chunks);
        printf("*****************\n");
        printf("end aggregation (end list)\n");
        printf("*****************\n");
    }

//    netcdf_debug_level=LOG_ALL;
    chunks_iterator_t dst_iter=var_details->chunks->begin();
    for(;dst_iter != var_details->chunks->end();dst_iter++) {
        chunk_details_iterator_t component_iter=(*dst_iter)->component_chunks.begin();
        if (((*dst_iter)->details->buf == NULL) && ((*dst_iter)->details->len > 0)) {
            (*dst_iter)->details->buf = (char *)malloc((*dst_iter)->details->len);
//            if (DEBUG > 3) printf("allocated dst_iter->details->buf(%p), len(%ld)\n",
//                    (*dst_iter)->details->buf,
//                    (*dst_iter)->details->len);
        } else {
//            if (DEBUG > 3) printf("did not allocate dst_iter->details->buf(%p)\n", (*dst_iter)->details->buf);
        }
        for(;component_iter != (*dst_iter)->component_chunks.end();component_iter++) {
            //if (DEBUG > 3) printf("copying component\n");
            copy_chunk(*component_iter, (*dst_iter)->details);
            //if (DEBUG > 3) printf("destroying component\n");
            destroy_chunk(*component_iter);
        }
        (*dst_iter)->component_chunks.clear();
    }
//    netcdf_debug_level=old;

//    netcdf_debug_level=LOG_ALL;
    //if (DEBUG > 3) printf("*****************\n");
    //if (DEBUG > 3) printf("chunks after aggregation\n");
    //if (DEBUG > 3) printf("*****************\n");
//    base_iter = var_details->chunks->begin();
//    for (;base_iter != var_details->chunks->end(); ++base_iter) {
//        base_chunk = *base_iter;
//        if (base_chunk != NULL)
//            print_chunk(base_chunk);
//    }
//    netcdf_debug_level=old;

    return(aggregation_success);
}

/*
 * Aggregate all variables in the file.
 *
 */
int try_aggregation(const int fd)
{
    int aggregation_success=FALSE;

    file_details_t    *file_details=NULL;
    var_map_iterator_t var_iter;
    per_var_details_t *var_details=NULL;

    if (DEBUG > 3) printf("entered try_aggregation - fd(%d)\n", fd);

    file_details = open_file_map[fd];
    if (file_details == NULL) {
        return(aggregation_success);
    }
    var_iter = file_details->vars.begin();
    for (; var_iter != file_details->vars.end(); var_iter++) {
        var_details = var_iter->second;
        if (var_details == NULL) {
//            if (DEBUG > 3) printf("var_details==NULL.  continuing\n");
            continue;
        } else {
//            if (DEBUG > 3) printf("aggregating var_name(%s)\n", var_details->var_name);
            while(try_aggregation(fd, var_details->var_name) == TRUE);
        }
    }

    aggregation_success = TRUE;

    return(aggregation_success);
}

int aggregate_data_ready_to_write(const int fd, const char *var_name)
{
//    file_details_t *details = open_file_map[fd];
//    int chunks_needed=0;

//    if (details->num_participants > 0) {
//        chunks_needed = details->num_participants;
//    } else {
//        chunks_needed = details->participants->size();
//    }
//
//    if (details->vars[var_name]->chunks_received == chunks_needed) {
//        return TRUE;
//    }

    return FALSE;
}

int cache_data_ready_to_write(const int fd, const char *var_name)
{
//    file_details_t *details = open_file_map[fd];
//    int chunks_needed=0;

//    if (details->num_participants > 0) {
//        chunks_needed = details->num_participants;
//    } else {
//        chunks_needed = details->participants->size();
//    }
//
//    if (details->vars[var_name]->chunks_received == chunks_needed) {
//        return TRUE;
//    }

    return FALSE;
}

aggregation_chunk_details_t **get_chunks(const int fd, const char *var_name, int *chunk_count)
{
    file_details_t *details=NULL;
    per_var_details_t *var_details=NULL;
    aggregation_chunk_details_t **chunks=NULL;
    chunks_iterator_t iter;

    if (DEBUG > 3) printf("entered get_chunks - fd(%d) var_name(%s)\n", fd, var_name);

    *chunk_count=0;

    details = open_file_map[fd];
    if (details == NULL) {
        return(NULL);
    }
    var_details = details->vars[var_name];
    if (var_details == NULL) {
        return(NULL);
    }

    *chunk_count = details->vars[var_name]->chunks->size();

    if (DEBUG > 3) printf("found %d chunks to return\n", *chunk_count);

    if (*chunk_count == 0) {
        return(NULL);
    }
    chunks = (aggregation_chunk_details_t **)malloc(*chunk_count*sizeof(aggregation_chunk_details_t *));

    var_details->chunks->sort(compare_chunks_for_caching);

    iter = var_details->chunks->begin();
    for (int i=0;iter != var_details->chunks->end(); ++iter,i++) {
        chunks[i] = (*iter)->details;
//        print_chunk(chunks[i]);
    }

    //if (DEBUG > 3) printf("finished\n");

    return(chunks);
}

aggregation_chunk_details_t **get_chunks(const int fd, int *chunk_count)
{
    file_details_t  *details=NULL;
    var_map_iterator_t var_iter;
    per_var_details_t *var_details=NULL;
    aggregation_chunk_details_t **chunks=NULL;
    chunks_iterator_t chunks_iter;

    if (DEBUG > 3) printf("entered get_chunks - fd(%d)\n", fd);

    *chunk_count=0;

    details = open_file_map[fd];
    if (details == NULL) {
        return(NULL);
    }
    var_iter = details->vars.begin();
    for (; var_iter != details->vars.end(); ++var_iter) {
        var_details = var_iter->second;
        *chunk_count += var_details->chunks->size();
    }

    if (DEBUG > 3) printf("found %d chunks to return\n", *chunk_count);

    if (*chunk_count == 0) {
        return(NULL);
    }
    chunks = (aggregation_chunk_details_t **)malloc(*chunk_count*sizeof(aggregation_chunk_details_t *));

    int i=0;
    var_iter = details->vars.begin();
    for (; var_iter != details->vars.end(); var_iter++) {
        var_details = var_iter->second;
        var_details->chunks->sort(compare_chunks_for_caching);
        chunks_iter = var_details->chunks->begin();
        for (;chunks_iter != var_details->chunks->end(); ++chunks_iter,i++) {
            chunks[i] = (*chunks_iter)->details;
//            print_chunk(chunks[i]);
        }
    }

    //if (DEBUG > 3) printf("finished\n");

    return(chunks);
}
