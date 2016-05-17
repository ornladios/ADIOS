/*
 * in the COPYING file in the top level directory of this source distribution.
 *
 * Copyright (c) 2008 - 2009.  UT-BATTELLE, LLC. All rights reserved.
 */

/* ADIOS list_methods utility
 *  list available 
 *    write methods
 *    read methods
 *    transform methods
 *
 * This is a sequential program but compiled with MPI to use libadios.a
 */

#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <inttypes.h>
#include <unistd.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <getopt.h>
#include <errno.h>

#include <math.h>     // NAN
#include <libgen.h>   // basename
#include <regex.h>    // regular expression matching
#include <fnmatch.h>  // shell pattern matching

#include "core/adios_internals.h" // write hooks and adios_transport_struct
#include "core/adios_read_hooks.h" // read hooks and adios_read_hooks_struct
#include "core/transforms/adios_transforms_hooks.h" 
#include "core/transforms/adios_transforms_hooks_read.h"
#include "core/transforms/adios_transforms_read.h"
#include "query/adios_query_hooks.h"


int print_data(void *data, int item, enum ADIOS_DATATYPES adiosvartype);


#ifdef WRITE
static struct adios_transport_struct * adios_transports = 0;
#endif
static struct adios_read_hooks_struct * adios_read_hooks = 0;
static struct adios_query_hooks_struct * adios_query_hooks = 0;

int main (int argc, char ** argv) {
    int  rank, size, i;
    MPI_Comm comm = MPI_COMM_WORLD;

    MPI_Init(&argc,&argv);
    MPI_Comm_rank(comm,&rank);
    MPI_Comm_size(comm,&size);


#ifdef WRITE
    adios_init_transports (&adios_transports);
#endif
    adios_read_hooks_init (&adios_read_hooks);
    adios_transform_read_init();
    adios_query_hooks_init(&adios_query_hooks);

    if(rank==0) {

#ifdef WRITE
        // print all write methods
        printf ("Available write methods (in XML <method> element or in adios_select_method()):\n");
        for (i = 0; i < ADIOS_METHOD_COUNT; i++) {    
            if (adios_transports[i].method_name) {
                printf("    \"%s\"\n", adios_transports[i].method_name);
            }
        }
#endif

        printf ("Available read methods (constants after #include \"adios_read.h\"):\n");
        for (i = 0; i < ADIOS_READ_METHOD_COUNT; i++) {    
            if (adios_read_hooks[i].method_name) {
                printf("    %s (=%d)\n", adios_read_hooks[i].method_name, i);
            }
        }

        printf ("Available data transformation methods (in XML transform tags in <var> elements):\n");
        for (i = (int)adios_transform_none; i < num_adios_transform_types; i++) {    
            if (adios_transform_is_implemented((enum ADIOS_TRANSFORM_TYPE)i)) {
            printf("    \"%s\"\t: %s\n", 
                    adios_transform_plugin_primary_xml_alias((enum ADIOS_TRANSFORM_TYPE)i),
                    adios_transform_plugin_desc((enum ADIOS_TRANSFORM_TYPE)i));
            }
        }

        printf ("Available query methods (in adios_query_set_method()):\n");
        for (i = 0; i < ADIOS_QUERY_METHOD_COUNT; i++) {
        	const enum ADIOS_QUERY_METHOD method = (enum ADIOS_QUERY_METHOD)i;
            if (adios_query_hooks[method].method_name) {
            	printf("    %s (=%d)\n", adios_query_hooks[method].method_name, i);
            }
        }
    }

    MPI_Barrier(comm);
    MPI_Finalize();
    return(0);
}

//check whether or not s+c goes over the global bound of v
//input:
//loc (tell you whether the overflow occurs, var definition or read/write)
//v (var info)
//s (offsets array)
//c (chunk block / local bounds array)

void checkOverflow(int loc, ADIOS_VARINFO* v, uint64_t* s, uint64_t* c) {
    int j;
    for(j=0; j<v->ndim; j++){
        if(s[j]+c[j]>v->dims[j]){
            if(loc==0)
                printf("in define: ");
            else //loc == 1
                printf("in read/write: ");
            printf("bound overflow happened. use debug mode\n");
        }
    }
}

//tell you what the size per element is based on the type
//input:
//adiosvartype (variable type structure)
//elemsize (pointer to the element size that you should set)
//output:
//tells you whether or not the adiosvartype is known.

int getTypeInfo( enum ADIOS_DATATYPES adiosvartype, int* elemsize){

    switch(adiosvartype) {
    case adios_unsigned_byte:
        *elemsize = 1;
        break;
    case adios_byte:
        *elemsize = 1;
        break;
    case adios_string:
        *elemsize = 1;
        break;

    case adios_unsigned_short:
        *elemsize = 2;
        break;
    case adios_short:
        *elemsize = 2;
        break;

    case adios_unsigned_integer:
        *elemsize = 4;
        break;
    case adios_integer:
        *elemsize = 4;
        break;

    case adios_unsigned_long:
        *elemsize = 8;
        break;
    case adios_long:
        *elemsize = 8;
        break;

    case adios_real:
        *elemsize = 4;
        break;

    case adios_double:
        *elemsize = 8;
        break;

    case adios_complex:
        *elemsize = 8;
        break;

    case adios_double_complex:
        *elemsize = 16;
        break;

    case adios_long_double: // do not know how to print
        //*elemsize = 16;
    default:
        return 1;
    }
    return 0;
}


//advance s by "by" number of elements.
//NOTE: you have to first make sure "by" doesn't go over "s" yourself. The function doesn't check this. If not, it could lead to error.
//input:
//v (variable info pointer)
//s (offset array pointer to start from)
//by (by how much elements do you want to advance?)
//rank (rank of your process)

void rS(ADIOS_VARINFO* v, uint64_t* s, uint64_t by, int rank){

    int q;
    uint64_t bulk = 1;

    for(q=1; q<v->ndim; q++)
        bulk *= v->dims[q];

    for(q=0; q<v->ndim; q++){

        //if(bulk == 0)
        //break;

        if(by == 0)
            break;

        uint64_t inc = by/bulk;

        if(inc >= 1){

            if(s[q]+inc<v->dims[q]){
                s[q] += inc;
            }else{
                //s[q-1]++;

                uint64_t r = 1;
                while(1){
                    if(s[q-r]+1 < v->dims[q-r]){
                        s[q-r]++;
                        break;
                    }else{
                        s[q-r] = 0;
                        r++;
                    }
                }

                uint64_t uinc = v->dims[q]-s[q];
                uint64_t new_inc = inc - uinc;
                s[q] = new_inc;
            }
            by -= inc*bulk;

        }


        if(q+1<v->ndim)
            bulk /= v->dims[q+1];
    }

}


//set chunk array block "c" based on the advised chunk_size
//NOTE: c has to be all 1's; otherwise, there would be error.
//input:
//chunk_size (advised chunk_size)
//v (variable info pointer)
//c (chunk array block you want to set)

void calcC(uint64_t chunk_size, ADIOS_VARINFO* v, uint64_t* c){
    int i;
    uint64_t tot = 1;
    uint64_t t;
    for(i=v->ndim-1; i>=0; i--){

        if(v->dims[i]*tot<=chunk_size){
            c[i] = v->dims[i];
            tot *= v->dims[i];
        }else{
            t = chunk_size/tot;
            c[i] = t;
            break;
        }
    }

}


//Calculate rough best estimate for chunk_size
//Input:
//total_size (total number of elements in the array)
//mne (maximum number of elements per for the chunk_size)
//np (number of cores you have)
//Output:
//chunk_size (in terms of the number of elements)

uint64_t calcChunkSize(uint64_t total_size, uint64_t mne, int np){

    uint64_t chunk_size = 0;
    if((total_size/np) <= mne)
        chunk_size = total_size/np;
    else
        chunk_size = mne;

    if(chunk_size<1)
        chunk_size = 1;

    return chunk_size;

}

//based on the the theoretical c you want to use, this function gives the real c, "uc" that you can use for the local bounds,
//without going out of the global bounds, starting from offset s
//NOTE: uc has to be all 1's
//Input:
//v (variable info pointer)
//s (offset pointer)
//c (chunk array block you want to use)
//uc (real chunk array calculated that you can use as local bounds)
//chunk_size (basically the product of dimensions of c)
//Output:
//remain_chunk (the number of elements still left to process after you process uc block of elements)

uint64_t checkBound(ADIOS_VARINFO* v, uint64_t* s, uint64_t* c, uint64_t* uc, uint64_t chunk_size){
    int i;
    uint64_t remain_chunk = chunk_size;
    int used_chunk = 1;
    for(i=v->ndim-1;i>=0;i--){
        if(s[i]+c[i]-1>=v->dims[i]){
            uc[i] = v->dims[i]-s[i];
            break;
        }
        uc[i] = c[i];
    }

    int j;

    for(j=0; j<v->ndim; j++)
        used_chunk *= uc[j];
    remain_chunk -= used_chunk;

    return remain_chunk;

}


//copy elements from "from" to "to"

void arrCopy(uint64_t* from, uint64_t* to){
    int i;
    for(i=0; i<10; i++)
        to[i] = from[i];
}


void getbasename (char *path, char **dirname, char **basename)
{
    char *work, *ptr;

    work = strdup(path);
    if (work[strlen(work)-1] == '/' && strlen(work)>1)
        work[strlen(work)-1] = '\0';

    ptr = rindex(work, '/');
    if (ptr && ptr != work) {
        // found a / and but not the beginning / of a full path
        ptr[0] = 0;
        *dirname = strdup(work);
        ptr[0] = '/';
        *basename = strdup(ptr+1);
    } else if (ptr == work) {
        // found / as the first character 
        *dirname = strdup("");
        *basename = strdup(ptr+1);
    } else {
        *dirname = NULL; //strdup(".");
        *basename = strdup(work);
    }
    free(work);
}


int print_data(void *data, int item, enum ADIOS_DATATYPES adiosvartype)
{
    if (data == NULL) {
        printf ( "null ");
        return 0;
    }
    // print next data item into vstr
    switch(adiosvartype) {
        case adios_unsigned_byte:
            printf ("%hhu", ((unsigned char *) data)[item]);
            break;
        case adios_byte:
            printf ("%hhd", ((signed char *) data)[item]);
            break;
        case adios_string:
            printf ("\"%s\"", ((char *) data)+item);
            break;
        case adios_string_array:
            // we expect one elemet of the array here
            printf("\"%s\"", *((char **)data+item));
            break;

        case adios_unsigned_short:
            printf ("%hu", ((unsigned short *) data)[item]);
            break;
        case adios_short:
            printf ("%hd", ((signed short *) data)[item]);
            break;

        case adios_unsigned_integer: 
            printf ("%u", ((unsigned int *) data)[item]);
            break;
        case adios_integer:    
            printf ("%d", ((signed int *) data)[item]);
            break;

        case adios_unsigned_long:
            printf ("%" PRIu64, ((uint64_t *) data)[item]);
            break;
        case adios_long:        
            printf ("%" PRId64, ((int64_t *) data)[item]);
            break;

        case adios_real:
            printf ("%g", ((float *) data)[item]);
            break;

        case adios_double:
            printf ("%g", ((double *) data)[item]);
            break;


        case adios_long_double:
            //printf ("%g ", ((double *) data)[item]);
            printf ("????????");
            break;


        case adios_complex:
            printf ("(%g,i%g)", ((float *) data)[2*item], ((float *) data)[2*item+1]);
            break;

        case adios_double_complex:
            printf ("(%g,i%g)", ((double *) data)[2*item], ((double *) data)[2*item+1]);
            break;

        case adios_unknown:
            break;
    } // end switch
    return 0;
}

