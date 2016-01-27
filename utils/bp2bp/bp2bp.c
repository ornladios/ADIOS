/*
 * in the COPYING file in the top level directory of this source distribution.
 *
 * Copyright (c) 2008 - 2009.  UT-BATTELLE, LLC. All rights reserved.
 */

/* ADIOS bp2bp utility
 *  read all variables and attributes from
 *    all groups in a BP file and output this to a adios file
 *
 * This is a sequential program.
 */


/* Now we have to get the last plane section working, to divide up any way we want.
   Then we can divide up the pieces with not just the last dimension, but also the
   dimension before that.... I.e. we can can make larger chunks....
   The idea is that we want as few chunks as possible, given the total number of procs
   which will read in the the data
 */

#ifndef _GNU_SOURCE
#   define _GNU_SOURCE
#endif


#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <unistd.h>
#include <inttypes.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <getopt.h>
#include <errno.h>
#include <limits.h>   // LONG_MAX36

#include <math.h>     // NAN
#include <libgen.h>   // basename
#include <regex.h>    // regular expression matching
#include <fnmatch.h>  // shell pattern matching

#include "mpi.h"
#include "adios_read.h"
#include "adios_types.h"
#include "adios.h"


#ifdef DMALLOC
#include "dmalloc.h"
#endif

MPI_Comm   comm = MPI_COMM_WORLD;
int        gflag = 10;


#define MAX_DIMS 100
#define DEBUG 0
#define PERFORMANCE_CHECK 0
#define TIMING 100

void checkOverflow(int loc, ADIOS_VARINFO* v, uint64_t* s, uint64_t* c);
int getTypeInfo( enum ADIOS_DATATYPES adiosvartype, int* elemsize);
void rS(ADIOS_VARINFO* v, uint64_t* s, uint64_t by, int rank);
void calcC(uint64_t chunk_size, ADIOS_VARINFO* v, uint64_t* c);
uint64_t calcChunkSize(uint64_t total_size, uint64_t mne, int np);
uint64_t checkBound(ADIOS_VARINFO* v, uint64_t* s, uint64_t* c, uint64_t* uc, uint64_t chunk_size);
void arrCopy(uint64_t* from, uint64_t* to);
void getbasename (char *path, char **dirname, char **basename);
int print_data(void *data, int item, enum ADIOS_DATATYPES adiosvartype);


int main (int argc, char ** argv) {
    //For varriable definitions:
    //gbounds = global bounds string, lbounds = local bounds string, offs = offset string, tstring = temp string to hold temperary stuff
    char       gbounds[1007], lbounds[1007], offs[1007],tstring[100];
    //size = number of cores,  gidx = adios group index
    int        rank, size, gidx, i, j, k, ii;
    //data = pointer to read-in data
    void       * data = NULL;
    uint64_t   s[] = {0,0,0,0,0,0,0,0,0,0};  //starting offset
    uint64_t   c[] = {1,1,1,1,1,1,1,1,1,1};  //chunk block array
    uint64_t   bytes_read = 0;
    int        element_size;
    int64_t    new_adios_group, m_adios_file;
    uint64_t   var_size;  //portion_bound,
    uint64_t   adios_groupsize, adios_totalsize;
    int        read_buffer;        //possible maximum size you the user would like for each chunk in MB
    int           write_buffer = 1536;  //actual buffer size you use in MB
    int        itime;
    int        WRITEME=1;
    uint64_t   chunk_size;   //chunk size in # of elements
    char      *var_path, *var_name; // full path cut into dir path and name
    MPI_Init(&argc,&argv);
    MPI_Comm_rank(comm,&rank);
    MPI_Comm_size(comm,&size);

    // timing numbers
    // we will time:
    // 0: adios_open, adios_group_size
    // 1: the total time to read in the data
    // 2: times around each write (will only work if we do NOT buffer....
    // 3: the time in the close
    // 4: fopen, fclose
    // 5: total time
    // timers: the total I/O time
    int        timers = 6;
    double     start_time[timers], end_time[timers], total_time[timers];

    if (TIMING==100) {
        for (itime=0;itime<timers;itime++) {
            start_time[itime] = 0;
            end_time[itime] = 0;
            total_time[itime]=0;
        }
        //MPI_Barrier(MPI_COMM_WORLD);
        start_time[5] = MPI_Wtime();
    }

    if(rank==0)
        printf("converting...\n");

    if (argc < 5) {
        if (rank==0) printf("Usage: %s <BP-file> <ADIOS-file> read_buffer(MB) write_buffer(MB) METHOD (LUSTRE_strip_count) (LUSTRE_strip_size) (LUSTRE_block_size)\n", argv[0]);
        return 1;
    }



    if(TIMING==100)
        start_time[4] = MPI_Wtime();
    ADIOS_FILE * f = adios_fopen (argv[1], MPI_COMM_SELF);
    if(TIMING==100){
        end_time[4] = MPI_Wtime();
        total_time[4] = end_time[4]-start_time[4];
    }
    adios_init_noxml(comm); // no xml will be used to write the new adios file
    read_buffer = atoi(argv[3]);
    write_buffer = atoi(argv[4]);
    adios_allocate_buffer (ADIOS_BUFFER_ALLOC_NOW, write_buffer); // allocate MB buffer



    if (f == NULL) {
        printf("rank=%d, file cant be opened\n", rank);
        if (DEBUG) printf ("%s\n", adios_errmsg());
        return -1;
    }


    for (gidx = 0; gidx < f->groups_count; gidx++) {    //group part
        adios_groupsize = 0;
        ADIOS_GROUP * g = adios_gopen (f, f->group_namelist[gidx]);


        if (g == NULL) {
            if (DEBUG) printf ("%s\n", adios_errmsg());
            printf("rank %d: group cannot be opened.\n", rank);
            return -1;
        }
        /* First create all of the groups */
        // now I need to create this group in the file that will be written

        adios_declare_group(&new_adios_group,f->group_namelist[gidx],"",adios_flag_yes);


        if(strcmp(argv[5],"MPI_LUSTRE")!=0)   //see whether or not the user uses MPI_LUSTRE method
            adios_select_method (new_adios_group, argv[5], "", "");  //non-MPI_LUSTRE methods... like MPI, POSIX....
        else{
            char lustre_pars[1000];
            strcpy(lustre_pars, "");
            strcat(lustre_pars, "stripe_count=");
            sprintf(tstring, "%d", atoi(argv[6]));
            strcat(lustre_pars, tstring);
            strcat(lustre_pars, ",stripe_size=");
            sprintf(tstring, "%d", atoi(argv[7]));
            strcat(lustre_pars, tstring);
            strcat(lustre_pars, ",block_size=");
            sprintf(tstring, "%d", atoi(argv[8]));
            strcat(lustre_pars, tstring);

            if(rank==0)
                printf("lustre_pars=%s\n", lustre_pars);

            adios_select_method (new_adios_group, argv[5], lustre_pars, "");  //Use MPI Lustre method

        }



        // variable definition part
        for (i = 0; i < g->vars_count; i++) {
            ADIOS_VARINFO * v = adios_inq_var_byid (g, i);
            getbasename (g->var_namelist[i], &var_path, &var_name);

            if (v->ndim == 0) 
            {   
                // scalars: every process does them the same.
                adios_define_var(new_adios_group,var_name,var_path,v->type,0,0,0);
                getTypeInfo( v->type, &element_size);    //element_size is size per element based on its type
                if (v->type == adios_string) {  //special case when the scalar is string.
                    adios_groupsize += strlen(v->value);
                } else {
                    adios_groupsize += element_size;
                }
            } 
            else 
            { 
                // vector variables
                getTypeInfo( v->type, &element_size);
                var_size=1;
                for (ii=0;ii<v->ndim;ii++) {
                    var_size*=v->dims[ii];
                }
                uint64_t total_size = var_size;  //total_size tells you the total number of elements in the current vector variable
                var_size*=element_size; //var_size tells you the size of the current vector variable in bytess

                //re-initialize the s and c variables
                for(j=0; j<v->ndim; j++){
                    s[j] = 0;
                    c[j] = 1;
                }

                //find the approximate chunk_size you would like to use.
                chunk_size = calcChunkSize(total_size, read_buffer*1024*1024/element_size, size);

                //set the chunk block array with the total size as close to chunk_size as possible
                calcC(chunk_size, v, c);
                strcpy(lbounds,"");
                for(j=0; j<v->ndim; j++){
                    sprintf(tstring, "%" PRId64 ",", c[j]);
                    strcat(lbounds, tstring);
                }
                printf("rank=%d, name=%s, chunk_size1=%" PRId64 " c[]=%s\n",rank,g->var_namelist[i],chunk_size,lbounds);


                chunk_size = 1;
                for(ii=0; ii<v->ndim; ii++)            //reset chunk_size based on the created c. Now the chunk_size is exact.
                    chunk_size *= c[ii];

                //current step points to where the process is in processing the vector. First sets with respect to rank.
                uint64_t current_step = rank*chunk_size;

                //First advance the starting point s by current_step. Of course, you don't do it if the current_step exceeds total_size.
                if(current_step<total_size)
                    rS(v, s, current_step, rank);

                uint64_t elements_defined = 0;  //First, the number of elements you have defined is 0.

                //You (the process) process your part of the vector when your current_step is smaller than the total_size
                while(current_step < total_size)
                {
                    //ts, temporary s, is introduced for the sake of the inner do while loop below. Copy s to ts.
                    uint64_t ts[] = {0,0,0,0,0,0,0,0,0,0};
                    arrCopy(s, ts);

                    //for every outer while iteration, you always have the whole chunk_size remained to process.
                    uint64_t remain_chunk = chunk_size;
                    if(current_step+chunk_size>total_size) //except when you are nearing the end of the vector....
                        remain_chunk = total_size-current_step;

                    //tc, temporary c, is introduced for the sake of the inner do while loop below. Copy s to tc.
                    uint64_t tc[] = {1,1,1,1,1,1,1,1,1,1};
                    arrCopy(c, tc);

                    do{
                        //how much of the remain chunk you wanna process? initially you think you can do all of it....
                        uint64_t used_chunk = remain_chunk;

                        //you feel like you should process the vector with tc block size, but given ts, you might go over bound.
                        uint64_t uc[] = {1,1,1,1,1,1,1,1,1,1};
                        //so you verify it by setting a new legit chunck block uc, and getting a new remain_chunk.
                        remain_chunk = checkBound(v, ts, tc, uc, remain_chunk);

                        //you check whether or not ts+uc goes over the bound. This is just checking to make sure there's no error.
                        //Thereotically, there should be no problem at all.
                        checkOverflow(0, v, ts, uc);


                        //the below code fragment simply calculates gbounds, and sets place holders for lbounds and offs.
                        strcpy(gbounds,"");
                        strcpy(lbounds,"");
                        strcpy(offs,"");

                        for(j=0; j<v->ndim-1; j++){
                            sprintf(tstring, "%d,", (int)v->dims[j]);
                            strcat(gbounds, tstring);
                            //sprintf(tstring, "ldim%d_%s,", j, var_name);
                            sprintf(tstring, "ldim%d,", j);
                            strcat(lbounds, tstring);
                            //sprintf(tstring, "offs%d_%s,", j, var_name);
                            sprintf(tstring, "offs%d,", j);
                            strcat(offs, tstring);
                        }

                        sprintf(tstring, "%d", (int)v->dims[v->ndim-1]);
                        strcat(gbounds, tstring);
                        //sprintf(tstring, "ldim%d_%s", v->ndim-1, var_name);
                        sprintf(tstring, "ldim%d", v->ndim-1);
                        strcat(lbounds, tstring);
                        //sprintf(tstring, "offs%d_%s", v->ndim-1, var_name);
                        sprintf(tstring, "offs%d", v->ndim-1);
                        strcat(offs, tstring);

                        //sprintf(tstring, "%d", v->ndim);
                        for(j=0; j<v->ndim; j++){
                            //sprintf(tstring, "ldim%d_%s", j, var_name);
                            sprintf(tstring, "ldim%d", j);
                            adios_define_var(new_adios_group, tstring, "bp2bp", adios_unsigned_long, 0, 0, 0);
                            //sprintf(tstring, "offs%d_%s", j, var_name);
                            sprintf(tstring, "offs%d", j);
                            adios_define_var(new_adios_group, tstring, "bp2bp", adios_unsigned_long, 0, 0, 0);
                        }

                        adios_define_var(new_adios_group,var_name,var_path,v->type,lbounds,gbounds,offs);


                        if (DEBUG){
                            strcpy(lbounds,"");
                            strcpy(offs,"");
                            for(j=0; j<v->ndim; j++){
                                sprintf(tstring, "%" PRId64 ",", ts[j]);
                                strcat(offs, tstring);
                                sprintf(tstring, "%" PRId64 ",", uc[j]);
                                strcat(lbounds, tstring);
                            }

                            printf("rank=%d, name=%s, gbounds=%s: lbounds=%s: offs=%s \n",rank,g->var_namelist[i],gbounds, lbounds, offs);
                        }

                        used_chunk -= remain_chunk; //you get the actual used_chunk here.
                        elements_defined += used_chunk;
                        if(remain_chunk!=0){
                            rS(v, ts, used_chunk, rank);  //advance ts by used_chunk.
                            for(k=0; k<10; k++)
                                tc[k] = 1;
                            calcC(remain_chunk, v, tc);   //based on the remain_chunk, calculate new tc chunk block remained to process.
                        }

                        adios_groupsize+= used_chunk*element_size+2*v->ndim*8;

                    }while(remain_chunk!=0);

                    current_step += size*chunk_size;  //once a whole chunk_size is processed, advance the current_step in roll-robin manner.

                    if(current_step<total_size){   //advance s in the same way.
                        rS(v, s, size*chunk_size, rank);
                    }
                }

                //beside checkOverflow above, here you check whether or not the total number of elements processed across processes matches
                //the total number of elements in the original vector.
                if(DEBUG){
                    uint64_t* sb = (uint64_t *) malloc(sizeof(uint64_t));
                    uint64_t* rb = (uint64_t *) malloc(sizeof(uint64_t));
                    sb[0] = elements_defined;
                    MPI_Reduce(sb,rb,1,MPI_UNSIGNED_LONG_LONG,MPI_SUM,0, comm);

                    if(rank==0 && rb[0]!=total_size)
                        printf("some array define mismatch. please use debug mode\n");
                    free(sb); free(rb);
                }
            }
            free (var_name);
            free (var_path);
        } // finished declaring all of the variables


        // Now we can define the attributes....
        for (i = 0; i < g->attrs_count; i++) {
            enum ADIOS_DATATYPES atype;
            int  asize;
            void *adata;
            adios_get_attr_byid (g, i, &atype, &asize, &adata);
            // if (DEBUG) printf("attribute name=%s\n",g->attr_namelist[i]);
            adios_define_attribute(new_adios_group,g->attr_namelist[i],"",atype,adata,0);
        }



        /*------------------------------ NOW WE WRITE -------------------------------------------- */
        // Now we have everything declared... now we need to write them out!!!!!!
        if (WRITEME==1) {
            // open up the file for writing....
            if (DEBUG) printf("rank=%d, opening file = %s, with group %s, size=%" PRId64 "\n",rank,argv[2],f->group_namelist[gidx],adios_groupsize);

            if(TIMING==100)
                start_time[0] = MPI_Wtime();

            adios_open(&m_adios_file, f->group_namelist[gidx],argv[2],"w",comm);
            adios_group_size( m_adios_file, adios_groupsize, &adios_totalsize);

            //get both the total adios_totalsize and total adios_groupsize summed across processes.
            uint64_t* sb = (uint64_t *) malloc(sizeof(uint64_t));;
            uint64_t* rb = (uint64_t *) malloc(sizeof(uint64_t));
            sb[0] = adios_groupsize;
            MPI_Reduce(sb,rb,1,MPI_UNSIGNED_LONG_LONG,MPI_SUM,0, comm);

            uint64_t* sb2 = (uint64_t *) malloc(sizeof(uint64_t));;
            uint64_t* rb2 = (uint64_t *) malloc(sizeof(uint64_t));
            sb2[0] = adios_totalsize;
            MPI_Reduce(sb2,rb2,1,MPI_UNSIGNED_LONG_LONG,MPI_SUM,0, comm);
            if(rank==0){
                printf("total adios_totalsize = %" PRId64 "\n", *rb2);
                printf("total adios_groupsize = %" PRId64 "\n", *rb);
            }
            free(sb); free(rb); free(sb2); free(rb2);

            if (TIMING==100) {
                end_time[0] = MPI_Wtime();
                total_time[0]+=end_time[0] - start_time[0];    //variable definition time taken
            }

            // now we have to write out the variables.... since they are all declared now
            // This will be the place we actually write out the data!!!!!!!!
            for (i = 0; i < g->vars_count; i++) {
                ADIOS_VARINFO * v = adios_inq_var_byid (g, i);
                getbasename (g->var_namelist[i], &var_path, &var_name);
                if (v->ndim == 0) 
                {
                    if (DEBUG) {
                        printf ("ADIOS WRITE SCALAR: rank=%d, name=%s value=",
                                rank,g->var_namelist[i]);
                        print_data (v->value, 0, v->type);
                        printf ("\n");
                    }
                    if (TIMING==100) {
                        start_time[2] = MPI_Wtime();
                    }
                    adios_write(m_adios_file,g->var_namelist[i],v->value);
                    if (TIMING==100) {
                        end_time[2] = MPI_Wtime();
                        total_time[2]+=end_time[2] - start_time[2];     //IO write time...
                    }
                } 
                else 
                {
                    for(j=0; j<v->ndim; j++){
                        s[j] = 0;
                        c[j] = 1;
                    }
                    getTypeInfo( v->type, &element_size);

                    uint64_t total_size = 1;
                    for (ii=0;ii<v->ndim;ii++)
                        total_size*=v->dims[ii];

                    chunk_size = calcChunkSize(total_size, read_buffer*1024*1024/element_size, size);
                    calcC(chunk_size, v, c);
                    chunk_size = 1;
                    for(ii=0; ii<v->ndim; ii++)
                        chunk_size *= c[ii];


                    uint64_t current_step = rank*chunk_size;
                    if(current_step<total_size)
                        rS(v, s, current_step, rank);

                    uint64_t elements_written = 0;

                    while(current_step < total_size)
                    {
                        uint64_t ts[] = {0,0,0,0,0,0,0,0,0,0};
                        arrCopy(s, ts);
                        uint64_t remain_chunk = chunk_size;
                        if(current_step+chunk_size>total_size)
                            remain_chunk = total_size-current_step;
                        uint64_t tc[] = {1,1,1,1,1,1,1,1,1,1};
                        arrCopy(c, tc);

                        do{
                            uint64_t uc[] = {1,1,1,1,1,1,1,1,1,1};
                            uint64_t used_chunk = remain_chunk;
                            remain_chunk = checkBound(v, ts, tc, uc, remain_chunk);

                            checkOverflow(1, v, ts, uc);

                            used_chunk -= remain_chunk;
                            elements_written += used_chunk;

                            //allocated space for data read-in
                            data = (void *) malloc(used_chunk*element_size);

                            if (TIMING==100) {
                                start_time[1] = MPI_Wtime();
                            }
                            if(PERFORMANCE_CHECK) printf("rank=%d, read start\n",rank);
                            bytes_read = adios_read_var_byid(g,v->varid,ts,uc,data);
                            if(PERFORMANCE_CHECK) printf("rank=%d, read end\n",rank);
                            if (TIMING==100) {
                                end_time[1] = MPI_Wtime();
                                total_time[1]+=end_time[1] -start_time[1];      //IO read time
                            }

                            if (DEBUG)
                                printf ("ADIOS WRITE: rank=%d, name=%s datasize=%" PRId64 "\n",rank,g->var_namelist[i],bytes_read);


                            if (TIMING==100) {
                                start_time[2] = MPI_Wtime();
                            }
                            if (DEBUG){
                                printf("rank=%d, write ts=",rank);
                                int k;
                                for(k=0; k<v->ndim; k++)
                                    printf("%" PRId64 ",", ts[k]);
                                printf("  uc=");
                                for(k=0; k<v->ndim; k++)
                                    printf("%" PRId64 ",", uc[k]);
                                printf("\n");
                            }

                            //local bounds and offets placeholders are not written out with actual values.
                            if(PERFORMANCE_CHECK) printf("rank=%d, adios write start\n", rank);
                            for(k=0; k<v->ndim; k++){
                                //sprintf(tstring, "ldim%d_%s", k, var_name);
                                sprintf(tstring, "ldim%d", k);
                                if (DEBUG) {
                                    printf ("ADIOS WRITE DIMENSION: rank=%d, name=%s value=",
                                            rank,tstring);
                                    print_data (&uc[k], 0, adios_unsigned_long);
                                    printf ("\n");
                                }
                                adios_write(m_adios_file, tstring, &uc[k]);

                                //sprintf(tstring, "offs%d_%s", k, var_name);
                                sprintf(tstring, "offs%d", k);
                                if (DEBUG) {
                                    printf ("ADIOS WRITE OFFSET: rank=%d, name=%s value=",
                                            rank,tstring);
                                    print_data (&ts[k], 0, adios_unsigned_long);
                                    printf ("\n");
                                }
                                adios_write(m_adios_file, tstring, &ts[k]);
                            }
                            adios_write(m_adios_file,g->var_namelist[i],data);
                            if(PERFORMANCE_CHECK) printf("rank=%d, adios write end\n", rank);


                            if (TIMING==100) {
                                end_time[2] = MPI_Wtime();
                                total_time[2]+=end_time[2] - start_time[2];   //IO write time
                            }

                            free(data);


                            if(remain_chunk!=0){
                                rS(v, ts, used_chunk, rank);
                                for(k=0; k<10; k++)
                                    tc[k] = 1;
                                calcC(remain_chunk, v, tc);
                            }

                        }while(remain_chunk!=0);

                        current_step += size*chunk_size;

                        if(current_step<total_size)
                            rS(v, s, size*chunk_size,rank);
                    }

                    if(DEBUG){
                        uint64_t* sb = (uint64_t *) malloc(sizeof(uint64_t));;
                        uint64_t* rb = (uint64_t *) malloc(sizeof(uint64_t));
                        sb[0] = elements_written;
                        MPI_Reduce(sb,rb,1,MPI_UNSIGNED_LONG_LONG,MPI_SUM,0, comm);
                        if(rank==0 && rb[0]!=total_size)
                            printf("some array read mismatch. please use debug mode\n");
                        free(sb); free(rb);
                    }
                }
                free (var_name);
                free (var_path);
            }// end of the writing of the variable..
            if (TIMING==100) {
                start_time[3] = MPI_Wtime();
            }
            if(PERFORMANCE_CHECK) printf("rank=%d, adios_close start\n", rank);
            adios_close(m_adios_file);
            if(PERFORMANCE_CHECK) printf("rank=%d, adios_close end\n", rank);
            if (TIMING==100) {
                end_time[3] = MPI_Wtime();
                total_time[3]+=end_time[3] - start_time[3];
            }
            adios_gclose(g);
        } //end of WRITEME
    } // end of all of the groups

    if(rank==0)
        printf("conversion done!\n");

    if(TIMING==100)
        start_time[4] = MPI_Wtime();
    adios_fclose(f);
    if(TIMING==100){
        end_time[4] = MPI_Wtime();
        total_time[4] = total_time[4]+end_time[4]-start_time[4];
    }
    adios_finalize(rank);


    // now, we write out the timing data, for each category, we give max, min, avg, std, all in seconds, across all processes.
    if(TIMING==100){

        // 0: adios_open, adios_group_size
        // 1: the total time to read in the data
        // 2: times around each write (will only work if we do NOT buffer....
        // 3: the time in the close
        // 4: fopen, fclose
        // 5: total time
        end_time[5] = MPI_Wtime();
        total_time[5] = end_time[5] - start_time[5];

        double sb[7];
        sb[0] = total_time[1]; sb[1] = total_time[4];   //read_var, fopen+fclose
        sb[2] = sb[0]+sb[1];
        sb[3] = total_time[0]; sb[4] = total_time[2]+total_time[3]; //adios_open+adios_group_size, write+close
        sb[5] = sb[3]+sb[4];
        sb[6] = total_time[5]; //total

        double * rb = NULL;

        if(rank==0)
            rb = (double *)malloc(size*7*sizeof(double));
        //MPI_Barrier(comm);
        MPI_Gather(sb, 7, MPI_DOUBLE, rb, 7, MPI_DOUBLE, 0, comm);

        if(rank==0){
            double read_avg1 = 0;
            double read_avg2 = 0;
            double tread_avg = 0;
            double write_avg1 = 0;
            double write_avg2 = 0;
            double twrite_avg = 0;
            double total_avg = 0;
            for(j=0; j<size; j++){
                read_avg1 += rb[7*j];
                read_avg2 += rb[7*j+1];
                tread_avg += rb[7*j+2];
                write_avg1 += rb[7*j+3];
                write_avg2 += rb[7*j+4];
                twrite_avg += rb[7*j+5];
                total_avg += rb[7*j+6];
            }
            read_avg1 /= size;
            read_avg2 /= size;
            tread_avg /= size;
            write_avg1 /= size;
            write_avg2 /= size;
            twrite_avg /= size;
            total_avg /= size;

            double read1_max = rb[0];
            double read1_min = rb[0];
            double read1_std = rb[0]-read_avg1; read1_std *= read1_std;

            double read2_max = rb[1];
            double read2_min = rb[1];
            double read2_std = rb[1]-read_avg2; read2_std *= read2_std;

            double tread_max = rb[2];
            double tread_min = rb[2];
            double tread_std = rb[2]-tread_avg; tread_std *= tread_std;

            double write1_max = rb[3];
            double write1_min = rb[3];
            double write1_std = rb[3]-write_avg1; write1_std *= write1_std;

            double write2_max = rb[4];
            double write2_min = rb[4];
            double write2_std = rb[4]-write_avg2; write2_std *= write2_std;

            double twrite_max = rb[5];
            double twrite_min = rb[5];
            double twrite_std = rb[5]-twrite_avg; twrite_std *= twrite_std;

            double total_max = rb[6];
            double total_min = rb[6];
            double total_std = rb[6]-total_avg; total_std *= total_std;

            for(j=1; j<size; j++){
                if(rb[7*j]>read1_max)
                    read1_max = rb[7*j];
                else if(rb[7*j]<read1_min)
                    read1_min = rb[7*j];
                double std = rb[7*j]-read_avg1; std *= std;
                read1_std += std;

                if(rb[7*j+1]>read2_max)
                    read2_max = rb[7*j+1];
                else if(rb[7*j+1]<read2_min)
                    read2_min = rb[7*j+1];
                std = rb[7*j+1]-read_avg2; std *= std;
                read2_std += std;

                if(rb[7*j+2]>tread_max)
                    tread_max = rb[7*j+2];
                else if(rb[7*j+2]<tread_min)
                    tread_min = rb[7*j+2];
                std = rb[7*j+2]-tread_avg; std *= std;
                tread_std += std;

                if(rb[7*j+3]>write1_max)
                    write1_max = rb[7*j+3];
                else if(rb[7*j+3]<write1_min)
                    write1_min = rb[7*j+3];
                std = rb[7*j+3]-write_avg1; std *= std;
                write1_std += std;

                if(rb[7*j+4]>write2_max)
                    write2_max = rb[7*j+4];
                else if(rb[7*j+4]<write2_min)
                    write2_min = rb[7*j+4];
                std = rb[7*j+4]-write_avg2; std *= std;
                write2_std += std;

                if(rb[7*j+5]>twrite_max)
                    twrite_max = rb[7*j+5];
                else if(rb[7*j+5]<twrite_min)
                    twrite_min = rb[7*j+5];
                std = rb[7*j+5]-twrite_avg; std *= std;
                twrite_std += std;

                if(rb[7*j+6]>total_max)
                    total_max = rb[7*j+6];
                else if(rb[7*j+6]<total_min)
                    total_min = rb[7*j+6];
                std = rb[7*j+6]-total_avg; std *= std;
                total_std += std;
            }
            read1_std /= size;  read1_std = sqrt(read1_std);
            read2_std /= size;  read2_std = sqrt(read2_std);
            tread_std /= size;    tread_std = sqrt(tread_std);
            write1_std /= size; write1_std = sqrt(write1_std);
            write2_std /= size; write2_std = sqrt(write2_std);
            twrite_std /= size;    twrite_std = sqrt(twrite_std);
            total_std /= size; total_std = sqrt(total_std);

            printf("---type---                       max\tmin\tavg\tstd\n");
            printf("---read_var---                   %lf\t%lf\t%lf\t%lf\n", read1_max, read1_min, read_avg1, read1_std);
            printf("---fopen+fclose---               %lf\t%lf\t%lf\t%lf\n", read2_max, read2_min, read_avg2, read2_std);
            printf("---total_read---                 %lf\t%lf\t%lf\t%lf\n", tread_max, tread_min, tread_avg, tread_std);
            printf("---adios_open+adios_groupsize--- %lf\t%lf\t%lf\t%lf\n", write1_max, write1_min, write_avg1, write1_std);
            printf("---write+close---                %lf\t%lf\t%lf\t%lf\n", write2_max, write2_min, write_avg2, write2_std);
            printf("---total_write---                %lf\t%lf\t%lf\t%lf\n", twrite_max, twrite_min, twrite_avg, twrite_std);
            printf("---total---                      %lf\t%lf\t%lf\t%lf\n", total_max, total_min, total_avg, total_std);
            free(rb);

        }

    }

    //    if (TIMING==100 && rank==0) {
    //        printf("------------------------------------------------------------------\n");
    //        printf("Define variables     = %lf\n",total_time[0]);
    //        printf("Read   variables     = %lf\n",total_time[1]);
    //        printf("Write  variables     = %lf\n",total_time[2]);
    //        printf("Close File for write = %lf\n",total_time[3]);
    //        printf("Total write time     = %lf\n",total_time[2] + total_time[3]);
    //        for (itime=0;itime<timers-1;itime++)
    //            total_time[timers-1]+=total_time[itime];
    //        printf("Total I/O time       = %lf\n",total_time[timers-1]);
    //    }
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

