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


#define MAX_DIMS 100
#define DEBUG 1
#define verbose 1
#define TIMING 100

int getTypeInfo( enum ADIOS_DATATYPES adiosvartype, int* elemsize);



int main (int argc, char ** argv)  
{
    char       filename [256]; 
    char       gbounds[64], lbounds[64], offs[64],tstring[32], tstring2[32];
    int        rank, size, gidx, i, j, k,l, ii;
    enum       ADIOS_DATATYPES attr_type;
    void       * data = NULL;
    uint64_t   start[] = {0,0,0,0,0,0,0,0,0,0};
    uint64_t   s[] = {0,0,0,0,0,0,0,0,0,0};
    uint64_t   c[] = {0,0,0,0,0,0,0,0,0,0};
    uint64_t   bytes_read = 0;
    int64_t    dims [MAX_DIMS];
    int        st, ct,out_size;
    char       ** grp_name;
    int64_t    m_adios_group, m_adios_file;
    int64_t    msize,msize2, var_size;
    uint64_t   adios_groupsize, adios_totalsize;
    int        err;
    int64_t    readsize;
    int        *read_planes,*read_planes_end;
    int        buff_size=768;
    int        flag,itime;
    int        WRITEME=1, WRITEMELAST=0;
    MPI_Init(&argc,&argv);
    MPI_Comm_rank(comm,&rank);
    MPI_Comm_size(comm,&size);

    // timing numbers
    // we will time:
    // 0: defining variables, after reading in the header.
    // 1: the total time to read in the data
    // 2: times around each write (will only work if we do NOT buffer....
    // 3: the time in the close
    // 4: the time spent in malloc
    // timers: the total I/O time
    int        timers = 6;
    double     start_time[timers], end_time[timers], total_time[timers];

    if (argc < 5) {
      if (rank==0) printf("Usage: %s <BP-file> <ADIOS-file> buffer_size(MB) METHOD \n", argv[0]);
      return 1;
    }

    ADIOS_FILE * f = adios_fopen (argv[1], comm);
    adios_init_noxml(); // no xml will be used to write the new adios file
    adios_allocate_buffer (ADIOS_BUFFER_ALLOC_NOW, buff_size); // allocate MB buffer

    if (f == NULL) {
        if (DEBUG) printf ("%s\n", adios_errmsg());
	return -1;
    }
    buff_size = atoi(argv[3]);
    adios_groupsize = 0;
    /* For all groups */

    if (TIMING==100) {
      for (itime=0;itime<timers;itime++) {
	start_time[itime] = 0;
	end_time[itime] = 0;
	total_time[itime]=0;
      }
      MPI_Barrier(MPI_COMM_WORLD);
      start_time[0] = MPI_Wtime();
    }
    for (gidx = 0; gidx < f->groups_count; gidx++) { 
      if (DEBUG) printf("Group %s:\n", f->group_namelist[gidx]);
      ADIOS_GROUP * g = adios_gopen (f, f->group_namelist[gidx]);
      
      if (g == NULL) {
	if (DEBUG) printf ("%s\n", adios_errmsg());
	return -1;
      }
      /* First create all of the groups */
      // now I need to create this group in the file that will be written
      adios_declare_group(&m_adios_group,f->group_namelist[gidx],"",adios_flag_yes);
      adios_select_method (m_adios_group, argv[4], "", "");
      
      // for each variable, I need to know how  much to read in... I have a buffer, so I can
      // think that I can read in this much data..
      read_planes = (int *) malloc(sizeof(int)*g->vars_count);
      read_planes_end = (int *) malloc(sizeof(int)*g->vars_count);
      //for (i = 0; i < g->vars_count; i++) {  
      for (i = 13; i < 14; i++) { // SAK --- just for testing 1 variable output....
	ADIOS_VARINFO * v = adios_inq_var_byid (g, i);
	// now I can declare the variable...
	// if it's a scalar then we declare like:
	if (v->ndim == 0) {
	  adios_define_var(m_adios_group,g->var_namelist[i],"",v->type,0,0,0);
	  err = getTypeInfo( v->type, &out_size);
	  adios_groupsize+= out_size;
	} else {
	  // we will do some string maniupulation to set the global bounds, offsets, and local bounds... 
	  j = 0 ;
	  
	  // find out how many planes to read in... we can do this easily
	  err = getTypeInfo( v->type, &out_size); // this is the first multiplier....
	  var_size=1;
	  for (ii=1;ii<v->ndim;ii++) { //figure out the size of just 1 plane we will read... 
	    var_size*=v->dims[ii];
	  }
	  var_size*=out_size; // this is the size for 1 plane... now we can figure out
	  bytes_read = 0;
	  flag = 0;
	  ii=0;
	  while (bytes_read < 1024*1024*buff_size && ii<v->dims[0])  {
	    bytes_read+=var_size;
	    ii++;
	  }
	    
	  // try just having this as our parmeter
	  
	  read_planes[i] = (int) buff_size;

	  // for now, we need to make sure that the first dimension	    
	  msize =1;
	  for (ii =1;ii<v->ndim;ii++) {
	    msize = msize * v->dims[ii]; // we are calculating the size of the group to output
	  }
	  flag =0;
	  msize2 = 0;
	  for ( ii=rank*read_planes[i];ii<v->dims[0];ii+=read_planes[i]*size ) { //split the first dimension...
	    if (read_planes[i] + ii > v->dims[0]) {
	      read_planes_end[i] = v->dims[0] - ii;
	      if (DEBUG==1 && i==13) printf("rank=%d, read_planes_end = %d %d %d: %d\n",rank,v->dims[0],ii,read_planes[i],read_planes_end[i]);
	      flag = 1;
	    }
	    //  if (read_planes[i] + ii <= v->dims[0]) {
	    strcpy(gbounds,"");
	    strcpy(lbounds,"");
	    strcpy(offs,"");
	    // we divide the first dimension up.. not the last...
	    sprintf(tstring,"%d,",(int)v->dims[0]);
	    strcat(gbounds,tstring);
	    
	    if (flag==0) {
	      sprintf(tstring,"%d,",(int) read_planes[i]);
	      msize2 = msize2 + read_planes[i];
	    }
	    else if (flag==1) { 	      
	      sprintf(tstring,"%d,",(int) read_planes_end[i]);
	      msize2 = msize2 + read_planes_end[i];
	    }
	    strcat(lbounds,tstring);
	    
	    sprintf(tstring,"%d,",(int) ii);//the offset is the plane ii
	    strcat(offs,tstring);
	    
	    for (j=1;j<v->ndim-1;j++) {
	      sprintf(tstring,"%d,",(int)v->dims[j]);
	      strcat(gbounds,tstring);
	      
	      sprintf(tstring,"%d,",(int)v->dims[j]);
	      strcat(lbounds,tstring);
	      
	      sprintf(tstring,"%d,",(int) 0);
	      strcat(offs,tstring);
	    }
	    
	    sprintf(tstring,"%d\0",(int)v->dims[j]);
	    strcat(gbounds,tstring);
	    
	    sprintf(tstring,"%d\0",(int)v->dims[j]);
	    strcat(lbounds,tstring);
	    
	    sprintf(tstring,"%d\0",(int) 0);
	    strcat(offs,tstring);
	    
	    adios_define_var(m_adios_group,g->var_namelist[i],"",v->type,lbounds,gbounds,offs);
	    if (DEBUG&&i==13) printf("rank=%d, name=%s, gbounds=%s: lbounds=%s: offs=%s,flag=%d\n",rank,g->var_namelist[i],gbounds, lbounds, offs,flag);
	    
	    // we have a new variable for each sub block being written........
	}// end of looping for the splitting of the last dimension
	  
	  
	adios_groupsize+= out_size * msize * msize2;
	
	}
    } // finished declaring all of the variables
    // Now we can define the attributes....
    
    // if (DEBUG) printf("  Attributes=%d:\n", g->attrs_count);
    for (i = 0; i < g->attrs_count; i++) {
      enum ADIOS_DATATYPES atype;
      int  asize;
      void *adata;
      adios_get_attr_byid (g, i, &atype, &asize, &adata);
      // if (DEBUG) printf("attribute name=%s\n",g->attr_namelist[i]);
      adios_define_attribute(m_adios_group,g->attr_namelist[i],"",atype,adata,g->attr_namelist[i]);
    }
    if (TIMING==100) {
      MPI_Barrier(MPI_COMM_WORLD);
      
      end_time[0] = MPI_Wtime();
      total_time[0]+=end_time[0] - start_time[0];
    }
      /*------------------------------ NOW WE WRITE -------------------------------------------- */
      // Now we have everything declared... now we need to write them out!!!!!!
      
      
      if (WRITEME==1) {
	// open up the file for writing....
	if (DEBUG) printf("opening file = %s, with group %s, size=%lld\n",argv[2],f->group_namelist[gidx],adios_groupsize);
	adios_open(&m_adios_file, f->group_namelist[gidx],argv[2],"w",comm);
	adios_group_size( m_adios_file, adios_groupsize, &adios_totalsize);
	// now we have to write out the variables.... since they are all declared now
	// This will be the place we actually write out the data!!!!!!!!
        //for (i = 0; i < 14;i++) { 
	//for (i = 0; i < g->vars_count; i++) {
	for (i = 13; i < 14;i++) { 
	  ADIOS_VARINFO * v = adios_inq_var_byid (g, i);
	  // first we can write out the scalars very easily........
	  if (v->ndim == 0) {
	    
	    if (TIMING==100) {
	      MPI_Barrier(MPI_COMM_WORLD);
	      start_time[2] = MPI_Wtime();
	    }
	    adios_write(m_adios_file,g->var_namelist[i],v->value); //I think there will be a problem for strings.... Please check this out SAK
	    if (TIMING==100) {
	      MPI_Barrier(MPI_COMM_WORLD);
	      
	      end_time[2] = MPI_Wtime();
	      total_time[2]+=end_time[2] - start_time[2];
	    }
	  } else {
	    // now we will read in the variable, and then write it out....
	    // allocate the memory, read in the data, and then write it out....
	    for ( ii=rank*read_planes[i];ii<v->dims[0];ii+=read_planes[i]*size ) { //split the first dimension...
	      readsize = 1; 
	      for (j=1;j<(v->ndim);j++) {
		readsize = readsize * v->dims[j];
		s[j] = 0;
		c[j] = v->dims[j];
	      } // we are reading in with a 1 on the last dimension... so this will be the correct size...
	      readsize*= read_planes[i];
	      s[0] = ii; // this is the plane we are reading....
	      c[0] = read_planes[i]; // we are only reading read_plane planes....
	      flag = 0;
	      // Now we will care if we are on the last plane....
	      if (read_planes[i] + ii > v->dims[0]) {
		c[0] = read_planes_end[i];
	      }
	      
	      err = getTypeInfo( v->type, &out_size);
	      // in this version all chunks will be the same size, so only allocate 1 time
	      /*if (ii=rank*read_planes[i])*/ data = (void *) malloc(readsize*out_size);
	      if (TIMING==100) {
		MPI_Barrier(MPI_COMM_WORLD);
		start_time[1] = MPI_Wtime();
	      }
	      bytes_read = adios_read_var_byid(g,v->varid,s,c,data);
	      if (TIMING==100) {
		MPI_Barrier(MPI_COMM_WORLD);
		end_time[1] = MPI_Wtime();
		total_time[1]+=end_time[1] -start_time[1];
	      }
	      if (ii==0) 
printf("RW %f MB\n",(float) (bytes_read/1024./1024.));
	      // ok... now we write this out....
              if (DEBUG) printf ("ADIOS WRITE: rank=%d, name=%s datasize=%d\n",rank,g->var_namelist[i],bytes_read);
	      if (TIMING==100) {
		MPI_Barrier(MPI_COMM_WORLD);
		start_time[2] = MPI_Wtime();
	      }
	      if (DEBUG && i==13) printf("rank=%d, write s[0]=%d, c[0]=%d\n",rank,s[0],c[0]); 
	      adios_write(m_adios_file,g->var_namelist[i],data);
	      if (TIMING==100) {
		MPI_Barrier(MPI_COMM_WORLD);
		
		end_time[2] = MPI_Wtime();
		total_time[2]+=end_time[2] - start_time[2];
	      }
	      free(data);
	    }
	    // Now we are on the last plane
	    if (WRITEMELAST==1) {
	      readsize = 1;	      
	      for (j=1;j<(v->ndim);j++) {
		readsize = readsize * v->dims[j];
		s[j] = 0;
		c[j] = v->dims[j];
	      } // we are reading in with a 1 on the last dimension... so this will be the correct size...
	      readsize*= read_planes_end[i];
	      s[0] = ii; // this is the plane we are reading....
	      c[0] = read_planes_end[i]; // we are only reading read_plane planes....
	      err = getTypeInfo( v->type, &out_size);
	      data = (void *) malloc(readsize*out_size);
	      
	      if (TIMING==100) {
		MPI_Barrier(MPI_COMM_WORLD);
		start_time[1] = MPI_Wtime();
	      }
	      bytes_read = adios_read_var_byid(g,v->varid,s,c,data);
	      if (TIMING==100) {
		MPI_Barrier(MPI_COMM_WORLD);
		end_time[1] = MPI_Wtime();
		total_time[1]+=end_time[1] -start_time[1];
	      }
	      
	      printf("Last plane %d %d\n",s[0],c[0]);
	      // ok... now we write this out....
	      if (DEBUG) printf ("%s %d\n",g->var_namelist[i],bytes_read);


	      
	      if (TIMING==100) {
		MPI_Barrier(MPI_COMM_WORLD);
		start_time[2] = MPI_Wtime();
	      }
	      adios_write(m_adios_file,g->var_namelist[i],data);
	      if (TIMING==100) {
		MPI_Barrier(MPI_COMM_WORLD);

		end_time[2] = MPI_Wtime();
		total_time[2]+=end_time[2] - start_time[2];
	      }
	      free(data);
	    }
	  }
	}// end of the writing of the variable..
	if (TIMING==100) {
	  MPI_Barrier(MPI_COMM_WORLD);
	  start_time[3] = MPI_Wtime();
	}
	adios_close(m_adios_file);
	if (TIMING==100) {
	  MPI_Barrier(MPI_COMM_WORLD);	  
	  end_time[3] = MPI_Wtime();
	  total_time[3]+=end_time[3] - start_time[3];
	}
	adios_gclose(g);
      } //WRITEME
    } // end of all of the groups
    adios_fclose(f);
    adios_finalize(rank);

    if (TIMING==100 && rank==0) {
      printf("------------------------------------------------------------------\n");
      printf("Define variables     = %lf\n",total_time[0]);
      printf("Read   variables     = %lf\n",total_time[1]);
      printf("Write  variables     = %lf\n",total_time[2]);
      printf("Close File for write = %lf\n",total_time[3]);
      printf("Total write time     = %lf\n",total_time[2] + total_time[3]);
      for (itime=0;itime<timers-1;itime++) 
	total_time[timers-1]+=total_time[itime];
      printf("Total I/O time       = %lf\n",total_time[timers-1]);	     
    }
    MPI_Finalize();
    return(0);
}

int getTypeInfo( enum ADIOS_DATATYPES adiosvartype, int* elemsize)
{
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
