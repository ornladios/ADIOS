/* 
 * Staged write of ADIOS files using a staging method
 *
 * Copyright (c) 2008 - 2012.  UT-BATTELLE, LLC. All rights reserved.
 */


/* Staged write example code.
   Assumptions:
     - one output step fits into the memory of the staged writer.
       Actually, this means, even more memory is needed than the size of output.
       We need to read each variable while also buffering all of them for output.
     - output steps contain the same variable set (no changes in variables)
     - attributes are the same for all steps (will write only once here)
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <inttypes.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <errno.h>
#include <math.h>
#include "mpi.h"
#include "utils.h"
#include "decompose.h"
#include "adios.h"
#include "adios_read.h"
#include "adios_error.h"

#define MIN(a,b) ((a) < (b) ? a : b)

static enum ADIOS_READ_METHOD read_method = ADIOS_READ_METHOD_BP;
//static enum ADIOS_READ_METHOD read_method = ADIOS_READ_METHOD_DATASPACES;

double FUZZ_FACTOR = 0.0;


// Input arguments
char   infilename1[256];    // File1
char   infilename2[256];    // File2
char   fuzzfactor[256];     //fuzz factor
char   outfilename[256];   // File to write
char   methodname[16];     // ADIOS write method
char   methodparams[256];  // ADIOS write method


// Global variables
int         rank, numproc;
MPI_Comm    comm; 
ADIOS_FILE *f1, *f2;      // stream for reading
int64_t    fh;     // ADIOS output file handle
int64_t     gh;     // ADIOS group for output definitions
//uint64_t    write_total1, write_total2; // data size written by one processor
uint64_t    largest_block1; // the largest variable block one process reads
char     ** group_namelist1; // name of ADIOS group
char     ** group_namelist2;
char       *readbuf1, *readbuf2; // read buffer
double     fuzz_factor; //fuzz factor
int        verbose = 0;


int process_metadata();
int diff();
int compare_buffer(char * variable_name, void *data1, void *data2, int total, enum ADIOS_DATATYPES adiosvartype);
int compare_data(char * variable_name, void *data1, void *data2, int item, enum ADIOS_DATATYPES adiosvartype);

void printUsage(char *fname)
{
    print0("Usage: %s file1 file2 [-f fuzz_factor] [-v]\n"
           "    file1	Input file1 path\n"
           "    file2	Input file2 path\n"
           "    fuzz factor The difference cutoff for float/double\n",
           fname);
}


int processArgs(int argc, char ** argv)
{
    if (argc < 3) {
        printUsage(argv[0]);
        return 1;
    }
    strncpy(infilename1,     argv[1], sizeof(infilename1));
    strncpy(infilename2,    argv[2], sizeof(infilename2));
    //strncpy(fuzzfactor,    argv[3], sizeof(fuzzfactor));

    int option = 0;
    while ((option = getopt (argc, argv, "vf:")) != -1){
      switch (option){
        case 'v':
          verbose = 1;
          break;
        case 'f': 
          fuzz_factor = atof(optarg);
          break;
        case '?':
          return 1;
      }
    }

    return 0;
}


int main (int argc, char ** argv) 
{
    int         err;
    int         retval = 0;

    MPI_Init (&argc, &argv);
    comm = MPI_COMM_WORLD;
    MPI_Comm_rank (comm, &rank);
    MPI_Comm_size (comm, &numproc);

    if (processArgs(argc, argv)) {
        return 1;
    }
    
    print0("Input stream      = %s\n", infilename1);
    print0("Input stream      = %s\n", infilename2);
    print0("Fuzz factor       = %f\n", fuzz_factor);
    
    err = adios_read_init_method(read_method, comm, 
                                 "max_chunk_size=100; "
                                 "app_id =32767; \n"
                                 "verbose= 3;"
                                 "poll_interval  =  100;"
                                );

    if (!err) {
        print0 ("%s\n", adios_errmsg());
    }

    //adios_init_noxml(comm);

    print0 ("Waiting to open stream %s...\n", infilename1);
    f1 = adios_read_open_file(infilename1, read_method, comm);

    f2 = adios_read_open_file(infilename2, read_method, comm);

    if (adios_errno == err_file_not_found) 
    {
        print ("rank %d: File not found: %s\n", rank, adios_errmsg());
        retval = adios_errno;
    } 
    else if (f1 == NULL) {
        print ("rank %d: Error at opening file: %s\n", rank, adios_errmsg());
        retval = adios_errno;
    } 
    else 
    {

       print0 ("File info:\n");
       print0 ("  # of variables: %d	%d:\n", f1->nvars, f2->nvars);

       retval = process_metadata();
       if (retval) return retval;

       retval = diff();
       if (retval) return retval;


       adios_read_close (f1);
       adios_read_close (f2);
    }

    adios_read_finalize_method (read_method);
    MPI_Finalize ();
    return retval;
}


typedef struct {
    ADIOS_VARINFO * v;
    uint64_t        start[10];
    uint64_t        count[10];
    uint64_t        writesize; // size of subset this process writes, 0: do not write
    int		    cross_ref;
} VarInfo;

VarInfo * varinfo1, *varinfo2;
uint64_t sum_count1, sum_count2;

int process_metadata()
{
    int i, j;
    ADIOS_VARINFO *v1, *v2; // shortcut pointer


    varinfo1 = (VarInfo *) malloc (sizeof(VarInfo) * f1->nvars);
    if (!varinfo1) {
        print("ERROR: rank %d cannot allocate %zu bytes\n", rank, sizeof(VarInfo)*f1->nvars);
        return 1;
    }
    varinfo2 = (VarInfo *) malloc (sizeof(VarInfo) * f2->nvars);
    if (!varinfo2) {
        print("ERROR: rank %d cannot allocate %zu bytes\n", rank, sizeof(VarInfo)*f2->nvars);
        return 1;
    }

    //get varinfo of file1
    print0 ("--------------- metadata of %s ----------------------\n", infilename1); 
    for (i=0; i<f1->nvars; i++) 
    {
        //print0 ("%s Get info on variable %d: %s\n", infilename1, i, f1->var_namelist[i]); 
        varinfo1[i].v = adios_inq_var_byid (f1, i);
	      varinfo1[i].cross_ref= -1;
        v1 = varinfo1[i].v; // just a shortcut
        if (v1 == NULL) {
            print ("%s rank %d: ERROR: Variable %s inquiry failed: %s\n", 
                   infilename1, rank, f1->var_namelist[i], adios_errmsg());
            return 1;
        }

        // print variable type and dimensions
        print0("    %-9s  %s", adios_type_to_string(v1->type), f1->var_namelist[i]);
        if (v1->ndim > 0) {
            print0("[%" PRIu64, v1->dims[0]);
            for (j = 1; j < v1->ndim; j++)
                print0(", %" PRIu64, v1->dims[j]);
            print0("] :\n");
        } else {
            print0("\tscalar\n");
        }

     }

    print0 ("--------------- metadata of %s ----------------------\n", infilename2); 
    for (i=0; i<f2->nvars; i++) {
        //print0 ("%s Get info on variable %d: %s\n", infilename2, i, f2->var_namelist[i]); 
        varinfo2[i].v = adios_inq_var_byid (f2, i);
        varinfo2[i].cross_ref = -1;
        v2 = varinfo2[i].v; // just a shortcut
        if (v2 == NULL) {
            print ("%s rank %d: ERROR: Variable %s inquiry failed: %s\n", 
                   infilename2, rank, f2->var_namelist[i], adios_errmsg());
            return 1;
        }
        
        // print variable type and dimensions
        print0("    %-9s  %s", adios_type_to_string(v2->type), f2->var_namelist[i]);
        if (v2->ndim > 0) {
            print0("[%" PRIu64, v2->dims[0]);
            for (j = 1; j < v2->ndim; j++)
                print0(", %" PRIu64, v2->dims[j]);
            print0("] :\n");
        } else {
            print0("\tscalar\n");
        }
     }

    print0 ("--------------- end of metadata ----------------------\n"); 
     
    for(i=0; i<f1->nvars; i++){
	    v1 = varinfo1[i].v;
	    //for(j=i; j!=i-1; ((j++)%(f2->nvars))){
	    for(j=i; j!=i-1; j = (j % f2->nvars)+1){
	      v2 = varinfo2[j].v;
	      if(strcmp(f1->var_namelist[i], f2->var_namelist[j]) == 0 &&
	        strcmp(adios_type_to_string(v1->type), adios_type_to_string(v2->type))==0){
          if(v1->ndim == v2->ndim){
	          varinfo1[i].cross_ref = j;
	          varinfo2[j].cross_ref = i;
          }
	        break;
	      }
	    }
    }
     
   return 0;
}

int diff(){
  int i;
  ADIOS_VARINFO *v1, *v2; // shortcut pointer

  //process variables in one file but not in the other
  for(i=0; i<f1->nvars; i++){
    if(varinfo1[i].cross_ref == -1){
      print0("%s is in %s, not %s, or the datatype or dimension of the variable is different in two files\n", f1->var_namelist[i],infilename1, infilename2);
    }
  }

  for(i=0; i<f2->nvars; i++){
    if(varinfo2[i].cross_ref == -1){
      print0("%s is in %s, not %s, or the datatype or dimension of the variable is different in two files\n", f2->var_namelist[i], infilename2, infilename1);
    }
  }

  for (i=0; i<f1->nvars; i++){
	  v1 = varinfo1[i].v;
	  int cross_ref = varinfo1[i].cross_ref;
	  v2 = varinfo2[cross_ref].v;
	  
    if(varinfo1[i].cross_ref >=0 && v1->ndim == 0) //scalar 1 value already in mem of every processes
    {
	    if(!rank){
	     int ret = compare_data(f1->var_namelist[i], v1->value, v2->value, 0, v1->type);
       if(!ret) print0("%s is the same in %s and %s\n", f1->var_namelist[i], infilename1, infilename2);
	    }
      continue;
    }

    if(varinfo1[i].cross_ref >= 0 && v1->ndim >0){
      //decompose
      int *decomp_values = (int *)malloc(sizeof(int)*v1->ndim);
      int x;for(x=0; x<v1->ndim; x++) {decomp_values[x] = numproc;} 
      decompose (numproc, rank, v1->ndim, v1->dims, decomp_values,
                varinfo1[i].count, varinfo1[i].start, &sum_count1);
      decompose (numproc, rank, v2->ndim, v2->dims, decomp_values,
                varinfo2[cross_ref].count, varinfo2[cross_ref].start, &sum_count2);
      /*for(x=0; x<v1->ndim; x++) print("%d %d %d %d %d %d %d\n",rank, varinfo1[i].count[x], varinfo1[i].start[x], sum_count1,\
                                                                     varinfo2[cross_ref].start[x], varinfo2[cross_ref].count[x], sum_count2); */

      //read from files
      uint64_t size = adios_type_size (v1->type, v1->value);
      char *readbuf1 = malloc(sum_count1 * size);
      ADIOS_SELECTION *sel1 = adios_selection_boundingbox (v1->ndim,
                                                           varinfo1[i].start, 
                                                           varinfo1[i].count);
      adios_schedule_read_byid (f1, sel1, i, 0, 1, readbuf1);
      adios_perform_reads (f1, 1);   

      char *readbuf2 = malloc(sum_count1 * size);
      ADIOS_SELECTION *sel2 = adios_selection_boundingbox (v2->ndim,
                                                           varinfo2[cross_ref].start, 
                                                           varinfo2[cross_ref].count);
      adios_schedule_read_byid (f2, sel2, i, 0, 1, readbuf2);
      adios_perform_reads (f2, 1);   

      int sum_count = MIN(sum_count1, sum_count2);
      int ret = compare_data(f1->var_namelist[i], readbuf1, readbuf2, sum_count, v1->type);
      int allret;
      MPI_Reduce(&ret, &allret, 1, MPI_INT, MPI_SUM, 0, comm);
      if(allret > 0){
        if(!verbose){
          print0("%s has different values in %s and %s\n", f1->var_namelist[i], infilename1, infilename2);
        } else {
          print0("%s has %d different values in %s and %s\n", f1->var_namelist[i], allret, infilename1, infilename2);
        }
      } else {
        print0("%s is the same in %s and %s\n", f1->var_namelist[i], infilename1, infilename2);
      }

      free(decomp_values);
      free(readbuf1);
      free(readbuf2);
      adios_selection_delete(sel1);
      adios_selection_delete(sel2);
    }
  }//for
  return 0;
}

int compare_buffer(char * variable_name, void *data1, void *data2, int total, enum ADIOS_DATATYPES adiosvartype){
  int i;
  int total_diff = 0;
  for(i=0; i<total; i++){
    total_diff += compare_data(variable_name, data1, data2, i, adiosvartype);
  }
  return total_diff;
}

int compare_data(char * variable_name, void *data1, void *data2, int item, enum ADIOS_DATATYPES adiosvartype)
{
    int ret = 0;
    if (data1 == NULL || data2 == NULL) {
        print("null data");
        return ret;
    }


    // print next data item into vstr
    switch(adiosvartype) {
        case adios_unsigned_byte:
            if(((unsigned char *) data1)[item] !=  ((unsigned char *) data2)[item] )//not identical
            {
                print("%s : %hhu in %s | %hhu in %s\n", variable_name, ((unsigned char *) data1)[item], infilename1, ((unsigned char *) data2)[item], infilename2);
                ret++;
            }
            break;
        case adios_byte:
            if(((signed char *) data1)[item] != ((signed char *) data2)[item])//not identical
            {
                print("%s : %hhd in %s | %hhd in %s\n", variable_name, ((signed char *) data1)[item], infilename1, ((signed char *) data2)[item], infilename2);
                ret++;
            }
            break;
        case adios_string:
            if(strcmp(((char *) data1), ((char *) data2))!= 0 )//not identical
            {
                print("%s : %s in %s | %s in %s\n", variable_name, (char *) data1, infilename1, (char *) data2, infilename2);
                ret++;
            }
            break;
        case adios_unsigned_short:
            if(((unsigned short*) data1)[item] != ((unsigned short *) data2)[item])//not identical
            {
                print("%s : %hu in %s | %hu in %s\n", variable_name, ((unsigned short *) data1)[item], infilename1, ((unsigned short *) data2)[item], infilename2);
                ret++;
            }
            break;
        case adios_short:
            if(((signed short*) data1)[item] != ((signed short *) data2)[item])//not identical
            {
                print("%s : %hd in %s | %hd in %s\n", variable_name, ((signed short *) data1)[item], infilename1, ((signed short *) data2)[item], infilename2);
                ret++;
            }
            break;
        case adios_unsigned_integer:
            if(((unsigned int*) data1)[item] != ((unsigned int*) data2)[item])//not identical
            {
                print("%s : %u in %s | %u in %s\n", variable_name, ((unsigned int*) data1)[item], infilename1, ((unsigned int*) data2)[item], infilename2);
                ret++;
            }
            break;
        case adios_integer:
            if(((signed int*) data1)[item] != ((signed int*) data2)[item])//not identical
            {
                print("%s : %d in %s | %d in %s\n", variable_name, ((signed int*) data1)[item], infilename1, ((signed int*) data2)[item], infilename2);
                ret++;
            }
            break;
        case adios_unsigned_long:
            if(((unsigned long long*) data1)[item] != ((unsigned long long*) data2)[item])//not identical
            {
                print("%s : %llu in %s | %llu in %s\n", variable_name, ((unsigned long long*) data1)[item], infilename1, ((unsigned long long*) data2)[item], infilename2);
                ret++;
            }
            break;
        case adios_long:
            if(((unsigned long long*) data1)[item] != ((unsigned long long*) data2)[item])//not identical
            {
                print("%s : %lld in %s | %lld in %s\n", variable_name, ((signed long long*) data1)[item], infilename1, ((signed long long*) data2)[item], infilename2);
                ret++;
            }
            break;
        case adios_real:
            {
                float a, b;
                a = ((float *) data1)[item];
                b = ((float *) data2)[item];
                if(fabs(a-b)> fuzz_factor){
                    print("%s : %g in %s | %g in %s\n", variable_name, a, infilename1, b, infilename2);
                    ret++;
                }
                break;
            }
        case adios_double:
            {
                double aa, bb;
                aa = ((double *) data1)[item];
                bb = ((double *) data2)[item];
                if(fabs(aa-bb)> fuzz_factor){
                    print("%s : %g in %s | %g in %s\n", variable_name, aa, infilename1, bb, infilename2);
                    ret++;
                }
                break;
            }
        case adios_long_double:
            //fprintf(outf,(f ? format : "%g "), ((double *) data1)[item]);
            //                fprintf(outf,(f ? format : "????????"));
            break;
        case adios_complex:
            {
                float a11, a12, b11, b12;
                a11 = ((float *) data1)[2*item];
                a12 = ((float *) data1)[2*item+1];
                b11 = ((float *) data2)[2*item];
                b12 = ((float *) data2)[2*item+1];
                if(fabs(a11-b11)> fuzz_factor || fabs(a12-b12)>fuzz_factor){
                    print("%s : %g i%g in %s | %g i%g in %s\n", variable_name, a11, b11, infilename1, a12, b12, infilename2);
                    ret++;
                }

                break;
            }
        case adios_double_complex:
            {
                double a21, a22, b21, b22;
                a21 = ((float *) data1)[2*item];
                a22 = ((float *) data1)[2*item+1];
                b21 = ((float *) data2)[2*item];
                b22 = ((float *) data2)[2*item+1];
                if(fabs(a21-b21)> fuzz_factor || fabs(a22-b22)>fuzz_factor){
                    print("%s : %g i%g in %s | %g i%g in %s\n", variable_name, a21, b21, infilename1, a22, b22, infilename2);
                    ret++;
                }
                break;
            }
        default:
            break;
    } // end switch
    return ret;
}

