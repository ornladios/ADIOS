#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "public/adios.h"
#include "public/adios_read.h"
#include <iapi.h>

#include "../../src/query/fastbit_adios.h"

  char gVarNameFastbitIdxKey[10] = "key";
  char gVarNameFastbitIdxOffset[10] = "offsets";
  char gVarNameFastbitIdxBms[10] = "bms";

  char gGroupNameFastbitIdx[20] = "notNamed";

  char* gBinningOption = 0;

  int64_t       gAdios_group;
  int64_t       gAdios_write_file;

  time_t        lastMeasured;
  time_t        indexRefresh;
  time_t        fileStarted;

  long   _lastMeasuredMillis;
  long   _indexRefreshMillis;
  long   _fileStartedMillis;

  uint64_t sum_nb=0, sum_nk=0, sum_no=0;

//void printData(void* data, enum ADIOS_DATATYPES type, uint64_t size);

void defineFastbitVar(int nblocks, const char* name, int64_t* ids, int adiosType, uint64_t* localDim, const char* globalStr, uint64_t* offset)
{
  int i=0;
  for (i = 0; i < nblocks; i++) 
  {
    //int offset= i*5;
    char offsetStr[100] = "";
    if (offset != NULL) {
      sprintf(offsetStr, "%llu", offset[i]);
    } 

    char dimStr[100];
    sprintf(dimStr, "%llu", localDim[i]);
    ids[i] = adios_define_var (gAdios_group, name, "", adiosType, dimStr, globalStr, offsetStr);
    //adios_set_transform (ids[i], "identity");
  }
}


// k is numbered from 1 to sum_nblocks
void verifyData(ADIOS_FILE* f, ADIOS_VARINFO* v, int k, int timestep) 
{
  uint64_t blockBytes = adios_type_size (v->type, v->value);
  int j=0;

  if (v->ndim <= 0) {
    return;
  }

  //printf("verify block[%d]: ", k);
  
      for (j=0; j<v->ndim; j++) 
      {  
	  blockBytes *= v->blockinfo[k].count[j];
	  //printf("%llu:%llu ", v->blockinfo[k].start[j], v->blockinfo[k].count[j]);
      }

      void* data = NULL;
      data = malloc(blockBytes);  
      ADIOS_SELECTION* sel =  adios_selection_boundingbox (v->ndim, v->blockinfo[k].start, v->blockinfo[k].count);
      int err = adios_schedule_read_byid(f, sel, v->varid, timestep, 1, data);      
      if (!err) {	
	   err = adios_perform_reads(f, 1);
      }
      //fastbit_adios_util_printData(data, v->type, blockBytes/adios_type_size(v->type, v->value));
      adios_selection_delete(sel);
      free(data);	 
      data = NULL;
}



void assertErr(long int errCode, const char* exp, const char* refName) 
{
  if (errCode < 0) {
    printf("errCode=%ld %s %s\n", errCode, exp, refName);
    exit(EXIT_FAILURE);
  } 
}

void fastbitIndex(const char* datasetName, void* data, uint64_t blockSize, FastBitDataType ft, 
		  double**keys, uint64_t*nk, int64_t**offsets, uint64_t*no,		  
		  uint32_t**bms, uint64_t*nb)
{
  long int fastbitErr;

  fastbitErr = fastbit_iapi_register_array(datasetName, ft, data, blockSize);
  assertErr(fastbitErr, "failed to register array with", datasetName);
  
  fastbitErr = fastbit_iapi_build_index(datasetName, (const char*)gBinningOption);
  assertErr(fastbitErr, "failed to build idx with ", datasetName);
  
  
  fastbitErr = fastbit_iapi_deconstruct_index(datasetName, keys, nk, offsets, no, bms, nb);
  assertErr(fastbitErr, "failed with fastbit_iapi_deconstruct on ", datasetName);
  
  //printf("nk/no/nb %lld %lld %lld\n", *nk, *no, *nb);
  
  fastbit_iapi_free_all();

  /*free(offsets);

  free(keys);

  free(bms);
  */
}


void logTimeMillis(const char* notes) 
{
#ifndef SHOW_TIMESTUDY
    return;
#endif
  long ms = fastbit_adios_getCurrentTimeMillis();

  if (notes == NULL) {
     printf("\n");
     _indexRefreshMillis = ms;
  } else {
    long d = ms - _lastMeasuredMillis;
    printf("   ELAPSED millis: %ld \t%s\n", d, notes);
  }
  _lastMeasuredMillis = ms;
}


void sumLogTimeMillis(int stage) {
#ifndef SHOW_TIMESTUDY
    return;
#endif

  long ms = fastbit_adios_getCurrentTimeMillis();

  if (stage == -1) { // init
    _lastMeasuredMillis = ms;
    _indexRefreshMillis = ms;
    _fileStartedMillis = ms;
  } else if (stage == 0) { // whole program    
    printf("\n==> Total time spent to process this FILE: %ld millis.\n", ms - _fileStartedMillis);    
  } else if (stage == 1) { // block
    printf("==>  Total time spent to process this block: %ld millis. \n", ms - _indexRefreshMillis);
    _indexRefreshMillis = ms;
  } else if (stage == 2) { // variable
    //printf("==>  total time spent to process this variable: %ld millis. \n", ms - _indexRefreshMillis);
    //_indexRefreshMillis = ms;
  }
}


void logTime(const char* notes) {
  /*
  char buff[300];
  struct tm * timeinfo;
  time_t now = time(0);
  timeinfo = localtime (&now);
  strftime(buff, 30, "%b %d %H:%M:%S", timeinfo); 

  if (notes == NULL) {
    //printf("\n  ====== %s\n", buff); 
     printf("\n"); 
     indexRefresh = now;
  } else {
     double d = difftime(now, lastMeasured);
     printf("%s elapsed (sec) %4.0f %s\n", buff, d, notes); 
  }
  lastMeasured = now;
  */
}

void sumLogTime(int stage) {
  /*
  time_t now = time(0);
  struct tm* timeinfo = localtime(&now);
  if (stage == -1) { // init
    lastMeasured = now;
    indexRefresh = now;
    fileStarted = now;
  } else if (stage == 0) { // whole program    
    double d = difftime(now, fileStarted);
    printf("==> total time spent to process this FILE: %4.0f seconds.\n", d);    
  } else if (stage == 1) { // block
    double d = difftime(now, indexRefresh);
    printf("==>  total time spent to process this block: %4.0f seconds. \n", d);
    indexRefresh = now;
  } else if (stage == 2) { // variable
    //double d = difftime(now, indexRefresh);
    //printf("==>  total time spent to process this variable: %4.0f seconds. \n", d);
    //indexRefresh = now;
  }
  */
}

void onBlock(int rank, ADIOS_FILE* f, ADIOS_VARINFO* v, int i, int j, int blockCounter, FastBitDataType ft)
{
      char bmsVarName[100];
      char keyVarName[100];
      char offsetName[100];

       int64_t       var_ids_bms[v->nblocks[i]];
       int64_t       var_ids_key[v->nblocks[i]];
       int64_t       var_ids_offset[v->nblocks[i]];


       
	 sprintf(bmsVarName, "bms-%d-%d-%d", v->varid, i, j);
	 sprintf(keyVarName, "key-%d-%d-%d", v->varid, i, j);
	 sprintf(offsetName, "offset-%d-%d-%d", v->varid, i, j);

	 //	 blockCounter++;
	 uint64_t blockSize = fastbit_adios_util_getBlockSize(v, blockCounter);
	 uint64_t blockDataByteSize = adios_type_size (v->type, v->value) * blockSize; 

	 char notes[100];
	 logTime(NULL); logTimeMillis(NULL);

	 sprintf(notes, "  reading data from adios  on varid=%d, time=%d, block: %d, size=%ld bytes=%ld", v->varid, i, j, blockSize, blockDataByteSize);

	 logTime(notes); logTimeMillis(notes);
	 localtime(&indexRefresh);

	    //printf("   %d th block / (%d), size= %llu bytes=%llu", j, blockSize, blockCounter, blockDataByteSize);
	    
	    void* data = malloc (blockDataByteSize);
	    ADIOS_SELECTION* blockSel = adios_selection_writeblock(j);

	    //adios_selcton_writeblock(num),  0 <= num <  nblocks[timestep]
	    //ADIOS_SELECTION* blockSel = adios_selection_writeblock(blockCounter);
	    int err = adios_schedule_read_byid(f, blockSel, v->varid, i, 1, data);
	    if (!err) {	
	      err = adios_perform_reads(f, 1);
	    } else {
	      printf("Unable to read block %d at timestep: %d \n", j, i);
	      return;
	      //break;
	    }
	    //fastbit_adios_util_printData(data, v->type, blockSize);
	   
	    double* keys = NULL;
	    int64_t *offsets = NULL;
	    uint32_t *bms = NULL;
	    uint64_t nk=0, no=0, nb=0;
	    
	    const char* datasetName = "test";
	    logTime("  data collected, fastbit start indexing"); 
	    logTimeMillis("  data collected, fastbit start indexing"); 
	    fastbitIndex(datasetName, data, blockSize, ft, &keys, &nk, &offsets, &no, &bms, &nb);
	    logTime("  indexed on block");
	    logTimeMillis("  indexed on block");

	    printf("    rank:%d, index created =  %llu, %llu, %llu, on var:%d, timestep: %d, block %d\n", rank, nb, nk, no, v->varid, i, j);
	    sum_nb += nb; sum_nk += nk, sum_no += no;

	    defineFastbitVar(1,bmsVarName, &var_ids_bms[j], adios_unsigned_integer, &nb,0,0);    			    
	    defineFastbitVar(1, keyVarName, &var_ids_key[j], adios_double, &nk, 0, 0);
	    defineFastbitVar(1, offsetName, &var_ids_offset[j], adios_long, &no, 0, 0); 

	    logTime("  write starts");
	    logTimeMillis("  write starts");
	    adios_write_byid(gAdios_write_file, var_ids_bms[j], bms);
	    adios_write_byid(gAdios_write_file, var_ids_key[j], keys);
	    adios_write_byid(gAdios_write_file, var_ids_offset[j], offsets);
	    logTime("  write ends");
	    logTimeMillis("  write ends");
	    sumLogTime(1);
	    sumLogTimeMillis(1);

	    adios_selection_delete(blockSel);
	    free(data);	 

	    verifyData(f, v, blockCounter, i);
}

void buildIndex_mpi(ADIOS_FILE* f, ADIOS_VARINFO* v, int rank, int size)
{

  adios_inq_var_blockinfo(f, v);

  int i=0;
  int j=0;
  int k=0;

  printf("building index on variable %d, binning op=%s\n", v->varid, gBinningOption);

  int blockCounter = -1;
  FastBitDataType ft = fastbit_adios_util_getFastbitDataType(v->type);

  for (i=0; i < v->nsteps; i++) {
       int nblocks = v->nblocks[i];

       for (j=0; j < v->nblocks[i]; j++) {
	 blockCounter++;
	 if (blockCounter % size == rank) {
	   onBlock(rank, f, v, i, j, blockCounter, ft);
	 }
       }
       printf(" rank %d, varid %d, timestep %d  total index created =  %llu, %llu, %llu, \n", rank, v->varid, i, sum_nb, sum_nk, sum_no);
       printf(" rank %d, varid %d, timestep %d  total bytes         =  %llu, %llu, %llu, \n", rank, v->varid, i,
	      adios_type_size(adios_unsigned_integer , NULL)*sum_nb,
	      adios_type_size(adios_double, NULL)*sum_nk, 
	      adios_type_size(adios_long, NULL)*sum_no);
       fastbit_cleanup();
       printf("\n");
  }
}

void buildIndexOnAllVar(ADIOS_FILE* f, int rank, int size) 
{
  int numVars = f->nvars;
  
  int i=0;
  for (i=0; i<numVars; i++) {
    char* varName = f->var_namelist[i];
    ADIOS_VARINFO* v = adios_inq_var(f, varName);

    if (rank == 0) {
      printf("\n==> building fastbit index on  %dth variable: %s, %s ", i, varName, gBinningOption);
    }

    if (v->ndim > 0) {
      buildIndex_mpi(f, v, rank, size);
    } else {
      if (rank == 0) {
	printf("\t ... skipping scalar ...\n");
      }
    }
    adios_free_varinfo(v);
  }
}

uint64_t estimateBytesOnVar(ADIOS_FILE* f, ADIOS_VARINFO* v) 
{
    if (v->ndim == 0) {
      return 0;
    }

    adios_inq_var_blockinfo(f, v);

    uint64_t result = 0;
    int i = 0; 
    int j = 0;
    for (i=0; i<v->sum_nblocks; i++) {
      ADIOS_VARBLOCK curr = v->blockinfo[i];

      uint64_t blockSize = fastbit_adios_util_getBlockSize(v, i);
      result += blockSize;
      //uint64_t blockDataByteSize = adios_type_size (v->type, v->value) * blockSize0; 

      /*uint64_t blockSize = 1;
      for (j=0; j<v->ndim; j++) {
	blockSize = blockSize* curr.count[j];
      }
      result += blockSize;
      */
    }

    return adios_type_size (v->type, v->value) * result; 
    //int typeSize = adios_type_size(v->type, NULL);
    //return result * typeSize;
}

int64_t getByteEstimationOnFile(ADIOS_FILE* f, int rank) 
{
  // check all vars
  int numVars = f->nvars;
  uint64_t bytes = 0;

  int i=0;
  for (i=0; i<numVars; i++) {
    char* varName = f->var_namelist[i];
    ADIOS_VARINFO* v = adios_inq_var(f, varName);

    uint64_t varSize = estimateBytesOnVar(f, v);
    if (rank == 0) {
      printf(" var: %s has size: %lld\n", varName, varSize);
    }
    bytes += varSize;

    adios_free_varinfo(v);
  }
  return bytes;

}
int64_t getByteEstimation(ADIOS_FILE* f, int rank, int argc, char** argv)
{
  uint64_t bytes = 0;

  if (argc >= 3) {
    int i=2;
    while (i<argc) {
      const char* varName = argv[i];      
      if(strstr(varName, "<binning prec") != NULL) {
	if (gBinningOption == NULL) {
	   gBinningOption = argv[i];
	}
	if (argc == 3) {
	  return getByteEstimationOnFile(f, rank);
	}
	i++;
	continue;
      }
      ADIOS_VARINFO * v = adios_inq_var(f, varName);
      if (v == NULL) {
	printf("Invalid variable: %s\n", varName);	
	return -1;
      }
      uint64_t varSize = estimateBytesOnVar(f, v);
      printf(" var: %s has size: %lld\n", varName, varSize);
      bytes += varSize;
      adios_free_varinfo(v);
      i++;
    } 
    return bytes;
  } else { 
    return getByteEstimationOnFile(f,rank);
  }
}

int main (int argc, char** argv) 
{
  fastbit_init(0);
  fastbit_set_verbose_level(0);

  ADIOS_FILE * f;
  //MPI_Comm    comm_dummy = 0;  // MPI_Comm is defined through adios_read.h 
  MPI_Comm comm_dummy = MPI_COMM_WORLD;


  int         rank, size;
  MPI_Init (&argc, &argv);			   
  MPI_Comm_rank (comm_dummy, &rank);
  MPI_Comm_size (comm_dummy, &size);

  adios_init_noxml (comm_dummy);

  char *idxFileName;
  
  if (argc < 2) {
    printf("Usage: index_fastbit fileName (attrName)");
    return 0;
  }

  f = adios_read_open_file (argv[1], ADIOS_READ_METHOD_BP, comm_dummy);
  if (f == NULL) {
    printf ("::%s\n", adios_errmsg());
    return -1;
  }

  adios_allocate_buffer (ADIOS_BUFFER_ALLOC_NOW, (f->file_size)*2/1048576 + 5); // +5MB for extra room in buffer
  adios_declare_group (&gAdios_group, gGroupNameFastbitIdx, "", adios_flag_yes);
  adios_select_method (gAdios_group, "MPI", "", "");

  
  idxFileName = fastbit_adios_util_getFastbitIndexFileName(argv[1]);

  unlink(idxFileName);
  adios_open (&gAdios_write_file, gGroupNameFastbitIdx, idxFileName, "w", comm_dummy);

  uint64_t adios_totalsize;
  
  uint64_t estimatedbytes = getByteEstimation(f, rank, argc, argv);
  // each processor is writing call this function, so the eventual size would add up 
  adios_group_size (gAdios_write_file, estimatedbytes/size+1048576, &adios_totalsize);     

  printf("=> adios open output file: %s, rank %d allocated %llu bytes... \n", idxFileName, rank, adios_totalsize); 

  sumLogTime(-1);
  sumLogTimeMillis(-1);

  if (argc >= 3) {
     int i=2;
     while (i<argc) {
        const char* varName = argv[i];
	if(strstr(varName, "<binning prec") != NULL) {
	  if (gBinningOption == NULL) {
	    gBinningOption = argv[i];
	  }
	  if (argc == 3) {
	    buildIndexOnAllVar(f, rank, size);
	    break;
	  }
	  i++;
	  continue;
	} else {
	  ADIOS_VARINFO * v = adios_inq_var(f, varName);
	  if (v == NULL) {
	     printf("No such variable: %s\n", varName);
	     return 0;
	   }	
	  printf("building fastbit index on  variable: %s\n", varName);
	  buildIndex_mpi(f, v, rank, size);
	  adios_free_varinfo(v);
	  i++;
	}
     }
  } else {
    buildIndexOnAllVar(f, rank, size);
  }

  sumLogTime(0);
  sumLogTimeMillis(0);

  adios_close(gAdios_write_file);
  adios_read_close(f);

  //
  // writing file clean up
  //


  // read back:
  f = adios_read_open_file (idxFileName, ADIOS_READ_METHOD_BP, comm_dummy);
  if (f == NULL) {
    printf("No such file: %s \n", idxFileName);
    return 0;
  }

  int numVars = f->nvars;
  
  int i=0;
  int k=0;
  int j=0;
  for (i=0; i<numVars; i++) {
      char* varName = f->var_namelist[i];
      ADIOS_VARINFO* v = adios_inq_var(f, varName);

       adios_inq_var_blockinfo(f,v);      
      int timestep = 0;
      for (k=0; k<v->sum_nblocks; k++) {
	  verifyData(f, v, k, timestep);
      }

      adios_free_varinfo(v);
  }

  adios_read_close(f);

  if (rank == 0) {
    printf(" ==>  index file is at: %s\n", idxFileName);
  }

  // clean up
  MPI_Barrier (comm_dummy);
  adios_finalize (rank);
  MPI_Finalize ();
  free (idxFileName);

  fastbit_cleanup();
  return 0;
}


