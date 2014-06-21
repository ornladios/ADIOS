/*
 * adios_read_ext_test.c
 *
 *  Created on: Jun 15, 2014
 *      Author: xczou
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <inttypes.h>
#include <assert.h>
#include "mpi.h"
#include "adios_read.h"
#include "adios_read_ext.h"

/*
 * RUN with one processor is enough
 */
int main (int argc, char ** argv)
{

    char        filename [256];
    int         i, j, datasize, if_any;
    MPI_Comm    comm = MPI_COMM_WORLD;
    enum ADIOS_READ_METHOD method = ADIOS_READ_METHOD_BP;
    ADIOS_SELECTION * sel1, * sel2;
    ADIOS_VARCHUNK * chunk = 0;
    double * data = NULL;
    uint64_t start[2], count[2], npoints, * points;

    MPI_Init (&argc, &argv);

    if (argc < 2 ){
    	printf(" usage: %s {input bp file} \n", argv[0]);
    	return 1;
    }
    adios_read_init_method (method, comm, NULL);

    ADIOS_FILE * fp = adios_read_open_file (argv[1], method, comm);

    if ( fp == NULL){
    	printf(" can not open file %s \n", argv[1]);
    	return 1;
    }

    printf("adios file : ");
    printf("nvar %d", fp->nvars);
    for(i=0; i < fp->nvars; i++) {
    	printf(" { %s }", fp->var_namelist[i]);
    }
    printf(" || nattrs %d", fp->nattrs);
    for(i=0; i < fp->nattrs; i++) {
       	printf(" { %s }", fp->attr_namelist[i]);
    }
    printf("\n");

    char varName[256] = "rdm";
	ADIOS_VARINFO* v = adios_inq_var(fp, varName);
	if (v == NULL) {
		printf(" Error! no such var:%s \n", varName);
		return 1;
	}

	printf("var info: ");
	printf("varId: %d, type: %d, ndim: %d, nsteps: %d, nblocks: %d, ",
			v->varid, v->type, v->ndim, v->nsteps, v->sum_nblocks);
	printf("dims: {");
	for(i = 0; i < v->ndim; i++){
		printf("%" PRIu64 ",", v->dims[i]);
	}
	printf("}\n");

	i = adios_inq_var_blockinfo(fp, v);
	if (i != 0){
		printf("error from adios_inq_var_blockinfo \n");
		return 1;
	}

	ADIOS_VARBLOCK * bi = v->blockinfo;
	printf("block info: ");
	for(i=0; i< v->sum_nblocks; i++){
		printf("{ ");
		for(j = 0; j < v->ndim; j ++){
			printf("%" PRIu64 ": %"PRIu64 " , ", bi[i].start[j], bi[i].count[j]);
		}
		printf(" }");
	}
	printf("\n");


    data_view_t  dv = PHYSICAL_DATA_VIEW;
    adios_read_set_data_view(fp, dv);
    int blockId = 0;
    //======For this example, it only has one block =========/
	ADIOS_VARTRANSFORM * tfv = adios_inq_var_transform(fp, v);
	ADIOS_TRANSFORM_METADATA * tmetas = tfv->transform_metadatas;
	ADIOS_TRANSFORM_METADATA tmeta = tmetas[blockId];
	assert(tmeta.length == 24);
	uint64_t * threeData = (uint64_t *) tmeta.content;
	printf("meta size: %" PRIu64 ", index size: %" PRIu64
			", data size: %" PRIu64 "\n"
			, threeData[0], threeData[1], threeData[2]);

    /*  ADIOS_VARINFO * varinfo = adios_inq_var (f, "temperature");
    if (varinfo)
    {
        int nranks;

        assert(varinfo->ndim == 2);

        nranks = varinfo->dims[0];
        assert(nranks % 4 == 0);
        assert(varinfo->dims[1] == 10);

        datasize = (nranks / 2) * varinfo->dims[1] * sizeof(double);
        data = malloc (datasize);

        start[0] = nranks / 4;
        start[1] = 2;
        count[0] = nranks / 2;
        count[1] = 6;

        sel1 = adios_selection_boundingbox (varinfo->ndim, start, count);

        adios_schedule_read (f, sel1, "temperature", 0, 1, data);
        adios_perform_reads (f, 1);

        printf("Subvolume at (%llu,%llu) of size (%llu,%llu):\n", start[0], start[1], count[0], count[1]);
        for (i = 0; i < count[0]; i++) {
            printf("[ ");
            for (j = 0; j < count[1]; j++) {
                printf("%.0lf ", data[i * count[1] + j]);
            }
            printf("]\n");
        }

        adios_selection_delete (sel1);
    }

    adios_free_varinfo (varinfo);*/
    adios_read_close (fp);

    adios_read_finalize_method (ADIOS_READ_METHOD_BP);


    MPI_Finalize ();
    return 0;
}

