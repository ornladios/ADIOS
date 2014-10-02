/*
 * ADIOS is freely available under the terms of the BSD license described
 * in the COPYING file in the top level directory of this source distribution.
 *
 * Copyright (c) 2008 - 2009.  UT-BATTELLE, LLC. All rights reserved.
 */

/* ADIOS C Example: read global arrays from a BP file
 *
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>
#include <string.h>
#include <assert.h>
#include "adios_selection.h"
#include "adios_query.h"
#include <mxml.h>
#include <sys/stat.h>

//#include "adios_types.h"
#define MAXDIM    10
#define MAXQUERY  1000

#define GET_ATTR2(n,attr,var,en)                                 \
    if (!strcasecmp (n, attr->name)) {                           \
        if (!var)                                                \
        {                                                        \
            var = attr->value;                                   \
            continue;                                            \
        }                                                        \
        else                                                     \
        {                                                        \
            printf ("xml: duplicate attribute %s on %s (ignored)",n,en); \
            continue;                                            \
        }                                                        \
    }

void tokenize_dimensions2 (const char * str, char *** tokens, int * count)
{
    if (!str)
    {
        *tokens = 0;
        *count = 0;

        return;
    }

    char * save_str = strdup (str);
    char * t = save_str;
    int i;

    trim_spaces (save_str);

    if (strlen (save_str) > 0)
        *count = 1;
    else
    {
        *tokens = 0;
        *count = 0;
        free (save_str);

        return;
    }

    while (*t)
    {
        if (*t == ',')
            (*count)++;
        t++;
    }

    *tokens = (char **) malloc (sizeof (char **) * *count);
    (*tokens) [0] = strdup (strtok (save_str, ","));
    for (i = 1; i < *count; i++)
    {
        (*tokens) [i] = strdup (strtok (NULL, ","));
    }

    free (save_str);
}
//end of stolen functions

// Stack for storing queries
typedef struct {
    int size;
    ADIOS_QUERY *stack[MAXQUERY];
} QueryStack;

// init query stack
void queryStackInit(QueryStack* queryStack)
{
    queryStack->size=0;
}

void queryPush(QueryStack* queryStack, ADIOS_QUERY *q)
{
    if (queryStack->size>=MAXQUERY) {
        printf("Query number exceeds MAXQUERY, exiting\n");
        exit(-1);
    }
    queryStack->stack[queryStack->size++] = q;

}

int queryStackSize(QueryStack* queryStack)
{
    return queryStack->size;
}

ADIOS_QUERY * queryPop(QueryStack* queryStack)
{
    if (queryStackSize(queryStack)==0) {
        printf("Error: popping empty query stack, exiting...\n");
        exit(-1);
    }
    return queryStack->stack[--queryStack->size];
}

// coordinates[2] + coordinates[1] * destcount[2] + coordinates[0]* destcount[1] * destcount[2];
 //coordinates[1] + coordinates[0] * destcount[1] ;
void printRids(const ADIOS_SELECTION_POINTS_STRUCT * pts,  uint64_t *deststart, uint64_t *destcount) {
	uint64_t i = 0, rid=0;
	if (pts->ndim == 3) {
		for (i = 0; i < pts->npoints; i++) {
			rid =  (pts->points[i * 3 + 2] - deststart[2]) + (pts->points[i * 3 + 1] - deststart[1])  * destcount[2] +  (pts->points[i * 3] - deststart[0]) * destcount[2] * destcount[1];
					printf("[ %"PRIu64" ] ,", rid);
		}
	}

	if (pts->ndim == 2) {
		for (i = 0; i < pts->npoints; i++) {
				rid =  (pts->points[i * 2 + 1] - deststart[1]) + (pts->points[i * 2 ] - deststart[0])  * destcount[1];
						printf("[ %"PRIu64" ] ,", rid);
		}
	}
	printf("\n");
}
void printPoints(const ADIOS_SELECTION_POINTS_STRUCT * pts) {
	uint64_t i = 0;
	if (pts->ndim == 3) {

		for (i = 0; i < pts->npoints; i++) {
			printf("[ %"PRIu64", %"PRIu64", %"PRIu64" ] ,", pts->points[i * 3],
					pts->points[i * 3 + 1], pts->points[i * 3 + 2]);
		}

	}
	if (pts->ndim == 2) {
		for (i = 0; i < pts->npoints; i++) {
			printf("[ %"PRIu64", %"PRIu64" ] ,", pts->points[i * 2],
					pts->points[i * 2 + 1]);
		}
	}
	printf("\n");

}

ADIOS_QUERY * createQueryConstraints(ADIOS_FILE* bf, const char* varName,
		ADIOS_SELECTION* box, const char * lb, const char * hb) {
	/* lb <= x <= hb */
	ADIOS_QUERY* q1 = adios_query_create(bf, varName, box, ADIOS_GTEQ, lb);
	ADIOS_QUERY* q2 = adios_query_create(bf, varName, box, ADIOS_LTEQ, hb);
	ADIOS_QUERY* q = adios_query_combine(q1, ADIOS_QUERY_OP_AND, q2);
	return q;

	/* x != lb
	 * return adios_query_create(bf, varName, box, ADIOS_NE, lb);*/
}

#define CHECK_ERROR_DATA(data, num, check) {                     \
		 uint64_t di = 0;                                        \
		 for(di = 0; di < (num); di++){                          \
				if (check)                                       \
					printf("error data: %f, ", (data)[di]);      \
		 }                                                       \
}

/*
 * bounding box selection & writeblock selection
 */
void multiSelection(ADIOS_FILE* bf, const char * b1, const char * b2,
		ADIOS_FILE *dataF, enum ADIOS_SELECTION_TYPE sel_type) {
	if (sel_type == ADIOS_SELECTION_BOUNDINGBOX) {
		printf("\n=============== testing multiple bound box  ===========\n");
	} else if (sel_type == ADIOS_SELECTION_WRITEBLOCK) {
		printf("\n=============== testing multiple writeblock selection ===========\n");
	}
	const char* varName1 = "temp";
	int ndim = 3;
	uint64_t start1[] = { 0, 0, 0 }; // block 0 -> 1st block
	uint64_t count1[] = { 64, 32, 32 };
	ADIOS_SELECTION *box1 = adios_selection_boundingbox(ndim, start1, count1);

	uint64_t start2[] = { 0, 0, 0 }; //block 2 -> 3rd block
	uint64_t count2[] = { 64, 32, 32 };
	ADIOS_SELECTION* box2 = adios_selection_boundingbox(ndim, start2, count2);

	/*
	 uint64_t start3[] = {64, 32, 0}; //block 4 -> 5th block
	 uint64_t count3[] = {64, 10,32};
	 ADIOS_SELECTION* outBox = adios_selection_boundingbox(ndim, start3, count3);
	 */

	ADIOS_VARINFO * dataV = adios_inq_var(dataF, varName1);
	double lb1v = atof(b1);
	double lb2v = atof(b2);
	int timestep = 0;
	adios_query_set_timestep(timestep);
	ADIOS_QUERY* q1, *q2;
	int64_t batchSize = 1000;
    if (sel_type == ADIOS_SELECTION_BOUNDINGBOX){
		q1 = adios_query_create(bf, varName1, box1, ADIOS_GT, b1); // > b1
		q2 = adios_query_create(bf, varName1, box2, ADIOS_LT, b2); // < b2

    }else if (sel_type == ADIOS_SELECTION_WRITEBLOCK){
    	ADIOS_SELECTION *block1= adios_selection_writeblock(0); // block 0 ->  1st block
    	ADIOS_SELECTION* block2 = adios_selection_writeblock(2); //block 2 -> 3rd block
		q1 = adios_query_create(bf, varName1, block1, ADIOS_GT, b1); // > b1
		q2 = adios_query_create(bf, varName1, block2, ADIOS_LT, b2); // < b2
    }
	// if box has a different shape, e.g. different count values, then combine() returns error
	ADIOS_QUERY* q = adios_query_combine(q1, ADIOS_QUERY_OP_AND, q2);

	// box3 is the same shape as other boxes , if it is in different shape from box1/box2, then error
	ADIOS_SELECTION* outBox = box1;
	while (1) {
		ADIOS_SELECTION* currBatch = NULL;
		int hasMore = adios_query_get_selection(q, batchSize, outBox,
				&currBatch);
		if (currBatch == NULL) // there is no results at all
			break;

		const ADIOS_SELECTION_POINTS_STRUCT * retrievedPts =
				&(currBatch->u.points);
		printf("retrieved points %" PRIu64 " \n", retrievedPts->npoints);
//		printPoints(retrievedPts);
//		printRids(retrievedPts, start1, count1);// outbox is box1
		/*double * data = (double *) malloc(
				retrievedPts->npoints * sizeof(double));
		adios_schedule_read_byid(dataF, currBatch, dataV->varid, 0, 1, data);
		adios_perform_reads(dataF, 1);
		CHECK_ERROR_DATA(data, retrievedPts->npoints, (data[di] <= lb1v ));
		free(data);*/
		adios_selection_delete(currBatch);
		if (hasMore == 0) { // there is no left results to retrieve
			break;
		}
	}

	//reset query engine
/*	q->_onTimeStep = -1; // NO_EVAL_BEFORE
	outBox = box2; // switch to different output box
	while (1) {
		ADIOS_SELECTION* currBatch = NULL;
		int hasMore = adios_query_get_selection(q, batchSize, outBox,
				&currBatch);
		if (currBatch == NULL) // there is no results at all
			break;

		const ADIOS_SELECTION_POINTS_STRUCT * retrievedPts =
				&(currBatch->u.points);
		printf("retrieved points %" PRIu64 " \n", retrievedPts->npoints);

		double * data = (double *) malloc(
				retrievedPts->npoints * sizeof(double));
		adios_schedule_read_byid(dataF, currBatch, dataV->varid, 0, 1, data);
		adios_perform_reads(dataF, 1);
		CHECK_ERROR_DATA(data, retrievedPts->npoints, (data[di] >= lb2v ));
		free(data);
		adios_selection_delete(currBatch);
		if (hasMore == 0) { // there is no left results to retrieve
			break;
		}
	}*/

	adios_query_free(q);

	adios_selection_delete(box1);
	adios_selection_delete(box2);
}


void oneDefinedBox(ADIOS_FILE* bf, const char * lb, const char * hb,
		ADIOS_FILE *dataF) {

	printf("\n=============== testing one single bounding box ===========\n");
	// temp data is in the range btw 100 and 200
	const char* varName1 = "temp";
	int ndim = 3;
	uint64_t start1[] = { 0, 0, 0 };
	uint64_t count1[] = { 64, 32, 32 };
	ADIOS_SELECTION *box1 = adios_selection_boundingbox(ndim, start1, count1);
	ADIOS_VARINFO * dataV = adios_inq_var(dataF, varName1);
	double lbv = atof(lb);
	double hbv = atof(hb);
	int timestep = 0;
	adios_query_set_timestep(timestep);
	int64_t batchSize = 10000;

	ADIOS_QUERY * q = createQueryConstraints(bf, varName1, box1, lb, hb);
	//      int64_t estimateResults = adios_query_estimate(q1);
	//      printf("estimated result number %"PRIu64 " \n", estimateResults);

	while (1) {
		ADIOS_SELECTION* currBatch = NULL;
		int hasMore = adios_query_get_selection(q, batchSize, box1, &currBatch);
		if (hasMore == 0) { // there is no left results to retrieve
			break;
		}
		assert(currBatch->type ==ADIOS_SELECTION_POINTS);
		const ADIOS_SELECTION_POINTS_STRUCT * retrievedPts =
				&(currBatch->u.points);
		printf("retrieved points %" PRIu64 " \n", retrievedPts->npoints);

	/*	double * data = (double *) malloc(
				retrievedPts->npoints * sizeof(double));
		adios_schedule_read_byid(dataF, currBatch, dataV->varid, 0, 1, data);
		adios_perform_reads(dataF, 1);
		// lb <= x <= lb
		CHECK_ERROR_DATA(data, retrievedPts->npoints, (data[di] > hbv || data[di] < lbv));*/
		//        CHECK_ERROR_DATA(data, retrievedPts->npoints, (data[di] == lbv));
//		free(data);
		//        printPoints(retrievedPts);
		adios_selection_delete(currBatch);

	}

	adios_query_free(q);
	adios_selection_delete(box1);

}

void retrieveAllValues(ADIOS_FILE* f) {
	uint64_t start[] = { 0 , 0};
	uint64_t count[] = { 8, 8 };
	int dim = 2;
	int c = 0;
	uint64_t totalElm = 1;
	for(c =0; c < dim ; c ++){
		totalElm *= count[c];
	}
	uint64_t* points = (uint64_t*) (malloc(dim * totalElm* sizeof(uint64_t)));
	uint64_t i = 0, j = 0, k = 0;
	for(i = 0; i < count[0]; i ++ ){
		for (j = 0 ; j < count[1]; j ++){
			points[k++] = i;
			points[k++] = j;
		}
	}

	ADIOS_SELECTION* currBatch  = adios_selection_points(dim, totalElm, points);

	const ADIOS_SELECTION_POINTS_STRUCT * retrievedPts =
						&(currBatch->u.points);

	printPoints(retrievedPts);
	void *data = malloc(totalElm * 4);

	int timestep = 1;
	adios_schedule_read (f, currBatch, "temp", timestep, 1, data);
	adios_perform_reads(f, 1);

	for (i = 0; i < totalElm; i++) {
		printf("%.6f\t", ((float*)data)[i]);
	}
	printf("\n");

	free(points);
}

void oneBoundingBoxForVars(ADIOS_FILE* f, ADIOS_FILE *dataF, const char * lbs, const char * hbs) {
	printf("\n=============== test oneBoundingBoxForVars ===========\n");
	uint64_t start[] = { 0, 0 };
	uint64_t count[] = { 8, 8 };
	//uint64_t start[] = { 0, 0, 0 };
	//uint64_t count[] = { 64, 32, 32 };
	//  uint64_t count[] = {32,16,16};
	ADIOS_SELECTION* box = adios_selection_boundingbox(2, start, count);

	const char* varName1 = "temp";
	double lb = atof(lbs);
	double hb = atof(hbs);

//	const char* value1 = "170.0";
//	const char* varName2 = "uvel";
//	const char* value2 = "0.5";
	//const char* value2 = "14";
//	double uvelConstraint = atof(value2);
//	ADIOS_QUERY* q1 = adios_query_create(f, varName1, box, ADIOS_LT, value1); // temp < 150.0
//	ADIOS_QUERY* q2 = adios_query_create(f, varName1, box, ADIOS_GT, value2); // uvel > 15
	//ADIOS_QUERY* q2 = adios_query_create(f, varName2, box, ADIOS_GT, value2); // uvel > 15

	ADIOS_QUERY* q1 = adios_query_create(f, varName1, box, ADIOS_LT, hbs); // temp < hb
	ADIOS_QUERY* q2 = adios_query_create(f, varName1, box, ADIOS_GT, lbs); // temp < lb

//	ADIOS_QUERY* q1 = adios_query_create(f, varName1, box, ADIOS_GT, lbs); // temp < hb
//	ADIOS_QUERY *q = q1;
	ADIOS_QUERY* q = adios_query_combine(q1, ADIOS_QUERY_OP_AND, q2);

	ADIOS_VARINFO * tempVar = adios_inq_var(dataF, varName1);
	//ADIOS_VARINFO * uvelVar = adios_inq_var(dataF, varName2);

	int64_t batchSize = 1220;

	int i = 0, timestep = 0 ;
	printf("times steps for variable is: %d \n", q1->_var->nsteps);
	for (timestep  = 0; timestep  < q1->_var->nsteps; timestep ++) {
		printf("querying on timestep %d \n", timestep );
		adios_query_set_timestep(timestep );

		ADIOS_SELECTION* currBatch = NULL;
		while ( adios_query_get_selection(q, batchSize, box,
				&currBatch)) {

			assert(currBatch->type ==ADIOS_SELECTION_POINTS);
			const ADIOS_SELECTION_POINTS_STRUCT * retrievedPts =
					&(currBatch->u.points);
			printf("retrieved points %" PRIu64 " \n", retrievedPts->npoints);

			printPoints(retrievedPts);

			int elmSize = adios_type_size(tempVar->type, NULL);
			void *data = malloc(retrievedPts->npoints * elmSize);

			// check returned temp data
//			adios_schedule_read_byid(dataF, currBatch, tempVar->varid, timestep , 1, data);
			adios_schedule_read (dataF, currBatch, varName1, timestep , 1, data);
			adios_perform_reads(dataF, 1);

			printf("Total data retrieved:%"PRIu64"\n", retrievedPts->npoints);
			if (tempVar->type == adios_double){
				for (i = 0; i < retrievedPts->npoints; i++) {
					printf("%.6f\t", ((double*)data)[i]);
				}
				printf("\n");
//				CHECK_ERROR_DATA(data, retrievedPts->npoints, (data[di] >= hb && data[di] <=lb));
			}else if (tempVar->type == adios_real){
				for (i = 0; i < retrievedPts->npoints; i++) {
					printf("%.6f\t", ((float*)data)[i]);
				}
				printf("\n");
			}



			// check return uvel data
			/*
			 adios_schedule_read_byid (dataF, currBatch, uvelVar->varid, 0, 1, data);
			 adios_perform_reads(dataF, 1);
			 CHECK_ERROR_DATA(data, retrievedPts->npoints, (data[di] <= uvelConstraint));
			 */
			free(data);
			adios_selection_delete(currBatch);
			currBatch = NULL;

		}

	}

	adios_query_free(q);
	adios_selection_delete(box);
}

int performQuery(ADIOS_FILE *f, ADIOS_QUERY* q, ADIOS_SELECTION* box, char* varname, int totalTS, uint64_t batchSize)
{

    int i = 0, timestep = 0 ;
    ADIOS_VARINFO * tempVar = adios_inq_var(f, varname);
    printf("times steps for variable is: %d, batch size is %llu\n", totalTS, batchSize);
    for (timestep  = 0; timestep  < totalTS; timestep ++) {
    	printf("querying on timestep %d \n", timestep );
    	adios_query_set_timestep(timestep);
    
    	ADIOS_SELECTION* currBatch = NULL;
    	while ( adios_query_get_selection(q, batchSize, box, &currBatch)) {
    
            assert(currBatch->type ==ADIOS_SELECTION_POINTS);
    	    const ADIOS_SELECTION_POINTS_STRUCT * retrievedPts = &(currBatch->u.points);
            printf("retrieved points %" PRIu64 " \n", retrievedPts->npoints);
    
    	    printPoints(retrievedPts);
    
    	    int elmSize = adios_type_size(tempVar->type, NULL);
    	    void *data = malloc(retrievedPts->npoints * elmSize);
    
    	    // check returned temp data
    	    adios_schedule_read_byid(f, currBatch, tempVar->varid, timestep , 1, data);
    	    adios_schedule_read (f, currBatch, varname, timestep , 1, data);
    	    adios_perform_reads(f, 1);
    
    	    printf("Total data retrieved:%"PRIu64"\n", retrievedPts->npoints);
    	    if (tempVar->type == adios_double){
    	    	for (i = 0; i < retrievedPts->npoints; i++) {
    	    		printf("%.6f\t", ((double*)data)[i]);
    	    	}
    	    	printf("\n");
    	    }
            else if (tempVar->type == adios_real){
    	    	for (i = 0; i < retrievedPts->npoints; i++) {
    	    		printf("%.6f\t", ((float*)data)[i]);
    	    	}
    	    	printf("\n");
    	    }
    
    
    	    free(data);
    	    adios_selection_delete(currBatch);
    	    currBatch = NULL;
    
        }

    }

    adios_query_free(q);
}

int parseXml(const char *inputxml, ADIOS_FILE* f)
{
    int i, j;
    FILE * fp = fopen (inputxml,"r");
    if (!fp){
        printf("missing xml input file %s \n", inputxml);
        exit(-1);
    }
    struct stat s;
    char * buffer = NULL;
    if (stat (inputxml, &s) == 0) {
        buffer = malloc (s.st_size + 1);
        buffer [s.st_size] = 0;
    }

    if (buffer)     {
        size_t bytes_read = fread (buffer, 1, s.st_size, fp);

        if (bytes_read != s.st_size) {
            printf("error reading input xml file: %s. Expected %ld Got %ld\n"
                            ,inputxml, s.st_size, bytes_read );
            fclose(fp);
            exit(-1);
        }
    }
    fclose (fp);
    mxml_node_t * doc = NULL;
    mxml_node_t * root = NULL;
    mxml_node_t * queryNode = NULL;
    doc = mxmlLoadString (NULL, buffer, MXML_TEXT_CALLBACK);
    free (buffer);
    buffer = NULL;
    root = doc;

    if (!doc) {
        printf( "unknown error parsing XML (probably structural)\n"
                "Did you remember to start the file with\n"
                "<?xml version=\"1.0\"?>\n");
        exit(-1);
    }
    if (strcasecmp (doc->value.element.name, "adios-alac-test-inputs")) {
        root = mxmlFindElement (doc, doc, "adios-alac-test-inputs", NULL, NULL, MXML_DESCEND_FIRST);
    }

    queryNode = mxmlFindElement(root, root, "query", NULL, NULL, MXML_DESCEND_FIRST);
    const char *numVarS=NULL;
    const char *timestepS=NULL;
    const char *batchsizeS=NULL;

    int numQuery = 0;
    int timestep = 1;
    uint64_t batchsize= 1;
    for (i = 0; i < queryNode->value.element.num_attrs; i++) {
        mxml_attr_t * attr = &queryNode->value.element.attrs [i];
        GET_ATTR2("num",attr,numVarS,"query");
        GET_ATTR2("timestep",attr,timestepS,"query");
        GET_ATTR2("batchsize",attr,batchsizeS,"query");
    }
    if ( !numVarS || !strcmp ( numVarS, "")) {
        printf("missing values for num attribute \n");
        mxmlRelease(doc);
        exit(-1);
    }
    else {
        numQuery  = atoi(numVarS);
        timestep  = atoi(timestepS);
        batchsize = strtoull(batchsizeS, NULL, 10);
    }

    mxml_node_t *outputNode     = NULL;
    const char *outputTypeS=NULL, *outputDimS=NULL, *outputStartS=NULL, *outputCountS=NULL, *outputWbIndexS=NULL;
    int outputDim;
    int outputWbIndex;
    int selType;
    uint64_t outputCount[MAXDIM];
    uint64_t outputStart[MAXDIM];
    char** outputCountTokens=NULL;
    char** outputStartTokens=NULL;
    ADIOS_SELECTION *outputBox;
    
    // Parse output selection info

    outputNode = mxmlFindElement(queryNode, queryNode, "output", NULL, NULL, MXML_DESCEND_FIRST);
    for (i = 0; i < outputNode->value.element.num_attrs; i++) {
        mxml_attr_t * attr = &outputNode->value.element.attrs [i];
        GET_ATTR2("type",attr,outputTypeS,"output");
        if ( strcmp(outputTypeS, "ADIOS_SELECTION_BOUNDINGBOX") == 0) {
            selType = ADIOS_SELECTION_BOUNDINGBOX;
            GET_ATTR2("dim",attr,outputDimS,"output");
            GET_ATTR2("start",attr,outputStartS,"output");
            GET_ATTR2("count",attr,outputCountS,"output");
        }
        else if ( strcmp(outputTypeS, "ADIOS_SELECTION_WRITEBLOCK") == 0) {
            selType = ADIOS_SELECTION_WRITEBLOCK;
            GET_ATTR2("index",attr,outputWbIndexS,"selection");
        }
    }
    if ( selType == ADIOS_SELECTION_BOUNDINGBOX ) {
        if ( !outputTypeS || !outputDimS || !outputStartS || !outputCountS || !strcmp (outputTypeS, "")|| !strcmp (outputDimS, "") || !strcmp (outputStartS, "") || !strcmp (outputCountS, "") ) {
            printf("missing values for output attribute \n");
            mxmlRelease(doc);
            exit(-1);
        }
        else {
            outputDim = atoi(outputDimS);
            if (outputDim > MAXDIM) {
                printf("QueryDim exceeds 10, readjust MAXDIM to larger value, exiting...\n");
                exit(-1);
            }
    
            tokenize_dimensions2(outputStartS, &outputStartTokens, &outputDim);
            tokenize_dimensions2(outputCountS, &outputCountTokens, &outputDim);

            for (j = 0; j < outputDim; j ++){
                outputStart[j] = atoi(outputStartTokens[j]);
                outputCount[j] = atoi(outputCountTokens[j]);
            }

            outputBox = adios_selection_boundingbox(outputDim, outputStart, outputCount);

            printf("Selected output boundingbox: dim:%d start:", outputDim);
            for (j = 0; j < outputDim; j ++){
                printf(" %d", outputStart[j]);
            }
            printf("\t count:");
            for (j = 0; j < outputDim; j ++){
                printf(" %d", outputCount[j]);
            }
            printf("\n");

        }

    }
    else if( selType == ADIOS_SELECTION_WRITEBLOCK ) {

        if ( !outputWbIndexS || !strcmp (outputWbIndexS, "") ) {
                printf("missing values for selection attribute \n");
                mxmlRelease(doc);
                exit(-1);
            }
            else {
                outputWbIndex = atoi(outputWbIndexS);
    	        outputBox = adios_selection_writeblock(outputWbIndex);

                printf("Selected output writeblock: %d\n", outputWbIndex);
            }
    }


    // Iterate all combine/entry nodes in <query>
    mxml_node_t *entryNode      = NULL;
    mxml_node_t *selectionNode  = NULL;
    const char *varNameS=NULL, *opS=NULL, *constraintS=NULL;
    const char *typeS=NULL, *dimS=NULL, *startS=NULL, *countS=NULL, *wbIndexS=NULL;
    int entryIter;
    int queryDim;
    int wbIndex;
    uint64_t queryCount[MAXDIM];
    uint64_t queryStart[MAXDIM];
    ADIOS_SELECTION *box, *block;
    ADIOS_QUERY *q, *q1, *q2, *qc;
    char** queryCountTokens=NULL;
    char** queryStartTokens=NULL;
    char *queryCombineOp=NULL;

    // init query stack
    QueryStack queryStack;
    queryStackInit(&queryStack);

    for (entryIter = 0; entryIter < (numQuery*2-1); entryIter++) {
        
        // Find entry node
        if (entryIter == 0) {
            entryNode = mxmlFindElement(queryNode, queryNode, "entry", NULL, NULL, MXML_DESCEND_FIRST);
        }
        else {
            // this is the only way I found for getting the next <entry> or <combine> node
            entryNode= mxmlWalkNext(entryNode, queryNode, MXML_NO_DESCEND);
            entryNode= mxmlWalkNext(entryNode, queryNode, MXML_NO_DESCEND);
        }
        
        // check if current node is <combine>
        if ( strcmp( (&(entryNode->value.element.attrs[0]))->name, "op") == 0 ) {
            queryCombineOp = (&(entryNode->value.element.attrs[0]))->value;
            printf("Found combine op %s\n", queryCombineOp);
            // pop up two query and perform the op
            if (queryStackSize(&queryStack)<2) {
                printf("Popping with less than 2 queries in query stack, exiting...\n");
                exit(-1);
            }

            q1 = queryPop(&queryStack);
            q2 = queryPop(&queryStack);
            if (strcmp(queryCombineOp, "AND") == 0 || strcmp(queryCombineOp, "and") == 0) {
                qc = adios_query_combine(q1, ADIOS_QUERY_OP_AND, q2);
            }
            else if (strcmp(queryCombineOp, "OR") == 0 || strcmp(queryCombineOp, "or") == 0) {
                qc = adios_query_combine(q1, ADIOS_QUERY_OP_OR, q2);
            }
            queryPush(&queryStack,qc);

            //adios_query_free(q1);
            //adios_query_free(q2);

            continue;
        }
        

        // Make sure all *S are NULL for verification
        varNameS=NULL, opS=NULL, constraintS=NULL;
        typeS=NULL, dimS=NULL, startS=NULL, countS=NULL, wbIndexS=NULL;

        for (i = 0; i < entryNode->value.element.num_attrs; i++) {
            mxml_attr_t * attr = &entryNode->value.element.attrs [i];
            GET_ATTR2("var",attr,varNameS,"entry");
            GET_ATTR2("op",attr,opS,"entry");
            GET_ATTR2("constraint",attr,constraintS,"entry");
        }
        if ( !varNameS || !opS || !constraintS || !strcmp (varNameS, "")|| !strcmp (opS, "") || !strcmp (constraintS, "") ) {
            printf("missing values for entry attribute \n");
            mxmlRelease(doc);
            exit(-1);
        }

        // Parse selection 
        selectionNode = mxmlFindElement(entryNode, entryNode, "selection", NULL, NULL, MXML_DESCEND_FIRST);
 
        for (i = 0; i < selectionNode->value.element.num_attrs; i++) {
            mxml_attr_t * attr = &selectionNode->value.element.attrs [i];
            GET_ATTR2("type",attr,typeS,"selection");
            if ( strcmp(typeS, "ADIOS_SELECTION_BOUNDINGBOX") == 0) {
                selType = ADIOS_SELECTION_BOUNDINGBOX;
                GET_ATTR2("dim",attr,dimS,"selection");
                GET_ATTR2("start",attr,startS,"selection");
                GET_ATTR2("count",attr,countS,"selection");
            }
            else if ( strcmp(typeS, "ADIOS_SELECTION_WRITEBLOCK") == 0) {
                selType = ADIOS_SELECTION_WRITEBLOCK;
                GET_ATTR2("index",attr,wbIndexS,"selection");
            }
        }
        if ( selType == ADIOS_SELECTION_BOUNDINGBOX ) {

            if ( !typeS || !dimS || !startS || !countS || !strcmp (typeS, "")|| !strcmp (dimS, "") || !strcmp (startS, "") || !strcmp (countS, "") ) {
                printf("missing values for selection attribute \n");
                mxmlRelease(doc);
                exit(-1);
            }
            else {
                queryDim = atoi(dimS);
                if (queryDim > MAXDIM) {
                    printf("QueryDim exceeds 10, readjust MAXDIM to larger value, exiting...\n");
                    exit(-1);
                }
        
                tokenize_dimensions2(startS, &queryStartTokens, &queryDim);
                tokenize_dimensions2(countS, &queryCountTokens, &queryDim);

                for (j = 0; j < queryDim; j ++){
                    queryStart[j] = atoi(queryStartTokens[j]);
                    queryCount[j] = atoi(queryCountTokens[j]);
                }

                box = adios_selection_boundingbox(queryDim, queryStart, queryCount);

                if( strcmp(opS, "<=") == 0 )
                    q = adios_query_create(f, varNameS, box, ADIOS_LTEQ, constraintS);
                else if( strcmp(opS, "<") == 0 )
                    q = adios_query_create(f, varNameS, box, ADIOS_LT, constraintS);
                else if( strcmp(opS, ">=") == 0 )
                    q = adios_query_create(f, varNameS, box, ADIOS_GTEQ, constraintS);
                else if( strcmp(opS, ">") == 0 )
                    q = adios_query_create(f, varNameS, box, ADIOS_GT, constraintS);
                else {
                    printf("Unsupported entry op %s\n", opS);
                    exit(-1);
                }

                queryPush(&queryStack,q);

                printf("Selected input bounding box:  dim:%d start:", queryDim);
                for (j = 0; j < queryDim; j ++){
                    printf(" %d", queryStart[j]);
                }
                printf("\t count:");
                for (j = 0; j < queryDim; j ++){
                    printf(" %d", queryCount[j]);
                }
                printf("\n");



            }
        } // selType == ADIOS_SELECTION_BOUNDINGBOX
        else {

            if ( !wbIndexS || !strcmp (wbIndexS, "") ) {
                printf("missing values for selection attribute \n");
                mxmlRelease(doc);
                exit(-1);
            }
            else {
                wbIndex = atoi(wbIndexS);

    	        block   = adios_selection_writeblock(wbIndex);
                if( strcmp(opS, "<=") == 0 )
                    q = adios_query_create(f, varNameS, block, ADIOS_LTEQ, constraintS);
                else if( strcmp(opS, "<") == 0 )
                    q = adios_query_create(f, varNameS, block, ADIOS_LT, constraintS);
                else if( strcmp(opS, ">=") == 0 )
                    q = adios_query_create(f, varNameS, block, ADIOS_GTEQ, constraintS);
                else if( strcmp(opS, ">") == 0 )
                    q = adios_query_create(f, varNameS, block, ADIOS_GT, constraintS);
                else {
                    printf("Unsupported entry op %s\n", opS);
                    exit(-1);
                }

                queryPush(&queryStack,q);
                printf("Selected input writeblock: %d\n", wbIndex);

            }

        }
       
   
        printf("Parsed entry: var=%s op=%s constraint=%s\n", varNameS, opS, constraintS);
    
    
    }

    // TODO: need some correct checking for box and varNameS
    performQuery(f, queryPop(&queryStack), outputBox, varNameS, timestep, batchsize);

}

int main(int argc, char ** argv) {

	char filename[256];
	char dataFileName[256];
        char xmlFileName[256];

	int i, j, datasize, if_any;
	MPI_Comm comm = MPI_COMM_WORLD;
	enum ADIOS_READ_METHOD method = ADIOS_READ_METHOD_BP;
	ADIOS_SELECTION * sel1, *sel2;
	ADIOS_VARCHUNK * chunk = 0;
	double * data = NULL;
	uint64_t start[2], count[2], npoints, *points;
	MPI_Init(&argc, &argv);
	if (argc != 3) {
		printf(" usage: %s {input bp file} {xml file}\n", argv[0]);
		//printf(" usage: %s {input bp file},  {bp file without transform} {lb} {hb} \n",	argv[0]);
		return 1;
	}
        else {
	    strcpy(xmlFileName,  argv[2]);
        }

	/* char lbstr[255], hbstr[255]; */

/* 	if (argc == 3) { */
/* 		strcpy(dataFileName, argv[2]); */
/* 		strcpy(xmlFileName,  argv[3]); */

/* 	} else { */
/* 		strcpy(dataFileName, "./xml/alacrity-2var-no-transform_524288.bp"); */
/* 	} */

	/* argc >= 4 ? strcpy(lbstr, argv[3]) : strcpy(lbstr, "0.0"); */
	/* argc >= 5 ? strcpy(hbstr, argv[4]) : strcpy(hbstr, "0.0"); */

	adios_read_init_method(method, comm, NULL);

	ADIOS_FILE * f = adios_read_open_file(argv[1], method, comm);

	if (f == NULL) {
		MPI_Finalize();
		printf(" can not open file %s \n", argv[1]);
		return 1;
	}

	/* ADIOS_FILE * dataF = adios_read_open_file(dataFileName, method, comm); */

	/* if (dataF == NULL) { */
	/* 	if (f) { */
	/* 		adios_read_close(f); */
	/* 		MPI_Finalize(); */
	/* 	} */
	/* 	printf(" can not open file %s \n", dataFileName); */
	/* 	return 1; */
	/* } */

	adios_query_init(ADIOS_QUERY_TOOL_ALACRITY);

	//====================== start to test ===================//
//	    oneDefinedBox(f , lbstr, hbstr, dataF);
//	multiSelection(f , lbstr, hbstr, dataF, ADIOS_SELECTION_BOUNDINGBOX);

//	multiSelection(f , lbstr, hbstr, dataF, ADIOS_SELECTION_WRITEBLOCK);

        parseXml(xmlFileName, f);
        //parseXml(xmlFileName, f, dataF);	

        //char *lbstr="0.9";
        //char *hbstr="0.9";
	//oneBoundingBoxForVars(f, dataF, lbstr, hbstr);
//	printf("retrieve all values \n");
//	retrieveAllValues(f);

	adios_read_close(f);
	/* adios_read_close(dataF); */

	adios_read_finalize_method(ADIOS_READ_METHOD_BP);

	MPI_Finalize();
	return 0;
}
