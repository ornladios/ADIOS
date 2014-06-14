#include "adios_query.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "adios_read.h"
#include <iapi.h>

void queryDetail(ADIOS_QUERY* q, int timeStep);

void getCoordinateFromPoints(uint64_t pos,
		const ADIOS_SELECTION_POINTS_STRUCT* sel, uint64_t* coordinates) {
	int i = 0;
	for (i = 0; i < sel->ndim; i++) {
		(coordinates)[i] = sel->points[pos * (sel->ndim) + i];
	}
}

//check point coordinates
//offset in the bounding box needs to be taken account
void getCoordinateFromBox(uint64_t pos,
		const ADIOS_SELECTION_BOUNDINGBOX_STRUCT* sel, int n,
		uint64_t* coordinates) {
	//printf("pos = %ld, n=%d \n", pos, n);
	if (n == 1) {
		coordinates[n - 1] = pos + sel->start[n - 1];
		return;
	}

	uint64_t lastDimSize = sel->count[n - 1];
	uint64_t res = pos % lastDimSize;

	//printf("      lastDim = %ld, res=%d \n", lastDimSize, res);
	coordinates[n - 1] = res + sel->start[n - 1];
	uint64_t stepUp = (pos - res) / lastDimSize;

	//printf("      coordinate[%d]=%ld\n", n-1, coordinates[n-1]);

	getCoordinateFromBox(stepUp, sel, n - 1, coordinates);
}

void getCoordinateFromVariable(uint64_t pos, const ADIOS_VARINFO* var, int n,
		uint64_t* coordinates) {
	if (n == 1) {
		coordinates[n - 1] = pos;
		return;
	}

	uint64_t lastDimSize = var->dims[n - 1];
	uint64_t res = pos % lastDimSize;

	coordinates[n - 1] = res;
	uint64_t stepUp = (pos - res) / lastDimSize;

	getCoordinateFromVariable(stepUp, var, n - 1, coordinates);
}

ADIOS_VARINFO* getAdiosVariable(ADIOS_FILE* f, const char* varName) {
	ADIOS_VARINFO * v = adios_inq_var(f, varName);

	if (v != NULL) {
		printf("  found variable [%s] in file\n", varName);
		return v;
	}

	return NULL;
}

FastBitDataType getFastbitDataType(enum ADIOS_DATATYPES type) {
	switch (type) {
	case adios_unsigned_byte:
		return FastBitDataTypeUByte;
		break;

	case adios_byte:
		return FastBitDataTypeByte;
		break;

	case adios_short:
		return FastBitDataTypeShort;
		break;

	case adios_unsigned_short:
		return FastBitDataTypeUShort;
		break;

	case adios_integer:
		return FastBitDataTypeInt;
		break;

	case adios_unsigned_integer:
		return FastBitDataTypeUInt;
		break;

	case adios_long:
		return FastBitDataTypeLong;
		break;

	case adios_unsigned_long:
		return FastBitDataTypeULong;
		break;

	case adios_string:
		return FastBitDataTypeUnknown;
		break;

	case adios_real:
		return FastBitDataTypeFloat;
		break;

	case adios_double:
		return FastBitDataTypeDouble;
		break;

	case adios_long_double:
		//sprintf (s, "%Lg", ((long double *) data)[idx]);
	case adios_complex:
		//sprintf (s, "(%g, %g)", ((float *) data)[2*idx], ((float *) data)[2*idx+1]);
	case adios_double_complex:
		//sprintf (s, "(%lg, %lg)", ((double *) data)[2*idx], ((double *) data)[2*idx+1]);
		return FastBitDataTypeDouble;
	}

}

FastBitCompareType getFastbitCompareType(enum ADIOS_PREDICATE_MODE op) {
	switch (op) {
	case ADIOS_LT:
		return FastBitCompareLess;
		break;
	case ADIOS_LTEQ:
		return FastBitCompareLessEqual;
		break;
	case ADIOS_GT:
		return FastBitCompareGreater;
		break;
	case ADIOS_GTEQ:
		return FastBitCompareGreaterEqual;
		break;
	case ADIOS_EQ:
		return FastBitCompareEqual;
		break;
	case ADIOS_NE:
		return FastBitCompareNotEqual;
		break;
	}
}

//
//
//
void adios_query_fastbit_init_method() {
	const char* conffile = 0;
#ifdef DEBUG
	int msglvl = 200;
#else
	int msglvl = 0;
#endif
	fastbit_init(conffile);
	fastbit_set_verbose_level(msglvl);

	printf("[fastbit has initialized with msglvl = %d]\n", msglvl);
}

void assertValue(char* input, char* endptr) {
	if (*endptr != '\0')
		if ((errno == ERANGE) || (errno != 0) || (endptr == input) || (*endptr
				!= '\0')) {
			//perror("strtol");
			printf("Exit due to :invalid integer value: %s\n", input);
			exit(EXIT_FAILURE);
		}
}

int readWithTimeStep(ADIOS_QUERY* q, int timeStep) {
	FastBitDataType dataType = getFastbitDataType(q->_var->type);
	FastBitCompareType compareOp = getFastbitCompareType(q->_op);

	int errorCode = adios_schedule_read_byid(q->_f, q->_sel, q->_var->varid,
			timeStep, 1, q->_dataSlice);
	printf("      schedule read error code = %d adios_error=%d \n", errorCode,
			adios_errno);
	if (errorCode != 0) {
		return errorCode;
	}

	adios_perform_reads(q->_f, 1); // return 0 regardless whether data is valid, so donnt need to check return value
	printf("      perfo read error code = %d adios_errno=%d\n", errorCode,
			adios_errno);
	if (adios_errno != 0) {
		return -1;
	}

	uint64_t dataSize = q->_rawDataSize;
	int j;
	printf("::\t %s At timestep: %llu datasize=%llu \n\t\t   data:  [",
			q->_condition, timeStep, dataSize);
	for (j = 0; j < dataSize; j++) {
		if (j < 64) {
			if ((j % 10) == 0) {
				printf(" \n\t\t           ");
			}
			if (dataType == FastBitDataTypeDouble) {
				printf(" %lg ", ((double *) (q->_dataSlice))[j]);
			} else if (dataType == FastBitDataTypeFloat) {
				printf(" %g ", ((float *) (q->_dataSlice))[j]);
			} else if (dataType == FastBitDataTypeUInt) {
				printf("%d ", ((uint32_t *) (q->_dataSlice))[j]);
			} else {
				//printf("\t%g ", ((uint32_t *)(q->_dataSlice))[j]);
				printf(" *  ");
			}
		} else {
			printf(" ... ");
			break;
		}
	}
	printf("]\n");

	//adios_free_varinfo(q->_var);

	//q->_queryInternalTimeStep = timeStep;
	char* endptr;
	if (dataType == FastBitDataTypeDouble) {
		double vv = strtod(q->_value, &endptr);
		assertValue(q->_value, endptr);
		q->_queryInternal = fastbit_selection_create(dataType, q->_dataSlice,
				dataSize, compareOp, &vv);
	} else if ((dataType == FastBitDataTypeInt) || (dataType
			== FastBitDataTypeLong) || (dataType == FastBitDataTypeUInt)
			|| (dataType == FastBitDataTypeULong)) {
		long vv = strtol(q->_value, &endptr, 10);

		assertValue(q->_value, endptr);
		q->_queryInternal = fastbit_selection_create(dataType, q->_dataSlice,
				dataSize, compareOp, &vv);
	} else if (dataType == FastBitDataTypeFloat) {
		float vv = strtof(q->_value, &endptr);
		assertValue(q->_value, endptr);
		q->_queryInternal = fastbit_selection_create(dataType, q->_dataSlice,
				dataSize, compareOp, &vv);
	} else {
		q->_queryInternal = NULL;
	}

	return 0;
}

int prepareData(ADIOS_QUERY* q, int timeStep) {
	if (q->_onTimeStep == timeStep) {
		printf("::\t query data has been read for timestep: %d\n", timeStep);
		return 0;
	}

	if (q->_var != NULL) {
		if (q->_queryInternal != NULL) {
			fastbit_selection_free(q->_queryInternal);
		}
		int errorCode = readWithTimeStep(q, timeStep);
		if (errorCode != 0) {
			return errorCode;
		}
	} else {
		int errorCode1 = prepareData((ADIOS_QUERY*) (q->_left), timeStep);
		if (errorCode1 != 0) {
			return errorCode1;
		}

		int errorCode2 = prepareData((ADIOS_QUERY*) (q->_right), timeStep);
		if (errorCode2 != 0) {
			return errorCode2;
		}

		ADIOS_QUERY* _left = (ADIOS_QUERY*) (q->_left);
		ADIOS_QUERY* _right = (ADIOS_QUERY*) (q->_right);

		if (q->_leftToRightOp == ADIOS_QUERY_OP_AND) {
			q->_queryInternal = fastbit_selection_combine(
					_left->_queryInternal, FastBitCombineAnd,
					_right->_queryInternal);
		} else {
			q->_queryInternal = fastbit_selection_combine(
					_left->_queryInternal, FastBitCombineOr,
					_right->_queryInternal);
		}
	}

	q->_onTimeStep = timeStep;
	q->_maxResultDesired = 0;
	q->_lastRead = 0;

	return 0;
}

int64_t adios_query_fastbit_estimate_method(ADIOS_QUERY* q, int timeStep) {
	// call fastbit_estimate_num_hits(selection)
	int errorCode = prepareData(q, timeStep);
	if (errorCode != 0) {
		return -1;
	}
	return fastbit_selection_estimate(q->_queryInternal);
}

int64_t adios_query_fastbit_evaluate_method(ADIOS_QUERY* q, int timeStep,
		uint64_t _maxResult) {
	int errorCode = prepareData(q, timeStep);
	if (errorCode != 0) {
		return -1;
	}

	if (q->_maxResultDesired > 0) { // evaluated already
		if (_maxResult <= 0) { // stay put
			return q->_maxResultDesired;
		}

		if (q->_maxResultDesired > _maxResult) {
			q->_maxResultDesired = _maxResult;
			return q->_maxResultDesired;
		}
		printf(":: user required more results. will evaluate again. \n");
	}

	int64_t numHits = fastbit_selection_evaluate(q->_queryInternal);
	printf(":: ==> fastbit_evaluate() num of hits found for [%s] = %llu\n",
			q->_condition, numHits);
	if (numHits <= _maxResult) {
		// take it as max
		q->_maxResultDesired = numHits;
	} else if (_maxResult > 0) {
		q->_maxResultDesired = _maxResult;
	} else {
		q->_maxResultDesired = numHits;
	}

	return q->_maxResultDesired;
}

void printOneSpatialCoordinate(int dim, uint64_t* spatialCoordinates) {
	int k;
	printf(" spatial = [");
	for (k = 0; k < dim; k++) {
		printf("%d ", spatialCoordinates[k]);
	}
	printf("]\n");

}

void fillUp(int dimSize, uint64_t* spatialCoordinates, uint64_t i,
		uint64_t* pointArray) {
	int k = 0;

	for (k = 0; k < dimSize; k++) {
		uint64_t idx = i * dimSize + k;
		printf(" points[%d] = %lld ", idx, spatialCoordinates[k]);
		pointArray[idx] = spatialCoordinates[k];
	}
	printf("\n");
}

ADIOS_SELECTION* getSpatialCoordinatesDefault(ADIOS_VARINFO* var,
		uint64_t* coordinates, uint64_t retrivalSize) {
	uint64_t arraySize = retrivalSize * (var->ndim);
	uint64_t* pointArray = (uint64_t*) (malloc(arraySize * sizeof(uint64_t)));

	int i;
	for (i = 0; i < retrivalSize; i++) {
		uint64_t spatialCoordinates[var->ndim];
		getCoordinateFromVariable(coordinates[i], var, var->ndim,
				spatialCoordinates);

		fillUp(var->ndim, spatialCoordinates, i, pointArray);
	}
	free(pointArray);
	return adios_selection_points(var->ndim, retrivalSize, pointArray);

}

ADIOS_SELECTION* getSpatialCoordinates(ADIOS_SELECTION* outputBoundry,
		uint64_t* coordinates, uint64_t retrivalSize) {
	int k = 0;
	uint64_t i = 0;

	switch (outputBoundry->type) {
	case ADIOS_SELECTION_BOUNDINGBOX: {
		const ADIOS_SELECTION_BOUNDINGBOX_STRUCT *bb = &(outputBoundry->u.bb);

		uint64_t arraySize = retrivalSize * (bb->ndim);
		uint64_t* pointArray =
				(uint64_t*) (malloc(arraySize * sizeof(uint64_t)));

		for (i = 0; i < retrivalSize; i++) {
			uint64_t spatialCoordinates[bb->ndim];
			getCoordinateFromBox(coordinates[i], bb, bb->ndim,
					spatialCoordinates);

			fillUp(bb->ndim, spatialCoordinates, i, pointArray);
		}
		free(pointArray);
		return adios_selection_points(bb->ndim, retrivalSize, pointArray);
	}
		break;
	case ADIOS_SELECTION_POINTS: {
		const ADIOS_SELECTION_POINTS_STRUCT *points =
				&(outputBoundry->u.points);

		uint64_t arraySize = retrivalSize * (points->ndim);
		uint64_t* pointArray =
				(uint64_t*) (malloc(arraySize * sizeof(uint64_t)));

		for (i = 0; i < retrivalSize; i++) {
			uint64_t spatialCoordinates[points->ndim];
			getCoordinateFromPoints(coordinates[i], points, spatialCoordinates);
			fillUp(points->ndim, spatialCoordinates, i, pointArray);
			/*
			 for (k=0; k<points->ndim; k++) {
			 uint64_t idx = i*(points->ndim)+k;
			 printf(" points[%d] = %lld \n", idx, spatialCoordinates[k]);
			 pointArray[idx] = spatialCoordinates[k];
			 }
			 */
		}
		free(pointArray);
		return adios_selection_points(points->ndim, retrivalSize, pointArray);

		//printOneSpatialCoordinate(points->ndim, spatialCoordinates);
	}
		break;
	default:
		printf("Error: Type of selection is not supported!");
	}
}

ADIOS_QUERY* getFirstLeaf(ADIOS_QUERY* q) {
	if (q == NULL) {
		return NULL;
	}

	if (q->_var != NULL) {
		return q;
	}
	return getFirstLeaf(q->_left);
}

int adios_query_fastbit_get_selection_method(ADIOS_QUERY* q,
		uint64_t batchSize, ADIOS_SELECTION* outputBoundry,
		ADIOS_SELECTION** result) {
	/*
	 if (q->_onTimeStep < 0) {
	 printf(":: Error: need to call evaluate first! Exit.\n");
	 return -1;
	 }
	 */
	adios_query_fastbit_evaluate_method(q, gCurrentTimeStep, 0);
	printf("::\t max=%llu _lastRead=%llu\n", q->_maxResultDesired, q->_lastRead);
	uint64_t retrivalSize = q->_maxResultDesired - q->_lastRead;
	if (retrivalSize == 0) {
		result = 0;
		printf(":: ==> no more results to fetch\n");
		return 0;
	}

	if (retrivalSize > batchSize) {
		retrivalSize = batchSize;
	}

	uint64_t coordinates[retrivalSize];

	fastbit_selection_get_coordinates(q->_queryInternal, coordinates,
			retrivalSize, q->_lastRead);

	q->_lastRead += retrivalSize;

	if (outputBoundry == 0) {
		if ((getFirstLeaf(q) == NULL) || (getFirstLeaf(q)->_var == NULL)) {
			printf(":: Error: unable to get a valid first leaf! Exit. \n");
			return -1;
		}
		*result = getSpatialCoordinatesDefault(getFirstLeaf(q)->_var,
				coordinates, retrivalSize);
	} else {
		*result = getSpatialCoordinates(outputBoundry, coordinates,
				retrivalSize);
	}
	// print results
	int i = 0;
	printf("\n:: coordinates: [\n");
	for (i = 0; i < retrivalSize; i++) {
		//printf("%lld ", coordinates[i]);
		if (i < 100) {
			printf("%lld ", coordinates[i]);
		} else {
			break;
		}
	}
	printf("]\n\n");

	if (q->_lastRead == q->_maxResultDesired) {
		return 0;
	} else {
		return 1;
	}
}

//void adios_fastbit_free_query(ADIOS_QUERY* query) 
int adios_query_fastbit_free_method(ADIOS_QUERY* query) {
	if (query == NULL) {
		return;
	}

	printf(":: free %s\n", query->_condition);
	free(query->_value);
	free(query->_dataSlice);
	free(query->_condition);

	//adios_selection_delete(query->_sel);
	adios_free_varinfo(query->_var);

	//fastbit_selection_free(query->_queryInternal);
	free(query);

}

void adios_query_fastbit_clean_method() {
	fastbit_iapi_free_all();
	fastbit_cleanup();
}

void assert(void* ptr, const char* notes) {
	if (ptr == NULL) {
		printf("::Error allocating memory: %s\n", notes);
		exit(EXIT_FAILURE);
	}
}

/*

 */

void usage(char* prog) {
	printf("Usage: %s <BP-file> [query]\n e.g. ./test_v1 my.bp \"x1 > 10\" \n",
			prog);
}

void getVarName(const char* sliceStr, char** varName, char** dimDef) {
	//char str[strlen(sliceStr)];
	//strncpy(str, sliceStr, strlen(sliceStr));
	char* str = strdup(sliceStr);
	//str[strlen(sliceStr)]=0;
	char* rest = strtok(str, "[");

	*varName = strdup(rest);

	rest = strtok(NULL, "]");

	if (rest == NULL) {
		//strcpy(dimDef, "");
		*dimDef = NULL;
	} else {
		*dimDef = strdup(rest);
	}

	free(str);
}

/*
 */

//
// dimDef is asusmed to be: start0 count0, start1 count1, .. , startN countN
//
void createBox(ADIOS_VARINFO* v, char* dimDef, uint64_t* start, uint64_t* count) {
	// assigns default:
	int i;

	for (i = 0; i < v->ndim; i++) {
		start[i] = 0;
		count[i] = v->dims[i];
	}

	if ((dimDef == NULL) || (strlen(dimDef) == 0)) {
		printf("::\t There is no restriction on dimention.\n");
		//use default
		return;
	}

	const char* comma = ",";

	// int column=-1, columnStarts=-1, columnEnds=0;

	char *dimSpecStart, *end;

	dimSpecStart = end = dimDef;
	i = 0;

	while (*end) {
		switch (*end) {
		case ',':
			*end = '\0';
			//*column = atoi(dimSpecStart);

			count[i] = atol(dimSpecStart);
			dimSpecStart = end + 1;
			*end = ',';
			i++;
			break;
		case ':':
			*end = '\0';
			//*columnStarts=atoi(dimSpecStart);
			start[i] = atol(dimSpecStart);
			*end = ':';
			dimSpecStart = end + 1;
			//break;
		}
		end++;
	}

	count[i] = atol(dimSpecStart);

	for (i = 0; i < v->ndim; i++) {
		if (count[i] <= 0) {
			count[i] = v->dims[i] - start[i];
		}
	}
	/*
	 if (*column ==-1) {
	 *column = atoi(dimSpecStart);
	 } else if (*columnStarts >= 0) {
	 *columnEnds = atoi(dimSpecStart);
	 }
	 */

}
/*
 void createBox(ADIOS_VARINFO* v, char* dimDef, uint64_t* start, uint64_t* count)
 {
 int column=-1, columnStarts=-1, columnEnds=0;
 parseSlice(dimDef, &column, &columnStarts, &columnEnds);
 printf(" column= %d, %d:%d \n\n", column, columnStarts, columnEnds);

 int j;
 for (j = 0; j < v->ndim; j++) {
 if (column == -1) {
 start[j] = 0;
 count[j] = v->dims[j];
 } else if (j == column) {

 if (columnStarts > 0) {
 start[j] = columnStarts;
 } else {
 start[j] = 0;
 }
 if (columnEnds > 0) {
 count[j] = columnEnds - start[j];
 } else {
 count[j] = v->dims[j] - start[j];
 }

 //start[j] = column;
 //count[j] = 1;
 } else {// not the column
 start[j] = 0;
 count[j] = v->dims[j];
 }
 }
 //start[v->ndim]=0;
 //count[v->n dim]=0;
 }
 */

enum ADIOS_PREDICATE_MODE getOp(const char* opStr) {
	if (strcmp(opStr, ">=") == 0) {
		return ADIOS_GTEQ;
	} else if (strcmp(opStr, "<=") == 0) {
		return ADIOS_LTEQ;
	} else if (strcmp(opStr, "<") == 0) {
		return ADIOS_LT;
	} else if (strcmp(opStr, ">") == 0) {
		return ADIOS_GT;
	} else if (strcmp(opStr, "=") == 0) {
		return ADIOS_EQ;
	} else { // if (strcmp(opStr, "!=") == 0) {
		return ADIOS_NE;
	}
}

ADIOS_QUERY* getQuery(const char* condition, ADIOS_FILE* f) {

	char* varStr;// = malloc(sizeof(char) * strlen(condition));

	//char opStr[5];
	char* opStr;// = malloc(sizeof(char)*5);
	//char value[strlen(condition)];
	char* valueStr; // = malloc(sizeof(char)*strlen(condition));

	//parseVar(condition, &varStr, &opStr, &valueStr); // these values are ok in parseVar but just valueStr became null after this call!!!??

	//char str[strlen(condition)];
	//strncpy(str, condition, strlen(condition));
	//str[strlen(condition)] = 0;
	char* str = strdup(condition);

	char * pch;
	pch = strtok(str, " ");

	varStr = strdup(pch);

	pch = strtok(NULL, " ");
	opStr = strdup(pch);

	pch = strtok(NULL, " ");
	valueStr = strdup(pch);

	char* varName; //[strlen(varStr)];
	char* dimDef; //[strlen(varStr)];
	getVarName(varStr, &varName, &dimDef);

	ADIOS_VARINFO* v = adios_inq_var(f, varName);
	//ADIOS_VARINFO* v = getAdiosVariable(f, varName);

	if (v == NULL) {
		free(valueStr);
		free(opStr);
		free(varStr);
		free(varName);
		free(str);
		if (dimDef != NULL) {
			free(dimDef);
		}
		return NULL;
	}

	uint64_t* start = malloc(sizeof(uint64_t) * v->ndim);
	uint64_t* count = malloc(sizeof(uint64_t) * v->ndim);

	createBox(v, dimDef, start, count);

	ADIOS_SELECTION* sel = adios_selection_boundingbox(v->ndim, start, count);

	adios_free_varinfo(v);

	ADIOS_QUERY* q =
			adios_query_create(f, varName, sel, getOp(opStr), valueStr);

	free(valueStr);
	free(opStr);
	free(varStr);
	free(varName);
	free(str);
	//free(start); free(count); // if deleted, then adios_sel values would be affected
	if (dimDef != NULL) {
		free(dimDef);
	}

	printf("::\t query created for: %s\n", condition);
	return q;
}

void queryDetail(ADIOS_QUERY* q, int timeStep) {
	int64_t estimated = adios_query_estimate(q);
	printf("::\t query estimated = %llu \n", estimated);

	int64_t numHits = adios_query_fastbit_evaluate_method(q, timeStep, 100000);
	printf("::\t query evaluated = %llu \n", numHits);

	uint64_t batchSize = 50;

	while (q->_maxResultDesired - q->_lastRead > 0) {
		//printf("::\t max=%llu _lastRead=%llu\n", q->_maxResultDesired, q->_lastRead);
		ADIOS_SELECTION* t;
		ADIOS_SELECTION* bound;
		adios_query_get_selection(q, batchSize, bound, &t);

		adios_selection_delete(t);
		//printf("::\t      max=%llu _lastRead=%llu\n", q->_maxResultDesired, q->_lastRead);
	}
}

/*
 void queryDetailOld(ADIOS_QUERY* q) {
 uint64_t numHits = adios_query_fastbit_evaluate_method(q, 0, 100000);

 printf(":: ==> Num of hits found for [%s] = %llu\n", q->condition, numHits);

 int64_t coordinates[numHits];
 fastbit_selection_get_coordinates(q->_queryInternal, coordinates, numHits, 0);

 int i=0;
 printf("\n:: coordinate: [");
 for (i=0; i<numHits; i++) {
 if (i< 100) {
 printf("%lld ", coordinates[i]);
 } else {
 printf(" ... ");
 break;
 }
 }
 printf("]\n\n");
 }
 */

/*
 int main (int argc, char ** argv)
 {
 if (argc <= 2) {
 usage(argv[0]);
 return 1;
 }

 if (argc == 4) {
 usage (argv[0]);
 return 1;
 }

 adios_query_init(ADIOS_QUERY_TOOL_FASTBIT);

 ADIOS_FILE * f;
 MPI_Comm    comm_dummy = 0;  // MPI_Comm is defined through adios_read.h

 f = adios_read_open_file (argv[1], ADIOS_READ_METHOD_BP, comm_dummy);
 if (f == NULL) {
 printf ("::%s\n", adios_errmsg());
 return -1;
 }

 const char* condition = argv[2];
 //printf(":: condition = %s \n", condition);

 ADIOS_QUERY* q = getQuery(condition, f);

 if (q != NULL) {
 if (argc > 3) {
 const char* opStr = argv[3];
 enum ADIOS_CLAUSE_OP_MODE op = ADIOS_QUERY_OP_AND;
 if ((strcmp(opStr, "OR") == 0) || (strcmp(opStr, "or") == 0)) {
 op = ADIOS_QUERY_OP_OR;
 }

 const char* condition2 = argv[4];
 ADIOS_QUERY* q2 = getQuery(condition2, f);

 if (q2 != NULL) {
 ADIOS_QUERY* combined = adios_query_combine(q, op, q2); //adios_query_combine_fastbit(q, op, q2);

 queryDetail(combined, 0);

 //adios_query_fastbit_free_method(q2);
 // looks like can not free q and q2 and combined. causes problem in fastbit
 fastbit_selection_free(combined->_queryInternal);

 adios_query_free(combined);
 adios_query_free(q);
 adios_query_free(q2);

 adios_query_clean();
 adios_read_close(f);

 return;
 }
 } else {
 queryDetail(q, 0);
 }

 fastbit_selection_free(q->_queryInternal);
 //adios_query_fastbit_free_method(q);
 adios_query_free(q);
 }

 adios_query_clean();
 adios_read_close(f);

 }
 */
