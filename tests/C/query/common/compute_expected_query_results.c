/*
 * compute_expected_query_results.c
 *
 *  Created on: Sep 30, 2014
 *      Author: David A. Boyuka II
 */

#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <assert.h>

#include <adios_read.h>
#include <adios_read_ext.h>
#include <adios_query.h>
#include <adios_logger.h>
#include "adios_query_xml_parse.h"

#define max(a,b) ((a)>(b)?(a):(b))

// Returned ADIOS selection must be freed after use, and will only be valid as long as the supplied varinfo struct is valid
static ADIOS_SELECTION * convertWBToBB(ADIOS_SELECTION *sel, int timestep, ADIOS_FILE *fp, ADIOS_VARINFO *varinfo) {
	assert(sel->type == ADIOS_SELECTION_WRITEBLOCK);

	const int wbindex = sel->u.block.index;
	const int abs_wbindex = adios_get_absolute_writeblock_index(varinfo, wbindex, timestep);

	if (!varinfo->blockinfo) {
		const data_view_t old_view = adios_read_set_data_view(fp, LOGICAL_DATA_VIEW);
		adios_read_bp_inq_var_blockinfo(fp, varinfo);
		adios_read_set_data_view(fp, old_view);
	}
	ADIOS_VARBLOCK *vb = &varinfo->blockinfo[abs_wbindex];

	return adios_selection_boundingbox(varinfo->ndim, vb->start, vb->count);
}

static uint64_t computeSelectionSizeInElements(ADIOS_SELECTION *sel) {
	switch (sel->type) {
	case ADIOS_SELECTION_BOUNDINGBOX: {
		int i;
		ADIOS_SELECTION_BOUNDINGBOX_STRUCT *bb = &sel->u.bb;
		uint64_t size = 1;
		for (i = 0; i < bb->ndim; ++i)
			size *= bb->count[i];
		return size;
	}
	default:
		fprintf(stderr, "Unsupported selection type %d at %s:%s", sel->type, __FILE__, __LINE__);
		exit(1);
		return 0;
	}
}

enum REDUCED_DATATYPE { RD_UNKNOWN = -1, RD_SIGNED_INT, RD_UNSIGNED_INT, RD_DOUBLE, RD_LONG_DOUBLE, RD_STRING, RD_DOUBLE_COMPLEX };
static enum REDUCED_DATATYPE getReducedDatatype(enum ADIOS_DATATYPES datatype) {
	switch (datatype) {
	case adios_byte:
	case adios_short:
	case adios_integer:
	case adios_long:
		return RD_SIGNED_INT;
	case adios_unsigned_byte:
	case adios_unsigned_short:
	case adios_unsigned_integer:
	case adios_unsigned_long:
		return RD_UNSIGNED_INT;
	case adios_real:
	case adios_double:
		return RD_DOUBLE;
	case adios_long_double:
		return RD_LONG_DOUBLE;
	case adios_string:
		return RD_STRING;
	case adios_complex:
	case adios_double_complex:
		return RD_DOUBLE_COMPLEX;
	}
	return RD_UNKNOWN;
}

#define ALLOCATE_REDUCED_DATATYPE(type) (type*)malloc(sizeof(type))
static void * allocateReducedDatatype(enum REDUCED_DATATYPE reduced_datatype, int maxstrlen) {
	switch (reduced_datatype) {
	case RD_SIGNED_INT:		return ALLOCATE_REDUCED_DATATYPE(int64_t);
	case RD_UNSIGNED_INT:	return ALLOCATE_REDUCED_DATATYPE(uint64_t);
	case RD_DOUBLE:			return ALLOCATE_REDUCED_DATATYPE(double);
	case RD_LONG_DOUBLE:	return ALLOCATE_REDUCED_DATATYPE(long double);
	case RD_STRING:			return (char*)malloc(maxstrlen + 1);
	case RD_DOUBLE_COMPLEX:
		fprintf(stderr, "Cannot handle complex or double complex datatypes (at %s:%s)\n", __FILE__, __LINE__);
		exit(1);
		return NULL;
	}
}

#define RETURN_REDUCED_DATATYPE(type, val) { type *rd = (type *)malloc(sizeof(type)); *rd = (val); return rd; }
static void * parseStringAsReducedDatatype(const char *str, enum REDUCED_DATATYPE reduced_datatype) {
	switch (reduced_datatype) {
	case RD_SIGNED_INT:		RETURN_REDUCED_DATATYPE(int64_t, strtoll(str, NULL, 0)); break;
	case RD_UNSIGNED_INT:	RETURN_REDUCED_DATATYPE(uint64_t, strtoull(str, NULL, 0)); break;
	case RD_DOUBLE:			RETURN_REDUCED_DATATYPE(double, strtod(str, NULL)); break;
	case RD_LONG_DOUBLE:	RETURN_REDUCED_DATATYPE(long double, strtold(str, NULL)); break;
	case RD_STRING:			return strdup(str);
	case RD_DOUBLE_COMPLEX:
		fprintf(stderr, "Cannot handle complex or double complex datatypes (at %s:%s)\n", __FILE__, __LINE__);
		exit(1);
		return NULL;
	}
}

static void castToReducedDatatype(const void *value, enum ADIOS_DATATYPES datatype, void *outValue) {
	switch (datatype) {
	case adios_byte:				*(int64_t*)outValue = *(const int8_t*)value; break;
	case adios_short:				*(int64_t*)outValue = *(const int16_t*)value; break;
	case adios_integer:				*(int64_t*)outValue = *(const int32_t*)value; break;
	case adios_long:				*(int64_t*)outValue = *(const int64_t*)value; break;
	case adios_unsigned_byte:		*(uint64_t*)outValue = *(const uint8_t*)value; break;
	case adios_unsigned_short:		*(uint64_t*)outValue = *(const uint16_t*)value; break;
	case adios_unsigned_integer:	*(uint64_t*)outValue = *(const uint32_t*)value; break;
	case adios_unsigned_long:		*(uint64_t*)outValue = *(const uint64_t*)value; break;
	case adios_real:				*(double*)outValue = *(const float*)value; break;
	case adios_double:				*(double*)outValue = *(const double*)value; break;
	case adios_long_double:			*(long double*)outValue = *(const long double*)value; break;
	case adios_string:				strcpy((char*)outValue, (const char*)value); break;
	case adios_complex:
	case adios_double_complex:
	default:
		fprintf(stderr, "Unsupported or invalid reduced datatype %d at %s:%s\n", datatype, __FILE__, __LINE__);
		exit(1);
	}
}

#define RETURN_COMPARE_PTRS_CAST_TO_TYPE(ptr1, ptr2, type) { const type __v1 = *(const type*)(ptr1), __v2 = *(const type*)(ptr2); return __v1 < __v2 ? -1 : __v1 > __v2 ? 1 : 0; }
static int compareReducedDatatypeValues(const void *v1, const void *v2, enum REDUCED_DATATYPE datatype) {
	switch (datatype) {
	case RD_SIGNED_INT:		RETURN_COMPARE_PTRS_CAST_TO_TYPE(v1, v2, int64_t); break;
	case RD_UNSIGNED_INT:	RETURN_COMPARE_PTRS_CAST_TO_TYPE(v1, v2, int64_t); break;
	case RD_DOUBLE:			RETURN_COMPARE_PTRS_CAST_TO_TYPE(v1, v2, double); break;
	case RD_LONG_DOUBLE:	RETURN_COMPARE_PTRS_CAST_TO_TYPE(v1, v2, long double); break;
	case RD_STRING:			return strcmp((const char*)v1, (const char*)v2);
	case RD_DOUBLE_COMPLEX:
	default:
		fprintf(stderr, "Unsupported or invalid reduced datatype %d at %s:%s\n", datatype, __FILE__, __LINE__);
		exit(1);
		return 0;
	}
}

// NOTE: both bound and value must be (the same) reduced datatype
static int compareConstraintBoundValue(const void *bound, const void *value, enum REDUCED_DATATYPE reducedDatatype, enum ADIOS_PREDICATE_MODE comparison) {
	int compare = compareReducedDatatypeValues(value, bound, reducedDatatype);

	switch (comparison) {
	case ADIOS_LT: return compare < 0;
	case ADIOS_LTEQ: return compare <= 0;
	case ADIOS_GT: return compare > 0;
	case ADIOS_GTEQ: return compare >= 0;
	case ADIOS_EQ: return compare == 0;
	case ADIOS_NE: return compare != 0;
	}
}

// Returns which points in the given buffer of data (buffer) for the given selection
// (insel) match the constraint in the given query (query), returning a list of points
// that are relative to (insel)
static ADIOS_SELECTION * scanBufferForMatchingPoints(const char *buffer, enum ADIOS_DATATYPES datatype, ADIOS_SELECTION *insel, ADIOS_QUERY *query) {
	assert(insel->type == ADIOS_SELECTION_BOUNDINGBOX); // For now, only support bounding boxes (and writeblocks, since they are converted to bounding boxes earlier)

	const int ndim = query->_var->ndim;
	const ADIOS_SELECTION_BOUNDINGBOX_STRUCT *bb = &insel->u.bb;

	const int datatypeSize = adios_type_size(datatype, NULL);
	const enum REDUCED_DATATYPE reducedDatatype = getReducedDatatype(datatype);
	const enum ADIOS_PREDICATE_MODE comparison = query->_op;

	const void *boundValue = parseStringAsReducedDatatype(query->_value, reducedDatatype);
	void *pointValue = allocateReducedDatatype(reducedDatatype, 0);

	uint64_t elemsRemaining = computeSelectionSizeInElements(insel);

	uint64_t npoints = 0;
	uint64_t pointsCapacity = 1;
	uint64_t *points = (uint64_t *)calloc(pointsCapacity, ndim * sizeof(uint64_t)); // First coordinate is at 0,0,0,...,0, since the results should be relative to the selection box
	assert(points);
	uint64_t *nextPoint = points;

	int i;
	while (elemsRemaining-- > 0) {
		// Cast the point's value to a reduced datatype
		castToReducedDatatype(buffer, datatype, pointValue);

		// Compare the point value to the bound value
		if (compareConstraintBoundValue(boundValue, pointValue, reducedDatatype, comparison)) {
			// Commit the current point as a real point
			++npoints;

			// Expand the point array if need be
			if (npoints == pointsCapacity) {
				const uint64_t next_point_offset = nextPoint - points;
				pointsCapacity *= 2;
				points = (uint64_t*)realloc(points, pointsCapacity * ndim * sizeof(uint64_t));
				assert(points);
				nextPoint = points + next_point_offset;
			}

			// Move to the next point
			const uint64_t *curPoint = nextPoint;
			nextPoint += ndim;

			// Copy the current point's coordinates to the next point's coordinates,
			// so they will continue to be incremented from here
			memcpy(nextPoint, curPoint, ndim * sizeof(uint64_t));
		}

		// Increment the next point's coordinates
		for (i = ndim - 1; i >= 0; --i) {
			++nextPoint[i];
			if (nextPoint[i] == bb->count[i]) {
				nextPoint[i] = 0;
			} else {
				break;
			}
		}

		// Advance in the input buffer
		buffer += datatypeSize;
	}

	free((void*)boundValue);
	free(pointValue);
	return adios_selection_points(ndim, npoints, points);
}

static int pointLexCompareNumDims;
static void setPointLexCompareNumDims(int ndim) {
	pointLexCompareNumDims = ndim;
}
static int pointLexCompare(const void *left, const void *right) {
	const int ndim = pointLexCompareNumDims;
	const uint64_t *leftpt = (const uint64_t *)left;
	const uint64_t *rightpt = (const uint64_t *)right;
	int i;
	for (i = 0; i < ndim; ++i) {
		if (*leftpt < *rightpt)
			return -1;
		else if (*leftpt > *rightpt)
			return 1;
		++leftpt;
		++rightpt;
	}	return 0;

}

static void sortPointsLexOrder(ADIOS_SELECTION *pointsel) {
	assert(pointsel->type == ADIOS_SELECTION_POINTS);
	setPointLexCompareNumDims(pointsel->u.points.npoints); // Sets the ndim for pointLexCompare via global variable, since we can't pass any parameters
	qsort(pointsel->u.points.points, pointsel->u.points.npoints, pointsel->u.points.ndim * sizeof(uint64_t), pointLexCompare);
}

// Returns a point selection with points in lexicographical order
static ADIOS_SELECTION * evaluateConstraint(ADIOS_QUERY *query, int timestep) {
	assert(!query->_left && !query->_right);
	assert(query->_var && query->_f && query->_sel);

	ADIOS_SELECTION *insel = query->_sel;
	int free_insel = 0;

	if (insel->type == ADIOS_SELECTION_WRITEBLOCK) {
		insel = convertWBToBB(insel, timestep, query->_f, query->_var);
		free_insel = 1;
	}

	const uint64_t buffersize = computeSelectionSizeInElements(insel) * adios_type_size(query->_var->type, NULL);
	char *buffer = (char *)malloc(buffersize);
	assert(buffer);

	adios_schedule_read_byid(query->_f, insel, query->_var->varid, timestep, 1, buffer);
	adios_perform_reads(query->_f, 1);

	ADIOS_SELECTION *results = scanBufferForMatchingPoints(buffer, query->_var->type, insel, query);
	sortPointsLexOrder(results); // Sort the matching points in lexicographical order

	if (free_insel)
		adios_selection_delete(insel);

	return results;
}

// Takes two point selections with points in lexicographical order, and
// returns a combined point selection with points in lexicographical order
static ADIOS_SELECTION * computePointListCombination(enum ADIOS_CLAUSE_OP_MODE op, ADIOS_SELECTION *leftsel, ADIOS_SELECTION *rightsel) {
	assert(leftsel->u.points.ndim == rightsel->u.points.ndim);
	const int ndim = leftsel->u.points.ndim;

	// Allocate an array for the combined points
	const uint64_t maxNewPoints =
			(op == ADIOS_QUERY_OP_AND) ?
				max(leftsel->u.points.npoints, rightsel->u.points.npoints) :
				leftsel->u.points.npoints + rightsel->u.points.npoints;
	uint64_t *newPoints = (uint64_t *)malloc(maxNewPoints * ndim * sizeof(uint64_t));
	assert(newPoints);

	// Set up iterator pointers
	const uint64_t *leftHeadPoint = leftsel->u.points.points;
	const uint64_t *rightHeadPoint = rightsel->u.points.points;
	const uint64_t *leftPointsEnd = leftHeadPoint + ndim * leftsel->u.points.npoints;
	const uint64_t *rightPointsEnd = rightHeadPoint + ndim * rightsel->u.points.npoints;

	uint64_t newNPoints = 0;
	uint64_t *curPoint = newPoints;

	// Perform a list conjunction/disjunction on the point lists
	setPointLexCompareNumDims(ndim);
	while (leftHeadPoint != leftPointsEnd && rightHeadPoint != rightPointsEnd) {
		// Compare the head points of both point lists
		const int compare = pointLexCompare(leftHeadPoint, rightHeadPoint);

		// If we are taking the OR of the lists, or we are taking the AND and
		// the head points are equal, copy the lesser of the head points to the output
		// (copying the lesser point works for both AND and OR cases)
		if (op == ADIOS_QUERY_OP_OR || (op == ADIOS_QUERY_OP_AND && compare == 0)) {
			memcpy(curPoint, (compare < 0) ? leftHeadPoint : rightHeadPoint, ndim * sizeof(uint64_t));
			++newNPoints;
			curPoint += ndim;
		}

		// Increment the point list with the lesser head point
		// (or both lists if the points are equal)
		if (compare <= 0)
			leftHeadPoint += ndim;
		if (compare >= 0)
			rightHeadPoint += ndim;
	}

	// If we are in OR mode, copy any remaining points in the unexhausted
	// list (if one of the two point lists is so)
	if (op == ADIOS_QUERY_OP_OR) {
		if (leftHeadPoint != leftPointsEnd) {
			const uint64_t coordsRemaining = (leftPointsEnd - leftHeadPoint);
			memcpy(curPoint, leftHeadPoint, coordsRemaining * sizeof(uint64_t));
			newNPoints += coordsRemaining / ndim;
		} else if (rightHeadPoint != rightPointsEnd) {
			const uint64_t coordsRemaining = (rightPointsEnd - rightHeadPoint);
			memcpy(curPoint, rightHeadPoint, coordsRemaining * sizeof(uint64_t));
			newNPoints += coordsRemaining / ndim;
		}
	}

	// Free the left and right point lists, since they are no longer needed
	free(leftsel->u.points.points);
	free(rightsel->u.points.points);
	adios_selection_delete(leftsel);
	adios_selection_delete(rightsel);

	// Return the combined point list
	return adios_selection_points(ndim, newNPoints, newPoints);
}

static ADIOS_SELECTION * evaluateQueryTree(ADIOS_QUERY *query, int timestep) {
	if (!query->_left && !query->_right) {
		return evaluateConstraint(query, timestep);
	} else if (query->_left && query->_right) {
		const enum ADIOS_CLAUSE_OP_MODE op = query->_leftToRightOp;
		ADIOS_SELECTION *leftsel = evaluateQueryTree(query->_left, timestep);
		ADIOS_SELECTION *rightsel = evaluateQueryTree(query->_right, timestep);

		ADIOS_SELECTION *combinedsel = computePointListCombination(op, leftsel, rightsel);
		return combinedsel;
	} else if (query->_left) {
		return evaluateQueryTree(query->_left, timestep);
	} else if (query->_right) {
		return evaluateQueryTree(query->_right, timestep);
	}
}

static ADIOS_SELECTION * derelativizePoints(ADIOS_SELECTION *inputPointsSel, ADIOS_SELECTION *outputSelection) {
	ADIOS_SELECTION_POINTS_STRUCT *inputPoints = &inputPointsSel->u.points;
	const int ndim = inputPoints->ndim;

	assert(outputSelection->type == ADIOS_SELECTION_BOUNDINGBOX);
	const uint64_t *outputOffset = outputSelection->u.bb.start;

	uint64_t i, j;
	uint64_t *curPoint = inputPoints->points;
	for (i = 0; i < inputPoints->npoints; ++i)
		for (j = 0; j < ndim; ++j)
			*curPoint++ += outputOffset[j];

	return inputPointsSel;
}

static ADIOS_SELECTION * computeExpectedQueryResults(ADIOS_QUERY *query, int timestep, ADIOS_SELECTION *outputSelection) {
	ADIOS_SELECTION *resultPointsSel = evaluateQueryTree(query, timestep);

	int freeOutputSelection = 0;
	if (outputSelection->type == ADIOS_SELECTION_WRITEBLOCK) {
		fprintf(stderr, "Writeblock output selections are currently not supported in compute_expected_query_results (at %s:%s)\n", __FILE__, __LINE__);
		abort();
		return NULL;
		//outputSelection = convertWBToBB(outputSelection, timestep, query->_fp);
	}

	derelativizePoints(resultPointsSel, outputSelection);

	if (freeOutputSelection)
		adios_selection_delete(outputSelection);

	return resultPointsSel;
}

static void printPointSelection(int timestep, ADIOS_SELECTION *sel) {
	assert(sel->type == ADIOS_SELECTION_POINTS);

	const ADIOS_SELECTION_POINTS_STRUCT *pstruct = &sel->u.points;
	const int ndim = pstruct->ndim;
	const uint64_t npoints = pstruct->npoints;
	const uint64_t *points = pstruct->points;

	uint64_t i;
	int j;
	for (i = 0; i < npoints; ++i) {
		printf("%d", timestep);
		for (j = 0; j < ndim; ++j) {
			printf(" %llu", *points++);
		}
		printf("\n");
	}
}

static void usage(const char *cmd) {
	fprintf(stderr, "Usage: %s <query XML file> <input BP file>", cmd);
}

#define SHIFT_N(n) { argc -= (n); argv += (n); }
#define SHIFT SHIFT_N(1)
int main(int argc, char **argv) {
	const char *cmd = *argv; SHIFT;
	if (argc != 2) {
		usage(cmd);
		exit(1);
	}

	const char *inputxml_filename = *argv; SHIFT;
	const char *bp_filename = *argv; SHIFT;

	const MPI_Comm comm = MPI_COMM_WORLD;

	MPI_Init(&argc, &argv);
	adios_read_init_method(ADIOS_READ_METHOD_BP, comm, "");
	adios_query_init(ADIOS_QUERY_TOOL_ALACRITY);

	ADIOS_FILE *bp_file = adios_read_open_file(bp_filename, ADIOS_READ_METHOD_BP, comm);
	if (bp_file == NULL) {
		log_error("Error: could not read input dataset %s\n", bp_filename);
		exit(1);
	}

	ADIOS_QUERY_TEST_INFO *testinfo = parseXml(inputxml_filename, bp_file);
	if (testinfo == NULL) {
		log_error("Error: could not read query XML file %s\n", inputxml_filename);
		exit(1);
	}

	int timestep;
	for (timestep = testinfo->fromStep; timestep < testinfo->fromStep + testinfo->numSteps; ++timestep) {
		ADIOS_SELECTION *result = computeExpectedQueryResults(testinfo->query, timestep, testinfo->outputSelection);
		printPointSelection(timestep, result);

		free(result->u.points.points);
		adios_selection_delete(result);
	}

	adios_selection_delete(testinfo->outputSelection); // TODO: leaks start[] and count[] if it's a BB
	adios_query_free(testinfo->query);
	free(testinfo);
	adios_read_close(bp_file);

	adios_query_clean();
	adios_read_finalize_method(ADIOS_READ_METHOD_BP);
	MPI_Finalize();
}

