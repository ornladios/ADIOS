#!/bin/bash
#
# Test if the ADIOS query framework can successfully perform a range of *serial* queries on simple, pre-defined datasets
# Parallel query tests are tested in *TODO*
#
# Uses:
# - tests/C/query/common/build_indexed_dataset (executable, produces the pre-defined datasets)
# - tests/C/query/common/xml-testcases/*/* (custom XML files describing interesting queries over the predefined datasets)
#   - e.g. tests/C/query/common/xml-testcases/DS1/simple-query.xml
# - tests/C/query/common/compute_expected_query_results (executable, runs queries via sequential scan over a BP file)
#
# Environment variables set by caller:
# MPIRUN        Run command
# NP_MPIRUN     Run commands option to set number of processes
# MAXPROCS      Max number of processes allowed
# HAVE_FORTRAN  yes or no
# SRCDIR        Test source dir (.. of this script)
# TRUNKDIR      ADIOS trunk dir

function die() {
  EC=${2-1}
  echo "$1" >&2
  exit $EC
}

QUERY_TEST_DIR="$TRUNKDIR/tests/C/query"
QUERY_COMMON_DIR="$QUERY_TEST_DIR/common"

# Some external tools to use
DATASET_BUILDER_EXE="$QUERY_COMMON_DIR/build_indexed_dataset"
QUERY_SEQSCAN_EXE="$QUERY_COMMON_DIR/compute_expected_query_results"
QUERY_EXE="$QUERY_COMMON_DIR/adios_query_test"
QUERY_XML_DIR="$QUERY_TEST_DIR/query-xmls/"

[ -x "$DATASET_BUILDER_EXE" ] || die "ERROR: $DATASET_BUILDER_EXE is not executable"
[ -x "$QUERY_SEQSCAN_EXE" ] || die "ERROR: $QUERY_SEQSCAN_EXE is not executable"
[ -x "$QUERY_EXE" ] || die "ERROR: $QUERY_EXE is not executable"
[ -d "$QUERY_XML_DIR" ] || die "ERROR: $QUERY_XML_DIR is not a directory"

# mpirun for serial command
MPIRUN_SERIAL="$MPIRUN $NP_MPIRUN 1 $EXEOPT"

# All pre-defined dataset IDs (which can be extracted from build_indexed_dataset)
ALL_DATASET_IDS="DS1 DS2 DS3"

# All query engine implementations to test
ALL_QUERY_ENGINES="alacrity" # TODO: Expand to FastBit once we understand it properly

for DSID in $ALL_DATASET_IDS; do
  # Build the pre-defined dataset with index (TODO: generalize for ALACRITY)
  $DATASET_BUILDER_EXE $DSID $DSID alacrity ||
    die "ERROR: $DATASET_BUILDER_EXE failed with exit code $?"
  
  INDEXED_DS="$DSID.bp"
  [ -f "$INDEXED_DS" ] ||
    die "ERROR: $DATASET_BUILDER_EXE did not produce expected output BP file \"$INDEXED_DS\""

  # Iterate over all interesting queries:
  for QUERY_XML in "$QUERY_XML_DIR/$DSID/"/*; do
    QUERY_NAME="$QUERY_XML"
    QUERY_NAME="${QUERY_NAME##*/}"
    QUERY_NAME="${QUERY_NAME%%.xml}"
    
    EXPECTED_POINTS_FILE="$QUERY_NAME.expected-points.txt"
    
    # Compute the expected results
    $MPIRUN_SERIAL "$QUERY_SEQSCAN_EXE" "$QUERY_XML" "$INDEXED_DS" > "$EXPECTED_POINTS_FILE" ||
      die "ERROR: $QUERY_SEQSCAN_EXE failed with exit code $?"

    # Run the query for each query engine implementation and compare to the expected results
    for QUERY_ENGINE in $ALL_QUERY_ENGINES; do
      OUTPUT_POINTS_FILE="$QUERY_NAME.$QUERY_ENGINE.output-points.txt"
    
      # Run the query through ADIOS Query to get actual results
      $MPIRUN_SERIAL "$QUERY_EXE" "$INDEXED_DS" "$QUERY_XML" "$QUERY_ENGINE" > "$OUTPUT_POINTS_FILE" ||
        die "ERROR: $QUERY_EXE failed with exit code $?"

      if ! diff -q "$EXPECTED_POINTS_FILE" "$OUTPUT_POINTS_FILE"; then
        echo "ERROR: ADIOS Query does not return the expected points matching query $QUERY_NAME on dataset $DSID using query engine $QUERY_ENGINE"
        echo "Compare \"$PWD/$EXPECTED_POINTS_FILE\" (expected points) vs. \"$PWD/$OUTPUT_POINTS_FILE\" (returned points)"
        echo "The BP file queried is at \"$INDEXED_DS\" and the query is specified by \"$QUERY_XML\""
        exit 1  
      fi
    done
  done
done
