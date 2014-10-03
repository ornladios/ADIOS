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

# Check for the executability of all executables and the "directory"-ness
# of all directories that we need
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
  
  # Rename the ADIOS XML used for writing so it's more clear what it is for
  mv "$DSID.xml" "$DSID.create.xml"
  
  # Ensure the dataset was actually produced
  INDEXED_DS="$DSID.bp"
  [ -f "$INDEXED_DS" ] ||
    die "ERROR: $DATASET_BUILDER_EXE did not produce expected output BP file \"$INDEXED_DS\""

  # Iterate over all pre-defined queries
  for QUERY_XML in "$QUERY_XML_DIR/$DSID/"/*; do
    QUERY_XML_BASENAME="${QUERY_XML##*/}"     # Strip the path
    QUERY_NAME="${QUERY_XML_BASENAME%%.xml}"  # Strip the .xml suffix
    
    # Copy the query XML into our working directory for easy access
    QUERY_XML_LOCAL="${QUERY_XML_BASENAME}"
    cp "$QUERY_XML" "$QUERY_XML_LOCAL" 
    
    # Compute the expected results
    EXPECTED_POINTS_FILE="$QUERY_NAME.expected-points.txt"
    $MPIRUN_SERIAL "$QUERY_SEQSCAN_EXE" "$INDEXED_DS" "$QUERY_XML" > "$EXPECTED_POINTS_FILE" ||
      die "ERROR: $QUERY_SEQSCAN_EXE failed with exit code $?"

    # NOTE: the sequential scan program produces a point list that is guaranteed to be sorted in C array order, so no need to sort it here

    # Run the query for each query engine implementation and compare to the expected results
    for QUERY_ENGINE in $ALL_QUERY_ENGINES; do
      OUTPUT_POINTS_FILE="$QUERY_NAME.$QUERY_ENGINE-points.txt"
    
      # Run the query through ADIOS Query to get actual results
      echo \
      $MPIRUN_SERIAL "$QUERY_EXE" "$INDEXED_DS" "$QUERY_XML_LOCAL" "$QUERY_ENGINE" \> "$OUTPUT_POINTS_FILE"
      $MPIRUN_SERIAL "$QUERY_EXE" "$INDEXED_DS" "$QUERY_XML_LOCAL" "$QUERY_ENGINE"  > "$OUTPUT_POINTS_FILE" ||
        die "ERROR: $QUERY_EXE failed with exit code $RET"

      # Sort the output points in C array order, since the query engine makes no guarantee as to the ordering of the results
      # Sort file in place (-o FILE) with numerical sort order (-n) on each of the first 9 fields (-k1,1 ...)
      # Assumes the output will have at most 8 dimensions (+ 1 timestep column == 9), add more if needed (or a generalized column counter)
      sort -n -k1,1 -k2,2 -k3,3 -k4,4 -k5,5 -k6,6 -k7,7 -k8,8 -k9,9 \
        "$OUTPUT_POINTS_FILE" -o "$OUTPUT_POINTS_FILE"

      # Compare the actual and expected results via diff (the matching points are sorted lexicographically
      if ! diff -q "$EXPECTED_POINTS_FILE" "$OUTPUT_POINTS_FILE"; then
        echo "ERROR: ADIOS Query does not return the expected points matching query $QUERY_NAME on dataset $DSID using query engine $QUERY_ENGINE"
        echo "Compare \"$PWD/$EXPECTED_POINTS_FILE\" (expected points) vs. \"$PWD/$OUTPUT_POINTS_FILE\" (returned points)"
        echo "The BP file queried is at \"$INDEXED_DS\" and the query is specified by \"$QUERY_XML\""
        exit 1  
      fi
    done
  done
done
