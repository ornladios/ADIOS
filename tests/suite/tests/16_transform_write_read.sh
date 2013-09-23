
PROCS=16
TEST_NAME=transforms_read_write
PROGRAM_NAME=adios_${TEST_NAME}
OUTPUT_FILENAME=${PROGRAM_NAME}.bp

cp ${SRCDIR}/programs/${PROGRAM_NAME} .
cp ${SRCDIR}/programs/adios_transforms.xml .

echo "Run C ${PROGRAM_NAME}"
$MPIRUN $NP_MPIRUN $PROCS ./${PROGRAM_NAME}
 
EX=$?
if [ ! -f ${OUTPUT_FILENAME} ]; then
    echo "ERROR: C version of ${PROGRAM_NAME} failed. No BP file is created. Exit code=$EX"
    exit 1
fi

echo "Check output with bpls"
$TRUNKDIR/utils/bpls/bpls -la ${OUTPUT_FILENAME} | grep -v -e endianness -e 'file size' > ${TEST_NAME}_bpls.txt
diff -q ${TEST_NAME}_bpls.txt $SRCDIR/reference/${TEST_NAME}_bpls.txt

if [ $? != 0 ]; then
    echo "ERROR: C version of ${PROGRAM_NAME} produced a file different from the reference."
    echo "Compare \"bpls -la $PWD/${OUTPUT_FILENAME} | grep -v -e endianness -e 'file size'\" to reference $SRCDIR/reference/${TEST_NAME}_bpls.txt"
    exit 1
fi

