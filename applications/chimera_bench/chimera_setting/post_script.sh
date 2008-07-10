#!/bin/sh

RESULT_DIR=wwww/chimera_results/$1/$2

cd $RESULT_DIR 

# extract timing info
echo "Open time goes to " $RESULT_DIR/open
grep "# Open        :" ./*/log|awk '{print NR "\t" $4 "\t" $5 "\t" $6 "\t" $7}' |tee $RESULT_DIR/open
echo "Write time goes to " $RESULT_DIR/write
grep "# Write       :" ./*/log|awk '{print NR "\t" $4 "\t" $5 "\t" $6 "\t" $7}' |tee $RESULT_DIR/write
echo "Close time goes to " $RESULT_DIR/close
grep "# Close       :" ./*/log|awk '{print NR "\t" $4 "\t" $5 "\t" $6 "\t" $7}' |tee $RESULT_DIR/close
echo "Total IO time goes to " $RESULT_DIR/total
grep "# Total       :" ./*/log|awk '{print NR "\t" $4 "\t" $5 "\t" $6 "\t" $7}' |tee $RESULT_DIR/total


