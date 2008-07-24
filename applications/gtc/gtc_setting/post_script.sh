#!/bin/sh

RESULT_DIR=wwww/gtc_results/$1/$2/$3
TEMP_DIR=wwww/temp/gtc/$1/$2/$3

cd $RESULT_DIR 
mkdir -p $TEMP_DIR

# extract timing info
for i in diagnosis.0 diagnosis.1 diagnosis.2 output3d.0 output3d.1 particles restart snapshot
do
echo "Open time goes to " $TEMP_DIR/$i.open
grep "# Open        :" ./*/$i.prof|awk '{print NR "\t" $4 "\t" $5 "\t" $6 "\t" $7}' |tee $TEMP_DIR/$i.open
echo "Write time goes to " $TEMP_DIR/$i.write
grep "# Write       :" ./*/$i.prof|awk '{print NR "\t" $4 "\t" $5 "\t" $6 "\t" $7}' |tee $TEMP_DIR/$i.write
echo "Close time goes to " $TEMP_DIR/$i.close
grep "# Close       :" ./*/$i.prof|awk '{print NR "\t" $4 "\t" $5 "\t" $6 "\t" $7}' |tee $TEMP_DIR/$i.close
echo "Total IO time goes to " $TEMP_DIR/$i.total
grep "# Total       :" ./*/$i.prof|awk '{print NR "\t" $4 "\t" $5 "\t" $6 "\t" $7}' |tee $TEMP_DIR/$i.total
done

