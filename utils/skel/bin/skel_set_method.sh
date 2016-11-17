#!/bin/bash

# Parameters are:
#       1) Method, in the correct case to fit the XML file format, e.g. POSIX or MPI_LUSTRE
#       2) Infile
#       3) Outfile
#       4) Method parameter string


# Two sed commands (-e): the first one replaces the XXX in each method="XXX" with the new method
# The second sed replaces all instances of the string '***skel-parameters***' with the new parameter string 
sed -e's+method="[a-zA-Z0-9\-\_]*"+method="'$1'"+g' \
    -e's+\*\*\*skel\-parameters\*\*\*+'"$4"'+g' < $2 > $3


