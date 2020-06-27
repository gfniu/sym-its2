#!/bin/bash

#
INPUT_PATH=''
OUTPUT_PATH=''

./cd-hit -c 1 -aS 1 -i $INPUT_PATH -G 0 -o $OUTPUT_PATH
./clstr_sort_by.pl no < ''$OUTPUT_PATH'.clstr' > ''$OUTPUT_PATH'.sorted'
