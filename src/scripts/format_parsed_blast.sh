#!/bin/bash
PARSED_BLAST=$1

sed -z 's/\n>/;>/g' $PARSED_BLAST | sed -z 's/\n/ /g' | sed -z 's/;/\n/g' | sed -z 's/\t/;/g' > "fix_${PARSED_BLAST}"
