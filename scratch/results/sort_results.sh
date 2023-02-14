#!/bin/bash

FILES=`ls *fil_0000000000.png | sed s/_0001\.fil_0000000000\.png//`

for FILE in $FILES
do
	mkdir -p $FILE
	mv $FILE* $FILE
done

