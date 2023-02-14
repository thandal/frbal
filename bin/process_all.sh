#!/bin/bash

while true;
do
	FIRST_MISSING_RESULT=`missing_results.sh | head --lines=1`
	if [ -z "$FIRST_MISSING_RESULT" ]
	then
		exit
	fi
	# Note: dist_fdmt_search.sh will exit if a search is already running.
	dist_fdmt_search.sh ${FIRST_MISSING_RESULT}_0001.fil
	sleep 60
done
