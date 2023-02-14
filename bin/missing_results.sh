#!/bin/bash

## This version does the symmetric difference (items that are unique in either)
#echo "$(ls Skynet*.fil | sed 's/_0001\.fil//'; ls ../results | grep Skynet)" | sort | uniq -u

#FIL_LIST=`ls -1 /users/nfairfie/scratch/fil/Skynet*.fil | xargs -L 1 basename | sed 's/_0001\.fil//' | sort`
#RESULTS_LIST=`ls -1 /users/nfairfie/scratch/results | xargs -L 1 basename | grep Skynet | sort`
#echo $FIL_LIST
#echo aaaaaaaaaaaaaaa
#echo $RESULTS_LIST

comm -23 <(ls -1 /users/nfairfie/scratch/fil/Skynet*.fil | xargs -L 1 basename | sed 's/_0001\.fil//' | sort) <(ls -1 /users/nfairfie/scratch/results | xargs -L 1 basename | grep Skynet | sort)


