#!/bin/bash
# This script searches Skynet data for long observations (up to *_0040.fits).

rm ~/long_observations.txt

cd ~/SkynetData
OBSERVATIONS=`ls -1`
for OBSERVATION in $OBSERVATIONS
do 
	NUM_FITS=`find $OBSERVATION -name \*_\?\?\?\?.fits | wc -l`
	if [ "$NUM_FITS" != "0" ]; then
		echo $NUM_FITS $OBSERVATION >> ~/long_observations.txt
	fi 
done

sort -r -n ~/long_observations.txt > ~/long_observations_sorted.txt 
