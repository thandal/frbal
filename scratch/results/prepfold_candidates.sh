#!/bin/bash

# source /home/pulsar64/guppi/guppi.bash

#DATA=/home/than/Google Drive/Astronomy/Radio/data
DATA=/users/nfairfie/scratch/fil
DIRS=`ls -d *_heimdall`

for DIR in $DIRS
do
	cd $DIR	
	cat 2*.cand | sort -r -n > sorted.cand
	for N in {1..5};
	do
        	prepfold -noxwin -noclip -accelfile sorted.cand -accelcand $N -filterbank ${DATA}/${DIR:0:37}
		exit
   	done
	cd ..
done

