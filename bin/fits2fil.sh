#!/bin/bash

# This script converts all PSRFITS .fits files from a given Skynet observation
# of the form:
#  Skynet_59590_M82_71548
# into .fil format in the ~/scratch/fits directory.
OBSERVATION=$1

OBSERVATION_PATH=~/SkynetData/$OBSERVATION
echo "OBSERVATION_PATH=$OBSERVATION_PATH"
if [ ! -d $OBSERVATION_PATH ];
then
	echo "Error: ~/SkynetData/$OBSERVATION_PATH not found, exiting."
	exit
fi

SCRATCH_PATH=~/scratch/fil
echo "SCRATCH_PATH=$SCRATCH_PATH"
mkdir -p $SCRATCH_PATH
cd $SCRATCH_PATH

# NOTE: there may be multiple collections of fits files, meaning multiple
# _0001.fits files!
FIRST_FITS=`ls $OBSERVATION_PATH/Cyborg/*_0001.fits`
echo "FIRST_FITS '$FIRST_FITS'"

for fits in $FIRST_FITS
do 
	echo "Processing $fits"
	if [ ! -f $fits ];
	then
		echo "Error: first .fits file not found, exiting."
		#echo `ls $OBSERVATION_PATH/Cyborg`
		exit
	fi
	
	#nice ~/bin/psrfits2fil $fits
	ionice --class 3 nice ~/bin/psrfits2fil $fits 695 950
	#nice ~/bin/psrfits2fil $fits 695 950
	#nice ~/bin/psrfits2fil $fits 747 910
	#nice ~/bin/psrfits2fil $fits 747 911
	#nice ~/bin/psrfits2fil $fits 747 911 flip 
done

#nice gzip *.fil
