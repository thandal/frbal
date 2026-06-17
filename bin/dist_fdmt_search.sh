#!/bin/bash
# Setup:
# * Install blimpy for your username:
#   $ pip install --user blimpy   # Seems to pick python3.7, and installs numpy, matplotlib, etc.
# * Install an ssh key
# For each observation:
# * On 20m-data, convert the .fits files to the much more compact .fil file with:
#   $ fits2fil.sh Skynet_XXXXX_OBSNAME_XXXXX
# * Run this script to dump the results plots:
#   $ dist_fdmt_search.sh Skynet_XXXXX_OBSNAME_XXXXX_0001.fil

# Note: some machines may be repeated, will double-up on worker tasks!
# 2021-10-15: removed euclid, its processes seem to hang...
# 2022-02-08: added euclid back in.
# 2022-02-08: added titania
MACHINES="euclid euclid euclid thales thales thales fourier fourier planck planck newton newton"
UNIQUE_MACHINES=`echo $MACHINES | tr ' ' '\n' | sort -u`

if [ 1 -ne $# ]
then
	echo "Argument required: 'kill', 'ps', or .fil filename"
	exit
fi

if [ "kill" == $1 ]
then
	for machine in $MACHINES
	do
		ssh -q $machine killall fdmt_search.py &
	done
	exit
fi

if [ "ps" == $1 ]
then
	for machine in $UNIQUE_MACHINES
	do
		echo "==== $machine ===="
		ssh -q $machine ps ux | grep fdmt_search | grep --invert-match grep
	done
	exit
fi

# When run on the host machine, there will be ssh sessions running fdmt_search.
DIST_PS=`ps u | grep ssh | grep fdmt_search | grep --invert-match grep | head --lines=1`

if [ "$DIST_PS" ]
then
	echo "A dist_fdmt_search appears to be already running: exiting."
	echo $DIST_PS
	exit
fi

# Cleanup
#find ~/scratch/results -size 0 | xargs rm

# Start the worker tasks.
for machine in $MACHINES
do
	echo "==== $machine ===="
	ssh $machine /users/nfairfie/bin/fdmt_search.py $1 &
done

# Machines -- https://www.gb.nrao.edu/pubcomputing/public.shtml
# Now -- https://greenbankobservatory.org/portal/gbt/processing/

#  euclid	Redhat EL8	128GB 2.1GHz 16 cores		Remote access only; lustre client, 2 x GPU
#  thales	Redhat EL8	128GB 2.1GHz 16 cores		Remote access only; lustre client, 1 x GPU
#  fourier	Redhat EL8	192GB 2.4GHz  8 cores		Remote access only; lustre client
#  planck	Redhat EL8	192GB 2.4GHz  8 cores		Remote access only; lustre client
#  newton	Redhat EL8	192GB 2.4GHz  8 cores		Remote access only; lustre client

# Others
#  titania                      128 GB, 2.5GHz 4 cores 		Observing machine! 1x Quadro P600
#  euler                        256 GB, 2.6GHz 32 cores  	"Connection closed by UNKNOWN port 65535"
#  maxwell                       64 GB, 8 cores
#  bratac                         8 GB, 3.4GHz 8 cores
#  kepler                        16 GB, 4 cores
#  20m-data                     32GB, 24 cores                  Also captures data from the 20m
#  sokar                        16 GB, 4 cores
#  korra                        16 GB, 4 cores
#  egret                         7 GB, 4 cores
#  belinda                       7 GB, 4 cores
#  cardano	Redhat EL8	24 GB, 8 cores			Remote access only, /home/scratch is local, no BeeGFS
# Unknown
#  blh0
#  mandc8

# NOT FOR GENERAL PURPOSE COMPUTING
#  arcturus	Redhat EL8	128GB, 2.1Ghz cores, 32 cores   Remote access only; lustre client for GBT pipeline runs
#  kepler	Redhat EL6	16GB, 2.4Ghz, 8 cores		Remote access only, /home/scratch is local, no lustre
#  galileo 64 GB, 24 cores Exclusive for Software Development Division

