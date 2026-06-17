#!/bin/sh
LD_LIBRARY_PATH=/users/nfairfie/heimdall_build/dedisp/lib:/users/nfairfie/cuda-10.1/lib64
mkdir -p ${1}_candidates
/users/nfairfie/heimdall_build/heimdall/Applications/heimdall \
	-v -V \
	-nsamps_gulp 4194304 \
	-output_dir ${1}_candidates \
	-f $1


# Machines with GPUs
#  euclid
#  thales
#  20m-data
