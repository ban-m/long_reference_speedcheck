#!/bin/bash

# This is a local runner for generating data for figure in paper.
# Note that you should correctly put the data before running this script.

# Constant settings. You may want to modify them to adjust your environment.
QUERY_SIZE=250

## Where your reference sequence is.
ECOLIREF=/media/ban-m/hdd/linux-folder/readuntils/data/ecolik12.fa

## Whare the model file is.
MODELPATH=/media/ban-m/hdd/linux-folder/readuntils/kmer_models/r9.2_180mv_250bps_6mer/template_median68pA.model

## Where the .fast5 file is.(More precisely, the runner script searches with the gblob: $QUERIES/*.fast5)
QUERIES=/data/speedcheck/

## Reference size length in kbp. For example, 200 means to use first 200 kbp region.
REFERENCE_SIZE=200


## The directry where training dataset is in. In fact, traning data is just a sequence of DTW score separated by newline. So you can easily create it with any method you like.
## Please see README.md for more detail.
TRAINING_DATA_ROOT="./data/data/scouting_training/200Kscouting"

set -eu

echo "refsize,num_scouts,num_packs,result,power" > proposed.csv
for num_packs in $(seq 2 1 3)
do
    for num_scouts in $(seq 14 2 20)
    do
	    for power in $(seq 35 1 38)
	    do
	        cargo run --release --bin scouting_threshold -- $QUERIES $MODELPATH $ECOLIREF $REFERENCE_SIZE $QUERY_SIZE ${num_scouts} ${num_packs} ${power} ${TRAINING_DATA_ROOT}${num_scouts}${num_packs}${power}Positive.dat ${TRAINING_DATA_ROOT}${num_scouts}${num_packs}${power}Negative.dat >> proposed.csv 2>> log 
	    done
    done
done
