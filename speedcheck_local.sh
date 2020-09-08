#!/bin/bash

# This is a local runner for generating data for figure in paper.
# Note that you should correctly put the data on previous to run this script.

# Constant settings.
WO_CHUNKING_QUERYSIZE=500
QUERY_SIZE=250
ECOLIREF=/media/ban-m/hdd/linux-folder/readuntils/data/ecolik12.fa
MODELPATH=/media/ban-m/hdd/linux-folder/readuntils/kmer_models/r9.2_180mv_250bps_6mer/template_median68pA.model
QUERIES=/data/speedcheck/
REFERENCE_SIZE=200
TRAINING_DATA_ROOT="./data/200KScouting"

set -eu

### baseline implimentation
### Baseline(Sub dynamic time warping without chunking)
echo "refsize,bandwidth,time,power" > result_local/baseline_sub_without_chunking.csv
cargo run --release --bin speedcheck  -- $QUERIES $MODELPATH $ECOLIREF $REFERENCE_SIZE $WO_CHUNKING_QUERYSIZE 0 0 >> result_local/baseline_sub_without_chunking.csv 2> log

### Baseline(Sub dynamic time warping with chunking)
echo "refsize,bandwidth,time,power" > result_local/baseline_sub_with_chunking.csv
for power in $(seq 35 1 40)
do
    cargo run --release --bin speedcheck $QUERIES $MODELPATH $ECOLIREF $REFERENCE_SIZE $QUERY_SIZE 0 ${power} >> result_local/baseline_sub_with_chunking.csv 2> log 
done


### Baseline v1.(Chunking and scouting)
echo "refsize,num_scouts,num_packs,time,power" > result_local/baseline_scouting_with_chunking.csv
for num_packs in $(seq 2 1 3)
do
    for num_scouts in $(seq 14 2 20)
    do
	    for power in $(seq 35 1 38)
	    do
	        cargo run --release  --bin sub_scouting $QUERIES $MODELPATH $ECOLIREF $REFERENCE_SIZE $QUERY_SIZE ${num_scouts} ${num_packs} ${power} >> result_local/baseline_scouting_with_chunking.csv 2> log
	    done
    done
done


### Proposed method(Threshold + Scouting + Chunking)
echo "refsize,num_scouts,num_packs,result,power" > result_local/proposed.csv
for num_packs in $(seq 2 1 3)
do
    for num_scouts in $(seq 14 2 20)
    do
	    for power in $(seq 35 1 38)
	    do
	        cargo run --release --bin scouting_threshold -- $QUERIES $MODELPATH $ECOLIREF $REFERENCE_SIZE $QUERY_SIZE ${num_scouts} ${num_packs} ${power} ${TRAINING_DATA_ROOT}${num_scouts}${num_packs}${power}Positive.dat ${TRAINING_DATA_ROOT}${num_scouts}${num_packs}${power}Negative.dat >> result_local/proposed.csv 2>> log 
	    done
    done
done


### Baseline(Sakoe-Chiba)
echo "refsize,bandwidth,time,power" > result_local/baseline_sakoechiba_with_chunking.csv
for power in $(seq 35 1 40)
do
    for bandwidth in $(seq 11 5 91)
    do
        cargo run --release --bin speedcheck $QUERIES $MODELPATH $ECOLIREF $REFERENCE_SIZE $QUERY_SIZE ${bandwidth} ${power} >>result_local/baseline_sakoechiba_with_chunking.csv 2>> log
    done
done



### many times
echo "refsize,num_scouts,num_packs,result,power" > result_local/proposed_refsize_varies.csv
echo "refsize,bandwidth,time,power" > result_local/baseline_sub_wo_chunking_refsize_varies.csv
for num_packs in $(seq 2 1 3)
do
    for num_scouts in $(seq 14 2 20)
    do
	    
	for refsize in 2 4 8 16 32 64 128 256
	do
	    root="./data/"${refsize}"KScouting"
	    cargo run --release --bin speedcheck  -- $QUERIES $MODELPATH $ECOLIREF ${refsize} $QUERY_SIZE 0 0 >> result_local/baseline_sub_wo_chunking_refsize_varies.csv
	    for power in $(seq 35 1 38)
	    do
		cargo run --release --bin scouting_threshold -- $QUERIES $MODELPATH $ECOLIREF ${refsize} $QUERY_SIZE ${num_scouts} ${num_packs} ${power} ${root}${num_scouts}${num_packs}${power}Positive.dat ${root}${num_scouts}${num_packs}${power}Negative.dat >> result_local/proposed_refsize_varies.csv 2>> log 
	    done
	done
    done
done
