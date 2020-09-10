#!/bin/bash
#$ -S /bin/bash
#$ -N time_check
#$ -cwd
#$ -pe smp 24
#$ -V


#### ================ Prepare dataset ==================


DATA=${PWD}/data
WO_CHUNKING_QUERYSIZE=500
QUERY_SIZE=250

## Download the reference
ECOLIREF=${DATA}/EColi_k12.fasta
if ! [ -e ${ECOLIREF} ]
then
    wget http://togows.dbcls.jp/entry/nucleotide/U00096.3.fasta -O ${ECOLIREF}
fi
REFERENCE_SIZE=200

## Donwload models.
if ! [ -d ${PWD}/kmer_models ]
then
    git clone https://github.com/nanoporetech/kmer_models.git
fi
MODELPATH=${PWD}/kmer_models/r9.2_180mv_250bps_6mer/template_median68pA.model


## Download signal data
## These are events extracted by Python script named `${PWD}/../score_calculate/scripts/extract.py`.
## It is obtained by
## ```bash
## wget https://s3.climb.ac.uk/nanopore/E_coli_K12_1D_R9.2_SpotON_2.tgz
## tar -xvf E_coli_K12_1D_R9.2_SpotON_2.tgz
## python3 ${PWD}/../score_calculate/scripts/extract.py \ 
## ${PWD}/E_coli_K12_1D_R9.2_SpotON_2/downloads/pass/ 1200 ${QUERIES}/events.json
## ```
## It requires ONT's fast5 API packages.
QUERY=${DATA}/events.json
if ! [ -e ${QUERY} ]
then
    wget https://mlab.cb.k.u-tokyo.ac.jp/~ban-m/read_until_paper/events.json.xz -O ${QUERY}.xz
    unxz ${QUERY}.xz
fi

## Download reads and SAM File. Mapping from query reads -> ECOLIREF.
READS=${DATA}/query.fasta
SAM=${DATA}/mapping.sam
if ! [ -e ${READS} ]
then
    wget https://s3.climb.ac.uk/nanopore/E_coli_K12_1D_R9.2_SpotON_2.pass.fasta -O ${READS}
fi

if ! [ -e ${SAM} ]
then
    minimap2 -a -x map-ont ${ECOLIREF} ${READS} > ${SAM}
fi


##### ==================== RUN ==========================

set -eu
### baseline implimentation
### Baseline(Sub dynamic time warping without chunking)
echo "refsize,bandwidth,time,power" > result_local/baseline_sub_without_chunking.csv
cargo run --release --bin speedcheck  -- $SAM $QUERY $MODELPATH $ECOLIREF $REFERENCE_SIZE $WO_CHUNKING_QUERYSIZE 0 0 \
      >> result_local/baseline_sub_without_chunking.csv 2> log
### Baseline(Sub dynamic time warping with chunking)
echo "refsize,bandwidth,time,power" > result_local/baseline_sub_with_chunking.csv
for power in $(seq 35 1 40)
do
    cargo run --release --bin speedcheck -- ${SAM} $QUERY $MODELPATH $ECOLIREF $REFERENCE_SIZE $QUERY_SIZE \
          0 ${power} >> result_local/baseline_sub_with_chunking.csv 2> log 
done


### Baseline v1.(Chunking and scouting)
echo "refsize,num_scouts,num_packs,time,power" > result_local/baseline_scouting_with_chunking.csv
for num_packs in $(seq 2 1 3)
do
    for num_scouts in $(seq 14 2 20)
    do
	    for power in $(seq 35 1 38)
        do
	        cargo run --release  --bin sub_scouting -- $SAM $QUERY $MODELPATH $ECOLIREF $REFERENCE_SIZE $QUERY_SIZE ${num_scouts} ${num_packs} ${power} >> result_local/baseline_scouting_with_chunking.csv 2> log
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
	        cargo run --release --bin scouting_threshold -- $SAM $QUERY $MODELPATH $ECOLIREF $REFERENCE_SIZE $QUERY_SIZE ${num_scouts} ${num_packs} ${power} >> result_local/proposed.csv 2>> log 
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
