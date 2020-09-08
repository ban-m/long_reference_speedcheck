# Long reference speedcheck: a program to record the time needed to process a query.

- Author: Bansho Masutani<banmasutani@gmail.com>
- Data: 2017-11-19
- Lang: Rust


This is a repository to record the elapsed time to process a query by Dyss algorithm.

To replicat our result in the paper, please follow the instruction below:
0. Install python3 and Rust language. If you have not installed Rust yet, just type:`curl https://sh.rustup.rs -sSf | sh`.
1. Clone and build [Dyss](https://bitbucket.org/ban-m/dyss/src/default/).
```bash
hg clone https://ban-m@bitbucket.org/ban-m/dyss
bash setup.sh
python3 ./src/dyss_debug.py --reference ./data/lambda.fa --model ./kmer_models/r9.4_180mv_450bps_6mer/template_median68pA.model --param ./data/parameters.csv --test ./data/test_reads/
```
2. In the parent directly of Dyss(i.e. after cd ../ in ./dyss), clone this repository.
3. Build this repository with `cargo build --release`
4. Run DySS package.
```
cargo run --release --bin scouting_threshold \
                    -- [QUERIES] [MODELPATH] [REFERENCE] [REFERENCE_SIZE] \
                    [QUERY_SIZE] [num_scouts] [num_packs] [power] \
                    [Dataset for Positive Score] \
                    [Dataset for Negative Score] \
                    >> [Result]
```

4. Edit speedcheck_only_dyss.sh so that it can run in your environment. For example, you can modify the MODELPATH variable to the directly to ../dyss/kmer_models/r9.2_180mv_250bps_6mer/template_median68pA.model(6-mer model for R9.2 chemistry provided by ONT).
4. Run the script by `bash speedcheck_only_dyss.sh`.

If you want to experiment with various reference size, you should first create the training data set. This is just a file which containts DTW scores separated by newline.
You can give the dataset file via command line argument.


