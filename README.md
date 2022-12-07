# CoinFlipWalk

## Tested Environment
- Ubuntu 16.04.1 LTS (Xenial Xerus)
- g++ 7.5.0
- gcc 7.5.0


## Data
We provide the dataset threads-stack-overflow in the directory "./dataset/".

The other three datasets used in the paper can be downloaded from:
- bitcoin-temporal: https://www.cs.cornell.edu/~arb/data/temporal-bitcoin/
- indochina-2004: http://data.law.di.unimi.it/webdata/indochina-2004/indochina-2004.graph
- colisten-Spotify: https://www.cs.cornell.edu/~arb/data/colisten-Spotify/

Note that the dataset indochina-2004 downloaded from LAW website needs to be decompressed using the toolkit webgraph provided by the website(http://webgraph.di.unimi.it/).  

After downloading the indochina-2004, you need first convert it into motif-based weighted graphs using the MAPPR method given in [Yin et al, KDD2017]. The code of MAPPR is published at  http://snap.stanford.edu/mappr/ . Please put the motif-based weighted graph returned by MAPPR in the directory: `CoinFlipWalk/dataset/` and rename them as indochina-2004.txt. In addition, we provide a small dataset `grqc.txt` in `./dataset/grqc_weighted` directory for convenience.

We will convert the txt files into four files, `.attribute`, `.outEdges`, `.outPtr`, `.outWEdges`.

## Reproduce the experimental results
```bash
bash query.sh
```
Note that you may need to run powermethod first to generate groundtruths for performance evaluations (i.e. max-additive-error, precision@50, conductance). Because powermethod is a time-consuming algorithm, there may be a long time to wait in the beginning. 

A command for a method is as follows.
```bash
./<algo> -d </path/to/dataset/directory> -f <filelabel> -e <epsilon sequence> -qn <the number of query> -l <random walk step> [-u]
```

`-u` is optional. If `-u` is added, it will only run the updating part and exit. If not, it will only run the query.

For example, to run queries with `RejectionWalk` for `threads-stack-overflow`, you can run
```
./RejectionWalk -d ./dataset/threads-stack-overflow_weighted/ -f threads-stack-overflow -e 1e-1 1e-2 1e-3 1e-4 1e-5 -qn 10 -l 10
```
To get the performance evaluations, please run
```bash
bash metric.sh
```