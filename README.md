# NeurIPS#5719: Edge-based Local Push for Personalized PageRank



## Tested Environment:
- Ubuntu 16.04.10
- g++ 11
- GCC 5.4.0


## Data:
We provide the dataset Youtube in the directory "./dataset/".

The other three real datasets used in the paper can be downloaded from:
(1) soc-LiveJournal1: http://snap.stanford.edu/data/soc-LiveJournal1.txt.gz
(2) indochina-2004: http://data.law.di.unimi.it/webdata/indochina-2004/indochina-2004.graph
(3) orkut: https://snap.stanford.edu/data/bigdata/communities/com-orkut.ungraph.txt.gz

Note that the dataset indochina-2004 downloaded from LAW website needs to be decompressed using the toolkit webgraph provided by the website(http://webgraph.di.unimi.it/).  
Four the four synthetic affinity graphs (i.e. AG_1, AG_2, AG_3, AG_4), we provide the codes to generate them. You can uncomment the line 6 or line 7 in the script file "run_script.sh" to run experiments on affinity graphs. 

After downloading the datasets, we need first convert them into motif-based weighted graphs using the MAPPR method given in [Yin et al, KDD2017]. The code of MAPPR is published at " http://snap.stanford.edu/mappr/ ". Please put the motif-based weighted graphs returned by MAPPR in the directory: "EdgePush-master/dataset/" and rename them as soc-LiveJournal1.txt, indochina-2004.txt and orkut-links.txt. In addition, we provide a small dataset "grqc.txt" in "dataset" directory for convenience. 


[Yin et al, KDD2017] Yin H, Benson A R, Leskovec J, et al. Local higher-order graph clustering[C]//Proceedings of the 23rd ACM SIGKDD international conference on knowledge discovery and data mining. 2017: 555-564.




## Reproduce the experimental results:
```
bash run_script.sh
```
Note that we may need to run powermethod first to generate PPR groundtruths for performance evaluations (i.e. l1Error, precision@50, conductance). Because powermethod is a time-consuming algorithm, there may be a long time to wait in the beginning. Our method EdgePush comes after powermethod. We regard the MAPPR method [Yin et al, KDD2017] as LocalPush and run the code of MAPPR published at "http://snap.stanford.edu/mappr/ ". 




## Parameters:
For simplicity and readability, we set the error parameter epsilon as its minimum value tested in the paper, corresponding to the last node shown in each experimental figure. We set the number of query node as 10 to help you see the experimental results more quickly. Note that the experimental reuslts with 10 or 50 query nodes nearly remain unchanged.    

The experimental results under other parameters can be deirved by the following ways.  

- -d \<directory\> 
- -f \<filelabel\>
- -algo \<algorithm\>
- [-e \<epsilon\> (default 1e-05)]
- [-qn \<querynum\> (default 10)]
- [-a \<alpha\> (default 0.2)]



## Instructions:
(1) datatset/: the directory to store datasets 
(2) query/: the directory to store query nodes
(3) result/: the directory to store PPR approximation results. 


备注：
- 数据集可以指定位置
- query就在当前文件夹下