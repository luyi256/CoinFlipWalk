import matplotlib
import matplotlib.pyplot as plt
import os
from datetime import datetime
import pytz
#! args
# 选择提取哪些数据集，会往本地写数据，文件：'.datalog/[date]_[measurement].data'
datasets = ["colisten-Spotify",]
# "indochina-2004", "tags-stack-overflow", "threads-stack-overflow", "colisten-Spotify", "twitter-2010" 
# "gottron-trec", "reuters", "bag-nytimes" "wikipedia-discussions-de" "edit-dewiki" "bag-pubmed"
# "indochina-2004" "tags-stack-overflow" "threads-stack-overflow" "colisten-Spotify" "temporal-reddit-reply"  "twitter-2010" "sx-stackoverflow" "CollegeMsg" "email-Eu-core-temporal"
# "youtube" "soc-LiveJournal1" "indochina-2004" "orkut-links" "tags-stack-overflow" "threads-stack-overflow" "blockchair" "colisten-Spotify" "bitcoin-temporal"
# "colisten-Spotify" "indochina-2004" "bitcoin-temporal""block_ethereum"
# "threads-stack-overflow" "colisten-Spotify" "bitcoin-temporal" "indochina-2004" "twitter-2010" "coauth-MAG" "coauth-AMiner" "temporal-reddit-reply" "colisten-Spotify" "twitter-2010" "coauth-MAG-History" "coauth-DBLP" "tags-ask-ubuntu"
# "revision_affinity_10000_1"  "revision_affinity_v2_10000_1"  "revision_affinity_10000_13" "revision_affinity_10000_20"
path = './result_query/10-13/' #  10-10_2_7039154 10-09_9189d8f
measures = ['conductance', 'maxerr', 'precision']
#! 常量设置
# 避免生成的pdf在苹果系统中无法正确显示线形
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
datasets_alias = {
    "twitter-2010": 'TW',
    "threads-stack-overflow": 'TH',
    "colisten-Spotify": 'CS',
    "indochina-2004": 'IC',
    "bitcoin-temporal": 'BT',
    "coauth-MAG": 'MAG',
    "coauth-AMiner": 'AM',
    "tags-stack-overflow": "TAG",
    "youtube": "YT",
    "soc-LiveJournal1": "LJ",
    "orkut-links": 'OL',
    "block_ethereum": 'ET',
    "blockchair": 'BC',
    "clickstream": "CL",
    "temporal-reddit-reply": 'RR',
    "revision_affinity_10000_1":'AG-1',
    "revision_affinity_10000_13":'AF_1_13',
    "revision_affinity_10000_20":'AF_1_20',
    "revision_affinity_v2_10000_1":'AG-2',
     "tags-ask-ubuntu":"ASK",
     "coauth-MAG-History":"MH",
     "coauth-DBLP":"DBLP",
      "sx-mathoverflow" :'MOF',
      "sx-stackoverflow":'SOF',
      "CollegeMsg":'CM',"email-Eu-core-temporal":'EEC',
      "libimseti":'LI',
      "wikiconflict":'WC',
      "yahoo-song":'YS',
      "wikiconflict_v2":'WC2',
      "gottron-trec":"GT",
       "reuters" :'RT',
       "bag-nytimes":'BN',
       "wikipedia-discussions-de":'WD',
       "bag-pubmed":'BP',
       "edit-dewiki":'ED',
       "yahoo-song-v2":'YS2',
       "movielens-10m_rating":'MR',
       "movielens-10m_rating-v2":'MR2',
       "yahoo-song-v3":'YS3'
}
# algos = ['MCAM', 'MCSS', 'MCPS', 'MCAR']
# algos = ['AliasWalk', 'CoinFlipWalk', 'PrefixWalk', 'RejectionWalk']
# algo_alias = ['AliasWalk', 'CoinFlipWalk', 'PrefixWalk', 'RejectionWalk']
algos = ['CoinFlipWalk', 'PrefixWalk', 'RejectionWalk']
algo_alias = ['CoinFlipWalk', 'PrefixWalk', 'RejectionWalk']
ylabel = {
    'conductance': 'Conductance',
    'maxerr': 'MaxError',
    'precision': 'Precision@50'
}
measureIdx = {'precision': 1, 'conductance': 2, 'maxerr': 3}
# colors = ['#263FA9', '#447B48', '#E85472', '#A686EC']
colors = ['#1F1DBF', '#E85472', '#447B48']
# markers = ["x", "s", "o", "+"]
markers = ["x", "*", "+"]
tz = pytz.timezone('Asia/Shanghai')
dt = datetime.now(tz).strftime("%m-%d")
figpath = './result_query/' + dt + '/fig/'
if (os.path.exists(figpath) == False):
    os.makedirs(figpath)
for measure in measures:
    for dataset in datasets:
        fig, ax = plt.subplots(figsize=(10,8))
        for k, algo in enumerate(algos):
            print(dataset)
            print(algo)
            runtime_file = open(path + algo + '_' + dataset + '_runtime.csv',
                                'r')
            runerror_file = open(path + algo + '_' + dataset + '_runerror.csv',
                                 'r')
            runtime = [float(i) for i in runtime_file.read().split(',')[:-1]]
            runerror = [[float(j) for j in i.split(',')]
                        for i in runerror_file.read().split('\n')[:-1]]
            print(runtime)
            print(runerror)
            # print([i[measureIdx[measure]] for i in runerror])
            if dataset=='colisten-Spotify':
                idx=[0,4]
            else:
                idx=[0,5]
            plt.plot(runtime[idx[0]:idx[1]], [i[measureIdx[measure]] for i in runerror][idx[0]:idx[1]],
                     label=algo,
                     marker=markers[k],
                     color=colors[k],
                     markersize=25,
                     markerfacecolor='none',
                     linewidth=3,
                     mew=3)
        if measureIdx[measure] == 3:
            ax.set_yscale('log')
        ax.set_xscale('log')
        ax.set_yticklabels([],minor=True)
        plt.xlabel('query time (s) - ' + datasets_alias[dataset], fontsize=36)
        ax.tick_params(axis='both', which='major', labelsize=28)
        plt.ylabel(ylabel[measure] + ' - ' + datasets_alias[dataset], fontsize=36)
        plt.legend(prop={'size': 25},fontsize=35)
        plt.tight_layout()
        plt.savefig(
            figpath +
            "{}-query-{}_AddErr.jpg".format(measure, datasets_alias[dataset]))
        plt.savefig(
            figpath +
            "{}-query-{}_AddErr.pdf".format(measure, datasets_alias[dataset]))
