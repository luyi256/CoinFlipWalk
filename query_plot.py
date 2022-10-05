import matplotlib
import matplotlib.pyplot as plt
import os
from datetime import datetime
import pytz
#! args
# 选择提取哪些数据集，会往本地写数据，文件：'.datalog/[date]_[measurement].data'
datasets = [ "coauth-MAG-History", "coauth-DBLP"]
# "threads-stack-overflow" "colisten-Spotify" "bitcoin-temporal" "indochina-2004" "twitter-2010","orkut-links" ,"tags-stack-overflow"
path = './analysis/'
measures = ['conductance', 'maxerr', 'precision']
#!! 算法是必须跑四个的，不然要牵一发而动全身
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
     "coauth-DBLP":"DBLP"
}
# algos = ['MCAM', 'MCSS', 'MCPS', 'MCAR']
# algos = ['AliasWalk', 'CoinFlipWalk', 'PrefixWalk', 'RejectionWalk']
# algo_alias = ['AliasWalk', 'CoinFlipWalk', 'PrefixWalk', 'RejectionWalk']
algos = ['CoinFlipWalk', 'PrefixWalk', 'RejectionWalk']
algo_alias = ['CoinFlipWalk', 'PrefixWalk', 'RejectionWalk']
ylabel = {
    'conductance': 'conductance',
    'maxerr': 'MaxAddErr',
    'precision': 'precision@50'
}
measureIdx = {'precision': 1, 'conductance': 2, 'maxerr': 3}
# colors = ['#263FA9', '#447B48', '#E85472', '#A686EC']
colors = ['#447B48', '#E85472', '#A686EC']
# markers = ["x", "s", "o", "+"]
markers = ["s", "o", "+"]
tz = pytz.timezone('Asia/Shanghai')
dt = datetime.now(tz).strftime("%m-%d")
figpath = path + '../result_query/' + dt + '/fig/'
if (os.path.exists(figpath) == False):
    os.makedirs(figpath)
for measure in measures:
    for dataset in datasets:
        fig, ax = plt.subplots()
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
            plt.plot(runtime, [i[measureIdx[measure]] for i in runerror],
                     label=algo,
                     marker=markers[k],
                     color=colors[k],
                     markersize=14,
                     markerfacecolor='none',
                     linewidth=2,
                     mew=2)
        if measureIdx[measure] == 3:
            ax.set_yscale('log')
        ax.set_xscale('log')
        plt.xlabel('query time(s)-' + datasets_alias[dataset], fontsize=15)
        ax.tick_params(axis='both', which='major', labelsize=15)
        plt.ylabel(measure + '-' + datasets_alias[dataset], fontsize=15)
        plt.legend(prop={'size': 10})
        plt.tight_layout()
        plt.savefig(
            figpath +
            "{}-query-{}_AddErr.jpg".format(measure, datasets_alias[dataset]))
        plt.savefig(
            figpath +
            "{}-query-{}_AddErr.pdf".format(measure, datasets_alias[dataset]))
