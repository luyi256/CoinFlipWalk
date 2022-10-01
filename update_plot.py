from turtle import color
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from datetime import datetime
import pytz
import os
# args
# 选择提取哪些数据集，会往本地写数据，文件：'./result_update/[date]_[measurement].data'
datasets = ["twitter-2010"]
# "threads-stack-overflow" "colisten-Spotify" "bitcoin-temporal" "indochina-2004" "twitter-2010"
log_measures = ['memory', 'memory_overhead', 'time', 'time_overhead']
plot_measures = ['memory', 'memory_overhead', 'time', 'time_overhead']

algos = ['CoinFlipWalk', 'PrefixWalk', 'RejectionWalk']
algo_alias = ['CoinFlipWalk', 'PrefixWalk', 'RejectionWalk']
# algos = ['AliasWalk', 'CoinFlipWalk', 'PrefixWalk', 'RejectionWalk']
# algos = ['MCAM', 'MCSS', 'MCPS', 'MCAR']
# algo_alias = ['AliasWalk', 'CoinFlipWalk', 'PrefixWalk', 'RejectionWalk']
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
    "temporal-reddit-reply": 'RR'
}
tz = pytz.timezone('Asia/Shanghai')
# 避免生成的pdf在苹果系统中无法正确显示线形
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

# RejectionWalk的位置，求overhead的时候需要减去
to_sub = 0
for i, algo in enumerate(algos):
    if algo == 'MCAR' or algo == 'RejectionWalk':
        to_sub = i
        break
print("Rejection at",to_sub)
time = []
memory = []
for dataset in datasets:
    tmp_time = []
    tmp_memory = []
    for algo in algos:
        filename = './runlog/' + algo + '_' + dataset + '_update.log'
        tmp_del_time = 0
        tmp_add_time = 0
        with open(filename, 'r') as f:
            for i, line in enumerate(f.readlines()):
                if i == 4:  # del time
                    tmp_del_time = float(line.split()[5])
                if i == 5:
                    tmp_add_time = float(line.split()[5])
                if i == 6:
                    tmp_memory.append(float(line.split()[4]))
                    break
        tmp_time.append((tmp_del_time + tmp_add_time) / 2)
    time.append(tmp_time)
    memory.append(tmp_memory)
print(time)
print(memory)
time_overhead = [[j - i[to_sub] for j in i] for i in time]
memory_overhead = [[j - i[to_sub] for j in i] for i in memory]
print(time_overhead)
print(memory_overhead)

dt = datetime.now(tz).strftime("%m-%d")
path = "./result_update/" + dt + '/'
if (os.path.exists(path) == False):
    os.makedirs(path)
for measure in log_measures:
    isTime = False if measure.find("time") == -1 else True
    isOverhead = False if measure.find("overhead") == -1 else True
    if isTime and not isOverhead:
        filename = dt + "_" + "time.data"
        arr = time
    elif isTime and isOverhead:
        filename = dt + "_" + "time_overhead.data"
        arr = time_overhead
    elif not isTime and not isOverhead:
        filename = dt + "_" + "memory.data"
        arr = memory
    elif not isTime and isOverhead:
        filename = dt + "_" + "memory_overhead.data"
        arr = memory_overhead

    f = open(path + filename, 'w')
    for i, subarr in enumerate(arr):
        f.write(datasets[i] + ':')
        f.write(str(subarr))
        f.write('\n')
    f.close()

num_dataset = len(datasets)
# colors = [(0, 0.4470, 0.7410), (0.8500, 0.3250, 0.0980),
#           (0.4940, 0.1840, 0.5560), (0.6353, 0.07843, 0.1843)]
# colors = ['#263FA9', '#447B48', '#E85472', '#A686EC']
colors = [ '#447B48', '#E85472', '#A686EC']
# ylims = [50, 500, 40, 6]#8382EB
#EC7F4B
#263FA9
#447B48
for measure in plot_measures:
    isTime = False if measure.find("time") == -1 else True
    isOverhead = False if measure.find("overhead") == -1 else True
    if isTime and not isOverhead:
        arr = time
    elif isTime and isOverhead:
        arr = time_overhead
    elif not isTime and not isOverhead:
        arr = memory
    elif not isTime and isOverhead:
        arr = memory_overhead
    for k, dataset in enumerate(datasets):
        fig, ax = plt.subplots(figsize=(10, 5))
        ax.bar(algos,
               arr[k],
               tick_label=algo_alias,
               color='none',
               log=isTime,
               linewidth=3)
        ax.tick_params(axis='x', labelsize=25, rotation=15)
        plt.yticks(fontsize=18)
        if isTime and isOverhead:
            ax.set_ylabel("overhead time (s)", fontsize=25)
        elif isTime and not isOverhead:
            ax.set_ylabel("time (s)", fontsize=25)
        elif not isTime and isOverhead:
            ax.set_ylabel("overhead memory (GB)", fontsize=25)
        else:
            ax.set_ylabel("memory (GB)", fontsize=25)
        hatches = ["/", "x", ".", "\\"]
        ax.set_title(datasets_alias[dataset], fontsize=25)
        for i, patch in enumerate(ax.patches):
            patch.set_hatch(hatches[i])
            patch.set_edgecolor(colors[i])
        plt.rcParams['hatch.linewidth'] = 3.0
        plt.tight_layout()
        plt.savefig(path +
                    "{}_{}.png".format(measure, datasets_alias[dataset]))
        # plt.savefig(path + "{}_{}.pdf".format(measure, datasets_alias[dataset]))
