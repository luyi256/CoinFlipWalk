from tkinter.tix import Tree
from tkinter.ttk import LabeledScale
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
datasets = ["tags-stack-overflow", "threads-stack-overflow", "colisten-Spotify","temporal-reddit-reply",]
# "indochina-2004", "tags-stack-overflow", "threads-stack-overflow", "colisten-Spotify", "twitter-2010" , "temporal-reddit-reply",
precision_idx = {
    "temporal-reddit-reply":5,
    "threads-stack-overflow":4,
    "tags-stack-overflow":5,
    "colisten-Spotify":4,
}
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
xlim={
    "temporal-reddit-reply":[4,1e-6],
    "threads-stack-overflow": [0.1,1e-6],
    "colisten-Spotify":[0,0],
    "tags-stack-overflow":[0,0]
}
tz = pytz.timezone('Asia/Shanghai')
# 避免生成的pdf在苹果系统中无法正确显示线形
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
markers = ["x", "*", "+"]
colors = ['#1F1DBF', '#E85472', '#447B48']
dt = datetime.now(tz).strftime("%m-%d")
path = "./result_precision/" + dt + '/'
if (os.path.exists(path) == False):
    os.makedirs(path)

for dataset in datasets:
    fig, ax = plt.subplots(figsize=(10,8))
    filename = './result_precision/readgraph_' + dataset + '_update.log'
    with open(filename, 'r') as f:
        add_idx=4
        del_idx=5
        memory_idx=6
        for i, line in enumerate(f.readlines()):
            if i == del_idx:  # del time
                tmp_del_time = float(line.split()[5])
                break
            if i == add_idx:
                tmp_add_time = float(line.split()[5])
        sub_update_time=(tmp_del_time + tmp_add_time) / 2
    for k,algo in enumerate(algos):
        print(dataset)
        print(algo)
        # read update
        filename = './result_precision/' + algo + '_' + dataset + '_update.log'
        if algo=='CoinFlipWalk':
            add_idx=5
            del_idx=6
            memory_idx=7
        else:
            add_idx=4
            del_idx=5
            memory_idx=6
        with open(filename, 'r') as f:
            for i, line in enumerate(f.readlines()):
                if i == del_idx:  # del time
                    tmp_del_time = float(line.split()[5])
                    break
                if i == add_idx:
                    tmp_add_time = float(line.split()[5])
        update_time=(tmp_del_time + tmp_add_time) / 2-sub_update_time
        # read query time
        query_file=open('./result_precision/' + algo + '_' + dataset + '_runtime.csv','r')
        query_time =  float(query_file.read().split(',')[precision_idx[dataset]-1])
        error_file=open('./result_precision/' + algo + '_' + dataset + '_runerror.csv','r')
        query_error=float(error_file.readlines()[precision_idx[dataset]-1].split(',')[1])
        plt.plot(query_time,update_time,marker=markers[k],color=colors[k],markersize=30,mew=3,markerfacecolor='none',)
        # plt.plot([query_time,query_time],[update_time,0],color=colors[k],linewidth=2.0,linestyle='--')
        # plt.plot([query_time,0],[update_time,update_time],color=colors[k],linewidth=2.0,linestyle='--')
        if datasets_alias[dataset]=='RR':
            plt.xlim((4, 18))
            plt.ylim((1.5e-6, 5.0e-6))
            if algo=='PrefixWalk':
                plt.text(query_time-0.25,update_time-5e-8,'({algo},{precision:.3f})'.format(algo=algo,precision=query_error),fontsize=28,color=colors[k],verticalalignment='top',horizontalalignment='right')
            elif algo=='CoinFlipWalk':
                plt.text(query_time+0.5,update_time-5e-8,'({algo},{precision:.3f})'.format(algo=algo,precision=query_error),fontsize=28,color=colors[k],verticalalignment='top',horizontalalignment='left')
            else:
                plt.text(query_time-0.5,update_time+8e-8,'({algo},{precision:.3f})'.format(algo=algo,precision=query_error),fontsize=28,color=colors[k],verticalalignment='bottom',horizontalalignment='left')
        elif datasets_alias[dataset]=='TH':
            plt.xlim((0.3, 1.0))
            plt.ylim((1.0e-6, 6.0e-6))
            if algo=='PrefixWalk':
                plt.text(query_time-0.007,update_time-8e-8,'({algo},{precision:.3f})'.format(algo=algo,precision=query_error),fontsize=28,color=colors[k],verticalalignment='top',horizontalalignment='right')
            elif algo=='CoinFlipWalk':
                plt.text(query_time+0.02,update_time+1e-8,'({algo},{precision:.3f})'.format(algo=algo,precision=query_error),fontsize=28,color=colors[k],verticalalignment='top',horizontalalignment='left')
            else:
                plt.text(query_time+0.008,update_time+8e-8,'({algo},{precision:.3f})'.format(algo=algo,precision=query_error),fontsize=28,color=colors[k],verticalalignment='bottom',horizontalalignment='left')
        elif datasets_alias[dataset]=='TAG':
            plt.xlim((-10, 120))
            plt.ylim((1e-7, 4.0e-6))
            if algo=='PrefixWalk':
                plt.text(query_time+2.1,update_time-6e-8,'({algo},{precision:.3f})'.format(algo=algo,precision=query_error),fontsize=28,color=colors[k],verticalalignment='top',horizontalalignment='left')
            elif algo=='CoinFlipWalk':
                plt.text(query_time+2.1,update_time+2e-7,'({algo},{precision:.3f})'.format(algo=algo,precision=query_error),fontsize=28,color=colors[k],verticalalignment='bottom',horizontalalignment='left')
            else:
                plt.text(query_time-1.3,update_time-6e-8,'({algo},{precision:.3f})'.format(algo=algo,precision=query_error),fontsize=28,color=colors[k],verticalalignment='top',horizontalalignment='right')
        elif datasets_alias[dataset]=='CS':
            plt.xlim((-10, 100))
            plt.ylim((4e-7, 1.2e-5))
            if algo=='PrefixWalk':
                plt.text(query_time+2.1,update_time+6e-8,'({algo},{precision:.3f})'.format(algo=algo,precision=query_error),fontsize=28,color=colors[k],verticalalignment='top',horizontalalignment='left')
            elif algo=='CoinFlipWalk':
                plt.text(query_time+2.1,update_time+5e-7,'({algo},{precision:.3f})'.format(algo=algo,precision=query_error),fontsize=28,color=colors[k],verticalalignment='bottom',horizontalalignment='left')
            else:
                plt.text(query_time-1.3,update_time-6e-8,'({algo},{precision:.3f})'.format(algo=algo,precision=query_error),fontsize=28,color=colors[k],verticalalignment='top',horizontalalignment='right')
        print(query_time)
        print(update_time)
        plt.xlabel('query time (s) - ' + datasets_alias[dataset], fontsize=36)
        plt.ylabel('update overhead (s) - ' + datasets_alias[dataset], fontsize=36)
    plt.tight_layout()
    ax.tick_params(axis='both', which='major', labelsize=28)
    ax.yaxis.get_offset_text().set(size=25)
    # ax.set_xlim(left=xlim[dataset][0])
    # ax.set_ylim(bottom=xlim[dataset][1])
    plt.grid(linestyle=':',linewidth=0.3,alpha=0.8)
    plt.tight_layout()
    plt.savefig(path +
                "{}.png".format(dataset))
    plt.savefig(path + "{}.pdf".format( datasets_alias[dataset]))
