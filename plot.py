from turtle import color
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

# strs = ["memory_overhead", "memory", "time", 'time_overhead']
strs = ["time"]
dataset = ['BT', 'CS', 'IC', 'TH']
algo = ['AliasWalk', 'CoinFlipWalk', 'PrefixWalk', 'RejectionWalk']
num_dataset = len(dataset)
# colors = [(0, 0.4470, 0.7410), (0.8500, 0.3250, 0.0980),
#           (0.4940, 0.1840, 0.5560), (0.6353, 0.07843, 0.1843)]
colors = ['#263FA9', '#447B48', '#E85472', '#A686EC']
# ylims = [50, 500, 40, 6]#8382EB
#EC7F4B
#263FA9
#447B48
for str in strs:
    isTime = False if str.find("time") == -1 else True
    isOverhead = False if str.find("overhead") == -1 else True
    file = './fig/' + str + '.txt'
    table = pd.read_csv(file, header=None)
    array = np.array(table)
    newarr = np.array([array[2, :], array[3, :], array[1, :], array[0, :]])
    for k in range(num_dataset):
        fig, ax = plt.subplots(figsize=(10, 5))
        print(type(newarr[0, k]))
        ax.bar(algo,
               newarr[:, k],
               tick_label=algo,
               color='none',
               log=isTime,
               lineWidth=3)
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
        ax.set_title(dataset[k], fontsize=25)
        for i, patch in enumerate(ax.patches):
            patch.set_hatch(hatches[i])
            patch.set_edgecolor(colors[i])
        plt.rcParams['hatch.linewidth'] = 3.0
        plt.tight_layout()

        plt.savefig("./fig/{}_{}.png".format(str, dataset[k]))
        plt.savefig("./fig/{}_{}.pdf".format(str, dataset[k]))
