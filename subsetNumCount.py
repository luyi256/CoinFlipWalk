
import seaborn as sns
import matplotlib.pyplot as plt

# args
datasets = ["colisten-Spotify","threads-stack-overflow" ,"indochina-2004" ,"twitter-2010","orkut-links" ,"tags-stack-overflow" ,"temporal-reddit-reply"]

for dataset in datasets:
  counting={}
  with open("./runlog/"+dataset+".subsetNum",'r') as f:
    data=list(map(int,f.read().split()[:-1]))
  num_bin=range(min(data),max(data),1)
  print(num_bin)
  plt.hist((data))
  plt.xticks(num_bin)
  plt.grid(alpha=0.4)
  plt.savefig("./subsetNum/"+dataset+".jpg")
  # subsetList=[]
  # subsetNumList=[]
  # for key,value in counting.items():
  #   subsetList.append(key)
  #   subsetNumList.append(value)
  