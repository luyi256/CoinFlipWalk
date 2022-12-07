import os
#!multiedge
path="/data/lu_yi/dataset/"
dataset="edit-enwiki"
outpath="../dataset/"+dataset+"_weighted/"
file="out."+dataset
good_d={}
with open(path+file,'r') as f:
  line=f.readline()
  line=f.readline()
  while line:
    edge = [int(i) for i in line.split()[0:2]]
    if edge[1] in good_d.keys():
      good_d[edge[1]].append(edge[0])
    else:
      good_d[edge[1]]=[edge[0],]
    line=f.readline()

graph={}
dd=[good_d]
for d in dd:
  for song in d.keys():
    size=len(d[song])
    for i in range(size-1):
      if d[song][i] in graph.keys():
        if d[song][i+1] in graph[d[song][i]].keys():
          graph[d[song][i]][d[song][i+1]]+=1
        else:
          graph[d[song][i]][d[song][i+1]]=1
      else:
        graph[d[song][i]]={}
        graph[d[song][i]][d[song][i+1]]=1
      if d[song][i+1] in graph.keys():
        if d[song][i] in graph[d[song][i+1]].keys():
          graph[d[song][i+1]][d[song][i]]+=1
        else:
          graph[d[song][i+1]][d[song][i]]=1
      else:
        graph[d[song][i+1]]={}
        graph[d[song][i+1]][d[song][i]]=1

if (os.path.exists(outpath) == False):
    os.makedirs(outpath)
f=open(outpath+dataset+".txt",'w')
for u in graph.keys():
  for v in graph[u].keys():
    f.write("{s} {dest} {w}\n".format(s=u,dest=v,w=graph[u][v]))

  

#! for wikiconflict
# path="../dataset/"
# dataset="wikiconflict"
# file="out."+dataset
# d={}
# min=1000000
# with open(path+file,'r') as f:
#   line=f.readline()
#   line=f.readline()
#   line=f.readline()
#   i=0
#   while line:
#     edge = [float(i) for i in line.split(' ')[0:3]]
#     edge[0]=int(edge[0])
#     edge[1]=int(edge[1])
#     if edge[2]==0:
#       edge[2]=0.1
#     if edge[2]*edge[2]<min:
#       min=edge[2]*edge[2]
#     if edge[0] in d.keys():
#       if edge[1] in d[edge[0]]:
#         d[edge[0]][edge[1]]+=edge[2]*edge[2]
#       else:
#         d[edge[0]][edge[1]]=edge[2]*edge[2]
#     else:
#       d[edge[0]]={}
#       d[edge[0]][edge[1]]=edge[2]*edge[2]
#     line=f.readline()
# print(min)
# outfile="wikiconflict_weighted.txt"
# f=open(path+outfile,'w')
# for u in d.keys():
#   for v in d[u].keys():
#     f.write("{s} {dest} {w}\n".format(s=u,dest=v,w=d[u][v]/min))


# opfile=dataset+".op"
# with open(path+opfile,'w') as f:
#   for i in range(op_num,0,-1):
#     f.write("{source} {dest}\n".format(source=l[-i][0],dest=l[-i][1]))
