#include <algorithm>
#include <iostream>
#include <string>
#include <sstream>
#include <cstdlib>
#include <vector>
#include <fstream>
#include <cstring>
#include <unordered_set>
#include <unordered_map>
#include <sys/time.h>
#include <sys/resource.h>
#include <sys/stat.h>
#include <unistd.h>
#include <errno.h>
#include <math.h>
#include "alias.h"
#include "BIT.h"
#include "rbtree.h"
#include "utils.h"

using namespace std;

typedef unsigned int uint;

struct node
{
	uint id;
	double w;
	node(uint _id, double _w) : id(_id), w(_w) {}
	node(const node &tmp)
	{
		id = tmp.id;
		w = tmp.w;
	}
	node &operator=(const node &tmp)
	{
		id = tmp.id;
		w = tmp.w;
		return *this;
	}
};

class Graph
{
public:
	uint n;
	string filedir, filelabel;
	unordered_map<uint, vector<node>> neighborList;
	unordered_map<uint, uint> outSizeList;
	unordered_map<uint, double> outWeightList;
	unordered_map<uint, double> outMaxWeightList;
	double totdeg;
	Graph() {}
	Graph(const string &_filedir, const string &_filelabel)
	{
		totdeg = 0;
		filedir = _filedir;
		filelabel = _filelabel;
		string neiNode, neiWeight, neiNum, graphAttr;
		graphAttr = filedir + filelabel + ".initattribute";
		ifstream graphAttrIn(graphAttr.c_str());
		cout << "FilePath: " << graphAttr.c_str() << endl;
		if (!graphAttrIn)
		{
			graphAttrIn.close();
			cout << "NO init file. Start to construct..." << endl;
			neiNode = filedir + filelabel + ".outEdges";
			neiWeight = filedir + filelabel + ".outWEdges";
			neiNum = filedir + filelabel + ".outPtr";
			graphAttr = filedir + filelabel + ".attribute";
			readFile(graphAttr, neiNode, neiWeight, neiNum);
			graphAttr = filedir + filelabel + ".initattribute";
			ofstream graphAttrOut(graphAttr.c_str());
			graphAttrOut << "n " << n;
			graphAttrOut.close();
			del(10000);
		}
		else
		{
			graphAttrIn.close();
			neiNode = filedir + filelabel + ".initoutEdges";
			neiWeight = filedir + filelabel + ".initoutWEdges";
			neiNum = filedir + filelabel + ".initoutPtr";
			readFile(graphAttr, neiNode, neiWeight, neiNum);
		}
		double pkm = peak_mem() / 1024.0 / 1024.0;
		cout << "Graph: peak memory: " << pkm << " G" << endl;
		double pkrss = peak_rss() / 1024.0 / 1024.0;
		cout << ", peak rss: " << pkrss << " G" << endl;
	}

	// void getQuery(vector<uint> &query, int num)
	// {
	// 	vector<double> outSizes(n, 0);
	// 	for (uint i = 0; i < n; i++)
	// 	{
	// 		uint outSize = outSizeList[i];
	// 		for (uint j = 0; j < outSize; j++)
	// 		{
	// 			uint nei = neighborList[i][j].id;
	// 			uint neiOutSize = outSizeList[nei];
	// 			outSizes[i] += neiOutSize - 1;
	// 			// for (uint k = 0; k < neiOutSize; k++)
	// 			// {
	// 			// 	uint neisnei = getOutVert(nei, k);
	// 			// 	if (neisnei != i)
	// 			// 		outSizes[i] += getOutSize(neisnei) - 1;
	// 			// }
	// 		}
	// 		outSizes[i] /= double(outSize);
	// 		outSizes[i] += outSize * 0.1 + outSizes[i] * 0.9;
	// 	}
	// 	vector<pair<double, uint>> pairs;
	// 	for (uint i = 0; i < n; i++)
	// 		pairs.push_back(make_pair(outSizes[i], i));
	// 	sort(pairs.begin(), pairs.end(), greater<pair<double, uint>>());
	// 	for (int i = 0; i < num; i++)
	// 	{
	// 		cout << pairs.at(i).first << endl;
	// 		query.push_back(pairs.at(i).second);
	// 	}
	// }

	void update()
	{
		ifstream opfile(filedir + "/" + filelabel + ".op", ios::in);

		double tottime = 0;
		int opnum = 0;
		uint s, t;
		double w;
		while (opfile >> s >> t >> w)
		{
			clock_t t0 = clock();
			add(s, t, w);
			clock_t t1 = clock();
			tottime += (t1 - t0) / (double)CLOCKS_PER_SEC;
			opnum++;
		}
		cout << filelabel << " avg update time: " << tottime / opnum << endl;
		opfile.close();
	}

	~Graph()
	{
	}

	void del(long num)
	{
		Random R;
		ofstream output(filedir + "/" + filelabel + ".op", ios::out);
		for (int i = 0; i < num; i++)
		{
			uint s = ceil(R.drand() * n);
			uint outSize = outSizeList[s];
			double outWeight = outWeightList[s];
			while (outSize <= 1)
			{
				s = ceil(R.drand() * n);
				outSize = outSizeList[s];
			}
			auto tmp = neighborList[s].end() - 1;
			neighborList[s].pop_back();
			outSizeList[s]--;
			outWeightList[s] -= tmp->w;
			output << s << " " << tmp->id << " " << tmp->w << endl;
		}
		output.close();
		ofstream outedges(filedir + filelabel + ".initoutEdges", ios::out | ios::binary);
		ofstream outwedges(filedir + filelabel + ".initoutWEdges", ios::out | ios::binary);
		ofstream outptr(filedir + filelabel + ".initoutPtr", ios::out | ios::binary);

		uint preSum = 0;
		for (uint i = 0; i < n; i++)
		{
			for (vector<node>::iterator it = neighborList[i].begin(); it != neighborList[i].end(); ++it)
			{
				outedges.write(reinterpret_cast<char *>(&(it->id)), sizeof(uint));
				outwedges.write(reinterpret_cast<char *>(&(it->w)), sizeof(double));
			}
			outptr.write(reinterpret_cast<char *>(&(preSum)), sizeof(uint));
			preSum += outSizeList[i];
		}
		outptr.write(reinterpret_cast<char *>(&(preSum)), sizeof(uint));
		outedges.close();
		outwedges.close();
		outptr.close();
	}

	void readFile(const string &graphAttr, const string &neiNode, const string &neiWeight, const string &neiNum)
	{
		cout << "Read graph attributes..." << endl;
		string tmp;
		ifstream graphAttrIn(graphAttr.c_str());
		graphAttrIn >> tmp >> n;
		cout << "n=" << n << endl;
		graphAttrIn.close();
		cout << "Read graph ..." << endl;
		//每个节点的出节点下标从哪里开始
		ifstream neiNumIn(neiNum.c_str(), ios::in | ios::binary);
		ifstream neiWeightIn(neiWeight.c_str(), ios::in | ios::binary);
		ifstream neiNodeIn(neiNode.c_str(), ios::in | ios::binary);
		uint outSize;
		uint outSizeSum = 0, preOutSizeSum = 0;
		neiNumIn.read(reinterpret_cast<char *>(&outSizeSum), sizeof(uint));
		for (uint i = 0; i < n; i++)
		{
			double outWeight = 0;
			neiNumIn.read(reinterpret_cast<char *>(&outSizeSum), sizeof(uint));
			outSize = outSizeSum - preOutSizeSum;
			preOutSizeSum = outSizeSum;
			double maxW = 0;
			for (uint j = 0; j < outSize; j++)
			{
				uint id;
				double w;
				neiNodeIn.read(reinterpret_cast<char *>(&id), sizeof(uint));
				neiWeightIn.read(reinterpret_cast<char *>(&w), sizeof(double));
				neighborList[i].push_back(node(id, w));
				outWeight += w;
				if (w > maxW)
					maxW = w;
			}
			outSizeList[i] = outSize;
			outWeightList[i] = outWeight;
			totdeg += outWeight;
			outMaxWeightList[i] = maxW;
		}
		neiNumIn.close();
		neiWeightIn.close();
		neiNodeIn.close();
	}

	virtual void add(uint s, uint t, double w)
	{
		neighborList[s].push_back(node{t, w});
		if (outMaxWeightList[s] < w)
			outMaxWeightList[s] = w;
		if (outSizeList.find(s) != outSizeList.end())
		{
			outSizeList[s]++;
			outWeightList[s] += w;
		}
		else
		{
			outSizeList[s] = 1;
			outWeightList[s] = w;
		}
	}
};

class BITPrefixSumGraph : public Graph
{
public:
	unordered_map<uint, BIT> BITList;
	void add(uint s, uint t, double w)
	{
		Graph::add(s, t, w);
		BITList[s].add(w);
	}
	BITPrefixSumGraph(const string &_filedir, const string &_filelabel) : Graph(_filedir, _filelabel)
	{
		for (uint i = 0; i < n; i++)
		{
			uint outSize = outSizeList[i];
			double *arr = new double[outSize];
			for (uint j = 0; j < outSize; j++)
				arr[j] = neighborList[i][j].w;
			BITList[i] = BIT(outSize, arr);
		}
	}
	BITPrefixSumGraph() {}
	~BITPrefixSumGraph() {}
};

class BSTPrefixSumGraph : public Graph
{
public:
	unordered_map<uint, bst> BSTList;
	void add(uint s, uint t, double w)
	{
		Graph::add(s, t, w);
		BSTList[s].insert(outWeightList[s]);
	}
	BSTPrefixSumGraph(const string &_filedir, const string &_filelabel) : Graph(_filedir, _filelabel)
	{
		for (uint i = 0; i < n; i++)
		{
			uint outSize = outSizeList[i];
			double prefixsum = 0;
			for (uint j = 0; j < outSize; j++)
			{
				prefixsum += neighborList[i][j].w;
				BSTList[i].insert(prefixsum);
			}
		}
	}
	BSTPrefixSumGraph() {}
	~BSTPrefixSumGraph() {}
};

class AliasMethodGraph : public Graph
{
public:
	unordered_map<uint, Alias> aliasList;
	void add(uint s, uint t, double w)
	{
		Graph::add(s, t, w);
		uint outSize = outSizeList[s];
		pair<int, double> *pi = new pair<int, double>[outSize];
		for (uint j = 0; j < outSize; j++)
		{
			pi[j] = make_pair(neighborList[s][j].id, neighborList[s][j].w);
		}
		aliasList[s] = Alias(pi, outSize);
	}
	AliasMethodGraph(const string &_filedir, const string &_filelabel) : Graph(_filedir, _filelabel)
	{
		for (uint i = 0; i < n; i++)
		{
			uint outSize = outSizeList[i];
			pair<int, double> *pi = new pair<int, double>[outSize];
			for (uint j = 0; j < outSize; j++)
			{
				pi[j] = make_pair(neighborList[i][j].id, neighborList[i][j].w);
			}
			aliasList[i] = Alias(pi, outSize);
		}
		cout << "sizeof Alias:" << sizeof(aliasList[0]) << endl;
	}
	AliasMethodGraph() {}
	~AliasMethodGraph() {}
};

class subsetGraph
{
public:
	uint n;
	string filedir, filelabel;
	unordered_map<uint, unordered_map<int, vector<node>>> neighborList;
	unordered_map<uint, uint> outSizeList;
	unordered_map<uint, double> outWeightList;
	subsetGraph() {}
	subsetGraph(const string &_filedir, const string &_filelabel)
	{
		filedir = _filedir;
		filelabel = _filelabel;
		string neiNode, neiWeight, neiNum, graphAttr;
		graphAttr = filedir + filelabel + ".initattribute";
		ifstream graphAttrIn(graphAttr.c_str());
		cout << "FilePath: " << graphAttr.c_str() << endl;
		if (!graphAttrIn)
		{
			graphAttrIn.close();
			cout << "NO init file. Start to construct..." << endl;
			neiNode = filedir + filelabel + ".outEdges";
			neiWeight = filedir + filelabel + ".outWEdges";
			neiNum = filedir + filelabel + ".outPtr";
			graphAttr = filedir + filelabel + ".attribute";
			readFile(graphAttr, neiNode, neiWeight, neiNum);
			graphAttr = filedir + filelabel + ".initattribute";
			ofstream graphAttrOut(graphAttr.c_str());
			graphAttrOut << "n " << n;
			graphAttrOut.close();
			del(10000);
		}
		else
		{
			graphAttrIn.close();
			neiNode = filedir + filelabel + ".initoutEdges";
			neiWeight = filedir + filelabel + ".initoutWEdges";
			neiNum = filedir + filelabel + ".initoutPtr";
			readFile(graphAttr, neiNode, neiWeight, neiNum);
		}
		double pkm = peak_mem() / 1024.0 / 1024.0;
		cout << "Graph: peak memory: " << pkm << " G" << endl;
		double pkrss = peak_rss() / 1024.0 / 1024.0;
		cout << ", peak rss: " << pkrss << " G" << endl;
	}

	~subsetGraph()
	{
	}

	void del(long num)
	{
		//是为了构造初始化的图，这里不实现了，先运行别的算法得到初始图文件即可。
	}

	void readFile(const string &graphAttr, const string &neiNode, const string &neiWeight, const string &neiNum)
	{
		cout << "Read graph attributes..." << endl;
		string tmp;
		ifstream graphAttrIn(graphAttr.c_str());
		graphAttrIn >> tmp >> n;
		cout << "n=" << n << endl;
		graphAttrIn.close();
		cout << "Read graph ..." << endl;
		//每个节点的出节点下标从哪里开始
		ifstream neiNumIn(neiNum.c_str(), ios::in | ios::binary);
		ifstream neiWeightIn(neiWeight.c_str(), ios::in | ios::binary);
		ifstream neiNodeIn(neiNode.c_str(), ios::in | ios::binary);
		uint outSize;
		uint outSizeSum = 0, preOutSizeSum = 0;
		neiNumIn.read(reinterpret_cast<char *>(&outSizeSum), sizeof(uint));
		for (uint i = 0; i < n; i++)
		{
			double outWeight = 0;
			neiNumIn.read(reinterpret_cast<char *>(&outSizeSum), sizeof(uint));
			outSize = outSizeSum - preOutSizeSum;
			preOutSizeSum = outSizeSum;
			vector<node> neighbor;
			for (uint j = 0; j < outSize; j++)
			{
				uint id;
				double w;
				neiNodeIn.read(reinterpret_cast<char *>(&id), sizeof(uint));
				neiWeightIn.read(reinterpret_cast<char *>(&w), sizeof(double));
				neighbor.push_back(node(id, w));
				outWeight += w;
			}
			outSizeList[i] = outSize;
			outWeightList[i] = outWeight;
			double gap = outWeight / outSize;
			for (uint j = 0; j < outSize; j++)
			{
				int subset = floor(log2(neighbor[j].w)) + 1;
				neighborList[i][subset].push_back(neighbor[j]);
			}
		}
		neiNumIn.close();
		neiWeightIn.close();
		neiNodeIn.close();
	}
	virtual void add(uint s, uint t, double w)
	{
		int subset = floor(log2(w)) + 1;
		neighborList[s][subset].push_back(node(t, w));
		if (outSizeList.find(s) != outSizeList.end())
		{
			outSizeList[s]++;
			outWeightList[s] += w;
		}
		else
		{
			outSizeList[s] = 1;
			outWeightList[s] = w;
		}
	}
	void update()
	{
		ifstream opfile(filedir + "/" + filelabel + ".op", ios::in);

		double tottime = 0;
		int opnum = 0;
		uint s, t;
		double w;
		while (opfile >> s >> t >> w)
		{
			clock_t t0 = clock();
			add(s, t, w);
			clock_t t1 = clock();
			tottime += (t1 - t0) / (double)CLOCKS_PER_SEC;
			opnum++;
		}
		cout << filelabel << " avg update time: " << tottime / opnum << endl;
		opfile.close();
	}
};
