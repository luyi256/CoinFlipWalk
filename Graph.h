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
#include "AVL.h"
#include "rbtree.h"
#include "utils.h"
#include <chrono>
using namespace std;

typedef unsigned int uint;
const int add_num = 100000;
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
	unordered_map<uint, unordered_map<uint, int>> adjList;
	unordered_map<uint, uint> outSizeList;
	unordered_map<uint, double> maxWeight;
	Graph() {}
	Graph(const string &_filedir, const string &_filelabel)
	{
		filedir = _filedir;
		filelabel = _filelabel;
		string neiNode, neiWeight, neiNum, graphAttr;
		neiNode = filedir + filelabel + ".outEdges";
		neiWeight = filedir + filelabel + ".outWEdges";
		neiNum = filedir + filelabel + ".outPtr";
		graphAttr = filedir + filelabel + ".attribute";
		cout << "FilePath: " << graphAttr.c_str() << endl;
		readFile(graphAttr, neiNode, neiWeight, neiNum);
		ifstream opfile(filedir + "/" + filelabel + ".op");
		if (!opfile)
			getAddEdge(add_num);
		opfile.close();
	}

	void getQuery(const string &queryname)
	{
		cout << "Generate query file..." << endl;
		pair<int, double> *aliasD = new pair<int, double>[n];
		for (uint i = 0; i < n; i++)
		{
			double outweight = 0.0;
			for (uint j = 0; j < outSizeList[i]; j++)
			{
				outweight += neighborList[i][j].w;
			}
			aliasD[i] = make_pair(i, outweight);
		}
		Alias alias = Alias(aliasD, n);
		Random R;
		ofstream query(queryname);
		for (int i = 0; i < 50; i++)
		{
			uint node = alias.generateRandom(R);
			query << node << endl;
		}
		query.close();
		delete[] aliasD;
	}

	void update()
	{
		ifstream opfile(filedir + "/" + filelabel + ".op", ios::in);
		uint sarr[add_num];
		uint tarr[add_num];
		double warr[add_num];
		int opnum = 0;
		while (opfile >> sarr[opnum] >> tarr[opnum] >> warr[opnum])
		{
			opnum++;
		}
		opfile.close();
		// delete
		auto begin = chrono::high_resolution_clock::now();
		for (int i = 0; i < add_num; i++)
		{
			uint s = sarr[i];
			uint t = tarr[i];
			if (adjList[s].find(t) == adjList[s].end())
			{
				cout << "delete a nonexistent neighbor" << endl;
				exit(-2);
			}
			uint neiidx = adjList[s][t];
			uint outsize = outSizeList[s];
			if (neiidx == outsize - 1)
				neighborList[s].pop_back();
			else
			{
				auto tmp = neighborList[s].end() - 1;
				neighborList[s][neiidx].id = tmp->id;
				neighborList[s][neiidx].w = tmp->w;
				adjList[s][tmp->id] = neiidx;
				neighborList[s].pop_back();
			}
			adjList[s].erase(t);
			outSizeList[s]--;
		}
		cout << filelabel << " new del update time: " << (chrono::duration_cast<chrono::nanoseconds>(chrono::high_resolution_clock::now() - begin).count() / 1000000000.0) / add_num << endl;
		// add
		begin = chrono::high_resolution_clock::now();
		for (int i = 0; i < add_num; i++)
		{
			neighborList[sarr[i]].push_back(node{tarr[i], warr[i]});
			outSizeList[sarr[i]]++;
			adjList[sarr[i]][tarr[i]] = neighborList[sarr[i]].size() - 1;
		}
		cout << filelabel << " new add update time: " << (chrono::duration_cast<chrono::nanoseconds>(chrono::high_resolution_clock::now() - begin).count() / 1000000000.0) / add_num << endl;
	}

	~Graph()
	{
	}

	void getAddEdge(long num)
	{
		pair<int, double> *aliasD = new pair<int, double>[n];
		unordered_map<uint, vector<uint>> hasAdded;
		for (uint i = 0; i < n; i++)
		{
			aliasD[i] = make_pair(i, outSizeList[i]);
		}
		Alias alias = Alias(aliasD, n);
		Random R;
		ofstream output(filedir + "/" + filelabel + ".op", ios::out);
		int i = num;
		while (i)
		{
			uint nodeidx = alias.generateRandom(R);
			uint outSize = outSizeList[nodeidx];
			while (outSize <= 1)
			{
				nodeidx = alias.generateRandom(R);
				outSize = outSizeList[nodeidx];
			}
			uint neiidx = floor(R.drand() * outSize);
			auto tmp = neighborList[nodeidx][neiidx];
			if (hasAdded[nodeidx].size() > 0)
			{
				int flag = 0;
				for (auto j : hasAdded[nodeidx])
					if (j == tmp.id)
					{
						flag = 1;
						break;
					}
				if (flag == 1)
					continue;
			}
			output << nodeidx << " " << tmp.id << " " << tmp.w << endl;
			hasAdded[nodeidx].push_back(tmp.id);
			i--;
		}
		delete[] aliasD;
		output.close();
		// ofstream outedges(filedir + filelabel + ".initoutEdges_distr", ios::out | ios::binary);
		// ofstream outwedges(filedir + filelabel + ".initoutWEdges_distr", ios::out | ios::binary);
		// ofstream outptr(filedir + filelabel + ".initoutPtr_distr", ios::out | ios::binary);

		// uint preSum = 0;
		// for (uint i = 0; i < n; i++)
		// {
		// 	for (vector<node>::iterator it = neighborList[i].begin(); it != neighborList[i].end(); ++it)
		// 	{
		// 		outedges.write(reinterpret_cast<char *>(&(it->id)), sizeof(uint));
		// 		outwedges.write(reinterpret_cast<char *>(&(it->w)), sizeof(double));
		// 	}
		// 	outptr.write(reinterpret_cast<char *>(&(preSum)), sizeof(uint));
		// 	preSum += outSizeList[i];
		// }
		// outptr.write(reinterpret_cast<char *>(&(preSum)), sizeof(uint));
		// outedges.close();
		// outwedges.close();
		// outptr.close();
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
		ifstream neiNumIn(neiNum.c_str(), ios::in | ios::binary);
		// ifstream neiWeightIn(neiWeight.c_str(), ios::in | ios::binary);
		ifstream neiNodeIn(neiNode.c_str(), ios::in | ios::binary);
		uint outSize;
		uint outSizeSum = 0, preOutSizeSum = 0;
		neiNumIn.read(reinterpret_cast<char *>(&outSizeSum), sizeof(uint));
		for (uint i = 0; i < n; i++)
		{
			neiNumIn.read(reinterpret_cast<char *>(&outSizeSum), sizeof(uint));
			outSize = outSizeSum - preOutSizeSum;
			preOutSizeSum = outSizeSum;
			double maxw = 0;
			for (uint j = 0; j < outSize; j++)
			{
				uint id;
				double w;
				neiNodeIn.read(reinterpret_cast<char *>(&id), sizeof(uint));
				// neiWeightIn.read(reinterpret_cast<char *>(&w), sizeof(double));
				w = 1;
				if (w > maxw)
					maxw = w;
				neighborList[i].push_back(node(id, w));
				adjList[i][id] = j;
			}
			outSizeList[i] = outSize;
			maxWeight[i] = maxw;
		}
		neiNumIn.close();
		// neiWeightIn.close();
		neiNodeIn.close();
	}
};

// class BITPrefixSumGraph : public Graph
// {
// public:
// 	unordered_map<uint, BIT> BITList;
// 	void add(uint s, uint t, double w)
// 	{
// 		Graph::add(s, t, w);
// 		BITList[s].add(w);
// 	}
// 	BITPrefixSumGraph(const string &_filedir, const string &_filelabel) : Graph(_filedir, _filelabel)
// 	{
// 		for (uint i = 0; i < n; i++)
// 		{
// 			uint outSize = outSizeList[i];
// 			double *arr = new double[outSize];
// 			for (uint j = 0; j < outSize; j++)
// 				arr[j] = neighborList[i][j].w;
// 			BITList[i] = BIT(outSize, arr);
// 		}
// 	}
// 	BITPrefixSumGraph() {}
// 	~BITPrefixSumGraph() {}
// };

class BSTPrefixSumGraph
{
public:
	unordered_map<uint, bst> BSTList;
	uint n;
	string filedir, filelabel;
	unordered_map<uint, vector<node>> neighborList;
	unordered_map<uint, unordered_map<uint, int>> adjList;
	unordered_map<uint, uint> outSizeList;
	unordered_map<uint, double> outWeightList;
	BSTPrefixSumGraph(const string &_filedir, const string &_filelabel)
	{
		filedir = _filedir;
		filelabel = _filelabel;
		ifstream opfile(filedir + "/" + filelabel + ".op");
		string neiNode, neiWeight, neiNum, graphAttr;
		neiNode = filedir + filelabel + ".outEdges";
		neiWeight = filedir + filelabel + ".outWEdges";
		neiNum = filedir + filelabel + ".outPtr";
		graphAttr = filedir + filelabel + ".attribute";
		cout << "FilePath: " << graphAttr.c_str() << endl;
		readFile(graphAttr, neiNode, neiWeight, neiNum);
		if (!opfile)
			getAddEdge(10000);
	}

	void update()
	{
		ifstream opfile(filedir + "/" + filelabel + ".op", ios::in);
		uint sarr[add_num];
		uint tarr[add_num];
		double warr[add_num];
		int opnum = 0;
		while (opfile >> sarr[opnum] >> tarr[opnum] >> warr[opnum])
		{
			opnum++;
		}
		opfile.close();
		// delete
		// auto begin = chrono::high_resolution_clock::now();
		// for (int i = 0; i < add_num; i++)
		// {
		// 	uint s = sarr[i];
		// 	uint t = tarr[i];
		// 	uint neiidx = adjList[s][t];
		// 	uint outsize = outSizeList[s];
		// 	if (neiidx == outsize - 1)
		// 		neighborList[s].pop_back();
		// 	else
		// 	{
		// 		auto tmp = neighborList[s].end() - 1;
		// 		neighborList[s][neiidx].id = tmp->id;
		// 		neighborList[s][neiidx].w = tmp->w;
		// 		adjList[s][tmp->id] = neiidx;
		// 		neighborList[s].pop_back();
		// 		adjList[s].erase(t);
		// 	}
		// 	outSizeList[s]--;
		// 	outWeightList[s] -= warr[i];
		// 	// BST无法处理删掉一个值的情况，因为存储的是前缀和，这样前面所有值都要改变，这里应该只能通过树状数组完成
		// 	BSTList[s].delete_value(outWeightList[s]);
		// }
		// cout << filelabel << "avg del time: " << (chrono::duration_cast<chrono::nanoseconds>(chrono::high_resolution_clock::now() - begin).count() / 1000000000.0) / add_num << endl;
		// add
		auto begin = chrono::high_resolution_clock::now();
		for (int i = 0; i < add_num; i++)
		{
			neighborList[sarr[i]].push_back(node{tarr[i], warr[i]});
			outSizeList[sarr[i]]++;
			outWeightList[sarr[i]] += warr[i];
			BSTList[sarr[i]].insert(outWeightList[sarr[i]]);
			adjList[sarr[i]][tarr[i]] = neighborList[sarr[i]].size() - 1;
		}
		cout << filelabel << " new add update time: " << (chrono::duration_cast<chrono::nanoseconds>(chrono::high_resolution_clock::now() - begin).count() / 1000000000.0) / add_num << endl;
	}

	void getAddEdge(long num)
	{
		cout << "No update file. Try MCAR first and then run MCPS" << endl;
		exit(-1);
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
		ifstream neiNumIn(neiNum.c_str(), ios::in | ios::binary);
		// ifstream neiWeightIn(neiWeight.c_str(), ios::in | ios::binary);
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
			double prefixsum = 0;
			for (uint j = 0; j < outSize; j++)
			{
				uint id;
				double w;
				neiNodeIn.read(reinterpret_cast<char *>(&id), sizeof(uint));
				// neiWeightIn.read(reinterpret_cast<char *>(&w), sizeof(double));
				w = 1;
				neighborList[i].push_back(node(id, w));
				outWeight += w;
				prefixsum += neighborList[i][j].w;
				adjList[i][id] = j;
				BSTList[i].insert(prefixsum);
			}
			outSizeList[i] = outSize;
			outWeightList[i] = outWeight;
		}
		neiNumIn.close();
		// neiWeightIn.close();
		neiNodeIn.close();
	}

	BSTPrefixSumGraph() {}
	~BSTPrefixSumGraph() {}
};

class AVLPrefixSumGraph
{
public:
	unordered_map<uint, AVLnode *> AVLList;
	uint n;
	string filedir, filelabel;
	unordered_map<uint, vector<node>> neighborList;
	unordered_map<uint, unordered_map<uint, int>> adjList;
	unordered_map<uint, uint> outSizeList;
	unordered_map<uint, double> outWeightList;
	AVLPrefixSumGraph(const string &_filedir, const string &_filelabel)
	{
		filedir = _filedir;
		filelabel = _filelabel;
		ifstream opfile(filedir + "/" + filelabel + ".op");
		string neiNode, neiWeight, neiNum, graphAttr;
		neiNode = filedir + filelabel + ".outEdges";
		neiWeight = filedir + filelabel + ".outWEdges";
		neiNum = filedir + filelabel + ".outPtr";
		graphAttr = filedir + filelabel + ".attribute";
		cout << "FilePath: " << graphAttr.c_str() << endl;
		readFile(graphAttr, neiNode, neiWeight, neiNum);
		if (!opfile)
			getAddEdge(10000);
	}

	void update()
	{
		ifstream opfile(filedir + "/" + filelabel + ".op", ios::in);
		uint sarr[add_num];
		uint tarr[add_num];
		double warr[add_num];
		int opnum = 0;
		while (opfile >> sarr[opnum] >> tarr[opnum] >> warr[opnum])
		{
			opnum++;
		}
		opfile.close();
		// delete
		auto begin = chrono::high_resolution_clock::now();
		for (int i = 0; i < add_num; i++)
		{
			uint s = sarr[i];
			uint t = tarr[i];
			double w = warr[i];
			if (adjList[s].find(t) == adjList[s].end())
			{
				cout << "delete a nonexistent neighbor" << endl;
				exit(-2);
			}
			uint neiidx = adjList[s][t];
			// if (AVLList[s]->leftCount + AVLList[s]->rightCount != outSizeList[s] || (AVLList[s]->leftSum + AVLList[s]->rightSum - outWeightList[s]) > 1e-6)
			// {
			// 	cout << "error" << endl;
			// 	exit(-1);
			// }
			uint outsize = outSizeList[s];
			if (neiidx == outsize - 1)
			{
				neighborList[s].pop_back();
				AVLList[s] = deleteLast(AVLList[s]);
			}
			else
			{
				auto tmp = neighborList[s].end() - 1;
				neighborList[s][neiidx].id = tmp->id;
				neighborList[s][neiidx].w = tmp->w;
				adjList[s][tmp->id] = neiidx;
				neighborList[s].pop_back();
				AVLList[s] = deleteLast(AVLList[s]);
				AVLList[s] = updateNode(AVLList[s], neiidx, tmp->w);
			}
			adjList[s].erase(t);
			outSizeList[s]--;
			outWeightList[s] -= w;
			// if (AVLList[s]->leftCount + AVLList[s]->rightCount != outSizeList[s] || (AVLList[s]->leftSum + AVLList[s]->rightSum - outWeightList[s]) > 1e-6)
			// {
			// 	cout << "error" << endl;
			// 	exit(-1);
			// }
		}
		cout << filelabel << " new add update time: " << (chrono::duration_cast<chrono::nanoseconds>(chrono::high_resolution_clock::now() - begin).count() / 1000000000.0) / add_num << endl;
		// add
		begin = chrono::high_resolution_clock::now();
		for (int i = 0; i < add_num; i++)
		{
			uint s = sarr[i];
			uint t = tarr[i];
			double w = warr[i];
			neighborList[s].push_back(node{t, w});
			outSizeList[s]++;
			outWeightList[s] += w;
			AVLList[s] = addNode(AVLList[s], outSizeList[s] - 1, new AVLnode(warr[i]));
			adjList[s][t] = neighborList[s].size() - 1;
			// if (AVLList[s]->leftCount + AVLList[s]->rightCount != outSizeList[s] || (AVLList[s]->leftSum + AVLList[s]->rightSum - outWeightList[s]) > 1e-6)
			// {
			// 	cout << "error" << endl;
			// 	exit(-1);
			// }
		}
		cout << filelabel << " new del update time: " << (chrono::duration_cast<chrono::nanoseconds>(chrono::high_resolution_clock::now() - begin).count() / 1000000000.0) / add_num << endl;
	}

	void getAddEdge(long num)
	{
		cout << "No update file. Try MCAR first and then run MCPS" << endl;
		exit(-1);
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
		ifstream neiNumIn(neiNum.c_str(), ios::in | ios::binary);
		// ifstream neiWeightIn(neiWeight.c_str(), ios::in | ios::binary);
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
			AVLList[i] = NULL;
			for (uint j = 0; j < outSize; j++)
			{
				uint id;
				double w;
				neiNodeIn.read(reinterpret_cast<char *>(&id), sizeof(uint));
				// neiWeightIn.read(reinterpret_cast<char *>(&w), sizeof(double));
				w = 1;
				neighborList[i].push_back(node(id, w));
				outWeight += w;
				adjList[i][id] = j;
				AVLList[i] = addNode(AVLList[i], j, new AVLnode(w));
			}
			outSizeList[i] = outSize;
			outWeightList[i] = outWeight;
		}
		neiNumIn.close();
		// neiWeightIn.close();
		neiNodeIn.close();
	}

	AVLPrefixSumGraph() {}
	~AVLPrefixSumGraph() {}
};

class AliasMethodGraph
{
public:
	unordered_map<uint, Alias> aliasList;
	uint n;
	string filedir, filelabel;
	unordered_map<uint, vector<node>> neighborList;
	unordered_map<uint, uint> outSizeList;
	unordered_map<uint, unordered_map<uint, int>> adjList;
	AliasMethodGraph(const string &_filedir, const string &_filelabel)
	{
		filedir = _filedir;
		filelabel = _filelabel;
		ifstream opfile(filedir + "/" + filelabel + ".op");
		string neiNode, neiWeight, neiNum, graphAttr;
		neiNode = filedir + filelabel + ".outEdges";
		neiWeight = filedir + filelabel + ".outWEdges";
		neiNum = filedir + filelabel + ".outPtr";
		graphAttr = filedir + filelabel + ".attribute";
		cout << "FilePath: " << graphAttr.c_str() << endl;
		readFile(graphAttr, neiNode, neiWeight, neiNum);
		if (!opfile)
			getAddEdge(10000);
	}

	// void alias(double *p, uint *h, uint outsize, double sum, const vector<node> &neighbors)
	// {
	// 	double avg = outsize / sum;
	// 	uint *small = new uint[outsize];
	// 	uint *big = new uint[outsize];
	// 	int big_cnt = 0, small_cnt = 0;
	// 	for (int j = 0; j < outsize; j++)
	// 	{
	// 		p[j] = neighbors[j].w * avg;
	// 		if (p[j] > 1)
	// 			big[big_cnt++] = j;
	// 		else
	// 			small[small_cnt++] = j;
	// 	}
	// 	small_cnt--;
	// 	big_cnt--;
	// 	while (small_cnt >= 0 && big_cnt >= 0)
	// 	{
	// 		uint smallId = small[small_cnt];
	// 		uint bigId = big[big_cnt];
	// 		h[smallId] = bigId;
	// 		p[bigId] -= (1 - p[smallId]);
	// 		if (p[bigId] < 1)
	// 		{
	// 			small[small_cnt] = bigId;
	// 			big_cnt--;
	// 		}
	// 		else
	// 			small_cnt--;
	// 	}
	// 	delete[] small;
	// 	delete[] big;
	// }

	void update()
	{
		ifstream opfile(filedir + "/" + filelabel + ".op", ios::in);
		uint sarr[add_num];
		uint tarr[add_num];
		double warr[add_num];
		int opnum = 0;
		while (opfile >> sarr[opnum] >> tarr[opnum] >> warr[opnum])
		{
			opnum++;
		}
		opfile.close();
		// delete
		auto begin = chrono::high_resolution_clock::now();
		for (int i = 0; i < add_num; i++)
		{
			uint s = sarr[i];
			uint t = tarr[i];
			if (adjList[s].find(t) == adjList[s].end())
			{
				cout << "delete a nonexistent neighbor" << endl;
				exit(-2);
			}
			uint neiidx = adjList[s][t];
			uint outsize = outSizeList[s];
			if (neiidx == outsize - 1)
			{
				neighborList[s].pop_back();
			}
			else
			{
				auto tmp = neighborList[s].end() - 1;
				neighborList[s][neiidx].id = tmp->id;
				neighborList[s][neiidx].w = tmp->w;
				adjList[s][tmp->id] = neiidx;
				neighborList[s].pop_back();
			}
			adjList[s].erase(t);
			outSizeList[s]--;
			uint outSize = outSizeList[s];
			pair<int, double> *pi = new pair<int, double>[outSize];
			for (uint j = 0; j < outSize; j++)
			{
				pi[j] = make_pair(neighborList[s][j].id, neighborList[s][j].w);
			}
			aliasList[s] = Alias(pi, outSize);
			delete[] pi;
		}
		cout << filelabel << " new del update time: " << (chrono::duration_cast<chrono::nanoseconds>(chrono::high_resolution_clock::now() - begin).count() / 1000000000.0) / add_num << endl;
		// add
		begin = chrono::high_resolution_clock::now();
		for (int i = 0; i < add_num; i++)
		{
			uint s = sarr[i];
			uint t = tarr[i];
			neighborList[s].push_back(node{tarr[i], warr[i]});
			outSizeList[s]++;
			uint outSize = outSizeList[s];
			// delete[] aliasList[s].first;
			// delete[] aliasList[s].second;
			pair<int, double> *pi = new pair<int, double>[outSize];
			for (uint j = 0; j < outSize; j++)
			{
				pi[j] = make_pair(neighborList[s][j].id, neighborList[s][j].w);
			}
			aliasList[s] = Alias(pi, outSize);
			delete[] pi;
		}
		cout << filelabel << " new add update time: " << (chrono::duration_cast<chrono::nanoseconds>(chrono::high_resolution_clock::now() - begin).count() / 1000000000.0) / add_num << endl;
	}

	void getAddEdge(long num)
	{
		cout << "No update file. Try MCAR first and then run MCAM" << endl;
		exit(-1);
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
		// ifstream neiWeightIn(neiWeight.c_str(), ios::in | ios::binary);
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
			pair<int, double> *pi = new pair<int, double>[outSize];
			for (uint j = 0; j < outSize; j++)
			{
				uint id;
				double w;
				neiNodeIn.read(reinterpret_cast<char *>(&id), sizeof(uint));
				// neiWeightIn.read(reinterpret_cast<char *>(&w), sizeof(double));
				w = 1;
				neighborList[i].push_back(node(id, w));
				outWeight += w;
				pi[j] = make_pair(id, w);
				adjList[i][id] = j;
			}
			aliasList[i] = Alias(pi, outSize);
			delete[] pi;
			// uint *h;
			// double *p;
			// p = new double[outSize];
			// h = new uint[outSize];
			// alias(p, h, outSize, outWeight, neighborList[i]);
			// aliasList[i] = make_pair(h, p);
			outSizeList[i] = outSize;
		}
		neiNumIn.close();
		// neiWeightIn.close();
		neiNodeIn.close();
	}
	AliasMethodGraph() {}
	~AliasMethodGraph() {}
};

// struct smallSetID
// {
// 	int id;
// 	smallSetID *next, *pre;
// };
class subsetGraph
{
public:
	uint n;
	string filedir, filelabel;
	unordered_map<uint, unordered_map<int, vector<node>>> neighborList;
	unordered_map<uint, uint> outSizeList;
	unordered_map<uint, double> outWeightList;
	unordered_map<uint, unordered_map<uint, pair<int, int>>> adjList;
	// unordered_map<uint, unordered_map<int, smallSetID *>> smallSetMap;
	subsetGraph() {}
	subsetGraph(const string &_filedir, const string &_filelabel)
	{
		filedir = _filedir;
		filelabel = _filelabel;
		ifstream opfile(filedir + "/" + filelabel + ".op");
		string neiNode, neiWeight, neiNum, graphAttr;
		neiNode = filedir + filelabel + ".outEdges";
		neiWeight = filedir + filelabel + ".outWEdges";
		neiNum = filedir + filelabel + ".outPtr";
		graphAttr = filedir + filelabel + ".attribute";
		cout << "FilePath: " << graphAttr.c_str() << endl;
		readFile(graphAttr, neiNode, neiWeight, neiNum);
		if (!opfile)
			getAddEdge(10000);
	}

	~subsetGraph()
	{
	}

	void getAddEdge(long num)
	{
		cout << "No update file. Try MCAR first and then run MCSS" << endl;
		exit(-1);
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
		ifstream neiNumIn(neiNum.c_str(), ios::in | ios::binary);
		// ifstream neiWeightIn(neiWeight.c_str(), ios::in | ios::binary);
		ifstream neiNodeIn(neiNode.c_str(), ios::in | ios::binary);
		uint outSize;
		uint outSizeSum = 0, preOutSizeSum = 0;
		neiNumIn.read(reinterpret_cast<char *>(&outSizeSum), sizeof(uint));
		for (uint i = 0; i < n; i++)
		{
			// smallSetHead[i] = new smallSetID();
			// smallSetHead[i]->next = NULL;
			// smallSetHead[i]->pre = NULL;
			// smallSetHead[i]->id = -1;
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
				// neiWeightIn.read(reinterpret_cast<char *>(&w), sizeof(double));
				w = 1;
				neighbor.push_back(node(id, w));
				outWeight += w;
			}
			outSizeList[i] = outSize;
			outWeightList[i] = outWeight;
			// int gap = floor(log2(outWeight / (outSize * outSize))) + 1;
			for (uint j = 0; j < outSize; j++)
			{
				int subset = pow(2, ceil(log2(neighbor[j].w)));
				// if (subset < gap && smallSetMap[i].find(subset) == smallSetMap[i].end())
				// {
				// 	smallSetID *newSmallSet = new smallSetID();
				// 	newSmallSet->id = subset;
				// 	newSmallSet->pre = smallSetHead[i];
				// 	newSmallSet->next = smallSetHead[i]->next;
				// 	smallSetHead[i]->next = newSmallSet;
				// 	smallSetMap[i][subset] = newSmallSet;
				// }
				neighborList[i][subset].push_back(neighbor[j]);
				adjList[i][neighbor[j].id] = make_pair(subset, neighborList[i][subset].size() - 1);
			}
		}
		neiNumIn.close();
		// neiWeightIn.close();
		neiNodeIn.close();
	}

	// void doGap(int gap, int oldGap, uint s)
	// {
	// 	if (gap == oldGap)
	// 		return;
	// 	if (gap > oldGap)
	// 	{
	// 		if (oldGap < 1)
	// 			oldGap = 1;
	// 		while (oldGap < gap)
	// 		{
	// 			smallSetID *newSmallSet = new smallSetID();
	// 			newSmallSet->id = oldGap;
	// 			newSmallSet->pre = smallSetHead[s];
	// 			newSmallSet->next = smallSetHead[s]->next;
	// 			smallSetHead[s]->next = newSmallSet;
	// 			smallSetMap[s][oldGap] = newSmallSet;
	// 			oldGap++;
	// 		}
	// 	}
	// 	else
	// 	{
	// 		oldGap--;
	// 		while (oldGap >= gap && oldGap >= 1)
	// 		{
	// 			smallSetID *oldSmallSet = smallSetMap[s][oldGap];
	// 			if (oldSmallSet->pre != NULL)
	// 				oldSmallSet->pre->next = oldSmallSet->next;
	// 			if (oldSmallSet->next != NULL)
	// 				oldSmallSet->next->pre = oldSmallSet->pre;
	// 			delete oldSmallSet;
	// 			smallSetMap[s].erase(oldGap);
	// 			oldGap--;
	// 		}
	// 	}
	// 	return;
	// }

	void update()
	{
		ifstream opfile(filedir + "/" + filelabel + ".op", ios::in);
		uint sarr[add_num];
		uint tarr[add_num];
		double warr[add_num];
		int opnum = 0;
		while (opfile >> sarr[opnum] >> tarr[opnum] >> warr[opnum])
		{
			opnum++;
		}
		opfile.close();
		// delete
		auto begin = chrono::high_resolution_clock::now();
		for (int i = 0; i < add_num; i++)
		{
			uint s = sarr[i];
			uint t = tarr[i];
			if (adjList[s].find(t) == adjList[s].end())
			{
				cout << "delete a nonexistent neighbor" << endl;
				exit(-2);
			}
			int subsetID = adjList[s][t].first;
			int idx = adjList[s][t].second;
			uint subsetSize = neighborList[s][subsetID].size();
			// double oldGap = floor(log2(outWeightList[s] / (outSizeList[s] * outSizeList[s]))) + 1;
			// cout << outSizeList[s] << endl;
			// cout << outWeightList[s] << endl;
			if (idx == subsetSize - 1)
				neighborList[s][subsetID].pop_back();
			else
			{
				int tmpsize = neighborList[s][subsetID].size();
				auto tmp = neighborList[s][subsetID].end() - 1;
				neighborList[s][subsetID][idx].id = tmp->id;
				outWeightList[s] -= neighborList[s][subsetID][idx].w;
				neighborList[s][subsetID][idx].w = tmp->w;
				adjList[s][tmp->id].second = idx;
				neighborList[s][subsetID].pop_back();
			}
			adjList[s].erase(t);
			outSizeList[s]--;
			// cout << outSizeList[s] << endl;
			// cout << outWeightList[s] << endl;
			// double gap = floor(log2(outWeightList[s] / (outSizeList[s] * outSizeList[s]))) + 1;
			// if (gap > oldGap && gap > 1)
			// 	doGap(gap, oldGap, s);
			// if (gap < oldGap && oldGap > 1)
			// 	doGap(gap, oldGap, s);
		}
		cout << filelabel << " new del update time: " << (chrono::duration_cast<chrono::nanoseconds>(chrono::high_resolution_clock::now() - begin).count() / 1000000000.0) / add_num << endl;
		// add
		begin = chrono::high_resolution_clock::now();
		for (int i = 0; i < add_num; i++)
		{
			uint s = sarr[i];
			uint t = tarr[i];
			int subset = pow(2, floor(log2(warr[i])) + 1);
			// double oldGap = floor(log2(outWeightList[s] / (outSizeList[s] * outSizeList[s]))) + 1;
			neighborList[s][subset].push_back(node(t, warr[i]));
			outSizeList[s]++;
			outWeightList[s] += warr[i];
			adjList[s][t] = make_pair(subset, neighborList[s][subset].size() - 1);
			// double gap = floor(log2(outWeightList[s] / (outSizeList[s] * outSizeList[s]))) + 1;
			// if (gap > oldGap && gap > 1)
			// 	doGap(gap, oldGap, s);
			// if (gap < oldGap && oldGap > 1)
			// 	doGap(gap, oldGap, s);
		}
		cout << filelabel << " new add update time: " << (chrono::duration_cast<chrono::nanoseconds>(chrono::high_resolution_clock::now() - begin).count() / 1000000000.0) / add_num << endl;
	}
};

class BITPrefixSumGraph
{
public:
	unordered_map<uint, BIT> BITList;
	uint n;
	string filedir, filelabel;
	unordered_map<uint, vector<node>> neighborList;
	unordered_map<uint, uint> outSizeList;
	unordered_map<uint, double> outWeightList;
	BITPrefixSumGraph(const string &_filedir, const string &_filelabel)
	{
		filedir = _filedir;
		filelabel = _filelabel;
		string neiNode, neiWeight, neiNum, graphAttr;
		graphAttr = filedir + filelabel + ".initattribute_distr";
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
			graphAttr = filedir + filelabel + ".initattribute_distr";
			ofstream graphAttrOut(graphAttr.c_str());
			graphAttrOut << "n " << n;
			graphAttrOut.close();
			getAddEdge(10000);
		}
		else
		{
			graphAttrIn.close();
			neiNode = filedir + filelabel + ".initoutEdges_distr";
			neiWeight = filedir + filelabel + ".initoutWEdges_distr";
			neiNum = filedir + filelabel + ".initoutPtr_distr";
			readFile(graphAttr, neiNode, neiWeight, neiNum);
		}
	}

	void update()
	{
		ifstream opfile(filedir + "/" + filelabel + ".op", ios::in);
		uint sarr[add_num];
		uint tarr[add_num];
		double warr[add_num];
		int opnum = 0;
		while (opfile >> sarr[opnum] >> tarr[opnum] >> warr[opnum])
		{
			opnum++;
		}
		opfile.close();
		auto begin = chrono::high_resolution_clock::now();
		for (int i = 0; i < add_num; i++)
		{
			neighborList[sarr[i]].push_back(node{tarr[i], warr[i]});
			outSizeList[sarr[i]]++;
			outWeightList[sarr[i]] += warr[i];
			BITList[sarr[i]].add(warr[i]);
		}
		cout << filelabel << " new add update time: " << (chrono::duration_cast<chrono::nanoseconds>(chrono::high_resolution_clock::now() - begin).count() / 1000000000.0) / add_num << endl;
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
			double *arr = new double[outSize];

			for (uint j = 0; j < outSize; j++)
			{
				uint id;
				double w;
				neiNodeIn.read(reinterpret_cast<char *>(&id), sizeof(uint));
				neiWeightIn.read(reinterpret_cast<char *>(&w), sizeof(double));
				neighborList[i].push_back(node(id, w));
				outWeight += w;
				arr[j] = neighborList[i][j].w;
			}
			BITList[i] = BIT(outSize, arr);
			outSizeList[i] = outSize;
			outWeightList[i] = outWeight;
		}
		neiNumIn.close();
		neiWeightIn.close();
		neiNodeIn.close();
	}
	void getAddEdge(long num)
	{
		cout << "No update file. Try MCAR first and then run MCPS_BIT" << endl;
		exit(-1);
	}
	BITPrefixSumGraph() {}
	~BITPrefixSumGraph() {}
};