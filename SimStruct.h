#ifndef SIMSTRUCT_H
#define SIMSTRUCT_H

#define INT_MAX 32767

#include <vector>
#include <algorithm>
#include <queue>
#include <functional>
#include <iostream>
#include <thread>
#include <string>
#include <sstream>
#include <fstream>
//#include "Graph.h"
#include "Random.h"
#include "alias.h"
#include "priority.h"
#include <unordered_map>
#include <unordered_set>
#include <sys/time.h>
#include <sys/stat.h>
#include <sys/resource.h>
#include <unistd.h>
#include <errno.h>
#include <queue>
#include <cmath>
#include <random>
#include <ctime>
#include <set>
#include <queue>
#include <math.h>
#include <map>
#include <chrono>

typedef unsigned int uint;
enum sorted
{
	UNSORTED,
	SORTED
};
class SimStruct
{
public:
	Graph g;	//Class Graph
	Random R;	//Class Random
	uint vert;	//# of vertice
	uint nedge; //# of edges
	string filelabel;
	uint candidate_count[2];
	uint *cs_exist[2];
	uint *candidate_set[2];
	double *prob[2];
	uint final_count;
	uint *final_exist;
	double *final_p;
	uint *final_node;
	uint L;
	double eps;
	double avg_time;
	uint actual_sample;
	// uint tot2push, pushNum;
	double pushRate;
	SimStruct(string name, string file_label, double epsilon, uint step, sorted isSorted)
	{
		filelabel = file_label;
		g.inputGraph(name, file_label, epsilon, isSorted);
		R = Random();
		vert = g.n;
		nedge = g.m;
		//是否存在在当前层中的candidate_set中
		cs_exist[0] = new uint[vert];
		cs_exist[1] = new uint[vert];
		//当前层candidate_set的点
		candidate_set[0] = new uint[vert];
		candidate_set[1] = new uint[vert];
		candidate_count[0] = 0;
		candidate_count[1] = 0;
		//当前层该点的probability
		prob[0] = new double[vert];
		prob[1] = new double[vert];

		final_count = 0;
		for (uint i = 0; i < vert; i++)
		{
			cs_exist[0][i] = 0;
			cs_exist[1][i] = 0;
			candidate_set[0][i] = 0;
			candidate_set[1][i] = 0;
			prob[0][i] = 0;
			prob[1][i] = 0;
		}
		L = step;
		eps = epsilon;
	}

	~SimStruct()
	{
		delete[] cs_exist[0];
		delete[] cs_exist[1];
		delete[] candidate_set[0];
		delete[] candidate_set[1];
		delete[] prob[0];
		delete[] prob[1];
	}

	virtual void query(uint u) {}
	uint random_walk_alias(uint u, uint len)
	{

		if (g.getOutSize(u) == 0)
			return u;
		uint curt;
		curt = u;
		uint i = 0;
		while (i++ < len)
		{
			if (g.getOutSize(curt) == 0)
				break;
			curt = g.alias[curt].generateRandom(R);
		}
		return curt;
	}
	uint random_walk_accept_reject(uint u, uint len)
	{
		if (g.getOutSize(u) == 0)
			return u;
		uint curt;
		curt = u;
		uint i = 0;
		while (i++ < len)
		{
			uint outSize = g.getOutSize(curt);
			if (outSize == 0)
				break;
			while (true)
			{
				double j = ceil(R.drand() * outSize);
				double r = R.drand();
				if (r < g.getOutEdgeWeight(curt, j) / g.getOutVertWeight(curt))
				{
					curt = g.getOutVert(curt, j);
					break;
				}
			}
		}
		return curt;
	}
	uint random_walk_prefixsum(uint u, uint len)
	{
		if (g.getOutSize(u) == 0)
			return u;
		uint curt;
		curt = u;
		uint i = 0;
		while (i++ < len)
		{
			uint outSize = g.getOutSize(curt);
			if (outSize == 0)
				break;
			double prefixsum = 0;
			double r = R.drand() * g.getOutVertWeight(curt);
			for (uint j = 0; j < outSize; j++)
			{
				prefixsum += g.getOutEdgeWeight(curt, j);
				if (r < prefixsum)
				{
					curt = g.getOutVert(curt, j);
					break;
				}
			}
		}
		return curt;
	}
};

class powermethod : public SimStruct
{
public:
	powermethod(string name, string file_label, double epsilon, uint step) : SimStruct(name, file_label, epsilon, step, UNSORTED)
	{
	}
	~powermethod()
	{
	}

	void query(uint u)
	{
		uint levelID = L % 2;
		for (uint j = 0; j < final_count; j++)
		{
			uint node = candidate_set[levelID][j];
			prob[levelID][node] = 0;
			cs_exist[levelID][node] = 0;
		}
		final_count = 0;
		uint tempLevel = 0;
		candidate_set[0][0] = u;
		candidate_count[0] = 1;
		candidate_count[1] = 0;
		prob[tempLevel][u] = 1.0;

		while (tempLevel < L)
		{
			uint tempLevelID = tempLevel % 2;
			uint newLevelID = (tempLevel + 1) % 2;
			uint candidateCnt = candidate_count[tempLevelID];
			//cout<<"Iteration "<<tempLevel<<": candidateCnt="<<candidateCnt<<endl;
			if (candidateCnt == 0)
			{
				cout << "candidateCnt=0 tempLevel=" << tempLevel << endl;
				break;
			}
			candidate_count[tempLevelID] = 0;
			for (uint j = 0; j < candidateCnt; j++)
			{
				uint tempNode = candidate_set[tempLevelID][j];
				double tempP = prob[tempLevelID][tempNode];
				cs_exist[tempLevelID][tempNode] = 0;
				prob[tempLevelID][tempNode] = 0;
				uint outSize = g.getOutSize(tempNode);
				double outVertWt = g.getOutVertWeight(tempNode);
				double incre = tempP / outVertWt;

				for (uint k = 0; k < outSize; k++)
				{
					uint newNode = g.getOutVert(tempNode, k);
					prob[newLevelID][newNode] += incre * g.getOutEdgeWeight(tempNode, k);
					if (cs_exist[newLevelID][newNode] == 0)
					{
						cs_exist[newLevelID][newNode] = 1;
						candidate_set[newLevelID][candidate_count[newLevelID]++] = newNode;
					}
				}
			}
			tempLevel++;
		}
		final_count = candidate_count[levelID];
	}
};

class localPush : public SimStruct
{
public:
	localPush(string name, string file_label, double epsilon, uint step) : SimStruct(name, file_label, epsilon, step, UNSORTED)
	{
	}
	~localPush()
	{
	}

	void query(uint u)
	{
		uint levelID = L % 2;
		for (uint j = 0; j < final_count; j++)
		{
			uint node = candidate_set[levelID][j];
			prob[levelID][node] = 0;
			cs_exist[levelID][node] = 0;
		}
		final_count = 0;
		uint tempLevel = 0;
		candidate_set[0][0] = u;
		candidate_count[0] = 1;
		candidate_count[1] = 0;
		prob[tempLevel][u] = 1.0;

		while (tempLevel < L)
		{
			uint tempLevelID = tempLevel % 2;
			uint newLevelID = (tempLevel + 1) % 2;
			uint candidateCnt = candidate_count[tempLevelID];
			if (candidateCnt == 0)
			{
				cout << "candidateCnt=0 tempLevel=" << tempLevel << endl;
				break;
			}
			candidate_count[tempLevelID] = 0;
			for (uint j = 0; j < candidateCnt; j++)
			{
				uint tempNode = candidate_set[tempLevelID][j];
				double tempP = prob[tempLevelID][tempNode];
				cs_exist[tempLevelID][tempNode] = 0;
				prob[tempLevelID][tempNode] = 0;
				if (tempP < eps)
					continue;
				uint outSize = g.getOutSize(tempNode);
				double outVertWt = g.getOutVertWeight(tempNode);
				double incre = tempP / outVertWt;

				for (uint k = 0; k < outSize; k++)
				{
					uint newNode = g.getOutVert(tempNode, k);
					prob[newLevelID][newNode] += incre * g.getOutEdgeWeight(tempNode, k);
					if (cs_exist[newLevelID][newNode] == 0)
					{
						cs_exist[newLevelID][newNode] = 1;
						candidate_set[newLevelID][candidate_count[newLevelID]++] = newNode;
					}
				}
			}
			tempLevel++;
		}
		final_count = candidate_count[levelID];
	}
};

class MCPS : public SimStruct
{
public:
	MCPS(string name, string file_label, double epsilon, uint step) : SimStruct(name, file_label, epsilon, step, UNSORTED)
	{
		final_p = new double[vert];
		final_node = new uint[vert];
		// final_c = new uint[vert];
		for (uint i = 0; i < vert; i++)
		{
			final_p[i] = 0;
			final_node[i] = 0;
		}
	}
	~MCPS()
	{
		delete[] final_p;
		delete[] final_node;
		// delete[] final_c;
	}

	void query(uint u)
	{
		for (uint j = 0; j < final_count; j++)
		{
			final_p[final_node[j]] = 0;
			cs_exist[1][final_node[j]] = 0;
		}
		final_count = 0;
		unsigned long long w = 1 / eps / 0.25;
		cout << "w=" << w << endl;
		for (unsigned long long ni = 0; ni < w; ni++)
		{
			uint tempj = random_walk_prefixsum(u, L);
			final_p[tempj] += 1.0 / w;
			if (cs_exist[1][tempj] == 0)
			{
				cs_exist[1][tempj] = 1;
				final_node[final_count++] = tempj;
				// final_c[tempj] = 0;
			}
			// final_c[tempj] += 1;
			// if (ni % 100000 == 0)
			// 	cout << "the " << ni << "th..." << endl;
		}
	}
};

class MCAM : public SimStruct
{
public:
	MCAM(string name, string file_label, double epsilon, uint step) : SimStruct(name, file_label, epsilon, step, UNSORTED)
	{
		final_p = new double[vert];
		final_node = new uint[vert];
		// final_c = new uint[vert];
		for (uint i = 0; i < vert; i++)
		{
			final_p[i] = 0;
			final_node[i] = 0;
		}
	}
	~MCAM()
	{
		delete[] final_p;
		delete[] final_node;
		// delete[] final_c;
	}

	void query(uint u)
	{
		for (uint j = 0; j < final_count; j++)
		{
			final_p[final_node[j]] = 0;
			cs_exist[1][final_node[j]] = 0;
		}
		final_count = 0;
		unsigned long long w = 1 / eps / 0.25;
		cout << "w=" << w << endl;
		for (unsigned long long ni = 0; ni < w; ni++)
		{
			uint tempj = random_walk_alias(u, L);
			final_p[tempj] += 1.0 / w;
			if (cs_exist[1][tempj] == 0)
			{
				cs_exist[1][tempj] = 1;
				final_node[final_count++] = tempj;
				// final_c[tempj] = 0;
			}
			// final_c[tempj] += 1;
			// if (ni % 100000 == 0)
			// 	cout << "the " << ni << "th..." << endl;
		}
	}
};

class MCAR : public SimStruct
{
public:
	MCAR(string name, string file_label, double epsilon, uint step) : SimStruct(name, file_label, epsilon, step, UNSORTED)
	{
		final_p = new double[vert];
		final_node = new uint[vert];
		// final_c = new uint[vert];
		for (uint i = 0; i < vert; i++)
		{
			final_p[i] = 0;
			final_node[i] = 0;
		}
	}
	~MCAR()
	{
		delete[] final_p;
		delete[] final_node;
		// delete[] final_c;
	}

	void query(uint u)
	{
		for (uint j = 0; j < final_count; j++)
		{
			final_p[final_node[j]] = 0;
			cs_exist[1][final_node[j]] = 0;
		}
		final_count = 0;
		unsigned long long w = 1 / eps / 0.25;
		cout << "w=" << w << endl;
		for (unsigned long long ni = 0; ni < w; ni++)
		{
			uint tempj = random_walk_accept_reject(u, L);
			final_p[tempj] += 1.0 / w;
			if (cs_exist[1][tempj] == 0)
			{
				cs_exist[1][tempj] = 1;
				final_node[final_count++] = tempj;
				// final_c[tempj] = 0;
			}
			// final_c[tempj] += 1;
			// if (ni % 100000 == 0)
			// 	cout << "the " << ni << "th..." << endl;
		}
	}
};

class edgepushByLevel : public SimStruct
{
public:
	edgepushByLevel(string name, string file_label, double epsilon, uint step) : SimStruct(name, file_label, epsilon, step, UNSORTED)
	{
	}
	~edgepushByLevel()
	{
	}

	void query(uint u)
	{
		uint levelID = L % 2;
		for (uint j = 0; j < final_count; j++)
		{
			uint node = candidate_set[levelID][j];
			prob[levelID][node] = 0;
			cs_exist[levelID][node] = 0;
		}
		final_count = 0;
		if (g.getOutSize(u) == 0)
			return;
		double Gku = g.initsort[u][0] - 1.0;
		if (Gku > 0)
			return;
		uint tempLevel = 0;
		candidate_set[0][0] = u;
		candidate_count[0] = 1;
		candidate_count[1] = 0;
		prob[tempLevel][u] = 1.0;
		while (tempLevel < L)
		{
			uint tempLevelID = tempLevel % 2;
			uint newLevelID = (tempLevel + 1) % 2;
			uint candidateCnt = candidate_count[tempLevelID];
			if (candidateCnt == 0)
			{
				cout << "candidateCnt=0 tempLevel=" << tempLevel << endl;
				break;
			}
			candidate_count[tempLevelID] = 0;
			for (uint j = 0; j < candidateCnt; j++)
			{
				uint tempNode = candidate_set[tempLevelID][j];
				double tempP = prob[tempLevelID][tempNode];
				cs_exist[tempLevelID][tempNode] = 0;
				prob[tempLevelID][tempNode] = 0;
				uint outSize = g.getOutSize(tempNode);
				if (outSize == 0)
					continue;
				double Gku = g.initsort[tempNode][0] - tempP;
				if (Gku >= 0)
					continue;
				double outVertWt = g.getOutVertWeight(tempNode);
				double incre = tempP / outVertWt;
				uint idx = 0;
				while (Gku < 0)
				{
					uint newNode = g.getPriorOutVert(tempNode, idx);
					prob[newLevelID][newNode] += incre * g.getPriorOutEdgeWeight(tempNode, idx);
					if (cs_exist[newLevelID][newNode] == 0)
					{
						cs_exist[newLevelID][newNode] = 1;
						candidate_set[newLevelID][candidate_count[newLevelID]++] = newNode;
					}
					if (idx + 1 < outSize)
						Gku = g.initsort[tempNode][++idx] - tempP;
					else
						break;
				}
				// if (Gku < 0)
				// 	pushNum += outSize;
				// else
				// 	pushNum += idx;
				// tot2push += outSize;
			}
			tempLevel++;
		}
		final_count = candidate_count[levelID];
	}
};

class MCSS : public SimStruct
{
public:
	MCSS(string name, string file_label, double epsilon, uint step) : SimStruct(name, file_label, epsilon, step, SORTED)
	{
		final_p = new double[vert];
		final_node = new uint[vert];
		final_exist = new uint[vert];
		for (uint i = 0; i < vert; i++)
		{
			final_p[i] = 0;
			final_node[i] = 0;
			final_exist[i] = 0;
		}
	}
	~MCSS()
	{
		delete[] final_p;
		delete[] final_node;
		delete[] final_exist;
	}

	void query(uint u)
	{
		// actual_sample = 0;
		for (uint j = 0; j < final_count; j++)
		{
			final_p[final_node[j]] = 0;
			final_exist[final_node[j]] = 0;
		}
		final_count = 0;
		uint levelID = L % 2;
		uint nr = L / eps * 4;
		cout << "samples=" << nr << endl;
		for (uint i = 0; i < nr; i++)
		{
			uint tempLevel = 0;
			candidate_set[0][0] = u;
			candidate_count[0] = 1;
			candidate_count[1] = 0;
			prob[tempLevel][u] = 1.0;

			while (tempLevel < L)
			{
				uint tempLevelID = tempLevel % 2;
				uint newLevelID = (tempLevel + 1) % 2;
				uint candidateCnt = candidate_count[tempLevelID];
				// if (candidateCnt == 0)
				// {
				// 	cout << "candidateCnt=0 tempLevel=" << tempLevel << endl;
				// 	break;
				// }
				if (candidateCnt == 0)
					break;
				candidate_count[tempLevelID] = 0;
				for (uint j = 0; j < candidateCnt; j++)
				{
					uint tempNode = candidate_set[tempLevelID][j];
					double tempP = prob[tempLevelID][tempNode];
					cs_exist[tempLevelID][tempNode] = 0;
					prob[tempLevelID][tempNode] = 0;
					uint outSize = g.getOutSize(tempNode);
					double outVertWt = g.getOutVertWeight(tempNode);
					double incre = tempP / outVertWt;
					uint k;
					for (k = 0; k < outSize; k++)
					{
						uint newNode = g.getOutVert(tempNode, k);
						double edgeIncre = incre * g.getOutEdgeWeight(tempNode, k);
						if (edgeIncre >= 1)
						{
							prob[newLevelID][newNode] += edgeIncre;
							if (cs_exist[newLevelID][newNode] == 0)
							{
								cs_exist[newLevelID][newNode] = 1;
								candidate_set[newLevelID][candidate_count[newLevelID]++] = newNode;
							}
						}
						else
							break;
					}
					while (k < outSize)
					{
						double edgeWt = g.getOutEdgeWeight(tempNode, k);
						double edgeIncre = incre * edgeWt;
						int seed = chrono::system_clock::now().time_since_epoch().count();
						default_random_engine generator(seed);
						geometric_distribution<int> distribution(edgeIncre);
						int rgeo = distribution(generator);
						// int rgeo=floor(-1.0/edgeIncre*log(R.drand()));
						uint next = k + rgeo;
						if (next >= outSize)
							break;
						double nextEdgeWt = g.getOutEdgeWeight(tempNode, next);
						double r = R.drand();
						double comp = nextEdgeWt / edgeWt;
						if (r < nextEdgeWt / edgeWt)
						{
							uint newNode = g.getOutVert(tempNode, next);
							prob[newLevelID][newNode] += 1;
							if (cs_exist[newLevelID][newNode] == 0)
							{
								cs_exist[newLevelID][newNode] = 1;
								candidate_set[newLevelID][candidate_count[newLevelID]++] = newNode;
							}
						}
						k = next + 1;
					}
				}
				tempLevel++;
			}
			uint cnt = candidate_count[levelID];
			for (uint j = 0; j < cnt; j++)
			{
				uint node = candidate_set[levelID][j];
				final_p[node] += (1.0 / nr) * prob[levelID][node];
				if (final_exist[node] == 0)
				{
					final_exist[node] = 1;
					final_node[final_count++] = node;
				}
				prob[levelID][node] = 0;
				cs_exist[levelID][node] = 0;
			}
			// actual_sample += cnt;
		}
	}
};

#endif