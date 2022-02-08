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
#include "alias.h"
#include "Graphother.h"
#include <unordered_map>
#include <unordered_set>
#include <sys/time.h>
#include <sys/stat.h>
#include <sys/resource.h>
#include <unistd.h>
#include <errno.h>
#include <queue>
#include <cmath>
#include <ctime>
#include <set>
#include <queue>
#include <math.h>
#include <map>
#include <chrono>

typedef unsigned int uint;

class SimStruct
{
public:
	uint L;
	double eps;
	double avg_time;
	uint final_count;
	uint *final_exist;
	double *final_p;
	uint *final_node;
	uint candidate_count[2];
	uint *cs_exist[2];
	uint *candidate_set[2];
	double *prob[2];
	SimStruct(string filedir, string filelabel, uint step)
	{
		L = step;
	}
	~SimStruct()
	{
	}

	void setEps(double epsilon)
	{
		eps = epsilon;
		avg_time = 0;
	}

	virtual void query(uint u) {}
	virtual void update() {}
};

class powermethod : public SimStruct
{
public:
	Graph g;
	uint candidate_count[2];
	uint *cs_exist[2];
	uint *candidate_set[2];
	double *prob[2];
	powermethod(string filedir, string filelabel, uint step) : SimStruct(filedir, filelabel, step)
	{
		g = Graph(filedir, filelabel);
		uint vert = g.n;
		final_p = new double[vert];
		final_node = new uint[vert];
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
			final_p[i] = 0;
			final_node[i] = 0;
			cs_exist[0][i] = 0;
			cs_exist[1][i] = 0;
			candidate_set[0][i] = 0;
			candidate_set[1][i] = 0;
			prob[0][i] = 0;
			prob[1][i] = 0;
		}
	}
	~powermethod()
	{
	}

	void query(uint u)
	{
		uint levelID = L % 2;
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
			// cout<<"Iteration "<<tempLevel<<": candidateCnt="<<candidateCnt<<endl;
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
				uint outSize = g.outSizeList[tempNode];
				double outVertWt = g.outWeightList[tempNode];
				double incre = tempP / outVertWt;

				for (uint k = 0; k < outSize; k++)
				{
					uint newNode = g.neighborList[tempNode][k].id;
					prob[newLevelID][newNode] += incre * g.neighborList[tempNode][k].w;
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
		for (uint j = 0; j < final_count; j++)
		{
			uint node = candidate_set[levelID][j];
			final_p[node] = prob[levelID][node];
			final_node[j] = node;
			prob[levelID][node] = 0;
			cs_exist[levelID][node] = 0;
		}
	}
	void update()
	{
		ifstream opfile(g.filedir + "/" + g.filelabel + ".op", ios::in);
		double tottime = 0;
		int opnum = 0;
		uint s, t;
		double w;
		while (opfile >> s >> t >> w)
		{
			clock_t t0 = clock();
			g.add(s, t, w);
			clock_t t1 = clock();
			tottime += (t1 - t0) / (double)CLOCKS_PER_SEC;
			opnum++;
		}
		cout << g.filelabel << " avg update time: " << tottime / opnum << endl;
		opfile.close();
	}
};

class MCPS : public SimStruct
{
public:
	Graph g;

	MCPS(string filedir, string filelabel, uint step) : SimStruct(filedir, filelabel, step)
	{
		g = Graph(filedir, filelabel);
		uint vert = g.n;
		final_p = new double[vert];
		final_node = new uint[vert];
		final_exist = new uint[vert];
		final_count = 0;
		for (uint i = 0; i < vert; i++)
		{
			final_p[i] = 0;
			final_node[i] = 0;
			final_exist[i] = 0;
		}
	}
	~MCPS()
	{
		delete[] final_p;
		delete[] final_node;
		delete[] final_exist;
	}

	void query(uint u)
	{
		for (uint j = 0; j < final_count; j++)
		{
			final_p[final_node[j]] = 0;
			final_exist[final_node[j]] = 0;
		}
		final_count = 0;
		unsigned long long w = 1 / eps / 0.25;
		cout << "w=" << w << endl;
		for (unsigned long long ni = 0; ni < w; ni++)
		{
			uint tempj = random_walk_prefixsum(u, L);
			final_p[tempj] += 1.0 / w;
			if (final_exist[tempj] == 0)
			{
				final_exist[tempj] = 1;
				final_node[final_count++] = tempj;
			}
		}
	}
	uint random_walk_prefixsum(uint u, uint len)
	{
		if (g.outSizeList[u] == 0)
			return u;
		uint curt;
		curt = u;
		uint i = 0;
		while (i++ < len)
		{
			uint outSize = g.outSizeList[curt];
			if (outSize == 0)
				break;
			double prefixsum = 0;
			double r = g.R.drand() * g.outWeightList[curt];
			for (uint j = 0; j < outSize; j++)
			{
				prefixsum += g.neighborList[curt][j].w;
				if (r < prefixsum)
				{
					curt = g.neighborList[curt][j].id;
					break;
				}
			}
		}
		return curt;
	}
	void update()
	{
		ifstream opfile(g.filedir + "/" + g.filelabel + ".op", ios::in);
		double tottime = 0;
		int opnum = 0;
		uint s, t;
		double w;
		while (opfile >> s >> t >> w)
		{
			clock_t t0 = clock();
			g.add(s, t, w);
			clock_t t1 = clock();
			tottime += (t1 - t0) / (double)CLOCKS_PER_SEC;
			opnum++;
		}
		cout << g.filelabel << " avg update time: " << tottime / opnum << endl;
		opfile.close();
	}
};

class MCPS_BIT : public SimStruct
{
public:
	BITPrefixSumGraph g;

	MCPS_BIT(string filedir, string filelabel, uint step) : SimStruct(filedir, filelabel, step)
	{
		g = BITPrefixSumGraph(filedir, filelabel);
		uint vert = g.n;
		final_p = new double[vert];
		final_node = new uint[vert];
		final_exist = new uint[vert];
		final_count = 0;
		for (uint i = 0; i < vert; i++)
		{
			final_p[i] = 0;
			final_node[i] = 0;
			final_exist[i] = 0;
		}
	}
	~MCPS_BIT()
	{
		delete[] final_p;
		delete[] final_node;
		delete[] final_exist;
	}

	void query(uint u)
	{
		for (uint j = 0; j < final_count; j++)
		{
			final_p[final_node[j]] = 0;
			final_exist[final_node[j]] = 0;
		}
		final_count = 0;
		unsigned long long w = 1 / eps / 0.25;
		cout << "w=" << w << endl;
		for (unsigned long long ni = 0; ni < w; ni++)
		{
			uint tempj = random_walk_prefixsum(u, L);
			final_p[tempj] += 1.0 / w;
			if (final_exist[tempj] == 0)
			{
				final_exist[tempj] = 1;
				final_node[final_count++] = tempj;
			}
		}
	}
	uint random_walk_prefixsum(uint u, uint len)
	{
		if (g.outSizeList[u] == 0)
			return u;
		uint curt;
		curt = u;
		uint i = 0;
		while (i++ < len)
		{
			uint outSize = g.outSizeList[curt];
			if (outSize == 0)
				break;
			double prefixsum = 0;
			double r = g.R.drand() * g.outWeightList[curt];
			for (uint j = 0; j < outSize; j++)
			{
				prefixsum += g.neighborList[curt][j].w;
				if (r < prefixsum)
				{
					curt = g.neighborList[curt][j].id;
					break;
				}
			}
		}
		return curt;
	}
	void update()
	{
		ifstream opfile(g.filedir + "/" + g.filelabel + ".op", ios::in);
		double tottime = 0;
		int opnum = 0;
		uint s, t;
		double w;
		while (opfile >> s >> t >> w)
		{
			clock_t t0 = clock();
			g.add(s, t, w);
			clock_t t1 = clock();
			tottime += (t1 - t0) / (double)CLOCKS_PER_SEC;
			opnum++;
		}
		cout << g.filelabel << " avg update time: " << tottime / opnum << endl;
		opfile.close();
	}
};

class MCAM : public SimStruct
{
public:
	AliasMethodGraph g;
	MCAM(string filedir, string filelabel, uint step) : SimStruct(filedir, filelabel, step)
	{
		g = AliasMethodGraph(filedir, filelabel);
		uint vert = g.n;
		final_p = new double[vert];
		final_node = new uint[vert];
		final_exist = new uint[vert];
		final_count = 0;
		for (uint i = 0; i < vert; i++)
		{
			final_p[i] = 0;
			final_node[i] = 0;
			final_exist[i] = 0;
		}
	}
	~MCAM()
	{
		delete[] final_p;
		delete[] final_node;
		delete[] final_exist;
	}

	void query(uint u)
	{
		for (uint j = 0; j < final_count; j++)
		{
			final_p[final_node[j]] = 0;
			final_exist[final_node[j]] = 0;
		}
		final_count = 0;
		unsigned long long w = 1 / eps / 0.25;
		cout << "w=" << w << endl;
		for (unsigned long long ni = 0; ni < w; ni++)
		{
			uint tempj = random_walk_alias(u, L);
			final_p[tempj] += 1.0 / w;
			if (final_exist[tempj] == 0)
			{
				final_exist[tempj] = 1;
				final_node[final_count++] = tempj;
			}
		}
	}
	uint random_walk_alias(uint u, uint len)
	{

		if (g.outSizeList[u] == 0)
			return u;
		uint curt;
		curt = u;
		uint i = 0;
		while (i++ < len)
		{
			if (g.outSizeList[u] == 0)
				break;
			curt = g.aliasList[curt].generateRandom(g.R);
		}
		return curt;
	}
	void update()
	{
		ifstream opfile(g.filedir + "/" + g.filelabel + ".op", ios::in);
		double tottime = 0;
		int opnum = 0;
		uint s, t;
		double w;
		while (opfile >> s >> t >> w)
		{
			clock_t t0 = clock();
			g.add(s, t, w);
			clock_t t1 = clock();
			tottime += (t1 - t0) / (double)CLOCKS_PER_SEC;
			opnum++;
		}
		cout << g.filelabel << " avg update time: " << tottime / opnum << endl;
		opfile.close();
		cout << "update result." << endl;
	}
};

class MCAR : public SimStruct
{
public:
	Graph g;

	MCAR(string filedir, string filelabel, uint step) : SimStruct(filedir, filelabel, step)
	{
		g = Graph(filedir, filelabel);
		uint vert = g.n;
		final_p = new double[vert];
		final_node = new uint[vert];
		final_exist = new uint[vert];
		final_count = 0;
		for (uint i = 0; i < vert; i++)
		{
			final_p[i] = 0;
			final_node[i] = 0;
			final_exist[i] = 0;
		}
	}
	~MCAR()
	{
		delete[] final_p;
		delete[] final_node;
		delete[] final_exist;
	}

	void query(uint u)
	{
		for (uint j = 0; j < final_count; j++)
		{
			final_p[final_node[j]] = 0;
			final_exist[final_node[j]] = 0;
		}
		final_count = 0;
		unsigned long long w = 1 / eps / 0.25;
		cout << "w=" << w << endl;
		for (unsigned long long ni = 0; ni < w; ni++)
		{
			uint tempj = random_walk_accept_reject(u, L);
			final_p[tempj] += 1.0 / w;
			if (final_exist[tempj] == 0)
			{
				final_exist[tempj] = 1;
				final_node[final_count++] = tempj;
			}
		}
	}
	uint random_walk_accept_reject(uint u, uint len)
	{
		if (g.outSizeList[u] == 0)
			return u;
		uint curt;
		curt = u;
		uint i = 0;
		while (i++ < len)
		{
			uint outSize = g.outSizeList[curt];
			if (outSize == 0)
				break;
			while (true)
			{
				double j = ceil(g.R.drand() * outSize);
				double r = g.R.drand();
				if (r < g.neighborList[curt][j].w / g.outWeightList[curt])
				{
					curt = g.neighborList[curt][j].id;
					break;
				}
			}
		}
		return curt;
	}
	void update()
	{
		ifstream opfile(g.filedir + "/" + g.filelabel + ".op", ios::in);
		double tottime = 0;
		int opnum = 0;
		uint s, t;
		double w;
		while (opfile >> s >> t >> w)
		{
			clock_t t0 = clock();
			g.add(s, t, w);
			clock_t t1 = clock();
			tottime += (t1 - t0) / (double)CLOCKS_PER_SEC;
			opnum++;
		}
		cout << g.filelabel << " avg update time: " << tottime / opnum << endl;
		opfile.close();
	}
};

class MCSS : public SimStruct
{
public:
	subsetGraph g;
	MCSS(string filedir, string filelabel, uint step) : SimStruct(filedir, filelabel, step)
	{
		g = subsetGraph(filedir, filelabel);
		uint vert = g.n;
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

		final_p = new double[vert];
		final_node = new uint[vert];
		final_exist = new uint[vert];
		final_count = 0;

		for (uint i = 0; i < vert; i++)
		{
			cs_exist[0][i] = 0;
			cs_exist[1][i] = 0;
			candidate_set[0][i] = 0;
			candidate_set[1][i] = 0;
			prob[0][i] = 0;
			prob[1][i] = 0;
			final_p[i] = 0;
			final_node[i] = 0;
			final_exist[i] = 0;
		}
	}
	~MCSS()
	{
		delete[] cs_exist[0];
		delete[] cs_exist[1];
		delete[] candidate_set[0];
		delete[] candidate_set[1];
		delete[] prob[0];
		delete[] prob[1];
		delete[] final_p;
		delete[] final_node;
		delete[] final_exist;
	}

	void query(uint u)
	{
		uint levelID = L % 2;
		for (uint j = 0; j < final_count; j++)
		{
			final_p[final_node[j]] = 0;
			final_exist[final_node[j]] = 0;
		}
		final_count = 0;
		uint nr = 0.1 * L / eps * 4;
		cout << "samples=" << nr << endl;
		uint tempLevel = 0;
		candidate_set[0][0] = u;
		candidate_count[0] = 1;
		candidate_count[1] = 0;
		prob[0][u] = 1.0;
		while (tempLevel <= L)
		{

			uint tempLevelID = tempLevel % 2;
			uint newLevelID = (tempLevel + 1) % 2;
			uint candidateCnt = candidate_count[tempLevelID];
			if (candidateCnt == 0)
				break;
			candidate_count[tempLevelID] = 0;
			for (uint j = 0; j < candidateCnt; j++)
			{
				uint tempNode = candidate_set[tempLevelID][j];
				double tempP = prob[tempLevelID][tempNode];
				cs_exist[tempLevelID][tempNode] = 0;
				prob[tempLevelID][tempNode] = 0;

				if (tempLevel == L)
				{
					if (final_exist[tempNode] == 0)
					{
						final_exist[tempNode] = 1;
						final_node[final_count++] = tempNode;
					}
					// cout<<"tempLevel="<<tempLevel<<" tempP="<<tempP<<" tempNode="<<tempNode<<endl;
					final_p[tempNode] = tempP;
					continue;
				}
				uint outSize = g.outSizeList[tempNode];
				double outVertWt = g.outWeightList[tempNode];
				double incre = tempP / outVertWt;
				for (auto setIt = g.neighborList[tempNode].begin(); setIt != g.neighborList[tempNode].end(); setIt++)
				{
					int setID = setIt->first;
					double increMax = incre * (pow(2, setID));

					if (increMax >= 1.0)
					{
						for (auto nodeIt = setIt->second.begin(); nodeIt != setIt->second.end(); nodeIt++)
						{
							uint newNode = nodeIt->id;
							prob[newLevelID][newNode] += incre * nodeIt->w;
							if (cs_exist[newLevelID][newNode] == 0)
							{
								cs_exist[newLevelID][newNode] = 1;
								candidate_set[newLevelID][candidate_count[newLevelID]++] = newNode;
							}
						}
					}
					else
					{
						double pmax = pow(2, setID) / outVertWt; // the maximum sampling probability in this subset;
						uint subsetSize = setIt->second.size();
						int seed = chrono::system_clock::now().time_since_epoch().count();
						std::default_random_engine generator(seed);
						std::geometric_distribution<int> distribution(pmax);

						for (uint k = 0; k < nr; k++)
						{
							int rgeo = distribution(generator);
							int curti = rgeo;
							while (curti < subsetSize)
							{
								node tmp = setIt->second[curti];
								double r = g.R.drand();
								if (r < tmp.w / pow(2, setID))
								{
									uint nodeidx = tmp.id;
									prob[newLevelID][nodeidx] += tempP;
									if (cs_exist[newLevelID][nodeidx] == 0)
									{
										cs_exist[newLevelID][nodeidx] = 1;
										candidate_set[newLevelID][candidate_count[newLevelID]++] = nodeidx;
									}
								}
								curti += 1;
								rgeo = distribution(generator);
								curti += rgeo;
							}
						}
					}
				}
			}
			tempLevel++;
		}
	}

	void update()
	{
		ifstream opfile(g.filedir + "/" + g.filelabel + ".op", ios::in);
		double tottime = 0;
		int opnum = 0;
		uint s, t;
		double w;
		while (opfile >> s >> t >> w)
		{
			clock_t t0 = clock();
			g.add(s, t, w);
			clock_t t1 = clock();
			tottime += (t1 - t0) / (double)CLOCKS_PER_SEC;
			opnum++;
		}
		cout << g.filelabel << " avg update time: " << tottime / opnum << endl;
		opfile.close();
	}
};

#endif
