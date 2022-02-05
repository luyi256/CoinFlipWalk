#ifndef ALIAS_H
#define ALIAS_H

#include <random>
#include <algorithm>
#include <stack>
#include "Random.h"
using namespace std;

class Alias
{
public:
	double *p;
	int *h;
	int *map1;
	int n;
	Alias() {}
	void input(pair<int, double> *pi, uint size)
	{
		n = size;
		double sum = 0;
		// stack<int> small;
		// stack<int> big;
		//初始化代价比stack高
		int *small = new int[n];
		int *big = new int[n];
		int small_cnt = 0, big_cnt = 0;
		p = new double[n];
		h = new int[n];
		map1 = new int[n];
		for (int i = 0; i < n; i++)
		{
			sum += pi[i].second;
			map1[i] = pi[i].first;
		}
		for (int i = 0; i < n; i++)
		{
			p[i] = pi[i].second * n / sum;
			if (p[i] > 1)
				// big.push(i);
				big[big_cnt++] = i;
			else
				// small.push(i);
				small[small_cnt++] = i;
		}
		//cnt的意义变成指示
		small_cnt--;
		big_cnt--;
		while (small_cnt >= 0 && big_cnt >= 0)
		{
			int smallId = small[small_cnt];
			int bigId = big[big_cnt];
			h[smallId] = bigId;
			p[bigId] -= (1 - p[smallId]);
			if (p[bigId] < 1)
			{
				small[small_cnt] = bigId;
				big_cnt--;
			}
			else
				small_cnt--;
		}
		delete[] small;
		delete[] big;
	}
	Alias(pair<int, double> *pi, uint size)
	{
		n = size;
		double sum = 0;
		// stack<int> small;
		// stack<int> big;
		//初始化代价比stack高
		int *small = new int[n];
		int *big = new int[n];
		int small_cnt = 0, big_cnt = 0;
		p = new double[n];
		h = new int[n];
		map1 = new int[n];
		for (int i = 0; i < n; i++)
		{
			sum += pi[i].second;
			map1[i] = pi[i].first;
		}
		for (int i = 0; i < n; i++)
		{
			p[i] = pi[i].second * n / sum;
			if (p[i] > 1)
				// big.push(i);
				big[big_cnt++] = i;
			else
				// small.push(i);
				small[small_cnt++] = i;
		}
		//cnt的意义变成指示
		small_cnt--;
		big_cnt--;
		while (small_cnt >= 0 && big_cnt >= 0)
		{
			int smallId = small[small_cnt];
			int bigId = big[big_cnt];
			h[smallId] = bigId;
			p[bigId] -= (1 - p[smallId]);
			if (p[bigId] < 1)
			{
				small[small_cnt] = bigId;
				big_cnt--;
			}
			else
				small_cnt--;
		}
		delete[] small;
		delete[] big;
	}

	Alias(const Alias &tmp)
	{
		n = tmp.n;
		p = new double[n];
		h = new int[n];
		map1 = new int[n];
		memcpy(p, tmp.p, sizeof(double) * n);
		memcpy(h, tmp.h, sizeof(int) * n);
		memcpy(map1, tmp.map1, sizeof(int) * n);
	}

	Alias &operator=(const Alias &tmp)
	{
		if (&tmp != this)
		{
			n = tmp.n;
			p = new double[n];
			h = new int[n];
			map1 = new int[n];
			memcpy(p, tmp.p, sizeof(double) * n);
			memcpy(h, tmp.h, sizeof(int) * n);
			memcpy(map1, tmp.map1, sizeof(int) * n);
		}

		return *this;
	}

	~Alias()
	{
		delete[] p;
		delete[] h;
		delete[] map1;
	}
	int generateRandom(Random &R)
	{
		int firstId = R.drand() * n;
		int answer = R.drand() < p[firstId] ? map1[firstId] : map1[h[firstId]];
		return answer;
	}
};

#endif
