#include <algorithm>
#include <iostream>
#include <string>
#include <sstream>
#include <cstdlib>
#include <vector>
#include "SimStruct.h"
#include <fstream>
#include <cstring>
#include <unordered_set>
#include <sys/time.h>
#include <sys/resource.h>
#include <unistd.h>
#include <errno.h>
//#include "alias.h"

void usage()
{
	cerr << "./edgepush" << endl
		 << "-d <directory>" << endl
		 << "-f <filelabel>" << endl
		 << "-algo <algorithm>" << endl
		 << "[-e <epsilon> (default 1e-05)]" << endl
		 << "[-a <teleport probability> (default 0.2)]" << endl
		 << "[-qn <querynum> (default 10)]" << endl;
}

typedef unsigned int uint;

long check_inc(long i, long max)
{
	if (i == max)
	{
		usage();
		exit(1);
	}
	return i + 1;
}

bool cmp(const pair<uint, double> &a, const pair<uint, double> &b)
{
	return a.second > b.second;
}

Graph g;

void gen_affinity(string filelabel, uint k)
{
	srand(time(NULL));
	uint vert = 10000;
	uint **pos = new uint *[vert];

	Random Ran = Random();
	for (uint u = 0; u < vert; u++)
	{
		pos[u] = new uint[k];
		for (uint v = 0; v < k; v++)
		{
			uint tmppos = Ran.generateRandom() % 10000;
			pos[u][v] = tmppos;
		}
	}

	stringstream ss_dir, ss_out;
	ss_dir << "./dataset/";
	ss_out << ss_dir.str() << filelabel << ".txt";
	cout << ss_out.str() << endl;

	mkpath(ss_dir.str());
	ofstream fout;
	fout.open(ss_out.str());
	if (!fout)
	{
		cout << "ERROR to open affinity graphs" << endl;
		return;
	}
	fout << vert << "\n";
	for (uint i = 0; i < vert; i++)
	{
		for (uint j = 0; j < vert; j++)
		{
			if (i == j)
			{
				continue;
			}
			double dist = 0;
			for (uint t = 0; t < k; t++)
			{
				dist += ((pos[i][t] - pos[j][t]) * (pos[i][t] - pos[j][t]) / (double)vert / (double)k / (double)k);
			}
			dist = exp((-1) * dist);
			dist = (double)vert * dist;
			if (dist < 1e-12)
			{
				dist = 1e-12;
			}
			fout << i << " " << j << " " << dist << "\n";
		}
	}

	for (uint i = 0; i < vert; i++)
	{
		delete[] pos[i];
	}
	fout.close();
}

int main(int argc, char **argv)
{
	long i = 1;
	char *endptr;
	string filedir = "./dataset/youtube/";
	string filelabel = "youtube";
	string algo = "EdgePush";
	long querynum = 10;
	double eps = 0.1;
	uint L = 10;
	while (i < argc)
	{
		if (!strcmp(argv[i], "-d"))
		{
			i = check_inc(i, argc);
			filedir = argv[i];
		}
		else if (!strcmp(argv[i], "-f"))
		{
			i = check_inc(i, argc);
			filelabel = argv[i];
		}
		else if (!strcmp(argv[i], "-algo"))
		{
			i = check_inc(i, argc);
			algo = argv[i];
		}
		else if (!strcmp(argv[i], "-e"))
		{
			i = check_inc(i, argc);
			eps = strtod(argv[i], &endptr);
			if ((eps == 0 || eps > 1) && endptr)
			{
				cerr << "Invalid eps argument" << endl;
				exit(1);
			}
		}
		else if (!strcmp(argv[i], "-qn"))
		{
			i = check_inc(i, argc);
			querynum = strtod(argv[i], &endptr);
			if ((querynum < -2) && endptr)
			{
				cerr << "Invalid querynum argument" << endl;
				exit(1);
			}
		}
		else if (!strcmp(argv[i], "-l"))
		{
			i = check_inc(i, argc);
			int tempL = int(strtod(argv[i], &endptr));
			if ((tempL < 0) && endptr)
			{
				cerr << "Invalid hop number" << endl;
				exit(1);
			}
			L = uint(tempL);
		}
		else
		{
			usage();
			exit(1);
		}
		i++;
	}

	cout << "=========" << endl;
	cout << "Dataset: " << filelabel << endl;
	cout << "Algorithm: " << algo << endl;
	cout << "L=" << L << endl;
	cout << "eps=" << eps << endl;
	cout << endl;

	if (filelabel == "AG_1" || filelabel == "AG_2" || filelabel == "AG_3" || filelabel == "AG_4")
	{
		ifstream f_ag;
		string ag_name;
		ag_name = "./dataset/" + filelabel + ".txt";
		f_ag.open(ag_name);
		if (!f_ag)
		{
			cout << "Generate affinity graphs..." << endl;
			uint tmpd = 1;
			if (filelabel == "AG_1")
			{
				tmpd = 1;
			}
			else if (filelabel == "AG_2")
			{
				tmpd = 2;
			}
			else if (filelabel == "AG_3")
			{
				tmpd = 3;
			}
			else if (filelabel == "AG_4")
			{
				tmpd = 4;
			}
			gen_affinity(filelabel, tmpd);
		}
		f_ag.close();
	}

	SimStruct *sim = NULL;

	if (algo == "powermethod")
		sim = new powermethod(filedir, filelabel, eps, L);
	else if (algo == "edgepushByLevel")
		sim = new edgepushByLevel(filedir, filelabel, eps, L);
	else if (algo == "localPush")
		sim = new localPush(filedir, filelabel, eps, L);
	else if (algo == "MCSS") // subset sampling
		sim = new MCSS(filedir, filelabel, eps, L);
	else if (algo == "MCPS") // prefix sum
		sim = new MCPS(filedir, filelabel, eps, L);
	else if (algo == "MCAR") // accept reject
		sim = new MCAR(filedir, filelabel, eps, L);
	else if (algo == "MCAM") // alias method
		sim = new MCAM(filedir, filelabel, eps, L);
	g = sim->g;

	if (querynum > sim->vert)
	{
		querynum = sim->vert;
	}
	cout << endl;
	cout << "querynum=" << querynum << endl;
	string queryname;
	queryname = "./query/" + filelabel + ".query";
	ifstream query;
	query.open(queryname);
	if (!query)
	{
		cout << "Generate query file..." << endl;
		pair<int, double> *aliasD = new pair<int, double>[sim->vert];
		uint *check = new uint[sim->vert]();
		for (int i = 0; i < sim->vert; i++)
		{
			aliasD[i] = make_pair(i, ((sim->g).getOutVertWeight(i) / (double)(sim->g).totaldeg));
		}
		Alias alias = Alias();
		alias.input(aliasD, sim->vert);
		ofstream data_idx("./query/" + filelabel + ".query");
		for (int i = 0; i < querynum; i++)
		{
			int tmpnode = alias.generateRandom(sim->R);
			while (((sim->g).getOutSize(tmpnode) == 0) || check[tmpnode] == 1)
			{
				tmpnode = alias.generateRandom(sim->R);
			}
			check[tmpnode] = 1;
			data_idx << tmpnode << "\n";
		}
		data_idx.close();
		query.open(queryname);
		if (!query)
		{
			cout << "ERROR:input query file:" << queryname << endl;
			return 0;
		}
	}
	exit(-1);
	cout << "Input query file from: " << queryname << endl;

	// stringstream ss_sample;
	// ss_sample << "./analysis/" << algo << "_" << filelabel << "_sample.csv";
	// ofstream writesample;
	// writesample.open(ss_sample.str(), ios::app);
	// double sample;
	// if (algo == "MCSS")
	// {
	// 	sample = L / eps * 4;
	// 	writesample << sample << ',';
	// }
	// else
	// {
	// 	sample = 1 / eps * 4;
	// 	writesample << sample << endl;
	// }
	// double samplerate = 0;
	for (uint i = 0; i < querynum; i++)
	{
		uint nodeId;
		query >> nodeId;
		cout << i << ": " << nodeId << endl;

		clock_t t0 = clock();
		sim->query(nodeId);
		clock_t t1 = clock();
		sim->avg_time += (t1 - t0) / (double)CLOCKS_PER_SEC;
		cout << "Query time for node " << nodeId << ": " << (t1 - t0) / (double)CLOCKS_PER_SEC << " s" << endl;

		stringstream ss_dir, ss;
		if (algo == "powermethod")
		{
			ss_dir << "./result/" << algo << "/" << filelabel << "/" << L << "/";
		}
		else
		{
			ss_dir << "./result/" << algo << "/" << filelabel << "/" << L << "/" << eps << "/";
		}
		mkpath(ss_dir.str());
		if (algo == "powermethod")
		{
			ss << ss_dir.str() << nodeId << "_gt.txt";
		}
		else
		{
			ss << ss_dir.str() << nodeId << ".txt";
		}
		cout << "Write query results in file: " << ss.str() << endl;
		ofstream fout;
		fout.open(ss.str());
		fout.setf(ios::fixed, ios::floatfield);
		fout.precision(15);
		if (!fout)
			cout << "Fail to open the writed file" << endl;
		if (algo == "localPush" || algo == "edgepushByLevel" || algo == "powermethod")
		{
			uint levelID = L % 2;
			for (uint j = 0; j < sim->final_count; j++)
			{
				uint node = sim->candidate_set[levelID][j];
				fout << node << " " << sim->prob[levelID][node] << endl;
			}
		}
		else
		{
			for (uint j = 0; j < sim->final_count; j++)
			{
				fout << sim->final_node[j] << " " << sim->final_p[sim->final_node[j]] << endl;
			}
		}
		fout.close();
		// samplerate += sim->actual_sample / sample;
	}
	// if (algo == "MCSS")
	// 	writesample << samplerate / querynum << endl;
	// writesample.close();
	cout << endl;
	cout << "query time: " << sim->avg_time / (double)querynum << " s" << endl;
	cout << "==== " << algo << " on " << filelabel << " done!====" << endl;
	cout << endl
		 << endl
		 << endl;
	query.close();
	if (algo == "powermethod")
		return 0;
	stringstream ss_run;
	ss_run << "./analysis/" << algo << "_" << filelabel << "_runtime.csv";
	ofstream writecsv;
	writecsv.open(ss_run.str(), ios::app);
	writecsv << sim->avg_time / (double)querynum << ',';
	writecsv.close();
	// if (algo != "edgepushByLevel")
	// 	return 0;
	// stringstream ss_push;
	// ss_push << "./pushrate/" << filelabel << "_" << algo << "_pushrate.csv";
	// writecsv.open(ss_push.str(), ios::app);
	// writecsv << "pushNum:" << sim->pushNum << " tot2push:" << sim->tot2push << " pushRate:" << double(sim->pushNum) / (double)sim->tot2push << ',' << endl;
	// writecsv.close();
	return 0;
}
