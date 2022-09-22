#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <cstring>
#include <cmath>
#include <sys/time.h>
#include <unordered_map>
#include <boost/math/distributions/binomial.hpp>
#include <boost/random.hpp>
#include "Graph.h"
using namespace std;
typedef unsigned int uint;

int main(int argc, char **argv)
{
    long querynum = 10;
    vector<double> epss;
    uint L = 10;
    string filedir, filelabel;
    int is_update = 0;
    argParser(argc, argv, filedir, filelabel, querynum, epss, L, is_update);
    subsetGraph g(filedir, filelabel);
    if (is_update)
    {
        g.update();
        double pkm = peak_mem() / 1024.0 / 1024.0;
        cout << "Total graph: peak memory: " << pkm << " G" << endl;
        exit(0);
    }
    string queryname;
    queryname = "./query/" + filelabel + ".query";
    ifstream query;
    cout << "Input query file from: " << queryname << endl;
    double *final_p = new double[g.n];
    uint *final_node = new uint[g.n];
    uint *final_exist = new uint[g.n];
    uint final_count = 0;
    uint *cs_exist[2];
    cs_exist[0] = new uint[g.n];
    cs_exist[1] = new uint[g.n];
    //当前层candidate_set的点
    uint *candidate_set[2];
    candidate_set[0] = new uint[g.n];
    candidate_set[1] = new uint[g.n];
    uint candidate_count[2];
    candidate_count[0] = 0;
    candidate_count[1] = 0;
    //当前层该点的probability
    double *prob[2];
    prob[0] = new double[g.n];
    prob[1] = new double[g.n];
    for (uint i = 0; i < g.n; i++)
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
    int seed = chrono::system_clock::now().time_since_epoch().count();
    boost::mt19937 rng;
    stringstream ss_run;
    ss_run << "./analysis/CoinFlipWalk_" << filelabel << "_runtime.csv";
    ofstream writecsv;
    writecsv.open(ss_run.str(), ios::app);
    for (auto epsIt = epss.begin(); epsIt != epss.end(); epsIt++)
    {
        double eps = *epsIt;
        query.open(queryname);
        double avg_time = 0;
        for (uint i = 0; i < querynum; i++)
        {
            Random R;
            uint u;
            query >> u;
            cout << i << ": " << u << endl;
            clock_t t0 = clock();
            uint levelID = L % 2;
            final_count = 0;
            uint nr = 0.1 * L / eps * 4;
            cout << "samples=" << nr << endl;
            for (uint k = 0; k < nr; k++)
            {
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
                            final_p[tempNode] += tempP / (double)nr;
                            continue;
                        }
                        uint outSize = g.outSizeList[tempNode];
                        double outVertWt = g.outWeightList[tempNode];
                        double incre = tempP / outVertWt;
                        for (auto setIt = g.neighborList[tempNode].begin(); setIt != g.neighborList[tempNode].end(); setIt++)
                        {
                            int setID = setIt->first;
                            double increMax = incre * setID;
                            uint subsetSize = setIt->second.size();
                            if (increMax >= 1)
                            {
                                for (uint setidx = 0; setidx < subsetSize; setidx++)
                                {
                                    node &tmpnode = setIt->second[setidx];
                                    uint newNode = tmpnode.id;
                                    prob[newLevelID][newNode] += incre * tmpnode.w;
                                    if (cs_exist[newLevelID][newNode] == 0)
                                    {
                                        cs_exist[newLevelID][newNode] = 1;
                                        candidate_set[newLevelID][candidate_count[newLevelID]++] = newNode;
                                    }
                                }
                            }
                            else
                            {
                                boost::binomial_distribution<> bio(subsetSize, increMax);
                                int rbio = bio(rng);
                                for (uint j = 0; j < rbio; j++)
                                {
                                    double r1 = floor(R.drand() * subsetSize);
                                    node tmp = setIt->second[r1];
                                    double r2 = R.drand();
                                    if (r2 < tmp.w / setID)
                                    {
                                        prob[newLevelID][tmp.id] += 1;
                                        if (cs_exist[newLevelID][tmp.id] == 0)
                                        {
                                            cs_exist[newLevelID][tmp.id] = 1;
                                            candidate_set[newLevelID][candidate_count[newLevelID]++] = tmp.id;
                                        }
                                    }
                                }
                            }
                        }
                    }
                    tempLevel++;
                }
            }

            clock_t t1 = clock();
            avg_time += (t1 - t0) / (double)CLOCKS_PER_SEC;
            cout << "Query time for node " << u << ": " << (t1 - t0) / (double)CLOCKS_PER_SEC << " s";
            stringstream ss_dir, ss;
            ss_dir << "./result/CoinFlipWalk/" << filelabel << "/" << L << "/" << eps << "/";
            ss << ss_dir.str() << u << ".txt";
            cout << "Write query results in file: " << ss.str() << endl;
            mkpath(ss_dir.str());
            ofstream fout;
            fout.open(ss.str());
            fout.setf(ios::fixed, ios::floatfield);
            fout.precision(15);
            if (!fout)
                cout << "Fail to open the writed file" << endl;
            cout << "final_count=" << final_count << endl;
            for (uint j = 0; j < final_count; j++)
            {
                fout << final_node[j] << " " << final_p[final_node[j]] << endl;
                final_p[final_node[j]] = 0;
                final_exist[final_node[j]] = 0;
            }
            fout.close();
        }
        query.close();
        cout << endl;
        cout << "query time: " << avg_time / (double)querynum << " s" << endl;
        cout << "==== "
             << "CoinFlipWalk"
             << " with " << eps << " on " << filelabel << " done!====" << endl;
        writecsv << avg_time / (double)querynum << ',';
    }

    delete[] cs_exist[0];
    delete[] cs_exist[1];
    delete[] candidate_set[0];
    delete[] candidate_set[1];
    delete[] prob[0];
    delete[] prob[1];
    delete[] final_p;
    delete[] final_node;
    delete[] final_exist;
    cout << endl
         << endl
         << endl;

    writecsv.close();
}
