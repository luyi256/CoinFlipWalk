
#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <cstring>
#include <cmath>
#include <sys/time.h>
#include <unordered_map>
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
    Graph g(filedir, filelabel);
    if (is_update)
    {
        g.update();
        double pkm = peak_mem() / 1024.0 / 1024.0;
        cout << "Total graph: peak memory: " << pkm << " G" << endl;
        exit(0);
    }
    double pkm = peak_mem() / 1024.0 / 1024.0;
    cout << "Total graph: peak memory: " << pkm << " G" << endl;
    string queryname;
    queryname = "./query/" + filelabel + ".query";
    ifstream query;
    cout << "Input query file from: " << queryname << endl;
    uint final_count = 0;
    uint *cs_exist[2];
    cs_exist[0] = new uint[g.n];
    cs_exist[1] = new uint[g.n];
    uint *candidate_set[2];
    candidate_set[0] = new uint[g.n];
    candidate_set[1] = new uint[g.n];
    uint candidate_count[2];
    candidate_count[0] = 0;
    candidate_count[1] = 0;
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
    }
    stringstream ss_run;
    ss_run << "./analysis/powermethod_" << filelabel << "_runtime.csv";
    ofstream writecsv;
    writecsv.open(ss_run.str(), ios::app);
    query.open(queryname);
    if (!query)
    {
        query.close();
        g.getQuery(queryname);
        query.open(queryname);
    }
    double avg_time = 0;
    for (uint i = 0; i < querynum; i++)
    {
        uint nodeId;
        query >> nodeId;
        cout << i << ": " << nodeId << endl;
        clock_t t0 = clock();
        final_count = 0;
        uint tempLevel = 0;
        candidate_set[0][0] = nodeId;
        candidate_count[0] = 1;
        candidate_count[1] = 0;
        prob[0][nodeId] = 1.0;
        clock_t t1 = clock();
        int levelID = L % 2;
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
                uint outSize = g.outSizeList[tempNode];
                double outVertWt = 0;
                for (auto iter = g.neighborList[tempNode].begin(); iter != g.neighborList[tempNode].end(); iter++)
                {
                    outVertWt += iter->w;
                }
                double incre = tempP / outVertWt;
                for (auto iter = g.neighborList[tempNode].begin(); iter != g.neighborList[tempNode].end(); iter++)
                {
                    auto newNode = iter->id;
                    auto newWeight = iter->w;

                    prob[newLevelID][newNode] += incre * newWeight;
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
        avg_time += (t1 - t0) / (double)CLOCKS_PER_SEC;
        // ofstream memout("mem_powermethod_" + filelabel + ".txt", ios::app);
        double pkm = peak_mem() / 1024.0 / 1024.0;
        cout << "Total process: peak memory: " << pkm << " G" << endl;
        // double pkrss = peak_rss() / 1024.0 / 1024.0;
        // memout << ", peak rss: " << pkrss << " G" << endl;
        cout << "Query time for node " << nodeId << ": " << (t1 - t0) / (double)CLOCKS_PER_SEC << " s";
        stringstream ss_dir, ss;
        ss_dir << "./result/powermethod/" << filelabel << "/" << L << "/";
        ss << ss_dir.str() << nodeId << "_gt.txt";
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
            uint nodeidx = candidate_set[levelID][j];
            fout << nodeidx << " " << prob[levelID][nodeidx] << endl;
            prob[levelID][nodeidx] = 0;
            cs_exist[levelID][nodeidx] = 0;
        }
        fout.close();
    }
    query.close();
    cout << endl;
    cout << "query time: " << avg_time / (double)querynum << " s" << endl;
    cout << "==== "
         << "powermethod"
         << " on " << filelabel << " done!====" << endl;
    writecsv << avg_time / (double)querynum << ',';

    delete[] cs_exist[0];
    delete[] cs_exist[1];
    delete[] candidate_set[0];
    delete[] candidate_set[1];
    delete[] prob[0];
    delete[] prob[1];
    cout << endl
         << endl
         << endl;

    writecsv.close();
}
