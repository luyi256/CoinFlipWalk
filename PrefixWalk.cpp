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

int main(int argc, char** argv)
{
    long i = 1;
    long querynum = 10;
    vector<double> epss;
    uint L = 10;
    string filedir, filelabel;
    int is_update = 0;
    argParser(argc, argv, filedir, filelabel, querynum, epss, L, is_update);
    AVLPrefixSumGraph g(filedir, filelabel);
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
    double* final_p = new double[g.n];
    uint* final_node = new uint[g.n];
    uint* final_exist = new uint[g.n];
    uint final_count = 0;
    for (uint i = 0; i < g.n; i++)
    {
        final_p[i] = 0;
        final_node[i] = 0;
        final_exist[i] = 0;
    }
    stringstream ss_run;
    ss_run << "./analysis/PrefixWalk_" << filelabel << "_" << L << "_runtime.csv";

    for (auto epsIt = epss.begin(); epsIt != epss.end(); epsIt++)
    {
        ofstream writecsv;
        writecsv.open(ss_run.str(), ios::app);
        double eps = *epsIt;
        query.open(queryname);
        double avg_time = 0;
        for (uint i = 0; i < querynum; i++)
        {
            uint nodeId;
            query >> nodeId;
            cout << i << ": " << nodeId << endl;
            clock_t t0 = clock();
            for (uint j = 0; j < final_count; j++)
            {
                final_p[final_node[j]] = 0;
                final_exist[final_node[j]] = 0;
            }
            final_count = 0;
            unsigned long long w = 1 / eps / 0.25;
            cout << "w=" << w << endl;
            Random R;
            for (unsigned long long ni = 0; ni < w; ni++)
            {
                uint u = nodeId;
                if (g.outSizeList[u] == 0)
                    return u;
                uint i = 0;
                double r;
                int nodeno;
                while (i++ < L)
                {
                    if (g.outSizeList[u] == 0)
                        break;
                    r = R.drand() * g.outWeightList[u];
                    uint tmpsize = g.outSizeList[u];
                    nodeno = prefixSumIndex(g.AVLList[u], r, 0) - 1;
                    if (nodeno >= tmpsize) {
                        cout << "error, out of outsize" << endl;
                        exit(-1);
                    }
                    u = g.neighborList[u][nodeno].id;
                }
                if (g.outSizeList[u] == 0 && i <= L)
                    continue;
                final_p[u] += 1.0 / w;
                if (final_exist[u] == 0)
                {
                    final_exist[u] = 1;
                    final_node[final_count++] = u;
                }
            }
            clock_t t1 = clock();
            avg_time += (t1 - t0) / (double)CLOCKS_PER_SEC;
            // ofstream memout("mem_PrefixWalk_" + filelabel + ".txt", ios::app);
            double pkm = peak_mem() / 1024.0 / 1024.0;
            cout << "Total process: peak memory: " << pkm << " G" << endl;
            // double pkrss = peak_rss() / 1024.0 / 1024.0;
            // memout << ", peak rss: " << pkrss << " G" << endl;
            cout << "Query time for node " << nodeId << ": " << (t1 - t0) / (double)CLOCKS_PER_SEC << " s";
            stringstream ss_dir, ss;
            ss_dir << "./result/PrefixWalk/" << filelabel << "/" << L << "/" << eps << "/";
            ss << ss_dir.str() << nodeId << ".txt";
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
            }
            fout.close();
        }
        query.close();
        cout << endl;
        cout << "query time: " << avg_time / (double)querynum << " s" << endl;
        cout << "==== "
            << "PrefixWalk"
            << " with " << eps << " on " << filelabel << " done!====" << endl;
        writecsv << avg_time / (double)querynum << ',';
        writecsv.close();
    }
    delete[] final_p;
    delete[] final_node;
    delete[] final_exist;
    cout << endl
        << endl
        << endl;
}
