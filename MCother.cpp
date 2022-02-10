#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <cstring>
#include <cmath>
#include <sys/time.h>
#include <unordered_map>
#include "SimStruct.h"
using namespace std;
typedef unsigned int uint;
int mkpath(string s, mode_t mode = 0755)
{
    size_t pre = 0, pos;
    string dir;
    int mdret;
    if (s[s.size() - 1] != '/')
    {
        s += '/';
    }
    while ((pos = s.find_first_of('/', pre)) != string::npos)
    {
        dir = s.substr(0, pos++);
        pre = pos;
        if (dir.size() == 0)
            continue;
        if ((mdret = ::mkdir(dir.c_str(), mode)) && errno != EEXIST)
        {
            return mdret;
        }
    }
    return mdret;
}
int main(int argc, char **argv)
{
    long i = 1;
    string filedir, filelabel, algoname;
    char *endptr;
    long querynum = 10;
    // long querynum=1;
    vector<double> epss;
    // uint L;//hanzhi
    uint L = 10; // hanzhi
    while (i < argc)
    {
        if (!strcmp(argv[i], "-d"))
        {
            if (++i < argc)
                filedir = argv[i];
        }
        else if (!strcmp(argv[i], "-f"))
        {
            if (++i < argc)
                filelabel = argv[i];
        }
        else if (!strcmp(argv[i], "-algo"))
        {
            if (++i < argc)
                algoname = argv[i];
        }
        else if (!strcmp(argv[i], "-e"))
        {
            while (++i < argc && argv[i][0] != '-')
            {
                double eps = strtod(argv[i], &endptr);
                if ((eps == 0 || eps > 1) && endptr)
                {
                    cerr << "Invalid eps argument" << endl;
                    exit(1);
                }
                epss.push_back(eps);
            }
            i--;
        }
        else if (!strcmp(argv[i], "-qn"))
        {
            if (++i < argc)
                querynum = strtod(argv[i], &endptr);
            if ((querynum < -2) && endptr)
            {
                cerr << "Invalid querynum argument:" << querynum << endl;
                exit(1);
            }
        }
        else if (!strcmp(argv[i], "-l"))
        {
            int tempL;
            if (++i < argc)
                tempL = int(strtod(argv[i], &endptr));
            if ((tempL < 0) && endptr)
            {
                cerr << "Invalid hop number" << endl;
                exit(1);
            }
            L = uint(tempL);
        }
        else
        {
            cerr << "Err args:" << argv[i] << endl;
            exit(1);
        }
        i++;
    }
    cout << "=========" << endl;
    cout << "Dataset: " << filelabel << endl;
    cout << "Algorithm: " << algoname << endl;
    cout << "L=" << L << endl;
    cout << endl;
    SimStruct *sim = NULL;
    if (algoname == "powermethod")
        sim = new powermethod(filedir, filelabel, L);
    else if (algoname == "MCPS") // prefix sum
        sim = new MCPS(filedir, filelabel, L);
    else if (algoname == "MCPS_BIT")
        sim = new MCPS_BIT(filedir, filelabel, L);
    else if (algoname == "MCAR") // accept reject
        sim = new MCAR(filedir, filelabel, L);
    else if (algoname == "MCAM") // alias method
        sim = new MCAM(filedir, filelabel, L);
    else if (algoname == "MCSS")
        sim = new MCSS(filedir, filelabel, L);
    cout << "Graph init done." << endl;
    sim->update();
    cout << endl;
    cout << "querynum=" << querynum << endl;
    string queryname;
    queryname = "./query/" + filelabel + ".query";
    ifstream query;
    cout << "Input query file from: " << queryname << endl;
    stringstream ss_run;
    ss_run << "./analysis/" << algoname << "_" << filelabel << "_runtime.csv";
    ofstream writecsv;
    writecsv.open(ss_run.str(), ios::app);
    for (auto epsIt = epss.begin(); epsIt != epss.end(); epsIt++)
    {
        double eps = *epsIt;
        sim->setEps(eps);
        query.open(queryname);
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
            if (algoname == "powermethod")
            {
                ss_dir << "./result/" << algoname << "/" << filelabel << "/" << L << "/";
            }
            else
            {
                ss_dir << "./result/" << algoname << "/" << filelabel << "/" << L << "/" << eps << "/";
            }
            if (algoname == "powermethod")
            {
                ss << ss_dir.str() << nodeId << "_gt.txt";
            }
            else
            {
                ss << ss_dir.str() << nodeId << ".txt";
            }
            cout << "Write query results in file: " << ss.str() << endl;
            mkpath(ss_dir.str());
            ofstream fout;
            fout.open(ss.str());
            fout.setf(ios::fixed, ios::floatfield);
            fout.precision(15);
            if (!fout)
                cout << "Fail to open the writed file" << endl;

            // if (false)
            // {
            //     uint levelID = L % 2;
            //     for (uint j = 0; j < sim->final_count; j++)
            //     {
            //         uint node = sim->candidate_set[levelID][j];
            //         fout << node << " " << sim->prob[levelID][node] << endl;
            //     }
            // }
            // else
            // {
            cout << "sim->final_count=" << sim->final_count << endl;
            for (uint j = 0; j < sim->final_count; j++)
            {
                fout << sim->final_node[j] << " " << sim->final_p[sim->final_node[j]] << endl;
            }
            // }
            fout.close();
            // samplerate += sim->actual_sample / sample;
        }
        query.close();
        cout << endl;
        cout << "query time: " << sim->avg_time / (double)querynum << " s" << endl;
        cout << "==== " << algoname << " with " << eps << " on " << filelabel << " done!====" << endl;
        if (algoname != "powermethod")
        {
            writecsv << sim->avg_time / (double)querynum << ',';
        }
    }
    cout << endl
         << endl
         << endl;

    writecsv.close();
}
