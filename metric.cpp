#include <algorithm>
#include <iostream>
#include <string>
#include <sstream>
#include <cstdlib>
#include <vector>
#include <fstream>
#include <cstring>
#include <unordered_set>
#include "Graphother.h"
#include "Random.h"
#include <random>
using namespace std;
Graph g;
vector<uint> query_set;
class metric
{
private:
    uint vert;
    double *gtvalues;
    vector<pair<double, uint> > unnm_algoanswers;
    vector<pair<double, int> > unnm_gtanswers;
    uint k;
    double *algovalues;
    bool *algocheck;
    vector<uint> topk_algo_Nodes;

public:
    double cal_maxAE();
    double cal_l1Error();
    double calPrecision();
    double conductance();
    metric(string filedir, string filelabel, string algoname, long querynum, uint L, double eps);
    ~metric();
};

double metric::calPrecision()
{
    if (unnm_gtanswers.size() < k)
    {
        return 1;
    }
    uint hitCount = 0;
    uint real_gtsize = k;
    for (uint i = k; i < unnm_gtanswers.size(); i++)
    {
        if (unnm_gtanswers[i].first != unnm_gtanswers[k - 1].first)
        {
            break;
        }
        else
        {
            real_gtsize += 1;
        }
    }
    for (uint i = 0; i < k; i++)
    {
        for (uint j = 0; j < real_gtsize; j++)
        {
            if (topk_algo_Nodes[i] == unnm_gtanswers[j].second)
            {
                hitCount += 1;
                break;
            }
        }
    }
    return (hitCount / (double)k); //?
}

double metric::conductance()
{
    uint cvert = g.n;
    double ctotald = g.totdeg;

    double numcut = 0;
    double crt_vol = 0;
    double crt_Phi = 1.0;
    uint nonzero = unnm_algoanswers.size();
    uint *cutcheck = new uint[cvert]();

    for (uint i = 0; i < nonzero; i++)
    {
        if (unnm_algoanswers[i].first == 0)
        {
            break;
        }
        uint crt_node = unnm_algoanswers[i].second;
        uint crt_OutSize = g.outSizeList[crt_node];
        double crt_outdeg = g.outWeightList[crt_node];

        if (crt_OutSize == 0)
        {
            continue;
        }
        crt_vol += crt_outdeg;
        cutcheck[crt_node] = 1;
        for (uint j = 0; j < crt_OutSize; j++)
        {
            uint tmpOutV = g.neighborList[crt_node][j].id; // tmp out-neighbor node
            if (cutcheck[tmpOutV] == 1)
            {
                numcut -= g.neighborList[crt_node][j].w;
            }
            else
            {
                numcut += g.neighborList[crt_node][j].w;
            }
        }

        if (numcut <= 0)
        {
            numcut = 0;
            continue;
        }
        double min_vol = crt_vol;
        if (crt_vol > (ctotald - crt_vol))
        {
            min_vol = ctotald - crt_vol;
        }

        if (min_vol == 0.0)
        {
            continue;
        }

        double tmp_Phi = numcut / min_vol;
        if (tmp_Phi <= crt_Phi)
        {
            crt_Phi = tmp_Phi;
        }
    }

    delete[] cutcheck;
    return crt_Phi;
}

double metric::cal_l1Error()
{
    double l1Err = 0;
    uint *nnzarr = new uint[vert]();
    double tmp_err;
    for (uint j = 0; j < unnm_algoanswers.size(); j++)
    {
        uint tmpnode = unnm_algoanswers[j].second;
        tmp_err = abs(unnm_algoanswers[j].first - gtvalues[tmpnode]);
        l1Err += tmp_err;
        nnzarr[tmpnode] = 1;
    }
    for (uint j = 0; j < vert; j++)
    {
        if (nnzarr[j] == 0)
        {
            tmp_err = gtvalues[j];
            l1Err += tmp_err;
        }
    }
    delete[] nnzarr;
    return l1Err;
}

double metric::cal_maxAE()
{
    double err = 0;
    uint *nnzarr = new uint[vert]();
    double tmp_err, max_err = 0, algoans;
    uint max_node = -1;
    for (uint j = 0; j < unnm_algoanswers.size(); j++)
    {
        uint tmpnode = unnm_algoanswers[j].second;
        tmp_err = abs(unnm_algoanswers[j].first - gtvalues[tmpnode]);
        if (tmp_err > max_err)
        {
            max_err = tmp_err;
            max_node = tmpnode;
            algoans = unnm_algoanswers[j].first;
        }
        nnzarr[tmpnode] = 1;
    }
    for (uint j = 0; j < vert; j++)
    {
        if (nnzarr[j] == 0)
        {
            tmp_err = gtvalues[j];
            if (tmp_err > max_err)
            {
                max_err = tmp_err;
                max_node = j;
                algoans = 0;
            }
        }
    }
    delete[] nnzarr;
    return max_err;
}

metric::metric(string filedir, string filelabel, string algoname, long querynum, uint L, double eps)
{
    k = 50;
    double avg_l1_error = 0, avg_pre50 = 0, avg_conductance = 0, avg_max_additive_error = 0;
    vert = g.n;

    for (uint i = 0; i < querynum; ++i)
    {
        uint u = query_set[i];
        stringstream ss_gt, ss_gt_dir;
        ss_gt << "./result/powermethod/" << filelabel << "/" << L << "/" << u << "_gt.txt";
        ifstream gtin;
        gtin.open(ss_gt.str());
        // cout << "query " << i << ":" << u << endl;
        // cout << "powermethod path:" << ss_gt.str() << endl;
        if (!gtin)
        {
            cout << "ERROR:unable to open groundtruth file " << ss_gt.str() << endl;
            return;
        }

        stringstream ss_algo;
        ss_algo << "./result/" << algoname << "/" << filelabel << "/" << L << "/" << eps << "/" << u << ".txt";
        ifstream algoin(ss_algo.str());
        if (!algoin)
        {
            cout << "ERROR:unable to open result file " << ss_algo.str() << endl;
            return;
        }

        gtvalues = new double[vert]();

        uint gtCnt = 0;
        uint gt_tempNode;
        double gt_tempSim;
        while (gtin >> gt_tempNode >> gt_tempSim)
        {
            if (gt_tempSim > 0.0)
            {
                gtvalues[gt_tempNode] = gt_tempSim;
                unnm_gtanswers.push_back(make_pair(gt_tempSim, gt_tempNode));
                gtCnt++;
            }
        }

        sort(unnm_gtanswers.begin(), unnm_gtanswers.end(), greater<pair<double, uint> >());
        uint precision_num = 50;
        if (gtCnt < 50)
        {
            precision_num = gtCnt;
        }

        // cout << "algo path" << ss_algo.str() << endl;
        vector<uint> algoNodes;
        algovalues = new double[vert]();
        algocheck = new bool[vert]();
        uint realCnt = 0;
        uint algo_tempNode;
        double algo_tempSim;
        while (algoin >> algo_tempNode >> algo_tempSim)
        {
            algoNodes.push_back(algo_tempNode);
            if (algo_tempSim > 0.0)
            {
                algovalues[algo_tempNode] = algo_tempSim;
                unnm_algoanswers.push_back(make_pair(algo_tempSim, algo_tempNode));
                algocheck[algo_tempNode] = true;
                realCnt++;
            }
        }
        sort(unnm_algoanswers.begin(), unnm_algoanswers.end(), greater<pair<double, uint> >());
        uint topknum = precision_num;
        if ((uint)unnm_algoanswers.size() < topknum)
        {
            topknum = (uint)unnm_algoanswers.size();
        }

        for (uint x = 0; x < topknum; x++)
        {
            topk_algo_Nodes.push_back(unnm_algoanswers[x].second);
        }
        //补完剩下不足的topknum
        if (topknum < precision_num)
        {
            uint tmpCnt = topknum;
            for (uint supid = 0; supid < vert; supid++)
            {
                if (tmpCnt >= precision_num)
                {
                    break;
                }
                else if (algocheck[supid] == false)
                {
                    topk_algo_Nodes.push_back(supid);
                    algocheck[supid] = true;
                    tmpCnt += 1;
                }
            }
        }
        avg_conductance += conductance();
        avg_l1_error += cal_l1Error();
        avg_pre50 += calPrecision();
        // printf("Ans for query no.%u\n", i);
        avg_max_additive_error += cal_maxAE();
        delete[] algovalues;
        delete[] gtvalues;
        delete[] algocheck;
        unnm_gtanswers.clear();
        unnm_algoanswers.clear();
        topk_algo_Nodes.clear();
    }
    avg_l1_error /= (double)querynum;
    avg_pre50 /= (double)querynum;
    avg_conductance /= (double)querynum;
    avg_max_additive_error /= (double)querynum;
    cout << "eps = " << eps << endl;
    cout << "Avg l1-error = " << avg_l1_error << endl;
    cout << "Avg normalized precision@50 = " << avg_pre50 << endl;
    cout << "Avg conductance = " << avg_conductance << endl;
    cout << "Avg max-additive-error = " << avg_max_additive_error << endl;
    stringstream ss_run;
    ss_run << "./analysis/" << algoname << "_" << filelabel << "_runerror.csv";
    ofstream writecsv;
    writecsv.open(ss_run.str(), ios::app);
    writecsv << avg_l1_error << ',' << avg_pre50 << ',' << avg_conductance << ',' << avg_max_additive_error << endl;
}

metric::~metric()
{
}

int main(int argc, char **argv)
{
    long i = 1;
    char *endptr;
    string filedir, filelabel, algoname;
    int querynum;
    double eps;
    uint L;
    vector<double> epss;
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
    g = Graph(filedir, filelabel);
    uint query_node;
    ifstream query_file;
    query_file.open("./query/" + filelabel + ".query");

    while (query_file >> query_node)
    {
        query_set.push_back(query_node);
    }
    for (vector<double>::iterator it = epss.begin(); it != epss.end(); it++)
        metric m(filedir, filelabel, algoname, querynum, L, *it);
    cout << "finish at ";
    time_t t = time(0); // get time now
    tm *now = localtime(&t);
    cout << (now->tm_year + 1900) << '-'
         << (now->tm_mon + 1) << '-'
         << now->tm_mday
         << "\n";
    return 0;
}
