#include <algorithm>
#include <iostream>
#include <string>
#include <sstream>
#include <cstdlib>
#include <vector>
#include <fstream>
#include <cstring>
#include <unordered_map>
#include "Random.h"
#include <random>
using namespace std;

vector<uint> query_set;

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
    Random R;
    unordered_map<uint, vector<node>> neighborList;
    unordered_map<uint, uint> outSizeList;
    unordered_map<uint, double> outWeightList;
    double totdeg;
    Graph() {}
    Graph(const string &_filedir, const string &_filelabel)
    {
        R = Random();
        totdeg = 0;
        filedir = _filedir;
        filelabel = _filelabel;
        string neiNode, neiWeight, neiNum, graphAttr;
        neiNode = filedir + filelabel + ".outEdges";
        neiWeight = filedir + filelabel + ".outWEdges";
        neiNum = filedir + filelabel + ".outPtr";
        graphAttr = filedir + filelabel + ".attribute";
        readFile(graphAttr, neiNode, neiWeight, neiNum);
    }

    ~Graph()
    {
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
            for (uint j = 0; j < outSize; j++)
            {
                uint id;
                double w;
                neiNodeIn.read(reinterpret_cast<char *>(&id), sizeof(uint));
                neiWeightIn.read(reinterpret_cast<char *>(&w), sizeof(double));
                neighborList[i].push_back(node(id, w));
                outWeight += w;
            }
            outSizeList[i] = outSize;
            outWeightList[i] = outWeight;
            totdeg += outWeight;
        }
        neiNumIn.close();
        neiWeightIn.close();
        neiNodeIn.close();
    }
};
Graph g;

class metric
{
private:
    uint vert;
    double *gtvalues;
    vector<pair<double, uint>> unnm_algoanswers;
    vector<pair<double, int>> unnm_gtanswers;
    uint k;
    double *algovalues;
    bool *algocheck;
    vector<uint> topk_algo_Nodes;

public:
    double cal_maxAE();
    double cal_l1Error();
    double calPrecision();
    double conductance();
    metric(string filedir, string filelabel, string algoname, long querynum, uint L, double eps, double thetad);
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
    return (hitCount / (double)k);
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
metric::metric(string filedir, string filelabel, string algoname, long querynum, uint L, double eps, double thetad)
{
    k = 50;
    double avg_l1_error = 0, avg_pre50 = 0, avg_conductance = 0, avg_max_additive_error = 0;
    vert = g.n;

    for (uint i = 0; i < querynum; ++i)
    {
        uint u = query_set[i];
        stringstream ss_gt, ss_gt_dir;
        ss_gt << "/home/lu_yi/result/powermethod/" << filelabel << "/" << L << "/" << u << "_gt.txt";
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
        cout << ss_algo.str() << endl; // hanzhi
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

        sort(unnm_gtanswers.begin(), unnm_gtanswers.end(), greater<pair<double, uint>>());
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
        sort(unnm_algoanswers.begin(), unnm_algoanswers.end(), greater<pair<double, uint>>());
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
    cout << "Avg precision@50 = " << avg_pre50 << endl;
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
    double thetad;
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
                epss.push_back(eps);
            }
            i--;
        }
        else if (!strcmp(argv[i], "-t"))
        {
            if (++i < argc)
            {
                thetad = strtod(argv[i], &endptr);
            }
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
        metric m(filedir, filelabel, algoname, querynum, L, *it, thetad);
    cout << "finish at ";
    time_t t = time(0); // get time now
    tm *now = localtime(&t);
    cout << (now->tm_year + 1900) << '-'
         << (now->tm_mon + 1) << '-'
         << now->tm_mday
         << "\n";
    return 0;
}
