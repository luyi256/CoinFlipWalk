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
    string neiNode, neiWeight, neiNum, graphAttr;
    neiNode = filedir + filelabel + ".outEdges";
    neiWeight = filedir + filelabel + ".outWEdges";
    neiNum = filedir + filelabel + ".outPtr";
    graphAttr = filedir + filelabel + ".attribute";
    cout << "FilePath: " << graphAttr.c_str() << endl;
    cout << "Read graph attributes..." << endl;
    string tmp;
    uint n;
    ifstream graphAttrIn(graphAttr.c_str());
    graphAttrIn >> tmp >> n;
    cout << "n=" << n << endl;
    graphAttrIn.close();
    cout << "Read graph ..." << endl;
    ifstream neiNumIn(neiNum.c_str(), ios::in | ios::binary);
    ifstream neiWeightIn(neiWeight.c_str(), ios::in | ios::binary);
    ifstream neiNodeIn(neiNode.c_str(), ios::in | ios::binary);
    uint outSize;
    uint outSizeSum = 0, preOutSizeSum = 0;
    neiNumIn.read(reinterpret_cast<char *>(&outSizeSum), sizeof(uint));
    uint totalEdge = 0;
    double totalDeg = 0;
    double product = 0;
    for (uint i = 0; i < n; i++)
    {
        neiNumIn.read(reinterpret_cast<char *>(&outSizeSum), sizeof(uint));
        outSize = outSizeSum - preOutSizeSum;
        preOutSizeSum = outSizeSum;
        double maxw = 0;
        for (uint j = 0; j < outSize; j++)
        {
            uint id;
            double w;
            neiNodeIn.read(reinterpret_cast<char *>(&id), sizeof(uint));
            neiWeightIn.read(reinterpret_cast<char *>(&w), sizeof(double));
            if (w > maxw)
                maxw = w;
            totalEdge++;
            totalDeg += w;
            product += sqrt(w);
        }
    }
    double x_norm = sqrt(totalEdge);
    double y_norm = sqrt(totalDeg);
    cout << "cos^2\psi: " << (product / (x_norm * y_norm)) * (product / (x_norm * y_norm));
    neiNumIn.close();
    neiWeightIn.close();
    neiNodeIn.close();
}
