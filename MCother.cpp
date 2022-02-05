#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <cstring>
#include <cmath>
#include <sys/time.h>
#include <unordered_map>
#include "Graphother.h"
using namespace std;
typedef unsigned int uint;

int main(int argc, char **argv)
{
    long i = 1;
    string filedir, filelabel;
    char *endptr;
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
        else
        {
            cerr << "Err args:" << argv[i] << endl;
            exit(1);
        }
        i++;
    }
    AliasMethodGraph g = AliasMethodGraph(filedir, filelabel);
    ifstream opfile(filedir + "/" + filelabel + ".op", ios::in);
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
    cout << filelabel << " avg update time: " << tottime / opnum << endl;
    opfile.close();
}
